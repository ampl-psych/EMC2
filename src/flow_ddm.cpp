// Derived from ~/Downloads/julian/nle/handover/nle-methods/port/src/flow_ddm.cpp
// (sha256 0a9ea58e9416b40699271628b29be2472bbf86647d2fdadafc27fa9882bd560b),
// the handover's source of truth, copied verbatim at d2202060. dev-nle Phase 2
// replaced its per-row conditioning with the batched evaluator in nle_flow.h
// (Armadillo GEMM per layer, distinct-row index, knot blocks); the model
// definition, the (SZ0, SZ, SZ1) reparameterisation and the per-trial spline
// arithmetic are unchanged. Agreement with the R port (port/R/flow_ddm.R) is
// tested to 1e-10 and with the golden vectors to 1e-9.
//
// C++ (Rcpp) evaluator for the exported DDM flow + choice classifier
// (DDM-character models: spline flow over rt conditioned on
// [standardized params, raw response code], plus a classifier MLP giving
// the logit of P(R=2); see flow/src/flow/config.py `joint_log_prob`).
//
// API:
//   ddm_build(fl)                        -> external pointer
//   ddm_eval_cpp(ptr, theta, rt, R)      -> one parameter vector, many
//                                           (rt, R) pairs
//   ddm_eval_trials_cpp(ptr, Theta, rt, R)
//                                        -> trial-wise rows; each distinct row
//                                           is conditioned once, all distinct
//                                           rows in batched network passes
// Ensemble (seed-ensembled production flows; equal-weight density mixture):
//   ddm_build_ensemble(list of fl)       -> external pointer to k members
//   ddm_ens_eval_trials_cpp(ptr, Theta, rt, R)
//                                        -> trial-wise mixture pdf/cdf
//                                           (mean over members); the batches
//                                           run member by member
//
// Outputs: joint (defective) pdf and cdf, log pdf, and P(R | params).
// Out-of-box parameters: pdf = 0, cdf = 0, log_pdf = -Inf, p_R = NA
// (rejection sentinel, not extrapolation).

#include "nle_flow.h"
using namespace Rcpp;

struct DdmModel {
  int n_ctx = 0;                          // number of parameters (no R)
  nle::Spline sp;
  std::vector<double> scaler_mean, scaler_scale;
  // Under the (SZ0, SZ, SZ1) reparameterisation the flow and the frozen
  // classifier take different encodings of the same draw, so the classifier
  // carries its own scaler. theta is ALWAYS supplied on EMC2's scale.
  bool triple = false;
  std::vector<double> clf_scaler_mean, clf_scaler_scale;
  std::vector<double> lower, upper;       // bounding box (sampled scale)
  nle::Mlp flow;                          // input dim n_ctx + 1 (R appended)
  nle::Mlp clf;                           // input dim n_ctx, output 1 logit
};

struct DdmEnsemble { std::vector<DdmModel> members; };

// EMC2 (v, log a, log t0, log s, qnorm Z, qnorm SZ, log sv, log st0) ->
// (v, log SZ0, log t0, log s, log SZ, log SZ1, log sv, log st0). Mirrors
// .emc2_to_triple() in port/R/flow_ddm.R and emc2_to_triple() in
// flow/scripts/export_flow_ddm.py. SZ0, SZ1 >= 0 always holds because
// sw <= 2*min(w, 1-w); the clamp only guards float cancellation.
static inline double std_norm_cdf(double x) {
  return 0.5 * std::erfc(-x * M_SQRT1_2);
}
static void emc2_to_triple(const double* theta, int n, double* out) {
  std::copy(theta, theta + n, out);
  const double a  = std::exp(theta[1]);
  const double w  = std_norm_cdf(theta[4]);
  const double sw = 2.0 * std_norm_cdf(theta[5]) * std::min(w, 1.0 - w);
  const double SZ  = sw * a;
  const double SZ0 = a * w - SZ / 2.0;
  const double SZ1 = a - SZ0 - SZ;
  const double eps = 1e-12;
  out[1] = std::log(std::max(SZ0, eps));
  out[4] = std::log(std::max(SZ,  eps));
  out[5] = std::log(std::max(SZ1, eps));
}

struct DdmOut { double *pdf, *cdf, *log_pdf, *p_R; };

// Trials t = 0..n-1 with parameter row uid[t] of Theta_u (n_ctx x U, one
// distinct row per column), rt[t] and R[t] in {1, 2}; equal-weight mixture over
// `mem` (one member = the flow itself). Every in-box distinct row runs the
// classifier once and the flow once per response it is observed with; both
// networks run on blocks of rows.
static void ddm_eval_core(const std::vector<const DdmModel*>& mem,
                          const double* Theta_u, int U, const int* uid,
                          const double* rt, const int* R, int n, DdmOut out) {
  const DdmModel& m0 = *mem[0];
  const int nc = m0.n_ctx;
  const int n_mem = (int)mem.size();

  // distinct rows in the box (members share the box)
  std::vector<int> inb_id(U, -1), inb;              // row u -> in-box index
  for (int u = 0; u < U; ++u)
    if (nle::in_box(m0.lower, m0.upper, Theta_u + (size_t)u * nc)) {
      inb_id[u] = (int)inb.size(); inb.push_back(u);
    }
  const int Uin = (int)inb.size();

  // flow columns: the (row, response) pairs that occur
  std::vector<int> fkey((size_t)2 * U, -1);
  for (int t = 0; t < n; ++t)
    if (inb_id[uid[t]] >= 0) fkey[(size_t)2 * uid[t] + R[t] - 1] = 0;
  std::vector<int> f_u, f_r;
  for (int u = 0; u < U; ++u)
    for (int r = 0; r < 2; ++r)
      if (fkey[(size_t)2 * u + r] == 0) {
        fkey[(size_t)2 * u + r] = (int)f_u.size(); f_u.push_back(u); f_r.push_back(r);
      }
  const int F = (int)f_u.size();
  std::vector<int> trial_f(n, -1);
  for (int t = 0; t < n; ++t) {
    if (inb_id[uid[t]] >= 0) trial_f[t] = fkey[(size_t)2 * uid[t] + R[t] - 1];
    else { out.pdf[t] = 0.0; out.cdf[t] = 0.0; out.log_pdf[t] = R_NegInf; out.p_R[t] = NA_REAL; }
  }
  const nle::Groups grp = nle::group_by(trial_f, F);

  std::vector<double> lp_mem;                       // n x n_mem (mixtures only)
  if (n_mem > 1) {
    lp_mem.assign((size_t)n * n_mem, R_NegInf);
    for (int t = 0; t < n; ++t)
      if (trial_f[t] >= 0) { out.pdf[t] = 0.0; out.cdf[t] = 0.0; out.p_R[t] = 0.0; }
  }

  std::vector<double> lpR((size_t)2 * Uin);         // log P(R = 1), log P(R = 2)
  std::vector<double> tf(nc);
  nle::KnotBlock kb;
  for (int k = 0; k < n_mem; ++k) {
    const DdmModel& m = *mem[k];
    const std::vector<double>& cmean = m.triple ? m.clf_scaler_mean : m.scaler_mean;
    const std::vector<double>& cscale = m.triple ? m.clf_scaler_scale : m.scaler_scale;

    // classifier: one pass per block of in-box rows
    for (int c0 = 0; c0 < Uin; c0 += nle::NLE_BLOCK) {
      const int nb = std::min(nle::NLE_BLOCK, Uin - c0);
      arma::mat X(nc, nb);
      for (int c = 0; c < nb; ++c) {
        const double* th = Theta_u + (size_t)inb[c0 + c] * nc;
        for (int j = 0; j < nc; ++j) X(j, c) = (th[j] - cmean[j]) / cscale[j];
      }
      const arma::mat logit = nle::mlp_forward_batch(m.clf, X);
      for (int c = 0; c < nb; ++c) {
        lpR[(size_t)2 * (c0 + c)]     = -nle::softplus(logit(0, c));    // log P(R = 1)
        lpR[(size_t)2 * (c0 + c) + 1] = -nle::softplus(-logit(0, c));   // log P(R = 2)
      }
    }

    // flow: one pass per block of (row, response) columns, then the trials
    const int K = m.sp.num_bins;
    for (int c0 = 0; c0 < F; c0 += nle::NLE_BLOCK) {
      const int nb = std::min(nle::NLE_BLOCK, F - c0);
      arma::mat X(nc + 1, nb);
      for (int c = 0; c < nb; ++c) {
        const double* th = Theta_u + (size_t)f_u[c0 + c] * nc;
        if (m.triple) { emc2_to_triple(th, nc, tf.data()); th = tf.data(); }
        for (int j = 0; j < nc; ++j) X(j, c) = (th[j] - m.scaler_mean[j]) / m.scaler_scale[j];
        X(nc, c) = (double)(f_r[c0 + c] + 1);       // raw response code
      }
      const arma::mat raw = nle::mlp_forward_batch(m.flow, X);
      kb.reset(K, nb);
      for (int c = 0; c < nb; ++c) nle::build_knots(m.sp, raw.colptr(c), kb, c);
      for (int c = 0; c < nb; ++c) {
        const int f = c0 + c;
        const double lp_R = lpR[(size_t)2 * inb_id[f_u[f]] + f_r[f]];
        const double p_R = std::exp(lp_R);
        for (int gi = grp.start[f]; gi < grp.start[f + 1]; ++gi) {
          const int t = grp.order[gi];
          const double u = std::log(rt[t]);
          double z, logdet;
          nle::rqs_inverse_one(kb.xc(c), kb.yc(c), kb.dc(c), K, u, z, logdet);
          const double lp = R::dnorm(z, 0.0, 1.0, 1) + logdet - u + lp_R;
          const double cdf = R::pnorm(z, 0.0, 1.0, 1, 0) * p_R;
          if (n_mem == 1) {
            out.log_pdf[t] = lp; out.pdf[t] = std::exp(lp); out.cdf[t] = cdf; out.p_R[t] = p_R;
          } else {
            lp_mem[(size_t)t * n_mem + k] = lp;
            out.pdf[t] += std::exp(lp); out.cdf[t] += cdf; out.p_R[t] += p_R;
          }
        }
      }
    }
  }

  if (n_mem > 1) {                                  // mean density; log via log-sum-exp
    const double log_k = std::log((double)n_mem);
    for (int t = 0; t < n; ++t) {
      if (trial_f[t] < 0) continue;
      out.pdf[t] /= n_mem; out.cdf[t] /= n_mem; out.p_R[t] /= n_mem;
      const double* lp = lp_mem.data() + (size_t)t * n_mem;
      double mx = lp[0];
      for (int k = 1; k < n_mem; ++k) mx = std::fmax(mx, lp[k]);
      if (!std::isfinite(mx)) { out.log_pdf[t] = std::log(out.pdf[t]); continue; }
      double s = 0.0;
      for (int k = 0; k < n_mem; ++k) s += std::exp(lp[k] - mx);
      out.log_pdf[t] = mx + std::log(s) - log_k;
    }
  }
}

// Trial-wise entry point shared by the single-model and ensemble exports.
static List ddm_eval_trials(const std::vector<const DdmModel*>& mem,
                            NumericMatrix theta, NumericVector rt, IntegerVector R) {
  const int nc = mem[0]->n_ctx;
  const int n = rt.size();
  if (theta.nrow() != n || theta.ncol() != nc)
    stop("theta must be length(rt) x %d.", nc);
  if (R.size() != n) stop("rt and R must have equal length.");
  for (int t = 0; t < n; ++t)
    if (R[t] != 1 && R[t] != 2) stop("R must be 1 or 2.");
  NumericVector pdf(n), cdf(n), log_pdf(n), p_R(n);
  const nle::RowIndex ix = nle::index_rows(theta.begin(), n, nc);
  const int U = ix.U();
  std::vector<double> Theta_u((size_t)U * nc);
  for (int u = 0; u < U; ++u)
    for (int j = 0; j < nc; ++j) Theta_u[(size_t)u * nc + j] = theta(ix.first[u], j);
  ddm_eval_core(mem, Theta_u.data(), U, ix.uid.data(), rt.begin(), R.begin(), n,
                DdmOut{pdf.begin(), cdf.begin(), log_pdf.begin(), p_R.begin()});
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                      _["log_pdf"] = log_pdf, _["p_R"] = p_R);
}

// ---------------------------------------------------------------------------
// R interface
// ---------------------------------------------------------------------------

static void load_ddm_from_list(List fl, DdmModel& m) {
  nle::load_spline(fl["spline"], m.sp, "ddm_build");
  List scaler = fl["scaler"];
  m.scaler_mean = as<std::vector<double>>(scaler["mean"]);
  m.scaler_scale = as<std::vector<double>>(scaler["scale"]);
  m.n_ctx = (int)m.scaler_mean.size();
  List bounds = fl["bounds_sampled"];
  m.lower = as<std::vector<double>>(bounds["lower"]);
  m.upper = as<std::vector<double>>(bounds["upper"]);
  if ((int)m.lower.size() != m.n_ctx || (int)m.upper.size() != m.n_ctx)
    stop("ddm_build: bounds_sampled length != n_params.");
  nle::load_mlp(fl["flow_mlp"], m.flow, "ddm_build");
  if (fl.containsElementNamed("context_encoding")) {
    std::string enc = as<std::string>(fl["context_encoding"]);
    m.triple = (enc == "triple");
  }
  if (m.triple) {
    if (!fl.containsElementNamed("classifier_scaler"))
      stop("ddm_build: context_encoding is 'triple' but classifier_scaler is absent.");
    List cs = fl["classifier_scaler"];
    m.clf_scaler_mean  = as<std::vector<double>>(cs["mean"]);
    m.clf_scaler_scale = as<std::vector<double>>(cs["scale"]);
    if ((int)m.clf_scaler_mean.size() != m.n_ctx)
      stop("ddm_build: classifier_scaler length != n_params.");
  }
  nle::load_mlp(fl["classifier_mlp"], m.clf, "ddm_build");
  if (m.flow.n_in() != m.n_ctx + 1)
    stop("ddm_build: flow input dim must be n_params + 1 (response).");
  if (m.flow.n_out() != 3 * m.sp.num_bins + 1)
    stop("ddm_build: flow output dim must be 3 * num_bins + 1.");
  if (m.clf.n_in() != m.n_ctx)
    stop("ddm_build: classifier input dim must be n_params.");
  if (m.clf.n_out() != 1)
    stop("ddm_build: classifier must have one output (the logit of P(R = 2)).");
}

// [[Rcpp::export]]
SEXP ddm_build(List fl) {
  DdmModel* m = new DdmModel();
  try { load_ddm_from_list(fl, *m); } catch (...) { delete m; throw; }
  XPtr<DdmModel> ptr(m, true);
  ptr.attr("class") = "ddm_flow_model";
  if (fl.containsElementNamed("model"))       // a label only; cards without an
    ptr.attr("model") = as<std::string>(fl["model"]);   // analytic model omit it
  return ptr;
}

// Equal-weight ensemble of DDM flows (same context/bounds, e.g. training
// seeds). fls: an R list of weight lists as accepted by ddm_build.
// [[Rcpp::export]]
SEXP ddm_build_ensemble(List fls) {
  if (fls.size() < 1) stop("ddm_build_ensemble: need at least one member.");
  DdmEnsemble* e = new DdmEnsemble();
  try {
    e->members.resize(fls.size());
    for (int k = 0; k < fls.size(); ++k) {
      load_ddm_from_list(fls[k], e->members[k]);
      if (e->members[k].n_ctx != e->members[0].n_ctx ||
          e->members[k].lower != e->members[0].lower ||
          e->members[k].upper != e->members[0].upper)
        stop("ddm_build_ensemble: members disagree on context dim or bounds.");
    }
  } catch (...) { delete e; throw; }
  XPtr<DdmEnsemble> ptr(e, true);
  ptr.attr("class") = "ddm_flow_ensemble";
  ptr.attr("n_members") = (int)e->members.size();
  return ptr;
}

// [[Rcpp::export]]
bool ddm_ptr_valid(SEXP ptr_) {
  return TYPEOF(ptr_) == EXTPTRSXP && R_ExternalPtrAddr(ptr_) != nullptr;
}

// One parameter vector (sampled scale), many (rt, R) pairs. Both networks
// run exactly once — the amortized path for MCMC likelihoods.
// [[Rcpp::export]]
List ddm_eval_cpp(SEXP ptr_, NumericVector theta, NumericVector rt,
                  IntegerVector R) {
  XPtr<DdmModel> ptr(ptr_);
  const DdmModel& m = *ptr;
  if ((int)theta.size() != m.n_ctx)
    stop("theta has %d elements; model expects %d.", theta.size(), m.n_ctx);
  const int n = rt.size();
  if (R.size() != n) stop("rt and R must have equal length.");
  for (int t = 0; t < n; ++t)
    if (R[t] != 1 && R[t] != 2) stop("R must be 1 or 2.");
  NumericVector pdf(n), cdf(n), log_pdf(n), p_R(n);
  std::vector<int> uid(n, 0);
  ddm_eval_core({&m}, theta.begin(), 1, uid.data(), rt.begin(), R.begin(), n,
                DdmOut{pdf.begin(), cdf.begin(), log_pdf.begin(), p_R.begin()});
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                      _["log_pdf"] = log_pdf, _["p_R"] = p_R,
                      _["in_box"] = nle::in_box(m.lower, m.upper, theta.begin()));
}

// Trial-wise parameter rows (n x K_ctx) with rt (n) and R (n). Each distinct
// row is conditioned once (both networks, batched over the distinct rows).
// [[Rcpp::export]]
List ddm_eval_trials_cpp(SEXP ptr_, NumericMatrix theta, NumericVector rt,
                         IntegerVector R) {
  XPtr<DdmModel> ptr(ptr_);
  return ddm_eval_trials({ptr.get()}, theta, rt, R);
}

// Trial-wise ensemble evaluation: equal-weight mixture over members.
// pdf/cdf/p_R are means over members; log_pdf = log(mean pdf), computed by
// log-sum-exp over the members' log densities.
// [[Rcpp::export]]
List ddm_ens_eval_trials_cpp(SEXP ptr_, NumericMatrix theta, NumericVector rt,
                             IntegerVector R) {
  XPtr<DdmEnsemble> ptr(ptr_);
  std::vector<const DdmModel*> mem;
  for (const DdmModel& m : ptr->members) mem.push_back(&m);
  return ddm_eval_trials(mem, theta, rt, R);
}
