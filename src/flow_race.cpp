// Derived from ~/Downloads/julian/nle/handover/nle-methods/port/src/flow_race.cpp
// (sha256 67701a892a5e20cf22cda0359c955cb74d324b526682d2971e9830804a11eda1),
// the handover's source of truth, copied verbatim at d2202060. dev-nle Phase 2
// replaced its per-row conditioning with the batched evaluator in nle_flow.h
// (Armadillo GEMM per layer, distinct-row index, knot blocks); the per-trial
// spline arithmetic is unchanged. Agreement with the R port
// (port/R/flow_race.R) is tested to 1e-10 and with the golden vectors to 1e-9.
//
// C++ (Rcpp) evaluator for exported race-model flows (LNR, RDM, ...).
// Model-agnostic: any context dimension, any MLP depth/widths, optional
// LayerNorm, any bin count. Only single-spline chains with an exp output
// transform and a standard-normal base are supported (what
// flow.config.spline_flow builds).
//
// Design for EMC2 race models:
//   flow_build(fl)                 -> external pointer, built once per session
//   flow_eval_cpp(ptr, theta, rt)  -> one parameter vector, many rts
//                                     (conditioner runs once: the MCMC case)
//   flow_eval_trials_cpp(ptr, Theta, rt)
//                                  -> trial-wise parameter rows; each distinct
//                                     row is conditioned once, all distinct
//                                     rows in one batched conditioner pass
//   nle_flow_native()              -> native calc_ll branch (nle_native.h,
//                                     model_NN.h): natural-scale rows; log pdf
//                                     for winners, log survivor for losers
//
// Optional card features (cards converted from Malte Lueken's flows; a card
// without them evaluates exactly as before):
//   input_scaling  the network sees (log(x + eps) - loc) / scale per input, in
//                  place of the scaler; with neither it sees the inputs as
//                  they are
//   affine         two more conditioner outputs after the 3K + 1 spline
//                  outputs, loc and raw_scale; with scale = softplus(raw_scale
//                  + offset) + min_scale and y = (log t - loc) / scale the
//                  spline is inverted at y (identity tails in y) and
//                  log pdf = dnorm(z, log) + log|dz/dy| - log(scale) - log t
//
// Out-of-box parameters (training-region bounding box, sampled scale):
// pdf = 0, cdf = 0 immediately (log_pdf = -Inf, log_sf = 0). Rejection
// sentinel, not extrapolation.
//
// Uses Rmath dnorm/pnorm so results are comparable with the R reference.

#include "nle_flow.h"
#include "nle_native.h"
#include "nle_wald.h"
using namespace Rcpp;

struct FlowModel {
  int n_ctx = 0;
  nle::Spline sp;
  int input_mode = 0;                                  // 0 as is, 1 scaler, 2 input_scaling
  std::vector<double> scaler_mean, scaler_scale;       // 1: (x - mean) / scale
  std::vector<double> in_eps, in_loc, in_scale;        // 2: (log(x + eps) - loc) / scale
  bool affine = false;                                 // affine stage after the spline
  double aff_min_scale = 0, aff_offset = 0;
  std::vector<double> lower, upper;                    // bounding box (sampled scale)
  nle::Mlp mlp;
};

// Any pointer may be null (that output is not computed). With `want`, trial
// t gets only its density outputs (pdf, log_pdf) if want[t] == 1 and only its
// distribution outputs (cdf, log_sf) if want[t] == 2: race winners need the
// density, losers the survivor function.
struct RaceOut { double *pdf, *cdf, *log_pdf, *log_sf; const unsigned char* want; };

// Trials t = 0..n-1 with parameter row uid[t] of Theta_u (n_ctx x U, one
// distinct row per column) and rt[t]. The conditioner runs once per in-box
// distinct row, on blocks of rows.
static void flow_eval_core(const FlowModel& m, const double* Theta_u, int U,
                           const int* uid, const double* rt, int n, RaceOut out) {
  auto want_d = [&](int t) { return out.want == nullptr || out.want[t] == 1; };
  auto want_p = [&](int t) { return out.want == nullptr || out.want[t] == 2; };
  const int nc = m.n_ctx;
  const int K = m.sp.num_bins;
  std::vector<int> inb_id(U, -1), inb;
  for (int u = 0; u < U; ++u)
    if (nle::in_box(m.lower, m.upper, Theta_u + (size_t)u * nc)) {
      inb_id[u] = (int)inb.size(); inb.push_back(u);
    }
  const int Uin = (int)inb.size();
  std::vector<int> trial_c(n, -1);
  for (int t = 0; t < n; ++t) {
    trial_c[t] = inb_id[uid[t]];
    if (trial_c[t] >= 0) continue;
    if (want_d(t)) {
      if (out.pdf) out.pdf[t] = 0.0;
      if (out.log_pdf) out.log_pdf[t] = R_NegInf;
    }
    if (want_p(t)) {
      if (out.cdf) out.cdf[t] = 0.0;
      if (out.log_sf) out.log_sf[t] = 0.0;
    }
  }
  const nle::Groups grp = nle::group_by(trial_c, Uin);

  nle::KnotBlock kb;
  for (int c0 = 0; c0 < Uin; c0 += nle::NLE_BLOCK) {
    const int nb = std::min(nle::NLE_BLOCK, Uin - c0);
    arma::mat X(nc, nb);
    for (int c = 0; c < nb; ++c) {
      const double* th = Theta_u + (size_t)inb[c0 + c] * nc;
      for (int j = 0; j < nc; ++j)
        X(j, c) = m.input_mode == 1 ? (th[j] - m.scaler_mean[j]) / m.scaler_scale[j] :
                  m.input_mode == 2 ? (std::log(th[j] + m.in_eps[j]) - m.in_loc[j]) / m.in_scale[j] : th[j];
    }
    const arma::mat raw = nle::mlp_forward_batch(m.mlp, X);
    kb.reset(K, nb);
    for (int c = 0; c < nb; ++c) nle::build_knots(m.sp, raw.colptr(c), kb, c);
    for (int c = 0; c < nb; ++c) {
      const int q = c0 + c;
      double loc = 0.0, scale = 1.0, log_scale = 0.0;
      if (m.affine) {
        const double* r = raw.colptr(c);
        loc = r[3 * K + 1];
        scale = nle::softplus(r[3 * K + 2] + m.aff_offset) + m.aff_min_scale;
        log_scale = std::log(scale);
      }
      for (int gi = grp.start[q]; gi < grp.start[q + 1]; ++gi) {
        const int t = grp.order[gi];
        const double u = std::log(rt[t]);
        double z, logdet;
        if (m.affine) {
          nle::rqs_inverse_one(kb.xc(c), kb.yc(c), kb.dc(c), K, (u - loc) / scale, z, logdet);
          logdet -= log_scale;
        } else nle::rqs_inverse_one(kb.xc(c), kb.yc(c), kb.dc(c), K, u, z, logdet);
        if (want_d(t)) {
          const double lp = R::dnorm(z, 0.0, 1.0, 1) + logdet - u;
          if (out.log_pdf) out.log_pdf[t] = lp;
          if (out.pdf) out.pdf[t] = std::exp(lp);
        }
        if (want_p(t)) {
          if (out.cdf) out.cdf[t] = R::pnorm(z, 0.0, 1.0, 1, 0);
          if (out.log_sf) out.log_sf[t] = R::pnorm(z, 0.0, 1.0, 0, 1);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// R interface
// ---------------------------------------------------------------------------

// TRUE if the external pointer still holds a live model (XPtrs read back
// from a saved session deserialize as NULL and must be rebuilt).
// [[Rcpp::export]]
bool flow_ptr_valid(SEXP ptr_) {
  return TYPEOF(ptr_) == EXTPTRSXP && R_ExternalPtrAddr(ptr_) != nullptr;
}

// [[Rcpp::export]]
SEXP flow_build(List fl) {
  FlowModel* m = new FlowModel();
  try {
    nle::load_spline(fl["spline"], m->sp, "flow_build");
    nle::load_mlp(fl["mlp"], m->mlp, "flow_build");
    m->n_ctx = m->mlp.n_in();
    auto has = [&](const char* f) { return fl.containsElementNamed(f) && !Rf_isNull(fl[f]); };
    if (has("scaler") && has("input_scaling"))
      stop("flow_build: the card has both a scaler and input_scaling.");
    if (has("scaler")) {
      List scaler = fl["scaler"];
      m->scaler_mean = as<std::vector<double>>(scaler["mean"]);
      m->scaler_scale = as<std::vector<double>>(scaler["scale"]);
      if ((int)m->scaler_mean.size() != m->n_ctx || (int)m->scaler_scale.size() != m->n_ctx)
        stop("flow_build: the scaler needs one mean and scale per MLP input.");
      m->input_mode = 1;
    } else if (has("input_scaling")) {
      List isc = fl["input_scaling"];
      m->in_eps = as<std::vector<double>>(isc["eps"]);
      m->in_loc = as<std::vector<double>>(isc["loc"]);
      m->in_scale = as<std::vector<double>>(isc["scale"]);
      if ((int)m->in_eps.size() != m->n_ctx || (int)m->in_loc.size() != m->n_ctx ||
          (int)m->in_scale.size() != m->n_ctx)
        stop("flow_build: input_scaling needs one eps, loc and scale per MLP input.");
      m->input_mode = 2;
    }
    m->affine = has("affine") && as<bool>(fl["affine"]);
    if (m->affine) {
      if (!has("affine_scale")) stop("flow_build: an affine card needs affine_scale (min_scale, offset).");
      List as_ = fl["affine_scale"];
      m->aff_min_scale = as<double>(as_["min_scale"]);
      m->aff_offset = as<double>(as_["offset"]);
    }
    if (!has("bounds_sampled")) stop("flow_build: the card has no bounds_sampled.");
    List bounds = fl["bounds_sampled"];
    m->lower = as<std::vector<double>>(bounds["lower"]);
    m->upper = as<std::vector<double>>(bounds["upper"]);
    if ((int)m->lower.size() != m->n_ctx || (int)m->upper.size() != m->n_ctx)
      stop("flow_build: bounds_sampled length != n_params.");
    if (m->mlp.n_out() != 3 * m->sp.num_bins + (m->affine ? 3 : 1))
      stop("flow_build: MLP output dim must be 3 * num_bins + %d.", m->affine ? 3 : 1);
  } catch (...) { delete m; throw; }
  XPtr<FlowModel> ptr(m, true);
  ptr.attr("class") = "flow_model";
  if (fl.containsElementNamed("model"))       // a label only; cards without an
    ptr.attr("model") = as<std::string>(fl["model"]);   // analytic model omit it
  return ptr;
}

// One parameter vector (sampled scale), many rts. The conditioner (MLP)
// runs exactly once — the amortized path for MCMC likelihoods.
// [[Rcpp::export]]
List flow_eval_cpp(SEXP ptr_, NumericVector theta, NumericVector rt) {
  XPtr<FlowModel> ptr(ptr_);
  const FlowModel& m = *ptr;
  if ((int)theta.size() != m.n_ctx)
    stop("theta has %d elements; model expects %d.", theta.size(), m.n_ctx);
  const int n = rt.size();
  NumericVector pdf(n), cdf(n), log_pdf(n), log_sf(n);
  std::vector<int> uid(n, 0);
  flow_eval_core(m, theta.begin(), 1, uid.data(), rt.begin(), n,
                 RaceOut{pdf.begin(), cdf.begin(), log_pdf.begin(), log_sf.begin(), nullptr});
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                      _["log_pdf"] = log_pdf, _["log_sf"] = log_sf,
                      _["in_box"] = nle::in_box(m.lower, m.upper, theta.begin()));
}

// Trial-wise parameter rows (n x K_ctx) paired with rt (n) — the shape of
// EMC2 dfun/pfun inputs. Each distinct row is conditioned once, whether or not
// its trials are adjacent (race winners/losers interleave accumulators).
// [[Rcpp::export]]
List flow_eval_trials_cpp(SEXP ptr_, NumericMatrix theta, NumericVector rt) {
  XPtr<FlowModel> ptr(ptr_);
  const FlowModel& m = *ptr;
  const int nc = m.n_ctx;
  const int n = rt.size();
  if (theta.nrow() != n || theta.ncol() != nc)
    stop("theta must be length(rt) x %d.", nc);
  NumericVector pdf(n), cdf(n), log_pdf(n), log_sf(n);
  const nle::RowIndex ix = nle::index_rows(theta.begin(), n, nc);
  const int U = ix.U();
  std::vector<double> Theta_u((size_t)U * nc);
  for (int u = 0; u < U; ++u)
    for (int j = 0; j < nc; ++j) Theta_u[(size_t)u * nc + j] = theta(ix.first[u], j);
  flow_eval_core(m, Theta_u.data(), U, ix.uid.data(), rt.begin(), n,
                 RaceOut{pdf.begin(), cdf.begin(), log_pdf.begin(), log_sf.begin(), nullptr});
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                      _["log_pdf"] = log_pdf, _["log_sf"] = log_sf);
}

// Batched forward pass of a bare MLP (the conditioner layout above): X has one
// input vector per ROW; returns one output vector per row. Internal: tests the
// batched GEMM path (incl. LayerNorm, which no shipped artefact uses) against
// plain R matrix algebra, and is the building block for plain-MLP likelihoods.
// [[Rcpp::export]]
NumericMatrix nle_mlp_forward(List mlp, NumericMatrix X) {
  nle::Mlp m;
  nle::load_mlp(mlp, m, "nle_mlp_forward", true);
  if (X.ncol() != m.n_in()) stop("X must have %d columns.", m.n_in());
  const arma::mat Xt = arma::mat(X.begin(), X.nrow(), X.ncol(), false, true).t();
  const arma::mat Y = nle::mlp_forward_batch(m, Xt);
  return wrap(arma::mat(Y.t()));
}

// The Wald rows of a hybrid neural likelihood (nle_wald.h), R path: density and
// CDF at decision times t > 0, one parameter row per time.
// [[Rcpp::export]]
List nle_wald_cpp(NumericVector t, NumericVector v, NumericVector b, NumericVector s) {
  const int n = t.size();
  if (v.size() != n || b.size() != n || s.size() != n) stop("v, b and s need one value per t.");
  NumericVector pdf(n), cdf(n);
  for (int i = 0; i < n; ++i) {
    if (!(t[i] > 0.0)) continue;
    pdf[i] = std::exp(nle::wald_log_pdf(t[i], v[i], b[i], s[i]));
    cdf[i] = nle::wald_cdf(t[i], v[i], b[i], s[i]);
  }
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf);
}

// ---------------------------------------------------------------------------
// Native calc_ll branch (nle_native.h)
// ---------------------------------------------------------------------------

const FlowModel* nle_flow_handle(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || !Rf_inherits(ptr, "flow_model") ||
      R_ExternalPtrAddr(ptr) == nullptr)
    stop("Neural likelihood: the evaluator is not a live race flow.");
  return static_cast<const FlowModel*>(R_ExternalPtrAddr(ptr));
}

int nle_flow_n_ctx(const FlowModel* f) { return f->n_ctx; }

void nle_flow_native(const FlowModel* f, const double* th, int m, const int* tf,
                     const double* tfc, const double* tn, const unsigned char* want,
                     double* log_pdf, double* log_sf) {
  std::vector<double> Theta_u;
  const nle::RowIndex ix = nle::index_net_rows(th, m, f->n_ctx, tf, tfc, Theta_u);
  flow_eval_core(*f, Theta_u.data(), ix.U(), ix.uid.data(), tn, m,
                 RaceOut{nullptr, nullptr, log_pdf, log_sf, want});
}
