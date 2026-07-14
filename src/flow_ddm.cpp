// C++ (Rcpp) evaluator for the exported DDM flow + choice classifier
// (DDM-character models: spline flow over rt conditioned on
// [standardized params, raw response code], plus a classifier MLP giving
// the logit of P(R=2); see flow/src/flow/config.py `joint_log_prob`).
//
// Direct translation of the validated R reference (port/R/flow_ddm.R).
// Shares its numerical core with flow_race.cpp (duplicated here so each
// file stays self-contained for sourceCpp; keep the helpers in sync).
//
// API, mirroring flow_race.cpp:
//   ddm_build(fl)                        -> external pointer
//   ddm_eval_cpp(ptr, theta, rt, R)      -> one parameter vector, many
//                                           (rt, R) pairs; the two networks
//                                           run once (knots for both
//                                           responses + one classifier pass)
//   ddm_eval_trials_cpp(ptr, Theta, rt, R)
//                                        -> trial-wise rows with
//                                           consecutive-duplicate caching
// Ensemble (seed-ensembled production flows; equal-weight density mixture):
//   ddm_build_ensemble(list of fl)       -> external pointer to k members
//   ddm_ens_eval_trials_cpp(ptr, Theta, rt, R)
//                                        -> trial-wise mixture pdf/cdf
//                                           (mean over members), duplicate
//                                           caching of all k conditioned
//                                           states
//
// Outputs: joint (defective) pdf and cdf, their logs, and P(R | params).
// Out-of-box parameters: pdf = 0, cdf = 0, log_pdf = -Inf (rejection
// sentinel, not extrapolation).

#include <Rcpp.h>
#include <cmath>
#include <cfloat>
#include <vector>
using namespace Rcpp;

struct Mlp {
  struct Layer { int nin, nout; std::vector<double> W, b; };
  std::vector<Layer> layers;
  bool use_norm = false;
  struct Norm { std::vector<double> scale, bias; double eps; };
  std::vector<Norm> norms;
};

struct DdmModel {
  int n_ctx = 0;                          // number of parameters (no R)
  int num_bins = 0;
  double range_min = 0, range_max = 0, min_bin_size = 0, min_knot_slope = 0;
  std::vector<double> scaler_mean, scaler_scale;
  std::vector<double> lower, upper;       // bounding box (sampled scale)
  Mlp flow;                                // input dim n_ctx + 1 (R appended)
  Mlp clf;                                 // input dim n_ctx, output 1 logit
};

struct Knots { std::vector<double> x, y, d; };

// per-parameter-vector amortized state: knots for R=1/R=2 + log P(R)
struct DdmCond {
  Knots kn[2];
  double log_p_R[2];
};

static inline double gelu_tanh(double x) {
  const double c = 0.7978845608028654;  // sqrt(2/pi)
  return 0.5 * x * (1.0 + std::tanh(c * (x + 0.044715 * x * x * x)));
}

static inline double softplus(double x) {
  return std::fmax(x, 0.0) + std::log1p(std::exp(-std::fabs(x)));
}

static void mlp_forward(const Mlp& m, std::vector<double> h,
                        std::vector<double>& out) {
  std::vector<double> h2;
  const int L = (int)m.layers.size();
  for (int li = 0; li < L; ++li) {
    const Mlp::Layer& lay = m.layers[li];
    h2.assign(lay.nout, 0.0);
    for (int i = 0; i < lay.nin; ++i) {
      const double hi = h[i];
      for (int j = 0; j < lay.nout; ++j)
        h2[j] += hi * lay.W[(size_t)j * lay.nin + i];
    }
    for (int j = 0; j < lay.nout; ++j) h2[j] += lay.b[j];
    if (li < L - 1) {
      if (m.use_norm) {
        const Mlp::Norm& nrm = m.norms[li];
        double mu = 0.0;
        for (double v : h2) mu += v;
        mu /= lay.nout;
        double var = 0.0;
        for (double v : h2) var += (v - mu) * (v - mu);
        var /= lay.nout;
        const double inv_sd = 1.0 / std::sqrt(var + nrm.eps);
        for (int j = 0; j < lay.nout; ++j)
          h2[j] = (h2[j] - mu) * inv_sd * nrm.scale[j] + nrm.bias[j];
      }
      for (int j = 0; j < lay.nout; ++j) h2[j] = gelu_tanh(h2[j]);
    }
    h.swap(h2);
  }
  out.swap(h);
}

static void build_knots(const DdmModel& m, const std::vector<double>& raw,
                        Knots& kn) {
  const int K = m.num_bins;
  const double range_size = m.range_max - m.range_min;
  const double budget = range_size - K * m.min_bin_size;
  kn.x.resize(K + 1); kn.y.resize(K + 1); kn.d.resize(K + 1);
  for (int part = 0; part < 2; ++part) {
    const double* u = raw.data() + part * K;
    double mx = u[0];
    for (int j = 1; j < K; ++j) mx = std::fmax(mx, u[j]);
    double sum = 0.0;
    std::vector<double> e(K);
    for (int j = 0; j < K; ++j) { e[j] = std::exp(u[j] - mx); sum += e[j]; }
    std::vector<double>& pos = part == 0 ? kn.x : kn.y;
    pos[0] = m.range_min;
    double acc = m.range_min;
    for (int j = 0; j < K - 1; ++j) {
      acc += (e[j] / sum) * budget + m.min_bin_size;
      pos[j + 1] = acc;
    }
    pos[K] = m.range_max;
  }
  const double offset = std::log(std::exp(1.0 - m.min_knot_slope) - 1.0);
  kn.d[0] = 1.0; kn.d[K] = 1.0;
  for (int j = 1; j < K; ++j)
    kn.d[j] = softplus(raw[2 * K + j] + offset) + m.min_knot_slope;
}

static inline void rqs_inverse_one(const Knots& kn, double u,
                                   double& z_out, double& logdet_out) {
  const int K = (int)kn.x.size() - 1;
  if (u <= kn.y[0]) { z_out = u - kn.y[0] + kn.x[0]; logdet_out = 0.0; return; }
  if (u >= kn.y[K]) { z_out = u - kn.y[K] + kn.x[K]; logdet_out = 0.0; return; }
  int lo = 0, hi = K;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (kn.y[mid] <= u) lo = mid; else hi = mid;
  }
  const int i = lo;
  const double bin_width  = kn.x[i + 1] - kn.x[i];
  const double bin_height = kn.y[i + 1] - kn.y[i];
  const double bin_slope  = bin_height / bin_width;
  const double d_lo = kn.d[i], d_hi = kn.d[i + 1];
  double w = (u - kn.y[i]) / bin_height;
  w = std::fmin(std::fmax(w, 0.0), 1.0);
  const double slopes_term = d_hi + d_lo - 2.0 * bin_slope;
  const double c = -bin_slope * w;
  const double b = d_lo - slopes_term * w;
  const double a = bin_slope - b;
  const double sqrt_diff = b * b - 4.0 * a * c;
  double safe_sqrt = std::sqrt(std::fmax(sqrt_diff, DBL_MIN));
  if (sqrt_diff <= 0.0) safe_sqrt = 0.0;
  const double num = (b >= 0.0) ? 2.0 * c : (-b + safe_sqrt);
  const double den = (b >= 0.0) ? (-b - safe_sqrt) : 2.0 * a;
  double zeta = num / den;
  zeta = std::fmin(std::fmax(zeta, 0.0), 1.0);
  z_out = kn.x[i] + bin_width * zeta;
  const double sq_z = zeta * zeta;
  const double z1mz = zeta - sq_z;
  const double sq_1mz = (1.0 - zeta) * (1.0 - zeta);
  const double denominator = bin_slope + slopes_term * z1mz;
  logdet_out = -2.0 * std::log(bin_slope)
      - std::log(d_hi * sq_z + 2.0 * bin_slope * z1mz + d_lo * sq_1mz)
      + 2.0 * std::log(denominator);
}

static inline bool in_box(const DdmModel& m, const double* theta) {
  for (int j = 0; j < m.n_ctx; ++j)
    if (std::isnan(theta[j]) || theta[j] < m.lower[j] || theta[j] > m.upper[j])
      return false;
  return true;
}

// run both networks once for a parameter vector
static void ddm_condition(const DdmModel& m, const double* theta, DdmCond& c) {
  std::vector<double> ctx(m.n_ctx + 1), ctx_clf(m.n_ctx), raw;
  for (int j = 0; j < m.n_ctx; ++j)
    ctx_clf[j] = ctx[j] = (theta[j] - m.scaler_mean[j]) / m.scaler_scale[j];
  for (int r = 0; r < 2; ++r) {
    ctx[m.n_ctx] = (double)(r + 1);      // raw response code
    mlp_forward(m.flow, ctx, raw);
    build_knots(m, raw, c.kn[r]);
  }
  std::vector<double> logit_v;
  mlp_forward(m.clf, ctx_clf, logit_v);
  const double logit = logit_v[0];
  c.log_p_R[0] = -softplus(logit);       // log P(R = 1)
  c.log_p_R[1] = -softplus(-logit);      // log P(R = 2)
}

static void eval_one(const DdmCond& c, double rt, int R, int o,
                     NumericVector& pdf, NumericVector& cdf,
                     NumericVector& log_pdf, NumericVector& p_R) {
  const int r = R - 1;
  const double u = std::log(rt);
  double z, logdet;
  rqs_inverse_one(c.kn[r], u, z, logdet);
  const double lp = R::dnorm(z, 0.0, 1.0, 1) + logdet - u + c.log_p_R[r];
  log_pdf[o] = lp;
  pdf[o] = std::exp(lp);
  cdf[o] = R::pnorm(z, 0.0, 1.0, 1, 0) * std::exp(c.log_p_R[r]);
  p_R[o] = std::exp(c.log_p_R[r]);
}

// ---------------------------------------------------------------------------
// R interface
// ---------------------------------------------------------------------------

static void load_mlp(List mlp_list, Mlp& m) {
  if (as<std::string>(mlp_list["activation"]) != "gelu_tanh")
    stop("ddm_build: unsupported activation.");
  List layers = mlp_list["layers"];
  for (int i = 0; i < layers.size(); ++i) {
    List l = layers[i];
    NumericMatrix W = l["W"];
    Mlp::Layer lay;
    lay.nin = W.nrow(); lay.nout = W.ncol();
    lay.W.assign(W.begin(), W.end());
    lay.b = as<std::vector<double>>(l["b"]);
    m.layers.push_back(std::move(lay));
  }
  m.use_norm = mlp_list.containsElementNamed("use_norm") &&
               as<bool>(mlp_list["use_norm"]);
  if (m.use_norm) {
    List norms = mlp_list["norms"];
    for (int i = 0; i < norms.size(); ++i) {
      List n = norms[i];
      Mlp::Norm nrm;
      nrm.scale = as<std::vector<double>>(n["scale"]);
      nrm.bias = as<std::vector<double>>(n["bias"]);
      nrm.eps = as<double>(n["eps"]);
      m.norms.push_back(std::move(nrm));
    }
  }
}

static void load_ddm_from_list(List fl, DdmModel& m) {
  List sp = fl["spline"];
  if (as<int>(sp["num_splines"]) != 1 ||
      as<std::string>(sp["output_transform"]) != "exp" ||
      as<std::string>(sp["base_distribution"]) != "standard_normal")
    stop("ddm_build: only single-spline / exp-transform / normal-base flows "
         "are supported.");
  m.num_bins = as<int>(sp["num_bins"]);
  m.range_min = as<double>(sp["range_min"]);
  m.range_max = as<double>(sp["range_max"]);
  m.min_bin_size = as<double>(sp["min_bin_size"]);
  m.min_knot_slope = as<double>(sp["min_knot_slope"]);
  List scaler = fl["scaler"];
  m.scaler_mean = as<std::vector<double>>(scaler["mean"]);
  m.scaler_scale = as<std::vector<double>>(scaler["scale"]);
  m.n_ctx = (int)m.scaler_mean.size();
  List bounds = fl["bounds_sampled"];
  m.lower = as<std::vector<double>>(bounds["lower"]);
  m.upper = as<std::vector<double>>(bounds["upper"]);
  load_mlp(fl["flow_mlp"], m.flow);
  load_mlp(fl["classifier_mlp"], m.clf);
  if (m.flow.layers[0].nin != m.n_ctx + 1)
    stop("ddm_build: flow input dim must be n_params + 1 (response).");
  if (m.clf.layers[0].nin != m.n_ctx)
    stop("ddm_build: classifier input dim must be n_params.");
}

// [[Rcpp::export]]
SEXP ddm_build(List fl) {
  DdmModel* m = new DdmModel();
  load_ddm_from_list(fl, *m);
  XPtr<DdmModel> ptr(m, true);
  ptr.attr("class") = "ddm_flow_model";
  ptr.attr("model") = as<std::string>(fl["model"]);
  return ptr;
}

struct DdmEnsemble { std::vector<DdmModel> members; };

// Equal-weight ensemble of DDM flows (same context/bounds, e.g. training
// seeds). fls: an R list of weight lists as accepted by ddm_build.
// [[Rcpp::export]]
SEXP ddm_build_ensemble(List fls) {
  if (fls.size() < 1) stop("ddm_build_ensemble: need at least one member.");
  DdmEnsemble* e = new DdmEnsemble();
  e->members.resize(fls.size());
  for (int k = 0; k < fls.size(); ++k) {
    load_ddm_from_list(fls[k], e->members[k]);
    if (e->members[k].n_ctx != e->members[0].n_ctx ||
        e->members[k].lower != e->members[0].lower ||
        e->members[k].upper != e->members[0].upper)
      stop("ddm_build_ensemble: members disagree on context dim or bounds.");
  }
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
  NumericVector pdf(n), cdf(n), log_pdf(n), p_R(n);
  if (!in_box(m, theta.begin())) {
    std::fill(log_pdf.begin(), log_pdf.end(), R_NegInf);
    std::fill(p_R.begin(), p_R.end(), NA_REAL);
    return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                        _["log_pdf"] = log_pdf, _["p_R"] = p_R,
                        _["in_box"] = false);
  }
  DdmCond c;
  ddm_condition(m, theta.begin(), c);
  for (int t = 0; t < n; ++t) {
    if (R[t] != 1 && R[t] != 2) stop("R must be 1 or 2.");
    eval_one(c, rt[t], R[t], t, pdf, cdf, log_pdf, p_R);
  }
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                      _["log_pdf"] = log_pdf, _["p_R"] = p_R,
                      _["in_box"] = true);
}

// Trial-wise parameter rows (n x K_ctx) with rt (n) and R (n). Consecutive
// duplicate rows reuse the conditioned state (both knot sets + classifier).
// [[Rcpp::export]]
List ddm_eval_trials_cpp(SEXP ptr_, NumericMatrix theta, NumericVector rt,
                         IntegerVector R) {
  XPtr<DdmModel> ptr(ptr_);
  const DdmModel& m = *ptr;
  const int n = rt.size();
  if (theta.nrow() != n || theta.ncol() != m.n_ctx)
    stop("theta must be length(rt) x %d.", m.n_ctx);
  if (R.size() != n) stop("rt and R must have equal length.");
  NumericVector pdf(n), cdf(n), log_pdf(n), p_R(n);
  std::vector<double> row(m.n_ctx), prev(m.n_ctx,
                                         std::numeric_limits<double>::quiet_NaN());
  DdmCond c;
  bool have = false, prev_in_box = false;
  for (int t = 0; t < n; ++t) {
    if (R[t] != 1 && R[t] != 2) stop("R must be 1 or 2.");
    for (int j = 0; j < m.n_ctx; ++j) row[j] = theta(t, j);
    const bool same = have && std::equal(row.begin(), row.end(), prev.begin());
    if (!same) {
      prev = row;
      prev_in_box = in_box(m, row.data());
      if (prev_in_box) ddm_condition(m, row.data(), c);
      have = true;
    }
    if (!prev_in_box) {
      pdf[t] = 0.0; cdf[t] = 0.0; log_pdf[t] = R_NegInf; p_R[t] = NA_REAL;
    } else {
      eval_one(c, rt[t], R[t], t, pdf, cdf, log_pdf, p_R);
    }
  }
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                      _["log_pdf"] = log_pdf, _["p_R"] = p_R);
}

// Trial-wise ensemble evaluation: equal-weight mixture over members.
// pdf/cdf/p_R are means over members; log_pdf = log(mean pdf). Consecutive
// duplicate parameter rows reuse all k conditioned states.
// [[Rcpp::export]]
List ddm_ens_eval_trials_cpp(SEXP ptr_, NumericMatrix theta, NumericVector rt,
                             IntegerVector R) {
  XPtr<DdmEnsemble> ptr(ptr_);
  const DdmEnsemble& e = *ptr;
  const int K = (int)e.members.size();
  const DdmModel& m0 = e.members[0];
  const int n = rt.size();
  if (theta.nrow() != n || theta.ncol() != m0.n_ctx)
    stop("theta must be length(rt) x %d.", m0.n_ctx);
  if (R.size() != n) stop("rt and R must have equal length.");

  NumericVector pdf(n), cdf(n), log_pdf(n), p_R(n);
  NumericVector pdf1(1), cdf1(1), lp1(1), pR1(1);   // per-member scratch
  std::vector<double> row(m0.n_ctx), prev(m0.n_ctx,
                                          std::numeric_limits<double>::quiet_NaN());
  std::vector<DdmCond> cond(K);
  bool have = false, prev_in_box = false;

  for (int t = 0; t < n; ++t) {
    if (R[t] != 1 && R[t] != 2) stop("R must be 1 or 2.");
    for (int j = 0; j < m0.n_ctx; ++j) row[j] = theta(t, j);
    const bool same = have && std::equal(row.begin(), row.end(), prev.begin());
    if (!same) {
      prev = row;
      prev_in_box = in_box(m0, row.data());   // members share bounds
      if (prev_in_box)
        for (int k = 0; k < K; ++k)
          ddm_condition(e.members[k], row.data(), cond[k]);
      have = true;
    }
    if (!prev_in_box) {
      pdf[t] = 0.0; cdf[t] = 0.0; log_pdf[t] = R_NegInf; p_R[t] = NA_REAL;
    } else {
      double s_pdf = 0.0, s_cdf = 0.0, s_pR = 0.0;
      for (int k = 0; k < K; ++k) {
        eval_one(cond[k], rt[t], R[t], 0, pdf1, cdf1, lp1, pR1);
        s_pdf += pdf1[0]; s_cdf += cdf1[0]; s_pR += pR1[0];
      }
      pdf[t] = s_pdf / K;
      cdf[t] = s_cdf / K;
      p_R[t] = s_pR / K;
      log_pdf[t] = std::log(pdf[t]);
    }
  }
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                      _["log_pdf"] = log_pdf, _["p_R"] = p_R);
}
