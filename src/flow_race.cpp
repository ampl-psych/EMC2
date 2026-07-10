// C++ (Rcpp) evaluator for exported race-model flows (LNR, RDM, ...),
// used by the LNRnn / RDMnn model types (R/model_NN.R).
//
// KEEP IN SYNC with the development copy in the likelihood-kde-normflow
// repo (port/src/flow_race.cpp), which is validated against float64 golden
// vectors from the Python training pipeline. Direct translation of the
// validated R reference (port/R/flow_race.R), which in turn mirrors distrax
// rational_quadratic_spline.py and flow/src/flow/config.py. Model-agnostic:
// any context dimension, any MLP depth/widths, optional LayerNorm, any bin
// count. Only single-spline chains with an exp output transform and a
// standard-normal base are supported (what flow.config.spline_flow builds).
//
// Design for EMC2 race models:
//   flow_build(fl)                 -> external pointer, built once per session
//   flow_eval_cpp(ptr, theta, rt)  -> one parameter vector, many rts
//                                     (conditioner runs once: the MCMC case)
//   flow_eval_trials_cpp(ptr, Theta, rt)
//                                  -> trial-wise parameter rows; consecutive
//                                     duplicate rows reuse the spline knots,
//                                     so cell-structured designs stay cheap
//
// Out-of-box parameters (training-region bounding box, sampled scale):
// pdf = 0, cdf = 0 immediately (log_pdf = -Inf, log_sf = 0). Rejection
// sentinel, not extrapolation.
//
// Uses Rmath dnorm/pnorm so results are bit-comparable with the R reference.

#include <Rcpp.h>
#include <cmath>
#include <cfloat>
#include <vector>
using namespace Rcpp;

struct FlowModel {
  int n_ctx = 0;
  int num_bins = 0;
  double range_min = 0, range_max = 0, min_bin_size = 0, min_knot_slope = 0;
  std::vector<double> scaler_mean, scaler_scale;
  std::vector<double> lower, upper;                    // bounding box (sampled scale)
  // MLP: layers[i] is column-major W (rows_in x cols_out) + bias
  struct Layer { int nin, nout; std::vector<double> W, b; };
  std::vector<Layer> layers;
  bool use_norm = false;
  struct Norm { std::vector<double> scale, bias; double eps; };
  std::vector<Norm> norms;
};

struct Knots {                 // amortized per-parameter-vector state
  std::vector<double> x, y, d; // each num_bins + 1
};

static inline double gelu_tanh(double x) {
  const double c = 0.7978845608028654;  // sqrt(2/pi)
  return 0.5 * x * (1.0 + std::tanh(c * (x + 0.044715 * x * x * x)));
}

static inline double softplus(double x) {
  return std::fmax(x, 0.0) + std::log1p(std::exp(-std::fabs(x)));
}

static void mlp_forward(const FlowModel& m, const double* theta,
                        std::vector<double>& raw) {
  std::vector<double> h(m.n_ctx), h2;
  for (int j = 0; j < m.n_ctx; ++j)
    h[j] = (theta[j] - m.scaler_mean[j]) / m.scaler_scale[j];

  const int L = (int)m.layers.size();
  for (int li = 0; li < L; ++li) {
    const FlowModel::Layer& lay = m.layers[li];
    h2.assign(lay.nout, 0.0);
    for (int i = 0; i < lay.nin; ++i) {
      const double hi = h[i];
      const double* Wcol = &lay.W[(size_t)i];  // row i of W, strided by nin
      for (int j = 0; j < lay.nout; ++j)
        h2[j] += hi * lay.W[(size_t)j * lay.nin + i];
      (void)Wcol;
    }
    for (int j = 0; j < lay.nout; ++j) h2[j] += lay.b[j];

    if (li < L - 1) {
      if (m.use_norm) {
        const FlowModel::Norm& nrm = m.norms[li];
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
  raw.swap(h);
}

// raw (3K+1) -> knot positions and slopes; mirrors distrax + "identity" bounds
static void build_knots(const FlowModel& m, const std::vector<double>& raw,
                        Knots& kn) {
  const int K = m.num_bins;
  const double range_size = m.range_max - m.range_min;
  const double budget = range_size - K * m.min_bin_size;

  kn.x.resize(K + 1); kn.y.resize(K + 1); kn.d.resize(K + 1);

  for (int part = 0; part < 2; ++part) {           // 0: widths/x, 1: heights/y
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
  kn.d[0] = 1.0; kn.d[K] = 1.0;                    // boundary_slopes = "identity"
  for (int j = 1; j < K; ++j)
    kn.d[j] = softplus(raw[2 * K + j] + offset) + m.min_knot_slope;
}

// Inverse spline at u = log(rt); writes z and log|dz/du|.
// Mirrors distrax _rational_quadratic_spline_inv (stable quadratic root,
// identity tails).
static inline void rqs_inverse_one(const Knots& kn, double u,
                                   double& z_out, double& logdet_out) {
  const int K = (int)kn.x.size() - 1;
  if (u <= kn.y[0]) { z_out = u - kn.y[0] + kn.x[0]; logdet_out = 0.0; return; }
  if (u >= kn.y[K]) { z_out = u - kn.y[K] + kn.x[K]; logdet_out = 0.0; return; }

  // binary search: largest idx with y[idx] <= u
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

static inline bool in_box(const FlowModel& m, const double* theta) {
  for (int j = 0; j < m.n_ctx; ++j)
    if (std::isnan(theta[j]) || theta[j] < m.lower[j] || theta[j] > m.upper[j])
      return false;
  return true;
}

// Evaluate all outputs for one knot set and a block of rts.
static void eval_block(const Knots& kn, const double* rt, int n, int offset,
                       NumericVector& pdf, NumericVector& cdf,
                       NumericVector& log_pdf, NumericVector& log_sf) {
  for (int t = 0; t < n; ++t) {
    const double u = std::log(rt[t]);
    double z, logdet;
    rqs_inverse_one(kn, u, z, logdet);
    const double lp = R::dnorm(z, 0.0, 1.0, 1) + logdet - u;
    const int o = offset + t;
    log_pdf[o] = lp;
    pdf[o] = std::exp(lp);
    cdf[o] = R::pnorm(z, 0.0, 1.0, 1, 0);
    log_sf[o] = R::pnorm(z, 0.0, 1.0, 0, 1);
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
  List sp = fl["spline"];
  if (as<int>(sp["num_splines"]) != 1 ||
      as<std::string>(sp["output_transform"]) != "exp" ||
      as<std::string>(sp["base_distribution"]) != "standard_normal")
    stop("flow_build: only single-spline / exp-transform / normal-base flows "
         "are supported (what flow.config.spline_flow builds).");

  FlowModel* m = new FlowModel();
  m->num_bins = as<int>(sp["num_bins"]);
  m->range_min = as<double>(sp["range_min"]);
  m->range_max = as<double>(sp["range_max"]);
  m->min_bin_size = as<double>(sp["min_bin_size"]);
  m->min_knot_slope = as<double>(sp["min_knot_slope"]);

  List scaler = fl["scaler"];
  m->scaler_mean = as<std::vector<double>>(scaler["mean"]);
  m->scaler_scale = as<std::vector<double>>(scaler["scale"]);
  m->n_ctx = (int)m->scaler_mean.size();

  List bounds = fl["bounds_sampled"];
  m->lower = as<std::vector<double>>(bounds["lower"]);
  m->upper = as<std::vector<double>>(bounds["upper"]);

  List mlp = fl["mlp"];
  if (as<std::string>(mlp["activation"]) != "gelu_tanh")
    stop("flow_build: unsupported activation.");
  List layers = mlp["layers"];
  for (int i = 0; i < layers.size(); ++i) {
    List l = layers[i];
    NumericMatrix W = l["W"];
    FlowModel::Layer lay;
    lay.nin = W.nrow(); lay.nout = W.ncol();
    lay.W.assign(W.begin(), W.end());              // column-major
    lay.b = as<std::vector<double>>(l["b"]);
    m->layers.push_back(std::move(lay));
  }
  m->use_norm = as<bool>(mlp["use_norm"]);
  if (m->use_norm) {
    List norms = mlp["norms"];
    for (int i = 0; i < norms.size(); ++i) {
      List n = norms[i];
      FlowModel::Norm nrm;
      nrm.scale = as<std::vector<double>>(n["scale"]);
      nrm.bias = as<std::vector<double>>(n["bias"]);
      nrm.eps = as<double>(n["eps"]);
      m->norms.push_back(std::move(nrm));
    }
  }

  XPtr<FlowModel> ptr(m, true);
  ptr.attr("class") = "flow_model";
  ptr.attr("model") = as<std::string>(fl["model"]);
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

  if (!in_box(m, theta.begin())) {
    std::fill(log_pdf.begin(), log_pdf.end(), R_NegInf);
    return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                        _["log_pdf"] = log_pdf, _["log_sf"] = log_sf,
                        _["in_box"] = false);
  }
  std::vector<double> raw;
  mlp_forward(m, theta.begin(), raw);
  Knots kn;
  build_knots(m, raw, kn);
  eval_block(kn, rt.begin(), n, 0, pdf, cdf, log_pdf, log_sf);
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                      _["log_pdf"] = log_pdf, _["log_sf"] = log_sf,
                      _["in_box"] = true);
}

// Trial-wise parameter rows (n x K_ctx) paired with rt (n) — the shape of
// EMC2 dfun/pfun inputs. Consecutive duplicate rows reuse the spline knots,
// so designs with cell-constant parameters pay the MLP once per cell.
// [[Rcpp::export]]
List flow_eval_trials_cpp(SEXP ptr_, NumericMatrix theta, NumericVector rt) {
  XPtr<FlowModel> ptr(ptr_);
  const FlowModel& m = *ptr;
  const int n = rt.size();
  if (theta.nrow() != n || theta.ncol() != m.n_ctx)
    stop("theta must be length(rt) x %d.", m.n_ctx);

  NumericVector pdf(n), cdf(n), log_pdf(n), log_sf(n);
  std::vector<double> row(m.n_ctx), prev(m.n_ctx,
                                         std::numeric_limits<double>::quiet_NaN());
  std::vector<double> raw;
  Knots kn;
  bool have_kn = false, prev_in_box = false;

  for (int t = 0; t < n; ++t) {
    for (int j = 0; j < m.n_ctx; ++j) row[j] = theta(t, j);
    const bool same = have_kn && std::equal(row.begin(), row.end(), prev.begin());
    if (!same) {
      prev = row;
      prev_in_box = in_box(m, row.data());
      if (prev_in_box) {
        mlp_forward(m, row.data(), raw);
        build_knots(m, raw, kn);
      }
      have_kn = true;
    }
    if (!prev_in_box) {
      pdf[t] = 0.0; cdf[t] = 0.0; log_pdf[t] = R_NegInf; log_sf[t] = 0.0;
    } else {
      eval_block(kn, &rt[t], 1, t, pdf, cdf, log_pdf, log_sf);
    }
  }
  return List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                      _["log_pdf"] = log_pdf, _["log_sf"] = log_sf);
}
