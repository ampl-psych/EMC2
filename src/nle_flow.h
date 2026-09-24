// Shared numerical core of the neural-likelihood evaluators (flow_ddm.cpp,
// flow_race.cpp; the plain-MLP evaluators planned for LANs reuse the MLP part).
//
// Origin: the per-row helpers of the NLE handover port
// (~/Downloads/julian/nle/handover/nle-methods/port/src/flow_{ddm,race}.cpp),
// which kept two in-sync copies. The spline arithmetic (build_knots,
// rqs_inverse_one) is unchanged formula for formula (build_knots takes its exps
// from the vectorised vec_exp); what is new is
//   * mlp_forward_batch: the conditioner runs on a block of U inputs at once,
//     one Armadillo GEMM per layer (GEMV when U = 1), instead of a hand loop
//     per parameter vector; GELU is evaluated in the equivalent logistic form
//     with the package's vectorised exp (math_utils.h) instead of std::tanh;
//   * RowIndex: distinct parameter rows are found with a hash index (plus a
//     consecutive-duplicate fast path), replacing the old consecutive-only
//     cache, so interleaved designs (race winners/losers, per-trial covariates
//     with repeats) also condition each distinct row once;
//   * KnotBlock: knot sets for a block of conditioned rows in flat storage.
// Log densities agree with the R port to ~1e-14 typically and ~5e-11 at worst
// (ill-conditioned spline points, where the handover's per-row C++ shows the
// same spread); tests/testthat/test-nn-port.R holds them to 1e-10.
//
// Everything below the R-facing loaders is reentrant (no R API calls, no
// shared mutable state), so an evaluator core can run inside OpenMP threads.

#ifndef EMC2_NLE_FLOW_H
#define EMC2_NLE_FLOW_H

#include <RcppArmadillo.h>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include "math_utils.h"

namespace nle {

// Columns conditioned per GEMM block: bounds the activation memory
// (widest layer x NLE_BLOCK doubles) while keeping the products large.
constexpr int NLE_BLOCK = 256;

inline double softplus(double x) {
  return std::fmax(x, 0.0) + std::log1p(std::exp(-std::fabs(x)));
}

// ---------------------------------------------------------------------------
// MLP. The artefact stores W as n_in x n_out; Wt[l] keeps its transpose
// (n_out x n_in) so that each layer is a plain product Wt * H: reference BLAS
// runs that form as column updates (vectorisable), ~2.4x faster than the
// dot-product form W.t() * H, and the per-element summation order is the same.
// Hidden layers: affine -> optional LayerNorm -> GELU(tanh); last layer affine.
// ---------------------------------------------------------------------------
struct Mlp {
  std::vector<arma::mat> Wt;
  std::vector<arma::vec> b;
  bool use_norm = false;
  std::vector<arma::vec> norm_scale, norm_bias;
  std::vector<double> norm_eps;
  int n_in()  const { return (int)Wt.front().n_cols; }
  int n_out() const { return (int)Wt.back().n_rows; }
};

inline void load_mlp(Rcpp::List mlp_list, Mlp& m, const char* who) {
  if (Rcpp::as<std::string>(mlp_list["activation"]) != "gelu_tanh")
    Rcpp::stop("%s: unsupported activation.", who);
  Rcpp::List layers = mlp_list["layers"];
  if (layers.size() < 1) Rcpp::stop("%s: MLP has no layers.", who);
  for (int i = 0; i < layers.size(); ++i) {
    Rcpp::List l = layers[i];
    Rcpp::NumericMatrix W = l["W"];
    Rcpp::NumericVector b = l["b"];
    if (b.size() != W.ncol()) Rcpp::stop("%s: layer %d bias length != ncol(W).", who, i + 1);
    if (i > 0 && W.nrow() != (int)m.Wt.back().n_rows)
      Rcpp::stop("%s: layer %d input dim does not match layer %d output dim.", who, i + 1, i);
    m.Wt.push_back(arma::mat(W.begin(), W.nrow(), W.ncol()).t());
    m.b.emplace_back(b.begin(), b.size());
  }
  m.use_norm = mlp_list.containsElementNamed("use_norm") &&
               !Rf_isNull(mlp_list["use_norm"]) && Rcpp::as<bool>(mlp_list["use_norm"]);
  if (m.use_norm) {
    Rcpp::List norms = mlp_list["norms"];
    if (norms.size() < (int)m.Wt.size() - 1)
      Rcpp::stop("%s: use_norm is set but fewer norms than hidden layers.", who);
    for (int i = 0; i < (int)m.Wt.size() - 1; ++i) {
      Rcpp::List n = norms[i];
      Rcpp::NumericVector sc = n["scale"], bi = n["bias"];
      m.norm_scale.emplace_back(sc.begin(), sc.size());
      m.norm_bias.emplace_back(bi.begin(), bi.size());
      m.norm_eps.push_back(Rcpp::as<double>(n["eps"]));
    }
  }
}

// GELU, tanh approximation (jax.nn.gelu(approximate=True)), over a block:
// 0.5 x (1 + tanh(y)) = x / (1 + exp(-2y)), y = sqrt(2/pi) (x + 0.044715 x^3).
// The logistic form is exact algebra, avoids the cancellation in 1 + tanh(y)
// for y << 0, and needs one exp per element, done by vec_exp (Accelerate's
// vvexp on macOS). x -> -Inf gives x / Inf = -0, as the tanh form does.
inline void gelu_block(double* z, arma::uword N, std::vector<double>& e) {
  const double c2 = -2.0 * 0.7978845608028654;       // -2 sqrt(2/pi)
  e.resize(N);
  for (arma::uword k = 0; k < N; ++k) e[k] = c2 * (z[k] + 0.044715 * z[k] * z[k] * z[k]);
  vec_exp(e.data(), (int)N);
  for (arma::uword k = 0; k < N; ++k) z[k] = z[k] / (1.0 + e[k]);
}

// X: n_in x U, one input vector per column -> n_out x U. One matrix product per
// layer (GEMV when U = 1); columns never interact.
inline arma::mat mlp_forward_batch(const Mlp& m, const arma::mat& X) {
  arma::mat H;
  std::vector<double> e;
  const int L = (int)m.Wt.size();
  for (int l = 0; l < L; ++l) {
    arma::mat Z = m.Wt[l] * (l == 0 ? X : H);
    Z.each_col() += m.b[l];
    if (l < L - 1) {
      const arma::uword nout = Z.n_rows, U = Z.n_cols;
      if (m.use_norm) {
        const double* sc = m.norm_scale[l].memptr();
        const double* bi = m.norm_bias[l].memptr();
        for (arma::uword c = 0; c < U; ++c) {
          double* z = Z.colptr(c);
          double mu = 0.0;
          for (arma::uword j = 0; j < nout; ++j) mu += z[j];
          mu /= nout;
          double var = 0.0;
          for (arma::uword j = 0; j < nout; ++j) var += (z[j] - mu) * (z[j] - mu);
          var /= nout;
          const double inv_sd = 1.0 / std::sqrt(var + m.norm_eps[l]);
          for (arma::uword j = 0; j < nout; ++j)
            z[j] = (z[j] - mu) * inv_sd * sc[j] + bi[j];
        }
      }
      gelu_block(Z.memptr(), Z.n_elem, e);
    }
    H.steal_mem(Z);
  }
  return H;
}

// ---------------------------------------------------------------------------
// Rational-quadratic spline (distrax conventions, identity tails/boundary
// slopes), exp output transform, standard-normal base.
// ---------------------------------------------------------------------------
struct Spline {
  int num_bins = 0;
  double range_min = 0, range_max = 0, min_bin_size = 0, min_knot_slope = 0;
};

inline void load_spline(Rcpp::List sp, Spline& s, const char* who) {
  if (Rcpp::as<int>(sp["num_splines"]) != 1 ||
      Rcpp::as<std::string>(sp["output_transform"]) != "exp" ||
      Rcpp::as<std::string>(sp["base_distribution"]) != "standard_normal")
    Rcpp::stop("%s: only single-spline / exp-transform / normal-base flows "
               "are supported (what flow.config.spline_flow builds).", who);
  s.num_bins = Rcpp::as<int>(sp["num_bins"]);
  s.range_min = Rcpp::as<double>(sp["range_min"]);
  s.range_max = Rcpp::as<double>(sp["range_max"]);
  s.min_bin_size = Rcpp::as<double>(sp["min_bin_size"]);
  s.min_knot_slope = Rcpp::as<double>(sp["min_knot_slope"]);
}

// Knot sets for a block of columns: x, y, d each (K + 1) per column.
struct KnotBlock {
  int K1 = 0;
  std::vector<double> x, y, d, e;   // e: softmax scratch (K)
  void reset(int num_bins, int ncol) {
    K1 = num_bins + 1;
    const size_t n = (size_t)K1 * ncol;
    if (x.size() < n) { x.resize(n); y.resize(n); d.resize(n); }
    e.resize(num_bins);
  }
  const double* xc(int c) const { return x.data() + (size_t)c * K1; }
  const double* yc(int c) const { return y.data() + (size_t)c * K1; }
  const double* dc(int c) const { return d.data() + (size_t)c * K1; }
};

// raw (3K + 1 conditioner outputs) -> knot positions and slopes of column c.
inline void build_knots(const Spline& s, const double* raw, KnotBlock& kb, int c) {
  const int K = s.num_bins;
  const double range_size = s.range_max - s.range_min;
  const double budget = range_size - K * s.min_bin_size;
  double* e = kb.e.data();
  for (int part = 0; part < 2; ++part) {           // 0: widths/x, 1: heights/y
    const double* u = raw + part * K;
    double mx = u[0];
    for (int j = 1; j < K; ++j) mx = std::fmax(mx, u[j]);
    for (int j = 0; j < K; ++j) e[j] = u[j] - mx;
    vec_exp(e, K);
    double sum = 0.0;
    for (int j = 0; j < K; ++j) sum += e[j];
    double* pos = (part == 0 ? kb.x.data() : kb.y.data()) + (size_t)c * kb.K1;
    pos[0] = s.range_min;
    double acc = s.range_min;
    for (int j = 0; j < K - 1; ++j) {
      acc += (e[j] / sum) * budget + s.min_bin_size;
      pos[j + 1] = acc;
    }
    pos[K] = s.range_max;
  }
  const double offset = std::log(std::exp(1.0 - s.min_knot_slope) - 1.0);
  double* d = kb.d.data() + (size_t)c * kb.K1;
  d[0] = 1.0; d[K] = 1.0;                          // boundary_slopes = "identity"
  // softplus(v) = max(v, 0) + log1p(exp(-|v|)), exps vectorised
  for (int j = 1; j < K; ++j) e[j - 1] = -std::fabs(raw[2 * K + j] + offset);
  vec_exp(e, K - 1);
  for (int j = 1; j < K; ++j)
    d[j] = (std::fmax(raw[2 * K + j] + offset, 0.0) + std::log1p(e[j - 1])) + s.min_knot_slope;
}

// Inverse spline at u = log(rt); writes z and log|dz/du|. Mirrors distrax
// _rational_quadratic_spline_inv (stable quadratic root, identity tails).
inline void rqs_inverse_one(const double* x, const double* y, const double* d,
                            int K, double u, double& z_out, double& logdet_out) {
  if (u <= y[0]) { z_out = u - y[0] + x[0]; logdet_out = 0.0; return; }
  if (u >= y[K]) { z_out = u - y[K] + x[K]; logdet_out = 0.0; return; }
  int lo = 0, hi = K;                  // largest idx with y[idx] <= u
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (y[mid] <= u) lo = mid; else hi = mid;
  }
  const int i = lo;
  const double bin_width  = x[i + 1] - x[i];
  const double bin_height = y[i + 1] - y[i];
  const double bin_slope  = bin_height / bin_width;
  const double d_lo = d[i], d_hi = d[i + 1];
  double w = (u - y[i]) / bin_height;
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
  z_out = x[i] + bin_width * zeta;
  const double sq_z = zeta * zeta;
  const double z1mz = zeta - sq_z;
  const double sq_1mz = (1.0 - zeta) * (1.0 - zeta);
  const double denominator = bin_slope + slopes_term * z1mz;
  logdet_out = -2.0 * std::log(bin_slope)
      - std::log(d_hi * sq_z + 2.0 * bin_slope * z1mz + d_lo * sq_1mz)
      + 2.0 * std::log(denominator);
}

// Training-region box on the sampled scale; NaN is out of box.
inline bool in_box(const std::vector<double>& lower, const std::vector<double>& upper,
                   const double* theta) {
  for (size_t j = 0; j < lower.size(); ++j)
    if (std::isnan(theta[j]) || theta[j] < lower[j] || theta[j] > upper[j]) return false;
  return true;
}

// ---------------------------------------------------------------------------
// Distinct rows of a column-major n x k matrix, in first-appearance order.
// Rows are equal when every entry is equal as a double (0 == -0; NaN rows are
// out of box whatever they are matched to).
// ---------------------------------------------------------------------------
struct RowIndex {
  std::vector<int> uid;     // length n: row -> distinct-row id
  std::vector<int> first;   // length U: first row carrying each id
  int U() const { return (int)first.size(); }
};

inline uint64_t canon_bits(double v) {
  if (v == 0.0) v = 0.0;                            // fold -0 onto +0
  uint64_t b;
  std::memcpy(&b, &v, sizeof b);
  return b;
}

inline RowIndex index_rows(const double* M, int n, int k) {
  RowIndex ix;
  ix.uid.resize(n);
  if (n == 0) return ix;
  auto row_eq = [&](int a, int b) {
    for (int j = 0; j < k; ++j)
      if (canon_bits(M[a + (size_t)j * n]) != canon_bits(M[b + (size_t)j * n])) return false;
    return true;
  };
  auto row_hash = [&](int a) {
    uint64_t h = 0x9E3779B97F4A7C15ULL;
    for (int j = 0; j < k; ++j) {                   // splitmix64-style mixing
      uint64_t z = canon_bits(M[a + (size_t)j * n]) + h + 0x9E3779B97F4A7C15ULL;
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
      h = z ^ (z >> 31);
    }
    return h;
  };
  size_t cap = 16;
  while (cap < 2 * (size_t)n) cap <<= 1;
  std::vector<int> table;                           // built lazily: many inputs are one run
  for (int t = 0; t < n; ++t) {
    if (t > 0 && row_eq(t, t - 1)) { ix.uid[t] = ix.uid[t - 1]; continue; }
    if (table.empty()) {
      table.assign(cap, -1);
      for (int u = 0; u < ix.U(); ++u) {            // (only ever the first row)
        size_t h = row_hash(ix.first[u]) & (cap - 1);
        while (table[h] >= 0) h = (h + 1) & (cap - 1);
        table[h] = u;
      }
    }
    size_t h = row_hash(t) & (cap - 1);
    for (;;) {
      const int u = table[h];
      if (u < 0) {
        ix.uid[t] = ix.U(); table[h] = ix.U(); ix.first.push_back(t);
        break;
      }
      if (row_eq(t, ix.first[u])) { ix.uid[t] = u; break; }
      h = (h + 1) & (cap - 1);
    }
  }
  return ix;
}

// Group items 0..n-1 by key in [0, n_keys) (key < 0 = skip): counting sort.
// Items of key q are order[start[q] .. start[q + 1]).
struct Groups {
  std::vector<int> start, order;
};

inline Groups group_by(const std::vector<int>& key, int n_keys) {
  Groups g;
  g.start.assign(n_keys + 1, 0);
  for (int q : key) if (q >= 0) ++g.start[q + 1];
  for (int q = 0; q < n_keys; ++q) g.start[q + 1] += g.start[q];
  g.order.resize(g.start[n_keys]);
  std::vector<int> pos(g.start.begin(), g.start.end() - 1);
  for (size_t t = 0; t < key.size(); ++t)
    if (key[t] >= 0) g.order[pos[key[t]]++] = (int)t;
  return g;
}

}  // namespace nle

#endif
