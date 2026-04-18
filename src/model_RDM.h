#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "RaceSpec.h"
#include "math_utils.h"  // must be before Rcpp
#include "pnorm_utils.h"
#include "ParamTable.h"

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Numerical stability clamps — applied in the fast path only.
// Guarantees a >= A_EPS and |l| >= L_EPS, eliminating degenerate branches
// in the core functions. The scalar wrappers use the asymptotic fallbacks
// instead, preserving numerical accuracy for small a.
// ---------------------------------------------------------------------------

constexpr double A_EPS = 1e-4;
constexpr double L_EPS = 1e-4;

inline void clamp_a_l(double& a, double& l)
{
  a = (a < A_EPS) ? A_EPS : a;
  l = (l > -L_EPS && l < L_EPS) ? (l >= 0.0 ? L_EPS : -L_EPS) : l;
}

// ---------------------------------------------------------------------------
// Asymptotic formulas for a -> 0 (point-mass starting point)
// ---------------------------------------------------------------------------

inline double pigt0(double t, double k, double l)
{
  const double mu     = k / l;
  const double lambda = k * k;
  const double z1     = std::sqrt(lambda / t) * (1.0 + t / mu);
  const double z2     = std::sqrt(lambda / t) * (1.0 - t / mu);
  return std::exp(2.0 * lambda / mu + std::log(PNORM_STD(z1, false, false))) + PNORM_STD(z2, false, false);
}

inline double digt0(double t, double k, double l)
{
  const double lambda = k * k;
  const double e = (l == 0.0) ? -0.5 * lambda / t : -(lambda / (2.0 * t)) * ((t * t * l * l) / (k * k) - 2.0 * t * l / k + 1.0);
  return std::exp(e + 0.5 * std::log(lambda) - 0.5 * std::log(2.0 * t * t * t * M_PI));
}

// ---------------------------------------------------------------------------
// Core scalar functions — assume t > 0, a >= A_EPS, |l| >= L_EPS
// ---------------------------------------------------------------------------

inline double digt_core(double t, double k, double l, double a)
{
  if (t <= 0.0) {
    return 0.0;
  }

  const double sqt      = std::sqrt(t);
  const double inv_sqt  = 1.0 / sqt;
  const double inv_t    = 1.0 / t;
  const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);

  // t1 part – same structure/order as in the old code
  const double temp1 = a - k + t * l;
  const double temp2 = a + k - t * l;
  const double t1a   = -0.5 * temp1 * temp1 * inv_t;
  const double t1b   = -0.5 * temp2 * temp2 * inv_t;
  const double t1    = inv_sqrt_2pi * (std::exp(t1a) - std::exp(t1b)) * inv_sqt;

  // t2 part – same structure/order as in the old code
  const double arg1 = (-k + a) * inv_sqt + sqt * l;
  const double arg2 = ( k + a) * inv_sqt - sqt * l;

  const double t2a = 2.0 * PNORM_STD(arg1, /*lower=*/true, /*logp=*/false) - 1.0;
  const double t2b = 2.0 * PNORM_STD(arg2, /*lower=*/true, /*logp=*/false) - 1.0;
  // const double t2a = std::erf(arg1 * M_SQRT1_2);
  // const double t2b = std::erf(arg2 * M_SQRT1_2);
  const double t2  = 0.5 * l * (t2a + t2b);

  const double sum = t1 + t2;

  if (sum <= 0.0 || !std::isfinite(sum)) {
    return 0.0;
  }

  double pdf = sum / (2.0 * a);

  if (!std::isfinite(pdf) || pdf < 0.0) {
    return 0.0;
  }
  return pdf;
}

inline double pigt_core(double t, double k, double l, double a)
{
  if (t <= 0.0) {
    return 0.0;
  }

  const double sqt      = std::sqrt(t);
  const double inv_sqt  = 1.0 / sqt;
  const double inv_t    = 1.0 / t;
  const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);

  // t1 term: sqt / sqrt(2π) * (exp(...) - exp(...)) – same order as old code
  const double tmp1 = k - a - t * l;
  const double tmp2 = a + k - t * l;
  const double t1a  = std::exp(-0.5 * tmp1 * tmp1 * inv_t);
  const double t1b  = std::exp(-0.5 * tmp2 * tmp2 * inv_t);
  const double t1   = sqt * inv_sqrt_2pi * (t1a - t1b);

  // t2 term – same structure/order as in the old code
  const double argA = -(k - a + t * l) * inv_sqt;
  const double argB = -(k + a + t * l) * inv_sqt;

  const double t2a = std::exp(2.0 * l * (k - a) +
                              PNORM_STD(argA, /*lower=*/true, /*logp=*/true));
  const double t2b = std::exp(2.0 * l * (k + a) +
                              PNORM_STD(argB, /*lower=*/true, /*logp=*/true));
  const double t2  = a + (t2b - t2a) / (2.0 * l);

  // t4 term – same structure/order as in the old code
  const double t4a = 2.0 * PNORM_STD((k + a) * inv_sqt - sqt * l,
                                     /*lower=*/true, /*logp=*/false) - 1.0;
  const double t4b = 2.0 * PNORM_STD((k - a) * inv_sqt - sqt * l,
                                     /*lower=*/true, /*logp=*/false) - 1.0;
  //  equivalent but no pnorm
  // const double t4a = std::erf((k + a - t * l) / (sqt * M_SQRT2));
  // const double t4b = std::erf((k - a - t * l) / (sqt * M_SQRT2));
  const double t4  = 0.5 * (t * l - a - k + 0.5 / l) * t4a + 0.5 * (k - a - t * l - 0.5 / l) * t4b;

  double cdf = 0.5 * (t4 + t2 + t1) / a;

  if (!std::isfinite(cdf) || cdf < 0.0) {
    return 0.0;
  }
  if (cdf > 1.0) {
    return 1.0;
  }
  return cdf;
}

// ---------------------------------------------------------------------------
// Scalar wrappers — used by legacy matrix-based functions and R exports.
// Use asymptotic fallback for small a, matching original behaviour exactly.
// ---------------------------------------------------------------------------

constexpr double A_ASYMPTOTIC = 1e-10;

inline double digt(double t, double k, double l, double a)
{
  if (t <= 0.0) return 0.0;
  if (a < A_ASYMPTOTIC) return digt0(t, k, l);
  clamp_a_l(a, l);
  return digt_core(t, k, l, a);
}

inline double pigt(double t, double k, double l, double a)
{
  if (t <= 0.0) return 0.0;
  if (a < A_ASYMPTOTIC) return pigt0(t, k, l);
  clamp_a_l(a, l);
  return pigt_core(t, k, l, a);
}

// ---------------------------------------------------------------------------
// Legacy matrix-based functions
// ---------------------------------------------------------------------------

NumericVector drdm_c(NumericVector rts, NumericMatrix pars,
                     LogicalVector idx, double min_ll, LogicalVector is_ok)
{
  // v=0, B=1, A=2, t0=3, s=4
  NumericVector out(sum(idx));
  int k = 0;
  for (int i = 0; i < rts.length(); i++) {
    if (idx[i]) {
      if (NumericVector::is_na(pars(i, 0))) {
        out[k] = 0.0;
      } else if ((rts[i] - pars(i, 3) > 0) && is_ok[i]) {
        double inv_s = 1.0 / pars(i, 4);
        out[k] = digt(rts[i] - pars(i, 3),
                      pars(i, 1) * inv_s + 0.5 * pars(i, 2) * inv_s,
                      pars(i, 0) * inv_s,
                      0.5 * pars(i, 2) * inv_s);
      } else {
        out[k] = min_ll;
      }
      k++;
    }
  }
  return out;
}

NumericVector prdm_c(NumericVector rts, NumericMatrix pars,
                     LogicalVector idx, double min_ll, LogicalVector is_ok)
{
  // v=0, B=1, A=2, t0=3, s=4
  NumericVector out(sum(idx));
  int k = 0;
  for (int i = 0; i < rts.length(); i++) {
    if (idx[i]) {
      if (NumericVector::is_na(pars(i, 0))) {
        out[k] = 0.0;
      } else if ((rts[i] - pars(i, 3) > 0) && is_ok[i]) {
        double inv_s = 1.0 / pars(i, 4);
        out[k] = pigt(rts[i] - pars(i, 3),
                      pars(i, 1) * inv_s + 0.5 * pars(i, 2) * inv_s,
                      pars(i, 0) * inv_s,
                      0.5 * pars(i, 2) * inv_s);
      } else {
        out[k] = min_ll;
      }
      k++;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// R-exported scalar functions
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0)
{
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    double t_eff = t[i] - t0[i];
    pdf[i] = (t_eff <= 0.0) ? 0.0 : digt(t_eff, B[i] + 0.5 * A[i], v[i], 0.5 * A[i]);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0)
{
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++) {
    double t_eff = t[i] - t0[i];
    cdf[i] = (t_eff <= 0.0) ? 0.0 : pigt(t_eff, B[i] + 0.5 * A[i], v[i], 0.5 * A[i]);
  }
  return cdf;
}

// ---------------------------------------------------------------------------
// Fast ParamTable-based functions
// ---------------------------------------------------------------------------

void drdm_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* raw)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* v  = &pt.base(0, spec.col_v);
  const double* B  = &pt.base(0, spec.col_B);
  const double* A  = &pt.base(0, spec.col_A);
  const double* t0 = &pt.base(0, spec.col_t0);
  const double* s  = &pt.base(0, spec.col_s);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (!win_ptr[i]) continue;

    const double t_eff = rt[i] - t0[i];
    const double inv_s = 1.0 / s[i];
    double a = 0.5 * A[i] * inv_s;
    double l = v[i] * inv_s;
    double k = B[i] * inv_s + a;

    clamp_a_l(a, l);
    raw[i] = digt_core(t_eff, k, l, a);
  }
}

void prdm_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* raw)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* v  = &pt.base(0, spec.col_v);
  const double* B  = &pt.base(0, spec.col_B);
  const double* A  = &pt.base(0, spec.col_A);
  const double* t0 = &pt.base(0, spec.col_t0);
  const double* s  = &pt.base(0, spec.col_s);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (win_ptr[i]) continue;

    const double t_eff = rt[i] - t0[i];
    const double inv_s = 1.0 / s[i];
    double a = 0.5 * A[i] * inv_s;
    double l = v[i] * inv_s;
    double k = B[i] * inv_s + a;

    clamp_a_l(a, l);
    raw[i] = pigt_core(t_eff, k, l, a);
  }
}


void drdm_prdm_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ raw,
                    RaceScratch& scratch)
{
  const double* __restrict__ rt = rts.begin();
  const double* __restrict__ v  = &pt.base(0, spec.col_v);
  const double* __restrict__ B  = &pt.base(0, spec.col_B);
  const double* __restrict__ A  = &pt.base(0, spec.col_A);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
  const double* __restrict__ s  = &pt.base(0, spec.col_s);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // Restrict-qualified pointers into scratch — lets clang prove no aliasing
  // with the input arrays during gather.
  double* __restrict__ sc_teff = scratch.t_eff.data();
  double* __restrict__ sc_v    = scratch.v.data();
  double* __restrict__ sc_B    = scratch.B.data();
  double* __restrict__ sc_A    = scratch.A.data();
  double* __restrict__ sc_s    = scratch.s.data();
  double* __restrict__ sc_out  = scratch.out.data();

  // --- Winners: gather ---
  for (int j = 0; j < n_win; ++j) {
    const int i  = idx_win[j];
    sc_teff[j]   = rt[i] - t0[i];
    sc_v[j]      = v[i];
    sc_B[j]      = B[i];
    sc_A[j]      = A[i];
    sc_s[j]      = s[i];
  }

  // --- Winners: compute (contiguous — cost-model override via simd) ---
#pragma omp simd
  for (int j = 0; j < n_win; ++j) {
    const double inv_s = 1.0 / sc_s[j];
    double a = 0.5 * sc_A[j] * inv_s;
    double l = sc_v[j]       * inv_s;
    double k = sc_B[j]       * inv_s + a;
    clamp_a_l(a, l);
    sc_out[j] = digt_core(sc_teff[j], k, l, a);
  }

  // --- Winners: scatter ---
  for (int j = 0; j < n_win; ++j) raw[idx_win[j]] = sc_out[j];

  // --- Losers: gather ---
  for (int j = 0; j < n_los; ++j) {
    const int i  = idx_los[j];
    sc_teff[j]   = rt[i] - t0[i];
    sc_v[j]      = v[i];
    sc_B[j]      = B[i];
    sc_A[j]      = A[i];
    sc_s[j]      = s[i];
  }

  // --- Losers: compute (contiguous — cost-model override via simd) ---
#pragma omp simd
  for (int j = 0; j < n_los; ++j) {
    const double inv_s = 1.0 / sc_s[j];
    double a = 0.5 * sc_A[j] * inv_s;
    double l = sc_v[j]       * inv_s;
    double k = sc_B[j]       * inv_s + a;
    clamp_a_l(a, l);
    sc_out[j] = pigt_core(sc_teff[j], k, l, a);
  }

  // --- Losers: scatter ---
  // fill in 1-CDF - survival!
  for (int j = 0; j < n_los; ++j) raw[idx_los[j]] = 1-sc_out[j];
}

#endif // rdm_h
