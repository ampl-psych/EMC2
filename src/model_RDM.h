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

[[gnu::always_inline]] inline double pigt0(double t, double k, double l)
{
  // CDF when A==0 (no between trial variability in start point)
  const double lambda = k * k;
  const double mu     = k / l;
  const double sqlt   = std::sqrt(lambda / t);   // shared
  const double tmu    = t / mu;                  // shared
  const double z1     = sqlt * (1.0 + tmu);
  const double z2     = sqlt * (1.0 - tmu);
  return pnorm_upper(z1) * std::exp(2.0 * lambda / mu) + PNORM_STD(z2, false, false);
}

[[gnu::always_inline]] inline double digt0(double t, double k, double l)
{
  // PDF when A==0 (no between trial variability in start point)
  const double lambda = k * k;
  const double tl_k = t * l / k;
  const double e = -0.5 * (lambda / t) * (tl_k - 1.0) * (tl_k - 1.0);
  return std::exp(e) * std::sqrt(lambda / (2.0 * M_PI * t * t * t));
}

// ---------------------------------------------------------------------------
// Core scalar functions — assume t > 0, a >= A_EPS, |l| >= L_EPS
// ---------------------------------------------------------------------------

[[gnu::always_inline]] inline double digt_core(double t, double k, double l, double a)
{
  // PDF when A>0
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

  double pdf = sum / (2.0 * a);
  return pdf;
}

[[gnu::always_inline]] inline double pigt_core(double t, double k, double l, double a)
{
  // CDF when A > 0
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

  return cdf;
}


// ---------------------------------------------------------------------------
// Scalar wrappers — used by R exports.
// Use asymptotic fallback for small a.
// ---------------------------------------------------------------------------

inline double digt(double t, double k, double l, double a)
{
  if (t <= 0.0) return 0.0;
  if (a < A_EPS) return digt0(t, k, l);
  clamp_a_l(a, l);
  double pdf = digt_core(t, k, l, a);
  return (std::isfinite(pdf) && pdf >= 0.0) ? pdf : 0.0;
}

inline double pigt(double t, double k, double l, double a)
{
  if (t <= 0.0) return 0.0;
  if (a < A_EPS) return pigt0(t, k, l);
  clamp_a_l(a, l);
  double cdf = pigt_core(t, k, l, a);
  if (!std::isfinite(cdf) || cdf < 0.0) return 0.0;
  if (cdf > 1.0) return 1.0;
  return cdf;
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
               double* ll_row)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* v  = &pt.base(0, spec.col_v);
  const double* B  = &pt.base(0, spec.col_B);
  const double* A  = &pt.base(0, spec.col_A);
  const double* t0 = &pt.base(0, spec.col_t0);
  const double* s  = &pt.base(0, spec.col_s);

  int* win_ptr = LOGICAL(winner);

  for (int i = 0; i < N; ++i) {
    if (!win_ptr[i]) continue;

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0) { ll_row[i] = 0.0; continue; }

    const double inv_s = 1.0 / s[i];
    double pdf;
    if (A[i] < A_EPS) {
      double l = v[i] * inv_s;
      double k = B[i] * inv_s;
      if (l > -L_EPS && l < L_EPS) l = (l >= 0.0 ? L_EPS : -L_EPS);
      pdf = digt0(t_eff, k, l);
    } else {
      double a = 0.5 * A[i] * inv_s;
      double l = v[i] * inv_s;
      double k = B[i] * inv_s + a;
      clamp_a_l(a, l);
      pdf = digt_core(t_eff, k, l, a);
    }
    ll_row[i] = (std::isfinite(pdf) && pdf >= 0.0) ? pdf : 0.0;
  }
}

void prdm_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* ll_row)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* v  = &pt.base(0, spec.col_v);
  const double* B  = &pt.base(0, spec.col_B);
  const double* A  = &pt.base(0, spec.col_A);
  const double* t0 = &pt.base(0, spec.col_t0);
  const double* s  = &pt.base(0, spec.col_s);

  int* win_ptr = LOGICAL(winner);

  for (int i = 0; i < N; ++i) {
    if (win_ptr[i]) continue;

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0) { ll_row[i] = 0.0; continue; }

    const double inv_s = 1.0 / s[i];
    double cdf;
    if (A[i] < A_EPS) {
      double l = v[i] * inv_s;
      double k = B[i] * inv_s;
      if (l > -L_EPS && l < L_EPS) l = (l >= 0.0 ? L_EPS : -L_EPS);
      cdf = pigt0(t_eff, k, l);
    } else {
      double a = 0.5 * A[i] * inv_s;
      double l = v[i] * inv_s;
      double k = B[i] * inv_s + a;
      clamp_a_l(a, l);
      cdf = pigt_core(t_eff, k, l, a);
    }
    if (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
    else if (cdf > 1.0)                    cdf = 1.0;
    ll_row[i] = cdf;
  }
}


// This new filling function checks whether A==0, if so --> runs digt0 and pigt0
void drdm_prdm_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch)
{
  // Note that for the losers, the *SURVIVAL* probability is filled, *NOT* the CDF
  const double* __restrict__ rt = rts.begin();
  const double* __restrict__ v  = &pt.base(0, spec.col_v);
  const double* __restrict__ B  = &pt.base(0, spec.col_B);
  const double* __restrict__ A  = &pt.base(0, spec.col_A);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
  const double* __restrict__ s  = &pt.base(0, spec.col_s);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // noA scratch (primary)
  double* __restrict__ sc_teff = scratch.t_eff.data();
  double* __restrict__ sc_v    = scratch.v.data();
  double* __restrict__ sc_B    = scratch.B.data();
  double* __restrict__ sc_s    = scratch.s.data();
  double* __restrict__ sc_out  = scratch.out.data();
  int*    __restrict__ sc_idx  = scratch.idx_win0.data();

  // core scratch (secondary)
  double* __restrict__ sc_teff_c = scratch.t_eff_c.data();
  double* __restrict__ sc_v_c    = scratch.v_c.data();
  double* __restrict__ sc_B_c    = scratch.B_c.data();
  double* __restrict__ sc_A_c    = scratch.A_c.data();
  double* __restrict__ sc_s_c    = scratch.s_c.data();
  double* __restrict__ sc_out_c  = scratch.out_c.data();
  int*    __restrict__ sc_idx_c  = scratch.idx_win_c.data();

  // =========================================================================
  // WINNERS
  // =========================================================================

  // --- One-pass gather: split into noA (primary, most likely) and core (secondary, most likely) ---
  int n_win_noA = 0, n_win_core = 0;
  for (int j = 0; j < n_win; ++j) {
    const int i       = idx_win[j];

    // first guard against teff < 0 (non-decision time can't be > rt)
    const double teff = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = 0.0; continue; }

    // Check whether A equals 0
    if (A[i] < A_EPS) {
      sc_teff[n_win_noA] = teff;
      sc_v   [n_win_noA] = v[i];
      sc_B   [n_win_noA] = B[i];
      sc_s   [n_win_noA] = s[i];
      sc_idx [n_win_noA] = i;
      n_win_noA++;
    } else {
      sc_teff_c[n_win_core] = teff;
      sc_v_c   [n_win_core] = v[i];
      sc_B_c   [n_win_core] = B[i];
      sc_A_c   [n_win_core] = A[i];
      sc_s_c   [n_win_core] = s[i];
      sc_idx_c [n_win_core] = i;
      n_win_core++;
    }
  }

  // --- Compute: digt0 winners ---
#pragma omp simd
  for (int j = 0; j < n_win_noA; ++j) {
    const double inv_s = 1.0 / sc_s[j];
    double l = sc_v[j] * inv_s;
    double k = sc_B[j] * inv_s;
    if (l > -L_EPS && l < L_EPS) l = (l >= 0.0 ? L_EPS : -L_EPS);  // still guard l
    sc_out[j] = digt0(sc_teff[j], k, l);
  }

  // --- Compute: digt_core winners ---
#pragma omp simd
  for (int j = 0; j < n_win_core; ++j) {
    const double inv_s = 1.0 / sc_s_c[j];
    double a = 0.5 * sc_A_c[j] * inv_s;
    double l = sc_v_c[j]       * inv_s;
    double k = sc_B_c[j]       * inv_s + a;
    clamp_a_l(a, l); // clamp a and l
    sc_out_c[j] = digt_core(sc_teff_c[j], k, l, a);
  }

  // --- Scatter: noA winners ---
  for (int j = 0; j < n_win_noA; ++j) {
    const double val = sc_out[j];
    ll_row[sc_idx[j]] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;    // guard against pdf < 0
  }

  // --- Scatter: core winners ---
  for (int j = 0; j < n_win_core; ++j) {
    const double val = sc_out_c[j];
    ll_row[sc_idx_c[j]] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;  // guard against pdf < 0
  }

  // =========================================================================
  // LOSERS
  // =========================================================================

  sc_idx  = scratch.idx_los0.data();
  sc_idx_c = scratch.idx_los_c.data();

  // --- One-pass gather: split into noA (primary, most likely) and core (secondary, most likely) ---
  int n_los_noA = 0, n_los_core = 0;
  for (int j = 0; j < n_los; ++j) {
    const int i       = idx_los[j];

    // first guard against teff < 0 (non-decision time can't be > rt)
    const double teff = rt[i] - t0[i];
    if (teff <= 0.0) { ll_row[i] = 1.0; continue; }

    // fill scratch depending on whether we have A tiny or not
    if (A[i] < A_EPS) {
      sc_teff[n_los_noA] = teff;
      sc_v   [n_los_noA] = v[i];
      sc_B   [n_los_noA] = B[i];
      sc_s   [n_los_noA] = s[i];
      sc_idx [n_los_noA] = i;
      n_los_noA++;
    } else {
      sc_teff_c[n_los_core] = teff;
      sc_v_c   [n_los_core] = v[i];
      sc_B_c   [n_los_core] = B[i];
      sc_A_c   [n_los_core] = A[i];
      sc_s_c   [n_los_core] = s[i];
      sc_idx_c [n_los_core] = i;
      n_los_core++;
    }
  }

  // --- Compute: pigt0 losers ---
#pragma omp simd
  for (int j = 0; j < n_los_noA; ++j) {
    const double inv_s = 1.0 / sc_s[j];
    double l = sc_v[j] * inv_s;
    double k = sc_B[j] * inv_s;
    if (l > -L_EPS && l < L_EPS) l = (l >= 0.0 ? L_EPS : -L_EPS);
    sc_out[j] = pigt0(sc_teff[j], k, l);
  }

  // --- Compute: pigt_core losers ---
#pragma omp simd
  for (int j = 0; j < n_los_core; ++j) {
    const double inv_s = 1.0 / sc_s_c[j];
    double a = 0.5 * sc_A_c[j] * inv_s;
    double l = sc_v_c[j]       * inv_s;
    double k = sc_B_c[j]       * inv_s + a;
    clamp_a_l(a, l);
    sc_out_c[j] = pigt_core(sc_teff_c[j], k, l, a);
  }

  // --- Scatter: noA losers ---
  for (int j = 0; j < n_los_noA; ++j) {
    double val = sc_out[j];
    if (!std::isfinite(val) || val < 0.0) val = 0.0;    // guard against cdf < 0
    else if (val > 1.0)                    val = 1.0;   // guard against cdf > 1
    ll_row[sc_idx[j]] = 1.0 - val;                      // fill in survival (not CDF)
  }

  // --- Scatter: core losers ---
  for (int j = 0; j < n_los_core; ++j) {
    double val = sc_out_c[j];
    if (!std::isfinite(val) || val < 0.0) val = 0.0;    // guard against cdf < 0
    else if (val > 1.0)                    val = 1.0;   // guard against cdf > 1
    ll_row[sc_idx_c[j]] = 1.0 - val;                    // fill in survival (not CDF)
  }
}



// void drdm_prdm_fast(const NumericVector& rts,
//                     const ParamTable& pt,
//                     const RaceSpec& spec,
//                     const std::vector<int>& idx_win,
//                     const std::vector<int>& idx_los,
//                     double* __restrict__ ll_row,
//                     RaceScratch& scratch)
// {
//   const double* __restrict__ rt = rts.begin();
//   const double* __restrict__ v  = &pt.base(0, spec.col_v);
//   const double* __restrict__ B  = &pt.base(0, spec.col_B);
//   const double* __restrict__ A  = &pt.base(0, spec.col_A);
//   const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
//   const double* __restrict__ s  = &pt.base(0, spec.col_s);
//
//   const int n_win = (int)idx_win.size();
//   const int n_los = (int)idx_los.size();
//
//   double* __restrict__ sc_teff = scratch.t_eff.data();
//   double* __restrict__ sc_v    = scratch.v.data();
//   double* __restrict__ sc_B    = scratch.B.data();
//   double* __restrict__ sc_A    = scratch.A.data();
//   double* __restrict__ sc_s    = scratch.s.data();
//   double* __restrict__ sc_out  = scratch.out.data();
//   int*    __restrict__ sc_idx  = scratch.idx_win0.data();
//
//   // --- Winners: gather (filter t_eff <= 0) ---
//   int n_win_valid = 0;
//   for (int j = 0; j < n_win; ++j) {
//     const int i       = idx_win[j];
//     const double teff = rt[i] - t0[i];
//     if (teff <= 0.0) { ll_row[i] = 0.0; continue; }
//     sc_teff[n_win_valid] = teff;
//     sc_v   [n_win_valid] = v[i];
//     sc_B   [n_win_valid] = B[i];
//     sc_A   [n_win_valid] = A[i];
//     sc_s   [n_win_valid] = s[i];
//     sc_idx [n_win_valid] = i;
//     n_win_valid++;
//   }
//
//   // --- Winners: compute ---
// #pragma omp simd
//   for (int j = 0; j < n_win_valid; ++j) {
//     const double inv_s = 1.0 / sc_s[j];
//     double a = 0.5 * sc_A[j] * inv_s;
//     double l = sc_v[j]       * inv_s;
//     double k = sc_B[j]       * inv_s + a;
//     clamp_a_l(a, l);
//     sc_out[j] = digt_core(sc_teff[j], k, l, a);
//   }
//
//   // --- Winners: scatter and guard---
//   for (int j = 0; j < n_win_valid; ++j) {
//     const double v = sc_out[j];
//     ll_row[sc_idx[j]] = (std::isfinite(v) && v >= 0.0) ? v : 0.0;   // guard against division by 0
//   }
//
//   sc_idx = scratch.idx_los0.data();
//
//   // --- Losers: gather (filter t_eff <= 0) ---
//   int n_los_valid = 0;
//   for (int j = 0; j < n_los; ++j) {
//     const int i       = idx_los[j];
//     const double teff = rt[i] - t0[i];
//     if (teff <= 0.0) { ll_row[i] = 1.0; continue; }
//     sc_teff[n_los_valid] = teff;
//     sc_v   [n_los_valid] = v[i];
//     sc_B   [n_los_valid] = B[i];
//     sc_A   [n_los_valid] = A[i];
//     sc_s   [n_los_valid] = s[i];
//     sc_idx [n_los_valid] = i;
//     n_los_valid++;
//   }
//
//   // --- Losers: compute ---
// #pragma omp simd
//   for (int j = 0; j < n_los_valid; ++j) {
//     const double inv_s = 1.0 / sc_s[j];
//     double a = 0.5 * sc_A[j] * inv_s;
//     double l = sc_v[j]       * inv_s;
//     double k = sc_B[j]       * inv_s + a;
//     clamp_a_l(a, l);
//     sc_out[j] = pigt_core(sc_teff[j], k, l, a);
//   }
//
//   // --- Losers: scatter ---
//   for (int j = 0; j < n_los_valid; ++j) {
//     double v = sc_out[j];
//     if (!std::isfinite(v) || v < 0.0) v = 0.0;
//     else if (v > 1.0)                  v = 1.0;
//     ll_row[sc_idx[j]] = 1.0 - v;
//   }
// }
//
// void drdm_prdm_noA_fast(const NumericVector& rts,
//                         const ParamTable& pt,
//                         const RaceSpec& spec,
//                         const std::vector<int>& idx_win,
//                         const std::vector<int>& idx_los,
//                         double* __restrict__ ll_row,
//                         RaceScratch& scratch)
// {
//   const double* __restrict__ rt = rts.begin();
//   const double* __restrict__ v  = &pt.base(0, spec.col_v);
//   const double* __restrict__ B  = &pt.base(0, spec.col_B);
//   const double* __restrict__ t0 = &pt.base(0, spec.col_t0);
//   const double* __restrict__ s  = &pt.base(0, spec.col_s);
//
//   const int n_win = (int)idx_win.size();
//   const int n_los = (int)idx_los.size();
//
//   double* __restrict__ sc_teff = scratch.t_eff.data();
//   double* __restrict__ sc_v    = scratch.v.data();
//   double* __restrict__ sc_B    = scratch.B.data();
//   double* __restrict__ sc_s    = scratch.s.data();
//   double* __restrict__ sc_out  = scratch.out.data();
//   int*    __restrict__ sc_idx  = scratch.idx_win0.data();
//
//   // --- Winners: gather (filter t_eff <= 0) ---
//   int n_win_valid = 0;
//   for (int j = 0; j < n_win; ++j) {
//     const int i       = idx_win[j];
//     const double teff = rt[i] - t0[i];
//     if (teff <= 0.0) { ll_row[i] = 0.0; continue; }
//     sc_teff[n_win_valid] = teff;
//     sc_v   [n_win_valid] = v[i];
//     sc_B   [n_win_valid] = B[i];
//     sc_s   [n_win_valid] = s[i];
//     sc_idx [n_win_valid] = i;
//     n_win_valid++;
//   }
//
//   // --- Winners: compute ---
// #pragma omp simd
//   for (int j = 0; j < n_win_valid; ++j) {
//     const double inv_s = 1.0 / sc_s[j];
//     double l = sc_v[j] * inv_s;
//     double k = sc_B[j] * inv_s;
//     if (l > -L_EPS && l < L_EPS) l = (l >= 0.0 ? L_EPS : -L_EPS);   // guard against division by 0
//     sc_out[j] = digt0(sc_teff[j], k, l);
//   }
//
//   // --- Winners: scatter ---
//   for (int j = 0; j < n_win_valid; ++j) {
//     const double v = sc_out[j];
//     ll_row[sc_idx[j]] = (std::isfinite(v) && v >= 0.0) ? v : 0.0;
//   }
//
//   sc_idx = scratch.idx_los0.data();
//
//   // --- Losers: gather (filter t_eff <= 0) ---
//   int n_los_valid = 0;
//   for (int j = 0; j < n_los; ++j) {
//     const int i       = idx_los[j];
//     const double teff = rt[i] - t0[i];
//     if (teff <= 0.0) { ll_row[i] = 1.0; continue; }
//     sc_teff[n_los_valid] = teff;
//     sc_v   [n_los_valid] = v[i];
//     sc_B   [n_los_valid] = B[i];
//     sc_s   [n_los_valid] = s[i];
//     sc_idx [n_los_valid] = i;
//     n_los_valid++;
//   }
//
//   // --- Losers: compute ---
// #pragma omp simd
//   for (int j = 0; j < n_los_valid; ++j) {
//     const double inv_s = 1.0 / sc_s[j];
//     double l = sc_v[j] * inv_s;
//     double k = sc_B[j] * inv_s;
//     if (l > -L_EPS && l < L_EPS) l = (l >= 0.0 ? L_EPS : -L_EPS);   // guard against division by 0
//     sc_out[j] = pigt0(sc_teff[j], k, l);
//   }
//
//   // --- Losers: scatter ---
//   for (int j = 0; j < n_los_valid; ++j) {
//     double v = sc_out[j];
//     if (!std::isfinite(v) || v < 0.0) v = 0.0;
//     else if (v > 1.0)                 v = 1.0;
//     ll_row[sc_idx[j]] = 1.0 - v;
//   }
// }





// // [[Rcpp::export]]
// double bench_digt_core_vec(Rcpp::NumericVector t,
//                            Rcpp::NumericVector k,
//                            Rcpp::NumericVector l,
//                            Rcpp::NumericVector a)
// {
//   const int n = t.size();
//   // For safety, use the smallest length of the inputs
//   const int N = std::min(std::min(k.size(), l.size()), a.size());
//
//   double acc = 0.0;
//   for (int i = 0; i < N; ++i) {
//     acc += digt_core(t[i], k[i], l[i], a[i]);
//   }
//   return acc;
// }
//
// // [[Rcpp::export]]
// double bench_digt0_vec(Rcpp::NumericVector t,
//                        Rcpp::NumericVector k,
//                        Rcpp::NumericVector l)
// {
//   const int n = t.size();
//   const int N = std::min(k.size(), l.size());
//
//   double acc = 0.0;
//   for (int i = 0; i < N; ++i) {
//     acc += digt0(t[i], k[i], l[i]);
//   }
//   return acc;
// }
//
// // [[Rcpp::export]]
// double bench_pigt_core_vec(Rcpp::NumericVector t,
//                            Rcpp::NumericVector k,
//                            Rcpp::NumericVector l,
//                            Rcpp::NumericVector a)
// {
//   const int n = t.size();
//   const int N = std::min(std::min(k.size(), l.size()), a.size());
//
//   double acc = 0.0;
//   for (int i = 0; i < N; ++i) {
//     acc += pigt_core(t[i], k[i], l[i], a[i]);
//   }
//   return acc;
// }
//
// // [[Rcpp::export]]
// double bench_pigt0_vec(Rcpp::NumericVector t,
//                        Rcpp::NumericVector k,
//                        Rcpp::NumericVector l)
// {
//   const int n = t.size();
//   const int N = std::min(k.size(), l.size());
//
//   double acc = 0.0;
//   for (int i = 0; i < N; ++i) {
//     acc += pigt0(t[i], k[i], l[i]);
//   }
//   return acc;
// }

#endif // rdm_h
