#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "RaceSpec.h"
#include "math_utils.h"  // must be before Rcpp
#include "wald_functions.h"
#include "ParamTable.h"

using namespace Rcpp;

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
      clamp_l(l);
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
      clamp_l(l);
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
      const double k = B[i] / s[i];
      if (k > K_MAX || k < 0.0) { ll_row[i] = 0.0; continue; }  // degenerate
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
    clamp_l(l);
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
      const double k = B[i] / s[i];
      if (k > K_MAX || k < 0.0) { ll_row[i] = 1.0; continue; }  // survival = 1
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
    clamp_l(l);
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
