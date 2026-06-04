#ifndef lba_h
#define lba_h

#include <Rcpp.h>
#include "RaceSpec.h"
#include "utility_functions.h"
#include "math_utils.h"
#include "pnorm_utils.h"
#include "ParamTable.h"

using namespace Rcpp;

constexpr double LBA_A_ASYMPTOTIC = 1e-10;


// ---------------------------------------------------------------------------
// Split compute functions — guards removed, called from compute loops
// ---------------------------------------------------------------------------

[[gnu::always_inline]] inline double dlba_core(double t, double A, double b,
                                               double v, double sv, double denom)
{
  const double zs     = t * sv;
  const double cmz    = b - t * v;
  const double cz     = cmz / zs;
  const double cz_max = (cmz - A) / zs;

  return (v * (PNORM_STD(cz,     /*lower=*/true, /*logp=*/false)
                 - PNORM_STD(cz_max, /*lower=*/true, /*logp=*/false))
            + sv * (DNORM_STD(cz_max) - DNORM_STD(cz))) / (A * denom);
}

[[gnu::always_inline]] inline double dlba_noA(double t, double b,
                                              double v, double sv, double denom)
{
  return DNORM(b / t, v, sv) * b / (t * t * denom);
}

[[gnu::always_inline]] inline double plba_core(double t, double A, double b,
                                               double v, double sv, double denom)
{
  const double zs     = t * sv;
  const double cmz    = b - t * v;
  const double xx     = cmz - A;
  const double cz     = cmz / zs;
  const double cz_max = xx / zs;

  return (1.0 + (zs * (DNORM_STD(cz_max) - DNORM_STD(cz))
                   + xx  * PNORM_STD(cz_max, /*lower=*/true, /*logp=*/false)
                   - cmz * PNORM_STD(cz,     /*lower=*/true, /*logp=*/false)) / A) / denom;
}

[[gnu::always_inline]] inline double plba_noA(double t, double b,
                                              double v, double sv, double denom)
{
  return PNORM_STD((b / t - v) / sv, /*lower=*/false, /*logp=*/false) / denom;
}

// Old pipeline, not sure whether to keep this (probably useful for censoring/truncation)
void plba_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* ll_row)
{
  const int N = rts.size();

  const double* rt  = rts.begin();
  const double* v   = &pt.base(0, spec.col_v);
  const double* sv  = &pt.base(0, spec.col_sv);
  const double* B   = &pt.base(0, spec.col_B);
  const double* A   = &pt.base(0, spec.col_A);
  const double* t0  = &pt.base(0, spec.col_t0);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (win_ptr[i]) continue;

    if (std::isnan(v[i])) { ll_row[i] = 0.0; continue; }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0) { ll_row[i] = 0.0; continue; }

    double denom = PNORM_STD(v[i] / sv[i], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;

    const double A_i = A[i];
    const double b_i = B[i] + A_i;

    // use no-A branch or A branch?
    double cdf = (A_i > LBA_A_ASYMPTOTIC) ? plba_core(t_eff, A_i, b_i, v[i], sv[i], denom) : plba_noA (t_eff,      b_i, v[i], sv[i], denom);

    if      (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
    else if (cdf > 1.0)                         cdf = 1.0;

    ll_row[i] = cdf;
  }
}

// Old pipeline, not sure whether to keep this (probably useful for censoring/truncation)
void dlba_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* ll_row)
{
  const int N = rts.size();

  const double* rt  = rts.begin();
  const double* v   = &pt.base(0, spec.col_v);
  const double* sv  = &pt.base(0, spec.col_sv);
  const double* B   = &pt.base(0, spec.col_B);
  const double* A   = &pt.base(0, spec.col_A);
  const double* t0  = &pt.base(0, spec.col_t0);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (!win_ptr[i]) continue;

    if (std::isnan(v[i])) { ll_row[i] = 0.0; continue; }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0) { ll_row[i] = 0.0; continue; }

    double denom = PNORM_STD(v[i] / sv[i], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;

    const double A_i = A[i];
    const double b_i = B[i] + A_i;

    // use no-A branch or A branch?
    double pdf = (A_i > LBA_A_ASYMPTOTIC) ? dlba_core(t_eff, A_i, b_i, v[i], sv[i], denom) : dlba_noA (t_eff,      b_i, v[i], sv[i], denom);

    ll_row[i] = (std::isfinite(pdf) && pdf >= 0.0) ? pdf : 0.0;
  }
}



// ---------------------------------------------------------------------------
// Hot path: one-pass gather split on A, core dominates -> primary scratch
// ---------------------------------------------------------------------------
void dlba_plba_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch)
{
  // Note that for the losers, the *SURVIVAL* probability is filled, *NOT* the CDF
  const double* __restrict__ rt  = rts.begin();
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ B   = &pt.base(0, spec.col_B);
  const double* __restrict__ A   = &pt.base(0, spec.col_A);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // primary scratch — core entries (A > LBA_A_ASYMPTOTIC)
  double* __restrict__ sc_teff  = scratch.t_eff.data();
  double* __restrict__ sc_v     = scratch.v.data();
  double* __restrict__ sc_sv    = scratch.sv.data();
  double* __restrict__ sc_B     = scratch.B.data();
  double* __restrict__ sc_A     = scratch.A.data();
  double* __restrict__ sc_out   = scratch.out.data();
  int*    __restrict__ sc_idx   = scratch.idx_win0.data();

  // secondary scratch — noA entries (A <= LBA_A_ASYMPTOTIC)
  // note: sc_sv_c reuses s_c since noA path needs sv but not s
  double* __restrict__ sc_teff_c = scratch.t_eff_c.data();
  double* __restrict__ sc_v_c    = scratch.v_c.data();
  double* __restrict__ sc_sv_c   = scratch.s_c.data();
  double* __restrict__ sc_B_c    = scratch.B_c.data();
  double* __restrict__ sc_out_c  = scratch.out_c.data();
  int*    __restrict__ sc_idx_c  = scratch.idx_win_c.data();

  // =========================================================================
  // WINNERS
  // =========================================================================

  // --- One-pass gather: core -> primary, noA -> secondary ---
  int n_win_core = 0, n_win_noA = 0;
  for (int j = 0; j < n_win; ++j) {
    const int i       = idx_win[j];
    const double teff = rt[i] - t0[i];
    if (std::isnan(v[i]) || teff <= 0.0) { ll_row[i] = 0.0; continue; }
    if (A[i] > LBA_A_ASYMPTOTIC) {
      sc_teff[n_win_core] = teff;
      sc_v   [n_win_core] = v[i];
      sc_sv  [n_win_core] = sv[i];
      sc_B   [n_win_core] = B[i];
      sc_A   [n_win_core] = A[i];
      sc_idx [n_win_core] = i;
      n_win_core++;
    } else {
      sc_teff_c[n_win_noA] = teff;
      sc_v_c   [n_win_noA] = v[i];
      sc_sv_c  [n_win_noA] = sv[i];
      sc_B_c   [n_win_noA] = B[i];
      sc_idx_c [n_win_noA] = i;
      n_win_noA++;
    }
  }

  // --- Compute: dlba_core winners ---
  for (int j = 0; j < n_win_core; ++j) {
    double denom = PNORM_STD(sc_v[j] / sc_sv[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = dlba_core(sc_teff[j], sc_A[j], sc_B[j] + sc_A[j],
                           sc_v[j], sc_sv[j], denom);
    sc_out[j] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;
  }

  // --- Compute: dlba_noA winners ---
  for (int j = 0; j < n_win_noA; ++j) {
    double denom = PNORM_STD(sc_v_c[j] / sc_sv_c[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = dlba_noA(sc_teff_c[j], sc_B_c[j], sc_v_c[j], sc_sv_c[j], denom);
    sc_out_c[j] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;
  }

  // --- Scatter ---
  for (int j = 0; j < n_win_core; ++j) ll_row[sc_idx  [j]] = sc_out  [j];
  for (int j = 0; j < n_win_noA;  ++j) ll_row[sc_idx_c[j]] = sc_out_c[j];

  // =========================================================================
  // LOSERS
  // =========================================================================

  sc_idx   = scratch.idx_los0.data();
  sc_idx_c = scratch.idx_los_c.data();

  // --- One-pass gather: core -> primary, noA -> secondary ---
  int n_los_core = 0, n_los_noA = 0;
  for (int j = 0; j < n_los; ++j) {
    const int i       = idx_los[j];
    const double teff = rt[i] - t0[i];
    if (std::isnan(v[i]) || teff <= 0.0) { ll_row[i] = 1.0; continue; }
    if (A[i] > LBA_A_ASYMPTOTIC) {
      sc_teff[n_los_core] = teff;
      sc_v   [n_los_core] = v[i];
      sc_sv  [n_los_core] = sv[i];
      sc_B   [n_los_core] = B[i];
      sc_A   [n_los_core] = A[i];
      sc_idx [n_los_core] = i;
      n_los_core++;
    } else {
      sc_teff_c[n_los_noA] = teff;
      sc_v_c   [n_los_noA] = v[i];
      sc_sv_c  [n_los_noA] = sv[i];
      sc_B_c   [n_los_noA] = B[i];
      sc_idx_c [n_los_noA] = i;
      n_los_noA++;
    }
  }

  // --- Compute: plba_core losers ---
  for (int j = 0; j < n_los_core; ++j) {
    double denom = PNORM_STD(sc_v[j] / sc_sv[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = plba_core(sc_teff[j], sc_A[j], sc_B[j] + sc_A[j],
                           sc_v[j], sc_sv[j], denom);
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    sc_out[j] = 1.0 - val;
  }

  // --- Compute: plba_noA losers ---
  for (int j = 0; j < n_los_noA; ++j) {
    double denom = PNORM_STD(sc_v_c[j] / sc_sv_c[j], /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
    double val = plba_noA(sc_teff_c[j], sc_B_c[j], sc_v_c[j], sc_sv_c[j], denom);
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    sc_out_c[j] = 1.0 - val;
  }

  // --- Scatter ---
  for (int j = 0; j < n_los_core; ++j) ll_row[sc_idx  [j]] = sc_out  [j];
  for (int j = 0; j < n_los_noA;  ++j) ll_row[sc_idx_c[j]] = sc_out_c[j];
}



// R exports

// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true)
{
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    double denom = 1.0;
    if (posdrift) {
      denom = PNORM_STD(v[i] / sv[i], /*lower=*/true, /*logp=*/false);
      if (denom < 1e-10) denom = 1e-10;
    }
    double val = (A[i] > LBA_A_ASYMPTOTIC)
      ? dlba_core(t[i], A[i], b[i], v[i], sv[i], denom)
        : dlba_noA (t[i],       b[i], v[i], sv[i], denom);
    if (val < 0.0) val = 0.0;
    pdf[i] = val;
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector plba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true)
{
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++) {
    double denom = 1.0;
    if (posdrift) {
      denom = PNORM_STD(v[i] / sv[i], /*lower=*/true, /*logp=*/false);
      if (denom < 1e-10) denom = 1e-10;
    }
    double val = (A[i] > LBA_A_ASYMPTOTIC)
      ? plba_core(t[i], A[i], b[i], v[i], sv[i], denom)
        : plba_noA (t[i],       b[i], v[i], sv[i], denom);
    if      (val < 0.0) val = 0.0;
    else if (val > 1.0) val = 1.0;
    cdf[i] = val;
  }
  return cdf;
}



#endif


