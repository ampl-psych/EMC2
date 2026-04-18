#ifndef lnr_h
#define lnr_h

#include <Rcpp.h>
#include "RaceSpec.h"
#include "utility_functions.h"
#include "math_utils.h"  // must be before Rcpp
#include "pnorm_utils.h"
#include "ParamTable.h"

using namespace Rcpp;

NumericVector plnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  // 0 = m, 1 = s, 2 = t0
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rts[i] - pars(i,2) > 0) && (is_ok[i] == TRUE)){
        out[k] = PLNORM(rts[i] - pars(i,2), pars(i,0), pars(i,1));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return out;
}

NumericVector dlnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  // 0 = m, 1 = s, 2 = t0
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rts[i] - pars(i,2) > 0) && (is_ok[i] == TRUE)){
        out[k] = DLNORM(rts[i] - pars(i,2), pars(i,0), pars(i,1));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return out;
}

void dlnr_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* raw)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* m  = &pt.base(0, spec.col_m);
  const double* s  = &pt.base(0, spec.col_s);
  const double* t0 = &pt.base(0, spec.col_t0);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (!win_ptr[i]) continue;

    if (std::isnan(m[i])) { raw[i] = 0.0; continue; }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0)    { raw[i] = 0.0; continue; }

    double pdf = DLNORM(t_eff, m[i], s[i]);
    raw[i] = (std::isfinite(pdf) && pdf >= 0.0) ? pdf : 0.0;
  }
}

void plnr_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* raw)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* m  = &pt.base(0, spec.col_m);
  const double* s  = &pt.base(0, spec.col_s);
  const double* t0 = &pt.base(0, spec.col_t0);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (win_ptr[i]) continue;

    if (std::isnan(m[i])) { raw[i] = 0.0; continue; }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0)    { raw[i] = 0.0; continue; }

    double cdf = PLNORM(t_eff, m[i], s[i]);
    if      (!std::isfinite(cdf) || cdf < 0.0) cdf = 0.0;
    else if (cdf > 1.0)                         cdf = 1.0;
    raw[i] = cdf;
  }
}


// Hot path: gather → compute → scatter
// Uses scratch.v for m, scratch.s for s — B, A, sv left unused.
void dlnr_plnr_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ raw,
                    RaceScratch& scratch)
{
  const double* __restrict__ rt = rts.begin();
  const double* __restrict__ m  = &pt.base(0, spec.col_m);
  const double* __restrict__ s  = &pt.base(0, spec.col_s);
  const double* __restrict__ t0 = &pt.base(0, spec.col_t0);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // --- Winners: gather (m → scratch.v, s → scratch.s) ---
  for (int j = 0; j < n_win; ++j) {
    const int i      = idx_win[j];
    scratch.t_eff[j] = rt[i] - t0[i];
    scratch.v[j]     = m[i];
    scratch.s[j]     = s[i];
  }

  // --- Winners: compute (contiguous — vectorisable) ---
  for (int j = 0; j < n_win; ++j) {
    if (std::isnan(scratch.v[j]) || scratch.t_eff[j] <= 0.0) {
      scratch.out[j] = 0.0;
      continue;
    }
    double val = DLNORM(scratch.t_eff[j], scratch.v[j], scratch.s[j]);
    scratch.out[j] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;
  }

  // --- Winners: scatter ---
  for (int j = 0; j < n_win; ++j) raw[idx_win[j]] = scratch.out[j];

  // --- Losers: gather ---
  for (int j = 0; j < n_los; ++j) {
    const int i      = idx_los[j];
    scratch.t_eff[j] = rt[i] - t0[i];
    scratch.v[j]     = m[i];
    scratch.s[j]     = s[i];
  }

  // --- Losers: compute (contiguous — vectorisable) ---
  for (int j = 0; j < n_los; ++j) {
    if (std::isnan(scratch.v[j]) || scratch.t_eff[j] <= 0.0) {
      scratch.out[j] = 0.0;
      continue;
    }
    double val = PLNORM(scratch.t_eff[j], scratch.v[j], scratch.s[j]);
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    scratch.out[j] = val;
  }

  // --- Losers: scatter ---
  // SURVIVOR! 1-CDF, not CDF
  for (int j = 0; j < n_los; ++j) raw[idx_los[j]] = 1-scratch.out[j];
}

#endif
