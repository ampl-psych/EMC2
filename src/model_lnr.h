#ifndef lnr_h
#define lnr_h

#include <Rcpp.h>
#include "utility_functions.h"
#include "ParamTable.h"
#include "RaceSpec.h"
#include "pnorm_utils.h"

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

#endif
