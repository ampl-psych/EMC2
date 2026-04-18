#ifndef lba_h
#define lba_h

#include <Rcpp.h>
#include "RaceSpec.h"
#include "utility_functions.h"
#include "math_utils.h"
#include "pnorm_utils.h"
#include "ParamTable.h"

using namespace Rcpp;


// double pnormP(double q, double mean = 0.0, double sd = 1.0,
//               bool lower = true, bool log = false){
//   return R::pnorm(q, mean, sd, lower, log);
// }
//
double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log = false){
  return R::dnorm(x, mean, sd, log);
}

double plba_norm(double t, double A, double b, double v, double sv,
                 bool posdrift = true)
{
  double denom = 1.0;
  if (posdrift) {
    denom = PNORM_STD(v / sv, /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
  }

  double cdf;

  if (A > 1e-10) {
    double zs    = t * sv;
    double cmz   = b - t * v;
    double xx    = cmz - A;
    double cz    = cmz / zs;       // = (cmz - mean) / sd  with mean=0, sd=zs...
    double cz_max = xx / zs;       // standardised arguments for N(0,1)

    cdf = (1.0 + (zs * (dnormP(cz_max) - dnormP(cz))
                    + xx  * PNORM_STD(cz_max, /*lower=*/true, /*logp=*/false)
                    - cmz * PNORM_STD(cz,     /*lower=*/true, /*logp=*/false)) / A) / denom;
  } else {
    // A ~ 0: starting point is fixed at 0, first-passage is just a
    // threshold crossing at b with drift v and noise sv
    // P(T < t) = P(b/t < v + sv*Z) = P(Z > (b/t - v)/sv) = 1 - Φ((b/t - v)/sv)
    cdf = PNORM_STD((b / t - v) / sv, /*lower=*/false, /*logp=*/false) / denom;
  }

  if (cdf < 0.0) return 0.0;
  if (cdf > 1.0) return 1.0;
  return cdf;
}

double dlba_norm(double t, double A, double b, double v, double sv,
                 bool posdrift = true)
{
  double denom = 1.0;
  if (posdrift) {
    denom = PNORM_STD(v / sv, /*lower=*/true, /*logp=*/false);
    if (denom < 1e-10) denom = 1e-10;
  }

  double pdf;

  if (A > 1e-10) {
    double zs    = t * sv;
    double cmz   = b - t * v;
    double cz    = cmz / zs;
    double cz_max = (cmz - A) / zs;

    pdf = (v * (PNORM_STD(cz,     /*lower=*/true, /*logp=*/false)
                  - PNORM_STD(cz_max, /*lower=*/true, /*logp=*/false))
             + sv * (dnormP(cz_max) - dnormP(cz))) / (A * denom);
  } else {
    pdf = dnormP(b / t, v, sv) * b / (t * t * denom);
  }

  if (pdf < 0.0) return 0.0;
  return pdf;
}


NumericVector dlba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = dlba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}

NumericVector plba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = plba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true)

{
  int n = t.size();
  NumericVector pdf(n);

  for (int i = 0; i < n; i++){
    pdf[i] = dlba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift);
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

  for (int i = 0; i < n; i++){
    cdf[i] = plba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift);
  }
  return cdf;
}


void dlba_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* raw)
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
    if (!win_ptr[i])
      continue;  // only winners get pdf

    if (std::isnan(v[i])) {
      raw[i] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0) {
      raw[i] = 0.0;
      continue;
    }

    const double A_i = A[i];
    const double B_i = B[i] + A_i;  // b = B + A, as in dlba_c_pt

    double pdf = dlba_norm(t_eff, A_i, B_i, v[i], sv[i], /*posdrift=*/true);
    if (!std::isfinite(pdf) || pdf < 0.0)
      pdf = 0.0;

    raw[i] = pdf;
  }
}

void plba_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* raw)
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
    if (win_ptr[i])
      continue;  // only losers get cdf

    if (std::isnan(v[i])) {
      raw[i] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0) {
      raw[i] = 0.0;
      continue;
    }

    const double A_i = A[i];
    const double B_i = B[i] + A_i;

    double cdf = plba_norm(t_eff, A_i, B_i, v[i], sv[i], /*posdrift=*/true);
    if (!std::isfinite(cdf) || cdf < 0.0)
      cdf = 0.0;
    else if (cdf > 1.0)
      cdf = 1.0;

    raw[i] = cdf;
  }
}

// Hot path: gather → compute → scatter
void dlba_plba_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ raw,
                    RaceScratch& scratch)
{
  const double* __restrict__ rt  = rts.begin();
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ B   = &pt.base(0, spec.col_B);
  const double* __restrict__ A   = &pt.base(0, spec.col_A);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);

  const int n_win = (int)idx_win.size();
  const int n_los = (int)idx_los.size();

  // --- Winners: gather ---
  for (int j = 0; j < n_win; ++j) {
    const int i      = idx_win[j];
    scratch.t_eff[j] = rt[i] - t0[i];
    scratch.v[j]     = v[i];
    scratch.sv[j]    = sv[i];
    scratch.B[j]     = B[i];
    scratch.A[j]     = A[i];
  }

  // --- Winners: compute (contiguous — vectorisable) ---
  for (int j = 0; j < n_win; ++j) {
    if (std::isnan(scratch.v[j]) || scratch.t_eff[j] <= 0.0) {
      scratch.out[j] = 0.0;
      continue;
    }
    double val = dlba_norm(scratch.t_eff[j],
                           scratch.A[j],
                                    scratch.B[j] + scratch.A[j],
                                                            scratch.v[j], scratch.sv[j], true);
    scratch.out[j] = (std::isfinite(val) && val >= 0.0) ? val : 0.0;
  }

  // --- Winners: scatter ---
  for (int j = 0; j < n_win; ++j) raw[idx_win[j]] = scratch.out[j];

  // --- Losers: gather ---
  for (int j = 0; j < n_los; ++j) {
    const int i      = idx_los[j];
    scratch.t_eff[j] = rt[i] - t0[i];
    scratch.v[j]     = v[i];
    scratch.sv[j]    = sv[i];
    scratch.B[j]     = B[i];
    scratch.A[j]     = A[i];
  }

  // --- Losers: compute (contiguous — vectorisable) ---
  for (int j = 0; j < n_los; ++j) {
    if (std::isnan(scratch.v[j]) || scratch.t_eff[j] <= 0.0) {
      scratch.out[j] = 0.0;
      continue;
    }
    double val = plba_norm(scratch.t_eff[j],
                           scratch.A[j],
                                    scratch.B[j] + scratch.A[j],
                                                            scratch.v[j], scratch.sv[j], true);
    if      (!std::isfinite(val) || val < 0.0) val = 0.0;
    else if (val > 1.0)                         val = 1.0;
    scratch.out[j] = 1-val; // SURVIVAL!
  }

  // --- Losers: scatter ---
  for (int j = 0; j < n_los; ++j) raw[idx_los[j]] = scratch.out[j];
}


#endif


