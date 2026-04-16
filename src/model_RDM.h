#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "ParamTable.h"
#include "RaceSpec.h"
#include "pnorm_utils.h"

using namespace Rcpp;

const double L_PI = 1.1447298858494001741434;  // std::log(M_PI)


// RDM
inline double pigt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.){
  //  return 0.;
  //}
  double mu = k / l;
  double lambda = k * k;

  double z1 = std::sqrt(lambda / t) * (1.0 + t / mu);
  double z2 = std::sqrt(lambda / t) * (1.0 - t / mu);

  double p1 = PNORM_STD(z1, /*lower=*/false, /*logp=*/false);
  double p2 = PNORM_STD(z2, /*lower=*/false, /*logp=*/false);

  return std::exp(std::exp(std::log(2. * lambda) - std::log(mu)) + std::log(p1)) + p2;
}

inline double digt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.) {
  //  return 0.;
  //}
  double lambda = k * k;
  double e;
  if (l == 0.) {
    e = -.5 * lambda / t;
  } else {
    double mu = k / l;
    e = - (lambda / (2. * t)) * ((t * t) / (mu * mu) - 2. * t / mu + 1.);
  }
  return std::exp(e + .5 * std::log(lambda) - .5 * std::log(2. * t * t * t * M_PI));
}

double pigt(double t, double k = 1., double l = 1., double a = .1, double threshold = 1e-10)
{
  if (t <= 0.0) {
    return 0.0;
  }
  if (a < threshold) {
    return pigt0(t, k, l);   // can be dropped if a is bounded away from 0
  }

  const double sqt     = std::sqrt(t);
  const double inv_sqt = 1.0 / sqt;
  const double lgt     = std::log(t);
  double cdf;

  if (std::fabs(l) < threshold) {
    // keep your existing l≈0 branch for now
    double t5a = 2.0 * R::pnorm((k + a) / sqt, 0.0, 1.0, true, false) - 1.0;
    double t5b = 2.0 * R::pnorm((-k - a) / sqt, 0.0, 1.0, true, false) - 1.0;

    double t6a = -0.5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
    double t6b = -0.5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);

    cdf = 1.0 + std::exp(t6a) - std::exp(t6b) + ((-k + a) * t5a - (k - a) * t5b) / (2.0 * a);
  } else {
    const double inv_t = 1.0 / t;

    // t1: simplified as sqrt(t)/(sqrt(2π)) * (exp(...) - exp(...))
    const double tmp1 = k - a - t * l;
    const double tmp2 = a + k - t * l;
    const double t1a  = std::exp(-0.5 * tmp1 * tmp1 * inv_t);
    const double t1b  = std::exp(-0.5 * tmp2 * tmp2 * inv_t);
    const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
    const double t1 = sqt * inv_sqrt_2pi * (t1a - t1b);

    // t2: unchanged algebraically
    const double argA = -(k - a + t * l) * inv_sqt;
    const double argB = -(k + a + t * l) * inv_sqt;
    const double t2a = std::exp(2.0 * l * (k - a) +
                                R::pnorm(argA, 0.0, 1.0, true, true));
    const double t2b = std::exp(2.0 * l * (k + a) +
                                R::pnorm(argB, 0.0, 1.0, true, true));
    const double t2 = a + (t2b - t2a) / (2.0 * l);

    // t4: unchanged
    const double t4a = 2.0 * R::pnorm((k + a) * inv_sqt - sqt * l,
                                      0.0, 1.0, true, false) - 1.0;
    const double t4b = 2.0 * R::pnorm((k - a) * inv_sqt - sqt * l,
                                      0.0, 1.0, true, false) - 1.0;
    const double t4 = 0.5 * (t * l - a - k + 0.5 / l) * t4a
    + 0.5 * (k - a - t * l - 0.5 / l) * t4b;

    cdf = 0.5 * (t4 + t2 + t1) / a;
  }

  if (cdf < 0.0 || std::isnan(cdf)) {
    return 0.0;
  }
  if (cdf > 1.0) {
    return 1.0;
  }
  return cdf;
}


double digt(double t, double k = 1., double l = 1., double a = .1, double threshold = 1e-10)
{
  if (t <= 0.) {
    return 0.0;  // you already ensure t>0 in the caller; this is just extra safety
  }

  // Small-a asymptotic: if your bounds guarantee a >> threshold, you can drop this.
  if (a < threshold) {
    return digt0(t, k, l);
  }

  double pdf;

  if (std::fabs(l) < threshold) {
    // small-l branch; again, can be dropped if |l| is bounded away from 0
    const double inv_two_t = 0.5 / t;
    const double km_a = k - a;
    const double kp_a = k + a;
    const double term =
      std::exp(-km_a * km_a * inv_two_t) -
      std::exp(-kp_a * kp_a * inv_two_t);

    // same expression as before, just with minor micro-optimisations
    pdf = std::exp(
      -0.5 * (M_LN2 + L_PI + std::log(t))  // normalising bits
      + std::log(term)                     // log(term)
      - M_LN2
      - std::log(a)
    );
  } else {
    // main branch
    const double sqt    = std::sqrt(t);
    const double inv_sqt = 1.0 / sqt;
    const double inv_t   = 1.0 / t;

    // t1: simplified and avoiding pow
    const double temp1 = a - k + t * l;
    const double temp2 = a + k - t * l;
    const double t1a   = -0.5 * temp1 * temp1 * inv_t;
    const double t1b   = -0.5 * temp2 * temp2 * inv_t;

    // 1/sqrt(2π)
    const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
    const double t1 = inv_sqrt_2pi * (std::exp(t1a) - std::exp(t1b)) * inv_sqt;

    // t2: no log/exp, just 0.5 * l
    const double arg1 = (-k + a) * inv_sqt + sqt * l;
    const double arg2 = ( k + a) * inv_sqt - sqt * l;
    const double t2a = 2.0 * R::pnorm(arg1, 0.0, 1.0, true, false) - 1.0;
    const double t2b = 2.0 * R::pnorm(arg2, 0.0, 1.0, true, false) - 1.0;
    const double t2  = 0.5 * l * (t2a + t2b);

    const double sum = t1 + t2;

    // pdf = exp(log(t1 + t2) - M_LN2 - log(a)) == (t1 + t2) / (2a)
    if (sum <= 0.0 || !std::isfinite(sum)) {
      pdf = 0.0;
    } else {
      pdf = sum / (2.0 * a);
    }
  }

  if (pdf < 0.0 || std::isnan(pdf)) {
    return 0.0;
  }
  return pdf;
}

NumericVector drdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = digt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

NumericVector prdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = pigt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}


// [[Rcpp::export]]
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt(t[i], B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt(t[i], B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return cdf;
}

#endif



// RDM model
// Minimal values of A and v/s to avoid numerical instability.
constexpr double A_EPS = 1e-4;   // or 1e-3
constexpr double L_EPS = 1e-4;   // or 1e-3

inline void clamp_a_l(double& a, double& l) {
  if (a < A_EPS) a = A_EPS;
  if (l > -L_EPS && l < L_EPS) l = (l >= 0.0 ? L_EPS : -L_EPS);
}

// Core digt and pigt
inline double digt_core(double t, double k, double l, double a)
{
  if (t <= 0.0) {
    return 0.0;
  }

  const double sqt      = std::sqrt(t);
  const double inv_sqt  = 1.0 / sqt;
  const double inv_t    = 1.0 / t;
  const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);

  // t1 part
  const double temp1 = a - k + t * l;
  const double temp2 = a + k - t * l;
  const double t1a   = -0.5 * temp1 * temp1 * inv_t;
  const double t1b   = -0.5 * temp2 * temp2 * inv_t;
  const double t1    = inv_sqrt_2pi * (std::exp(t1a) - std::exp(t1b)) * inv_sqt;

  // t2 part
  const double arg1 = (-k + a) * inv_sqt + sqt * l;
  const double arg2 = ( k + a) * inv_sqt - sqt * l;

  const double t2a = 2.0 * PNORM_STD(arg1, /*lower=*/true, /*logp=*/false) - 1.0;
  const double t2b = 2.0 * PNORM_STD(arg2, /*lower=*/true, /*logp=*/false) - 1.0;
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

// pigt core
inline double pigt_core(double t, double k, double l, double a)
{
  if (t <= 0.0) {
    return 0.0;
  }

  const double sqt      = std::sqrt(t);
  const double inv_sqt  = 1.0 / sqt;
  const double inv_t    = 1.0 / t;
  const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);

  // t1 term: sqt / sqrt(2π) * (exp(...) - exp(...))
  const double tmp1 = k - a - t * l;
  const double tmp2 = a + k - t * l;
  const double t1a  = std::exp(-0.5 * tmp1 * tmp1 * inv_t);
  const double t1b  = std::exp(-0.5 * tmp2 * tmp2 * inv_t);
  const double t1   = sqt * inv_sqrt_2pi * (t1a - t1b);

  // t2 term
  const double argA = -(k - a + t * l) * inv_sqt;
  const double argB = -(k + a + t * l) * inv_sqt;

  const double t2a = std::exp(2.0 * l * (k - a) +
                              PNORM_STD(argA, /*lower=*/true, /*logp=*/true));
  const double t2b = std::exp(2.0 * l * (k + a) +
                              PNORM_STD(argB, /*lower=*/true, /*logp=*/true));
  const double t2  = a + (t2b - t2a) / (2.0 * l);

  // t4 term
  const double t4a = 2.0 * PNORM_STD((k + a) * inv_sqt - sqt * l,
                                      /*lower=*/true, /*logp=*/false) - 1.0;
  const double t4b = 2.0 * PNORM_STD((k - a) * inv_sqt - sqt * l,
                                      /*lower=*/true, /*logp=*/false) - 1.0;
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

// pdf for RDM, winners only (contiguous + mask)
void drdm_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* raw /* length n_trials */)
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
    if (!win_ptr[i])
      continue;  // only winners get pdf; others left unchanged

    double t_eff = rt[i] - t0[i];

    double inv_s = 1.0 / s[i];
    double a     = 0.5 * A[i] * inv_s;
    double l     = v[i]   * inv_s;
    double k     = B[i]   * inv_s + a;

    clamp_a_l(a, l);

    double pdf = digt_core(t_eff, k, l, a);   // handles t<=0 internally

    raw[i] = pdf;  // may be 0, NaN, etc.; interpreted later
  }
}

// cdf for RDM, losers only (contiguous + mask)
void prdm_fast(const NumericVector& rts,
               const ParamTable& pt,
               const RaceSpec& spec,
               const LogicalVector& winner,
               double* raw /* length n_trials */)
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
    if (win_ptr[i])
      continue;  // only losers get cdf

    double t_eff = rt[i] - t0[i];

    double inv_s = 1.0 / s[i];
    double a     = 0.5 * A[i] * inv_s;
    double l     = v[i]   * inv_s;
    double k     = B[i]   * inv_s + a;

    clamp_a_l(a, l);

    double cdf = pigt_core(t_eff, k, l, a);   // handles t<=0 internally

    raw[i] = cdf;  // may be NaN, 0, 1, etc.; interpreted later
  }
}
