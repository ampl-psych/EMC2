#ifndef wald_functions_h
#define wald_functions_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "pnorm_utils.h"

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Numerical stability clamps — applied in the fast path only.
// Guarantees a >= A_EPS and |l| >= L_EPS, eliminating degenerate branches
// in the core functions. The scalar wrappers use the asymptotic fallbacks
// instead, preserving numerical accuracy for small a.
// ---------------------------------------------------------------------------

constexpr double A_EPS = 1e-4;
constexpr double L_EPS = 1e-4;
constexpr double K_MAX = 1e6;  // upper limit for threshold/noise ratio.

inline void clamp_l(double& l)
{
  l = (l > -L_EPS && l < L_EPS) ? (l >= 0.0 ? L_EPS : -L_EPS) : l;
}

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

// const double L_PI = 1.1447298858494001741434;  // std::log(M_PI)
//
// double pigt0(double t, double k = 1., double l = 1.){
//   //if (t <= 0.){
//   //  return 0.;
//   //}
//   double mu = k / l;
//   double lambda = k * k;
//
//   double p1 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. + t/mu), 0., 1., true, false);
//   double p2 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. - t/mu), 0., 1., true, false);
//
//   return std::exp(std::exp(std::log(2. * lambda) - std::log(mu)) + std::log(p1)) + p2;
// }
//
// double digt0(double t, double k = 1., double l = 1.){
//   //if (t <= 0.) {
//   //  return 0.;
//   //}
//   double lambda = k * k;
//   double e;
//   if (l == 0.) {
//     e = -.5 * lambda / t;
//   } else {
//     double mu = k / l;
//     e = - (lambda / (2. * t)) * ((t * t) / (mu * mu) - 2. * t / mu + 1.);
//   }
//   return std::exp(e + .5 * std::log(lambda) - .5 * std::log(2. * t * t * t * M_PI));
// }
//
// double pigt(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10){
//   if (t <= 0.){
//     return 0.;
//   }
//   if (a < threshold){
//     return pigt0(t, k, l);
//   }
//
//   double sqt = std::sqrt(t);
//   double lgt = std::log(t);
//   double cdf;
//
//   if (l < threshold){
//     double t5a = 2. * R::pnorm((k + a) / sqt, 0., 1., true, false) - 1;
//     double t5b = 2. * R::pnorm((- k - a) / sqt, 0., 1., true, false) - 1;
//
//     double t6a = - .5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
//     double t6b = - .5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);
//
//     cdf = 1. + std::exp(t6a) - std::exp(t6b) + ((- k + a) * t5a - (k - a) * t5b) / (2. * a);
//   } else {
//     double t1a = std::exp(- .5 * std::pow(k - a - t * l, 2) / t);
//     double t1b = std::exp(- .5 * std::pow(a + k - t * l, 2) / t);
//     double t1 = std::exp(.5* (lgt - M_LN2 - L_PI)) * (t1a - t1b);
//
//     double t2a = std::exp(2. * l * (k - a) + R::pnorm(- (k - a + t * l) / sqt, 0., 1., true, true));
//     double t2b = std::exp(2. * l * (k + a) + R::pnorm(- (k + a + t * l) / sqt, 0., 1., true, true));
//     double t2 = a + (t2b - t2a) / (2. * l);
//
//     double t4a = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
//     double t4b = 2. * R::pnorm((k - a) / sqt - sqt * l, 0., 1., true, false) - 1.;
//     double t4 = .5 * (t * l - a - k + .5 / l) * t4a + .5 * (k - a - t * l - .5 / l) * t4b;
//
//     cdf = .5 * (t4 + t2 + t1) / a;
//   }
//   if (cdf < 0. || std::isnan(cdf)) {
//     return 0.;
//   }
//   return cdf;
// }
//
// double digt(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10){
//   if (t <= 0.){
//     return 0.;
//   }
//   if (a < threshold){
//     return digt0(t, k, l);
//   }
//   double pdf;
//   if (l < threshold){
//     double term = std::exp(- (k - a) * (k - a) / (2. * t)) - std::exp(- (k + a) * (k + a) / (2. * t));
//     pdf = std::exp(-.5 * (M_LN2 + L_PI + std::log(t)) + std::log(term) - M_LN2 - std::log(a));
//   } else {
//     double sqt = std::sqrt(t);
//
//     double t1a = - std::pow(a - k + t * l, 2) / (2. * t);
//     double t1b = - std::pow(a + k - t * l, 2) / (2. * t);
//     double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);
//
//     double t2a = 2. * R::pnorm((- k + a) / sqt + sqt * l, 0., 1., true, false) - 1.;
//     double t2b = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
//     double t2 = std::exp(std::log(.5) + std::log(l)) * (t2a + t2b);
//
//     pdf = std::exp(std::log(t1 + t2) - M_LN2 - std::log(a));
//   }
//   if (pdf < 0. || std::isnan(pdf)) {
//     return 0.;
//   }
//   return pdf;
// }

#endif
