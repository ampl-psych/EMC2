#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "RaceSpec.h"
#include "math_utils.h"  // must be before Rcpp
#include "pnorm_utils.h"
#include "ParamTable.h"
// #include "CensorSpec.h"

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
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0);

NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0);

// ---------------------------------------------------------------------------
// Fast ParamTable-based functions
// ---------------------------------------------------------------------------
// void drdm_fast(const NumericVector& rts,
//                const ParamTable& pt,
//                const RaceSpec& spec,
//                const std::vector<int>& idx,
//                double* __restrict__ ll_row);

// void prdm_fast(const NumericVector& rts,
//                const ParamTable& pt,
//                const RaceSpec& spec,
//                const std::vector<int>& idx,
//                double* __restrict__ ll_row);

void drdm_prdm_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RaceSpec& spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__ ll_row,
                    RaceScratch& scratch);

// void rdm_censor(const CensorSpec&    censor,
//                 const ParamTable&    pt,
//                 const RaceSpec&      spec,
//                 double* __restrict__ ll_row,
//                 RaceScratch&         scratch);

void rdm_truncate(const std::vector<int>& idx,
                  const std::vector<double>& bound,
                  const ParamTable& pt,
                  const RaceSpec& spec,
                  double* __restrict__ out,
                  RaceScratch& scratch);

#endif // rdm_h
