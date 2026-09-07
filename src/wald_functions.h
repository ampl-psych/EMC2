#ifndef WALD_FUNCTIONS_H
#define WALD_FUNCTIONS_H

// pigt0 / digt0 — Wald (inverse Gaussian) PDF/survival function.
// Thread-safe (via pnorm_utils.h),
// parameterised as threshold k, drift l (mu = k/l, lambda = k*k), matching
// the original hot-loop implementation.
//
// Assumes upstream validation: t > 0 (zero-likelihood short-circuit happens
// before this is called), l already passed through clamp_l() so it's
// bounded away from 0. No NaN/limit/spike guards in the default path -- see
// WALD_RACE_DEFENSIVE at the bottom if you ever need this called on
// unclamped/unvalidated input.

#include <cmath>
#include <Rcpp.h>
#include "nan_check.h"
#include "r_constants.h"
#include "composite_functions.h" // log1p_exp
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
// digt0 — Wald PDF.
[[gnu::always_inline]] inline double digt0(double t, double k, double l)
{
  const double lambda = k * k;
  const double tl_k = t * l / k;
  const double e = -0.5 * (lambda / t) * (tl_k - 1.0) * (tl_k - 1.0);
  return std::exp(e) * std::sqrt(lambda / (2.0 * M_PI * t * t * t));
}

// pigt0 - Wald CDF.
[[gnu::always_inline]] inline double pigt0(double t, double k, double l)
{
  const double lambda = k * k;
  const double mu     = k / l;
  const double sqlt   = std::sqrt(lambda / t);   // shared
  const double tmu    = t / mu;                  // shared
  const double z1     = sqlt * (1.0 + tmu);
  const double z2     = sqlt * (1.0 - tmu);
  return pnorm_upper(z1) * std::exp(2.0 * lambda / mu) + PNORM_STD(z2, false, false);
}


// ---------------------------------------------------------------------------
// pigt0 — Wald survival function P(T > t), overflow-fixed.
//
// Original (naive) formula:
//   S(t) = (1 - Phi(z1)) * exp(2*lambda/mu) + (1 - Phi(z2))
// i.e. pnorm_upper(z1)*exp(2*lambda/mu) + pnorm_upper(z2), with
//   z1 = sqrt(lambda/t)*(1 + t/mu),  z2 = sqrt(lambda/t)*(1 - t/mu).
//
// exp(2*lambda/mu) = exp(2*k*l) overflows double precision once k*l
// exceeds ~354. Fix: keep that term in log-space and combine via
// log-sum-exp, exponentiating only once at the end -- same trick pwald
// uses, but simpler here: both summands in this formula are already
// non-negative (they're upper-tail probabilities, or an upper-tail
// probability scaled by a positive constant), so a plain log1p_exp
// suffices. pwald needs log1m_exp *and* a separate Giner & Smyth (2016)
// asymptotic branch because it expresses the CCDF as 1 - F(t), which
// creates a cancellation problem this formulation never has since it's
// already S(t) directly, not 1-F(t).
//
// Residual caveat: PNORM_STD(..., logp=true) computes log(Phi) as
// std::log(cdf) after computing cdf itself (see pnorm_utils.h), so it
// underflows to the PNORM_LOG_ZERO floor once z gtr ~38, independent of
// this fix. That only bites in the truly extreme tail (z1 gtr ~38, i.e.
// very large k*l together with very small t) and pwald doesn't fully
// solve it either without its own asymptotic branch -- flagging it here
// rather than building that branch, since it's a smaller, second-order
// concern than the overflow this change fixes, and your upstream checks
// already bound the parameter ranges.
// ---------------------------------------------------------------------------
// [[gnu::always_inline]] inline double pigt0(double t, double k, double l)
// {
//   const double lambda = k * k;
//   const double mu     = k / l;
//   const double sqlt   = std::sqrt(lambda / t);
//   const double tmu    = t / mu;
//   const double z1     = sqlt * (1.0 + tmu);
//   const double z2     = sqlt * (1.0 - tmu);
//
//   // log of first summand: log(1 - Phi(z1)) + 2*lambda/mu
//   const double log_a = PNORM_STD(z1, false, true) + 2.0 * lambda / mu;
//   // log of second summand: log(1 - Phi(z2))
//   const double log_b = PNORM_STD(z2, false, true);
//
//   const double m = (log_a > log_b) ? log_a : log_b;
//   return std::exp(m + log1p_exp(-std::fabs(log_a - log_b)));
// }


// ===========================================================================
// WALD_DEFENSIVE — optional, off by default.
//
// Adds guards around pigt0/digt0 (NaN / negative params, x<=0, mu==Inf, x==mu
// spike), at the cost of extra branches per call. Define WALD_DEFENSIVE before
// including this header to get pigt0_safe / digt0_safe alongside the fast
// versions above. Given you already validate t and l upstream, you probably
// don't need this in the hot loop -- it's here mainly for cases outside that
// path (public-facing R function, prior-predictive sweeps, unconstrained
// optimizer excursions) where inputs aren't guaranteed sane.
// ===========================================================================
#ifdef WALD_DEFENSIVE

namespace wald_detail {

// mu = k/l < 0.0 is treated as invalid, same as pwald's bad_params. Note
// lambda = k*k can never be negative for real k, so unlike pwald (which
// takes mu/lambda as free parameters) there's no separate lambda<0 check
// needed here -- it's structurally unreachable.
[[gnu::always_inline]] inline bool bad_params(double t, double k, double l)
{
  if (is_nan(t) || is_nan(k) || is_nan(l)) return true;
  return (k / l) < 0.0;
}

// t<0 / t==0 dropped throughout (guaranteed handled upstream).
[[gnu::always_inline]] inline bool lower_limit(double t, double k, double l)
{
  const double lambda = k * k;
  const double mu = k / l;
  return (t < mu) && is_inf(lambda);
}

// isinf(t) dropped (t guaranteed finite upstream).
[[gnu::always_inline]] inline bool upper_limit(double t, double k, double l)
{
  const double lambda = k * k;
  const double mu = k / l;
  return (t > 0.0 && lambda == 0.0) ||
    (t > mu && (mu == 0.0 || is_inf(lambda)));
}

[[gnu::always_inline]] inline bool is_spike(double t, double k, double l)
{
  const double mu = k / l;
  const double lambda = k * k;
  return (t == mu) && (mu == 0.0 || is_inf(lambda));
}

} // namespace wald_detail

[[gnu::always_inline]] inline double digt0_safe(double t, double k, double l)
{
  if (wald_detail::bad_params(t, k, l)) return na_real();
  if (wald_detail::lower_limit(t, k, l) || wald_detail::upper_limit(t, k, l)) return 0.0;
  if (wald_detail::is_spike(t, k, l)) return pos_inf();
  return digt0(t, k, l);
}

// pigt0 is only ever a survival function P(T>t) (no lower_tail toggle,
// unlike pwald), so the limit/spike return values below follow the
// mathematically consistent convention for S(t) directly, rather than
// literally porting pwald's branches (which return CDF-style values and,
// at the exact spike, return 1.0 unconditionally regardless of tail
// direction -- worth double-checking against your intent if you rely on
// that behaviour elsewhere).
[[gnu::always_inline]] inline double pigt0_safe(double t, double k, double l)
{
  if (wald_detail::bad_params(t, k, l)) return na_real();
  if (wald_detail::lower_limit(t, k, l)) return 1.0;  // below the point mass: S(t)=1
  if (wald_detail::upper_limit(t, k, l)) return 0.0;  // past the point mass: S(t)=0
  if (wald_detail::is_spike(t, k, l))    return 0.0;  // at the point mass itself: S(t)=0
  return pigt0(t, k, l);
}

#endif // WALD_DEFENSIVE

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
  return (is_finite(pdf) && pdf >= 0.0) ? pdf : 0.0;
}

inline double pigt(double t, double k, double l, double a)
{
  if (t <= 0.0) return 0.0;
  if (is_pos_inf(t)) return 1.0;
  if (a < A_EPS) return pigt0(t, k, l);
  clamp_a_l(a, l);
  double cdf = pigt_core(t, k, l, a);
  if (!is_finite(cdf) || cdf < 0.0) return 0.0;
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

#endif // WALD_FUNCTIONS_H
