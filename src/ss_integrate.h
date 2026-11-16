#ifndef ss_integrate_h
#define ss_integrate_h

// ----------------------------------------------------------------------------
// Stop-signal integration and prepared stop distributions.
//
// All stop-signal integrals (stop-success probability, lower-censoring
// response mass) go through the package's own adaptive cubature
// (`hcubature` in hcubature.h / hcubature.cpp, as used by the DDM and the
// ex-Gaussian race). Integrands use hcubature's callback signature:
//   int f(unsigned dim, const double* x, void* params, unsigned fdim, double* out)
// evaluating at x[0] and writing the value to out[0].
//
// hcubature has no infinite-range support, so callers integrate over a finite
// window. For the stop process (ex-Gaussian with mean muS and scale
// sigmaS/tauS) the window is
//   [max(lb, muS - K_SIGMA*sigmaS), min(upper, muS + K_SIGMA*sigmaS + K_TAU*tauS)]
// which contains all but a negligible part of the mass; an infinite lower
// bound (exgS_lb = -Inf, the "no truncation" exception) is thereby clipped.
// ----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include "pnorm_utils.h"
#include "hcubature.h"   // the package's adaptive cubature

constexpr double SS_WINDOW_K_SIGMA = 8.0;
constexpr double SS_WINDOW_K_TAU   = 16.0;

using ss_integrand_fn = int (*)(unsigned, const double*, void*, unsigned, double*);

// Integrate f over [lo, hi]. Returns 0 for an empty/invalid window.
inline double ss_integrate(ss_integrand_fn f, void* params,
                           double lo, double hi,
                           double abs_tol, double rel_tol,
                           std::size_t max_eval) {
  if (!(hi > lo) || !is_finite(lo) || !is_finite(hi)) return 0.0;
  double a = lo, b = hi, val = 0.0, err = 0.0;
  hcubature(f, params, 1, &a, &b, max_eval, abs_tol, rel_tol, &val, &err);
  return val;
}

// Finite stop-process window (see header comment). `upper` is in stop-relative
// time (i.e. time since SSD); pass R_PosInf when unrestricted.
inline double ss_stop_window_lo(double lb, double muS, double sigS,
                                double k_sigma = SS_WINDOW_K_SIGMA) {
  double lo = muS - k_sigma * sigS;
  return (is_finite(lb) && lb > lo) ? lb : lo;
}
inline double ss_stop_window_hi(double upper, double muS, double sigS, double tauS,
                                double k_sigma = SS_WINDOW_K_SIGMA,
                                double k_tau = SS_WINDOW_K_TAU) {
  double hi = muS + k_sigma * sigS + k_tau * tauS;
  return (is_finite(upper) && upper < hi) ? upper : hi;
}

enum class SSStopFamily { ExGaussian, Normal, Lognormal, Weibull };

// Prepared non-exGaussian stop distributions. Integration uses standardized
// normal/lognormal time or Weibull cumulative hazard, including the Jacobian.
// ExGaussian stops use the existing ExGaussian class in exgaussian_functions.h.
struct SSStop {
  SSStopFamily family;
  double location, scale, lower, inv_scale, log_scale, alpha, log_tail_factor;

  // (a, b) = (mu, sigma), (meanlog, sdlog), or (shape, scale).
  // lb is the normal lower truncation bound; alpha is its standardized value.
  SSStop(SSStopFamily f, double a, double b, double lb=0.0)
    : family(f), location(a), scale(b), lower(lb), inv_scale(1.0/b),
      log_scale(std::log(b)), alpha((lb-a)/b), log_tail_factor(0.0) {
    if (family==SSStopFamily::Normal)
      log_tail_factor=alpha>0 ? std::log(pnorm_detail::normal_tail_factor(alpha))
                             : PNORM_STD(alpha,false,true);
  }

  // log P(T_stop > t), conditional on lower truncation for normal stops.
  double log_survivor(double t) const {
    if (family==SSStopFamily::Normal) {
      if (t<=lower) return 0.0;
      if (is_pos_inf(t)) return neg_inf();
      if (alpha>0) {
        // Use scaled normal tails to avoid underflow in the truncation ratio.
        const double d=(t-lower)*inv_scale;
        return -d*(alpha+0.5*d) +
          std::log(pnorm_detail::normal_tail_factor(alpha+d))-log_tail_factor;
      }
      return PNORM_STD((t-location)*inv_scale,false,true)-log_tail_factor;
    }
    if (t<=0) return 0.0;
    if (family==SSStopFamily::Lognormal)
      return PNORM_STD((std::log(t)-location)*inv_scale,false,true);
    return -std::exp(location*(std::log(t)-log_scale));
  }

  // Map the integration coordinate to stop-relative time.
  // For alpha > 0, normal time is scaled by the tail width sigma/(1+alpha).
  double time(double x) const {
    if (family==SSStopFamily::Weibull)
      return x<=0 ? 0.0 : std::exp(log_scale+std::log(x)/location);
    if (family==SSStopFamily::Lognormal) return std::exp(location+scale*x);
    return alpha>0 ? lower+scale*(x/(1.0+alpha)) : location+scale*x;
  }

  double coordinate(double t) const {
    if (family==SSStopFamily::Weibull)
      return t<=0 ? 0.0 : std::exp(location*(std::log(t)-log_scale));
    if (family==SSStopFamily::Lognormal)
      return t<=0 ? neg_inf() : (std::log(t)-location)*inv_scale;
    return alpha>0 ? ((t-lower)*inv_scale)*(1.0+alpha) : (t-location)*inv_scale;
  }

  // Stop density in integration coordinates: log(f_stop(t(x)) * dt/dx).
  double log_weight(double x) const {
    if (family==SSStopFamily::Weibull) return -x;
    if (family==SSStopFamily::Normal && alpha>0) {
      const double d=x/(1.0+alpha);
      return -d*(alpha+0.5*d)-pnorm_detail::LOG_SQRT_2PI-
        log_tail_factor-std::log1p(alpha);
    }
    return -0.5*x*x-pnorm_detail::LOG_SQRT_2PI-
      (family==SSStopFamily::Normal ? log_tail_factor : 0.0);
  }

  void bounds(double upper, double& lo, double& hi) const {
    // Omitted stop mass is below 3e-20 (10 normal SDs or exp(-45)).
    // Positive normal truncation is scaled by its tail width, sigma/(1+alpha).
    lo=family==SSStopFamily::Weibull ? 0.0 :
       family==SSStopFamily::Normal ? (alpha>0 ? 0.0 : std::max(alpha,-10.0)) : -10.0;
    hi=family==SSStopFamily::Weibull ? 45.0 :
       family==SSStopFamily::Normal && alpha>0 ? 48.0 : 10.0;
    hi=std::min(hi,coordinate(upper));
  }
};

#endif
