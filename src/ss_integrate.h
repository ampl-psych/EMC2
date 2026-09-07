#ifndef ss_integrate_h
#define ss_integrate_h

// ----------------------------------------------------------------------------
// Stop-signal 1-D integration helper.
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

#include <cmath>
#include <cstddef>
#include <Rcpp.h>    // R_FINITE, R_PosInf
#include "hcubature.h"   // the package's adaptive cubature

constexpr double SS_WINDOW_K_SIGMA = 8.0;
constexpr double SS_WINDOW_K_TAU   = 16.0;

using ss_integrand_fn = int (*)(unsigned, const double*, void*, unsigned, double*);

// Integrate f over [lo, hi]. Returns 0 for an empty/invalid window.
inline double ss_integrate(ss_integrand_fn f, void* params,
                           double lo, double hi,
                           double abs_tol, double rel_tol,
                           std::size_t max_eval) {
  if (!(hi > lo) || !R_FINITE(lo) || !R_FINITE(hi)) return 0.0;
  double a = lo, b = hi, val = 0.0, err = 0.0;
  hcubature(f, params, 1, &a, &b, max_eval, abs_tol, rel_tol, &val, &err);
  return val;
}

// Finite stop-process window (see header comment). `upper` is in stop-relative
// time (i.e. time since SSD); pass R_PosInf when unrestricted.
inline double ss_stop_window_lo(double lb, double muS, double sigS,
                                double k_sigma = SS_WINDOW_K_SIGMA) {
  double lo = muS - k_sigma * sigS;
  return (R_FINITE(lb) && lb > lo) ? lb : lo;
}
inline double ss_stop_window_hi(double upper, double muS, double sigS, double tauS,
                                double k_sigma = SS_WINDOW_K_SIGMA,
                                double k_tau = SS_WINDOW_K_TAU) {
  double hi = muS + k_sigma * sigS + k_tau * tauS;
  return (R_FINITE(upper) && upper < hi) ? upper : hi;
}

#endif
