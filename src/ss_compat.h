#ifndef ss_compat_h
#define ss_compat_h

// ---------------------------------------------------------------------------
// SS-port compatibility shim.
//
// Zach's stop-signal headers (exgaussian_functions.h, ss_exg_analytic.h,
// model_SS_EXG.h, model_SS_RDEX.h) were written against that branch's helper
// layer (composite_functions.h, wald_functions.h, a 700-line
// utility_functions.h). The cens branch has a different, leaner helper layout:
//   * no log1m / log_diff_exp,
//   * no emc2_isnan / emc2_isinf / emc2_isfinite,
//   * pnorm provided as the PNORM_STD(x, lower, logp) macro (pnorm_utils.h),
//     not a pnorm_std() function.
//
// This header provides exactly the small, self-contained helpers the ported SS
// code needs, routing pnorm through cens's own PNORM_STD so the SS models use
// the package's configured normal-CDF implementation. Definitions are copied
// verbatim from Zach's composite_functions.h / utility_functions.h (the generic
// numerics are branch-independent); nothing here pulls in GSL.
// ---------------------------------------------------------------------------

#include <Rcpp.h>
#include <cmath>
#include "pnorm_utils.h"   // PNORM_STD(x, lower, logp)

// --- finiteness predicates (verbatim from Zach utility_functions.h) ---------
// NB: use R_FINITE/ISNAN, not -ffast-math-foldable intrinsics, so the analytic
// overflow/finiteness guards in ss_exg_analytic.h stay live.
static inline bool emc2_isfinite(double x) { return (bool)R_FINITE(x); }
static inline bool emc2_isinf(double x)    { return !R_FINITE(x) && !ISNAN(x); }
static inline bool emc2_isnan(double x)    { return (bool)ISNAN(x); }

// --- pnorm_std: route Zach's call signature onto cens's PNORM_STD macro ------
inline double pnorm_std(double x, bool lower = true, bool log_p = false) {
  return PNORM_STD(x, lower, log_p);
}

// --- log-space helpers (verbatim from Zach composite_functions.h) -----------
inline double log1m(double x) {
  if (ISNAN(x)) return NA_REAL;
  if (x >= 1.0) return R_NegInf;          // log(1 - x) = log(0 or negative)
  if (x == R_NegInf) return 0.0;          // defensive clamp
  return std::log1p(-x);
}

inline double log_diff_exp(double a, double b) {
  if (ISNAN(a) || ISNAN(b)) return NA_REAL;
  if (a == R_PosInf) return NA_REAL;      // +Inf - anything undefined in log space
  if (a < b) return NA_REAL;              // log of negative number is undefined
  if (a == b) return R_NegInf;            // log(exp(a) - exp(a)) = log(0)
  if (b == R_NegInf) return a;

  double diff = b - a;
  if (diff > -0.693147) {                 // diff > -log(2): use expm1 (avoid cancellation)
    return a + std::log(-std::expm1(diff));
  } else {                                // exp(diff) < 0.5: log1m(exp(diff)) is stable
    return a + log1m(std::exp(diff));
  }
}

// log(1 - exp(x)) for x <= 0 (verbatim from Zach composite_functions.h)
inline double log1m_exp(double x) {
  if (ISNAN(x)) return NA_REAL;
  if (x > 0.0) return NA_REAL;            // 1 - exp(x) < 0 for x > 0
  if (x == 0.0) return R_NegInf;          // log(1 - 1) = log(0)
  if (x == R_NegInf) return 0.0;          // log(1 - 0) = 0
  if (x > -0.693147) return std::log(-std::expm1(x));
  return log1m(std::exp(x));
}

// log(exp(a) + exp(b)) (verbatim from Zach composite_functions.h)
inline double log_sum_exp(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  if (a > b) {
    double diff = b - a;
    if (!(diff >= -37.0)) return a;       // exp(-37) < double eps; NaN diffs -> NA
    return a + std::log1p(std::exp(diff));
  } else {
    double diff = a - b;
    if (!(diff >= -37.0)) return b;
    return b + std::log1p(std::exp(diff));
  }
}

// log(theta*exp(lambda1) + (1-theta)*exp(lambda2)) (verbatim from Zach composite_functions.h)
inline double log_mix(double theta, double lambda1, double lambda2) {
  if (ISNAN(theta) || ISNAN(lambda1) || ISNAN(lambda2)) return NA_REAL;
  if (theta < 0.0 || theta > 1.0) return NA_REAL;
  if (theta == 0.0) return lambda2;
  if (theta == 1.0) return lambda1;
  if (lambda1 == R_NegInf && lambda2 == R_NegInf) return R_NegInf;
  if (lambda1 == R_NegInf) return log1m(theta) + lambda2;
  if (lambda2 == R_NegInf) return std::log(theta) + lambda1;
  return log_sum_exp(std::log(theta) + lambda1, log1m(theta) + lambda2);
}

#endif // ss_compat_h
