// pnorm_utils.h
// Unified pnorm/plnorm/dlnorm implementation for all race models.
//
// To use the fast approximation, define USE_FAST_PNORM before including
// this header, or pass -DUSE_FAST_PNORM to the compiler (e.g. in Makevars).
//
// Default: exact R::pnorm / R::plnorm.

#ifndef PNORM_UTILS_H
#define PNORM_UTILS_H

#include <cmath>
#include <Rcpp.h>

// ---------------------------------------------------------------------------
// Implementation detail — keep out of global namespace
// ---------------------------------------------------------------------------

namespace pnorm_detail {

// Fast approximation — rational approximation of Φ(x)
// Source: https://stackoverflow.com/a/23119456 (CC BY-SA 3.0)

// static const double RT2PI = 2.506628274631000502415765284811;
// static const double SPLIT  = 7.07106781186547;

// static const double N0 = 220.206867912376;
// static const double N1 = 221.213596169931;
// static const double N2 = 112.079291497871;
// static const double N3 = 33.912866078383;
// static const double N4 = 6.37396220353165;
// static const double N5 = 0.700383064443688;
// static const double N6 = 3.52624965998911e-02;
//
// static const double M0 = 440.413735824752;
// static const double M1 = 793.826512519948;
// static const double M2 = 637.333633378831;
// static const double M3 = 296.564248779674;
// static const double M4 = 86.7807322029461;
// static const double M5 = 16.064177579207;
// static const double M6 = 1.75566716318264;
// static const double M7 = 8.83883476483184e-02;
//
// inline double phi(double x)
// {
//   const double z = std::fabs(x);
//   double c = 0.0;
//
//   if (z <= 37.0) {
//     const double e = std::exp(-z * z / 2.0);
//     if (z < SPLIT) {
//       const double n = (((((N6*z + N5)*z + N4)*z + N3)*z + N2)*z + N1)*z + N0;
//       const double d = ((((((M7*z + M6)*z + M5)*z + M4)*z + M3)*z + M2)*z + M1)*z + M0;
//       c = e * n / d;
//     } else {
//       const double f = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
//       c = e / (RT2PI * f);
//     }
//   }
//   return x <= 0.0 ? c : 1.0 - c;
// }

inline double phi(double x)
{
  return 0.5 * std::erfc(-x * M_SQRT1_2);
}

inline double fast_pnorm_std(double x, bool lower, bool logp)
{
  double cdf = phi(x);
  if (!lower) cdf = 1.0 - cdf;
  if (logp) {
    if (cdf <= 0.0) return R_NegInf;
    return std::log(cdf);
  }
  return cdf;
}

inline double exact_pnorm_std(double x, bool lower, bool logp)
{
  return R::pnorm(x, 0.0, 1.0, lower, logp);
}

} // namespace pnorm_detail

// ---------------------------------------------------------------------------
// PNORM_STD(x, lower, logp)
// Evaluates the standard normal CDF at x.
// ---------------------------------------------------------------------------

#ifdef USE_FAST_PNORM
#define PNORM_STD(x, lower, logp) pnorm_detail::fast_pnorm_std((x), (lower), (logp))
#else
#define PNORM_STD(x, lower, logp) pnorm_detail::exact_pnorm_std((x), (lower), (logp))
#endif

// ---------------------------------------------------------------------------
// PLNORM(t, m, s) / DLNORM(t, m, s)
// Lognormal CDF and PDF, always lower-tail, non-log scale.
// DLNORM is always the direct formula (no pnorm involved, always fast).
// PLNORM uses PNORM_STD so it respects the USE_FAST_PNORM switch.
// ---------------------------------------------------------------------------

namespace pnorm_detail {

inline double fast_plnorm(double t, double m, double s)
{
  return PNORM_STD((std::log(t) - m) / s, /*lower=*/true, /*logp=*/false);
}

inline double fast_dlnorm(double t, double m, double s)
{
  const double log_t = std::log(t);
  const double z     = (log_t - m) / s;
  return std::exp(-0.5 * z * z) / (t * s * std::sqrt(2.0 * M_PI));
}

} // namespace pnorm_detail

#ifdef USE_FAST_PNORM
#define PLNORM(t, m, s) pnorm_detail::fast_plnorm((t), (m), (s))
#else
#define PLNORM(t, m, s) R::plnorm((t), (m), (s), true, false)
#endif

// DLNORM: always use the direct formula — no #ifdef needed
#define DLNORM(t, m, s) pnorm_detail::fast_dlnorm((t), (m), (s))

#endif // PNORM_UTILS_H
