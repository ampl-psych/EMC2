// pnorm_utils.h
// Unified pnorm/plnorm/dlnorm implementation for all race models.
//
// Compile-time precision modes via PNORM_MODE (set in src/Makevars):
//
//   PNORM_MODE 0  — R::pnorm, exact double precision
//                   Requires Rcpp / R API. NOT thread-safe.
//   PNORM_MODE 1  (default) — Hart (1968) rational approximation,
//                   near-double precision (~1.2e-15 max abs error).
//                   Two-branch: rational poly (|x| < 7.07) /
//                   continued fraction (|x| >= 7.07).
//                   Pure C++; thread-safe.
//   PNORM_MODE 2  — A&S 7.1.26 rational approximation,
//                   single-precision quality (~1.5e-7 max abs error).
//                   Branchless polynomial; best vectorisability.
//                   Pure C++; thread-safe.
//
// Example Makevars lines:
//   PKG_CXXFLAGS = -DPNORM_MODE=0   # exact R::pnorm (not thread-safe)
//   PKG_CXXFLAGS = -DPNORM_MODE=1   # fast, near-double (default)
//   PKG_CXXFLAGS = -DPNORM_MODE=2   # fastest, single-precision quality
//
// Public macros (all modes):
//   PNORM_STD(x, lower, logp)   — standard normal CDF
//   PLNORM(t, m, s)             — log-normal CDF (lower tail, no log)
//   DLNORM(t, m, s)             — log-normal PDF (always exact formula)

#ifndef PNORM_UTILS_H
#define PNORM_UTILS_H

#include <cmath>

#ifndef PNORM_MODE
#define PNORM_MODE 1  // default: Hart near-double, pure C++, thread-safe
#endif

// Rcpp is only needed for mode 0 (R::pnorm). Modes 1 and 2 are pure C++
// and safe to call from worker threads.
#if PNORM_MODE == 0
#include <Rcpp.h>
#endif

// Large finite sentinel used in place of -inf for logp underflow in modes 1
// and 2. Safe under -ffinite-math-only / -Ofast. Any log-likelihood floor
// will treat this identically to -inf in practice.
#define PNORM_LOG_ZERO (-1e30)

namespace pnorm_detail {

// ===========================================================================
// MODE 0 — exact R::pnorm
// Requires R API. NOT safe to call from non-main threads.
// ===========================================================================

#if PNORM_MODE == 0
inline double exact_pnorm_std(double x, bool lower, bool logp)
{
  return R::pnorm(x, 0.0, 1.0, lower, logp);
}
#endif

// ===========================================================================
// MODE 1 — Hart (1968) rational approximation
// Source: https://stackoverflow.com/a/23119456 (CC BY-SA 3.0)
// Max absolute error: ~1.2e-15 over the full real line.
// Two-branch: rational poly for |x| < SPLIT, continued fraction otherwise.
// Pure C++; thread-safe.
// ===========================================================================

#if PNORM_MODE == 1
namespace hart {

constexpr double RT2PI = 2.506628274631000502415765284811;
constexpr double SPLIT = 7.07106781186547;

constexpr double N0 = 220.206867912376;
constexpr double N1 = 221.213596169931;
constexpr double N2 = 112.079291497871;
constexpr double N3 = 33.912866078383;
constexpr double N4 = 6.37396220353165;
constexpr double N5 = 0.700383064443688;
constexpr double N6 = 3.52624965998911e-02;

constexpr double M0 = 440.413735824752;
constexpr double M1 = 793.826512519948;
constexpr double M2 = 637.333633378831;
constexpr double M3 = 296.564248779674;
constexpr double M4 = 86.7807322029461;
constexpr double M5 = 16.064177579207;
constexpr double M6 = 1.75566716318264;
constexpr double M7 = 8.83883476483184e-02;

inline double phi(double x)
{
  const double z = std::fabs(x);
  double c = 0.0;

  if (z <= 37.0) {
    const double e = std::exp(-z * z / 2.0);
    if (z < SPLIT) {
      const double n = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
      const double d = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
      c = e * n / d;
    } else {
      const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
      c = e / (RT2PI * f);
    }
  }
  return x <= 0.0 ? c : 1.0 - c;
}

} // namespace hart

inline double hart_pnorm_std(double x, bool lower, bool logp)
{
  double cdf = hart::phi(x);
  if (!lower) cdf = 1.0 - cdf;
  if (logp) return (cdf > 0.0) ? std::log(cdf) : PNORM_LOG_ZERO;
  return cdf;
}
#endif // PNORM_MODE == 1

// ===========================================================================
// MODE 2 — A&S 7.1.26 rational approximation
// Max absolute error: ~1.5e-7 over the full real line.
// Branchless polynomial path; best auto-vectorisability.
// Pure C++; thread-safe.
// ===========================================================================

#if PNORM_MODE == 2
namespace as7126 {

inline double fast_erf(double x)
{
  constexpr double p  = 0.3275911;
  constexpr double a1 =  0.254829592;
  constexpr double a2 = -0.284496736;
  constexpr double a3 =  1.421413741;
  constexpr double a4 = -1.453152027;
  constexpr double a5 =  1.061405429;

  const double ax   = std::fabs(x);
  const double t    = 1.0 / (1.0 + p * ax);
  const double poly = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
  const double e    = std::exp(-ax * ax);
  return std::copysign(1.0 - poly * e, x);
}

inline double phi(double x)
{
  constexpr double inv_sqrt2 = 0.70710678118654752440;
  return 0.5 * (1.0 + fast_erf(x * inv_sqrt2));
}

} // namespace as7126

inline double as_pnorm_std(double x, bool lower, bool logp)
{
  double cdf = as7126::phi(x);
  if (!lower) cdf = 1.0 - cdf;
  if (logp) return (cdf > 0.0) ? std::log(cdf) : PNORM_LOG_ZERO;
  return cdf;
}
#endif // PNORM_MODE == 2

// ===========================================================================
// Shared helpers
// ===========================================================================

inline double fast_dlnorm(double t, double m, double s)
{
  const double log_t = std::log(t);
  const double z     = (log_t - m) / s;
  return std::exp(-0.5 * z * z) / (t * s * std::sqrt(2.0 * M_PI));
}

} // namespace pnorm_detail

// ===========================================================================
// PNORM_STD(x, lower, logp)
//
// All four (lower, logp) combinations used in the codebase are handled.
// The compiler folds dead branches away when lower/logp are literal bools.
//
//   (true,  false) — Φ(x)            — digt_core t2, pigt_core t4
//   (true,  true)  — log Φ(x)        — pigt_core t2
//   (false, false) — 1 − Φ(x)        — pigt0
//   (false, true)  — log(1 − Φ(x))   — safe fallback
// ===========================================================================

#if   PNORM_MODE == 2
#define PNORM_STD(x, lower, logp) pnorm_detail::as_pnorm_std((x), (lower), (logp))
#elif PNORM_MODE == 1
#define PNORM_STD(x, lower, logp) pnorm_detail::hart_pnorm_std((x), (lower), (logp))
#else
#define PNORM_STD(x, lower, logp) pnorm_detail::exact_pnorm_std((x), (lower), (logp))
#endif

// ===========================================================================
// PLNORM(t, m, s) — log-normal CDF, lower tail, no log
// DLNORM(t, m, s) — log-normal PDF (always exact; no R call needed)
// ===========================================================================

#define PLNORM(t, m, s) PNORM_STD((std::log(t) - (m)) / (s), true, false)
#define DLNORM(t, m, s) pnorm_detail::fast_dlnorm((t), (m), (s))

#endif // PNORM_UTILS_H
