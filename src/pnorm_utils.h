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
// qnorm (all modes except 0): Wichura AS241 (Appl. Stat. 37:477-484, 1988).
//   Three-region rational approximation; full double precision in a single
//   pass (~1e-15 max absolute error). No refinement step needed.
//   Pure C++; thread-safe.
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
inline double exact_qnorm_std(double p)
{
  return R::qnorm(p, 0.0, 1.0, true, false);
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

[[gnu::always_inline]] inline double phi(double x)
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

[[gnu::always_inline]] inline double fast_erf(double x)
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

[[gnu::always_inline]] inline double phi(double x)
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
// qnorm — Wichura AS241 (Appl. Stat. 37:477-484, 1988)
//
// Three-region rational approximation giving full double precision (~1e-15
// max absolute error) in a single pass. No refinement step needed.
// This is the same algorithm used internally by GSL and R's qnorm().
//
// Domain: p in (0, 1). Returns ±Inf at boundaries.
// ===========================================================================

#if PNORM_MODE != 0
namespace as241 {

// Coefficients for the central region (|p - 0.5| <= 0.425)
constexpr double A0 =  3.3871328727963666080e0;
constexpr double A1 =  1.3314166789178437745e+2;
constexpr double A2 =  1.9715909503065514427e+3;
constexpr double A3 =  1.3731693765509461125e+4;
constexpr double A4 =  4.5921953931549871457e+4;
constexpr double A5 =  6.7265770927008700853e+4;
constexpr double A6 =  3.3430575583588128105e+4;
constexpr double A7 =  2.5090809287301226727e+3;

constexpr double B1 =  4.2313330701600911252e+1;
constexpr double B2 =  6.8718700749205790830e+2;
constexpr double B3 =  5.3941960214247511077e+3;
constexpr double B4 =  2.1213794301586595867e+4;
constexpr double B5 =  3.9307895800092710610e+4;
constexpr double B6 =  2.8729085735721942674e+4;
constexpr double B7 =  5.2264952788528545610e+3;

// Coefficients for the tail regions (0.425 < |p - 0.5| < 0.5)
constexpr double C0 =  1.42343711074721209650e0;
constexpr double C1 =  4.63033784615654529590e0;
constexpr double C2 =  5.76949722146864628717e0;
constexpr double C3 =  3.64784832476320460504e0;
constexpr double C4 =  1.27045825245236838258e0;
constexpr double C5 =  2.41780725177450611770e-1;
constexpr double C6 =  2.27001535109994502416e-2;
constexpr double C7 =  7.74545433931955401503e-4;

constexpr double D1 =  2.05319162663775882187e0;
constexpr double D2 =  1.67638483950684182780e0;
constexpr double D3 =  6.89767334985100004550e-1;
constexpr double D4 =  1.48103976427480074590e-1;
constexpr double D5 =  1.51986665636164571966e-2;
constexpr double D6 =  5.47593808499534494600e-4;
constexpr double D7 =  1.05075007164441684324e-9;

// Coefficients for the far tail (|p - 0.5| >= 0.5, i.e. p < 0.02425 or p > 0.97575)
constexpr double E0 =  6.65790464350110377720e0;
constexpr double E1 =  5.46378491116411436990e0;
constexpr double E2 =  1.78482653991729133580e0;
constexpr double E3 =  2.96560571828504891230e-1;
constexpr double E4 =  2.65321895265761230930e-2;
constexpr double E5 =  1.24266094738807843860e-3;
constexpr double E6 =  2.71155556874348757815e-5;
constexpr double E7 =  2.01033439929228813265e-7;

constexpr double F1 =  5.99832206555887937690e-1;
constexpr double F2 =  1.36929880922735805310e-1;
constexpr double F3 =  1.48753612908506508940e-2;
constexpr double F4 =  7.86869131145613259100e-4;
constexpr double F5 =  1.84631831751005468180e-5;
constexpr double F6 =  1.42151175831644588870e-7;
constexpr double F7 =  2.04426310338993978564e-15;

constexpr double SPLIT1 = 0.425;
constexpr double SPLIT2 = 5.0;
constexpr double CONST1 = 0.180625;   // 0.5 - SPLIT1^2 / 2
constexpr double CONST2 = 1.6;

[[gnu::always_inline]] inline double qnorm(double p)
{
  if (p <= 0.0) return -INFINITY;
  if (p >= 1.0) return  INFINITY;

  const double q = p - 0.5;

  if (std::fabs(q) <= SPLIT1) {
    // Central region
    const double r = CONST1 - q * q;
    return q * (((((((A7*r+A6)*r+A5)*r+A4)*r+A3)*r+A2)*r+A1)*r+A0) /
      (((((((B7*r+B6)*r+B5)*r+B4)*r+B3)*r+B2)*r+B1)*r+1.0);
  }

  // Tail regions — work with the smaller of p and 1-p
  double r = std::sqrt(-std::log(q < 0.0 ? p : 1.0 - p));
  double x;

  if (r <= SPLIT2) {
    r -= CONST2;
    x = (((((((C7*r+C6)*r+C5)*r+C4)*r+C3)*r+C2)*r+C1)*r+C0) /
      (((((((D7*r+D6)*r+D5)*r+D4)*r+D3)*r+D2)*r+D1)*r+1.0);
  } else {
    r -= SPLIT2;
    x = (((((((E7*r+E6)*r+E5)*r+E4)*r+E3)*r+E2)*r+E1)*r+E0) /
      (((((((F7*r+F6)*r+F5)*r+F4)*r+F3)*r+F2)*r+F1)*r+1.0);
  }

  return q < 0.0 ? -x : x;
}

} // namespace as241

inline double cpp_qnorm_std(double p)
{
  return as241::qnorm(p);
}
#endif // PNORM_MODE != 0

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
#define QNORM_STD(p)              pnorm_detail::cpp_qnorm_std((p))
#elif PNORM_MODE == 1
#define PNORM_STD(x, lower, logp) pnorm_detail::hart_pnorm_std((x), (lower), (logp))
#define QNORM_STD(p)              pnorm_detail::cpp_qnorm_std((p))
#else
#define PNORM_STD(x, lower, logp) pnorm_detail::exact_pnorm_std((x), (lower), (logp))
#define QNORM_STD(p)              pnorm_detail::exact_qnorm_std((p))
#endif

// ===========================================================================
// PLNORM(t, m, s) — log-normal CDF, lower tail, no log
// DLNORM(t, m, s) — log-normal PDF (always exact; no R call needed)
// ===========================================================================

#define PLNORM(t, m, s) PNORM_STD((std::log(t) - (m)) / (s), true, false)
#define DLNORM(t, m, s) pnorm_detail::fast_dlnorm((t), (m), (s))


// Upper tail phi(x) = 1 - phi(x), without log. Used by pigt0 to avoid log(pnorm).
inline double pnorm_upper(double x)
{
#if   PNORM_MODE == 2
  return 1.0 - pnorm_detail::as7126::phi(x);
#elif PNORM_MODE == 1
  return pnorm_detail::hart::phi(-x);   // phi(-x) == upper tail
#else
  return R::pnorm(x, 0.0, 1.0, false, false);
#endif
}


// ---------------------------------------------------------------------------
// fast_dnorm — standard normal PDF and general normal PDF
// Exact formula; no R API needed.
// ---------------------------------------------------------------------------

namespace pnorm_detail {

[[gnu::always_inline]] inline double fast_dnorm_std(double x)
{
  constexpr double INV_SQRT_2PI = 0.3989422804014327;
  return INV_SQRT_2PI * std::exp(-0.5 * x * x);
}

[[gnu::always_inline]] inline double fast_dnorm(double x, double mean, double sd)
{
  return fast_dnorm_std((x - mean) / sd) / sd;
}

} // namespace pnorm_detail

#define DNORM_STD(x)          pnorm_detail::fast_dnorm_std((x))
#define DNORM(x, mean, sd)    pnorm_detail::fast_dnorm((x), (mean), (sd))

#endif // PNORM_UTILS_H
