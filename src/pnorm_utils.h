// pnorm_utils.h
// Unified pnorm/plnorm/dlnorm implementation for all race models.
//
// Define USE_FAST_PNORM (e.g. -DUSE_FAST_PNORM in Makevars) to use the
// branchless A&S approximation in the hot path.  Default: exact R::pnorm.
//
// The fast path provides two primitives:
//   fast_phi(x)      — Φ(x), lower-tail, no log  [vectorisable: pure arithmetic]
//   fast_log_phi(x)  — log Φ(x)                  [vectorisable: arithmetic + log]
//
// All PNORM_STD(x, lower, logp) call sites in digt_core / pigt_core use one
// of exactly four combinations; each is mapped explicitly below so that
// lower/logp are never silently ignored.

#ifndef PNORM_UTILS_H
#define PNORM_UTILS_H

#include <cmath>
#include <Rcpp.h>

namespace pnorm_detail {

// ---------------------------------------------------------------------------
// A&S approximation core — erf via rational polynomial (Horner form)
// No branches on the polynomial path; only one copysign + one exp.
// Vectorisable by Apple clang when inlined into a contiguous loop.
// Max absolute error: ~1.5e-7 over the full real line.
// ---------------------------------------------------------------------------

inline double fast_erf(double x)
{
  // A&S 7.1.26 — rational approximation to erf(x) for x >= 0,
  // extended to all x via symmetry.
  constexpr double p  = 0.3275911;
  constexpr double a1 =  0.254829592;
  constexpr double a2 = -0.284496736;
  constexpr double a3 =  1.421413741;
  constexpr double a4 = -1.453152027;
  constexpr double a5 =  1.061405429;

  const double ax = std::fabs(x);
  const double t  = 1.0 / (1.0 + p * ax);
  const double poly = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
  const double e    = std::exp(-ax * ax);
  // erf(ax) ≈ 1 - poly * e; extend to negative x via copysign
  return std::copysign(1.0 - poly * e, x);
}

// Φ(x) = 0.5 * (1 + erf(x / sqrt(2)))
inline double fast_phi(double x)
{
  constexpr double inv_sqrt2 = 0.70710678118654752440;
  return 0.5 * (1.0 + fast_erf(x * inv_sqrt2));
}

// log Φ(x) — used in pigt_core t2 term (logp=true call site)
// Delegates to std::log(fast_phi(x)); clamped to avoid -Inf in hot path.
inline double fast_log_phi(double x)
{
  const double p = fast_phi(x);
  // For x < -8 the approximation saturates to ~0; clamp to a finite floor.
  return (p > 0.0) ? std::log(p) : -1000.0;
}

// ---------------------------------------------------------------------------
// Exact fallbacks
// ---------------------------------------------------------------------------

inline double exact_pnorm_std(double x, bool lower, bool logp)
{
  return R::pnorm(x, 0.0, 1.0, lower, logp);
}

// ---------------------------------------------------------------------------
// Rational approximation for the upper tail — used by pigt0 (lower=false)
// Reuses fast_phi via symmetry: Φ_upper(x) = 1 - Φ(x) = Φ(-x)
// ---------------------------------------------------------------------------

inline double fast_phi_upper(double x)
{
  return fast_phi(-x);
}

} // namespace pnorm_detail

// ---------------------------------------------------------------------------
// PNORM_STD(x, lower, logp)
//
// All four (lower, logp) combinations used in the codebase are mapped
// explicitly so that the fast path never silently ignores an argument.
//
//   (true,  false) — Φ(x)          — digt_core t2, pigt_core t4
//   (true,  true)  — log Φ(x)      — pigt_core t2
//   (false, false) — 1 - Φ(x)      — pigt0
//   (false, true)  — log(1 - Φ(x)) — (not currently used; safe fallback)
// ---------------------------------------------------------------------------

#ifdef USE_FAST_PNORM

// Dispatch at compile time on the literal bool arguments.
// The compiler folds the dead branches away entirely.
#define PNORM_STD(x, lower, logp)                                         \
( (lower) && !(logp) ? pnorm_detail::fast_phi(x)                          \
    : (lower) &&  (logp) ? pnorm_detail::fast_log_phi(x)                  \
    : !(lower) && !(logp) ? pnorm_detail::fast_phi_upper(x)               \
    :                       std::log(pnorm_detail::fast_phi_upper(x)) )

#else

#define PNORM_STD(x, lower, logp) pnorm_detail::exact_pnorm_std((x), (lower), (logp))

#endif // USE_FAST_PNORM

  // ---------------------------------------------------------------------------
  // PLNORM(t, m, s) / DLNORM(t, m, s)
  // ---------------------------------------------------------------------------

  namespace pnorm_detail {

  inline double fast_plnorm(double t, double m, double s)
  {
    // Reuses PNORM_STD so it respects the USE_FAST_PNORM switch.
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

  // DLNORM: always the direct formula — no R call, no #ifdef needed
#define DLNORM(t, m, s) pnorm_detail::fast_dlnorm((t), (m), (s))

#endif // PNORM_UTILS_H
