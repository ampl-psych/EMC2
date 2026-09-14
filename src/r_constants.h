#ifndef R_CONSTANTS_H
#define R_CONSTANTS_H

/**
 * @file r_constants.h
 * @brief Native C++ replacements for selected R floating-point constants.
 *
 * Provides R-independent replacements for:
 *
 *   NA_REAL
 *   R_NegInf
 *   R_PosInf
 *
 * The definitions are intended for numerical C++ code that may also be
 * compiled as part of an R/Rcpp package, but should not otherwise depend
 * on the R API.
 */

#include <cstdint>
#include <cstring>
#include <limits>

// positive infinity
[[gnu::always_inline]] inline double pos_inf() {
  return std::numeric_limits<double>::infinity();
}

// negative infinity
[[gnu::always_inline]] inline double neg_inf() {
  return -std::numeric_limits<double>::infinity();
}

// R-compatible NA value for double precision.
// R represents NA_REAL using a particular IEEE-754 NaN payload. The following
// reproduces that bit pattern without depending on R headers.
// The representation used by R for NA_REAL is 0x7ff00000000007a2 on IEEE-754
// double-precision platforms.
[[gnu::always_inline]] inline double na_real() {
  constexpr std::uint64_t bits = UINT64_C(0x7ff00000000007a2);
  double value;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

#endif  // R_CONSTANTS_H
