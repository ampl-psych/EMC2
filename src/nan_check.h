#ifndef NAN_CHECK_H
#define NAN_CHECK_H

#include <cstdint>
#include <cstring>

// ---------------------------------------------------------------------------
// NaN / NA-safe floating-point checks
//
// R's NA_real_ is a specific NaN bit pattern (0x7FF00000000007A2), and
// ordinary NaN checks are unreliable when -ffast-math / -ffinite-math-only
// is enabled: the compiler assumes finite values and may optimise away
// x != x, std::isnan(x), and even __builtin_isnan(x) (Clang).
//
// These helpers operate on raw bit patterns via memcpy, which is fully
// immune to floating-point optimisation flags and catches both NaN and
// R's NA_real_ (since NA_real_ is itself a NaN).
// ---------------------------------------------------------------------------

/// Returns true if x is any NaN (including R's NA_real_).
/// Safe under -ffast-math / -ffinite-math-only.
inline bool is_nan(double x) {
  uint64_t bits;
  std::memcpy(&bits, &x, sizeof(bits));
  // NaN: exponent bits all 1, mantissa non-zero
  return (bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL
  && (bits & 0x000FFFFFFFFFFFFFULL) != 0ULL;
}

/// Returns true if x is finite (not NaN, not Inf, not -Inf).
/// Safe under -ffast-math / -ffinite-math-only.
inline bool is_finite(double x) {
  uint64_t bits;
  std::memcpy(&bits, &x, sizeof(bits));
  // Finite: exponent bits not all 1
  return (bits & 0x7FF0000000000000ULL) != 0x7FF0000000000000ULL;
}

/// Returns true if x is +Inf or -Inf.
/// Safe under -ffast-math / -ffinite-math-only.
inline bool is_inf(double x) {
  uint64_t bits;
  std::memcpy(&bits, &x, sizeof(bits));
  // Inf: exponent all 1, mantissa exactly 0
  return (bits & 0x7FFFFFFFFFFFFFFFULL) == 0x7FF0000000000000ULL;
}

#endif // NAN_CHECK_H
