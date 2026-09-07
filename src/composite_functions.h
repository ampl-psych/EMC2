#ifndef composite_functions_h
#define composite_functions_h

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include "r_constants.h"
#include "nan_check.h"

/**
 * @file composite_funs.h
 * @brief Numerically stable, thread-safe, -ffast-math-safe convenience
 *        functions for statistical computing
 *
 * This header provides composite functions that are more efficient and
 * numerically stable than their naive implementations. NaN/Inf detection
 * uses the project's is_nan / is_finite / is_inf helpers instead of
 * ISNAN / direct equality against R_NegInf / R_PosInf, since the latter
 * are not guaranteed to behave correctly under -ffast-math
 * (-ffinite-math-only assumes no operand is ever NaN/Inf and the
 * compiler is free to fold such checks away).
 *
 * The Rcpp::NumericVector overloads are convenience wrappers for use from
 * R and should not themselves be called from worker threads.
 */

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------

namespace composite_detail {
  constexpr double LOG2 = 0.693147180559945309417232121458176568;
  constexpr double quiet_nan = std::numeric_limits<double>::quiet_NaN();
}  // namespace composite_detail


// -----------------------------------------------------------------------------
// log(1 - x)
// -----------------------------------------------------------------------------

/**
 * @brief Compute log(1 - x) in a numerically stable way.
 *
 * @param x Input value.
 * @return log(1 - x), or NaN if x > 1 or x is NaN.
 */
[[gnu::always_inline]] inline double log1m(double x) {
  if (is_nan(x) || x > 1.0) { return na_real(); }
  if (is_inf(x)) { return (x > 0.0) ? na_real() : pos_inf(); }
  if (x == 1.0) { return neg_inf(); }
  return std::log1p(-x);
}


/**
 * @brief Compute log(1 - x) for a vector.
 *
 * Rcpp convenience wrapper. Do not call from worker threads.
 */
[[gnu::always_inline]] inline Rcpp::NumericVector log1m(const Rcpp::NumericVector& x) {
  const int n = x.size();
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = log1m(x[i]);
  }
  return result;
}


// -----------------------------------------------------------------------------
// log(1 + exp(x)) / softplus
// -----------------------------------------------------------------------------

/**
 * @brief Compute log(1 + exp(x)) in a numerically stable way.
 *
 * Uses a four-regime piecewise approximation (equivalent in accuracy to
 * the Stan Math implementation) rather than a single log1p(exp(x))
 * call, which loses precision or overflows outside a narrow range:
 *  - x <= -37:  exp(x) alone is accurate to double precision
 *  - x <=  18:  log1p(exp(x)) is accurate and exp(x) doesn't overflow
 *  - x <=  33.3: exp(-x) is small enough that x + exp(-x) is accurate
 *                and avoids overflow from exp(x)
 *  - x >   33.3: log(1+exp(x)) rounds to x itself at double precision
 *
 * @param x Input value.
 * @return log(1 + exp(x)).
 */
[[gnu::always_inline]] inline double log1p_exp(double x) {
  if (is_nan(x)) { return na_real(); }
  if (is_inf(x)) { return (x > 0.0) ? pos_inf() : 0.0; }
  if (x <= -37.0) { return std::exp(x) ; }
  if (x <=  18.0) { return std::log1p(std::exp(x)) ; }
  if (x <=  33.3) { return x + std::exp(-x); }
  return x;
}


/**
 * @brief Compute log(1 + exp(x)) for a vector.
 *
 * Rcpp convenience wrapper. Do not call from worker threads.
 */
[[gnu::always_inline]] inline Rcpp::NumericVector log1p_exp(const Rcpp::NumericVector& x) {
  const int n = x.size();
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = log1p_exp(x[i]);
  }
  return result;
}


// -----------------------------------------------------------------------------
// log(1 - exp(x))
// -----------------------------------------------------------------------------

/**
 * @brief Compute log(1 - exp(x)) in a numerically stable way.
 *
 * This function computes log(1 - exp(x)) for x <= 0 without intermediate
 * overflow or underflow issues. Uses expm1 (accurate for x near 0, i.e.
 * exp(x) near 1) when x is close to 0, and log1p (accurate for exp(x)
 * near 0) when x is very negative, following the crossover at -log(2)
 * where the two formulations have comparable accuracy.
 *
 * @param x Input value (must be <= 0).
 */
[[gnu::always_inline]] inline double log1m_exp(double x) {
  if (is_nan(x) || x > 0.0) { return na_real(); }
  if (is_inf(x)) { return 0.0; } // can only be neg inf, given previous line
  if (x == 0.0) { return neg_inf(); }
  // Near zero, expm1(x) gives much better accuracy than 1 - exp(x).
  if (x > -composite_detail::LOG2) { return std::log(-std::expm1(x)); }
  // Here exp(x) is sufficiently far from 1 that this is well conditioned.
  return log1m(std::exp(x));
  // return std::log1p(-std::exp(x));
}


/**
 * @brief Compute log(1 - exp(x)) for a vector.
 *
 * Rcpp convenience wrapper. Do not call from worker threads.
 */
[[gnu::always_inline]] inline Rcpp::NumericVector log1m_exp(const Rcpp::NumericVector& x) {
  const int n = x.size();
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = log1m_exp(x[i]);
  }
  return result;
}


// -----------------------------------------------------------------------------
// log(exp(a) + exp(b)) / log(sum(exp(x)))
// -----------------------------------------------------------------------------

/**
 * @brief Compute log(exp(a) + exp(b)) in a numerically stable way.
 *
 * @param a First log-scale value.
 * @param b Second log-scale value.
 * @return log(exp(a) + exp(b)).
 */
[[gnu::always_inline]] inline double log_sum_exp(double a, double b) {
  if (is_nan(a) || is_nan(b)) { return na_real(); }
  // These cases must be handled before subtracting a and b, since
  // +Inf - +Inf is NaN.
  if (is_inf(a) && a < 0.0) { return b; }
  if (is_inf(b) && b < 0.0) { return a; }
  if (is_inf(a) && a > 0.0 && is_inf(b) && b > 0.0) { return pos_inf(); }
  if (a > b) {
    return a + std::log1p(std::exp(b - a));
  }
  return b + std::log1p(std::exp(a - b));
}


/**
 * @brief Compute log(exp(a) + exp(b)) for two vectors.
 *
 * Rcpp convenience wrapper. Do not call from worker threads.
 */
[[gnu::always_inline]] inline Rcpp::NumericVector log_sum_exp(
    const Rcpp::NumericVector& a, const Rcpp::NumericVector& b
) {
  const int n = a.size();
  if (n != b.size()) { Rcpp::stop("Vectors must have the same length"); }
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = log_sum_exp(a[i], b[i]);
  }
  return result;
}


/**
 * @brief Compute log(sum(exp(x))) in a numerically stable way.
 *
 * @param x Vector of log-scale values.
 * @return log(sum(exp(x))), with -Inf for an empty vector.
 */
[[gnu::always_inline]] inline double log_sum_exp(const Rcpp::NumericVector& x) {
  const int n = x.size();
  if (n == 0) { return neg_inf(); }
  if (n == 1) { return x[0]; }
  double max_val = neg_inf();
  bool any_nan = false;
  for (double val : x) {
    if (is_nan(val)) { any_nan = true; break; }
    if (val > max_val) { max_val = val; }
  }
  if (any_nan) { return na_real(); }
  if (is_inf(max_val) && max_val > 0.0) { return pos_inf(); }
  if (is_inf(max_val) && max_val < 0.0) { return neg_inf(); }
  double sum = 0.0;
  for (double val : x) {
    if (!(is_inf(val) && val < 0.0)) {
      sum += std::exp(val - max_val);
    }
  }
  return max_val + std::log(sum);
}


// -----------------------------------------------------------------------------
// log(exp(a) - exp(b))
// -----------------------------------------------------------------------------

/**
 * @brief Compute log(exp(a) - exp(b)) in a numerically stable way.
 *
 * Requires a >= b.
 *
 * @param a First log-scale value.
 * @param b Second log-scale value.
 * @return
 *   log(exp(a) - exp(b)) if a > b,
 *   -Inf if a == b,
 *   na_real if a < b or either argument is NaN.
 */
[[gnu::always_inline]] inline double log_diff_exp(double a, double b) {
  if (is_nan(a) || is_nan(b) || a < b) { return na_real(); }
  if (is_inf(a) && a > 0.0) {
    if (is_inf(b) && b > 0.0) { return na_real(); }
    return pos_inf();
  }
  if (a == b) { return neg_inf(); }
  if (is_inf(b) && b < 0.0) { return a; }
  const double diff = b - a;
  // When b is close to a, exp(diff) is close to 1 and
  // log1m_exp() avoids catastrophic cancellation.
  if (diff > -composite_detail::LOG2) { return a + log1m_exp(diff); }
  // When exp(diff) is small, 1 - exp(diff) is well conditioned.
  return a + std::log1p(-std::exp(diff));
}


/**
 * @brief Compute log(exp(a) - exp(b)) for two vectors.
 *
 * Rcpp convenience wrapper. Do not call from worker threads.
 */
[[gnu::always_inline]] inline Rcpp::NumericVector log_diff_exp(
    const Rcpp::NumericVector& a, const Rcpp::NumericVector& b
) {
  const int n = a.size();
  if (n != b.size()) { Rcpp::stop("Vectors must have the same length"); }
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = log_diff_exp(a[i], b[i]);
  }
  return result;
}


// -----------------------------------------------------------------------------
// Log mixture
// -----------------------------------------------------------------------------

/**
 * @brief Compute the log density of a two-component mixture.
 *
 * Computes:
 *
 *   log(theta * exp(lambda1) + (1 - theta) * exp(lambda2))
 *
 * where theta must lie in [0, 1].
 *
 * @param theta Mixing proportion.
 * @param lambda1 Log density of first component.
 * @param lambda2 Log density of second component.
 * @return Log mixture density, or NaN for invalid input.
 */
[[gnu::always_inline]] inline double log_mix(double theta, double lambda1, double lambda2) {
  if (is_nan(theta) || is_nan(lambda1) || is_nan(lambda2)) { return na_real(); }
  if (theta < 0.0 || theta > 1.0) { return na_real(); }
  if (theta == 0.0) { return lambda2; }
  if (theta == 1.0) { return lambda1; }
  bool l1_neg_inf = is_inf(lambda1) && lambda1 < 0.0;
  bool l2_neg_inf = is_inf(lambda2) && lambda2 < 0.0;
  if (l1_neg_inf && l2_neg_inf) { return neg_inf(); }
  if (l1_neg_inf) { return log1m(theta) + lambda2; }
  if (l2_neg_inf) { return std::log(theta) + lambda1; }
  return log_sum_exp(std::log(theta) + lambda1, log1m(theta) + lambda2);
}


/**
 * @brief Compute log mixture density for vectors.
 *
 * Rcpp convenience wrapper. Do not call from worker threads.
 */
[[gnu::always_inline]] inline Rcpp::NumericVector log_mix(
    const Rcpp::NumericVector& theta,
    const Rcpp::NumericVector& lambda1,
    const Rcpp::NumericVector& lambda2
) {
  const int n = theta.size();
  if (n != lambda1.size() || n != lambda2.size()) {
    Rcpp::stop("All vectors must have the same length");
  }
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = log_mix(theta[i], lambda1[i], lambda2[i]);
  }
  return result;
}


// -----------------------------------------------------------------------------
// log(mean(exp(x)))
// -----------------------------------------------------------------------------

/**
 * @brief Compute log(mean(exp(x))) in a numerically stable way.
 *
 * @param x Vector of log-scale values.
 * @return log(mean(exp(x))), with -Inf for an empty vector.
 */
[[gnu::always_inline]] inline double log_mean_exp(const Rcpp::NumericVector& x) {
  const int n = x.size();
  if (n == 0) {return neg_inf();}
  return log_sum_exp(x) - std::log(static_cast<double>(n));
}

#endif // composite_functions_h
