#ifndef DBM_MATH_H
#define DBM_MATH_H

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "nan_check.h"

// Mean, mode, and log-precision of Beta distribution
// Used by Beta-binomial kernels
inline double beta_mode(double a, double b) {
  if (a == b)              return 0.5;
  if (a < 1.0 || b < 1.0) return (a < b) ? 0.0 : 1.0;
  return (a - 1.0) / (a + b - 2.0);
}

inline double beta_mean(double a, double b) {
  return a / (a + b);
}

// inline double beta_variance(double a, double b) {
//   const double ab_sum = a + b;
//   return (a * b) / (ab_sum * ab_sum * (ab_sum + 1.0));
// }

inline double beta_log_precision(double a, double b) {
  static const double ceiling = -std::log(std::numeric_limits<double>::epsilon());
  const double ab_sum = a + b;
  const double result = 2.0 * std::log(ab_sum) +
    std::log1p(ab_sum) - std::log(a) - std::log(b);
  if (!is_finite(result) || result > ceiling) {
    return ceiling;
  }
  return result;
}

// Functions for discretised Beta distribution, used by DBM kernel
inline double dbeta_log(double x, double a, double b) {
  if (x <= 0.0 || x >= 1.0) return -std::numeric_limits<double>::infinity();
  return (a - 1.0) * std::log(x)
    + (b - 1.0) * std::log(1.0 - x)
    - std::lgamma(a) - std::lgamma(b) + std::lgamma(a + b);
}

inline double dbeta_val(double x, double a, double b) {
  return std::exp(dbeta_log(x, a, b));
}

inline double normalise_inplace(std::vector<double>& v) {
  double s = 0.0;
  for (double x : v) s += x;
  if (s > 0.0) { double inv = 1.0 / s; for (double& x : v) x *= inv; }
  return s;
}

inline double mean_discrete(const std::vector<double>& x,
                            const std::vector<double>& w) {
  double m = 0.0;
  const int n = static_cast<int>(x.size());
  for (int i = 0; i < n; ++i) m += x[i] * w[i];
  return m;
}

inline double mode_discrete(const std::vector<double>& x,
                            const std::vector<double>& w) {
  int idx = static_cast<int>(
    std::max_element(w.begin(), w.end()) - w.begin());
  return x[idx];
}

inline double var_discrete(const std::vector<double>& x,
                           const std::vector<double>& w) {
  double m = 0.0;
  double sx2 = 0.0;
  const int n = static_cast<int>(x.size());
  for (int i = 0; i < n; ++i) {
    m += x[i] * w[i];
    sx2 += x[i] * x[i] * w[i];
  }
  const double variance = sx2 - (m * m);
  return std::max(0.0, variance);
}

inline double log_precision_discrete(const std::vector<double>& x,
                                     const std::vector<double>& w) {
  static const double ceiling = -std::log(std::numeric_limits<double>::epsilon());
  const double result = -std::log(var_discrete(x, w));
  if (!is_finite(result) || result > ceiling) {
    return ceiling;
  }
  return result;
}


// Shannon surprise for binary observations, used by Beta-binomial and DBM kernels
inline double shannon_surprise(double pred, double obs) {
  static const double inv_ln2 = 1.0 / std::log(2.0);
  static const double tiny    = std::numeric_limits<double>::min();
  const double pred_safe = std::min(std::max(pred, tiny), 1.0 - tiny);
  const double like      = (obs == 1.0) ? pred_safe : (1.0 - pred_safe);
  return -std::log(like) * inv_ln2;
}

#endif // DBM_MATH_H
