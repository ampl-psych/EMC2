#ifndef DBM_MATH_H
#define DBM_MATH_H

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

// For DBM kernels
inline double dbeta_log(double x, double a, double b) {
  if (x <= 0.0 || x >= 1.0) return -std::numeric_limits<double>::infinity();
  return (a - 1.0) * std::log(x)
    + (b - 1.0) * std::log(1.0 - x)
    - std::lgamma(a) - std::lgamma(b) + std::lgamma(a + b);
}

inline double dbeta_val(double x, double a, double b) {
  return std::exp(dbeta_log(x, a, b));
}

inline double beta_mode(double a, double b) {
  if (a == b)              return 0.5;
  if (a < 1.0 || b < 1.0) return (a < b) ? 0.0 : 1.0;
  return (a - 1.0) / (a + b - 2.0);
}

inline double beta_mean(double a, double b) {
  return a / (a + b);
}

inline double normalise_inplace(std::vector<double>& v) {
  double s = 0.0;
  for (double x : v) s += x;
  if (s > 0.0) { double inv = 1.0 / s; for (double& x : v) x *= inv; }
  return s;
}

inline double mean_discrete(const std::vector<double>& x,
                            const std::vector<double>& w) {
  double s = 0.0;
  const int n = static_cast<int>(x.size());
  for (int i = 0; i < n; ++i) s += x[i] * w[i];
  return s;
}

inline double mode_discrete(const std::vector<double>& x,
                            const std::vector<double>& w) {
  int idx = static_cast<int>(
    std::max_element(w.begin(), w.end()) - w.begin());
  return x[idx];
}

inline double shannon_surprise(double pred, double obs) {
  static const double inv_ln2 = 1.0 / std::log(2.0);
  static const double tiny    = std::numeric_limits<double>::min();
  const double pred_safe = std::min(std::max(pred, tiny), 1.0 - tiny);
  const double like      = (obs == 1.0) ? pred_safe : (1.0 - pred_safe);
  return -std::log(like) * inv_ln2;
}

#endif // DBM_MATH_H
