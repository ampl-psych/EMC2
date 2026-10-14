#ifndef DBM_MATH_H
#define DBM_MATH_H

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "nan_check.h"

// Mean, mode, and log-precision of Beta distribution
// overwrite inline
inline void beta_summary(double a, double b,
                         double& mean, double& mode,
                         double& log_precision) {
  static const double eps = std::numeric_limits<double>::epsilon();
  static const double ceiling = -std::log(eps);

  const double s = a + b;
  const double inv_s = 1.0 / s;
  const double p = a * inv_s;
  const double q = b * inv_s;

  // Mean
  mean = p;

  // Mode determination
  if(a == b) {
    mode = 0.5;
  } else if(a < 1.0 || b < 1.0) {
    mode = (a < b) ? 0.0 : 1.0;
  } else {
    mode = (a - 1.0) / (s - 2.0);
  }

  const double variance = (p * q) / (s + 1.0);
  log_precision = !(variance > eps) ? ceiling : -std::log(variance);
}

// Functions for discretised Beta distribution, used by DBM kernel
inline double dbeta_log(double x, double a, double b) {
  if (x <= 0.0 || x >= 1.0) return -std::numeric_limits<double>::infinity();
  return (a - 1.0) * std::log(x) + (b - 1.0) * std::log(1.0 - x) - std::lgamma(a) - std::lgamma(b) + std::lgamma(a + b);
}

inline double dbeta_val(double x, double a, double b) {
  return std::exp(dbeta_log(x, a, b));
}

inline double normalise_inplace(std::vector<double>& v) {
  // vectorise both passes
  const std::size_t n = v.size();
  double* const data = v.data();

  double s = 0.0;

#pragma omp simd reduction(+:s)
  for (std::size_t i = 0; i < n; ++i) {
    s += data[i];
  }

  if (s > 0.0) {
    const double inv = 1.0 / s;

#pragma omp simd
    for (std::size_t i = 0; i < n; ++i) {
      data[i] *= inv;
    }
  }

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
  const int n = static_cast<int>(x.size());
  double m = 0.0;
  for (int i = 0; i < n; ++i) {
    m += x[i] * w[i];
  }
  double v = 0.0;
  for (int i = 0; i < n; ++i) {
    const double d = x[i] - m;
    v += d * d * w[i];
  }
  return std::max(0.0, v);
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

// inlined inplace calculation of discrete summary stats
inline void discrete_summary(const std::vector<double>& x,
                             const std::vector<double>& w,
                             double& mean,
                             double& mode,
                             double& log_precision) {
  static const double ceiling = -std::log(std::numeric_limits<double>::epsilon());

    const int n = static_cast<int>(x.size());

    double mean_sum = 0.0;
    double max_weight = w[0];

// This is vectorisable
#pragma omp simd reduction(+:mean_sum) reduction(max:max_weight)
    for (int i = 0; i < n; ++i) {
      mean_sum += x[i] * w[i];
      max_weight = std::max(max_weight, w[i]);
    }

    mean = mean_sum;

    // Recover the corresponding mode value
    const auto it = std::find(w.begin(), w.end(), max_weight);
    mode = x[static_cast<std::size_t>(it - w.begin())];


    // Second pass: centered variance.
    double variance = 0.0;

// This is vectorisable as well
#pragma omp simd reduction(+:variance)
    for (int i = 0; i < n; ++i) {
      const double d = x[i] - mean;
      variance += d * d * w[i];
    }
    const double result = -std::log(variance);
    log_precision = (!is_finite(result) || result > ceiling) ? ceiling : result;
}


// Shannon surprise for binary observations, used by Beta-binomial and DBM kernels
inline double shannon_surprise(double pred, double obs) {
  static const double inv_ln2 = 1.0 / std::log(2.0);
  static const double tiny    = std::numeric_limits<double>::min();
  const double pred_safe = std::min(std::max(pred, tiny), 1.0 - tiny);
  const double like      = (obs == 1.0) ? pred_safe : (1.0 - pred_safe);
  return -std::log(like) * inv_ln2;
}

// inline double shannon_entropy(double pred) {
//   static const double inv_ln2 = 1.0 / std::log(2.0);
//   static const double tiny    = std::numeric_limits<double>::min();
//   const double pred_safe = std::min(std::max(pred, tiny), 1.0 - tiny);
//   const double entropy = pred_safe * std::log(pred_safe)
//     + (1.0 - pred_safe) * std::log(1.0 - pred_safe);
//   return -entropy * inv_ln2;
// }


// Normalised discretised Beta(a, b) on a grid over [0, 1].
// log_x[i] = log(x_i), log_1mx[i] = log1p(-x_i)   (precomputed once per grid)
// Endpoints get zero mass, as dbeta_log returned -Inf at x <= 0 and x >= 1.
inline void beta_prior_grid(const std::vector<double>& log_x,
                            const std::vector<double>& log_1mx,
                            double a, double b,
                            std::vector<double>& out) {
  const int gs = static_cast<int>(out.size());

  // Invalid shape parameters: make it visible rather than silently returning zeros
  if (!(a > 0.0) || !(b > 0.0) || !is_finite(a) || !is_finite(b)) {
    std::fill(out.begin(), out.end(), std::numeric_limits<double>::quiet_NaN());
    return;
  }

  double mx = -std::numeric_limits<double>::infinity();
  for (int i = 1; i < gs - 1; ++i) {
    out[i] = (a - 1.0) * log_x[i] + (b - 1.0) * log_1mx[i];
    mx = std::max(mx, out[i]);
  }
  out[0] = out[gs - 1] = 0.0;
  for (int i = 1; i < gs - 1; ++i) out[i] = std::exp(out[i] - mx);
  normalise_inplace(out);
}

#endif // DBM_MATH_H
