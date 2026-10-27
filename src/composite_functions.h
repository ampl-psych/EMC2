#ifndef composite_functions_h
#define composite_functions_h

#include <Rcpp.h>
#include <cmath>
#include "r_constants.h"
#include "nan_check.h"
#include "hcubature.h"

namespace composite_detail {
  constexpr double LOG2 = 0.693147180559945309417232121458176568;
}

/**
 * @file composite_funs.h
 * @brief Numerically stable convenience functions for statistical computing
 *
 * This header provides composite functions that are more efficient and
 * numerically stable than their naive implementations.
 */

/**
 * @brief Compute log(1 - x) in a numerically stable way
 *
 * This function computes log(1 - x) without precision loss that would
 * occur with the naive computation when x is close to 0.
 *
 * @param x Input value (should be < 1 for real result)
 * @return log(1 - x)
 */
inline double log1m(double x) {
  if (ISNAN(x)) return NA_REAL;
  if (x >= 1.0) return R_NegInf;  // log(1 - x) = log(0 or negative)
  if (x == R_NegInf) return 0.0;  // log(1 - (-Inf)) = log(Inf) = Inf, but we return 0 for log(1)

  return std::log1p(-x);
}

/**
 * @brief Compute log(1 - x) for vectors in a numerically stable way
 *
 * Vectorized version of log1m for element-wise operations.
 *
 * @param x Vector of input values
 * @return Vector of log(1 - x[i])
 */
inline Rcpp::NumericVector log1m(const Rcpp::NumericVector& x) {
  int n = x.size();
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = log1m(x[i]);
  }
  return result;
}

/**
 * @brief Compute log(1 + exp(x)) in a numerically stable way
 *
 * This function computes log(1 + exp(x)) without intermediate overflow
 * issues that would occur with the naive computation. Also known as
 * the "softplus" function in machine learning.
 *
 * @param x Input value
 * @return log(1 + exp(x))
 */
inline double log1p_exp(double x) {
  if (ISNAN(x)) return NA_REAL;
  if (x == R_PosInf) return R_PosInf;
  if (x == R_NegInf) return 0.0;  // log(1 + 0) = 0

  if (x > 0.0) {
    return x + std::log1p(std::exp(-x));
  } else {
    return std::log1p(std::exp(x));
  }
}

/**
 * @brief Compute log(1 + exp(x)) for vectors in a numerically stable way
 *
 * Vectorized version of log1p_exp for element-wise operations.
 *
 * @param x Vector of input values
 * @return Vector of log(1 + exp(x[i]))
 */
inline Rcpp::NumericVector log1p_exp(const Rcpp::NumericVector& x) {
  int n = x.size();
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = log1p_exp(x[i]);
  }
  return result;
}

/**
 * @brief Compute log(1 - exp(x)) in a numerically stable way
 *
 * This function computes log(1 - exp(x)) for x <= 0 without intermediate
 * overflow or underflow issues. Uses expm1 near zero and log1p for
 * well-separated terms. Returns NA_REAL for x > 0.
 *
 * @param x Input value (must be <= 0)
 * @return log(1 - exp(x)) if x < 0, log(0) = -Inf if x = 0, NA_REAL if x > 0
 */
[[gnu::always_inline]] inline double log1m_exp(double x) {
  if (is_nan(x) || x > 0.0) { return na_real(); }
  if (is_neg_inf(x)) { return 0.0; }
  if (x == 0.0) { return neg_inf(); }
  if (x > -composite_detail::LOG2) { return std::log(-std::expm1(x)); }
  return log1m(std::exp(x));
}

/**
 * @brief Compute log(1 - exp(x)) for vectors in a numerically stable way
 *
 * Vectorized version of log1m_exp for element-wise operations.
 *
 * @param x Vector of input values
 * @return Vector of log(1 - exp(x[i]))
 */
inline Rcpp::NumericVector log1m_exp(const Rcpp::NumericVector& x) {
  int n = x.size();
  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = log1m_exp(x[i]);
  }
  return result;
}

/**
 * @brief Compute log(exp(a) + exp(b)) in a numerically stable way
 *
 * This function computes the logarithm of the sum of exponentials without
 * intermediate overflow or underflow issues that would occur with the
 * naive computation log(exp(a) + exp(b)).
 *
 * @param a First log-scale value
 * @param b Second log-scale value
 * @return log(exp(a) + exp(b))
 */
inline double log_sum_exp(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  if (a > b) {
    return a + std::log1p(std::exp(b - a));
  } else {
    return b + std::log1p(std::exp(a - b));
  }
}

/**
 * @brief Compute log(exp(a) + exp(b)) for vectors in a numerically stable way
 *
 * Vectorized version of log_sum_exp for element-wise operations.
 *
 * @param a First vector of log-scale values
 * @param b Second vector of log-scale values
 * @return Vector of log(exp(a[i]) + exp(b[i]))
 */
inline Rcpp::NumericVector log_sum_exp(const Rcpp::NumericVector& a,
                                       const Rcpp::NumericVector& b) {
  int n = a.size();
  if (n != b.size()) {
    Rcpp::stop("Vectors must have the same length");
  }

  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = log_sum_exp(a[i], b[i]);
  }
  return result;
}

/**
 * @brief Compute log(sum(exp(x))) for a vector in a numerically stable way
 *
 * This function computes the log-sum-exp of a vector, which is equivalent
 * to log(sum(exp(x))) but avoids overflow/underflow issues.
 *
 * @param x Vector of log-scale values
 * @return log(sum(exp(x[i])))
 */
inline double log_sum_exp(const Rcpp::NumericVector& x) {
  if (x.size() == 0) return R_NegInf;
  if (x.size() == 1) return x[0];

  // following line could replace the for-loop, but would require #include <algorithm>
  // double max_val = *std::max_element(x.begin(), x.end());
  double max_val = R_NegInf;
  for (double val : x) {
    if (val > max_val) {
      max_val = val;
    }
  }
  if (max_val == R_NegInf) return R_NegInf;

  double sum = 0.0;
  for (double val : x) {
    if (val != R_NegInf) {
      sum += std::exp(val - max_val);
    }
  }

  return max_val + std::log(sum);
}

/**
 * @brief Compute log(exp(a) - exp(b)) in a numerically stable way
 *
 * This function computes the logarithm of the difference of exponentials
 * without intermediate overflow or underflow issues. Requires a >= b.
 *
 * @param a First log-scale value (must be >= b)
 * @param b Second log-scale value (must be <= a)
 * @return log(exp(a) - exp(b)) if a > b, -Inf if a == b,
 * NA_REAL if a < b or either input is NaN; +Inf for a=+Inf and b<+Inf.
 */
[[gnu::always_inline]] inline double log_diff_exp(double a, double b) {
  if (is_nan(a) || is_nan(b) || a < b) { return na_real(); }
  if (is_pos_inf(a)) {
    if (is_pos_inf(b)) { return na_real(); }
    return pos_inf();
  }
  if (a == b) { return neg_inf(); }
  if (is_neg_inf(b)) { return a; }
  const double diff = b - a;
  // Close arguments require expm1 to retain their difference.
  if (diff > -composite_detail::LOG2) { return a + log1m_exp(diff); }
  return a + std::log1p(-std::exp(diff));
}

/**
 * @brief Compute log(exp(a) - exp(b)) for vectors in a numerically stable way
 *
 * Vectorized version of log_diff_exp for element-wise operations.
 *
 * @param a First vector of log-scale values
 * @param b Second vector of log-scale values
 * @return Vector of log(exp(a[i]) - exp(b[i]))
 */
inline Rcpp::NumericVector log_diff_exp(const Rcpp::NumericVector& a,
                                        const Rcpp::NumericVector& b) {
  int n = a.size();
  if (n != b.size()) {
    Rcpp::stop("Vectors must have the same length");
  }

  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = log_diff_exp(a[i], b[i]);
  }
  return result;
}

/**
 * @brief Compute log mixture density in a numerically stable way
 *
 * This function computes log(theta * exp(lambda1) + (1-theta) * exp(lambda2))
 * which is the log density of a mixture distribution with mixing proportion
 * theta and component log densities lambda1 and lambda2.
 *
 * @param theta Mixing proportion (must be in [0, 1])
 * @param lambda1 Log density of first component
 * @param lambda2 Log density of second component
 * @return log(theta * exp(lambda1) + (1-theta) * exp(lambda2))
 */
inline double log_mix(double theta, double lambda1, double lambda2) {
  // Input validation
  if (ISNAN(theta) || ISNAN(lambda1) || ISNAN(lambda2)) return NA_REAL;
  if (theta < 0.0 || theta > 1.0) return NA_REAL;

  // Handle boundary cases
  if (theta == 0.0) return lambda2;
  if (theta == 1.0) return lambda1;

  // Handle infinite log densities
  if (lambda1 == R_NegInf && lambda2 == R_NegInf) return R_NegInf;
  if (lambda1 == R_NegInf) return log1m(theta) + lambda2;  // Using log1m for log(1-theta)
  if (lambda2 == R_NegInf) return std::log(theta) + lambda1;

  // General case: log(theta * exp(lambda1) + (1-theta) * exp(lambda2))
  return log_sum_exp(std::log(theta) + lambda1, log1m(theta) + lambda2);  // Using log1m for log(1-theta)
}

/**
 * @brief Compute log mixture density for vectors in a numerically stable way
 *
 * Vectorized version of log_mix for element-wise operations.
 *
 * @param theta Vector of mixing proportions (each must be in [0, 1])
 * @param lambda1 Vector of log densities for first components
 * @param lambda2 Vector of log densities for second components
 * @return Vector of log mixture densities
 */
inline Rcpp::NumericVector log_mix(const Rcpp::NumericVector& theta,
                                   const Rcpp::NumericVector& lambda1,
                                   const Rcpp::NumericVector& lambda2) {
  int n = theta.size();
  if (n != lambda1.size() || n != lambda2.size()) {
    Rcpp::stop("All vectors must have the same length");
  }

  Rcpp::NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = log_mix(theta[i], lambda1[i], lambda2[i]);
  }
  return result;
}

/**
 * @brief Compute log(mean(exp(x))) for a vector in a numerically stable way
 *
 * This function computes the log of the mean of exponentials, which is
 * equivalent to log(mean(exp(x))) = log(sum(exp(x))/n) = log_sum_exp(x) - log(n)
 *
 * @param x Vector of log-scale values
 * @return log(mean(exp(x[i])))
 */
inline double log_mean_exp(const Rcpp::NumericVector& x) {
  if (x.size() == 0) return R_NegInf;
  return log_sum_exp(x) - std::log(static_cast<double>(x.size()));
}

// Integral on [0,1] of a positive unimodal function supplied on the log scale.
// Split outward from its mode so cubature cannot miss a narrow endpoint peak.
// Used for ill-conditioned closed forms and underflowing race integrals.
// The outer likelihood integrals target 1e-8 relative error. Reserve a
// smaller error budget for distribution evaluations inside those integrals.
constexpr double DISTRIBUTION_REL_TOL = 1e-10;
inline bool probability_roundoff_exceeds(double value,double magnitude) {
  return !(value>0.0) || 32*std::numeric_limits<double>::epsilon()*magnitude
                         > DISTRIBUTION_REL_TOL*value;
}

template<class F> inline double log_integral_unit(F log_value,double mode,
                                                 double rel_tol=DISTRIBUTION_REL_TOL) {
  if (mode<0.0) {
    // Only needed when an ordinary race integral underflows. The distribution
    // fallbacks supply their analytical mode directly.
    constexpr double ratio=0.3819660112501051518;
    double a=0,b=1,c=ratio,d=1-ratio;
    double fa=log_value(a),fb=log_value(b),fc=log_value(c),fd=log_value(d);
    for (int i=0;i<80;++i) {
      const double high=std::max(std::max(fa,fb),std::max(fc,fd));
      const double low=std::min(std::min(fa,fb),std::min(fc,fd));
      if (high-low<0.25 || !(a<c && c<d && d<b)) break;
      if (fc>fd) { b=d;fb=fd;d=c;fd=fc;c=a+ratio*(b-a);fc=log_value(c); }
      else { a=c;fa=fc;c=d;fc=fd;d=b-ratio*(b-a);fd=log_value(d); }
    }
    mode=fa>=fb && fa>=fc && fa>=fd ? a : (fb>=fc && fb>=fd ? b : (fc>=fd ? c : d));
  }
  const double scale=log_value(mode);
  if (!is_finite(scale)) return scale;
  struct Work { F& f; double scale; } work{log_value,scale};
  auto integrand=[](unsigned,const double* x,void* p,unsigned,double* out) -> int {
    const auto& w=*static_cast<Work*>(p);
    out[0]=std::exp(w.f(x[0])-w.scale);
    return 0;
  };
  double total=0.0;
  for (int direction : {-1,1}) {
    const double end=direction<0 ? 0.0 : 1.0;
    double step=std::fabs(end-mode), at=mode;
    if (step==0.0) continue;
    while (step>1e-15 && log_value(mode+direction*step)<scale-1.0) step*=0.5;
    while (at!=end) {
      const double next=direction<0 ? std::max(end,at-step) : std::min(end,at+step);
      double lo=std::min(at,next),hi=std::max(at,next),value=0.0,error=0.0;
      if (lo==hi) break;
      const double tolerance=std::max(rel_tol,16*std::numeric_limits<double>::epsilon()*std::fabs(scale));
      hcubature(integrand,&work,1,&lo,&hi,10000,tolerance*total/128,
                tolerance,&value,&error);
      total+=value;
      at=next;
      step*=2.0;
    }
  }
  return scale+std::log(total);
}

#endif // composite_functions_h
