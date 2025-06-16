#ifndef phi_functions_h
#define phi_functions_h

#include <Rcpp.h>
#include <cmath>
#include <limits>

const double LOG_STD_NORMAL_CONST = -0.91893853320467274178032973640562;
const double LOG_EPSILON = std::log(std::numeric_limits<double>::epsilon());

/**
 * Inverse logit function with numerical stability
 * @param x Input value
 * @return 1 / (1 + exp(-x))
 */
inline double inv_logit(double x) {
  if (x < 0.0) {
    double exp_x = std::exp(x);
    if (x < LOG_EPSILON) {
      return exp_x;
    } else {
      return exp_x / (1.0 + exp_x);
    }
  } else {
    return 1.0 / (1.0 + std::exp(-x));
  }
}

/**
 * Fast approximation of the standard normal CDF using logistic approximation
 * Based on Bowling et al. (2009) "A Logistic Approximation to the Cumulative Normal Distribution"
 * Journal of Industrial Engineering and Management 2(1): 114-27
 * @param x Input value
 * @return Approximate Phi(x)
 */
inline double phi_approx(double x) {
  return inv_logit(0.07056 * std::pow(x, 3.0) + 1.5976 * x);
}

/**
 * Safe computation of the standard normal CDF with Mills ratio approximation for extreme values
 * @param z Input value
 * @param lower_tail If true, return P(Z <= z), otherwise P(Z > z)
 * @param log_p If true, return log probability
 * @param lower_thresh Threshold below which to use Mills ratio approximation
 * @param upper_thresh Threshold above which to use Mills ratio approximation
 * @return Phi(z) or log(Phi(z))
 */
double phi_safe(
    double z,
    bool lower_tail = true,
    bool log_p = false,
    double lower_thresh = -15.0,
    double upper_thresh = 8.3
) {
  double log_out;

  if (z < lower_thresh) {
    if (lower_tail) {
      log_out = LOG_STD_NORMAL_CONST - 0.5 * z * z - std::log(-z);
    } else {
      log_out = 0.0;
    }
  } else if (z > upper_thresh) {
    if (lower_tail) {
      log_out = 0.0;
    } else {
      log_out = LOG_STD_NORMAL_CONST - 0.5 * z * z - std::log(z);
    }
  } else {
    log_out = R::pnorm(z, 0.0, 1.0, lower_tail, true);
  }

  return log_p ? log_out : std::exp(log_out);
}


#endif // phi_functions_h
