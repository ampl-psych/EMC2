#ifndef exgaussian_functions_h
#define exgaussian_functions_h

// Truncated ex-Gaussian density / CDF for the stop-signal models. The plain
// ex-Gaussian (dexg / pexg / SIG_TAU_EPS) and the log-space helpers come from
// the package's own ex-Gaussian race implementation (model_exgaussian.h).

#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include "model_exgaussian.h"     // dexg, pexg, SIG_TAU_EPS, log1m, log_diff_exp
#include "composite_functions.h"  // log1m_exp, log_sum_exp, log_mix

// probability density function of truncated ex-Gaussian distribution
inline double dtexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const double lower = R_NegInf,
    const double upper = R_PosInf,
    const bool log_d = false
) {

  if (lower == R_NegInf && upper == R_PosInf) {
    return dexg(x, mu, sigma, tau, log_d);
  }
  if (sigma <= 0. || tau <= 0.) return NA_REAL;
  if (lower >= upper) return NA_REAL;
  if (x <= lower || x >= upper) return log_d ? R_NegInf : 0.;

  double x_ld = dexg(x, mu, sigma, tau, true);
  if (x_ld == R_NegInf) return log_d ? R_NegInf : 0.;

  double lower_lcdf, upper_lcdf;
  if (lower == R_NegInf) {
    lower_lcdf = R_NegInf;
  } else {
    lower_lcdf = pexg(lower, mu, sigma, tau, true, true);
  }
  if (upper == R_PosInf) {
    upper_lcdf = 0.;
  } else {
    upper_lcdf = pexg(upper, mu, sigma, tau, true, true);
  }

  if (lower_lcdf == upper_lcdf) return log_d ? R_NegInf : 0.;

  double log_normaliser;
  if (lower_lcdf == R_NegInf) {
    log_normaliser = upper_lcdf;
  } else {
    log_normaliser = log_diff_exp(upper_lcdf, lower_lcdf);
  }
  if (log_normaliser == R_NegInf) return log_d ? R_NegInf : 0.;

  double log_out = x_ld - log_normaliser;
  return log_d ? log_out : std::exp(log_out);
}

// cumulative distribution function of truncated ex-Gaussian distribution
inline double ptexg(
    const double q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const double lower = R_NegInf,
    const double upper = R_PosInf,
    const bool lower_tail = true,
    const bool log_p = false
) {

  if (lower == R_NegInf && upper == R_PosInf) {
    return pexg(q, mu, sigma, tau, lower_tail, log_p);
  }
  if (sigma <= 0. || tau <= 0.) return NA_REAL;
  if (lower >= upper) return NA_REAL;
  if (q <= lower) {
    double out = lower_tail ? 0. : 1.;
    return log_p ? (out == 0. ? R_NegInf : 0.) : out;
  }
  if (q >= upper) {
    double out = lower_tail ? 1. : 0.;
    return log_p ? (out == 0. ? R_NegInf : 0.) : out;
  }

  double q_cdf = pexg(q, mu, sigma, tau);

  double lower_cdf, upper_cdf;
  if (lower == R_NegInf) {
    lower_cdf = 0.;
  } else {
    lower_cdf = pexg(lower, mu, sigma, tau);
  }
  if (upper == R_PosInf) {
    upper_cdf = 1.;
  } else {
    upper_cdf = pexg(upper, mu, sigma, tau);
  }

  double normaliser = upper_cdf - lower_cdf;
  if (normaliser <= 0.) return NA_REAL;

  double out;
  if (lower_tail) {
    out = (q_cdf - lower_cdf) / normaliser;
  } else {
    out = (upper_cdf - q_cdf) / normaliser;
  }
  out = std::max(0., std::min(1., out));

  return log_p ? std::log(out) : out;
}


#endif
