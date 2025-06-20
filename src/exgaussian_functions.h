#ifndef exgaussian_functions_h
#define exgaussian_functions_h

#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include "composite_functions.h"

const double SIG_TAU_EPS = 1e-12;

// probability density function of ex-Gaussian distribution
double dexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool log_d = false
) {

  // minimal parameter validation - assuming all inputs are finite
  if (sigma <= 0. || tau <= 0.) return NA_REAL;

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  // compute Phi term
  double z = (x - mu) / sig_p - sig_p / tau_p;
  double log_phi = R::pnorm(z, 0.0, 1.0, true, true);
  if (std::isnan(log_phi)) {
    return NA_REAL;
  }
  if (std::isinf(log_phi)) {
    return log_d ? log_phi : std::exp(log_phi);
  }

  // compute exp term
  double log_exp = (mu - x) / tau_p + (sig_p * sig_p) / (2. * tau_p * tau_p);

  // final output: log density of ex-Gaussian
  double log_out = -std::log(tau_p) + log_exp + log_phi;

  return log_d ? log_out : std::exp(log_out);
}

// cumulative distribution function of ex-Gaussian distribution
double pexg(
    const double q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool lower_tail = true,
    const bool log_p = false
) {

  // minimal parameter validation - assuming mu, sigma, and tau are all finite
  if (sigma <= 0. || tau <= 0.) return NA_REAL;

  // handle infinite q
  if (std::isinf(q)) {
    double out = (q < 0.) ? 0. : 1.;
    if (!lower_tail) out = 1. - out;
    return log_p ? std::log(out) : out;
  }

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  // compute the two Phi terms
  double log_phi_1 = R::pnorm((q - mu) / sig_p, 0.0, 1.0, true, true);
  double log_phi_2 = R::pnorm(
    (q - mu) / sig_p - sig_p / tau_p, 0.0, 1.0, true, true
  );

  // compute the exp term in log space
  double log_exp_term = (mu - q) / tau_p + (sig_p * sig_p) / (2. * tau_p * tau_p);

  // combined second term
  double log_second_term = log_exp_term + log_phi_2;

  // now obtain ex-Gaussian log CDF
  double log_cdf_lower;
  if (log_phi_1 > log_second_term) {
    log_cdf_lower = log_diff_exp(log_phi_1, log_second_term);
  } else {
    log_cdf_lower = R_NegInf;
  }

  double out;
  if (lower_tail) {
    out = log_p ? log_cdf_lower : std::exp(log_cdf_lower);
  } else {
    if (log_cdf_lower == R_NegInf) {
      out = log_p ? 0. : 1.;
    } else {
      double cdf_lower = std::exp(log_cdf_lower);
      if (cdf_lower >= 1. - 1e-15) {
        out = log_p ? R_NegInf : 0.;
      } else {
        out = log_p ? log1m(cdf_lower) : -std::expm1(log_cdf_lower);
      }
    }
  }

  return(out);
}

// probability density function of truncated ex-Gaussian distribution
double dtexg(
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
  // CHECK this might be too harsh; maybe x < lower || x > upper instead
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
double ptexg(
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
  // CHECK this might be too harsh; maybe q < lower, q > upper instead
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
