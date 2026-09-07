#ifndef exgaussian_functions_h
#define exgaussian_functions_h

#include <Rcpp.h>
#include <cmath>
#include "composite_functions.h"  // log1m, log1m_exp, log_sum_exp, log_diff_exp, log_mix
#include "pnorm_utils.h"
#include "r_constants.h"
#include "nan_check.h"
using namespace Rcpp;

inline constexpr double SIG_TAU_EPS = 1e-12;
inline constexpr double LOG_SQRT_2PI = 0.91893853320467274178;

// Laplace continued fraction denominator for the Mills ratio.
// For z > 0, log Phi(-z) = log phi(z) - log(mills_cf_denom(z)).
// Used by both fast_log_upper_tail and dexg's tail branch.
inline double mills_cf_denom(double z) {
  return z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
}

// probability density function of ex-Gaussian distribution
inline double dexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool log_d = false
) {

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  // Numerically stable branch for extreme tails where z = (x-mu)/sigma - sigma/tau is very negative.
  {
    double y = (x - mu) / sig_p;
    double a = sig_p / tau_p;
    double z = y - a;
    if (z < -8.0) {
      double log_out_stable = -std::log(tau_p) - LOG_SQRT_2PI
      - 0.5 * y * y - std::log(mills_cf_denom(-z));
      return log_d ? log_out_stable : std::exp(log_out_stable);
    }
  }

  // compute Phi term
  double z = (x - mu) / sig_p - sig_p / tau_p;
  double log_phi = PNORM_STD(z, true, true);

  // compute exp term
  double log_exp = (mu - x) / tau_p + (sig_p * sig_p) / (2. * tau_p * tau_p);

  // final output: log density of ex-Gaussian
  double log_out = -std::log(tau_p) + log_exp + log_phi;

  return log_d ? log_out : std::exp(log_out);
}

// cumulative distribution function of ex-Gaussian distribution
inline double pexg(
    const double q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool lower_tail = true,
    const bool log_p = false
) {

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  // compute the two Phi terms
  double log_phi_1 = PNORM_STD((q - mu) / sig_p, true, true);
  double log_phi_2 = PNORM_STD((q - mu) / sig_p - sig_p / tau_p, true, true);

  // compute the exp term in log space
  double log_exp_term = (mu - q) / tau_p + (sig_p * sig_p) / (2. * tau_p * tau_p);

  // combined second term
  double log_second_term = log_exp_term + log_phi_2;

  // now obtain ex-Gaussian log CDF
  double log_cdf_lower;
  if (log_phi_1 > log_second_term) {
    log_cdf_lower = log_diff_exp(log_phi_1, log_second_term);
  } else {
    log_cdf_lower = neg_inf();
  }

  double out;
  if (lower_tail) {
    out = log_p ? log_cdf_lower : std::exp(log_cdf_lower);
  } else {
    if (is_neg_inf(log_cdf_lower)) {
      out = log_p ? 0. : 1.;
    } else {
      double cdf_lower = std::exp(log_cdf_lower);
      if (cdf_lower >= 1. - 1e-15) {
        out = log_p ? neg_inf() : 0.;
      } else {
        out = log_p ? log1m(cdf_lower) : -std::expm1(log_cdf_lower);
      }
    }
  }

  return(out);
}


// probability density function of truncated ex-Gaussian distribution
inline double dtexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const double lower = neg_inf(),
    const double upper = pos_inf(),
    const bool log_d = false
) {

  if (is_neg_inf(lower) && is_pos_inf(upper)) {
    return dexg(x, mu, sigma, tau, log_d);
  }
  if (sigma <= 0. || tau <= 0. || lower >= upper) return na_real();
  if (x <= lower || x >= upper) return log_d ? neg_inf() : 0.;

  double x_ld = dexg(x, mu, sigma, tau, true);
  if (is_neg_inf(x_ld)) return log_d ? neg_inf() : 0.;

  double lower_lcdf, upper_lcdf;
  if (is_neg_inf(lower)) {
    lower_lcdf = neg_inf();
  } else {
    lower_lcdf = pexg(lower, mu, sigma, tau, true, true);
  }
  if (is_pos_inf(upper)) {
    upper_lcdf = 0.;
  } else {
    upper_lcdf = pexg(upper, mu, sigma, tau, true, true);
  }

  if (lower_lcdf == upper_lcdf) return log_d ? neg_inf() : 0.;

  double log_normaliser;
  if (is_neg_inf(lower_lcdf)) {
    log_normaliser = upper_lcdf;
  } else {
    log_normaliser = log_diff_exp(upper_lcdf, lower_lcdf);
  }
  if (is_neg_inf(log_normaliser)) return log_d ? neg_inf() : 0.;

  double log_out = x_ld - log_normaliser;
  return log_d ? log_out : std::exp(log_out);
}

// cumulative distribution function of truncated ex-Gaussian distribution
inline double ptexg(
    const double q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const double lower = neg_inf(),
    const double upper = pos_inf(),
    const bool lower_tail = true,
    const bool log_p = false
) {

  if (is_neg_inf(lower) && is_pos_inf(upper)) {
    return pexg(q, mu, sigma, tau, lower_tail, log_p);
  }
  if (sigma <= 0. || tau <= 0. || lower >= upper) return na_real();
  if (q <= lower) {
    double out = lower_tail ? 0. : 1.;
    return log_p ? (out == 0. ? neg_inf() : 0.) : out;
  }
  if (q >= upper) {
    double out = lower_tail ? 1. : 0.;
    return log_p ? (out == 0. ? neg_inf() : 0.) : out;
  }

  double q_cdf = pexg(q, mu, sigma, tau);

  double lower_cdf, upper_cdf;
  if (is_neg_inf(lower)) {
    lower_cdf = 0.;
  } else {
    lower_cdf = pexg(lower, mu, sigma, tau);
  }
  if (is_pos_inf(upper)) {
    upper_cdf = 1.;
  } else {
    upper_cdf = pexg(upper, mu, sigma, tau);
  }

  double normaliser = upper_cdf - lower_cdf;
  if (normaliser <= 0.) return na_real();

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
