#ifndef exgaussian_functions_h
#define exgaussian_functions_h

#include <Rcpp.h>
#include <cmath>
#include "composite_functions.h"
#include "pnorm_utils.h"
#include "r_constants.h"
#include "nan_check.h"
using namespace Rcpp;

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

static constexpr double SIG_TAU_EPS = 1e-12;

// probability density function of ex-Gaussian distribution
inline double dexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool log_d = false,
    const double normal_ratio = 0.,
    const double shifted_exp_sigma = 0.
) {

  // minimal parameter validation - assuming all inputs are finite
  if (sigma <= 0. || tau <= 0.) return na_real();

  if (shifted_exp_sigma > 0. && sigma < shifted_exp_sigma) {
    return R::dexp(x - mu, tau, log_d);
  }

  if (normal_ratio > 0. && tau < normal_ratio * sigma) {
    return R::dnorm(x, mu, sigma, log_d);
  }

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  // Numerically stable branch for extreme tails where z = (x-mu)/sigma - sigma/tau is very negative.
  // In this regime, naive log formulation suffers catastrophic cancellation between
  // +sigma^2/(2 tau^2) and log Phi(z) ≈ -(sigma/tau)^2/2 + ...
  // Using Mills ratio, one can simplify to a stable expression. For a = -z:
  // Phi(z) ≈ phi(z) * (1/a - 1/a^3 + 3/a^5 - 15/a^7 + 105/a^9).
  // The leading 1/a term gives the original first-order approximation; the
  // remaining terms substantially improve accuracy near the z < -8 cutoff.
  {
    double y = (x - mu) / sig_p;
    double a = sig_p / tau_p;
    double z = y - a;
    if (z < -8.0) {
      double mills_a = -z;
      double inv_a2 = 1.0 / (mills_a * mills_a);
      double mills_correction = 1.0 - inv_a2 + 3.0 * inv_a2 * inv_a2
                                - 15.0 * inv_a2 * inv_a2 * inv_a2
                                + 105.0 * inv_a2 * inv_a2 * inv_a2 * inv_a2;
      double log_out_stable = -std::log(tau_p) - 0.5*std::log(2.0 * M_PI)
                              - 0.5 * y * y - std::log(mills_a)
                              + std::log(mills_correction);
      return log_d ? log_out_stable : std::exp(log_out_stable);
    }
  }

  // compute Phi term
  double z = (x - mu) / sig_p - sig_p / tau_p;
  double log_phi = PNORM_STD(z, true, true);
  if (is_nan(log_phi)) {
    return na_real();
  }
  if (is_inf(log_phi)) {
    return log_d ? log_phi : std::exp(log_phi);
  }

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
    const bool log_p = false,
    const double normal_ratio = 0.,
    const double shifted_exp_sigma = 0.
) {

  if (is_nan(q) || is_nan(mu) || is_nan(sigma) || is_nan(tau)) return na_real();
  if (sigma <= 0. || tau <= 0.) return na_real();

  if (shifted_exp_sigma > 0. && sigma < shifted_exp_sigma) {
    return R::pexp(q - mu, tau, lower_tail, log_p);
  }

  if (normal_ratio > 0. && tau < normal_ratio * sigma) {
    return PNORM_STD((q - mu) / sigma, lower_tail, log_p);
  }

  // handle infinite q
  if (is_inf(q)) {
    double out = (q < 0.) ? 0. : 1.;
    if (!lower_tail) out = 1. - out;
    return log_p ? std::log(out) : out;
  }

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  const double z = (q - mu) / sig_p;
  const double a = sig_p / tau_p;
  const double w = z - a;
  double log_second_term;
  if (w < -40.0) {
    // Scale Phi(w) analytically before combining with its large exponential.
    // This CDF-only continued fraction does not change dexg's Mills polynomial.
    const double v = -w;
    double denom = v;
    for (int k = 12; k >= 1; --k) denom = v + k / denom;
    log_second_term = -(0.5 * z) * z - 0.5 * std::log(2.0 * M_PI) - std::log(denom);
  } else {
    log_second_term = a * (0.5 * a - z) + PNORM_STD(w, true, true);
  }
  if (lower_tail) {
    const double log_phi = PNORM_STD(z, true, true);
    const double out = log_phi > log_second_term ? log_diff_exp(log_phi, log_second_term) : neg_inf();
    return log_p ? out : std::exp(out);
  }
  // S(q) = Phi(-z) + exp(-a*z + a*a/2) * Phi(z-a), without 1-F(q).
  const double log_normal_survivor = PNORM_STD(z, false, true);
  double out;
  if (is_neg_inf(log_normal_survivor)) out = log_second_term;
  else if (is_neg_inf(log_second_term)) out = log_normal_survivor;
  else {
    const double hi = std::max(log_normal_survivor, log_second_term);
    const double lo = std::min(log_normal_survivor, log_second_term);
    out = hi + std::log1p(std::exp(lo - hi));
  }
  if (is_nan(out)) return na_real();
  out = std::min(0.0, out);
  return log_p ? out : std::exp(out);
}

// Shared interval mass for truncated densities, CDFs and survivors.
inline double exg_log_interval(double lower, double upper, double mu, double sigma,
                               double tau, double normal_ratio, double shifted_exp_sigma) {
  if (is_nan(lower) || is_nan(upper) || is_nan(mu) || is_nan(sigma) || is_nan(tau)) return na_real();
  if (sigma <= 0. || tau <= 0. || lower > upper) return na_real();
  if (lower == upper) return neg_inf();
  if (shifted_exp_sigma > 0. && sigma < shifted_exp_sigma) {
    const double lo = std::max(lower, mu);
    if (upper <= lo) return neg_inf();
    return -(lo - mu) / tau + R::pexp(upper - lo, tau, true, true);
  }
  if (is_neg_inf(lower)) return pexg(upper, mu, sigma, tau, true, true, normal_ratio, shifted_exp_sigma);
  if (is_pos_inf(upper)) return pexg(lower, mu, sigma, tau, false, true, normal_ratio, shifted_exp_sigma);

  const double width = upper - lower;
  const double mid = lower + 0.5 * width;
  const bool normal = normal_ratio > 0. && tau < normal_ratio * sigma;
  const double sig = normal ? sigma : std::max(sigma, SIG_TAU_EPS);
  if (width / sig * std::max(1.0, std::fabs((mid - mu) / sig)) < 1e-2) {
    // Simpson's rule in log space avoids cancellation at nearly equal endpoints.
    const double l = dexg(lower, mu, sigma, tau, true, normal_ratio, shifted_exp_sigma);
    const double m = dexg(mid, mu, sigma, tau, true, normal_ratio, shifted_exp_sigma);
    const double u = dexg(upper, mu, sigma, tau, true, normal_ratio, shifted_exp_sigma);
    if (is_nan(l) || is_nan(m) || is_nan(u)) return na_real();
    const double peak = std::max(l, std::max(m, u));
    if (is_neg_inf(peak)) return neg_inf();
    return std::log(width) + peak + std::log((std::exp(l-peak)+4.0*std::exp(m-peak)+std::exp(u-peak))/6.0);
  }
  const double uf = pexg(upper, mu, sigma, tau, true, true, normal_ratio, shifted_exp_sigma);
  if (uf <= -composite_detail::LOG2) {
    return log_diff_exp(uf, pexg(lower, mu, sigma, tau, true, true, normal_ratio, shifted_exp_sigma));
  }
  return log_diff_exp(pexg(lower, mu, sigma, tau, false, true, normal_ratio, shifted_exp_sigma),
                      pexg(upper, mu, sigma, tau, false, true, normal_ratio, shifted_exp_sigma));
}

// probability density function of truncated ex-Gaussian distribution
inline double dtexg(
    const double x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const double lower = neg_inf(),
    const double upper = pos_inf(),
    const bool log_d = false,
    const double normal_ratio = 0.,
    const double shifted_exp_sigma = 0.
) {

  if (is_neg_inf(lower) && is_pos_inf(upper)) {
    return dexg(x, mu, sigma, tau, log_d, normal_ratio, shifted_exp_sigma);
  }
  if (is_nan(x) || is_nan(lower) || is_nan(upper)) return na_real();
  if (sigma <= 0. || tau <= 0.) return na_real();
  if (lower >= upper) return na_real();
  if (x <= lower || x >= upper) return log_d ? neg_inf() : 0.;

  double x_ld = dexg(x, mu, sigma, tau, true, normal_ratio, shifted_exp_sigma);
  if (is_neg_inf(x_ld)) return log_d ? neg_inf() : 0.;

  const double log_normaliser = exg_log_interval(lower, upper, mu, sigma, tau, normal_ratio, shifted_exp_sigma);
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
    const bool log_p = false,
    const double normal_ratio = 0.,
    const double shifted_exp_sigma = 0.
) {

  if (is_neg_inf(lower) && is_pos_inf(upper)) {
    return pexg(q, mu, sigma, tau, lower_tail, log_p, normal_ratio, shifted_exp_sigma);
  }
  if (is_nan(q) || is_nan(lower) || is_nan(upper)) return na_real();
  if (sigma <= 0. || tau <= 0.) return na_real();
  if (lower >= upper) return na_real();
  if (q <= lower) {
    double out = lower_tail ? 0. : 1.;
    return log_p ? (out == 0. ? neg_inf() : 0.) : out;
  }
  if (q >= upper) {
    double out = lower_tail ? 1. : 0.;
    return log_p ? (out == 0. ? neg_inf() : 0.) : out;
  }

  const double norm = exg_log_interval(lower, upper, mu, sigma, tau, normal_ratio, shifted_exp_sigma);
  if (is_nan(norm) || is_neg_inf(norm)) return na_real();
  const double numerator = lower_tail ? exg_log_interval(lower, q, mu, sigma, tau, normal_ratio, shifted_exp_sigma)
                                     : exg_log_interval(q, upper, mu, sigma, tau, normal_ratio, shifted_exp_sigma);
  double out = numerator - norm;
  if (is_nan(out)) return na_real();
  // Derive the larger tail from the smaller interval to keep complements consistent.
  if (out > -composite_detail::LOG2) {
    const double other = lower_tail ? exg_log_interval(q, upper, mu, sigma, tau, normal_ratio, shifted_exp_sigma) - norm
                                   : exg_log_interval(lower, q, mu, sigma, tau, normal_ratio, shifted_exp_sigma) - norm;
    if (is_nan(other)) return na_real();
    out = log1m_exp(std::min(0.0, other));
  }
  out = std::min(0.0, out);
  return log_p ? out : std::exp(out);
}


#endif
