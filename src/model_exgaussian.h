#ifndef exgaussian_h
#define exgaussian_h

#include <Rcpp.h>
#include "RaceSpec.h"
#include "utility_functions.h"
#include "math_utils.h"  // must be before Rcpp
#include "pnorm_utils.h"
#include "ParamTable.h"
#include "hcubature.h"

using namespace Rcpp;

inline constexpr double SIG_TAU_EPS = 1e-12;
inline constexpr double LOG_SQRT_2PI = 0.91893853320467274178;

// Laplace continued fraction denominator for the Mills ratio.
// For z > 0, log Phi(-z) = log phi(z) - log(mills_cf_denom(z)).
// Used by both fast_log_upper_tail and dexg's tail branch.
inline double mills_cf_denom(double z) {
  return z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
}

inline double log1m(double x) {
  return std::log1p(-x);
}

inline double log_diff_exp(double a, double b) {
  double diff = b - a;

  // For numerical stability when a and b are close (mirrors log1m_exp branching on diff).
  // When diff is close to 0, exp(diff) ≈ 1 and log1m(exp(diff)) suffers catastrophic
  // cancellation; expm1 avoids this. When diff is far from 0, exp(diff) is small and
  // log1m(exp(diff)) is already stable.
  if (diff > -0.693147) {  // diff > -log(2), i.e., exp(diff) close to 1 → use expm1
    return a + std::log(-std::expm1(diff));
  } else {                 // exp(diff) < 0.5 → log1m(exp(diff)) is stable
    return a + log1m(std::exp(diff));
  }
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

void dexg_pexg_fast(const double*           rt,
                    const ParamTable&       pt,
                    const RaceSpec&         spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__    ll_row,
                    RaceScratch&            scratch);

void exg_survivor(const std::vector<int>&    idx,
                  const std::vector<double>& bound,
                  const ParamTable&          pt,
                  const RaceSpec&            spec,
                  double* __restrict__       out,
                  RaceScratch&               scratch);


void exg_survivor_with_response(const std::vector<int>&    idx,
                                const std::vector<int>&    winner,
                                const std::vector<double>& lower,
                                const std::vector<double>& upper,
                                int                        n_acc,
                                const ParamTable&          pt,
                                const RaceSpec&            spec,
                                double* __restrict__       out);

#endif
