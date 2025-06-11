#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
#include "race_integrate.h"

const double SIG_TAU_EPS = 1e-12;
const double LOG_STD_NORMAL_CONST = -0.91893853320467274178032973640562;

double phi_safe(
    double z,
    bool lower_tail = true,
    bool log_p = false,
    double lower_thresh = -15.,
    double upper_thresh = 8.3
) {
  double log_out;

  // use Mills ratio approximation for extreme values, otherwise R's pnorm
  if (z < lower_thresh) {
    if (lower_tail) {
      log_out = LOG_STD_NORMAL_CONST - 0.5 * z * z - std::log(-z);
    } else {
      log_out = 0.;
    }
  } else if (z > upper_thresh) {
    if (lower_tail) {
      log_out = 0.;
    } else {
      log_out = LOG_STD_NORMAL_CONST - 0.5 * z * z - std::log(z);
    }
  } else {
    log_out = R::pnorm(z, 0., 1., lower_tail, true);
  }

  return log_p ? log_out : std::exp(log_out);
}

double log_diff_exp(double a, double b) {
  if (a <= b) return R_NegInf;
  if (b == R_NegInf) return a;
  return a + std::log1p(-std::exp(b - a));
}

double log_sum_exp(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  if (a > b) {
    return a + std::log1p(std::exp(b - a));
  } else {
    return b + std::log1p(std::exp(a - b));
  }
}

// [[Rcpp::export]]
double dexg(
    double x,
    double mu = 5.,
    double sigma = 1.,
    double tau = 1.,
    bool log_d = false
) {

  // minimal parameter validation - assuming all inputs are finite
  if (sigma <= 0. || tau <= 0.) return NA_REAL;

  // protect against numerical issues due to extremely small sigma or tau values
  double tau_p = std::max(tau, SIG_TAU_EPS);
  double sig_p = std::max(sigma, SIG_TAU_EPS);

  // compute Phi term
  double z = (x - mu) / sig_p - sig_p / tau_p;
  double log_phi = phi_safe(z, true, true);
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


// [[Rcpp::export]]
double pexg(
    double q,
    double mu = 5.,
    double sigma = 1.,
    double tau = 1.,
    bool lower_tail = true,
    bool log_p = false
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
  double log_phi_1 = phi_safe((q - mu) / sig_p, true, true);
  double log_phi_2 = phi_safe((q - mu) / sig_p - sig_p / tau_p, true, true);

  // compute the exp term in log space
  double log_exp_term = (mu - q) / tau_p + (sig_p * sig_p) / (2. * tau_p * tau_p);

  // combined second term
  double log_second_term = log_exp_term + log_phi_2;

  // final output: ex-Gaussian CDF
  double log_out;
  if (log_phi_1 > log_second_term) {
    log_out = log_diff_exp(log_phi_1, log_second_term);
  } else {
    log_out = R_NegInf;
  }

  if (!lower_tail) {
    double temp_out = std::exp(log_out);
    if (log_out == R_NegInf || temp_out < 1e-15) {
      log_out = 0.;
    } else if (log_out == 0. || temp_out > 1. - 1e-15) {
      log_out = R_NegInf;
    } else {
      log_out = std::log1p(-temp_out);
    }
  }

  return log_p ? log_out : std::exp(log_out);
}

// [[Rcpp::export]]
double dtexg(
    double x,
    double mu = 5.,
    double sigma = 1.,
    double tau = 1.,
    double lower = R_NegInf,
    double upper = R_PosInf,
    bool log_d = false
) {

  if (lower >= upper) return NA_REAL;
  if (x <= lower || x >= upper) return log_d ? R_NegInf : 0.;
  if (lower == R_NegInf && upper == R_PosInf) {
    return dexg(x, mu, sigma, tau, log_d);
  }

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

// [[Rcpp::export]]
double ptexg(
    double q,
    double mu = 5.,
    double sigma = 1.,
    double tau = 1.,
    double lower = R_NegInf,
    double upper = R_PosInf,
    bool lower_tail = true,
    bool log_p = false
) {

  if (lower >= upper) return NA_REAL;
  if (q <= lower) {
    double out = lower_tail ? 0. : 1.;
    return log_p ? (out == 0. ? R_NegInf : 0.) : out;
  }
  if (q >= upper) {
    double out = lower_tail ? 1. : 0.;
    return log_p ? (out == 0. ? R_NegInf : 0.) : out;
  }
  if (lower == R_NegInf && upper == R_PosInf) {
    return pexg(q, mu, sigma, tau, lower_tail, log_p);
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
  out = std::clamp(out, 0., 1.);

  return log_p ? std::log(out) : out;
}




// // [[Rcpp::export]]
// NumericVector exg_lpdf(
//     NumericVector rt,
//     NumericMatrix pars,
//     LogicalVector idx,
//     double min_ll
// ) {
//   // pars columns: mu=0, sigma=1, tau=2, muS=3, sigmaS=4, tauS=5, tf=6, gf=7, SSD=8
//   int n_out = sum(idx);
//   NumericVector out(n_out);
//   int k = 0;
//
//   for (int i = 0; i < rt.size(); i++) {
//     if (!idx[i]) continue;
//
//     double log_d = dexg(rt[i], pars(i, 0), pars(i, 1), pars(i, 2), true);
//     if (!std::isfinite(log_d)) {
//       out[k] = min_ll;
//     } else {
//       out[k] = log_d;
//     }
//     k++;
//   }
//
//   return(out);
// }
//
//
// // [[Rcpp::export]]
// NumericVector exg_lccdf(
//     NumericVector rt,
//     NumericMatrix pars,
//     LogicalVector idx,
//     double min_ll
// ) {
//   // pars columns: mu=0, sigma=1, tau=2, muS=3, sigmaS=4, tauS=5, tf=6, gf=7, SSD=8
//   int n_out = sum(idx);
//   NumericVector out(n_out);
//   int k = 0;
//
//   for (int i = 0; i < rt.size(); i++) {
//     if (!idx[i]) continue;
//
//     double log_s = pexg(rt[i], pars(i, 0), pars(i, 1), pars(i, 2), false, true);
//     if (!std::isfinite(log_s)) {
//       out[k] = min_ll;
//     } else {
//       out[k] = log_s;
//     }
//     k++;
//   }
//
//   return(out);
// }


// NumericVector exg_go_race(
//     NumericVector rts,
//     NumericMatrix pars,
//     LogicalVector idx,
//     double min_ll,
//     LogicalVector is_ok
// ) {
// }
//
// NumericVector exg_stopfail_race(
//     NumericVector rts,
//     NumericMatrix pars,
//     LogicalVector idx,
//     double min_ll,
//     LogicalVector is_ok
// ) {
// }
//
// NumericVector exg_stopsuccess_race(
//     NumericMatrix pars,
//     LogicalVector idx,
//     double min_ll,
//     LogicalVector is_ok
// ) {
// }



// // [[Rcpp::export]]
// NumericVector dEXGrace(NumericMatrix dt,
//                        NumericVector mu, NumericVector sigma, NumericVector tau){
//   int n = mu.size();
//   NumericVector out(dt.nrow());
//   out = dEXG(dt(0, _), mu[0], sigma[0], tau[0], false);
//   for (int i = 1; i < n; i++){
//     out = out * pEXG(dt(i, _), mu[i], sigma[i], tau[i], false, false);
//   }
//   return out;
// }
//
// // [[Rcpp::export]]
// NumericVector stopfn_exg(NumericVector t,
//                          NumericVector mu, NumericVector sigma, NumericVector tau,
//                          double SSD){
//   NumericVector tmp(mu.size() * t.size());
//   tmp = rep_each(t, mu.size()) + SSD;
//   NumericMatrix dt(mu.size(), t.size(), tmp.begin());
//   dt(0, _) = dt(0, _) - SSD;
//   return dEXGrace(dt, mu, sigma, tau);
// }
