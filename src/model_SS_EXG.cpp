#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
#include "race_integrate.h"


// [[Rcpp::export]]
double dEXG(
    double x,
    double mu = 5.,
    double sigma = 1.,
    double tau = 1.,
    bool log_d = false
) {

  if (sigma <= 0 || tau <= 0) return NA_REAL;

  // if Gaussian SD is practically zero, treat as shifted exponential
  if (sigma < 1e-4) {
    return R::dexp(x - mu, tau, log_d);
  }

  // if exponential rate has very small value relative to Gaussian SD, treat as
  // Gaussian
  if (tau < .05 * sigma) {
    return R::dnorm(x, mu, sigma, log_d);
  }

  // main ex-Gaussian computation:
  double sig2 = sigma * sigma;
  double z = x - mu - sig2 / tau;
  double log_phi = std::log(R::pnorm(z / sigma, 0., 1., true, false));
  double log_exp = -std::log(tau) - (z + sig2 / (2. * tau)) / tau;
  double out = log_exp + log_phi;

  return log_d ? out : std::exp(out);
}


// [[Rcpp::export]]
double pEXG(
    double q,
    double mu = 5.,
    double sigma = 1.,
    double tau = 1.,
    bool lower_tail = true,
    bool log_p = false
) {

  if (sigma <= 0 || tau <= 0) return NA_REAL;

  // if Gaussian SD is practically zero, treat as shifted exponential
  if (sigma < 1e-4) {
    return R::pexp(q - mu, tau, lower_tail, log_p);
  }

  // if exponential rate has very small value relative to Gaussian SD, treat as
  // Gaussian
  if (tau < .05 * sigma) {
    return R::pnorm(q, mu, sigma, lower_tail, log_p);
  }

  // main ex-Gaussian computation:
  double out;
  if (traits::is_infinite<REALSXP>(q)) {
    out = (q < 0.) ? 0. : 1.;
  } else {
    double mu2 = mu * mu;
    double sig2 = sigma * sigma;
    double sig2_tau = sig2 / tau;
    double z = q - mu - sig2_tau;
    double norm1 = R::pnorm((q - mu) / sigma, 0., 1., true, false);
    double norm2 = R::pnorm(z / sigma, 0., 1., true, false);
    double exponent = (
      std::pow(mu + sig2_tau, 2.) - mu2 - 2. * q * sig2_tau
    ) / (2. * sig2);
    out = norm1 - std::exp(std::log(norm2) + exponent);
  }

  if (!lower_tail) out = 1. - out;

  return log_p ? std::log(out) : out;
}


// [[Rcpp::export]]
NumericVector exg_lpdf(
    NumericVector rts,
    NumericMatrix pars,
    LogicalVector idx,
    double min_ll,
    LogicalVector is_ok
) {
  // pars columns: mu=0, sigma=1, tau=2, muS=3, sigmaS=4, tauS=5, tf=6, gf=7, SSD=8
  int n_out = sum(idx);
  NumericVector out(n_out);
  int k = 0;

  for (int i = 0; i < rts.size(); i++) {
    if (!idx[i]) continue;

    if (NumericVector::is_na(pars(i, 0))) {
      out[k] = R_NegInf;
    } else if (is_ok[i] == TRUE) {
      out[k] = dEXG(rts[i], pars(i, 0), pars(i, 1), pars(i, 2), true);
    } else{
      out[k] = min_ll;
    }
    k++;
  }

  return(out);
}


// [[Rcpp::export]]
NumericVector exg_lccdf(
    NumericVector rts,
    NumericMatrix pars,
    LogicalVector idx,
    double min_ll,
    LogicalVector is_ok
) {
  // pars columns: mu=0, sigma=1, tau=2, muS=3, sigmaS=4, tauS=5, tf=6, gf=7, SSD=8
  int n_out = sum(idx);
  NumericVector out(n_out);
  int k = 0;

  for (int i = 0; i < rts.size(); i++) {
    if (!idx[i]) continue;

    if (NumericVector::is_na(pars(i, 0))) {
      out[k] = R_NegInf;
    } else if (is_ok[i] == TRUE) {
      out[k] = pEXG(rts[i], pars(i, 0), pars(i, 1), pars(i, 2), false, true);
    } else{
      out[k] = min_ll;
    }
    k++;
  }

  return(out);
}


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
