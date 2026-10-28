#ifndef ss_exg_h
#define ss_exg_h

#include <cmath>
#include <string>
#include <vector>
#include <Rcpp.h>
#include "utility_functions.h"
#include "exgaussian_functions.h"
#include "pnorm_utils.h"
using namespace Rcpp;

NumericVector dexg_c(
    const NumericVector x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool log_d = false
) {
  int n = x.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = dexg(x[i], mu, sigma, tau, log_d);
  }
  return(out);
}

NumericVector pexg_c(
    const NumericVector q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool lower_tail = true,
    const bool log_p = false
) {
  int n = q.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = pexg(q[i], mu, sigma, tau, lower_tail, log_p);
  }
  return(out);
}

NumericVector protect_finite(const NumericVector& x, double min_ll) {
  NumericVector out = clone(x);
  for (int i = 0; i < out.size(); i++) {
    if (!R_FINITE(out[i])) {
      out[i] = min_ll;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector dEXGrace(
    NumericMatrix dt, NumericVector mu, NumericVector sigma, NumericVector tau,
    double min_ll
){
  int n = mu.size();
  NumericVector log_out(dt.nrow());
  log_out = protect_finite(
    dexg_c(dt(0, _), mu[0], sigma[0], tau[0], true),
    min_ll
  );
  for (int i = 1; i < n; i++) {
    log_out += protect_finite(
      pexg_c(dt(i, _), mu[i], sigma[i], tau[i], false, true),
      min_ll
    );
  }
  return exp(log_out);
}

// [[Rcpp::export]]
NumericVector stopfn_exg(
    NumericVector t, NumericVector mu, NumericVector sigma, NumericVector tau,
    double SSD, double min_ll
){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dEXGrace(dt, mu, sigma, tau, min_ll);
}

// old

// [[Rcpp::export]]
NumericVector pEXG_old(NumericVector q,
                       double mu = 5., double sigma = 1., double tau = 1.,
                       bool lower_tail = true, bool log_p = false) {
  int n = q.size();
  if (tau <= 0 || sigma <= 0) {
    NumericVector cdf(n, NA_REAL);
    return cdf;
  }

  NumericVector cdf(n);
  if (sigma < 1e-4){
    for (int i = 0; i < n; i++){
      cdf[i] = R::pexp(q[i] - mu, tau, lower_tail, log_p);
    }
    return cdf;
  }

  for (int i = 0; i < n; i++){
    if (!traits::is_infinite<REALSXP>(q[i])){
      //      if (tau > .05 * sigma){
      double z_i = q[i] - mu - (sigma * sigma) / tau;
      cdf[i] = R::pnorm((q[i] - mu) / sigma, 0., 1., true, false) - std::exp(std::log(R::pnorm(z_i / sigma, 0., 1., true, false)) + (std::pow((mu + (sigma * sigma / tau)), 2) - mu * mu - 2. * q[i] * (sigma * sigma / tau)) / (2. * sigma * sigma));
      //      } else {
      //        cdf[i] = R::pnorm(q[i], mu, sigma, true, false);
      //      }
    } else {
      if (q[i] < 0) {
        cdf[i] = 0.;
      } else {
        cdf[i] = 1.;
      }
    }
  }
  if (!lower_tail){
    for(int i = 0; i < n; i++){
      cdf[i] = 1. - cdf[i];
    }
  }
  if (log_p){
    for(int i = 0; i < n; i++){
      cdf[i] = std::log(cdf[i]);
    }
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector dEXG_old(NumericVector x,
                       double mu = 5., double sigma = 1., double tau = 1.,
                       bool log_d = false) {
  int n = x.size();
  if (tau <= 0 || sigma <= 0) {
    NumericVector pdf(n, NA_REAL);
    return pdf;
  }

  NumericVector pdf(n);
  if (sigma < 1e-4){
    for (int i = 0; i < n; i++){
      pdf[i] = R::dexp(x[i] - mu, tau, log_d);
    }
    return pdf;
  }

  for (int i = 0; i < n; i++){
    //    if (tau > .05 * sigma){
    double z_i = x[i] - mu - (sigma * sigma) / tau;
    pdf[i] = - std::log(tau) - (z_i + (sigma * sigma)/(2. * tau)) / tau + std::log(R::pnorm(z_i / sigma, 0., 1., true, false));
    //    } else {
    //      pdf[i] = R::dnorm(x[i], mu, sigma, true);
    //    }
  }
  if (!log_d){
    for(int i = 0; i < n; i++){
      pdf[i] = std::exp(pdf[i]);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector dEXGrace_old(NumericMatrix dt,
                           NumericVector mu, NumericVector sigma, NumericVector tau){
  int n = mu.size();
  NumericVector out(dt.nrow());
  out = dEXG_old(dt(0, _), mu[0], sigma[0], tau[0], false);
  for (int i = 1; i < n; i++){
    out = out * pEXG_old(dt(i, _), mu[i], sigma[i], tau[i], false, false);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_exg_old(NumericVector t,
                             NumericVector mu, NumericVector sigma, NumericVector tau,
                             double SSD){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dEXGrace_old(dt, mu, sigma, tau);
}


// [[Rcpp::export]]
NumericVector pTEXG_vec(
    NumericVector q, double mu = 5., double sigma = 1., double tau = 1., double lb = .05,
    bool lower_tail = true, bool log_p = false
) {
  int n = q.size();
  if (tau <= 0. || sigma <= 0.) {
    NumericVector cdf(n, NA_REAL);
    return cdf;
  }
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    cdf[i] = ptexg(q[i], mu, sigma, tau, lb, R_PosInf, lower_tail, log_p);
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector dTEXG_vec(
    NumericVector x, double mu = 5., double sigma = 1., double tau = 1., double lb = .05,
    bool log_d = false
) {
  int n = x.size();
  if (tau <= 0. || sigma <= 0.) {
    NumericVector pdf(n, NA_REAL);
    return pdf;
  }
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    pdf[i] = dtexg(x[i], mu, sigma, tau, lb, R_PosInf, log_d);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector dTEXGrace(
    NumericMatrix dt,
    NumericVector mu, NumericVector sigma, NumericVector tau, NumericVector lb
){
  int n = mu.size();
  NumericVector out(dt.nrow());
  out = dTEXG_vec(dt(0, _), mu[0], sigma[0], tau[0], lb[0]);
  for (int i = 1; i < n; i++){
    out = out * pTEXG_vec(dt(i, _), mu[i], sigma[i], tau[i], lb[i], false);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_texg(
    NumericVector t,
    NumericVector mu, NumericVector sigma, NumericVector tau, NumericVector lb,
    double SSD
){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dTEXGrace(dt, mu, sigma, tau, lb);
}

#endif
