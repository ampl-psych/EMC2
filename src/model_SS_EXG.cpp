#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector pEXG(NumericVector q,
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
      if (tau > .05 * sigma){
        double z_i = q[i] - mu - (sigma * sigma) / tau;
        cdf[i] = R::pnorm((q[i] - mu) / sigma, 0., 1., true, false) - std::exp(std::log(R::pnorm(z_i / sigma, 0., 1., true, false)) + (std::pow((mu + (sigma * sigma / tau)), 2) - mu * mu - 2. * q[i] * (sigma * sigma / tau)) / (2. * sigma * sigma));
      } else {
        cdf[i] = R::pnorm(q[i], mu, sigma, true, false);
      }
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
NumericVector dEXG(NumericVector x,
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
    if (tau > .05 * sigma){
      double z_i = x[i] - mu - (sigma * sigma) / tau;
      pdf[i] = - std::log(tau) - (z_i + (sigma * sigma)/(2. * tau)) / tau + std::log(R::pnorm(z_i / sigma, 0., 1., true, false));
    } else {
      pdf[i] = R::dnorm(x[i], mu, sigma, true);
    }
  }
  if (!log_d){
    for(int i = 0; i < n; i++){
      pdf[i] = std::exp(pdf[i]);
    }
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector dEXGrace(NumericMatrix dt,
                       NumericVector mu, NumericVector sigma, NumericVector tau){
  int n = mu.size();
  NumericVector out(dt.nrow());
  out = dEXG(dt(0, _), mu[0], sigma[0], tau[0], false);
  for (int i = 1; i < n; i++){
    out = out * pEXG(dt(i, _), mu[i], sigma[i], tau[i], false, false);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_exg(NumericVector t,
                         NumericVector mu, NumericVector sigma, NumericVector tau,
                         double SSD){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dEXGrace(dt, mu, sigma, tau);
}
