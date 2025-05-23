#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector pEXG(
    NumericVector q,
    double mu = 5.,
    double sigma = 1.,
    double tau = 1.,
    bool lower_tail = true,
    bool log_p = false
) {

  int n = q.size();
  NumericVector cdf(n);

  // return NA for invalid parameter values
  if (sigma <= 0 || tau <= 0) {
    std::fill(cdf.begin(), cdf.end(), NA_REAL);
    return cdf;
  }

  // if Gaussian SD is practically zero, treat as shifted exponential
  if (sigma < 1e-4) {
    for (int i = 0; i < n; i++) {
      cdf[i] = R::pexp(q[i] - mu, tau, lower_tail, log_p);
    }
    return cdf;
  }

  for (int i = 0; i < n; i++) {
    if (!traits::is_infinite<REALSXP>(q[i])) {
      if (tau > .05 * sigma) {
        // if exponential rate has reasonable value, treat as ex-Gaussian
        double delta = sigma * sigma / tau;
        double z = q[i] - mu - delta;
        double norm1 = R::pnorm((q[i] - mu) / sigma, 0., 1., true, false);
        double norm2 = R::pnorm(z / sigma, 0., 1., true, false);
        double exponent = (
          std::pow(mu + delta, 2.) - mu * mu - 2. * q[i] * delta
        ) / (2. * sigma * sigma);
        cdf[i] = norm1 - std::exp(std::log(norm2) + exponent);
      } else {
        // if not, treat as Gaussian
        cdf[i] = R::pnorm(q[i], mu, sigma, true, false);
      }
    } else {
      // handle negative and positive infinity input q
      cdf[i] = (q[i] < 0.0) ? 0. : 1.;
    }
  }

  if (!lower_tail || log_p) {
    for (int i = 0; i < n; i++) {
      // complementary CDF
      if (!lower_tail) cdf[i] = 1. - cdf[i];
      // log transform
      if (log_p) cdf[i] = std::log(cdf[i]);
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
