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

  // if exponential rate has very small value relative to Gaussian SD, treat as
  // Gaussian
  if (tau < .05 * sigma) {
    for (int i = 0; i < n; i++) {
      cdf[i] = R::pnorm(q[i], mu, sigma, lower_tail, log_p);
    }
    return cdf;
  }

  // main ex-Gaussian computation:

  const double mu2 = mu * mu;
  const double sig2 = sigma * sigma;
  const double sig2_tau = sig2 / tau;

  // evaluate ex-Gaussian cumulative distribution function
  for (int i = 0; i < n; i++) {
    double q_i = q[i];
    if (!traits::is_infinite<REALSXP>(q_i)) {
      double z = q_i - mu - sig2_tau;
      double norm1 = R::pnorm((q_i - mu) / sigma, 0., 1., true, false);
      double norm2 = R::pnorm(z / sigma, 0., 1., true, false);
      double exponent = (
        std::pow(mu + sig2_tau, 2.) - mu2 - 2. * q_i * sig2_tau
      ) / (2. * sig2);
      cdf[i] = norm1 - std::exp(std::log(norm2) + exponent);
    } else {
      // handle negative and positive infinity input q
      cdf[i] = (q_i < 0.0) ? 0. : 1.;
    }
  }

  if (!lower_tail || log_p) {
    for (int i = 0; i < n; i++) {
      // convert to complementary CDF a.k.a. survival function
      if (!lower_tail) cdf[i] = 1. - cdf[i];
      // convert to log-probability
      if (log_p) cdf[i] = std::log(cdf[i]);
    }
  }

  return cdf;
}


// [[Rcpp::export]]
NumericVector dEXG(
    NumericVector x,
    double mu = 5.,
    double sigma = 1.,
    double tau = 1.,
    bool log_d = false
) {

  int n = x.size();
  NumericVector pdf(n);

  // return NA for invalid parameter values
  if (sigma <= 0 || tau <= 0) {
    std::fill(pdf.begin(), pdf.end(), NA_REAL);
    return pdf;
  }

  // if Gaussian SD is practically zero, treat as shifted exponential
  if (sigma < 1e-4) {
    for (int i = 0; i < n; i++) {
      pdf[i] = R::dexp(x[i] - mu, tau, log_d);
    }
    return pdf;
  }

  // if exponential rate has very small value relative to Gaussian SD, treat as
  // Gaussian
  if (tau < .05 * sigma) {
    for (int i = 0; i < n; i++) {
      pdf[i] = R::dnorm(x[i], mu, sigma, log_d);
    }
    return pdf;
  }

  // main ex-Gaussian computation:

  const double sig2 = sigma * sigma;

  // evaluate ex-Gaussian log probability density function
  for (int i = 0; i < n; i++) {
    double z = x[i] - mu - sig2 / tau;
    double log_phi = std::log(R::pnorm(z / sigma, 0., 1., true, false));
    double log_exp = -std::log(tau) - (z + sig2 / (2. * tau)) / tau;
    pdf[i] = log_exp + log_phi;
  }

  if (!log_d){
    // convert from log-probability to probability
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
