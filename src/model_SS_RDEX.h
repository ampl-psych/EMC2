#ifndef ss_rdex_h
#define ss_rdex_h

#include <cmath>
#include <Rcpp.h>
#include "wald_functions.h"
#include "exgaussian_functions.h"
#include "race_integrate.h"
using namespace Rcpp;

// ----------------------------------------------------------------------------
// TODO write new code, borrowing from model_SS_EXG.h
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// OLD STUFF BELOW
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector pEXG_RDEX(NumericVector q,
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
NumericVector dEXG_RDEX(NumericVector x,
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
NumericVector dWald_RDEX(NumericVector t, double v,
                         double B, double A, double t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0;
    if (t[i] <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt(t[i], B + .5 * A, v, .5 * A);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald_RDEX(NumericVector t, double v,
                         double B, double A, double t0){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0;
    if (t[i] <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt(t[i], B + .5 * A, v, .5 * A);
    }
  }
  return cdf;
}


// [[Rcpp::export]]
NumericVector dRDEXrace(NumericMatrix dt,
                        double mu, double sigma, double tau,
                        NumericVector v, NumericVector B, NumericVector A,
                        NumericVector t0, bool exgWinner = true){
  int n = v.size();
  NumericVector out(dt.nrow());
  if (exgWinner){
    out = dEXG_RDEX(dt(0, _), mu, sigma, tau, false);
    out = out * (1. - pWald_RDEX(dt(1, _), v[0], B[0], A[0], t0[0]));
  } else {
    out = dWald_RDEX(dt(0, _), v[0], B[0], A[0], t0[0]);
    out = out * (1. - pEXG_RDEX(dt(1, _), mu, sigma, tau));
  }
  for (int i = 1; i < n; i++){
    out = out * (1. - pWald_RDEX(dt(i + 1, _), v[i], B[i], A[i], t0[i]));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_rdex(NumericVector t, int n_acc,
                          double mu, double sigma, double tau,
                          NumericVector v, NumericVector B, NumericVector A,
                          NumericVector t0, double SSD){
  NumericVector tmp( (n_acc + 1) * t.size());
  tmp = rep_each(t, n_acc + 1) + SSD;
  NumericMatrix dt(n_acc + 1, t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dRDEXrace(dt, mu, sigma, tau, v, B, A, t0);
}

#endif
