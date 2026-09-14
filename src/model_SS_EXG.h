#ifndef ss_exg_h
#define ss_exg_h

#include <cmath>
#include <vector>
#include "nan_check.h" // is_finite
#include "r_constants.h" // pos / neg infinity; na_real
#include "exgaussian_functions.h"
using namespace Rcpp;

// R-facing helpers for the stop-signal ex-Gaussian race (used by the R
// reference likelihood, pstopTEXG in R/model_SS.R). The C++ likelihood lives
// in ss_fast.h.

// ----------------------------------------------------------------------------
// EXPORTS USED BY THE R REFERENCE LIKELIHOOD (pstopTEXG in R/model_SS.R)
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector pTEXG_vec(
    NumericVector q, double mu = 5., double sigma = 1., double tau = 1., double lb = .05,
    bool lower_tail = true, bool log_p = false
) {
  int n = q.size();
  if (tau <= 0. || sigma <= 0.) {
    NumericVector cdf(n, na_real());
    return cdf;
  }
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    cdf[i] = ptexg(q[i], mu, sigma, tau, lb, pos_inf(), lower_tail, log_p);
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
    NumericVector pdf(n, na_real());
    return pdf;
  }
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    pdf[i] = dtexg(x[i], mu, sigma, tau, lb, pos_inf(), log_d);
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

