#ifndef ss_rdex_h
#define ss_rdex_h

#include <cmath>
#include <string>
#include <vector>
#include <Rcpp.h>
#ifndef rdm_h
#include "wald_functions.h"
#endif
#include "exgaussian_functions.h"
#include "pnorm_utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dWald_RDEX(
    NumericVector t,
    double v, double B, double A, double t0, double s
) {
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    double t_eff = t[i] - t0;
    pdf[i] = 0.;
    if (t_eff > 0.) {
      pdf[i] = digt(t_eff, (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald_RDEX(
    NumericVector t,
    double v, double B, double A, double t0, double s
) {
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++) {
    double t_eff = t[i] - t0;
    cdf[i] = 0.;
    if (t_eff > 0.) {
      cdf[i] = pigt(t_eff, (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
    }
  }
  return cdf;
}


// [[Rcpp::export]]
NumericVector pTEXG_RDEX(
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
NumericVector dTEXG_RDEX(
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
NumericVector dRDEXrace(
    NumericMatrix dt,
    double mu, double sigma, double tau, double lb,
    NumericVector v, NumericVector B, NumericVector A, NumericVector t0, NumericVector s,
    bool exgWinner = true
) {
  int n = v.size();
  NumericVector out(dt.nrow());
  if (exgWinner) {
    out = dTEXG_RDEX(dt(0, _), mu, sigma, tau, lb);
    out = out * (1. - pWald_RDEX(dt(1, _), v[0], B[0], A[0], t0[0], s[0]));
  } else {
    out = dWald_RDEX(dt(0, _), v[0], B[0], A[0], t0[0], s[0]);
    out = out * (1. - pTEXG_RDEX(dt(1, _), mu, sigma, tau, lb));
  }
  for (int i = 1; i < n; i++){
    out = out * (1. - pWald_RDEX(dt(i + 1, _), v[i], B[i], A[i], t0[i], s[i]));
  }
  return out;
}


// [[Rcpp::export]]
NumericVector stopfn_rdex(
    NumericVector t, int n_acc,
    double mu, double sigma, double tau, double lb,
    NumericVector v, NumericVector B, NumericVector A, NumericVector t0, NumericVector s,
    double SSD
) {
  NumericVector tmp((n_acc + 1) * t.size());
  tmp = rep_each(t, n_acc + 1) + SSD;
  NumericMatrix dt(n_acc + 1, t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dRDEXrace(dt, mu, sigma, tau, lb, v, B, A, t0, s);
}

#endif
