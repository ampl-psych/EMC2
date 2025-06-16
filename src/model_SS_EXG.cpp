#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
#include "exgaussian_functions.h"
#include "race_integrate.h"

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

// [[Rcpp::export]]
NumericVector dEXGrace(NumericMatrix dt,
                       NumericVector mu, NumericVector sigma, NumericVector tau){
  int n = mu.size();
  NumericVector out(dt.nrow());
  out = dexg_c(dt(0, _), mu[0], sigma[0], tau[0], false);
  for (int i = 1; i < n; i++){
    out = out * pexg_c(dt(i, _), mu[i], sigma[i], tau[i], false, false);
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
