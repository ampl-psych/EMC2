#ifndef ss_exg_h
#define ss_exg_h

#include <cmath>
#include <vector>
#include <Rcpp.h>
#include "nan_check.h" // is_finite
#include "exgaussian_functions.h"
#include "ss_integrate.h"      // cens hcubature wrapper + finite window
using namespace Rcpp;

// ----------------------------------------------------------------------------
// TRUNCATED EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

// wrapper around truncated ex-Gaussian log PDF for race function
inline NumericVector texg_go_lpdf(
    // single RT, broadcast to repped vector (one per accumulator)
    NumericVector rt,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // accumulator index: for which accumulator(s) should log density be computed?
    LogicalVector idx,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {

  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;

  NumericVector out(n_acc_selected);
  int k = 0;

  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;

    // input args: x, mu, sigma, tau, exg_lb, upper = Inf, log_d = TRUE
    double log_d = dtexg(
      rt[i], pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, true
    );
    out[k] = is_finite(log_d) ? log_d : min_ll;

    k++;
  }

  return(out);
}

// wrapper around truncated ex-Gaussian log complementary CDF for race function
inline NumericVector texg_go_lccdf(
    // single RT, broadcast to repped vector (one per accumulator)
    NumericVector rt,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // accumulator index: for which accumulator(s) should survivor probability be computed?
    LogicalVector idx,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {

  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;

  NumericVector out(n_acc_selected);
  int k = 0;

  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;

    // input args: q, mu, sigma, tau, exg_lb, upper = Inf, lower_tail = FALSE, log_p = TRUE
    double log_s = ptexg(
      rt[i], pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, false, true
    );
    out[k] = is_finite(log_s) ? log_s : min_ll;

    k++;
  }

  return(out);
}

// ----------------------------------------------------------------------------
// STOP-SUCCESS INTEGRAL
// P(stop process finishes before every go accumulator | SSD), not accounting
// for trigger failure and go failure. Integrand in hcubature callback form
// (ss_integrate.h); x is stop-relative time, go survivors evaluated at x + SSD.
// ----------------------------------------------------------------------------

struct texg_stop_success_pars {
  double SSD;
  double min_ll;
  // stop params (truncated EXG)
  double muS, sigS, tauS, lbS;
  // go params per accumulator (truncated EXG)
  int n_go;
  std::vector<double> muG, sigG, tauG, lbG;

  texg_stop_success_pars(double SSD_, double min_ll_, const NumericMatrix& pars_)
    : SSD(SSD_),
      min_ll(min_ll_),
      muS(pars_(0, 3)), sigS(pars_(0, 4)), tauS(pars_(0, 5)), lbS(pars_(0, 9)),
      n_go(pars_.nrow()), muG(n_go), sigG(n_go), tauG(n_go), lbG(n_go)
  {
    for (int i = 0; i < n_go; ++i) {
      muG[i]  = pars_(i, 0);
      sigG[i] = pars_(i, 1);
      tauG[i] = pars_(i, 2);
      lbG[i]  = pars_(i, 8);
    }
  }
};

static int texg_stop_success_integrand(unsigned /*dim*/, const double* x, void* p,
                                       unsigned /*fdim*/, double* out) {
  const texg_stop_success_pars* w = static_cast<const texg_stop_success_pars*>(p);
  const double xx = x[0];
  // log density of stop process finishing at time xx
  double log_fS = dtexg(xx, w->muS, w->sigS, w->tauS, w->lbS, R_PosInf, true);
  if (!is_finite(log_fS)) { log_fS = w->min_ll; }
  // log probability that no go accumulator has finished by xx + SSD
  double log_S_go = 0.0;
  for (int i = 0; i < w->n_go; ++i) {
    double log_Si = ptexg(xx + w->SSD, w->muG[i], w->sigG[i], w->tauG[i], w->lbG[i],
                          R_PosInf, false, true);
    if (!is_finite(log_Si)) { log_Si = w->min_ll; }
    log_S_go += log_Si;
  }
  // output: sum of (1) log winner density (stop) and (2) sum of log survival
  // probabilities (go), exponentiated to put on likelihood scale
  out[0] = std::exp(log_fS + log_S_go);
  return 0;
}

// log P(stop wins), integrating over stop finish times in [lb, upper] (stop-
// relative time), clipped to the finite window of ss_integrate.h.
static inline double ss_texg_stop_success_lpdf(
    double SSD,
    NumericMatrix pars,
    double min_ll,
    double upper = R_PosInf,
    int max_subdiv = 30,
    double abs_tol = 1e-5,
    double rel_tol = 1e-4,
    double k_sigma = SS_WINDOW_K_SIGMA,
    double k_tau = SS_WINDOW_K_TAU
) {
  texg_stop_success_pars w(SSD, min_ll, pars);
  const double lo = ss_stop_window_lo(w.lbS, w.muS, w.sigS, k_sigma);
  const double hi = ss_stop_window_hi(upper, w.muS, w.sigS, w.tauS, k_sigma, k_tau);
  // max_subdiv is an evaluation-budget proxy (kept from the upstream API)
  const std::size_t max_eval = static_cast<std::size_t>(max_subdiv) * 64;
  double res = ss_integrate(texg_stop_success_integrand, &w, lo, hi,
                            abs_tol, rel_tol, max_eval);
  return (!is_finite(res) || res <= 0.0) ? min_ll : std::log(res);
}

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

