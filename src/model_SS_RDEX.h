#ifndef ss_rdex_h
#define ss_rdex_h

#include <cmath>
#include <vector>
#include <Rcpp.h>
#include "model_RDM.h"          // cens Wald functions (digt/pigt/digt0/pigt0)
#include "exgaussian_functions.h"
#include "ss_integrate.h"      // cens hcubature wrapper + finite window
using namespace Rcpp;

// ----------------------------------------------------------------------------
// HYBRID WALD / EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

// wrapper around Wald log PDF for race function
inline NumericVector rdex_go_lpdf(
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

    double dt_i = rt[i] - pars(i, 3);
    double log_d = R_NegInf;
    if (dt_i > 0.) {
      log_d = std::log(
        digt(
          dt_i,
          (pars(i, 1) / pars(i, 4)) + .5 * (pars(i, 2) / pars(i, 4)),
          pars(i, 0) / pars(i, 4),
          .5 * (pars(i, 2) / pars(i, 4))
        )
      );
    }

    out[k] = R_FINITE(log_d) ? log_d : min_ll;

    k++;
  }

  return(out);
}

// wrapper around Wald log complementary CDF for race function
inline NumericVector rdex_go_lccdf(
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

    double dt_i = rt[i] - pars(i, 3);
    double log_s = 0.; // log(1)
    if (dt_i > 0.) {
      log_s = log1m(
        pigt(
          dt_i,
          (pars(i, 1) / pars(i, 4)) + .5 * (pars(i, 2) / pars(i, 4)),
          pars(i, 0) / pars(i, 4),
          .5 * (pars(i, 2) / pars(i, 4))
        )
      );
    }

    out[k] = R_FINITE(log_s) ? log_s : min_ll;

    k++;
  }

  return(out);
}

// ----------------------------------------------------------------------------
// STOP-SUCCESS INTEGRAL
// P(stop process (truncated exG) finishes before every go accumulator (Wald) |
// SSD), not accounting for trigger failure and go failure. Integrand in
// hcubature callback form (ss_integrate.h); x is stop-relative time, go
// survivors evaluated at x + SSD.
// ----------------------------------------------------------------------------

struct rdex_stop_success_pars {
  double SSD;
  double min_ll;
  // stop params (truncated EXG): columns muS=5, sigmaS=6, tauS=7, exgS_lb=10
  double muS, sigS, tauS, lbS;
  // precomputed go Wald params: columns v=0, B=1, A=2, t0=3, s=4
  int n_go;
  std::vector<double> alpha, nu, gamma, t0;

  rdex_stop_success_pars(double SSD_, double min_ll_, const NumericMatrix& pars_)
    : SSD(SSD_),
      min_ll(min_ll_),
      muS(pars_(0, 5)), sigS(pars_(0, 6)), tauS(pars_(0, 7)), lbS(pars_(0, 10)),
      n_go(pars_.nrow()), alpha(n_go), nu(n_go), gamma(n_go), t0(n_go)
  {
    for (int i = 0; i < n_go; ++i) {
      double s = pars_(i, 4);
      alpha[i] = (pars_(i, 1) / s) + .5 * (pars_(i, 2) / s);
      nu[i]    =  pars_(i, 0) / s;
      gamma[i] = .5 * (pars_(i, 2) / s);
      t0[i]    =  pars_(i, 3);
    }
  }
};

static int rdex_stop_success_integrand(unsigned /*dim*/, const double* x, void* p,
                                       unsigned /*fdim*/, double* out) {
  const rdex_stop_success_pars* w = static_cast<const rdex_stop_success_pars*>(p);
  const double xx = x[0];
  // log density of stop process finishing at time xx
  double log_fS = dtexg(xx, w->muS, w->sigS, w->tauS, w->lbS, R_PosInf, true);
  if (!R_FINITE(log_fS)) { log_fS = w->min_ll; }
  // log probability that no go accumulator has finished by xx + SSD
  double log_S_go = 0.0;
  for (int i = 0; i < w->n_go; ++i) {
    double dt_i = (xx + w->SSD) - w->t0[i];
    if (dt_i > 0.0) {
      double log_Si = std::log(1.0 - pigt(dt_i, w->alpha[i], w->nu[i], w->gamma[i]));
      if (!R_FINITE(log_Si)) { log_Si = w->min_ll; }
      log_S_go += log_Si;
    }
  }
  // output: sum of (1) log winner density (stop) and (2) sum of log survival
  // probabilities (go), exponentiated to put on likelihood scale
  out[0] = std::exp(log_fS + log_S_go);
  return 0;
}

// log P(stop wins), integrating over stop finish times in [lb, upper] (stop-
// relative time), clipped to the finite window of ss_integrate.h.
static inline double ss_rdex_stop_success_lpdf(
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
  rdex_stop_success_pars w(SSD, min_ll, pars);
  const double lo = ss_stop_window_lo(w.lbS, w.muS, w.sigS, k_sigma);
  const double hi = ss_stop_window_hi(upper, w.muS, w.sigS, w.tauS, k_sigma, k_tau);
  const std::size_t max_eval = static_cast<std::size_t>(max_subdiv) * 64;
  double res = ss_integrate(rdex_stop_success_integrand, &w, lo, hi,
                            abs_tol, rel_tol, max_eval);
  return (!R_FINITE(res) || res <= 0.0) ? min_ll : std::log(res);
}

// ----------------------------------------------------------------------------
// EXPORTS USED BY THE R REFERENCE LIKELIHOOD (pstopHybrid in R/model_SS.R)
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dWald_RDEX(
    NumericVector t,
    double v, double B, double A, double t0, double s
) {
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    const double tt = t[i] - t0;   // do not mutate the input vector
    pdf[i] = 0.;
    if (tt > 0.) {
      pdf[i] = digt(tt, (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
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
    const double tt = t[i] - t0;   // do not mutate the input vector
    cdf[i] = 0.;
    if (tt > 0.) {
      cdf[i] = pigt(tt, (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
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

