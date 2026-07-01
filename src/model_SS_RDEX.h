#ifndef ss_rdex_h
#define ss_rdex_h

#include <cmath>
#include <vector>
#include <cstddef>
#include <Rcpp.h>
// cens port: Wald go-process uses cens's own digt/pigt(t,k,l,a) (model_RDM.h),
// which are the equivalents of Zach's digt_impl/pigt_impl(t,k,l,a) (the 4-arg
// (t,k,l,a) form, default a-threshold). The stop-success integral runs on cens's
// adaptive cubature (stop_quad.h). No GSL/GL.
#include "exgaussian_functions.h"
#include "stop_quad.h"
#include "model_RDM.h"   // digt(t,k,l,a), pigt(t,k,l,a)

using namespace Rcpp;

// Parameter-matrix column layout (rows = accumulators):
//   go (Wald): 0 v, 1 B, 2 A, 3 t0, 4 s
//   stop (truncated ex-Gaussian): 5 muS, 6 sigS, 7 tauS, 10 exgS_lb

// ----------------------------------------------------------------------------
// STOP-SUCCESS INTEGRAND (cens hcubature callback form)
// Integrate over the stop finishing time x: f_S(x) * prod_i S_Wald,i(x + SSD).
// ----------------------------------------------------------------------------
struct SSRDEXWrapper { NumericMatrix pars; double SSD; };

static int ss_rdex_integrand(unsigned /*dim*/, const double* x, void* p,
                             unsigned /*fdim*/, double* retval) {
  SSRDEXWrapper* wp = static_cast<SSRDEXWrapper*>(p);
  const double xx = x[0];
  const int n_acc = wp->pars.nrow();
  const double muS = wp->pars(0, 5), sigS = wp->pars(0, 6), tauS = wp->pars(0, 7);
  const double lbS = wp->pars(0, 10);
  double fS = dtexg(xx, muS, sigS, tauS, lbS, R_PosInf, false);
  if (fS <= 0.0) { retval[0] = 0.0; return 0; }
  double S_go_all = 1.0;
  for (int i = 0; i < n_acc; ++i) {
    const double s = wp->pars(i, 4);
    const double alpha = (wp->pars(i, 1) / s) + .5 * (wp->pars(i, 2) / s);
    const double nu    =  wp->pars(i, 0) / s;
    const double gamma = .5 * (wp->pars(i, 2) / s);
    const double t0    =  wp->pars(i, 3);
    const double dt_i  = (xx + wp->SSD) - t0;
    double Si = 1.0;
    if (dt_i > 0.) Si = 1.0 - pigt(dt_i, alpha, nu, gamma);
    S_go_all *= Si;
    if (S_go_all <= 0.0) { retval[0] = 0.0; return 0; }
  }
  retval[0] = fS * S_go_all;
  return 0;
}

// ----------------------------------------------------------------------------
// HYBRID WALD / EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

// wrapper around Wald log PDF for race function
NumericVector rdex_go_lpdf(NumericVector rt, NumericMatrix pars,
                           LogicalVector idx, double min_ll) {
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
      log_d = std::log(digt(
        dt_i,
        (pars(i, 1) / pars(i, 4)) + .5 * (pars(i, 2) / pars(i, 4)),
        pars(i, 0) / pars(i, 4),
        .5 * (pars(i, 2) / pars(i, 4))));
    }
    out[k] = R_FINITE(log_d) ? log_d : R_NegInf;
    k++;
  }
  return(out);
}

// wrapper around Wald log complementary CDF for race function
NumericVector rdex_go_lccdf(NumericVector rt, NumericMatrix pars,
                            LogicalVector idx, double min_ll) {
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
      log_s = log1m(pigt(
        dt_i,
        (pars(i, 1) / pars(i, 4)) + .5 * (pars(i, 2) / pars(i, 4)),
        pars(i, 0) / pars(i, 4),
        .5 * (pars(i, 2) / pars(i, 4))));
    }
    out[k] = emc2_isfinite(log_s) ? log_s : R_NegInf;
    k++;
  }
  return(out);
}

// go race log likelihood function, not accounting for go failure
double ss_rdex_go_lpdf(double RT, NumericMatrix pars,
                       LogicalVector winner, double min_ll) {
  const int n_acc = pars.nrow();
  double logpdf_winner = 0.0;
  double logsurv_losers = 0.0;
  bool winner_found = false;
  for (int i = 0; i < n_acc; ++i) {
    double s = pars(i, 4);
    double alpha = (pars(i, 1) / s) + .5 * (pars(i, 2) / s);
    double nu    =  pars(i, 0) / s;
    double gamma = .5 * (pars(i, 2) / s);
    double t0    =  pars(i, 3);
    double dt_i = RT - t0;
    if (winner[i]) {
      double ld = R_NegInf;
      if (dt_i > 0.) ld = std::log(digt(dt_i, alpha, nu, gamma));
      logpdf_winner += emc2_isfinite(ld) ? ld : min_ll;
      winner_found = true;
    } else {
      double ls = 0.; // log(1)
      if (dt_i > 0.) ls = log1m(pigt(dt_i, alpha, nu, gamma));
      logsurv_losers += emc2_isfinite(ls) ? ls : min_ll;
    }
  }
  if (!winner_found) return min_ll;
  double out = logpdf_winner + logsurv_losers;
  return emc2_isfinite(out) ? out : min_ll;
}

// go vs stop race log likelihood for failed-inhibition stop trials.
double ss_rdex_stop_fail_lpdf(double RT, double SSD, NumericMatrix pars,
                              LogicalVector winner, double min_ll) {
  double go_lprob = ss_rdex_go_lpdf(RT, pars, winner, min_ll);
  // stop pars: muS=5, sigS=6, tauS=7, lbS=10
  double stop_survivor_lprob = ptexg(
    RT - SSD, pars(0, 5), pars(0, 6), pars(0, 7), pars(0, 10), R_PosInf, false, true);
  if (!emc2_isfinite(stop_survivor_lprob)) stop_survivor_lprob = min_ll;
  return go_lprob + stop_survivor_lprob;
}

// stop-success integral, "integrate" route: cens adaptive hcubature over the
// finite window (was GSL qags/qagiu).
static inline double ss_rdex_stop_success_lpdf(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int max_subdiv = 100,
    double abs_tol = 1e-8, double rel_tol = 1e-6,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU) {
  SSRDEXWrapper w = {pars, SSD};
  const double muS = pars(0, 5), sigS = pars(0, 6), tauS = pars(0, 7);
  const double lbS = pars(0, 10);
  const double ub_heur = muS + k_sigma * sigS + k_tau * tauS;
  const double lo = emc2_isfinite(lbS) ? lbS
                                       : muS - k_sigma * sigS - k_tau * tauS;
  const double ub = emc2_isfinite(upper) ? upper : ub_heur;
  if (!(ub > lo)) return min_ll;
  const std::size_t max_eval = static_cast<std::size_t>(max_subdiv) * 1000u;
  double res = ss_adaptive_integrate(ss_rdex_integrand, &w, lo, ub,
                                     abs_tol, rel_tol, max_eval);
  if (!emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// --- fixed (GL-labelled) route: composite Gauss-Kronrod over [lo, ub] ---
static inline double ss_rdex_stop_success_lpdf_gl(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int n_nodes = 64,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU) {
  SSRDEXWrapper w = {pars, SSD};
  const double muS = pars(0, 5), sigS = pars(0, 6), tauS = pars(0, 7);
  const double lbS = pars(0, 10);
  const double ub_heur = muS + k_sigma * sigS + k_tau * tauS;
  const double lo = emc2_isfinite(lbS) ? lbS
                                       : muS - k_sigma * sigS - k_tau * tauS;
  double ub = emc2_isfinite(upper) ? upper : ub_heur;
  if (!(ub > lo)) ub = lo + 1e-12;
  double res = ss_fixed_integrate(ss_rdex_integrand, &w, lo, ub, n_nodes);
  if (!emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// --- "auto" dispatch: RDEX has no analytic form (Wald survivor arguments are
//     nonlinear in t), so auto = fixed route with the tight-stop-density bump. ---
static inline double ss_rdex_stop_success_lpdf_autodisp(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int n_nodes = 64,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU) {
  const double muS = pars(0, 5), sigS = pars(0, 6), tauS = pars(0, 7);
  const double lbS = pars(0, 10);
  const double lo = emc2_isfinite(lbS) ? lbS
                                       : muS - k_sigma * sigS - k_tau * tauS;
  const double ub = emc2_isfinite(upper) ? upper
                                         : muS + k_sigma * sigS + k_tau * tauS;
  const int n_eff = gl_auto_nodes(n_nodes, lo, ub, sigS);
  return ss_rdex_stop_success_lpdf_gl(SSD, pars, min_ll, upper, n_eff,
                                      k_sigma, k_tau);
}

// --- LIVE entry point: reads stop_method_config(), set by calc_ll_manager()
//     from SSRDEX(stop_method =, stop_n_nodes =) once per likelihood call. ---
static inline double ss_rdex_stop_success_lpdf_live(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int max_subdiv = 100,
    double abs_tol = 1e-8, double rel_tol = 1e-6,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU) {
  const StopMethodConfig& c = stop_method_config();
  switch (c.method) {
  case STOP_METHOD_INTEGRATE:
    return ss_rdex_stop_success_lpdf(SSD, pars, min_ll, upper, max_subdiv,
                                     abs_tol, rel_tol, k_sigma, k_tau);
  case STOP_METHOD_GL:
    return ss_rdex_stop_success_lpdf_gl(SSD, pars, min_ll, upper, c.n_nodes,
                                        k_sigma, k_tau);
  default:   // STOP_METHOD_AUTO / STOP_METHOD_ANALYTIC (no RDEX analytic form)
    return ss_rdex_stop_success_lpdf_autodisp(SSD, pars, min_ll, upper,
                                              c.n_nodes, k_sigma, k_tau);
  }
}

// --- Rcpp test wrapper: integral value (not log) for one trial. ---
// [[Rcpp::export]]
double ss_rdex_stop_success_value(
    double SSD, NumericMatrix pars, std::string method = "integrate",
    double upper = -1.0, int n_nodes = 64,
    double k_sigma = 8.0, double k_tau = 16.0,
    int max_subdiv = 100, double abs_tol = 1e-8, double rel_tol = 1e-6) {
  const double min_ll = -1e10;
  if (upper <= 0.0) upper = R_PosInf;   // sentinel: <=0 means auto/Inf
  double lp;
  if (method == "gl") {
    lp = ss_rdex_stop_success_lpdf_gl(SSD, pars, min_ll, upper, n_nodes,
                                      k_sigma, k_tau);
  } else if (method == "auto") {
    lp = ss_rdex_stop_success_lpdf_autodisp(SSD, pars, min_ll, upper, n_nodes,
                                            k_sigma, k_tau);
  } else if (method == "live") {
    lp = ss_rdex_stop_success_lpdf_live(SSD, pars, min_ll, upper, max_subdiv,
                                        abs_tol, rel_tol, k_sigma, k_tau);
  } else {
    lp = ss_rdex_stop_success_lpdf(SSD, pars, min_ll, upper, max_subdiv,
                                   abs_tol, rel_tol, k_sigma, k_tau);
  }
  return (lp <= min_ll) ? 0.0 : std::exp(lp);
}

// ----------------------------------------------------------------------------
// R-FACING VECTOR HELPERS (used by the R stop_success / plotting routes)
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector pEXG_RDEX(NumericVector q, double mu = 5., double sigma = 1.,
                        double tau = 1., bool lower_tail = true, bool log_p = false) {
  int n = q.size();
  if (tau <= 0 || sigma <= 0) { NumericVector cdf(n, NA_REAL); return cdf; }
  NumericVector cdf(n);
  if (sigma < 1e-4) {
    for (int i = 0; i < n; i++) cdf[i] = R::pexp(q[i] - mu, tau, lower_tail, log_p);
    return cdf;
  }
  for (int i = 0; i < n; i++) {
    if (!traits::is_infinite<REALSXP>(q[i])) {
      if (tau > .05 * sigma) {
        double z_i = q[i] - mu - (sigma * sigma) / tau;
        double mu_term = mu + (sigma * sigma / tau);
        cdf[i] = pnorm_std((q[i] - mu) / sigma) - std::exp(pnorm_std(z_i / sigma, true, true) + (mu_term * mu_term - mu * mu - 2. * q[i] * (sigma * sigma / tau)) / (2. * sigma * sigma));
      } else {
        cdf[i] = pnorm_std((q[i] - mu) / sigma);
      }
    } else {
      cdf[i] = (q[i] < 0) ? 0. : 1.;
    }
  }
  if (!lower_tail) for (int i = 0; i < n; i++) cdf[i] = 1. - cdf[i];
  if (log_p)       for (int i = 0; i < n; i++) cdf[i] = std::log(cdf[i]);
  return cdf;
}

// [[Rcpp::export]]
NumericVector dEXG_RDEX(NumericVector x, double mu = 5., double sigma = 1.,
                        double tau = 1., bool log_d = false) {
  int n = x.size();
  if (tau <= 0 || sigma <= 0) { NumericVector pdf(n, NA_REAL); return pdf; }
  NumericVector pdf(n);
  if (sigma < 1e-4) {
    for (int i = 0; i < n; i++) pdf[i] = R::dexp(x[i] - mu, tau, log_d);
    return pdf;
  }
  for (int i = 0; i < n; i++) {
    if (tau > .05 * sigma) {
      double z_i = x[i] - mu - (sigma * sigma) / tau;
      pdf[i] = - std::log(tau) - (z_i + (sigma * sigma)/(2. * tau)) / tau + pnorm_std(z_i / sigma, true, true);
    } else {
      pdf[i] = R::dnorm(x[i], mu, sigma, true);
    }
  }
  if (!log_d) for (int i = 0; i < n; i++) pdf[i] = std::exp(pdf[i]);
  return pdf;
}

// [[Rcpp::export]]
NumericVector dWald_RDEX_old(NumericVector t, double v, double B, double A, double t0) {
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    t[i] = t[i] - t0;
    pdf[i] = (t[i] <= 0) ? 0. : digt(t[i], B + .5 * A, v, .5 * A);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector dWald_RDEX(NumericVector t, double v, double B, double A, double t0, double s) {
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    t[i] = t[i] - t0;
    pdf[i] = 0.;
    if (t[i] > 0.) pdf[i] = digt(t[i], (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pWald_RDEX_old(NumericVector t, double v, double B, double A, double t0) {
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++) {
    t[i] = t[i] - t0;
    cdf[i] = (t[i] <= 0) ? 0. : pigt(t[i], B + .5 * A, v, .5 * A);
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector pWald_RDEX(NumericVector t, double v, double B, double A, double t0, double s) {
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++) {
    t[i] = t[i] - t0;
    cdf[i] = 0.;
    if (t[i] > 0.) cdf[i] = pigt(t[i], (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector pTEXG_RDEX(NumericVector q, double mu = 5., double sigma = 1.,
                         double tau = 1., double lb = .05,
                         bool lower_tail = true, bool log_p = false) {
  int n = q.size();
  if (tau <= 0. || sigma <= 0.) { NumericVector cdf(n, NA_REAL); return cdf; }
  NumericVector cdf(n);
  for (int i = 0; i < n; i++)
    cdf[i] = ptexg(q[i], mu, sigma, tau, lb, R_PosInf, lower_tail, log_p);
  return cdf;
}

// [[Rcpp::export]]
NumericVector dTEXG_RDEX(NumericVector x, double mu = 5., double sigma = 1.,
                         double tau = 1., double lb = .05, bool log_d = false) {
  int n = x.size();
  if (tau <= 0. || sigma <= 0.) { NumericVector pdf(n, NA_REAL); return pdf; }
  NumericVector pdf(n);
  for (int i = 0; i < n; i++)
    pdf[i] = dtexg(x[i], mu, sigma, tau, lb, R_PosInf, log_d);
  return pdf;
}

// [[Rcpp::export]]
NumericVector dRDEXrace_old(NumericMatrix dt, double mu, double sigma, double tau,
                            NumericVector v, NumericVector B, NumericVector A,
                            NumericVector t0, bool exgWinner = true) {
  int n = v.size();
  NumericVector out(dt.ncol());
  if (exgWinner) {
    out = dEXG_RDEX(dt(0, _), mu, sigma, tau, false);
    out = out * (1. - pWald_RDEX_old(dt(1, _), v[0], B[0], A[0], t0[0]));
  } else {
    out = dWald_RDEX_old(dt(0, _), v[0], B[0], A[0], t0[0]);
    out = out * (1. - pEXG_RDEX(dt(1, _), mu, sigma, tau));
  }
  for (int i = 1; i < n; i++)
    out = out * (1. - pWald_RDEX_old(dt(i + 1, _), v[i], B[i], A[i], t0[i]));
  return out;
}

// [[Rcpp::export]]
NumericVector dRDEXrace(NumericMatrix dt, double mu, double sigma, double tau, double lb,
                        NumericVector v, NumericVector B, NumericVector A,
                        NumericVector t0, NumericVector s, bool exgWinner = true) {
  int n = v.size();
  NumericVector out(dt.ncol());
  if (exgWinner) {
    out = dTEXG_RDEX(dt(0, _), mu, sigma, tau, lb);
    out = out * (1. - pWald_RDEX(dt(1, _), v[0], B[0], A[0], t0[0], s[0]));
  } else {
    out = dWald_RDEX(dt(0, _), v[0], B[0], A[0], t0[0], s[0]);
    out = out * (1. - pTEXG_RDEX(dt(1, _), mu, sigma, tau, lb));
  }
  for (int i = 1; i < n; i++)
    out = out * (1. - pWald_RDEX(dt(i + 1, _), v[i], B[i], A[i], t0[i], s[i]));
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_rdex_old(NumericVector t, int n_acc,
                              double mu, double sigma, double tau,
                              NumericVector v, NumericVector B, NumericVector A,
                              NumericVector t0, double SSD) {
  NumericVector tmp((n_acc + 1) * t.size());
  tmp = rep_each(t, n_acc + 1) + SSD;
  NumericMatrix dt(n_acc + 1, t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dRDEXrace_old(dt, mu, sigma, tau, v, B, A, t0);
}

// [[Rcpp::export]]
NumericVector stopfn_rdex(NumericVector t, int n_acc,
                          double mu, double sigma, double tau, double lb,
                          NumericVector v, NumericVector B, NumericVector A,
                          NumericVector t0, NumericVector s, double SSD) {
  NumericVector tmp((n_acc + 1) * t.size());
  tmp = rep_each(t, n_acc + 1) + SSD;
  NumericMatrix dt(n_acc + 1, t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dRDEXrace(dt, mu, sigma, tau, lb, v, B, A, t0, s);
}

#endif
