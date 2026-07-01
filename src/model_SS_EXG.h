#ifndef ss_exg_h
#define ss_exg_h

#include <cmath>
#include <vector>
#include <cstddef>
#include <Rcpp.h>
// cens port: Zach's header pulled utility_functions.h + wald_functions.h +
// gsl_utils.h + gl_quad.h + <gsl/...>. On cens the stop-success integral runs on
// the package's own adaptive cubature (stop_quad.h -> gauss.h/hcubature.cpp), so
// none of GSL/GL is needed; the ex-Gaussian numerics and the analytic closed
// form come from the ported headers below.
#include "exgaussian_functions.h"
#include "ss_exg_analytic.h"
#include "stop_quad.h"

using namespace Rcpp;

// Parameter-matrix column layout (rows = accumulators):
//   0 muG, 1 sigG, 2 tauG, 3 muS, 4 sigS, 5 tauS, 8 exg_lb (go), 9 exgS_lb (stop)

// ----------------------------------------------------------------------------
// STOP-SUCCESS INTEGRANDS (cens hcubature callback form)
//   int f(unsigned dim, const double* x, void* p, unsigned fdim, double* retval)
// Integrate over the stop finishing time x: f_S(x) * prod_i S_go,i(x + SSD).
// ----------------------------------------------------------------------------
struct SSEXGWrapper { NumericMatrix pars; double SSD; };

// truncated ex-Gaussian stop x truncated ex-Gaussian go survivors (LIVE path)
static int ss_texg_integrand(unsigned /*dim*/, const double* x, void* p,
                             unsigned /*fdim*/, double* retval) {
  SSEXGWrapper* wp = static_cast<SSEXGWrapper*>(p);
  const double xx = x[0];
  const int n_acc = wp->pars.nrow();
  const double muS = wp->pars(0, 3), sigS = wp->pars(0, 4), tauS = wp->pars(0, 5);
  const double lbS = wp->pars(0, 9);
  double fS = dtexg(xx, muS, sigS, tauS, lbS, R_PosInf, false);
  if (fS <= 0.0) { retval[0] = 0.0; return 0; }
  double S_go_all = 1.0;
  for (int i = 0; i < n_acc; ++i) {
    double Si = ptexg(xx + wp->SSD, wp->pars(i, 0), wp->pars(i, 1),
                      wp->pars(i, 2), wp->pars(i, 8), R_PosInf, false, false);
    S_go_all *= Si;
    if (S_go_all <= 0.0) { retval[0] = 0.0; return 0; }
  }
  retval[0] = fS * S_go_all;
  return 0;
}

// regular (untruncated) ex-Gaussian stop x ex-Gaussian go survivors
static int ss_exg_integrand(unsigned /*dim*/, const double* x, void* p,
                            unsigned /*fdim*/, double* retval) {
  SSEXGWrapper* wp = static_cast<SSEXGWrapper*>(p);
  const double xx = x[0];
  const int n_acc = wp->pars.nrow();
  const double muS = wp->pars(0, 3), sigS = wp->pars(0, 4), tauS = wp->pars(0, 5);
  double fS = dexg(xx, muS, sigS, tauS, false);
  if (fS <= 0.0) { retval[0] = 0.0; return 0; }
  double S_go_all = 1.0;
  for (int i = 0; i < n_acc; ++i) {
    double Si = pexg(xx + wp->SSD, wp->pars(i, 0), wp->pars(i, 1),
                     wp->pars(i, 2), false, false);
    S_go_all *= Si;
    if (S_go_all <= 0.0) { retval[0] = 0.0; return 0; }
  }
  retval[0] = fS * S_go_all;
  return 0;
}

// ----------------------------------------------------------------------------
// TRUNCATED EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

// wrapper around truncated ex-Gaussian log PDF for race function
NumericVector texg_go_lpdf(NumericVector rt, NumericMatrix pars,
                           LogicalVector idx, double min_ll) {
  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;
  NumericVector out(n_acc_selected);
  int k = 0;
  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;
    double log_d = dtexg(
      rt[i], pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, true);
    out[k] = emc2_isfinite(log_d) ? log_d : R_NegInf;
    k++;
  }
  return(out);
}

// wrapper around truncated ex-Gaussian log complementary CDF for race function
NumericVector texg_go_lccdf(NumericVector rt, NumericMatrix pars,
                            LogicalVector idx, double min_ll) {
  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;
  NumericVector out(n_acc_selected);
  int k = 0;
  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;
    double log_s = ptexg(
      rt[i], pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, false, true);
    out[k] = emc2_isfinite(log_s) ? log_s : R_NegInf;
    k++;
  }
  return(out);
}

// go race log likelihood function, not accounting for go failure
double ss_texg_go_lpdf(double RT, NumericMatrix pars,
                       LogicalVector winner, double min_ll) {
  const int n_acc = pars.nrow();
  double logpdf_winner = 0.0;
  double logsurv_losers = 0.0;
  bool winner_found = false;
  for (int i = 0; i < n_acc; ++i) {
    if (winner[i]) {
      double ld = dtexg(RT, pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, true);
      logpdf_winner += emc2_isfinite(ld) ? ld : min_ll;
      winner_found = true;
    } else {
      double ls = ptexg(RT, pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, false, true);
      logsurv_losers += emc2_isfinite(ls) ? ls : min_ll;
    }
  }
  if (!winner_found) return min_ll;
  double out = logpdf_winner + logsurv_losers;
  return emc2_isfinite(out) ? out : min_ll;
}

// go vs stop race log likelihood for failed-inhibition (stop-process-losing)
// stop trials, not accounting for trigger failure and go failure.
double ss_texg_stop_fail_lpdf(double RT, double SSD, NumericMatrix pars,
                              LogicalVector winner, double min_ll) {
  double go_lprob = ss_texg_go_lpdf(RT, pars, winner, min_ll);
  // NB SSD subtracted from observed RT to get stop finish time; stop pars 3..5,9
  double stop_survivor_lprob = ptexg(
    RT - SSD, pars(0, 3), pars(0, 4), pars(0, 5), pars(0, 9), R_PosInf, false, true);
  if (!emc2_isfinite(stop_survivor_lprob)) stop_survivor_lprob = min_ll;
  return go_lprob + stop_survivor_lprob;
}

// stop-success integral, "integrate" route: cens adaptive hcubature over the
// finite window. (Was GSL qags/qagiu; cens has no GSL, and hcubature needs
// finite bounds, so the infinite-upper case integrates over the same heuristic
// window as the deterministic route.)
static inline double ss_texg_stop_success_lpdf(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int max_subdiv = 100,
    double abs_tol = 1e-8, double rel_tol = 1e-6,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU) {
  SSEXGWrapper w = {pars, SSD};
  const double muS = pars(0, 3), sigS = pars(0, 4), tauS = pars(0, 5);
  const double lbS = pars(0, 9);
  const double ub_heur = muS + k_sigma * sigS + k_tau * tauS;
  const double lo = emc2_isfinite(lbS) ? lbS
                                       : muS - k_sigma * sigS - k_tau * tauS;
  const double ub = emc2_isfinite(upper) ? upper : ub_heur;
  if (!(ub > lo)) return min_ll;
  const std::size_t max_eval = static_cast<std::size_t>(max_subdiv) * 1000u;
  double res = ss_adaptive_integrate(ss_texg_integrand, &w, lo, ub,
                                     abs_tol, rel_tol, max_eval);
  if (!emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// ----------------------------------------------------------------------------
// REGULAR EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

NumericVector exg_go_lpdf(NumericVector rt, NumericMatrix pars,
                          LogicalVector idx, double min_ll) {
  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;
  NumericVector out(n_acc_selected);
  int k = 0;
  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;
    double log_d = dexg(rt[i], pars(i, 0), pars(i, 1), pars(i, 2), true);
    out[k] = emc2_isfinite(log_d) ? log_d : R_NegInf;
    k++;
  }
  return(out);
}

NumericVector exg_go_lccdf(NumericVector rt, NumericMatrix pars,
                           LogicalVector idx, double min_ll) {
  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;
  NumericVector out(n_acc_selected);
  int k = 0;
  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;
    double log_s = pexg(rt[i], pars(i, 0), pars(i, 1), pars(i, 2), false, true);
    out[k] = emc2_isfinite(log_s) ? log_s : R_NegInf;
    k++;
  }
  return(out);
}

double ss_exg_go_lpdf(double RT, NumericMatrix pars,
                      LogicalVector winner, double min_ll) {
  const int n_acc = pars.nrow();
  double logpdf_winner = 0.0;
  double logsurv_losers = 0.0;
  bool winner_found = false;
  for (int i = 0; i < n_acc; ++i) {
    if (winner[i]) {
      double ld = dexg(RT, pars(i, 0), pars(i, 1), pars(i, 2), true);
      logpdf_winner += emc2_isfinite(ld) ? ld : min_ll;
      winner_found = true;
    } else {
      double ls = pexg(RT, pars(i, 0), pars(i, 1), pars(i, 2), false, true);
      logsurv_losers += emc2_isfinite(ls) ? ls : min_ll;
    }
  }
  if (!winner_found) return min_ll;
  double out = logpdf_winner + logsurv_losers;
  return emc2_isfinite(out) ? out : min_ll;
}

double ss_exg_stop_fail_lpdf(double RT, double SSD, NumericMatrix pars,
                             LogicalVector winner, double min_ll) {
  double go_lprob = ss_exg_go_lpdf(RT, pars, winner, min_ll);
  double stop_survivor_lprob = pexg(
    RT - SSD, pars(0, 3), pars(0, 4), pars(0, 5), false, true);
  if (!emc2_isfinite(stop_survivor_lprob)) stop_survivor_lprob = min_ll;
  return go_lprob + stop_survivor_lprob;
}

// regular ex-Gaussian stop-success, "integrate" route over the finite window.
static inline double ss_exg_stop_success_lpdf(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int max_subdiv = 100,
    double abs_tol = 1e-8, double rel_tol = 1e-6,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU) {
  SSEXGWrapper w = {pars, SSD};
  const double muS = pars(0, 3), sigS = pars(0, 4), tauS = pars(0, 5);
  const double halfw = k_sigma * sigS + k_tau * tauS;
  const double lo = muS - halfw;
  const double ub = emc2_isfinite(upper) ? std::min(upper, muS + halfw)
                                         : muS + halfw;
  if (!(ub > lo)) return min_ll;
  const std::size_t max_eval = static_cast<std::size_t>(max_subdiv) * 1000u;
  double res = ss_adaptive_integrate(ss_exg_integrand, &w, lo, ub,
                                     abs_tol, rel_tol, max_eval);
  if (!emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// ----------------------------------------------------------------------------
// GAUSS-KRONROD (fixed) AND DISPATCH VARIANTS OF THE STOP-SUCCESS INTEGRAL
// "gl" route: deterministic composite Gauss-Kronrod over the SAME finite window
// (stop_quad.h ss_fixed_integrate). Kept the "gl" name for API/test back-compat;
// the engine is cens's gauss_kronrod, not Gauss-Legendre.
// ----------------------------------------------------------------------------

// --- truncated ex-Gaussian stop-success, fixed (GL-labelled) route ---
static inline double ss_texg_stop_success_lpdf_gl(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int n_nodes = 64,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU) {
  SSEXGWrapper w = {pars, SSD};
  const double muS = pars(0, 3), sigS = pars(0, 4), tauS = pars(0, 5);
  const double lbS = pars(0, 9);
  const double ub_heur = muS + k_sigma * sigS + k_tau * tauS;
  const double lo = emc2_isfinite(lbS) ? lbS
                                       : muS - k_sigma * sigS - k_tau * tauS;
  const double ub = emc2_isfinite(upper) ? upper : ub_heur;
  if (!(ub > lo)) return min_ll;
  double res = ss_fixed_integrate(ss_texg_integrand, &w, lo, ub, n_nodes);
  if (!emc2_isfinite(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// --- "auto"/"analytic" dispatch ---
// n_go == 1 with an infinite upper limit: exact closed form (full-line or
// truncated; see ss_exg_analytic.h). Guard trips, kinked domains, finite upper
// limits, and n_go >= 2: fixed route with an auto panel-density bump.
static inline double ss_texg_stop_success_lpdf_autodisp(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int n_nodes = 64,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU) {
  if (pars.nrow() == 1 && !emc2_isfinite(upper)) {
    bool ok = false;
    double p = ss_texg_stop_success_analytic1(
      SSD, pars(0, 0), pars(0, 1), pars(0, 2), pars(0, 8),
      pars(0, 3), pars(0, 4), pars(0, 5), pars(0, 9), ok);
    if (ok) return std::log(p);
  }
  const double muS = pars(0, 3), sigS = pars(0, 4), tauS = pars(0, 5);
  const double lbS = pars(0, 9);
  const double lo = emc2_isfinite(lbS) ? lbS
                                       : muS - k_sigma * sigS - k_tau * tauS;
  const double ub = emc2_isfinite(upper) ? upper
                                         : muS + k_sigma * sigS + k_tau * tauS;
  const int n_eff = gl_auto_nodes(n_nodes, lo, ub, sigS);
  return ss_texg_stop_success_lpdf_gl(SSD, pars, min_ll, upper, n_eff,
                                      k_sigma, k_tau);
}

// --- LIVE entry point: reads the process-global stop_method_config() that
//     calc_ll_manager() sets from SSEXG(stop_method =, stop_n_nodes =) once per
//     likelihood call. ---
static inline double ss_texg_stop_success_lpdf_live(
    double SSD, NumericMatrix pars, double min_ll,
    double upper = R_PosInf, int max_subdiv = 100,
    double abs_tol = 1e-8, double rel_tol = 1e-6,
    double k_sigma = SS_WINDOW_K_SIGMA, double k_tau = SS_WINDOW_K_TAU) {
  const StopMethodConfig& c = stop_method_config();
  switch (c.method) {
  case STOP_METHOD_INTEGRATE:
    return ss_texg_stop_success_lpdf(SSD, pars, min_ll, upper, max_subdiv,
                                     abs_tol, rel_tol, k_sigma, k_tau);
  case STOP_METHOD_GL:
    return ss_texg_stop_success_lpdf_gl(SSD, pars, min_ll, upper, c.n_nodes,
                                        k_sigma, k_tau);
  default:   // STOP_METHOD_AUTO / STOP_METHOD_ANALYTIC
    return ss_texg_stop_success_lpdf_autodisp(SSD, pars, min_ll, upper,
                                              c.n_nodes, k_sigma, k_tau);
  }
}

// --- stop_method config setter/getter (called from calc_ll_manager and tests) ---
// [[Rcpp::export]]
void emc2_set_stop_method(std::string method = "auto", int n_nodes = 64) {
  StopMethodConfig& c = stop_method_config();
  if (method == "auto")           c.method = STOP_METHOD_AUTO;
  else if (method == "integrate") c.method = STOP_METHOD_INTEGRATE;
  else if (method == "gl")        c.method = STOP_METHOD_GL;
  else if (method == "analytic")  c.method = STOP_METHOD_ANALYTIC;
  else Rcpp::stop("unknown stop_method '" + method + "'");
  if (n_nodes < 2) Rcpp::stop("stop_n_nodes must be >= 2");
  c.n_nodes = n_nodes;
}
// [[Rcpp::export]]
Rcpp::List emc2_get_stop_method() {
  const StopMethodConfig& c = stop_method_config();
  static const char* method_names[4] = {"auto", "integrate", "gl", "analytic"};
  return Rcpp::List::create(
    Rcpp::_["method"]  = std::string(method_names[c.method]),
    Rcpp::_["n_nodes"] = c.n_nodes);
}

// --- which branch would "auto" take for this trial? (test/benchmark aid) ---
// [[Rcpp::export]]
std::string ss_texg_stop_success_auto_branch(
    double SSD, NumericMatrix pars, double upper = -1.0) {
  if (upper <= 0.0) upper = R_PosInf;   // sentinel: <=0 means auto/Inf
  if (pars.nrow() != 1) return "gl_ngo2plus";
  if (emc2_isfinite(upper)) return "gl_finite_upper";
  bool ok = false;
  ss_texg_stop_success_analytic1(
    SSD, pars(0, 0), pars(0, 1), pars(0, 2), pars(0, 8),
    pars(0, 3), pars(0, 4), pars(0, 5), pars(0, 9), ok);
  if (!ok) return "gl_guard";
  const bool full_line = (pars(0, 9) == R_NegInf) && (pars(0, 8) == R_NegInf);
  return full_line ? "analytic_fullline" : "analytic_trunc";
}

// --- Rcpp test wrapper: returns the INTEGRAL VALUE (not log) for one trial.
//     method = "integrate" -> adaptive hcubature (config-independent),
//              "gl"        -> fixed composite Gauss-Kronrod with n_nodes,
//              "analytic"/"auto" -> the auto dispatch rule,
//              "live"      -> the live dispatcher (honours emc2_set_stop_method).
// NB k_sigma/k_tau defaults are literal copies of SS_WINDOW_K_SIGMA/K_TAU —
// compileAttributes copies defaults verbatim into the generated R wrapper.
// [[Rcpp::export]]
double ss_texg_stop_success_value(
    double SSD, NumericMatrix pars, std::string method = "integrate",
    double upper = -1.0, int n_nodes = 64,
    double k_sigma = 8.0, double k_tau = 16.0,
    int max_subdiv = 100, double abs_tol = 1e-8, double rel_tol = 1e-6) {
  const double min_ll = -1e10;
  if (upper <= 0.0) upper = R_PosInf;   // sentinel: <=0 means auto/Inf
  double lp;
  if (method == "gl") {
    lp = ss_texg_stop_success_lpdf_gl(SSD, pars, min_ll, upper, n_nodes,
                                      k_sigma, k_tau);
  } else if (method == "auto" || method == "analytic") {
    lp = ss_texg_stop_success_lpdf_autodisp(SSD, pars, min_ll, upper, n_nodes,
                                            k_sigma, k_tau);
  } else if (method == "live") {
    lp = ss_texg_stop_success_lpdf_live(SSD, pars, min_ll, upper, max_subdiv,
                                        abs_tol, rel_tol, k_sigma, k_tau);
  } else {
    lp = ss_texg_stop_success_lpdf(SSD, pars, min_ll, upper, max_subdiv,
                                   abs_tol, rel_tol, k_sigma, k_tau);
  }
  return (lp <= min_ll) ? 0.0 : std::exp(lp);
}

// ----------------------------------------------------------------------------
// R-FACING VECTOR HELPERS (used by the R stop_success / plotting routes)
// ----------------------------------------------------------------------------

NumericVector dexg_c(const NumericVector x, const double mu = 5.,
                     const double sigma = 1., const double tau = 1.,
                     const bool log_d = false) {
  int n = x.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = dexg(x[i], mu, sigma, tau, log_d);
  return(out);
}

NumericVector pexg_c(const NumericVector q, const double mu = 5.,
                     const double sigma = 1., const double tau = 1.,
                     const bool lower_tail = true, const bool log_p = false) {
  int n = q.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = pexg(q[i], mu, sigma, tau, lower_tail, log_p);
  return(out);
}

NumericVector protect_finite(const NumericVector& x, double min_ll) {
  NumericVector out = clone(x);
  for (int i = 0; i < out.size(); i++) if (!R_FINITE(out[i])) out[i] = min_ll;
  return out;
}

// [[Rcpp::export]]
NumericVector dEXGrace(NumericMatrix dt, NumericVector mu, NumericVector sigma,
                       NumericVector tau, double min_ll) {
  int n = mu.size();
  NumericVector log_out(dt.ncol());
  log_out = protect_finite(dexg_c(dt(0, _), mu[0], sigma[0], tau[0], true), min_ll);
  for (int i = 1; i < n; i++) {
    log_out += protect_finite(pexg_c(dt(i, _), mu[i], sigma[i], tau[i], false, true), min_ll);
  }
  return exp(log_out);
}

// [[Rcpp::export]]
NumericVector stopfn_exg(NumericVector t, NumericVector mu, NumericVector sigma,
                         NumericVector tau, double SSD, double min_ll) {
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dEXGrace(dt, mu, sigma, tau, min_ll);
}

// [[Rcpp::export]]
NumericVector pTEXG_vec(NumericVector q, double mu = 5., double sigma = 1.,
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
NumericVector dTEXG_vec(NumericVector x, double mu = 5., double sigma = 1.,
                        double tau = 1., double lb = .05, bool log_d = false) {
  int n = x.size();
  if (tau <= 0. || sigma <= 0.) { NumericVector pdf(n, NA_REAL); return pdf; }
  NumericVector pdf(n);
  for (int i = 0; i < n; i++)
    pdf[i] = dtexg(x[i], mu, sigma, tau, lb, R_PosInf, log_d);
  return pdf;
}

// [[Rcpp::export]]
NumericVector dTEXGrace(NumericMatrix dt, NumericVector mu, NumericVector sigma,
                        NumericVector tau, NumericVector lb) {
  int n = mu.size();
  int n_points = dt.ncol();
  NumericVector out(n_points);
  out = dTEXG_vec(dt(0, _), mu[0], sigma[0], tau[0], lb[0]);
  for (int i = 1; i < n; i++)
    out = out * pTEXG_vec(dt(i, _), mu[i], sigma[i], tau[i], lb[i], false);
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_texg(NumericVector t, NumericVector mu, NumericVector sigma,
                          NumericVector tau, NumericVector lb, double SSD) {
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dTEXGrace(dt, mu, sigma, tau, lb);
}

#endif
