#ifndef ss_likelihood_h
#define ss_likelihood_h

// ---------------------------------------------------------------------------
// Stop-signal trial-loop likelihood (ported from Zach's particle_ll.cpp).
//
// Computes the per-trial SS log-likelihood: go trials, stop-fail (signal-
// respond), stop-success (withheld), stop-triggered (ST) accumulators, trigger
// failure (tf), go failure (gf), deadline censoring (UC), and lower censoring
// (LC). The expensive stop-success integral is delegated to the model headers'
// _live dispatchers (cens hcubature, set via emc2_set_stop_method). The one
// remaining GSL integral on Zach's branch (the lower-censoring response-mass
// integral) is converted here to cens's adaptive hcubature (ss_adaptive_integrate).
//
// Must be included AFTER model_SS_EXG.h / model_SS_RDEX.h (function pointers
// reference texg_go_lpdf / rdex_go_lpdf / ss_*_stop_success_lpdf_live) and after
// utility_functions.h (submat_rcpp).
// ---------------------------------------------------------------------------

#include <cmath>
#include <string>
#include <Rcpp.h>
#include "ss_compat.h"               // log1m, log1m_exp, log_sum_exp, log_mix, emc2_isfinite
#include "stop_quad.h"               // ss_adaptive_integrate, SS_WINDOW_K_SIGMA/K_TAU
#include "exgaussian_functions.h"    // ptexg (stop survivor)
#include "model_SS_EXG.h"
#include "model_SS_RDEX.h"
#include "utility_functions.h"       // submat_rcpp(NumericMatrix, LogicalVector)

using namespace Rcpp;

// --- SS helper pointer types ------------------------------------------------
using ss_go_pdf_fn       = NumericVector (*)(NumericVector, NumericMatrix, LogicalVector, double);
using ss_stop_surv_fn    = double (*)(double, NumericMatrix);
using ss_stop_success_fn = double (*)(double, NumericMatrix, double, double, int, double, double, double, double);

// Model-specific stop-survivor wrappers (read the fixed stop columns).
static inline double stop_logsurv_texg_fn(double q, NumericMatrix P) {
  // EXG stop: muS=3, sigmaS=4, tauS=5, exgS_lb=9
  return ptexg(q, P(0, 3), P(0, 4), P(0, 5), P(0, 9), R_PosInf, false, true);
}
static inline double stop_logsurv_rdex_fn(double q, NumericMatrix P) {
  // RDEX stop: muS=5, sigmaS=6, tauS=7, exgS_lb=10
  return ptexg(q, P(0, 5), P(0, 6), P(0, 7), P(0, 10), R_PosInf, false, true);
}

// Fetch a data column, substituting default_val for absent column or NAs.
inline NumericVector ss_get_col_with_default(const DataFrame& df,
                                             const std::string& name,
                                             double default_val) {
  if (df.containsElementNamed(name.c_str())) {
    NumericVector col = df[name];
    bool has_na = false;
    for (int i = 0; i < col.size(); ++i)
      if (NumericVector::is_na(col[i])) { has_na = true; break; }
    if (!has_na) return col;
    NumericVector out = clone(col);
    for (int i = 0; i < out.size(); ++i)
      if (NumericVector::is_na(out[i])) out[i] = default_val;
    return out;
  }
  return NumericVector(df.nrows(), default_val);
}

static inline double sum_log_terms(const NumericVector& x) {
  double out = 0.0;
  for (int i = 0; i < x.size(); ++i) out += x[i];
  return out;
}

static inline double ss_log_go_density_for_winner(
    double rt, const NumericMatrix& P, const LogicalVector& is_go, int winner_idx,
    ss_go_pdf_fn go_lpdf_ptr, ss_go_pdf_fn go_lccdf_ptr, double min_ll) {
  const int n_acc = P.nrow();
  NumericVector rt_go(n_acc, rt);
  LogicalVector go_win_mask(n_acc, false);
  LogicalVector go_loss_mask(n_acc, false);
  int go_loss_count = 0;
  for (int i = 0; i < n_acc; ++i) {
    const bool go_i = (is_go[i] == TRUE);
    if (!go_i) continue;
    if (i == winner_idx) go_win_mask[i] = true;
    else { go_loss_mask[i] = true; ++go_loss_count; }
  }
  NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
  double out = (lw.size() > 0) ? sum_log_terms(lw) : R_NegInf;
  if (!R_FINITE(out)) return R_NegInf;
  if (go_loss_count == 0) return out;
  NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
  return out + sum_log_terms(ls);
}

static inline double ss_log_st_density_for_winner(
    double rt, double SSD, const NumericMatrix& P, const LogicalVector& is_st,
    int winner_idx, ss_go_pdf_fn go_lpdf_ptr, ss_go_pdf_fn go_lccdf_ptr, double min_ll) {
  const int n_acc = P.nrow();
  double rt_st_val = rt - SSD;
  if (rt_st_val < 0.0) rt_st_val = 0.0;
  NumericVector rt_st(n_acc, rt_st_val);
  LogicalVector st_win_mask(n_acc, false);
  LogicalVector st_loss_mask(n_acc, false);
  int st_loss_count = 0;
  for (int i = 0; i < n_acc; ++i) {
    const bool st_i = (is_st[i] == TRUE);
    if (!st_i) continue;
    if (i == winner_idx) st_win_mask[i] = true;
    else { st_loss_mask[i] = true; ++st_loss_count; }
  }
  NumericVector lw = go_lpdf_ptr(rt_st, P, st_win_mask, min_ll);
  double out = (lw.size() > 0) ? sum_log_terms(lw) : R_NegInf;
  if (!R_FINITE(out)) return R_NegInf;
  if (st_loss_count == 0) return out;
  NumericVector ls = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
  return out + sum_log_terms(ls);
}

static inline double ss_trial_log_response_density(
    double rt, const NumericMatrix& P, const NumericVector& lR_trial,
    const LogicalVector& is_go, const LogicalVector& is_st, double SSD,
    bool stop_signal_presented, double tf, double gf,
    ss_go_pdf_fn go_lpdf_ptr, ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr, ss_stop_success_fn stop_success_ptr,
    double min_ll, int response_code = NA_INTEGER) {
  const int n_acc = P.nrow();
  const bool response_known = response_code != NA_INTEGER;
  int n_accG = 0, n_accST = 0;
  for (int i = 0; i < n_acc; ++i) {
    if (is_go[i] == TRUE) ++n_accG;
    if (is_st[i] == TRUE) ++n_accST;
  }

  double log_go_any = R_NegInf;
  for (int i = 0; i < n_acc; ++i) {
    if (is_go[i] != TRUE) continue;
    if (response_known && lR_trial[i] != response_code) continue;
    const double go_lprob = ss_log_go_density_for_winner(
      rt, P, is_go, i, go_lpdf_ptr, go_lccdf_ptr, min_ll);
    if (!R_FINITE(go_lprob)) continue;

    double log_term = log1m(gf) + go_lprob;
    if (stop_signal_presented && rt > SSD) {
      double rt_eff = rt - SSD;
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = stop_logsurv_ptr(rt_eff, P);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = min_ll;

      double st_loss_sum = 0.0;
      if (n_accST > 0) {
        NumericVector rt_st(n_acc, rt_eff);
        LogicalVector st_loss_mask(n_acc, false);
        int st_loss_count = 0;
        for (int j = 0; j < n_acc; ++j) {
          st_loss_mask[j] = (is_st[j] == TRUE);
          if (st_loss_mask[j] == TRUE) ++st_loss_count;
        }
        if (st_loss_count > 0) {
          NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
          st_loss_sum = sum_log_terms(ls_st);
        }
      }
      const double comp_tf = go_lprob;
      const double comp_notf = go_lprob + log_stop_surv + st_loss_sum;
      log_term = log1m(gf) + log_mix(tf, comp_tf, comp_notf);
    }
    log_go_any = log_sum_exp(log_go_any, log_term);
  }

  if (!stop_signal_presented || n_accST == 0) return log_go_any;

  double log_st_any = R_NegInf;
  NumericMatrix P_go = submat_rcpp(P, is_go);
  for (int i = 0; i < n_acc; ++i) {
    if (is_st[i] != TRUE) continue;
    if (response_known && lR_trial[i] != response_code) continue;
    if (rt <= SSD) continue;
    const double st_base = ss_log_st_density_for_winner(
      rt, SSD, P, is_st, i, go_lpdf_ptr, go_lccdf_ptr, min_ll);
    if (!R_FINITE(st_base)) continue;

    double go_loss_sum = 0.0;
    if (n_accG > 0) {
      NumericVector rt_go(n_acc, rt);
      NumericVector ls_go = go_lccdf_ptr(rt_go, P, is_go, min_ll);
      go_loss_sum = sum_log_terms(ls_go);
    }
    double rt_eff = rt - SSD;
    if (rt_eff < 0.0) rt_eff = 0.0;
    double log_pstop = stop_success_ptr(SSD, P_go, min_ll, rt_eff,
                                        100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
    if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;

    const double term_gf = std::log(gf) + st_base;
    const double term_stop_win = log1m(gf) + log_pstop + st_base;
    const double term_stop_lose = log1m(gf) + log1m_exp(log_pstop) + st_base + go_loss_sum;
    const double log_term =
      log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));
    log_st_any = log_sum_exp(log_st_any, log_term);
  }
  return log_sum_exp(log_go_any, log_st_any);
}

// Lower-censoring response-mass integral. cens port: was a GSL qags call on
// Zach's branch; here it integrates the response density over [lower, upper]
// with cens's adaptive hcubature (stop_quad.h).
struct SsLcIntegrandData {
  NumericMatrix P;
  NumericVector lR_trial;
  LogicalVector is_go;
  LogicalVector is_st;
  double SSD;
  bool stop_signal_presented;
  double tf;
  double gf;
  ss_go_pdf_fn go_lpdf_ptr;
  ss_go_pdf_fn go_lccdf_ptr;
  ss_stop_surv_fn stop_logsurv_ptr;
  ss_stop_success_fn stop_success_ptr;
  double min_ll;
  int response_code;
};

static int ss_lc_integrand(unsigned /*dim*/, const double* x, void* params,
                           unsigned /*fdim*/, double* retval) {
  SsLcIntegrandData* d = static_cast<SsLcIntegrandData*>(params);
  const double logf = ss_trial_log_response_density(
    x[0], d->P, d->lR_trial, d->is_go, d->is_st, d->SSD, d->stop_signal_presented,
    d->tf, d->gf, d->go_lpdf_ptr, d->go_lccdf_ptr, d->stop_logsurv_ptr,
    d->stop_success_ptr, d->min_ll, d->response_code);
  retval[0] = R_FINITE(logf) ? std::exp(logf) : 0.0;
  return 0;
}

static inline double ss_integrate_lc_response_mass(
    double lower, double upper, const NumericMatrix& P, const NumericVector& lR_trial,
    const LogicalVector& is_go, const LogicalVector& is_st, double SSD,
    bool stop_signal_presented, double tf, double gf,
    ss_go_pdf_fn go_lpdf_ptr, ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr, ss_stop_success_fn stop_success_ptr,
    double min_ll, int response_code = NA_INTEGER) {
  if (!(upper > lower)) return min_ll;
  SsLcIntegrandData data{
    P, lR_trial, is_go, is_st, SSD, stop_signal_presented, tf, gf,
    go_lpdf_ptr, go_lccdf_ptr, stop_logsurv_ptr, stop_success_ptr, min_ll, response_code
  };
  double res = ss_adaptive_integrate(ss_lc_integrand, &data, lower, upper,
                                     1e-8, 1e-6, /*max_eval=*/200000u);
  if (!R_FINITE(res) || res <= 0.0) return min_ll;
  return std::log(res);
}

// ---------------------------------------------------------------------------
// Main SS trial-loop likelihood. pars: one row per accumulator per trial.
// ---------------------------------------------------------------------------
inline double c_log_likelihood_ss(
    NumericMatrix pars, DataFrame data, const int n_trials,
    IntegerVector expand, double min_ll, const std::vector<int>& is_ok,
    ss_go_pdf_fn go_lpdf_ptr, ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr, ss_stop_success_fn stop_success_ptr,
    int idx_tf, int idx_gf) {

  const int n_out = expand.length();
  bool any_ok = false;
  for (size_t i = 0; i < is_ok.size(); ++i) if (is_ok[i]) { any_ok = true; break; }
  if (!any_ok) return static_cast<double>(n_out) * min_ll;

  NumericVector lls(n_trials);
  NumericVector RT = data["rt"];
  IntegerVector R = data["R"];
  NumericVector SSD = data["SSD"];
  NumericVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  NumericVector LT = ss_get_col_with_default(data, "LT", 0.0);
  NumericVector UC = ss_get_col_with_default(data, "UC", R_PosInf);
  NumericVector LC = ss_get_col_with_default(data, "LC", 0.0);
  bool has_lI = data.containsElementNamed("lI");
  IntegerVector lI = has_lI ? as<IntegerVector>(data["lI"]) : IntegerVector(lR.size(), 2);
  // cens censoring convention: a per-trial `missingness` code identifies every
  // non-observed trial (1=lower, 2=upper/deadline, 3=both, 4=withheld/intrinsic
  // NR; NA=observed). rt is finite (observed) or NA — never Inf/-Inf. Absent
  // column => treat as observed, with a fallback below that maps NA-rt to code 4.
  bool has_missingness = data.containsElementNamed("missingness");
  IntegerVector missingness = has_missingness ? as<IntegerVector>(data["missingness"])
                                              : IntegerVector(lR.size(), NA_INTEGER);

  NumericVector unique_lR = unique(lR);
  const int n_acc = unique_lR.length();

  NumericVector tt(n_acc);
  auto log_surv_mask = [&](double t, const NumericMatrix& Pcur,
                           const LogicalVector& mask) -> double {
    tt.fill(t);
    NumericVector ls = go_lccdf_ptr(tt, Pcur, mask, min_ll);
    double out = 0.0;
    for (int i = 0; i < ls.size(); ++i) out += ls[i];
    return out;
  };

  for (int trial = 0; trial < n_trials; trial++) {
    if (is_ok[trial * n_acc] != 1) { lls[trial] = min_ll; continue; }
    int start_row = trial * n_acc;
    int end_row   = (trial + 1) * n_acc - 1;
    NumericMatrix P = pars(Range(start_row, end_row), _);
    NumericVector lR_trial = lR[Range(start_row, end_row)];
    IntegerVector lI_trial = lI[Range(start_row, end_row)];
    LogicalVector is_go(n_acc, false), is_st(n_acc, false);
    int n_accG = 0, n_accST = 0;
    if (has_lI) {
      int go_code = max(lI_trial);
      for (int i = 0; i < n_acc; i++) {
        const bool go_i = (lI_trial[i] == go_code);
        is_go[i] = go_i; is_st[i] = !go_i;
        if (go_i) ++n_accG; else ++n_accST;
      }
    } else {
      for (int i = 0; i < n_acc; i++) { is_go[i] = true; ++n_accG; }
    }

    double tf = P(0, idx_tf);
    double gf = P(0, idx_gf);
    double rt = RT[start_row];
    bool response_observed = R[start_row] != NA_INTEGER;
    bool stop_signal_presented = emc2_isfinite(SSD[start_row]);
    double uc = UC[start_row];
    bool response_is_go = false;
    if (response_observed) {
      int r_obs = R[start_row];
      for (int i = 0; i < n_acc; i++) {
        if (lR[start_row + i] == r_obs) { response_is_go = is_go[i]; break; }
      }
    }

    NumericVector rt_go(n_acc, rt);
    double rt_st_val = rt - SSD[start_row];
    if (rt_st_val < 0.0) rt_st_val = 0.0;
    NumericVector rt_st(n_acc, rt_st_val);

    LogicalVector win_mask = winner[Range(start_row, end_row)];
    LogicalVector go_win_mask(n_acc), go_loss_mask(n_acc);
    for (int i = 0; i < n_acc; i++) {
      go_win_mask[i] = (win_mask[i] && is_go[i]);
      go_loss_mask[i] = (!win_mask[i] && is_go[i]);
    }

    // --- classify the trial by its missingness code (cens convention) ---
    int miss = has_missingness ? missingness[start_row] : NA_INTEGER;
    // Fallback for data lacking an explicit code: a non-finite (NA) rt is an
    // intrinsic no-response (withheld), never an observed trial.
    if (miss == NA_INTEGER && !R_FINITE(rt)) miss = 4;

    // upper-censored / deadline log-mass, using UC[start_row] as the deadline.
    auto upper_deadline_lp = [&]() -> double {
      double logS_go = (n_accG > 0) ? log_surv_mask(uc, P, is_go) : 0.0;
      if (!stop_signal_presented)
        return log_sum_exp(std::log(gf), log1m(gf) + logS_go);
      double uc_eff = uc - SSD[start_row];
      if (!R_FINITE(uc_eff) || uc_eff <= 0.0) uc_eff = 0.0;
      bool stop_can_act = R_FINITE(uc_eff) && (uc_eff > 0.0);
      double log_pstop = R_NegInf, logS_stop = 0.0;
      NumericMatrix P_go = submat_rcpp(P, is_go);
      if (stop_can_act) {
        logS_stop = stop_logsurv_ptr(uc_eff, P);
        log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, uc_eff, 100, 1e-8, 1e-6,
                                     SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
        if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
      }
      double log_core_trig = log_sum_exp(log_pstop, logS_go + logS_stop);
      if (n_accST == 0) {
        double log_no_nogf = log_mix(tf, logS_go, log_core_trig);
        return log_sum_exp(std::log(gf), log1m(gf) + log_no_nogf);
      }
      double logS_st = log_surv_mask(uc_eff, P, is_st);
      double log_trig = logS_st + log_sum_exp(std::log(gf), log1m(gf) + log_core_trig);
      double log_tfbranch = log_sum_exp(std::log(gf), log1m(gf) + logS_go);
      return log_mix(tf, log_tfbranch, log_trig);
    };

    // lower-censored log-mass over [LT, LC].
    auto lower_censor_lp = [&]() -> double {
      return ss_integrate_lc_response_mass(
        LT[start_row], LC[start_row], P, lR_trial, is_go, is_st, SSD[start_row],
        stop_signal_presented, tf, gf, go_lpdf_ptr, go_lccdf_ptr, stop_logsurv_ptr,
        stop_success_ptr, min_ll, response_observed ? R[start_row] : NA_INTEGER);
    };

    // withheld / intrinsic no-response (stop-success or go-failure).
    auto withheld_lp = [&]() -> double {
      if (!stop_signal_presented) return std::log(gf);   // go failure
      if (n_accST == 0) {
        NumericMatrix P_go = submat_rcpp(P, is_go);
        double log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, R_PosInf,
                                            100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
        if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
        return log_sum_exp(std::log(gf), log1m(gf) + log1m(tf) + log_pstop);  // stop-success
      }
      return std::log(gf) + std::log(tf);                // ST intrinsic NR (out of scope)
    };

    if (miss == 1) { lls[trial] = lower_censor_lp(); continue; }
    if (miss == 2) { lls[trial] = upper_deadline_lp(); continue; }
    if (miss == 3) { lls[trial] = log_sum_exp(lower_censor_lp(), upper_deadline_lp()); continue; }
    if (miss == 4) { lls[trial] = withheld_lp(); continue; }
    // miss == NA: observed response (rt finite). Guard against inconsistent data.
    if (!R_FINITE(rt)) { lls[trial] = min_ll; continue; }

    // Response observed
    if (!stop_signal_presented) {
      double go_lprob = 0.0;
      NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
      go_lprob = (lw.size() > 0) ? sum(lw) : R_NegInf;
      if (!R_FINITE(go_lprob)) go_lprob = R_NegInf;
      if (n_accG > 1) {
        NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
        for (int i = 0; i < ls.size(); i++) go_lprob += ls[i];
      }
      lls[trial] = log1m(gf) + go_lprob;
      continue;
    }

    if (response_is_go) {
      double go_lprob = 0.0;
      NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
      go_lprob = (lw.size() > 0) ? sum(lw) : R_NegInf;
      if (!R_FINITE(go_lprob)) go_lprob = R_NegInf;
      if (n_accG > 1) {
        NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
        for (int i = 0; i < ls.size(); i++) go_lprob += ls[i];
      }
      double rt_eff = rt - SSD[start_row];
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_stop_surv = stop_logsurv_ptr(rt_eff, P);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = min_ll;
      double st_loss_sum = 0.0;
      if (n_accST > 0) {
        LogicalVector st_loss_mask(n_acc);
        for (int i = 0; i < n_acc; i++) st_loss_mask[i] = (is_st[i] && !win_mask[i]);
        NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
        for (int i = 0; i < ls_st.size(); i++) st_loss_sum += ls_st[i];
      }
      double comp_tf = go_lprob;
      double comp_notf = go_lprob + log_stop_surv + st_loss_sum;
      lls[trial] = log1m(gf) + log_mix(tf, comp_tf, comp_notf);
      continue;
    } else {
      LogicalVector st_win_mask(n_acc);
      for (int i = 0; i < n_acc; i++) st_win_mask[i] = (win_mask[i] && is_st[i]);
      NumericVector lw_st = go_lpdf_ptr(rt_st, P, st_win_mask, min_ll);
      double st_winner_logpdf = (lw_st.size() > 0) ? sum(lw_st) : R_NegInf;
      if (!R_FINITE(st_winner_logpdf)) st_winner_logpdf = R_NegInf;
      double st_loss_sum = 0.0;
      if (n_accST > 1) {
        LogicalVector st_loss_mask(n_acc);
        for (int i = 0; i < n_acc; i++) st_loss_mask[i] = (!win_mask[i] && is_st[i]);
        NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
        for (int i = 0; i < ls_st.size(); i++) st_loss_sum += ls_st[i];
      }
      double go_loss_sum = 0.0;
      if (n_accG > 0) {
        NumericVector ls_go = go_lccdf_ptr(rt_go, P, is_go, min_ll);
        for (int i = 0; i < ls_go.size(); i++) go_loss_sum += ls_go[i];
      }
      NumericMatrix P_go = submat_rcpp(P, is_go);
      double rt_eff = rt - SSD[start_row];
      if (rt_eff < 0.0) rt_eff = 0.0;
      double log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, rt_eff,
                                          100, 1e-8, 1e-6, SS_WINDOW_K_SIGMA, SS_WINDOW_K_TAU);
      if (!R_FINITE(log_pstop)) log_pstop = R_NegInf;
      double st_base = st_winner_logpdf + st_loss_sum;
      double term_gf = std::log(gf) + st_base;
      double term_stop_win = log1m(gf) + log_pstop + st_base;
      double term_stop_lose = log1m(gf) + log1m_exp(log_pstop) + st_base + go_loss_sum;
      lls[trial] = log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));
      continue;
    }
  }

  lls[is_na(lls)] = min_ll;
  lls[is_infinite(lls)] = min_ll;
  lls[lls < min_ll] = min_ll;
  double sum_ll = 0.0;
  for (int i = 0; i < n_out; ++i) sum_ll += lls[expand[i] - 1];  // expand is 1-based
  return sum_ll;
}

// Resolve function pointers + fixed column indices for SS (EXG or RDEX).
struct SSModelAdapter {
  ss_go_pdf_fn    go_lpdf_ptr;
  ss_go_pdf_fn    go_lccdf_ptr;
  ss_stop_surv_fn stop_logsurv_ptr;
  ss_stop_success_fn stop_success_ptr;
  int idx_tf;
  int idx_gf;
};

static inline SSModelAdapter resolve_ss_adapter(const std::string& type_std) {
  SSModelAdapter a;
  const bool is_exg = type_std.find("EXG") != std::string::npos;
  a.go_lpdf_ptr      = is_exg ? texg_go_lpdf  : rdex_go_lpdf;
  a.go_lccdf_ptr     = is_exg ? texg_go_lccdf : rdex_go_lccdf;
  a.stop_logsurv_ptr = is_exg ? stop_logsurv_texg_fn : stop_logsurv_rdex_fn;
  a.stop_success_ptr = is_exg ? ss_texg_stop_success_lpdf_live : ss_rdex_stop_success_lpdf_live;
  a.idx_tf           = is_exg ? 6 : 8;
  a.idx_gf           = is_exg ? 7 : 9;
  return a;
}

#endif // ss_likelihood_h
