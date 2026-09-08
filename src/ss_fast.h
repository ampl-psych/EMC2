#ifndef ss_fast_h
#define ss_fast_h

// ---------------------------------------------------------------------------
// Stop-signal likelihood (SSEXG / SSRDEX) in the race-model style:
//
//  * SSSpec: everything that depends on the DATA only (trial kinds, row masks,
//    winner / response indices, SSD, censoring bounds) is computed once per
//    dadm by make_ss_spec() at the R boundary (the only place Rcpp is used).
//  * ss_log_likelihood(): the per-particle work reads parameters straight from
//    the ParamTable columns and uses only plain C++ (no Rcpp objects, no R
//    API), so it is thread-safe: calc_ll_multithreaded runs it on per-thread
//    copies (SSSpec::rebind, as CensorSpec does).
//  * The stop-success integral (P(stop wins) over stop finish times) is
//    memoised per particle keyed on its inputs (SSD, upper bound, go/stop
//    parameter values) so trials sharing a design cell and SSD - most of them
//    without trends - integrate once.
//
// The likelihood terms (go / signal-respond / stop-triggered responses, the
// no-response and deadline masses, lower-censoring response mass) mirror the
// R reference log_likelihood_race_ss() in R/model_SS.R.
// ---------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <Rcpp.h>
#include "ParamTable.h"
#include "nan_check.h"            // is_finite
#include "r_constants.h"          // pos_inf, neg_inf
#include "composite_functions.h"  // log1m, log1m_exp, log_sum_exp, log_mix
#include "exgaussian_functions.h" // dtexg, ptexg
#include "wald_functions.h"       // digt, pigt
#include "ss_integrate.h"         // ss_integrate (hcubature), window helpers

// integration accuracy for the stop-success and lower-censoring integrals
// (matches the R reference: rel 1e-6 / abs 1e-8)
constexpr double SSF_ABS_TOL = 1e-8;
constexpr double SSF_REL_TOL = 1e-6;
constexpr std::size_t SSF_STOP_MAX_EVAL = 100u * 64u;
constexpr std::size_t SSF_LC_MAX_EVAL   = 200000u;

enum class SSKind : std::int8_t {
  ObsGo = 0,        // go trial, response observed
  ObsStopFailGo,    // stop trial, a go accumulator responded (signal-respond)
  ObsST,            // stop trial, a stop-triggered accumulator responded
  Withheld,         // no response, no deadline (intrinsic NR)
  UpperCens,        // missingness 2 with finite UC (no response by the deadline)
  LowerCens,        // missingness 1
  BothCens,         // missingness 3
  Invalid           // e.g. code NA with non-finite rt (should not occur)
};

struct SSSpec {
  bool is_exg = true;
  // parameter columns in the ParamTable (resolved once by name)
  int col_mu = -1, col_sigma = -1, col_tau = -1, col_lb = -1;   // SSEXG go
  int col_v = -1, col_B = -1, col_A = -1, col_t0 = -1, col_s = -1; // SSRDEX go
  int col_muS = -1, col_sigmaS = -1, col_tauS = -1, col_lbS = -1;
  int col_tf = -1, col_gf = -1;

  int n_trials = 0, n_rows = 0, n_acc = 0;

  // per trial (data only)
  std::vector<SSKind>  kind;
  std::vector<int>     base_row, n_accG, n_accST, winner_local, resp_local;
  std::vector<double>  rt, ssd, rt_eff, lt, lc, uc, uc_eff;
  std::vector<std::uint8_t> stop_presented, stop_can_act;
  // per row (data only)
  std::vector<std::uint8_t> is_go, is_st;

  // per-particle memo of the stop-success integral (see header comment):
  // hash -> entry; the full key is kept in memo_keys for verification
  struct MemoEntry { std::size_t off; std::size_t len; double val; };
  mutable std::unordered_map<std::uint64_t, MemoEntry> memo;
  mutable std::vector<double> memo_keys;

  const ParamTable* pt = nullptr;
  void rebind(const ParamTable& p) { pt = &p; }
};

// ---------------------------------------------------------------------------
// Parameter access (all inline; the branch on is_exg is uniform per dadm)
// ---------------------------------------------------------------------------
namespace ssf {

inline double par(const ParamTable& pt, int row, int col) { return pt.base(row, col); }

// log density of a go accumulator (row r) finishing at time t
inline double go_lpdf(const SSSpec& s, const ParamTable& pt, int r, double t) {
  if (s.is_exg) {
    double v = dtexg(t, par(pt, r, s.col_mu), par(pt, r, s.col_sigma), par(pt, r, s.col_tau),
                     par(pt, r, s.col_lb), pos_inf(), true);
    return is_finite(v) ? v : neg_inf();
  }
  const double sc = par(pt, r, s.col_s);
  const double dt = t - par(pt, r, s.col_t0);
  if (!(dt > 0.0)) return neg_inf();
  const double A = par(pt, r, s.col_A) / sc;
  double d = digt(dt, par(pt, r, s.col_B) / sc + 0.5 * A, par(pt, r, s.col_v) / sc, 0.5 * A);
  double v = std::log(d);
  return is_finite(v) ? v : neg_inf();
}

// log survivor of a go accumulator (row r) at time t
inline double go_lsurv(const SSSpec& s, const ParamTable& pt, int r, double t) {
  if (s.is_exg) {
    double v = ptexg(t, par(pt, r, s.col_mu), par(pt, r, s.col_sigma), par(pt, r, s.col_tau),
                     par(pt, r, s.col_lb), pos_inf(), false, true);
    return is_finite(v) ? v : neg_inf();
  }
  const double sc = par(pt, r, s.col_s);
  const double dt = t - par(pt, r, s.col_t0);
  if (!(dt > 0.0)) return 0.0;
  const double A = par(pt, r, s.col_A) / sc;
  double v = log1m(pigt(dt, par(pt, r, s.col_B) / sc + 0.5 * A, par(pt, r, s.col_v) / sc, 0.5 * A));
  return is_finite(v) ? v : neg_inf();
}

// log survivor of the stop process (parameters on row r0) at stop-relative time q
inline double stop_lsurv(const SSSpec& s, const ParamTable& pt, int r0, double q) {
  double v = ptexg(q, par(pt, r0, s.col_muS), par(pt, r0, s.col_sigmaS), par(pt, r0, s.col_tauS),
                   par(pt, r0, s.col_lbS), pos_inf(), false, true);
  return is_finite(v) ? v : neg_inf();
}

// ---------------------------------------------------------------------------
// Stop-success integral: P(stop finishes before every go accumulator | SSD),
// over stop-relative finish times in [lb, upper] clipped to the finite window.
// ---------------------------------------------------------------------------
struct StopIntegrand {
  const SSSpec* s; const ParamTable* pt;
  double ssd, muS, sigS, tauS, lbS;
  const int* go_rows; int n_go;
};

inline int stop_success_integrand(unsigned, const double* x, void* p, unsigned, double* out) {
  const StopIntegrand* w = static_cast<const StopIntegrand*>(p);
  const double xx = x[0];
  double lf = dtexg(xx, w->muS, w->sigS, w->tauS, w->lbS, pos_inf(), true);
  if (!is_finite(lf)) { out[0] = 0.0; return 0; }
  double ls = 0.0;
  for (int i = 0; i < w->n_go; ++i) {
    double v = go_lsurv(*w->s, *w->pt, w->go_rows[i], xx + w->ssd);
    if (!is_finite(v)) { out[0] = 0.0; return 0; }
    ls += v;
  }
  out[0] = std::exp(lf + ls);
  return 0;
}

// go_rows: rows of the go accumulators of the trial; r0: any row of the trial
// (stop parameters are trial-constant)
inline double stop_success_lp_raw(const SSSpec& s, const ParamTable& pt, int r0,
                                  const int* go_rows, int n_go, double ssd, double upper) {
  StopIntegrand w{&s, &pt, ssd, par(pt, r0, s.col_muS), par(pt, r0, s.col_sigmaS),
                  par(pt, r0, s.col_tauS), par(pt, r0, s.col_lbS), go_rows, n_go};
  const double lo = ss_stop_window_lo(w.lbS, w.muS, w.sigS);
  const double hi = ss_stop_window_hi(upper, w.muS, w.sigS, w.tauS);
  double res = ss_integrate(stop_success_integrand, &w, lo, hi, SSF_ABS_TOL, SSF_REL_TOL,
                            SSF_STOP_MAX_EVAL);
  return (!is_finite(res) || res <= 0.0) ? neg_inf() : std::log(res);
}

// 64-bit mixing hash over the bit patterns of the key doubles
inline std::uint64_t mix64(std::uint64_t h, double x) {
  std::uint64_t b; std::memcpy(&b, &x, sizeof b);
  h ^= b + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  h *= 0xff51afd7ed558ccdULL; h ^= h >> 33;
  return h;
}

// memoised stop-success integral; key = (ssd, upper, stop params, go params)
inline double stop_success_lp(const SSSpec& s, const ParamTable& pt, int trial, double upper) {
  const int r0 = s.base_row[trial];
  int go_rows[64]; int n_go = 0;
  for (int k = 0; k < s.n_acc && n_go < 64; ++k)
    if (s.is_go[r0 + k]) go_rows[n_go++] = r0 + k;
  // build the key at the end of memo_keys
  std::vector<double>& keys = s.memo_keys;
  const std::size_t k0 = keys.size();
  keys.push_back(s.ssd[trial]); keys.push_back(upper);
  keys.push_back(par(pt, r0, s.col_muS)); keys.push_back(par(pt, r0, s.col_sigmaS));
  keys.push_back(par(pt, r0, s.col_tauS)); keys.push_back(par(pt, r0, s.col_lbS));
  for (int i = 0; i < n_go; ++i) {
    const int r = go_rows[i];
    if (s.is_exg) {
      keys.push_back(par(pt, r, s.col_mu)); keys.push_back(par(pt, r, s.col_sigma));
      keys.push_back(par(pt, r, s.col_tau)); keys.push_back(par(pt, r, s.col_lb));
    } else {
      keys.push_back(par(pt, r, s.col_v)); keys.push_back(par(pt, r, s.col_B));
      keys.push_back(par(pt, r, s.col_A)); keys.push_back(par(pt, r, s.col_t0));
      keys.push_back(par(pt, r, s.col_s));
    }
  }
  const std::size_t klen = keys.size() - k0;
  std::uint64_t h = 0x243f6a8885a308d3ULL;
  for (std::size_t i = k0; i < keys.size(); ++i) h = mix64(h, keys[i]);
  auto it = s.memo.find(h);
  if (it != s.memo.end() && it->second.len == klen) {
    bool same = true;
    for (std::size_t i = 0; i < klen; ++i)
      if (keys[it->second.off + i] != keys[k0 + i]) { same = false; break; }
    if (same) { keys.resize(k0); return it->second.val; }
  }
  const double lp = stop_success_lp_raw(s, pt, r0, go_rows, n_go, s.ssd[trial], upper);
  if (it == s.memo.end()) s.memo.emplace(h, SSSpec::MemoEntry{k0, klen, lp});
  else keys.resize(k0);   // hash collision with a different key: do not memoise
  return lp;
}

// ---------------------------------------------------------------------------
// Trial response density at time t (go path; ST accumulators handled too),
// used for the lower-censoring response mass. Mirrors
// ss_trial_log_response_density in ss_likelihood.h.
// ---------------------------------------------------------------------------
inline double trial_log_response_density(const SSSpec& s, const ParamTable& pt, int trial,
                                         double t, int resp_local) {
  const int r0 = s.base_row[trial];
  const bool stop = s.stop_presented[trial];
  const double tf = par(pt, r0, s.col_tf), gf = par(pt, r0, s.col_gf);
  const double ssd = s.ssd[trial];
  double log_go_any = neg_inf();
  for (int w = 0; w < s.n_acc; ++w) {
    if (!s.is_go[r0 + w]) continue;
    if (resp_local >= 0 && w != resp_local) continue;
    double go_lprob = go_lpdf(s, pt, r0 + w, t);
    if (!is_finite(go_lprob)) continue;
    for (int k = 0; k < s.n_acc; ++k)
      if (s.is_go[r0 + k] && k != w) go_lprob += go_lsurv(s, pt, r0 + k, t);
    if (!is_finite(go_lprob)) continue;
    double term = log1m(gf) + go_lprob;
    if (stop && t > ssd) {
      const double te = t - ssd;
      double lss = stop_lsurv(s, pt, r0, te);
      double st_loss = 0.0;
      for (int k = 0; k < s.n_acc; ++k)
        if (s.is_st[r0 + k]) st_loss += go_lsurv(s, pt, r0 + k, te);
      term = log1m(gf) + log_mix(tf, go_lprob, go_lprob + lss + st_loss);
    }
    log_go_any = log_sum_exp(log_go_any, term);
  }
  if (!stop || s.n_accST[trial] == 0) return log_go_any;

  double log_st_any = neg_inf();
  for (int w = 0; w < s.n_acc; ++w) {
    if (!s.is_st[r0 + w]) continue;
    if (resp_local >= 0 && w != resp_local) continue;
    if (t <= ssd) continue;
    const double te = t - ssd;
    double st_base = go_lpdf(s, pt, r0 + w, te);
    if (!is_finite(st_base)) continue;
    for (int k = 0; k < s.n_acc; ++k)
      if (s.is_st[r0 + k] && k != w) st_base += go_lsurv(s, pt, r0 + k, te);
    double go_loss = 0.0;
    for (int k = 0; k < s.n_acc; ++k)
      if (s.is_go[r0 + k]) go_loss += go_lsurv(s, pt, r0 + k, t);
    double lps = stop_success_lp(s, pt, trial, te);
    const double term_gf   = std::log(gf) + st_base;
    const double term_win  = log1m(gf) + lps + st_base;
    const double term_lose = log1m(gf) + log1m_exp(lps) + st_base + go_loss;
    log_st_any = log_sum_exp(log_st_any,
                             log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_win, term_lose)));
  }
  return log_sum_exp(log_go_any, log_st_any);
}

struct LcIntegrand { const SSSpec* s; const ParamTable* pt; int trial; int resp_local; };

inline int lc_integrand(unsigned, const double* x, void* p, unsigned, double* out) {
  const LcIntegrand* d = static_cast<const LcIntegrand*>(p);
  double lf = trial_log_response_density(*d->s, *d->pt, d->trial, x[0], d->resp_local);
  out[0] = is_finite(lf) ? std::exp(lf) : 0.0;
  return 0;
}

inline double lower_censor_lp(const SSSpec& s, const ParamTable& pt, int trial, double min_ll) {
  const double lo = s.lt[trial], hi = s.lc[trial];
  if (!(hi > lo)) return min_ll;
  LcIntegrand d{&s, &pt, trial, s.resp_local[trial]};
  double res = ss_integrate(lc_integrand, &d, lo, hi, SSF_ABS_TOL, SSF_REL_TOL, SSF_LC_MAX_EVAL);
  return (!is_finite(res) || res <= 0.0) ? min_ll : std::log(res);
}

// log P(no response at all) on a trial without a deadline
inline double withheld_lp(const SSSpec& s, const ParamTable& pt, int trial) {
  const int r0 = s.base_row[trial];
  const double gf = par(pt, r0, s.col_gf);
  if (!s.stop_presented[trial]) return std::log(gf);                 // go failure
  const double tf = par(pt, r0, s.col_tf);
  if (s.n_accST[trial] == 0) {
    double lps = stop_success_lp(s, pt, trial, pos_inf());
    return log_sum_exp(std::log(gf), log1m(gf) + log1m(tf) + lps);   // stop success
  }
  return std::log(gf) + std::log(tf);                                // ST intrinsic NR
}

// log P(no response by the deadline UC)
inline double upper_deadline_lp(const SSSpec& s, const ParamTable& pt, int trial) {
  const double uc = s.uc[trial];
  if (!is_finite(uc)) return withheld_lp(s, pt, trial);
  const int r0 = s.base_row[trial];
  const double gf = par(pt, r0, s.col_gf), tf = par(pt, r0, s.col_tf);
  double logS_go = 0.0;
  for (int k = 0; k < s.n_acc; ++k)
    if (s.is_go[r0 + k]) logS_go += go_lsurv(s, pt, r0 + k, uc);
  if (!s.stop_presented[trial])
    return log_sum_exp(std::log(gf), log1m(gf) + logS_go);
  const double uc_eff = s.uc_eff[trial];
  double lps = neg_inf(), logS_stop = 0.0;
  if (s.stop_can_act[trial]) {
    logS_stop = stop_lsurv(s, pt, r0, uc_eff);
    lps = stop_success_lp(s, pt, trial, uc_eff);
  }
  const double log_core = log_sum_exp(lps, logS_go + logS_stop);
  if (s.n_accST[trial] == 0)
    return log_sum_exp(std::log(gf), log1m(gf) + log_mix(tf, logS_go, log_core));
  double logS_st = 0.0;
  for (int k = 0; k < s.n_acc; ++k)
    if (s.is_st[r0 + k]) logS_st += go_lsurv(s, pt, r0 + k, uc_eff);
  const double log_trig = logS_st + log_sum_exp(std::log(gf), log1m(gf) + log_core);
  const double log_tfb  = log_sum_exp(std::log(gf), log1m(gf) + logS_go);
  return log_mix(tf, log_tfb, log_trig);
}

// observed response: go accumulator won (go trial or signal-respond)
inline double observed_go_lp(const SSSpec& s, const ParamTable& pt, int trial) {
  const int r0 = s.base_row[trial], w = s.winner_local[trial];
  const double t = s.rt[trial];
  const double gf = par(pt, r0, s.col_gf);
  double go_lprob = go_lpdf(s, pt, r0 + w, t);
  if (!is_finite(go_lprob)) go_lprob = neg_inf();
  for (int k = 0; k < s.n_acc; ++k)
    if (s.is_go[r0 + k] && k != w) go_lprob += go_lsurv(s, pt, r0 + k, t);
  if (!s.stop_presented[trial]) return log1m(gf) + go_lprob;
  const double tf = par(pt, r0, s.col_tf), te = s.rt_eff[trial];
  double lss = stop_lsurv(s, pt, r0, te);
  double st_loss = 0.0;
  for (int k = 0; k < s.n_acc; ++k)
    if (s.is_st[r0 + k] && k != w) st_loss += go_lsurv(s, pt, r0 + k, te);
  return log1m(gf) + log_mix(tf, go_lprob, go_lprob + lss + st_loss);
}

// observed response: stop-triggered accumulator won
inline double observed_st_lp(const SSSpec& s, const ParamTable& pt, int trial) {
  const int r0 = s.base_row[trial], w = s.winner_local[trial];
  const double t = s.rt[trial], te = s.rt_eff[trial];
  const double gf = par(pt, r0, s.col_gf), tf = par(pt, r0, s.col_tf);
  double st_base = go_lpdf(s, pt, r0 + w, te);
  if (!is_finite(st_base)) st_base = neg_inf();
  for (int k = 0; k < s.n_acc; ++k)
    if (s.is_st[r0 + k] && k != w) st_base += go_lsurv(s, pt, r0 + k, te);
  double go_loss = 0.0;
  for (int k = 0; k < s.n_acc; ++k)
    if (s.is_go[r0 + k]) go_loss += go_lsurv(s, pt, r0 + k, t);
  double lps = stop_success_lp(s, pt, trial, te);
  const double term_gf   = std::log(gf) + st_base;
  const double term_win  = log1m(gf) + lps + st_base;
  const double term_lose = log1m(gf) + log1m_exp(lps) + st_base + go_loss;
  return log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_win, term_lose));
}

} // namespace ssf

// ---------------------------------------------------------------------------
// Per-particle likelihood: fills ll_buf[trial] (compressed trials). Bounds,
// clamping, expansion and summation are done by the caller (apply_bounds +
// expand_clamp_sum), as for the race models. Thread-safe given a per-thread
// SSSpec (rebind) and ParamTable.
// ---------------------------------------------------------------------------
inline void ss_log_likelihood(const SSSpec& s, double min_ll, double* ll_buf) {
  const ParamTable& pt = *s.pt;
  s.memo.clear();
  s.memo_keys.clear();
  for (int t = 0; t < s.n_trials; ++t) {
    double v;
    switch (s.kind[t]) {
      case SSKind::ObsGo:
      case SSKind::ObsStopFailGo: v = ssf::observed_go_lp(s, pt, t); break;
      case SSKind::ObsST:         v = ssf::observed_st_lp(s, pt, t); break;
      case SSKind::Withheld:      v = ssf::withheld_lp(s, pt, t); break;
      case SSKind::UpperCens:     v = ssf::upper_deadline_lp(s, pt, t); break;
      case SSKind::LowerCens:     v = ssf::lower_censor_lp(s, pt, t, min_ll); break;
      case SSKind::BothCens:      v = log_sum_exp(ssf::lower_censor_lp(s, pt, t, min_ll),
                                                  ssf::upper_deadline_lp(s, pt, t)); break;
      default:                    v = min_ll; break;
    }
    ll_buf[t] = v;
  }
}

// ---------------------------------------------------------------------------
// Data-only setup (R boundary; the only Rcpp use in this file)
// ---------------------------------------------------------------------------
inline Rcpp::NumericVector ss_col_or_default(const Rcpp::DataFrame& df, const char* name,
                                             int n, double def) {
  if (df.containsElementNamed(name)) {
    Rcpp::NumericVector col = df[name];
    Rcpp::NumericVector out = Rcpp::clone(col);
    for (int i = 0; i < out.size(); ++i) if (Rcpp::NumericVector::is_na(out[i])) out[i] = def;
    return out;
  }
  return Rcpp::NumericVector(n, def);
}

inline SSSpec make_ss_spec(const std::string& type, const Rcpp::DataFrame& data,
                           const ParamTable& pt) {
  using namespace Rcpp;
  SSSpec s;
  s.is_exg = (type == "SSEXG");
  if (s.is_exg) {
    s.col_mu = pt.base_index_for("mu"); s.col_sigma = pt.base_index_for("sigma");
    s.col_tau = pt.base_index_for("tau"); s.col_lb = pt.base_index_for("exg_lb");
  } else {
    s.col_v = pt.base_index_for("v"); s.col_B = pt.base_index_for("B");
    s.col_A = pt.base_index_for("A"); s.col_t0 = pt.base_index_for("t0");
    s.col_s = pt.base_index_for("s");
  }
  s.col_muS = pt.base_index_for("muS"); s.col_sigmaS = pt.base_index_for("sigmaS");
  s.col_tauS = pt.base_index_for("tauS"); s.col_lbS = pt.base_index_for("exgS_lb");
  s.col_tf = pt.base_index_for("tf"); s.col_gf = pt.base_index_for("gf");

  NumericVector RT = data["rt"];
  IntegerVector R  = data["R"];
  NumericVector SSD = data["SSD"];
  IntegerVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  const int n_rows = RT.size();
  IntegerVector lI = data.containsElementNamed("lI") ? as<IntegerVector>(data["lI"])
                                                     : IntegerVector(n_rows, 2);
  const bool has_miss = data.containsElementNamed("missingness");
  IntegerVector miss = has_miss ? as<IntegerVector>(data["missingness"])
                                : IntegerVector(n_rows, NA_INTEGER);
  NumericVector LT = ss_col_or_default(data, "LT", n_rows, 0.0);
  NumericVector LC = ss_col_or_default(data, "LC", n_rows, 0.0);
  NumericVector UC = ss_col_or_default(data, "UC", n_rows, pos_inf());

  s.n_acc = unique(lR).size();
  s.n_rows = n_rows;
  s.n_trials = n_rows / s.n_acc;
  const int n = s.n_trials, na = s.n_acc;
  s.kind.assign(n, SSKind::Invalid);
  s.base_row.assign(n, 0); s.n_accG.assign(n, 0); s.n_accST.assign(n, 0);
  s.winner_local.assign(n, -1); s.resp_local.assign(n, -1);
  s.rt.assign(n, 0.0); s.ssd.assign(n, 0.0); s.rt_eff.assign(n, 0.0);
  s.lt.assign(n, 0.0); s.lc.assign(n, 0.0); s.uc.assign(n, 0.0); s.uc_eff.assign(n, 0.0);
  s.stop_presented.assign(n, 0); s.stop_can_act.assign(n, 0);
  s.is_go.assign(n_rows, 1); s.is_st.assign(n_rows, 0);

  for (int t = 0; t < n; ++t) {
    const int r0 = t * na;
    s.base_row[t] = r0;
    // go / stop-triggered accumulators (lI: 2 = go, 1 = ST; as in the R reference
    // the go code is the maximum code present in the trial)
    int go_code = lI[r0];
    for (int k = 1; k < na; ++k) if (lI[r0 + k] > go_code) go_code = lI[r0 + k];
    for (int k = 0; k < na; ++k) {
      const bool go = (lI[r0 + k] == go_code);
      s.is_go[r0 + k] = go; s.is_st[r0 + k] = !go;
      if (go) ++s.n_accG[t]; else ++s.n_accST[t];
    }
    const double rt = RT[r0], ssd = SSD[r0], uc = UC[r0];
    const bool rt_ok = !NumericVector::is_na(rt) && is_finite(rt);
    s.rt[t] = rt_ok ? rt : neg_inf();
    s.ssd[t] = NumericVector::is_na(ssd) ? pos_inf() : ssd;
    s.stop_presented[t] = is_finite(s.ssd[t]);
    s.rt_eff[t] = rt_ok ? std::max(rt - s.ssd[t], 0.0) : 0.0;
    s.lt[t] = LT[r0]; s.lc[t] = LC[r0]; s.uc[t] = uc;
    double uce = uc - s.ssd[t];
    if (!is_finite(uce) || uce <= 0.0) uce = 0.0;
    s.uc_eff[t] = uce;
    s.stop_can_act[t] = is_finite(uce) && uce > 0.0;
    // winner / response
    const bool resp_obs = (R[r0] != NA_INTEGER);
    for (int k = 0; k < na; ++k) {
      if (winner[r0 + k]) s.winner_local[t] = k;
      if (resp_obs && lR[r0 + k] == R[r0]) s.resp_local[t] = k;
    }
    // classification
    int m = has_miss ? miss[r0] : NA_INTEGER;
    if (m == NA_INTEGER && !rt_ok) m = -1;
    if (m == 1)       s.kind[t] = SSKind::LowerCens;
    else if (m == 2)  s.kind[t] = is_finite(uc) ? SSKind::UpperCens : SSKind::Withheld;
    else if (m == 3)  s.kind[t] = SSKind::BothCens;
    else if (m == -1) s.kind[t] = SSKind::Withheld;
    else if (!rt_ok || s.winner_local[t] < 0) s.kind[t] = SSKind::Invalid;
    else if (!s.stop_presented[t]) s.kind[t] = SSKind::ObsGo;
    else s.kind[t] = s.is_go[r0 + s.winner_local[t]] ? SSKind::ObsStopFailGo : SSKind::ObsST;
  }
  s.rebind(pt);
  return s;
}

#endif // ss_fast_h
