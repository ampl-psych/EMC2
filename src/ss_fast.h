#ifndef ss_fast_h
#define ss_fast_h

// ---------------------------------------------------------------------------
// Stop-signal likelihood (all go/stop families) in the race-model style:
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

// Integration accuracy for exGaussian stop-success and lower-censoring integrals
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
  SSStopFamily stop_family = SSStopFamily::ExGaussian;
  // parameter columns in the ParamTable (resolved once by name)
  int col_mu = -1, col_sigma = -1, col_tau = -1, col_lb = -1;   // SSEXG go
  int col_v = -1, col_B = -1, col_A = -1, col_t0 = -1, col_s = -1; // SSRDEX go
  int col_stop1 = -1, col_stop2 = -1, col_tauS = -1, col_lbS = -1;
  int col_tf = -1, col_gf = -1;
  int col_normal_ratio = -1, col_shifted_exp_sigma = -1;

  int n_trials = 0, n_rows = 0, n_acc = 0;

  // per trial (data only)
  std::vector<SSKind>  kind;
  std::vector<int>     base_row, n_accG, n_accST, winner_local, resp_local;
  std::vector<int>     go_rows, go_offsets;
  std::vector<double>  rt, ssd, lt, lc, uc, uc_eff;
  std::vector<std::uint8_t> stop_presented, stop_can_act;
  // per row (data only)
  std::vector<std::uint8_t> is_go, is_st;

  // per-particle memo of the stop-success integral (see header comment):
  // hash -> entry; the full key is kept in memo_keys for verification
  struct MemoEntry { std::size_t off; std::size_t len; double val; };
  mutable std::unordered_map<std::uint64_t, MemoEntry> memo;
  mutable std::vector<double> memo_keys;

  // Prepared per particle: go distributions per row, stop distributions per trial.
  mutable std::vector<ExGaussian> go_exg, stop_exg;
  mutable std::vector<SSStop> stop_other;

  const ParamTable* pt = nullptr;
  void rebind(const ParamTable& p) { pt = &p; }
};

// ---------------------------------------------------------------------------
// Parameter access (all inline; the branch on is_exg is uniform per dadm)
// ---------------------------------------------------------------------------
namespace ssf {

inline double par(const ParamTable& pt, int row, int col) { return pt.base(row, col); }

inline double normal_ratio(const SSSpec& s, const ParamTable& pt, int row) {
  return s.col_normal_ratio < 0 ? 0.0 : par(pt,row,s.col_normal_ratio);
}
inline double shifted_exp_sigma(const SSSpec& s, const ParamTable& pt, int row) {
  return s.col_shifted_exp_sigma < 0 ? 0.0 : par(pt,row,s.col_shifted_exp_sigma);
}

// log density of a go accumulator (row r) finishing at time t
inline double go_lpdf(const SSSpec& s, const ParamTable& pt, int r, double t) {
  if (s.is_exg) {
    double v = s.go_exg[r].pdf(t, true);
    return is_finite(v) ? v : neg_inf();
  }
  const double sc = par(pt, r, s.col_s);
  const double dt = t - par(pt, r, s.col_t0);
  if (!(dt > 0.0)) return neg_inf();
  const double A = par(pt, r, s.col_A) / sc;
  double v = digt(dt, par(pt, r, s.col_B) / sc + 0.5 * A, par(pt, r, s.col_v) / sc, 0.5 * A, true);
  return is_finite(v) ? v : neg_inf();
}

// log survivor of a go accumulator (row r) at time t
inline double go_lsurv(const SSSpec& s, const ParamTable& pt, int r, double t) {
  if (s.is_exg) {
    double v = s.go_exg[r].prob(t, false, true);
    return is_finite(v) ? v : neg_inf();
  }
  const double sc = par(pt, r, s.col_s);
  const double dt = t - par(pt, r, s.col_t0);
  if (!(dt > 0.0)) return 0.0;
  const double A = par(pt, r, s.col_A) / sc;
  double v = pigt(dt, par(pt, r, s.col_B) / sc + 0.5 * A, par(pt, r, s.col_v) / sc, 0.5 * A, false, true);
  return is_finite(v) ? v : neg_inf();
}

// log survivor of the stop process (parameters on row r0) at stop-relative time q
inline double stop_lsurv(const SSSpec& s, int r0, double q) {
  double v = s.stop_family==SSStopFamily::ExGaussian
    ? s.stop_exg[r0 / s.n_acc].prob(q, false, true)
    : s.stop_other[r0 / s.n_acc].log_survivor(q);
  return is_finite(v) ? v : neg_inf();
}

// ---------------------------------------------------------------------------
// Stop-success integral: P(stop finishes before every go accumulator | SSD),
// over stop-relative finish times in [lb, upper] clipped to the finite window.
// ---------------------------------------------------------------------------
struct StopIntegrand {
  const SSSpec* s; const ParamTable* pt;
  double ssd;
  const int* go_rows; int n_go;
  const ExGaussian* exg;
  const SSStop* other;
};

// f_stop(t) * prod_i S_go,i(t + SSD), including the Jacobian when time is transformed.
inline int stop_success_integrand(unsigned, const double* x, void* p, unsigned, double* out) {
  const StopIntegrand* w = static_cast<const StopIntegrand*>(p);
  const double xx = w->other ? w->other->time(x[0]) : x[0];
  double lf = w->other ? w->other->log_weight(x[0]) : w->exg->pdf(xx, true);
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
  const bool exg=s.stop_family==SSStopFamily::ExGaussian;
  const int trial=r0/s.n_acc;
  StopIntegrand w{&s, &pt, ssd, go_rows, n_go,
                  exg ? &s.stop_exg[trial] : nullptr,
                  exg ? nullptr : &s.stop_other[trial]};
  double lo,hi;
  if (exg) {
    const auto& stop=*w.exg;
    lo=ss_stop_window_lo(stop.lower,stop.mu,stop.sigma);
    hi=ss_stop_window_hi(upper,stop.mu,stop.sigma,stop.tau);
  } else w.other->bounds(upper,lo,hi);
  if (!(hi > lo)) return neg_inf();
  // The stop window can extend well beyond the go finishing times. Clip hi
  // where the go survival product falls below exp(-30) for exGaussian stops
  // or exp(-45) for the other families. Bisection uses its monotone decrease.
  auto log_S_go = [&](double x) {
    const double t=w.other ? w.other->time(x) : x;
    double ls = 0.0;
    for (int i = 0; i < n_go; ++i) {
      double v = go_lsurv(s, pt, go_rows[i], t + ssd);
      if (!is_finite(v)) return neg_inf();
      ls += v;
    }
    return ls;
  };
  const double LS_MIN = exg ? -30.0 : -45.0;
  const double ls_lo=log_S_go(lo);
  if (ls_lo < LS_MIN) return neg_inf();
  auto locate = [&](double target) {
    double a=lo, b=hi;
    for (int it=0; it<(exg ? 30 : 48); ++it) {
      const double m=0.5*(a+b);
      if (m==a || m==b) break;
      if (log_S_go(m)<target) b=m; else a=m;
    }
    return b;
  };
  if (log_S_go(hi)<LS_MIN) hi=locate(LS_MIN);
  // A broad stop distribution can compress the go transition into a tiny
  // part of this interval. Split where go survival first changes appreciably,
  // so quadrature cannot miss it. The preceding part is flat to 1e-12 relative.
  double split=lo;
  if (!exg && log_S_go(hi)<ls_lo-1e-12) split=locate(ls_lo-1e-12);
  const double abs_tol=exg ? SSF_ABS_TOL : 1e-18;
  const double rel_tol=exg ? SSF_REL_TOL : 1e-8;
  double res=ss_integrate(stop_success_integrand,&w,lo,split,abs_tol,rel_tol,SSF_STOP_MAX_EVAL)+
    ss_integrate(stop_success_integrand,&w,split,hi,abs_tol,rel_tol,SSF_STOP_MAX_EVAL);
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
  const int* go_rows=s.go_rows.data()+s.go_offsets[trial];
  const int n_go=s.n_accG[trial];
  // build the key at the end of memo_keys
  std::vector<double>& keys = s.memo_keys;
  const std::size_t k0 = keys.size();
  keys.push_back(s.ssd[trial]); keys.push_back(upper);
  keys.push_back(par(pt, r0, s.col_stop1)); keys.push_back(par(pt, r0, s.col_stop2));
  if (s.col_tauS>=0) keys.push_back(par(pt,r0,s.col_tauS));
  if (s.col_lbS>=0) keys.push_back(par(pt,r0,s.col_lbS));
  keys.push_back(normal_ratio(s,pt,r0)); keys.push_back(shifted_exp_sigma(s,pt,r0));
  for (int i = 0; i < n_go; ++i) {
    const int r = go_rows[i];
    if (s.is_exg) {
      keys.push_back(par(pt, r, s.col_mu)); keys.push_back(par(pt, r, s.col_sigma));
      keys.push_back(par(pt, r, s.col_tau)); keys.push_back(par(pt, r, s.col_lb));
      keys.push_back(normal_ratio(s,pt,r)); keys.push_back(shifted_exp_sigma(s,pt,r));
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
// Density for a specified go or stop-triggered response at time t, shared by
// observed trials and the lower-censoring integral. ST densities use t - SSD.
// ---------------------------------------------------------------------------
inline double response_lp(const SSSpec& s, const ParamTable& pt, int trial,
                          double t, int w) {
  const int r0=s.base_row[trial];
  if (s.is_go[r0+w]) {
    const double gf = par(pt, r0, s.col_gf);
    double go_lprob = go_lpdf(s, pt, r0 + w, t);
    if (!is_finite(go_lprob)) return neg_inf();
    for (int k = 0; k < s.n_acc; ++k)
      if (s.is_go[r0 + k] && k != w) go_lprob += go_lsurv(s, pt, r0 + k, t);
    if (!s.stop_presented[trial] || t<=s.ssd[trial]) return log1m(gf) + go_lprob;
    const double tf = par(pt, r0, s.col_tf), te = t-s.ssd[trial];
    double lss = stop_lsurv(s, r0, te);
    double st_loss = 0.0;
    for (int k = 0; k < s.n_acc; ++k)
      if (s.is_st[r0 + k] && k != w) st_loss += go_lsurv(s, pt, r0 + k, te);
    // Trigger failure, or a triggered stop and ST racers still running at t.
    return log1m(gf) + log_mix(tf, go_lprob, go_lprob + lss + st_loss);
  }
  if (!s.stop_presented[trial] || t<=s.ssd[trial]) return neg_inf();
  const double te=t-s.ssd[trial];
  const double gf = par(pt, r0, s.col_gf), tf = par(pt, r0, s.col_tf);
  double st_base = go_lpdf(s, pt, r0 + w, te);
  if (!is_finite(st_base)) return neg_inf();
  for (int k = 0; k < s.n_acc; ++k)
    if (s.is_st[r0 + k] && k != w) st_base += go_lsurv(s, pt, r0 + k, te);
  double go_loss = 0.0;
  for (int k = 0; k < s.n_acc; ++k)
    if (s.is_go[r0 + k]) go_loss += go_lsurv(s, pt, r0 + k, t);
  double lps = stop_success_lp(s, pt, trial, te);
  // An ST response can follow go failure, a stop win, or go and stop still running.
  const double term_gf   = std::log(gf) + st_base;
  const double term_win  = log1m(gf) + lps + st_base;
  const double term_lose = log1m(gf) + stop_lsurv(s, r0, te) + st_base + go_loss;
  return log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_win, term_lose));
}

// Sum over all possible winners when the response identity is missing.
inline double trial_log_response_density(const SSSpec& s, const ParamTable& pt, int trial,
                                         double t, int resp_local) {
  if (resp_local>=0) return response_lp(s,pt,trial,t,resp_local);
  double out=neg_inf();
  for (int w=0; w<s.n_acc; ++w)
    out=log_sum_exp(out,response_lp(s,pt,trial,t,w));
  return out;
}

struct LcIntegrand { const SSSpec* s; const ParamTable* pt; int trial; int resp_local; };

inline int lc_integrand(unsigned, const double* x, void* p, unsigned, double* out) {
  const LcIntegrand* d = static_cast<const LcIntegrand*>(p);
  double lf = trial_log_response_density(*d->s, *d->pt, d->trial, x[0], d->resp_local);
  out[0] = is_finite(lf) ? std::exp(lf) : 0.0;
  return 0;
}

inline double lower_censor_lp(const SSSpec& s, const ParamTable& pt, int trial) {
  const double lo = s.lt[trial], hi = s.lc[trial];
  if (!(hi > lo)) return neg_inf();
  LcIntegrand d{&s, &pt, trial, s.resp_local[trial]};
  double res = ss_integrate(lc_integrand, &d, lo, hi, SSF_ABS_TOL, SSF_REL_TOL, SSF_LC_MAX_EVAL);
  return (!is_finite(res) || res <= 0.0) ? neg_inf() : std::log(res);
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
    logS_stop = stop_lsurv(s, r0, uc_eff);
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

} // namespace ssf

// Reset the integral cache and prepare distributions for this particle.
inline void ss_prepare_distributions(const SSSpec& s) {
  const ParamTable& pt = *s.pt;
  s.memo.clear();
  s.memo_keys.clear();
  // Parameter-dependent preparation lives in the worker's SSSpec and is
  // refreshed on every likelihood evaluation, including after rebind().
  s.go_exg.clear();
  if (s.is_exg) {
    s.go_exg.reserve(s.n_rows);
    for (int r=0; r<s.n_rows; ++r)
      s.go_exg.emplace_back(ssf::par(pt,r,s.col_mu), ssf::par(pt,r,s.col_sigma),
                           ssf::par(pt,r,s.col_tau), ssf::par(pt,r,s.col_lb), pos_inf(),
                           ssf::normal_ratio(s,pt,r), ssf::shifted_exp_sigma(s,pt,r));
  }
  s.stop_exg.clear();
  s.stop_other.clear();
  if (s.stop_family==SSStopFamily::ExGaussian) {
    s.stop_exg.reserve(s.n_trials);
    for (int t=0; t<s.n_trials; ++t) {
      const int r=s.base_row[t];
      s.stop_exg.emplace_back(ssf::par(pt,r,s.col_stop1), ssf::par(pt,r,s.col_stop2),
                             ssf::par(pt,r,s.col_tauS), ssf::par(pt,r,s.col_lbS), pos_inf(),
                             ssf::normal_ratio(s,pt,r), ssf::shifted_exp_sigma(s,pt,r));
    }
  } else {
    s.stop_other.reserve(s.n_trials);
    for (int t=0; t<s.n_trials; ++t) {
      const int r=s.base_row[t];
      s.stop_other.emplace_back(s.stop_family,ssf::par(pt,r,s.col_stop1),ssf::par(pt,r,s.col_stop2),
                                s.col_lbS>=0 ? ssf::par(pt,r,s.col_lbS) : 0.0);
    }
  }
}

// ---------------------------------------------------------------------------
// Per-particle likelihood: fills ll_buf[trial] (compressed trials). The caller
// checks bounds first and applies penalties, clamping, expansion and summation
// afterwards, as for the race models. Thread-safe given a per-thread
// SSSpec (rebind) and ParamTable.
// ---------------------------------------------------------------------------
inline void ss_log_likelihood(const SSSpec& s, double* ll_buf,
                              const std::vector<int>& is_ok) {
  ss_prepare_distributions(s);
  const ParamTable& pt=*s.pt;
  for (int t = 0; t < s.n_trials; ++t) {
    if (!is_ok[s.base_row[t]]) continue;
    double v;
    switch (s.kind[t]) {
      case SSKind::ObsGo:
      case SSKind::ObsStopFailGo:
      case SSKind::ObsST:         v = ssf::response_lp(s, pt, t, s.rt[t], s.winner_local[t]); break;
      case SSKind::Withheld:      v = ssf::withheld_lp(s, pt, t); break;
      case SSKind::UpperCens:     v = ssf::upper_deadline_lp(s, pt, t); break;
      case SSKind::LowerCens:     v = ssf::lower_censor_lp(s, pt, t); break;
      case SSKind::BothCens:      v = log_sum_exp(ssf::lower_censor_lp(s, pt, t),
                                                  ssf::upper_deadline_lp(s, pt, t)); break;
      default:                    v = neg_inf(); break;
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

inline void ss_set_parameter_columns(SSSpec& s, const std::string& type, const ParamTable& pt) {
  s.is_exg = type=="SSEXG" || type=="SSLNORM" || type=="SSNORM" || type=="SSWEIBULL";
  if (type=="SSLNORM" || type=="SSRDLNORM") s.stop_family=SSStopFamily::Lognormal;
  else if (type=="SSNORM" || type=="SSRDNORM") s.stop_family=SSStopFamily::Normal;
  else if (type=="SSWEIBULL" || type=="SSRDWEIBULL") s.stop_family=SSStopFamily::Weibull;
  if (s.is_exg) {
    s.col_mu = pt.base_index_for("mu"); s.col_sigma = pt.base_index_for("sigma");
    s.col_tau = pt.base_index_for("tau"); s.col_lb = pt.base_index_for("exg_lb");
  } else {
    s.col_v = pt.base_index_for("v"); s.col_B = pt.base_index_for("B");
    s.col_A = pt.base_index_for("A"); s.col_t0 = pt.base_index_for("t0");
    s.col_s = pt.base_index_for("s");
  }
  switch (s.stop_family) {
    case SSStopFamily::ExGaussian:
      s.col_tauS=pt.base_index_for("tauS");
      s.col_lbS=pt.base_index_for("exgS_lb");
      // Preserve approximation support for the existing EXG-stop models only.
      if (pt.name_to_base_idx.count("exg_normal_ratio"))
        s.col_normal_ratio=pt.base_index_for("exg_normal_ratio");
      if (pt.name_to_base_idx.count("exg_shifted_exp_sigma"))
        s.col_shifted_exp_sigma=pt.base_index_for("exg_shifted_exp_sigma");
      [[fallthrough]];
    case SSStopFamily::Normal:
      s.col_stop1=pt.base_index_for("muS"); s.col_stop2=pt.base_index_for("sigmaS");
      if (s.stop_family==SSStopFamily::Normal) s.col_lbS=pt.base_index_for("normS_lb");
      break;
    case SSStopFamily::Lognormal:
      s.col_stop1=pt.base_index_for("meanlogS"); s.col_stop2=pt.base_index_for("sdlogS");
      break;
    case SSStopFamily::Weibull:
      s.col_stop1=pt.base_index_for("shapeS"); s.col_stop2=pt.base_index_for("scaleS");
      break;
  }
}

inline SSSpec make_ss_spec(const std::string& type, const Rcpp::DataFrame& data,
                           const ParamTable& pt) {
  using namespace Rcpp;
  SSSpec s;
  ss_set_parameter_columns(s,type,pt);
  s.col_tf=pt.base_index_for("tf"); s.col_gf=pt.base_index_for("gf");

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
  s.rt.assign(n, 0.0); s.ssd.assign(n, 0.0);
  s.lt.assign(n, 0.0); s.lc.assign(n, 0.0); s.uc.assign(n, 0.0); s.uc_eff.assign(n, 0.0);
  s.stop_presented.assign(n, 0); s.stop_can_act.assign(n, 0);
  s.go_rows.reserve(n_rows); s.go_offsets.resize(n);
  s.is_go.assign(n_rows, 1); s.is_st.assign(n_rows, 0);

  for (int t = 0; t < n; ++t) {
    const int r0 = t * na;
    s.base_row[t] = r0;
    s.go_offsets[t] = s.go_rows.size();
    // go / stop-triggered accumulators (lI: 2 = go, 1 = ST; as in the R reference
    // the go code is the maximum code present in the trial)
    int go_code = lI[r0];
    for (int k = 1; k < na; ++k) if (lI[r0 + k] > go_code) go_code = lI[r0 + k];
    for (int k = 0; k < na; ++k) {
      const bool go = (lI[r0 + k] == go_code);
      s.is_go[r0 + k] = go; s.is_st[r0 + k] = !go;
      if (go) { ++s.n_accG[t]; s.go_rows.push_back(r0+k); } else ++s.n_accST[t];
    }
    const double rt = RT[r0], ssd = SSD[r0], uc = UC[r0];
    const bool rt_ok = !NumericVector::is_na(rt) && is_finite(rt);
    s.rt[t] = rt_ok ? rt : neg_inf();
    s.ssd[t] = NumericVector::is_na(ssd) ? pos_inf() : ssd;
    s.stop_presented[t] = is_finite(s.ssd[t]);
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
