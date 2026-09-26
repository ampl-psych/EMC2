// Neural-likelihood models in the compiled likelihood pipeline (calc_ll /
// calc_ll_multithreaded, type "NN"). Header-only; include from particle_ll.cpp.
//
// The model list's registration (R/nn_register.R, nn_native_args()) gives the
// evaluator, the network's inputs as EMC2 parameter names in the network's
// order with their scales, and an optional parameter subtracted from rt
// (`pre`). Per particle, after the ordinary parameter pipeline:
//   natural-scale columns of the ParamTable -> network inputs -> evaluator
//   flow_joint, regression_joint, mlp_joint: log p(rt - pre, R | theta) per trial
//   flow_race:  log pdf for the winning accumulator, log survivor for the
//               losers, summed per trial
// Bounds, expansion and clamping are the pipeline's, as for other models.
// This reproduces the R path (dfun/pfun via log_likelihood_ddm/race), which
// stays the reference; the one numerical difference is that the survivor is
// taken as log(1 - F) computed directly (pnorm upper tail) instead of
// log(1 - pfun()).

#ifndef EMC2_MODEL_NN_H
#define EMC2_MODEL_NN_H

#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include "ParamTable.h"
#include "nle_native.h"

struct NnCheck {        // a run-time refusal (see nn_check_pars() in R)
  int col;
  double value;
  bool fixed;           // true: any value != `value` refused; false: all == `value` refused
  std::string msg;
};

struct NnSpec {
  bool race = false;                       // flow_race; else a joint (two-response) network
  bool mlp_kind = false;                   // regression_joint / mlp_joint (plain MLP)
  const MlpLik* mlp = nullptr;
  const DdmEnsemble* ddm = nullptr;
  const FlowModel* flow = nullptr;
  int nc = 0;                              // network inputs
  std::vector<int> col;                    // ParamTable column per input
  std::vector<int> tf;                     // 0 identity, 1 log, 2 probit, 3 zerobox
  std::vector<double> tfc;                 // zerobox constant per input (NaN otherwise)
  int pre_col = -1;                        // parameter subtracted from rt
  int n_rows = 0, n_trials = 0, n_acc = 1;
  const double* rt = nullptr;
  Rcpp::NumericVector rt_keep;             // keeps rt alive
  std::vector<int> R;                      // flow_joint: response code (1, 2)
  std::vector<unsigned char> winner;       // flow_race
  std::vector<unsigned char> absent;       // flow_race: accumulator not in the race (RACE)
  std::vector<NnCheck> checks;
};

// Per-thread scratch, reused across particles.
struct NnScratch {
  std::vector<int> rows;
  std::vector<double> th, tn, lp, lsf, ll_row;
  std::vector<int> Rr;
  std::vector<unsigned char> want;
};

inline bool nn_has_col(const Rcpp::DataFrame& data, const char* nm) {
  Rcpp::CharacterVector nms = data.names();
  for (int i = 0; i < nms.size(); ++i)
    if (std::string(nms[i]) == nm) return true;
  return false;
}

inline NnSpec make_nn_spec(const Rcpp::List& nn, const Rcpp::DataFrame& data,
                           const ParamTable& pt, int n_lR) {
  using namespace Rcpp;
  NnSpec s;
  const std::string label = as<std::string>(nn["label"]);
  const std::string what = "Neural likelihood '" + label + "': ";
  const std::string kind = as<std::string>(nn["kind"]);
  if (kind == "flow_joint") s.race = false;
  else if (kind == "flow_race") s.race = true;
  else if (kind == "regression_joint" || kind == "mlp_joint") s.mlp_kind = true;
  else stop(what + "no compiled likelihood for kind '" + kind + "'");

  s.n_rows = data.nrow();
  s.n_acc = n_lR;
  s.n_trials = s.n_rows / n_lR;
  if (!s.race && n_lR != 1) stop(what + "a joint (DDM-type) network takes one data row per trial");

  // Censoring and truncation need the network's CDF inside the pipeline's
  // censor/truncation machinery, which is not wired for neural likelihoods.
  if (nn_has_col(data, "missingness")) {
    IntegerVector miss = data["missingness"];
    for (int i = 0; i < miss.size(); ++i)
      if (!IntegerVector::is_na(miss[i]))
        stop(what + "censored data (a non-missing 'missingness' code) is not supported");
  }
  if (nn_has_col(data, "LT")) {
    NumericVector lt = data["LT"];
    for (int i = 0; i < lt.size(); ++i)
      if (lt[i] > 0.0 && std::isfinite(lt[i])) stop(what + "truncated data (LT > 0) is not supported");
  }
  if (nn_has_col(data, "UT")) {
    NumericVector ut = data["UT"];
    for (int i = 0; i < ut.size(); ++i)
      if (std::isfinite(ut[i])) stop(what + "truncated data (finite UT) is not supported");
  }

  // evaluator and inputs
  SEXP ptr = nn["ptr"];
  if (s.race) { s.flow = nle_flow_handle(ptr); s.nc = nle_flow_n_ctx(s.flow); }
  else if (s.mlp_kind) { s.mlp = nle_mlp_handle(ptr); s.nc = nle_mlp_n_ctx(s.mlp); }
  else        { s.ddm = nle_ddm_handle(ptr);   s.nc = nle_ddm_n_ctx(s.ddm); }
  CharacterVector pars = nn["pars"];
  IntegerVector codes = nn["transform_codes"];
  NumericVector zc = nn.containsElementNamed("zerobox_c") ? NumericVector(nn["zerobox_c"])
                                                           : NumericVector(codes.size(), NA_REAL);
  if (zc.size() != codes.size()) stop(what + "zerobox_c needs one entry per input");
  if (pars.size() != s.nc || codes.size() != s.nc)
    stop(what + "the registration names " + std::to_string((int)pars.size()) +
         " inputs; the evaluator takes " + std::to_string(s.nc));
  for (int j = 0; j < s.nc; ++j) {
    const std::string p = as<std::string>(pars[j]);
    auto it = pt.name_to_base_idx.find(p);
    if (it == pt.name_to_base_idx.end()) stop(what + "the model does not supply input '" + p + "'");
    if (codes[j] < 0 || codes[j] > 3) stop(what + "unknown transform code for input '" + p + "'");
    if (codes[j] == 3 && !(zc[j] > 0.0 && std::isfinite(zc[j])))
      stop(what + "zerobox input '" + p + "' has no positive constant");
    s.col.push_back(it->second);
    s.tf.push_back(codes[j]);
    s.tfc.push_back(codes[j] == 3 ? zc[j] : NAN);
  }
  const std::string pre = as<std::string>(nn["pre"]);
  if (!pre.empty()) {
    auto it = pt.name_to_base_idx.find(pre);
    if (it == pt.name_to_base_idx.end()) stop(what + "the model does not supply pre parameter '" + pre + "'");
    s.pre_col = it->second;
  }

  // data
  s.rt_keep = data["rt"];
  s.rt = s.rt_keep.begin();
  if (!s.race) {
    IntegerVector R = data["R"];
    s.R.assign(R.begin(), R.end());
    for (int i = 0; i < s.n_rows; ++i)
      if (s.R[i] != 1 && s.R[i] != 2) stop(what + "a joint network needs a response (1 or 2) on every row");
  } else {
    LogicalVector w = data["winner"];
    s.winner.resize(s.n_rows);
    for (int i = 0; i < s.n_rows; ++i) s.winner[i] = (w[i] == TRUE);
    if (nn_has_col(data, "RACE")) {          // as log_likelihood_race: lR > RACE absent
      IntegerVector lR = data["lR"];
      RObject race_obj = data["RACE"];
      std::vector<double> race(s.n_rows);
      if (Rf_isFactor(race_obj)) {
        IntegerVector rc(race_obj);
        CharacterVector lev = rc.attr("levels");
        for (int i = 0; i < s.n_rows; ++i)
          race[i] = IntegerVector::is_na(rc[i]) ? NA_REAL : std::atof(CHAR(STRING_ELT(lev, rc[i] - 1)));
      } else {
        NumericVector rn(race_obj);
        std::copy(rn.begin(), rn.end(), race.begin());
      }
      s.absent.resize(s.n_rows);
      for (int i = 0; i < s.n_rows; ++i) s.absent[i] = (lR[i] > race[i]);
    }
  }

  // run-time refusals, as nn_check_pars() on the R path
  auto add_checks = [&](const char* vals, const char* msgs, bool fixed) {
    NumericVector v = nn[vals];
    CharacterVector m = nn[msgs];
    if (v.size() == 0) return;
    CharacterVector nms = v.names();
    for (int k = 0; k < v.size(); ++k) {
      auto it = pt.name_to_base_idx.find(as<std::string>(nms[k]));
      if (it == pt.name_to_base_idx.end()) continue;
      s.checks.push_back(NnCheck{it->second, v[k], fixed, as<std::string>(m[k])});
    }
  };
  add_checks("fixed", "fixed_msg", true);
  add_checks("refuse", "refuse_msg", false);
  return s;
}

// Main thread only (stops): the refusals the R path's Ttransform makes, on one
// mapped particle. design() refuses these already; this is the backstop.
inline void nn_check_refusals(const NnSpec& s, const ParamTable& pt) {
  for (const NnCheck& c : s.checks) {
    const double* v = pt.base.colptr(c.col);
    bool refuse = !c.fixed;
    for (int i = 0; i < pt.n_trials; ++i) {
      if (c.fixed && !std::isnan(v[i]) && v[i] != c.value) { refuse = true; break; }
      if (!c.fixed && !(v[i] == c.value)) { refuse = false; break; }
    }
    if (refuse) Rcpp::stop(c.msg);
  }
}

// Trial log-likelihoods of one mapped particle into ll_trial (n_trials),
// before bounds and clamping. Reentrant.
inline void nn_trial_ll(const NnSpec& s, const ParamTable& pt, NnScratch& w,
                        double min_ll, double* ll_trial) {
  const double* pre = s.pre_col >= 0 ? pt.base.colptr(s.pre_col) : nullptr;
  // rows the network evaluates: time > 0 (else zero density / CDF, as the
  // R path's dfun/pfun) and, in a race, accumulators taking part
  w.rows.clear(); w.tn.clear();
  for (int r = 0; r < s.n_rows; ++r) {
    if (!s.absent.empty() && s.absent[r]) continue;
    const double t = pre ? s.rt[r] - pre[r] : s.rt[r];
    if (!(t > 0.0)) continue;
    w.rows.push_back(r); w.tn.push_back(t);
  }
  const int m = (int)w.rows.size();
  w.th.resize((size_t)m * s.nc);
  for (int j = 0; j < s.nc; ++j) {
    const double* c = pt.base.colptr(s.col[j]);
    double* dst = w.th.data() + (size_t)j * m;
    for (int i = 0; i < m; ++i) dst[i] = c[w.rows[i]];
  }
  w.lp.resize(m);

  if (!s.race) {
    w.Rr.resize(m);
    for (int i = 0; i < m; ++i) w.Rr[i] = s.R[w.rows[i]];
    if (m > 0) {
      if (s.mlp_kind) nle_mlp_native(s.mlp, w.th.data(), m, s.tf.data(), s.tfc.data(), w.tn.data(), w.Rr.data(), w.lp.data());
      else nle_ddm_native(s.ddm, w.th.data(), m, s.tf.data(), s.tfc.data(), w.tn.data(), w.Rr.data(), w.lp.data());
    }
    std::fill(ll_trial, ll_trial + s.n_trials, R_NegInf);
    for (int i = 0; i < m; ++i) ll_trial[w.rows[i]] = w.lp[i];
    return;
  }

  w.lsf.resize(m);
  w.want.resize(m);
  for (int i = 0; i < m; ++i) w.want[i] = s.winner[w.rows[i]] ? 1 : 2;
  if (m > 0) nle_flow_native(s.flow, w.th.data(), m, s.tf.data(), s.tfc.data(), w.tn.data(), w.want.data(),
                             w.lp.data(), w.lsf.data());
  // not evaluated: a winner has no density, a loser has not finished
  w.ll_row.resize(s.n_rows);
  for (int r = 0; r < s.n_rows; ++r) w.ll_row[r] = s.winner[r] ? R_NegInf : 0.0;
  for (int i = 0; i < m; ++i) {
    const int r = w.rows[i];
    w.ll_row[r] = s.winner[r] ? w.lp[i] : w.lsf[i];
  }
  for (int t = 0; t < s.n_trials; ++t) {
    double ll = 0.0;
    for (int k = 0; k < s.n_acc; ++k) {
      const double v = w.ll_row[(size_t)t * s.n_acc + k];
      ll += std::isnan(v) ? min_ll : v;       // as log_likelihood_race
    }
    ll_trial[t] = ll;
  }
}

#endif
