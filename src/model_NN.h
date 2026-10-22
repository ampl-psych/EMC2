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
//   hybrid flow_race: a row whose hybrid parameter (e.g. the pulse height amp)
//               is exactly 0 is the closed-form Wald (nle_wald.h) and never
//               reaches the network
//   scaled flow_joint (DDMnle): the card models the unit-scale time u = (rt - t0)(s/a)^2 given
//               kappa = v a / s^2, eta = sv a / s^2, Z, SZ; c = (s/a)^2 and the uniform st0
//               window are exact: f = c g(u) (st0 = 0), [G(u1) - G(u0)] / st0 (st0 > 0),
//               as port/R/ddm_scaled.R's ddm_sc_density() (the R reference)
//   flow_cat:   K categories k = rating_fold(R, RR); the card's convolution
//               parameters (t0, st0) shift and smear the decision time exactly,
//               censored rows (missingness 1 / 2) are scored from LC / UC by the
//               evaluator (flow_cat.cpp; R reference nn_loglik_cat())
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
#include "nle_wald.h"

struct NnCheck {        // a run-time refusal (see nn_check_pars() in R)
  int col;
  double value;
  bool fixed;           // true: any value != `value` refused; false: all == `value` refused
  std::string msg;
};

struct NnSpec {
  bool race = false;                       // flow_race; else a joint (two-response) network
  bool mlp_kind = false;                   // regression_joint / mlp_joint (plain MLP)
  bool cat_kind = false;                   // flow_cat: K categories, uniform t0/st0 convolution, censoring
  bool scaled = false;                     // scaled DDM card (see the header)
  int sc_col[8] = {-1, -1, -1, -1, -1, -1, -1, -1};   // v, a, sv, t0, st0, s, Z, SZ
  const MlpLik* mlp = nullptr;
  const DdmEnsemble* ddm = nullptr;
  const FlowModel* flow = nullptr;
  const CatModel* cat = nullptr;
  int conv_mean_col = -1, conv_range_col = -1;   // flow_cat: t0 and st0 of the convolution
  std::vector<int> k;                      // flow_cat: category 0..K-1 (-1 on censored rows)
  std::vector<int> miss;                   // flow_cat: 0 observed, 1 lower, 2 upper censored
  std::vector<double> LC, UC;              // flow_cat: censoring bounds per row
  int nc = 0;                              // network inputs
  std::vector<int> col;                    // ParamTable column per input
  std::vector<int> tf;                     // 0 identity, 1 log, 2 probit, 3 zerobox
  std::vector<double> tfc;                 // zerobox constant per input (NaN otherwise)
  int pre_col = -1;                        // parameter subtracted from rt
  int hyb_col = -1;                        // hybrid race: rows with this parameter == 0 are Wald
  int wald_col[3] = {-1, -1, -1};          // the Wald's v, b, s
  double wald_lower[3] = {0, 0, 0}, wald_upper[3] = {0, 0, 0};   // their training box (closed)
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
  std::vector<int> rows, wrows;
  std::vector<double> th, tn, wtn, lp, lsf, ll_row;
  std::vector<int> Rr;
  std::vector<double> cdf;
  std::vector<int> sc_kind, sc_row;        // scaled: request kind (0 pdf, 1 G(u1), 2 G(u0)), trial
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
  else if (kind == "flow_cat") s.cat_kind = true;
  else stop(what + "no compiled likelihood for kind '" + kind + "'");

  s.n_rows = data.nrow();
  s.n_acc = n_lR;
  s.n_trials = s.n_rows / n_lR;
  if (!s.race && n_lR != 1) stop(what + "a joint (DDM-type) network takes one data row per trial");

  // Censoring and truncation need the network's CDF inside the pipeline's
  // censor/truncation machinery, which is not wired for neural likelihoods,
  // except the K-category flow, whose evaluator handles lower/upper censoring
  // (codes 1 and 2) itself through its CDF.
  if (nn_has_col(data, "missingness")) {
    IntegerVector miss = data["missingness"];
    for (int i = 0; i < miss.size(); ++i)
      if (!IntegerVector::is_na(miss[i])) {
        if (s.cat_kind && (miss[i] == 1 || miss[i] == 2)) continue;
        stop(what + (s.cat_kind ? "missingness code " + std::to_string((int)miss[i]) + " is not supported " +
                                  "(only 1, lower-censored, and 2, upper-censored)"
                                : "censored data (a non-missing 'missingness' code) is not supported"));
      }
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
  else if (s.cat_kind) { s.cat = nle_cat_handle(ptr); s.nc = nle_cat_n_ctx(s.cat); }
  else        { s.ddm = nle_ddm_handle(ptr);   s.nc = nle_ddm_n_ctx(s.ddm); }
  s.scaled = nn.containsElementNamed("scaled") && as<bool>(nn["scaled"]);
  if (s.scaled) {
    if (s.race || s.mlp_kind || s.cat_kind) stop(what + "only a joint flow can be a scaled DDM card");
    if (s.nc != 4) stop(what + "a scaled DDM card takes 4 inputs (kappa, eta, Z, SZ)");
    CharacterVector sp = nn["sc_pars"];
    if (sp.size() != 8) stop(what + "sc_pars must name v, a, sv, t0, st0, s, Z, SZ");
    for (int j = 0; j < 8; ++j) {
      auto it = pt.name_to_base_idx.find(as<std::string>(sp[j]));
      if (it == pt.name_to_base_idx.end())
        stop(what + "the model does not supply '" + as<std::string>(sp[j]) + "'");
      s.sc_col[j] = it->second;
    }
  }
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
    if (s.scaled) it = pt.name_to_base_idx.end();      // inputs are computed, not columns
    else if (it == pt.name_to_base_idx.end()) stop(what + "the model does not supply input '" + p + "'");
    if (codes[j] < 0 || codes[j] > 3) stop(what + "unknown transform code for input '" + p + "'");
    if (codes[j] == 3 && !(zc[j] > 0.0 && std::isfinite(zc[j])))
      stop(what + "zerobox input '" + p + "' has no positive constant");
    s.col.push_back(s.scaled ? -1 : it->second);
    s.tf.push_back(codes[j]);
    s.tfc.push_back(codes[j] == 3 ? zc[j] : NAN);
  }
  if (s.cat_kind) {                         // the convolution's mean and range
    for (int w = 0; w < 2; ++w) {
      const char* field = w == 0 ? "conv_mean" : "conv_range";
      if (!nn.containsElementNamed(field)) stop(what + "the registration names no " + field);
      const std::string nm = as<std::string>(nn[field]);
      auto it = pt.name_to_base_idx.find(nm);
      if (it == pt.name_to_base_idx.end()) stop(what + "the model does not supply convolution parameter '" + nm + "'");
      (w == 0 ? s.conv_mean_col : s.conv_range_col) = it->second;
    }
  }
  const std::string pre = as<std::string>(nn["pre"]);
  if (!pre.empty() && s.scaled) stop(what + "a scaled DDM card takes t0 itself; pre must be empty");
  if (!pre.empty()) {
    auto it = pt.name_to_base_idx.find(pre);
    if (it == pt.name_to_base_idx.end()) stop(what + "the model does not supply pre parameter '" + pre + "'");
    s.pre_col = it->second;
  }

  const std::string hybrid = nn.containsElementNamed("hybrid") ? as<std::string>(nn["hybrid"]) : "";
  if (!hybrid.empty()) {
    if (!s.race) stop(what + "only a race flow can be a hybrid");
    auto it = pt.name_to_base_idx.find(hybrid);
    if (it == pt.name_to_base_idx.end()) stop(what + "the model does not supply hybrid parameter '" + hybrid + "'");
    s.hyb_col = it->second;
    CharacterVector wp = nn["wald_pars"];
    NumericVector wl = nn["wald_lower"], wu = nn["wald_upper"];
    if (wp.size() != 3 || wl.size() != 3 || wu.size() != 3) stop(what + "a hybrid needs the Wald's v, b and s");
    for (int j = 0; j < 3; ++j) {
      auto iw = pt.name_to_base_idx.find(as<std::string>(wp[j]));
      if (iw == pt.name_to_base_idx.end())
        stop(what + "the model does not supply Wald parameter '" + as<std::string>(wp[j]) + "'");
      s.wald_col[j] = iw->second;
      s.wald_lower[j] = wl[j]; s.wald_upper[j] = wu[j];
    }
  }

  // data
  s.rt_keep = data["rt"];
  s.rt = s.rt_keep.begin();
  if (s.cat_kind) {
    // category k = rating_fold(R, RR, K / 2): the first response level counts
    // its ratings down (R = 1: k = Kr - RR), the second up (R = 2: k = Kr + RR - 1)
    const int K = nle_cat_n_categories(s.cat), Kr = K / 2;
    if (2 * Kr != K) stop(what + "rating_fold needs an even number of categories");
    IntegerVector R = data["R"];
    if (!nn_has_col(data, "RR")) stop(what + "a K-category network needs the rating column RR");
    NumericVector RR = data["RR"];
    s.k.assign(s.n_rows, -1); s.miss.assign(s.n_rows, 0);
    s.LC.assign(s.n_rows, 0.0); s.UC.assign(s.n_rows, R_PosInf);
    if (nn_has_col(data, "missingness")) {
      IntegerVector ms = data["missingness"];
      for (int i = 0; i < s.n_rows; ++i) s.miss[i] = IntegerVector::is_na(ms[i]) ? 0 : ms[i];
    }
    if (nn_has_col(data, "LC")) { NumericVector lc = data["LC"]; for (int i = 0; i < s.n_rows; ++i) s.LC[i] = lc[i]; }
    if (nn_has_col(data, "UC")) { NumericVector uc = data["UC"]; for (int i = 0; i < s.n_rows; ++i) s.UC[i] = uc[i]; }
    for (int i = 0; i < s.n_rows; ++i) {
      if (s.miss[i] != 0) continue;
      if (IntegerVector::is_na(R[i]) || (R[i] != 1 && R[i] != 2) || NumericVector::is_na(RR[i]))
        stop(what + "an observed row needs a response R (1 or 2) and a rating RR");
      const int rr = (int)RR[i];
      if (rr < 1 || rr > Kr) stop(what + "RR must lie in 1.." + std::to_string(Kr));
      s.k[i] = (R[i] == 1) ? (Kr - rr) : (Kr + rr - 1);
    }
  } else if (!s.race) {
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

// Scaled DDM card: one network call for all trials' densities (st0 = 0, or a window too narrow
// to difference) and CDF pairs (st0 > 0). Mirrors ddm_sc_density() in port/R/ddm_scaled.R.
inline void nn_trial_ll_scaled(const NnSpec& s, const ParamTable& pt, NnScratch& w, double* ll_trial) {
  const double *V = pt.base.colptr(s.sc_col[0]), *A = pt.base.colptr(s.sc_col[1]),
               *SV = pt.base.colptr(s.sc_col[2]), *T0 = pt.base.colptr(s.sc_col[3]),
               *ST0 = pt.base.colptr(s.sc_col[4]), *S = pt.base.colptr(s.sc_col[5]),
               *Z = pt.base.colptr(s.sc_col[6]), *SZ = pt.base.colptr(s.sc_col[7]);
  std::fill(ll_trial, ll_trial + s.n_trials, R_NegInf);
  w.rows.clear(); w.tn.clear(); w.Rr.clear(); w.sc_kind.clear(); w.sc_row.clear();
  std::vector<double>& th = w.th;
  th.clear();
  // requests in row order; theta kept row-major here, transposed below
  for (int r = 0; r < s.n_rows; ++r) {
    const double d = s.rt[r] - T0[r];
    if (!(d > 0.0)) continue;                         // before onset: zero density
    const double s2 = S[r] * S[r], c2 = s2 / (A[r] * A[r]);
    const double un[4] = {V[r] * A[r] / s2, SV[r] * A[r] / s2, Z[r], SZ[r] > 0.0 ? SZ[r] : 0.0};
    const double st0 = ST0[r];
    auto push = [&](int kind, double u) {
      th.insert(th.end(), un, un + 4); w.tn.push_back(u); w.Rr.push_back(s.R[r]);
      w.sc_kind.push_back(kind); w.sc_row.push_back(r);
    };
    if (st0 <= 0.0 || st0 < 1e-7 * d) push(0, (d - std::max(st0, 0.0) / 2.0) * c2);
    else {
      push(1, d * c2);
      const double u0 = (d - st0) * c2;
      if (u0 > 0.0) push(2, u0);
    }
  }
  const int m = (int)w.tn.size();
  if (m == 0) return;
  w.lp.resize(m); w.cdf.resize(m);
  std::vector<double>& thc = w.wtn;                  // column-major copy (network order)
  thc.resize((size_t)m * 4);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < 4; ++j) thc[(size_t)j * m + i] = th[(size_t)i * 4 + j];
  nle_ddm_native_full(s.ddm, thc.data(), m, s.tf.data(), s.tfc.data(), w.tn.data(), w.Rr.data(),
                      w.lp.data(), w.cdf.data());
  for (int i = 0; i < m; ++i) {
    const int r = w.sc_row[i];
    if (w.sc_kind[i] == 0) {
      const double s2 = S[r] * S[r];
      ll_trial[r] = std::log(s2 / (A[r] * A[r])) + w.lp[i];
    } else if (w.sc_kind[i] == 1) {
      const bool has0 = i + 1 < m && w.sc_kind[i + 1] == 2 && w.sc_row[i + 1] == r;
      const double G0 = has0 ? w.cdf[i + 1] : 0.0;
      ll_trial[r] = std::log(std::max(w.cdf[i] - G0, 0.0) / ST0[r]);
    }
  }
}

// Trial log-likelihoods of one mapped particle into ll_trial (n_trials),
// before bounds and clamping. Reentrant.
inline void nn_trial_ll(const NnSpec& s, const ParamTable& pt, NnScratch& w,
                        double min_ll, double* ll_trial) {
  if (s.scaled) { nn_trial_ll_scaled(s, pt, w, ll_trial); return; }
  if (s.cat_kind) {
    // every row goes to the evaluator (censored rows have no rt; the evaluator
    // scores them from LC / UC), the convolution parameters alongside
    const int m = s.n_rows;
    w.th.resize((size_t)m * s.nc);
    for (int j = 0; j < s.nc; ++j) {
      const double* c = pt.base.colptr(s.col[j]);
      std::copy(c, c + m, w.th.data() + (size_t)j * m);
    }
    w.lp.resize(m);
    nle_cat_native(s.cat, w.th.data(), m, s.tf.data(), s.tfc.data(), s.rt,
                   pt.base.colptr(s.conv_mean_col), pt.base.colptr(s.conv_range_col),
                   s.k.data(), s.miss.data(), s.LC.data(), s.UC.data(), w.lp.data());
    std::copy(w.lp.begin(), w.lp.end(), ll_trial);
    return;
  }
  const double* pre = s.pre_col >= 0 ? pt.base.colptr(s.pre_col) : nullptr;
  // rows the network evaluates: time > 0 (else zero density / CDF, as the
  // R path's dfun/pfun) and, in a race, accumulators taking part
  // (a hybrid's rows with the hybrid parameter exactly 0 go to the Wald)
  const double* hyb = s.hyb_col >= 0 ? pt.base.colptr(s.hyb_col) : nullptr;
  w.rows.clear(); w.tn.clear(); w.wrows.clear(); w.wtn.clear();
  for (int r = 0; r < s.n_rows; ++r) {
    if (!s.absent.empty() && s.absent[r]) continue;
    const double t = pre ? s.rt[r] - pre[r] : s.rt[r];
    if (!(t > 0.0)) continue;
    if (hyb && hyb[r] == 0.0) { w.wrows.push_back(r); w.wtn.push_back(t); }
    else { w.rows.push_back(r); w.tn.push_back(t); }
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
  if (!w.wrows.empty()) {
    const double* wv = pt.base.colptr(s.wald_col[0]);
    const double* wb = pt.base.colptr(s.wald_col[1]);
    const double* ws = pt.base.colptr(s.wald_col[2]);
    for (size_t i = 0; i < w.wrows.size(); ++i) {
      const int r = w.wrows[i];
      const double p[3] = {wv[r], wb[r], ws[r]};
      bool in = true;                         // outside the box: as the network's rows
      for (int j = 0; j < 3; ++j)
        if (std::isnan(p[j]) || p[j] < s.wald_lower[j] || p[j] > s.wald_upper[j]) in = false;
      if (!in) continue;
      w.ll_row[r] = s.winner[r] ? nle::wald_log_pdf(w.wtn[i], p[0], p[1], p[2])
                                : std::log1p(-nle::wald_cdf(w.wtn[i], p[0], p[1], p[2]));
    }
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
