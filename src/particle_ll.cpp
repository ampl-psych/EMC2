#include <Rcpp.h>
#include <unordered_map>

// Utilities first — no dependencies on model types
#include "utility_functions.h"
#include "transform_utils.h"
#include "ParamTable.h"
#include "TrendEngine.h"
#include "math_utils.h"
#include "composite_functions.h"

// Model headers — each includes RaceSpec.h themselves
#include "model_LBA.h"
#include "model_lnr.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_CDM.h"
#include "model_MRI.h"
#include "model_SS_EXG.h"
#include "model_SS_RDEX.h"

// RaceSetup last — references functions defined in model headers above
#include "RaceSetup.h"
using namespace Rcpp;


// =============================================================================
// PipelineCache — pre-computed specs and masks for the parameter pipeline
// =============================================================================

struct PipelineCache {
  std::unordered_set<std::string> postmap_param_set;
  std::vector<TransformSpec>      postmap_specs;
  std::vector<TransformSpec>      premap_specs;       // empty if no premap trend
  std::vector<TransformSpec>      pretransform_specs; // empty if no pretransform trend

  Rcpp::LogicalVector mask_premap;          // regular premap designs
  Rcpp::LogicalVector mask_premap_reparam;  // reparam targets that are premap
  Rcpp::LogicalVector mask_map;             // regular main designs
  Rcpp::LogicalVector mask_reparam;         // reparam in main step
};

PipelineCache make_pipeline_cache(
    const ParamTable& param_table,
    const Rcpp::List& designs,
    const std::vector<TransformSpec>& transform_specs,
    TrendRuntime* trend_runtime_ptr)
{
  static const std::unordered_set<std::string> empty_set;

  PipelineCache cache;

  const auto& premap_set       = trend_runtime_ptr ? trend_runtime_ptr->premap_trend_params()       : empty_set;
  const auto& pretransform_set = trend_runtime_ptr ? trend_runtime_ptr->pretransform_trend_params() : empty_set;

  cache.postmap_param_set = param_names_excluding(param_table, { &premap_set, &pretransform_set });
  cache.postmap_specs     = filter_specs_by_param_set(param_table, transform_specs, cache.postmap_param_set);

  if (trend_runtime_ptr && trend_runtime_ptr->has_premap()) {
    cache.premap_specs = filter_specs_by_param_set(param_table, transform_specs, premap_set);
  }
  if (trend_runtime_ptr && trend_runtime_ptr->has_pretransform()) {
    cache.pretransform_specs = filter_specs_by_param_set(param_table, transform_specs, pretransform_set);
  }

  Rcpp::CharacterVector dnames = designs.names();
  const int n_designs = dnames.size();

  // Figure out which parameters are *targets* for reparameterisations
  std::unordered_set<std::string> reparam_set;
  for (int i = 0; i < n_designs; ++i) {
    Rcpp::RObject dm = designs[i];
    Rcpp::RObject pd_attr = dm.attr("parameter_design");
    if (!Rf_isNull(pd_attr) && Rcpp::as<bool>(pd_attr)) {
      reparam_set.insert(Rcpp::as<std::string>(dnames[i]));
    }
  }

  cache.mask_premap         = Rcpp::LogicalVector(n_designs, false);
  cache.mask_premap_reparam = Rcpp::LogicalVector(n_designs, false);
  cache.mask_map            = Rcpp::LogicalVector(n_designs, false);
  cache.mask_reparam        = Rcpp::LogicalVector(n_designs, false);

  if (trend_runtime_ptr && trend_runtime_ptr->has_premap()) {
    Rcpp::LogicalVector base_premap = trend_runtime_ptr->premap_design_mask(designs);
    const auto& premap_pars = trend_runtime_ptr->premap_trend_params();

    for (int i = 0; i < n_designs; ++i) {
      std::string nm = Rcpp::as<std::string>(dnames[i]);
      bool is_rep = (reparam_set.count(nm) > 0);
      bool is_pre = base_premap[i] || (is_rep && premap_pars.count(nm) > 0);

      if      ( is_rep &&  is_pre) cache.mask_premap_reparam[i] = true;
      else if (!is_rep &&  is_pre) cache.mask_premap[i]         = true;
      else if ( is_rep && !is_pre) cache.mask_reparam[i]        = true;
      else                         cache.mask_map[i]            = true;
    }
  } else {
    for (int i = 0; i < n_designs; ++i) {
      std::string nm = Rcpp::as<std::string>(dnames[i]);
      bool is_rep = (reparam_set.count(nm) > 0);
      if (is_rep) cache.mask_reparam[i] = true;
      else        cache.mask_map[i]     = true;
    }
  }

  return cache;
}


// =============================================================================
// PipelineContext — live runtime state, owns objects for the particle loop lifetime
// =============================================================================

struct PipelineContext {
  Rcpp::NumericMatrix            particle_matrix;   // after pretransform + constants
  ParamTable                     param_table;
  std::vector<TransformSpec>     transform_specs;
  std::unique_ptr<TrendPlan>     trend_plan;
  std::unique_ptr<TrendRuntime>  trend_runtime;
  Rcpp::CharacterVector          keep_names;
  std::vector<int>               pm_col_to_base_idx;
};

PipelineContext make_pipeline_context(
    Rcpp::NumericMatrix particle_matrix,
    const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& constants,
    const Rcpp::List& designs,
    const Rcpp::List& transforms,
    const Rcpp::List& pretransforms,
    const Rcpp::Nullable<Rcpp::List>& trend)
{
  PipelineContext ctx;

  // 1. Pre-transform
  std::vector<TransformSpec> t_specs = make_transform_specs_matrix(particle_matrix, pretransforms);
  ctx.particle_matrix = c_do_transform_matrix(particle_matrix, t_specs);

  // 2. Append constants
  bool has_constants = !(constants.size() == 1 && Rcpp::NumericVector::is_na(constants[0]));
  if (has_constants) {
    ctx.particle_matrix = add_constants_columns(ctx.particle_matrix, constants);
  }

  // 3. Build ParamTable from first particle
  Rcpp::NumericVector p_vector = ctx.particle_matrix(0, Rcpp::_);
  p_vector.attr("names") = colnames(ctx.particle_matrix);
  ctx.param_table = ParamTable::from_p_vector_and_designs(p_vector, designs, data.nrow());

  // 4. Transform specs
  ctx.transform_specs = make_transform_specs_pt(ctx.param_table, transforms);

  // 5. Trend objects and keep_names
  if (!trend.isNull()) {
    ctx.trend_plan.reset(new TrendPlan(Rcpp::List(trend.get()), data));
    ctx.trend_runtime.reset(new TrendRuntime(*ctx.trend_plan));
    ctx.trend_runtime->bind_all_to_paramtable(ctx.param_table);

    Rcpp::CharacterVector dnames = designs.names();
    const auto& trend_params = ctx.trend_runtime->all_trend_params();
    ctx.keep_names = names_excluding(dnames, { &trend_params });
  } else {
    ctx.keep_names = designs.names();
  }

  // 6. Column-index lookup: particle matrix column -> ParamTable base index
  Rcpp::CharacterVector pm_names = colnames(ctx.particle_matrix);
  ctx.pm_col_to_base_idx.assign(pm_names.size(), -1);
  for (int j = 0; j < pm_names.size(); ++j) {
    std::string nm = Rcpp::as<std::string>(pm_names[j]);
    auto it = ctx.param_table.name_to_base_idx.find(nm);
    if (it != ctx.param_table.name_to_base_idx.end()) {
      ctx.pm_col_to_base_idx[j] = it->second;
    }
  }

  return ctx;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix do_transform(Rcpp::NumericMatrix pars, Rcpp::List transform) {
  // Build the specs for these parameters
  std::vector<TransformSpec> specs = make_transform_specs_matrix(pars, transform);
  // Apply transformation in place and return
  return c_do_transform_matrix(pars, specs);
}


// =============================================================================
// run_pars_pipeline — runs steps 3-7 in place on param_table
// =============================================================================

void run_pars_pipeline(ParamTable& param_table,
                       const Rcpp::List& designs,
                       TrendRuntime* trend_runtime,
                       const PipelineCache& cache)
{
  if (trend_runtime) {
    // 0) Ensure kernels are reset
    trend_runtime->reset_all_kernels();
  }

  // 1) Premap trends: MAP premap trend parameters, TRANSFORM them, RUN kernels+bases
  if (trend_runtime && trend_runtime->has_premap()) {
    param_table.map_from_designs(designs, cache.mask_premap);
    param_table.map_from_designs(designs, cache.mask_premap_reparam);
    if (!cache.premap_specs.empty()) {
      c_do_transform_pt(param_table, cache.premap_specs);
    }
    for (BaseRuntime& base : trend_runtime->premap_bases) {
      trend_runtime->apply_base(base, param_table);
    }
  }

  // 2) Map designs for remaining parameters
  param_table.map_from_designs(designs, cache.mask_map);
  param_table.map_from_designs(designs, cache.mask_reparam);

  // 3) Pretransform trends: TRANSFORM pretransform trend parameters, RUN kernels+bases
  if (trend_runtime && trend_runtime->has_pretransform()) {
    if (!cache.pretransform_specs.empty()) {
      c_do_transform_pt(param_table, cache.pretransform_specs);
    }
    for (BaseRuntime& base : trend_runtime->pretransform_bases) {
      trend_runtime->apply_base(base, param_table);
    }
  }

  // 4) Transforms for all parameters excluding trend pars used so far
  c_do_transform_pt(param_table, cache.postmap_specs);

  // 5) Posttransform trends
  if (trend_runtime && trend_runtime->has_posttransform()) {
    for (BaseRuntime& base : trend_runtime->posttransform_bases) {
      trend_runtime->apply_base(base, param_table);
    }
  }
}

// =============================================================================
// Extractors — call after run_pars_pipeline
// =============================================================================

NumericMatrix get_pars_matrix(ParamTable& param_table,
                              const Rcpp::CharacterVector& keep_names)
{
  return param_table.materialize_by_param_names(keep_names);
}

NumericMatrix get_all_pars(ParamTable& param_table)
{
  return param_table.materialize();
}

NumericMatrix get_covariate_matrix(ParamTable& param_table,
                                   TrendRuntime* trend_runtime,
                                   const std::vector<int>& kernel_output_codes)
{
  if (!trend_runtime) {
    Rcpp::stop("return_kernel_matrix/return_covariate_matrix requested but no trend was provided");
  }
  std::vector<int> codes = kernel_output_codes;
  if (codes.empty()) codes.push_back(1);  // default: main trajectory
  return trend_runtime->all_kernel_outputs(param_table, codes);
}


// =============================================================================
// Likelihood functions — Call within calc_ll branches
// =============================================================================

double c_log_likelihood_race(ParamTable& pt,
                             const RaceModelSetup& setup,
                             const NumericVector& rts,
                             const LogicalVector& winner,
                             const std::vector<int>& is_ok,
                             const std::vector<int>& idx_win,
                             const std::vector<int>& idx_los,
                             const IntegerVector& expand,
                             double min_ll,
                             int n_acc,
                             NumericVector& ll_row,
                             NumericVector& ll_trial,
                             RaceScratch& scratch)
{
  const int n_winners = (int)idx_win.size();

  double* ll_row_ptr = ll_row.begin();
  double* ll_ptr  = ll_trial.begin();
  const int* ok_ptr = is_ok.data();

  // 1) Fill log(pdf) for winners and log(1-cdf) for losers into ll_row.
  //    fill_both stores pdf / (1-cdf); vec_log transforms the whole array
  //    in one vectorised pass (vvlog on Apple, libmvec on Linux/x86).
  //    Invalid inputs (<=0, nan) produce -inf or nan, which the clamp below
  //    catches — no per-element branching needed.
  //
  //   // setup.fill_both() refers to gather-scatter implementations.
  //   // on linux/x86, this is significantly faster. macOS/arm64 doesn't care

  setup.fill_both(rts, pt, setup.spec, idx_win, idx_los, ll_row_ptr, scratch);
  vec_log(ll_row_ptr, ll_row.size());  // bulk log over entire ll_row buffer

  // 2) Per-trial log-likelihood into ll_trial.
  //    ll_row now contains log(pdf) at winner indices, log(1-cdf) at loser indices.
  //    Clamp to min_ll using !(v > min_ll) which catches -inf and nan.
  auto clamp = [min_ll](double v) {
    return (v > min_ll) ? v : min_ll;
  };

  if (n_acc == 1) {
    for (int t = 0; t < n_winners; ++t) {
      const int i_win = idx_win[t];
      ll_ptr[t] = ok_ptr[i_win] ? clamp(ll_row_ptr[i_win]) : min_ll;
    }
  } else {
    for (int t = 0; t < n_winners; ++t) {
      const int base = t * n_acc;

      if (!ok_ptr[idx_win[t]]) {
        // lr_all guarantees that ok_ptr are the same value for all accumulators in a trial
        // so only check here, no need to check for the other accumulators
        ll_ptr[t] = min_ll;
        continue;
      }

      // The current data format guarantees n_acc per trial, so we can just sum now
      // ll_row_ptr contains either the log-PDF (winners) or log(1-CDF) (losers)
      // Clamp here. There's a second clamp later on but not really needed probably
      double ll = 0.0;
      for (int k = 0; k < n_acc; ++k) {
        ll += clamp(ll_row_ptr[base + k]);
      }
      ll_ptr[t] = ll;
    }
  }

  // 3) Expand and sum
  const int  m       = expand.size();
  const int* exp_ptr = expand.begin();
  double sum_ll = 0.0;

#pragma omp simd reduction(+:sum_ll)
  for (int i = 0; i < m; ++i) {
    sum_ll += clamp(ll_ptr[exp_ptr[i] - 1]);
  }
  return sum_ll;
}


// Mirror of c_log_likelihood_DDM but for CDM (R is numeric angle)
double c_log_likelihood_CDM(NumericMatrix pars, DataFrame data,
                                   const int n_trials, IntegerVector expand,
                                   double min_ll, const std::vector<int>& is_ok){
  NumericVector rts = data["rt"]; // numeric
  NumericVector Rs  = data["R"];  // numeric angles
  CharacterVector dnames = data.names();
  const bool has_R2 = sum(contains(dnames, "R2")) == 1;
  const bool has_R3 = sum(contains(dnames, "R3")) == 1;
  LogicalVector is_ok_r(is_ok.size());
  std::copy(is_ok.begin(), is_ok.end(), is_ok_r.begin());
  NumericVector lls(n_trials);
  NumericVector lls_exp(expand.length());
  if (has_R2 && has_R3 && pars.ncol() >= 8) {
    NumericVector R2s = data["R2"];
    NumericVector R3s = data["R3"];
    lls = c_dHSDM(rts, Rs, R2s, R3s, pars, is_ok_r);
  } else if (has_R2 && pars.ncol() >= 7) {
    NumericVector R2s = data["R2"];
    lls = c_dSDM(rts, Rs, R2s, pars, is_ok_r);
  } else {
    lls = c_dCDM(rts, Rs, pars, is_ok_r);
  }
  lls_exp = c_expand(lls, expand); // decompress
  lls_exp[is_na(lls_exp)] = min_ll;
  lls_exp[is_infinite(lls_exp)] = min_ll;
  lls_exp[lls_exp < min_ll] = min_ll;
  return sum(lls_exp);
}

double c_log_likelihood_PSDM(NumericMatrix pars, DataFrame data,
                             const int n_trials, IntegerVector expand,
                             double min_ll, const std::vector<int>& is_ok){
  NumericVector rts = data["rt"];
  NumericVector Rs = data["R"];
  LogicalVector is_ok_r(is_ok.size());
  std::copy(is_ok.begin(), is_ok.end(), is_ok_r.begin());
  NumericVector lls(n_trials);
  NumericVector lls_exp(expand.length());
  lls = c_dPSDM(rts, Rs, pars, is_ok_r);
  lls_exp = c_expand(lls, expand);
  lls_exp[is_na(lls_exp)] = min_ll;
  lls_exp[is_infinite(lls_exp)] = min_ll;
  lls_exp[lls_exp < min_ll] = min_ll;
  return sum(lls_exp);
}

double c_log_likelihood_PHSDM(NumericMatrix pars, DataFrame data,
                              const int n_trials, IntegerVector expand,
                              double min_ll, const std::vector<int>& is_ok){
  NumericVector rts = data["rt"];
  NumericVector Rs = data["R"];
  NumericVector R2s = data["R2"];
  LogicalVector is_ok_r(is_ok.size());
  std::copy(is_ok.begin(), is_ok.end(), is_ok_r.begin());
  NumericVector lls(n_trials);
  NumericVector lls_exp(expand.length());
  lls = c_dPHSDM(rts, Rs, R2s, pars, is_ok_r);
  lls_exp = c_expand(lls, expand);
  lls_exp[is_na(lls_exp)] = min_ll;
  lls_exp[is_infinite(lls_exp)] = min_ll;
  lls_exp[lls_exp < min_ll] = min_ll;
  return sum(lls_exp);
}

// SS helper pointer types
using ss_go_pdf_fn = NumericVector (*)(NumericVector, NumericMatrix, LogicalVector, double);
using ss_stop_surv_fn = double (*)(double, NumericMatrix);
using ss_stop_success_fn = double (*)(double, NumericMatrix, double, double, int, double, double, double, double);

// Model-specific stop survivor wrappers (read fixed columns)
static inline double stop_logsurv_texg_fn(double q, NumericMatrix P) {
  // EXG stop: muS=3, sigmaS=4, tauS=5, exgS_lb=9
  return ptexg(q, P(0, 3), P(0, 4), P(0, 5), P(0, 9), R_PosInf, false, true);
}
static inline double stop_logsurv_rdex_fn(double q, NumericMatrix P) {
  // RDEX stop: muS=5, sigmaS=6, tauS=7, exgS_lb=10
  return ptexg(q, P(0, 5), P(0, 6), P(0, 7), P(0, 10), R_PosInf, false, true);
}

double c_log_likelihood_ss(
    NumericMatrix pars,
    DataFrame data,
    const int n_trials,
    IntegerVector expand,
    double min_ll,
    LogicalVector is_ok,
    ss_go_pdf_fn go_lpdf_ptr,
    ss_go_pdf_fn go_lccdf_ptr,
    ss_stop_surv_fn stop_logsurv_ptr,
    ss_stop_success_fn stop_success_ptr,
    int idx_tf,
    int idx_gf
) {
  // initialise local variables
  const int n_out = expand.length();
  if (is_true(all(!is_ok))) {
    NumericVector lls_expanded(n_out, min_ll);
    return(sum(lls_expanded));
  }
  NumericVector lls(n_trials);
  NumericVector lls_expanded(n_out);
  // extract data
  NumericVector RT = data["rt"];
  IntegerVector R = data["R"];
  NumericVector SSD = data["SSD"];
  NumericVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  bool has_lI = data.containsElementNamed("lI");
  IntegerVector lI = has_lI ? as<IntegerVector>(data["lI"]) : IntegerVector(lR.size(), 2);

  // dimensional expectations: pars has one row per accumulator per trial

  // compute log likelihoods (generalized, matching R's log_likelihood_race_ss)
  NumericVector unique_lR = unique(lR);
  const int n_acc = unique_lR.length();
  // n_trials equals data rows grouped by accumulators
  for (int trial = 0; trial < n_trials; trial++) {
    if (is_ok[trial] != 1) { lls[trial] = min_ll; continue; }

    int start_row = trial * n_acc;
    int end_row   = (trial + 1) * n_acc - 1;
    // basic bounds are guaranteed by correct n_trials passed into this function
    NumericMatrix P = pars(Range(start_row, end_row), _);
    IntegerVector lI_trial = lI[Range(start_row, end_row)];
    LogicalVector is_go(n_acc, true), is_st(n_acc, false);
    // determine go/ST accumulators if present
    if (has_lI) {
      int go_code = max(lI_trial);
      for (int i = 0; i < n_acc; i++) {
        is_go[i] = (lI_trial[i] == go_code);
        is_st[i] = !is_go[i];
      }
    } else {
      for (int i = 0; i < n_acc; i++) is_st[i] = false;
    }
    int n_accG = sum(is_go);
    int n_accST = sum(is_st);

    double tf = P(0, idx_tf);
    double gf = P(0, idx_gf);

    double rt = RT[start_row];
    bool response_observed = R[start_row] != NA_INTEGER;
    bool stop_signal_presented = std::isfinite(SSD[start_row]);

    // Identify whether observed response is GO or ST (when response observed)
    bool response_is_go = false;
    if (response_observed) {
      int r_obs = R[start_row];
      for (int i = 0; i < n_acc; i++) {
        if (lR[start_row + i] == r_obs) {
          response_is_go = is_go[i];
          break;
        }
      }
    }

    // Build rt vectors for go and st contexts
    NumericVector rt_go(n_acc, rt);
    NumericVector rt_st(n_acc, rt - SSD[start_row]);

    // GO masks for current trial
    LogicalVector win_mask = winner[Range(start_row, end_row)];
    LogicalVector go_win_mask(n_acc); // winner and go
    LogicalVector go_loss_mask(n_acc);
    for (int i = 0; i < n_acc; i++) {
      go_win_mask[i] = (win_mask[i] && is_go[i]);
      go_loss_mask[i] = (!win_mask[i] && is_go[i]);
    }

    if (!response_observed) {
      // No response
      if (!stop_signal_presented) {
        lls[trial] = std::log(gf);
      } else if (n_accST == 0) {
        // Stop trial, no ST accumulators: gf + (1-gf)*(1-tf)*pStop
        // Early-skip: if mixture weight (1-gf)*(1-tf) is tiny, skip integral
        double mix_w = (1.0 - gf) * (1.0 - tf);
        double comp1 = std::log(gf);
        if (mix_w <= 1e-6) {
          lls[trial] = comp1; // effectively just gf component
        } else {
          NumericMatrix P_go = submat_rcpp(P, is_go);
          double log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, R_PosInf,
                                              50, 1e-6, 1e-5, 8.0, 16.0);
          double comp2 = log1m(gf) + log1m(tf) + log_pstop;
          lls[trial] = log_sum_exp(comp1, comp2);
        }
      } else {
        // Not handled in R either; keep minimal ll
        lls[trial] = min_ll;
      }
      continue;
    }

    // Response observed
    if (!stop_signal_presented) {
      // GO trial with response: (1-gf) * GO race ll
      double go_lprob = 0.0;
      NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
      go_lprob = (lw.size() > 0) ? sum(lw) : min_ll;
      if (n_accG > 1) {
        NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
        for (int i = 0; i < ls.size(); i++) go_lprob += (R_FINITE(ls[i]) ? ls[i] : min_ll);
      }
      lls[trial] = log1m(gf) + go_lprob;
      continue;
    }

    // Stop trial with response
    if (response_is_go) {
      // GO wins on stop trial: (1-gf) * [ tf * go + (1-tf) * (go + stop_surv + st_loss) ]
      // go race ll at observed rt
      double go_lprob = 0.0;
      NumericVector lw = go_lpdf_ptr(rt_go, P, go_win_mask, min_ll);
      go_lprob = (lw.size() > 0) ? sum(lw) : min_ll;
      if (n_accG > 1) {
        NumericVector ls = go_lccdf_ptr(rt_go, P, go_loss_mask, min_ll);
        for (int i = 0; i < ls.size(); i++) go_lprob += (R_FINITE(ls[i]) ? ls[i] : min_ll);
      }
      // stop survivor at observed rt (rt - SSD)
      double log_stop_surv = stop_logsurv_ptr(rt - SSD[start_row], P);
      if (!R_FINITE(log_stop_surv)) log_stop_surv = min_ll;
      // ST losers survivors (if any)
      double st_loss_sum = 0.0;
      if (n_accST > 0) {
        LogicalVector st_loss_mask(n_acc);
        for (int i = 0; i < n_acc; i++) st_loss_mask[i] = (is_st[i] && !win_mask[i]);
        NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
        for (int i = 0; i < ls_st.size(); i++) st_loss_sum += (R_FINITE(ls_st[i]) ? ls_st[i] : min_ll);
      }
      double comp_tf = go_lprob; // only go
      double comp_notf = go_lprob + log_stop_surv + st_loss_sum; // fair race (stop loses)
      lls[trial] = log1m(gf) + log_mix(tf, comp_tf, comp_notf);
      continue;
    } else {
      // ST wins on stop trial
      // ST winner log pdf at rt - SSD
      LogicalVector st_win_mask(n_acc);
      for (int i = 0; i < n_acc; i++) st_win_mask[i] = (win_mask[i] && is_st[i]);
      NumericVector lw_st = go_lpdf_ptr(rt_st, P, st_win_mask, min_ll);
      double st_winner_logpdf = (lw_st.size() > 0) ? sum(lw_st) : min_ll;
      // ST losers survivors
      double st_loss_sum = 0.0;
      if (n_accST > 1) {
        LogicalVector st_loss_mask(n_acc);
        for (int i = 0; i < n_acc; i++) st_loss_mask[i] = (!win_mask[i] && is_st[i]);
        NumericVector ls_st = go_lccdf_ptr(rt_st, P, st_loss_mask, min_ll);
        for (int i = 0; i < ls_st.size(); i++) st_loss_sum += (R_FINITE(ls_st[i]) ? ls_st[i] : min_ll);
      }
      // GO losers survivors
      double go_loss_sum = 0.0;
      if (n_accG > 0) {
        NumericVector ls_go = go_lccdf_ptr(rt_go, P, is_go, min_ll);
        for (int i = 0; i < ls_go.size(); i++) go_loss_sum += (R_FINITE(ls_go[i]) ? ls_go[i] : min_ll);
      }
      // Stop success probability up to observed rt (only go racers influence integral)
      NumericMatrix P_go = submat_rcpp(P, is_go);
      double log_pstop = stop_success_ptr(SSD[start_row], P_go, min_ll, rt,
                                          50, 1e-6, 1e-5, 8.0, 16.0);

      double st_base = st_winner_logpdf + st_loss_sum;
      // mixture over gf and pStop, never tf when ST wins
      double term_gf = std::log(gf) + st_base; // go failure -> only ST race
      double term_stop_win = log1m(gf) + log_pstop + st_base; // stop beats go -> only ST race
      double term_stop_lose = log1m(gf) + log1m_exp(log_pstop) + st_base + go_loss_sum; // all race, no stop win
      lls[trial] = log1m(tf) + log_sum_exp(term_gf, log_sum_exp(term_stop_win, term_stop_lose));
      continue;
    }
  }
  lls[is_na(lls)] = min_ll;
  lls[is_infinite(lls)] = min_ll;
  lls[lls < min_ll] = min_ll;
  lls_expanded = c_expand(lls, expand); // decompress
  return(sum(lls_expanded));
}


double c_log_likelihood_DDM(NumericMatrix pars, DataFrame data,
                            const int n_trials, IntegerVector expand,
                            double min_ll, std::vector<int> is_ok){
  const int n_out = expand.length();
  NumericVector rts = data["rt"];
  IntegerVector R = data["R"];
  NumericVector lls(n_trials);
  lls = d_DDM_Wien(rts, R, pars, is_ok);

  // lls_exp = c_expand(lls, expand); // decompress
  // // lls_exp = lls;
  // lls_exp[is_na(lls_exp)] = min_ll;
  // lls_exp[is_infinite(lls_exp)] = min_ll;
  // lls_exp[lls_exp < min_ll] = min_ll;
  // return(sum(lls_exp));
  // More SIMD-friendly == faster
  // decompress

  const double* lls_ptr    = lls.begin();
  const int*    expand_ptr = expand.begin();

  double sum_ll = 0.0;

  // expand is 1-based, so subtract 1
#pragma omp simd reduction(+:sum_ll)
  for (int i = 0; i < n_out; ++i) {
    int idx = expand_ptr[i] - 1;
    double v = lls_ptr[idx];

    if (!std::isfinite(v) || v < min_ll) {
      v = min_ll;
    }
    sum_ll += v;
  }

  return sum_ll;
}


int c_col_index(const CharacterVector& names, const std::string& target) {
  for (int i = 0; i < names.size(); ++i) {
    if (Rcpp::as<std::string>(names[i]) == target) return i;
  }
  stop("Column not found: " + target);
}

NumericVector c_expand_ordered_cut(NumericVector raw_cut, int n_lR) {
  if (raw_cut.size() % n_lR != 0) {
    stop("cut vector length must be divisible by the number of response levels");
  }

  NumericVector cut = clone(raw_cut);
  const int n_trials = cut.size() / n_lR;

  for (int trial = 0; trial < n_trials; ++trial) {
    const int base = trial * n_lR;
    if (n_lR == 2) {
      cut[base + 1] = cut[base];
      continue;
    }

    double current = cut[base];
    for (int r = 1; r < n_lR - 1; ++r) {
      current += std::exp(raw_cut[base + r]);
      cut[base + r] = current;
    }
    cut[base + n_lR - 1] = cut[base + n_lR - 2];
  }

  return cut;
}

double c_ordered_cdf(double x, double location, double scale, bool probit) {
  if (x == R_NegInf) return 0.0;
  if (x == R_PosInf) return 1.0;
  if (probit) return R::pnorm(x, location, scale, true, false);
  return R::plogis(x, location, scale, true, false);
}

double c_log_likelihood_ordered(NumericMatrix pars, DataFrame data,
                                const int n_lR, IntegerVector expand,
                                double min_ll, std::vector<int> is_ok, bool probit) {
  const int loc_idx = c_col_index(colnames(pars), "location");
  const int scale_idx = c_col_index(colnames(pars), "scale");
  const int cut_idx = c_col_index(colnames(pars), "cut");
  const LogicalVector winner = data["winner"];
  const IntegerVector lR = data["lR"];
  const NumericVector cut = c_expand_ordered_cut(pars(_, cut_idx), n_lR);

  NumericVector ll_trial(sum(winner));
  int out_idx = 0;

  for (int i = 0; i < pars.nrow(); ++i) {
    if (!winner[i]) continue;

    if (is_ok[i] != 1) {
      ll_trial[out_idx++] = min_ll;
      continue;
    }

    const int level = lR[i];
    const double location = pars(i, loc_idx);
    const double scale = pars(i, scale_idx);
    const double upper = (level == n_lR) ? R_PosInf : cut[i];
    const double lower = (level == 1) ? R_NegInf : cut[i - 1];
    const double prob = c_ordered_cdf(upper, location, scale, probit) -
      c_ordered_cdf(lower, location, scale, probit);

    double ll = min_ll;
    if (R_FINITE(prob) && prob > 0) {
      ll = std::log(prob);
      if (!R_FINITE(ll) || ll < min_ll) ll = min_ll;
    }
    ll_trial[out_idx++] = ll;
  }

  NumericVector ll_exp = c_expand(ll_trial, expand);
  ll_exp[is_na(ll_exp)] = min_ll;
  ll_exp[is_infinite(ll_exp)] = min_ll;
  ll_exp[ll_exp < min_ll] = min_ll;
  return sum(ll_exp);
}

double c_log_likelihood_multinomial_logit(NumericMatrix pars, DataFrame data,
                                          const int n_lR, IntegerVector expand,
                                          double min_ll, std::vector<int> is_ok) {
  const int utility_idx = c_col_index(colnames(pars), "utility");
  const LogicalVector winner = data["winner"];
  const int n_trials = pars.nrow() / n_lR;
  NumericVector ll_trial(n_trials);

  for (int trial = 0; trial < n_trials; ++trial) {
    const int base = trial * n_lR;
    if (is_ok[base] != 1) {
      ll_trial[trial] = min_ll;
      continue;
    }

    double max_utility = pars(base, utility_idx);
    for (int r = 1; r < n_lR; ++r) {
      const double value = pars(base + r, utility_idx);
      if (value > max_utility) max_utility = value;
    }

    double denom = 0.0;
    double chosen = NA_REAL;
    for (int r = 0; r < n_lR; ++r) {
      const double value = std::exp(pars(base + r, utility_idx) - max_utility);
      denom += value;
      if (winner[base + r]) chosen = value;
    }

    double ll = min_ll;
    if (R_FINITE(denom) && denom > 0 && R_FINITE(chosen) && chosen > 0) {
      ll = std::log(chosen / denom);
      if (!R_FINITE(ll) || ll < min_ll) ll = min_ll;
    }
    ll_trial[trial] = ll;
  }

  NumericVector ll_exp = c_expand(ll_trial, expand);
  ll_exp[is_na(ll_exp)] = min_ll;
  ll_exp[is_infinite(ll_exp)] = min_ll;
  ll_exp[ll_exp < min_ll] = min_ll;
  return sum(ll_exp);
}



// [[Rcpp::export]]
NumericVector calc_ll(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                      List designs, String type, List bounds, List transforms, List pretransforms,
                      CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue) {
  const int n_particles = particle_matrix.nrow();
  const int n_trials    = data.nrow();

  NumericVector  lls(n_particles);
  std::vector<int> is_ok(n_trials, 1);

  // Shared setup -- context holds the param_table as well as designs, constants, trend etc
  PipelineContext ctx = make_pipeline_context(particle_matrix, data, constants,
                                              designs, transforms, pretransforms, trend);
  TrendRuntime* trend_runtime_ptr = ctx.trend_runtime ? ctx.trend_runtime.get() : nullptr;

  // Bounds — built once from structure, not values
  NumericMatrix   minmax   = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<BoundSpec> bound_specs = make_bound_specs_pt(minmax, mm_names, ctx.param_table, bounds);

  PipelineCache cache = make_pipeline_cache(ctx.param_table, designs,
                                            ctx.transform_specs, trend_runtime_ptr);


  // -----------------------------------------------------------------------
  // DDM
  // -----------------------------------------------------------------------
  if (type == "DDM") {
    IntegerVector expand = data.attr("expand");
    for (int i = 0; i < n_particles; ++i) {
      std::fill(is_ok.begin(), is_ok.end(), 1);
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok);
    }
  } else if(type == "CDM"){
    IntegerVector expand = data.attr("expand");
    for (int i = 0; i < n_particles; ++i) {
      std::fill(is_ok.begin(), is_ok.end(), 1);
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);
      lls[i] = c_log_likelihood_CDM(pars, data, n_trials, expand, min_ll, is_ok);
    }
  } else if(type == "PSDM"){
    IntegerVector expand = data.attr("expand");
    for (int i = 0; i < n_particles; ++i) {
      std::fill(is_ok.begin(), is_ok.end(), 1);
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);
      lls[i] = c_log_likelihood_PSDM(pars, data, n_trials, expand, min_ll, is_ok);
    }
  } else if(type == "PHSDM"){
    IntegerVector expand = data.attr("expand");
    for (int i = 0; i < n_particles; ++i) {
      std::fill(is_ok.begin(), is_ok.end(), 1);
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);
      lls[i] = c_log_likelihood_PHSDM(pars, data, n_trials, expand, min_ll, is_ok);
    }
  } else if(type == "ORDERED_PROBIT" || type == "ORDERED_LOGIT"){
    IntegerVector expand = data.attr("expand");
    IntegerVector lR = data["lR"];
    const int n_lR = unique(lR).length();
    const bool is_probit = (type == "ORDERED_PROBIT");
    for (int i = 0; i < n_particles; ++i) {
      std::fill(is_ok.begin(), is_ok.end(), 1);
      // Fill from particle row
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      // Run parameter mapping pipeline
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      // This one still requires a Rcpp::NumericMatrix, can't operate on the param_table directly
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // Do bound in-place
      lr_all(is_ok, n_lR);   // also in-place
      lls[i] = c_log_likelihood_ordered(pars, data, n_lR, expand, min_ll, is_ok, is_probit);
    }
  } else if(type == "MULTINOMIAL_LOGIT"){
    IntegerVector expand = data.attr("expand");
    IntegerVector lR = data["lR"];
    const int n_lR = unique(lR).length();
    for (int i = 0; i < n_particles; ++i) {
      std::fill(is_ok.begin(), is_ok.end(), 1);
      // Fill from particle row
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      // Run parameter mapping pipeline
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      // This one still requires a Rcpp::NumericMatrix, can't operate on the param_table directly
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // Do bound in-place
      lr_all(is_ok, n_lR);   // also in-place
      lls[i] = c_log_likelihood_multinomial_logit(pars, data, n_lR, expand, min_ll, is_ok);
    }
  // -----------------------------------------------------------------------
  // MRI / MRI_AR1
  // -----------------------------------------------------------------------
  } else if (type == "MRI" || type == "MRI_AR1") {
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    const bool is_ar1 = (type == "MRI_AR1");
    for (int i = 0; i < n_particles; ++i) {
      std::fill(is_ok.begin(), is_ok.end(), 1);
      // Fill from particle row
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      // Run parameter mapping pipeline
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      // This one still requires a Rcpp::NumericMatrix, can't operate on the param_table directly
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // Do bound in-place
      lls[i] = is_ar1 ? c_log_likelihood_MRI_white(pars, y, is_ok, n_trials, n_pars, min_ll)
        : c_log_likelihood_MRI(pars, y, is_ok, n_trials, n_pars, min_ll);
      }
  } else if(type == "SSEXG" || type == "SSRDEX"){
    IntegerVector expand = data.attr("expand");
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    int n_trials_ss = (n_lR > 0) ? (n_trials / n_lR) : n_trials;
    // Pick function pointers and indices based on type
    ss_go_pdf_fn go_lpdf_ptr = (type == "SSEXG") ? texg_go_lpdf : rdex_go_lpdf;
    ss_go_pdf_fn go_lccdf_ptr = (type == "SSEXG") ? texg_go_lccdf : rdex_go_lccdf;
    ss_stop_surv_fn stop_logsurv_ptr = (type == "SSEXG") ? stop_logsurv_texg_fn : stop_logsurv_rdex_fn;
    ss_stop_success_fn stop_success_ptr = (type == "SSEXG") ? ss_texg_stop_success_lpdf : ss_rdex_stop_success_lpdf;
    int idx_tf = (type == "SSEXG") ? 6 : 8;
    int idx_gf = (type == "SSEXG") ? 7 : 9;
    for (int i = 0; i < n_particles; ++i) {
      p_vector = p_matrix(i, _);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend, full_t_specs);
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR); // reduce to per-trial ok
      lls[i] = c_log_likelihood_ss(pars, data, n_trials_ss, expand, min_ll, is_ok,
                                   go_lpdf_ptr, go_lccdf_ptr,
                                   stop_logsurv_ptr, stop_success_ptr,
                                   idx_tf, idx_gf);
    }
  } else{
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    // Love me some good old ugly but fast c++ pointers
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    if(type == "LBA"){
      dfun = dlba_c;
      pfun = plba_c;
    } else if(type == "RDM"){
      dfun = drdm_c;
      pfun = prdm_c;
    } else{
      dfun = dlnr_c;
      pfun = plnr_c;
    }
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
  // -----------------------------------------------------------------------
  // Race models (RDM, LBA, LNR)
  // -----------------------------------------------------------------------
  } else {
    NumericVector lR     = data["lR"];
    IntegerVector expand = data.attr("expand");
    const int n_acc      = unique(lR).length();
    NumericVector rts    = data["rt"];
    LogicalVector winner = data["winner"];

    // Precompute winner/loser index lists (once, outside particle loop)
    std::vector<int> idx_win, idx_los;
    idx_win.reserve(n_trials);
    idx_los.reserve(n_trials);
    int* win_flag = LOGICAL(winner);

    // Pre-read RACE info needed for phantom filtering
    const bool has_race_col = (sum(contains(data.names(), "RACE")) == 1);
    NumericVector   NACC;
    CharacterVector vals_NACC;
    if (has_race_col) {
      NACC      = data["RACE"];
      vals_NACC = NACC.attr("levels");
    }

    // Identify which rows in the dadm correspond to winners, to losers, and which should be skipped entirely
    for (int i = 0; i < n_trials; ++i) {
      if (win_flag[i]) {
        idx_win.push_back(i);
      } else {
        // skip phantom accumulators — data-dependent, built once
        if (has_race_col && lR[i] > atoi(vals_NACC[NACC[i] - 1])) continue;
        idx_los.push_back(i);
      }
    }
    const int n_winners = (int)idx_win.size();

    // Scratch buffers (reused across particles)
    NumericVector ll_row(n_trials);  // stores (log)likelihood of row in dadm
    NumericVector ll_trial(n_winners); // stores (log)likelihood of trials

    RaceScratch scratch;
    scratch.reserve(std::max((int)idx_win.size(), (int)idx_los.size()));

    // Race model setup
    RaceModelSetup setup = make_race_setup(type, ctx.param_table);

    // Begin particle loop
    for (int i = 0; i < n_particles; ++i) {
      std::fill(is_ok.begin(), is_ok.end(), 1);
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);
      lr_all(is_ok, n_acc);
      std::fill(ll_row.begin(), ll_row.end(), 1.0);
      lls[i] = c_log_likelihood_race(
        ctx.param_table, setup,  // operates directly on param_table - no need for param extraction
        rts, winner, is_ok,
        idx_win, idx_los, expand,
        min_ll, n_acc, ll_row, ll_trial,
        scratch);
    }
  }

  return lls;
}


// [[Rcpp::export]]
List get_pars_c_wrapper(NumericMatrix particle_matrix,
                        DataFrame data,
                        NumericVector constants,
                        List designs,
                        List bounds,
                        List transforms,
                        List pretransforms,
                        Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                        bool return_kernel_matrix = false,
                        bool return_all_pars = false,
                        IntegerVector kernel_output_codes = 1)
{
  if (Rf_isNull(colnames(particle_matrix))) {
    stop("p_matrix must have column names for pretransforms/transform specs");
  }
  const int n_trials    = data.nrow();
  const int n_particles = particle_matrix.nrow();
  const bool has_lR     = (sum(contains(data.names(), "lR")) == 1);
  const int n_lR = has_lR ? unique(Rcpp::as<IntegerVector>(data["lR"])).length() : 1;

  // Shared setup
  PipelineContext ctx = make_pipeline_context(particle_matrix, data, constants,
                                              designs, transforms, pretransforms, trend);
  TrendRuntime* trend_runtime_ptr = ctx.trend_runtime ? ctx.trend_runtime.get() : nullptr;

  // Pipeline cache (built once, reused across particles)
  PipelineCache cache = make_pipeline_cache(ctx.param_table, designs,
                                            ctx.transform_specs, trend_runtime_ptr);

  // kernel_output_codes: IntegerVector -> std::vector<int>
  std::vector<int> kernel_codes(kernel_output_codes.begin(), kernel_output_codes.end());

  NumericMatrix   minmax   = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<BoundSpec> bound_specs = make_bound_specs_pt(minmax, mm_names, ctx.param_table, bounds);
  std::vector<int> is_ok(n_trials, 1);

  List result(n_particles);

  for (int i = 0; i < n_particles; ++i) {
    // Update param_table for this particle (skip refill on first particle)
    if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
    run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);

    std::fill(is_ok.begin(), is_ok.end(), 1);
    c_do_bound_pt(ctx.param_table, bound_specs, is_ok);
    if(n_lR > 1) lr_all(is_ok, n_lR);

    // Convert is_ok to LogicalVector for the attribute
    LogicalVector ok_attr(is_ok.begin(), is_ok.end());
    // Extract and store
    NumericMatrix mat;
    if (return_kernel_matrix) {
      mat = get_covariate_matrix(ctx.param_table, trend_runtime_ptr, kernel_codes);
    } else if (return_all_pars) {
      mat = get_all_pars(ctx.param_table);
    } else {
      mat = get_pars_matrix(ctx.param_table, ctx.keep_names);
    }

    mat.attr("ok") = ok_attr;
    result[i] = mat;
  }

  return result;
}
