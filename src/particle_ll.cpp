#include <Rcpp.h>
#include <unordered_map>

// Utilities first — no dependencies on model types
#include "utility_functions.h"
#include "transform_utils.h"
#include "ParamTable.h"
#include "TrendEngine.h"
#include "math_utils.h"

// Model headers — each includes RaceSpec.h themselves
#include "model_LBA.h"
#include "model_lnr.h"
#include "model_RDM.h"
#include "model_DDM.h"
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
    ctx.trend_plan.reset(new TrendPlan(trend, data));
    ctx.trend_runtime.reset(new TrendRuntime(*ctx.trend_plan));
    ctx.trend_runtime->bind_all_ops_to_paramtable(ctx.param_table);

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
    for (TrendOpRuntime& op : trend_runtime->premap_ops) {
      trend_runtime->apply_base_for_op(op, param_table);
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
    for (TrendOpRuntime& op : trend_runtime->pretransform_ops) {
      trend_runtime->apply_base_for_op(op, param_table);
    }
  }

  // 4) Transforms for all parameters excluding trend pars used so far
  c_do_transform_pt(param_table, cache.postmap_specs);

  // 5) Posttransform trends
  if (trend_runtime && trend_runtime->has_posttransform()) {
    for (TrendOpRuntime& op : trend_runtime->posttransform_ops) {
      trend_runtime->apply_base_for_op(op, param_table);
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

double c_log_likelihood_stop_signal(NumericMatrix pars, DataFrame data,
                                    IntegerVector expand, double min_ll,
                                    const std::vector<int>& is_ok, const String& type) {
  NumericVector RT = data["rt"];
  IntegerVector R = data["R"];
  NumericVector SSD = data["SSD"];
  NumericVector lR = data["lR"];
  LogicalVector winner = data["winner"];
  NumericVector unique_lR = unique(lR);
  const int n_acc = unique_lR.length();
  const int n_trials = lR.length() / n_acc;

  LogicalVector trial_ok(n_trials);
  for (int trial = 0; trial < n_trials; ++trial) {
    trial_ok[trial] = is_ok[trial * n_acc];
  }

  NumericVector ll_trial;
  if (type == "SSEXG") {
    ll_trial = ss_texg_lpdf(RT, R, SSD, lR, winner, pars, trial_ok, min_ll);
  } else {
    ll_trial = ss_rdex_lpdf(RT, R, SSD, lR, winner, pars, trial_ok, min_ll);
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
    const int n_lR = unique(lR).length();
    for (int i = 0; i < n_particles; ++i) {
      std::fill(is_ok.begin(), is_ok.end(), 1);
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);
      lr_all(is_ok, n_lR);
      lls[i] = c_log_likelihood_stop_signal(pars, data, expand, min_ll, is_ok, type);
    }
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
NumericMatrix get_pars_c_wrapper(NumericMatrix particle_matrix,
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

  // Shared setup
  PipelineContext ctx = make_pipeline_context(particle_matrix, data, constants,
                                              designs, transforms, pretransforms, trend);
  TrendRuntime* trend_runtime_ptr = ctx.trend_runtime ? ctx.trend_runtime.get() : nullptr;

  // Pipeline cache
  PipelineCache cache = make_pipeline_cache(ctx.param_table, designs, ctx.transform_specs, trend_runtime_ptr);

  // kernel_output_codes: IntegerVector -> std::vector<int>
  std::vector<int> kernel_codes(kernel_output_codes.begin(), kernel_output_codes.end());

  // Run pipeline (single particle — no loop needed)
  run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);

  // Extract and return
  if (return_kernel_matrix) {
    return get_covariate_matrix(ctx.param_table, trend_runtime_ptr, kernel_codes);
  } else if (return_all_pars) {
    return get_all_pars(ctx.param_table);
  } else {
    return get_pars_matrix(ctx.param_table, ctx.keep_names);
  }
}
