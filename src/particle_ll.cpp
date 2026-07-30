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


// Helper functions for ordered choice models
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

// =============================================================================
// expand_clamp_sum — hoisted final step for all models with expand
// =============================================================================

inline double expand_clamp_sum(const double* ll_ptr, const int* exp_ptr,
                               const int m, const double min_ll,
                               double* trialwise_out = nullptr)
{
  double sum_ll = 0.0;
  if (trialwise_out) {
    for (int i = 0; i < m; ++i) {
      const double v = ll_ptr[exp_ptr[i] - 1];
      const double c = (v > min_ll) ? v : min_ll;
      trialwise_out[i] = c;
      sum_ll += c;
    }
  } else {
#pragma omp simd reduction(+:sum_ll)
    for (int i = 0; i < m; ++i) {
      const double v = ll_ptr[exp_ptr[i] - 1];
      sum_ll += (v > min_ll) ? v : min_ll;
    }
  }
  return sum_ll;
}

// MRI variant — no expand
inline double clamp_sum(const double* ll_ptr, const int n, const double min_ll,
                        double* trialwise_out = nullptr)
{
  double sum_ll = 0.0;
  if (trialwise_out) {
    for (int i = 0; i < n; ++i) {
      const double c = (ll_ptr[i] > min_ll) ? ll_ptr[i] : min_ll;
      trialwise_out[i] = c;
      sum_ll += c;
    }
  } else {
#pragma omp simd reduction(+:sum_ll)
    for (int i = 0; i < n; ++i) {
      sum_ll += (ll_ptr[i] > min_ll) ? ll_ptr[i] : min_ll;
    }
  }
  return sum_ll;
}



// =============================================================================
// Likelihood functions — Call within calc_ll branches
// =============================================================================

void c_log_likelihood_race(ParamTable& pt,
                           const RaceModelSetup& setup,
                           const NumericVector& rts,
                           const LogicalVector& winner,
                           // const std::vector<int>& is_ok,
                           const std::vector<int>& idx_win,
                           const std::vector<int>& idx_los,
                           int n_acc,
                           NumericVector& ll_row,
                           double* ll_buf,
                           RaceScratch& scratch)
{
  const int n_winners = (int)idx_win.size();

  double* ll_row_ptr = ll_row.begin();
  // const int* ok_ptr = is_ok.data();

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
  if (n_acc == 1) {
    for (int t = 0; t < n_winners; ++t) {
      const int i_win = idx_win[t];
      ll_buf[t] = ll_row_ptr[i_win];
    }
  } else {
    for (int t = 0; t < n_winners; ++t) {
      const int base = t * n_acc;

      // The current data format guarantees n_acc per trial, so we can just sum now
      // ll_row_ptr contains either the log-PDF (winners) or log(1-CDF) (losers)
      // Clamp here. There's a second clamp later on but not really needed probably
      // SM changed his mind. do *not* clamp here. If we clamp here, a bad value (e.g., -inf) will be corrected and the whole trial might become min_ll + log(pdf) [winner!],
      // whereas arguably it should be min_ll
      // Only clamp at the end while reducing.
      double ll = 0.0;
      for (int k = 0; k < n_acc; ++k) {
        ll += ll_row_ptr[base + k];
      }
      ll_buf[t] = ll;
    }
  }
}


void c_log_likelihood_DDM(NumericMatrix pars, DataFrame data,
                            const int n_trials,
                            double* ll_buf) {
  NumericVector rts = data["rt"];
  IntegerVector R = data["R"];
  NumericVector lls = d_DDM_Wien(rts, R, pars);

  const double* src = lls.begin();
  for (int i = 0; i < n_trials; ++i) ll_buf[i] = src[i];
}


void c_log_likelihood_ordered(NumericMatrix pars, DataFrame data,
                              const int n_lR,    // std::vector<int>& is_ok,
                              bool probit,
                              double* ll_buf)    // size = n_winners
{
  const int loc_idx   = c_col_index(colnames(pars), "location");
  const int scale_idx = c_col_index(colnames(pars), "scale");
  const int cut_idx   = c_col_index(colnames(pars), "cut");
  const LogicalVector winner = data["winner"];
  const IntegerVector lR     = data["lR"];
  const NumericVector cut    = c_expand_ordered_cut(pars(_, cut_idx), n_lR);

  int out_idx = 0;
  for (int i = 0; i < pars.nrow(); ++i) {
    if (!winner[i]) continue;

    const int    level    = lR[i];
    const double location = pars(i, loc_idx);
    const double scale    = pars(i, scale_idx);
    const double upper    = (level == n_lR) ? R_PosInf : cut[i];
    const double lower    = (level == 1)    ? R_NegInf : cut[i - 1];
    const double prob     = c_ordered_cdf(upper, location, scale, probit) - c_ordered_cdf(lower, location, scale, probit);

    ll_buf[out_idx++] = std::log(prob);
  }
}

void c_log_likelihood_multinomial_logit(NumericMatrix pars, DataFrame data,
                                        const int n_lR,  // std::vector<int>& is_ok,
                                        double* ll_buf)  // size = n_choice_trials
{
  const int utility_idx  = c_col_index(colnames(pars), "utility");
  const LogicalVector winner = data["winner"];
  const int n_choice_trials  = pars.nrow() / n_lR;

  for (int trial = 0; trial < n_choice_trials; ++trial) {
    const int base = trial * n_lR;

    double max_utility = pars(base, utility_idx);
    for (int r = 1; r < n_lR; ++r) {
      const double v = pars(base + r, utility_idx);
      if (v > max_utility) max_utility = v;
    }

    double denom = 0.0, chosen = NA_REAL;
    for (int r = 0; r < n_lR; ++r) {
      const double v = std::exp(pars(base + r, utility_idx) - max_utility);
      denom += v;
      if (winner[base + r]) chosen = v;
    }

    ll_buf[trial] = std::log(chosen / denom);  // -inf/nan caught by expand_clamp_sum
  }
}

void c_log_likelihood_MRI(NumericMatrix pars, NumericVector y,
                          int n, int m,
                          double* ll_buf)             // size = n
{
  for (int i = 0; i < n; i++) {
    double s = 0.0;
    for (int j = 0; j < m - 1; j++) s += pars(i, j);
    ll_buf[i] = R::dnorm(y[i], s, pars(i, m - 1), true);
  }
}


void c_log_likelihood_MRI_white(NumericMatrix pars, NumericVector y,
                                int n, int m,
                                double* ll_buf)        // size = n
{
  // First observation: stationary variance
  double s = 0.0;
  for (int j = 0; j < m - 2; j++) s += pars(0, j);
  ll_buf[0] = R::dnorm(y[0], s, pars(0, m - 1), true);

  // t >= 1: AR(1) conditional likelihood
  for (int i = 1; i < n; i++) {
    double s = 0.0;
    for (int j = 0; j < m - 2; j++) s += pars(i, j);
    const double sigma_i   = pars(i, m - 1);
    const double rho_i     = pars(i, m - 2);
    const double cond_sd   = sigma_i * std::sqrt(1.0 - rho_i * rho_i);
    const double cond_mean = s + rho_i * (y[i - 1] - ([&]() {
      double s_prev = 0.0;
      for (int j = 0; j < m - 2; j++) s_prev += pars(i - 1, j);
      return s_prev;
    })());
    ll_buf[i] = R::dnorm(y[i], cond_mean, cond_sd, true);
  }
}


// [[Rcpp::export]]
NumericMatrix calc_ll(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                      List designs, String type, List bounds, List transforms, List pretransforms,
                      CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                      bool return_trialwise = false) {

  const int n_particles = particle_matrix.nrow();
  const int n_trials    = data.nrow();
  const bool has_lR         = (sum(contains(data.names(), "lR")) == 1);
  const int n_lR            = has_lR ? unique(IntegerVector(data["lR"])).length() : 1;
  const int n_choice_trials = n_trials / n_lR;
  const int out_rows        = return_trialwise ? n_choice_trials : 1;

  // Column i = particle i's trialwise buffer (contiguous)
  NumericMatrix result(out_rows, n_particles);  // column i = particle i, contiguous

  std::vector<int>     is_ok(n_trials, 1);
  std::vector<double>  ll_buf(n_choice_trials);  // compressed scratch, reused per particle


  // Shared setup -- context holds the param_table as well as designs, constants, trend etc
  PipelineContext ctx = make_pipeline_context(particle_matrix, data, constants,
                                              designs, transforms, pretransforms, trend);
  TrendRuntime* trend_runtime_ptr = ctx.trend_runtime ? ctx.trend_runtime.get() : nullptr;
  // Bounds — built once from structure
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
    const int     n_exp  = expand.size();
    for (int i = 0; i < n_particles; ++i) {
      // Map p_vector to trialwise parameters
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);

      // calculate raw (compressed) trialwise log-likelihoods (fills ll_buf)
      c_log_likelihood_DDM(pars, data, n_trials, ll_buf.data());

      // Check for bound violations (fills is_ok)
      std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
      for (int t = 0; t < n_choice_trials; ++t) if(!is_ok[t]) ll_buf[t] = min_ll;  // apply min_ll when is_ok = false (overwrite ll_buf)

      // Determine output location (tw is a pointer to the correct address in result) and protect via expand, clamp, sum
      double* tw = return_trialwise ? result.column(i).begin() : nullptr;
      const double sum = expand_clamp_sum(ll_buf.data(), expand.begin(), n_exp, min_ll, tw);
      if (!return_trialwise) result(0, i) = sum;
      }

  // -----------------------------------------------------------------------
  // Ordered probit / logit
  // -----------------------------------------------------------------------
  } else if (type == "ORDERED_PROBIT" || type == "ORDERED_LOGIT") {
    IntegerVector  expand    = data.attr("expand");
    const int      n_exp     = expand.size();
    const bool     is_probit = (type == "ORDERED_PROBIT");
    for (int i = 0; i < n_particles; ++i) {
      // Map p_vector to trialwise parameters
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);

      // calculate raw (compressed) trialwise log-likelihoods (fills ll_buf)
      c_log_likelihood_ordered(pars, data, n_lR, is_probit, ll_buf.data());

      // Handle not-ok parameter values (out of bound)
      std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok [parameter dependent]
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
      lr_all(is_ok, n_lR);                                 // make sure the ok value is shared across accumulators / levels or lR
      for (int t = 0; t < n_choice_trials; ++t) if(!is_ok[t * n_lR]) ll_buf[t] = min_ll;  // apply min_ll when is_ok = false (overwrite ll_buf)

      // Determine output location (tw is a pointer to the correct address in result) and protect via expand, clamp, sum
      double* tw = return_trialwise ? result.column(i).begin() : nullptr;
      const double sum = expand_clamp_sum(ll_buf.data(), expand.begin(), n_exp, min_ll, tw);
      if (!return_trialwise) result(0, i) = sum;
    }

  // -----------------------------------------------------------------------
  // Multinomial logit
  // -----------------------------------------------------------------------
  } else if (type == "MULTINOMIAL_LOGIT") {
    IntegerVector expand = data.attr("expand");
    const int     n_exp  = expand.size();
    for (int i = 0; i < n_particles; ++i) {
      // Map p_vector to trialwise parameters
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);

      // calculate raw (compressed) trialwise log-likelihoods (fills ll_buf)
      c_log_likelihood_multinomial_logit(pars, data, n_lR, ll_buf.data());

      // Handle not-ok parameter values (out of bound)
      std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok [parameter dependent]
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
      lr_all(is_ok, n_lR);                                 // make sure the ok value is shared across accumulators / levels or lR
      for (int t = 0; t < n_choice_trials; ++t) if(!is_ok[t * n_lR]) ll_buf[t] = min_ll;  // apply min_ll when is_ok = false (overwrite ll_buf)

      // Determine output location (tw is a pointer to the correct address in result) and protect via expand, clamp, sum
      double* tw = return_trialwise ? result.column(i).begin() : nullptr;
      const double sum = expand_clamp_sum(ll_buf.data(), expand.begin(), n_exp, min_ll, tw);
      if (!return_trialwise) result(0, i) = sum;
    }

  // -----------------------------------------------------------------------
  // MRI / MRI_AR1
  // -----------------------------------------------------------------------
  } else if (type == "MRI" || type == "MRI_AR1") {
    const int     n_pars = p_types.length();
    NumericVector y      = extract_y(data);
    const bool    is_ar1 = (type == "MRI_AR1");
    for (int i = 0; i < n_particles; ++i) {
      // Map p_vector to trialwise parameters
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);

      // Fill log-likelihood buffer
      if (is_ar1) c_log_likelihood_MRI(pars, y, n_trials, n_pars, ll_buf.data());
      else        c_log_likelihood_MRI_white(pars, y, n_trials, n_pars, ll_buf.data());

      // Check for bound violations (fills is_ok)
      std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
      for (int t = 0; t < n_choice_trials; ++t) if(!is_ok[t]) ll_buf[t] = min_ll;  // apply min_ll when is_ok = false (overwrite ll_buf)

      // Determine output location (tw is a pointer to the correct address in result) and protect via clamp, sum [[no expanding here, compression doesn't work for MRI]]
      double* tw = return_trialwise ? result.column(i).begin() : nullptr;
      const double sum = clamp_sum(ll_buf.data(), n_trials, min_ll, tw);
      if (!return_trialwise) result(0, i) = sum;
    }

  // -----------------------------------------------------------------------
  // Race models (RDM, LBA, LNR)
  // -----------------------------------------------------------------------
  } else {
    NumericVector lR     = data["lR"];
    IntegerVector expand = data.attr("expand");
    const int n_exp      = expand.size();
    NumericVector rts    = data["rt"];
    LogicalVector winner = data["winner"];
    int*           win_flag = LOGICAL(winner);

    // Pre-read RACE info needed for phantom filtering
    const bool has_race_col = (sum(contains(data.names(), "RACE")) == 1);
    NumericVector   NACC;     // to-do: reimplement this with missingness = 4 instead of assumig the value corresponds to NACC
    CharacterVector vals_NACC;
    if (has_race_col) {
      NACC      = data["RACE"];
      vals_NACC = NACC.attr("levels");
    }

    // Precompute winner/loser index lists (once, outside particle loop)
    std::vector<int> idx_win, idx_los;
    idx_win.reserve(n_trials);
    idx_los.reserve(n_trials);
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
    // const int n_winners = (int)idx_win.size();

    // Scratch buffers (reused across particles)
    NumericVector ll_row(n_trials);                 // stores (log)likelihood of row in dadm
    // ll_trial sized n_choice_trials; assumes n_winners <= n_choice_trials (one winner per trial at most)
    std::vector<double> ll_trial(n_choice_trials);  // compressed scratch for (log)likelihoods in race (compressed! so needs expanding)

    // Race model setup
    RaceScratch scratch;
    scratch.reserve(std::max((int)idx_win.size(), (int)idx_los.size()));
    RaceModelSetup setup = make_race_setup(type, ctx.param_table);

    // Begin particle loop
    for (int i = 0; i < n_particles; ++i) {
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, designs, trend_runtime_ptr, cache);

      std::fill(ll_row.begin(), ll_row.end(), 1.0); // re-fill row-wise ll -- helps with RACE functionality (all trials skipped get likelihood = 1)
      c_log_likelihood_race(ctx.param_table, setup, rts, winner,
                            idx_win, idx_los, n_lR, ll_row, ll_trial.data(), scratch);

      // Handle not-ok parameter values (out of bound)
      std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok [parameter dependent]
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
      lr_all(is_ok, n_lR);                                 // make sure the ok value is shared across accumulators / levels or lR
      for (int t = 0; t < n_choice_trials; ++t) if(!is_ok[t * n_lR]) ll_trial[t] = min_ll;  // apply min_ll when is_ok = false (overwrite ll_trial)

      double* tw = return_trialwise ? result.column(i).begin() : nullptr;
      const double sum = expand_clamp_sum(ll_trial.data(), expand.begin(), n_exp, min_ll, tw);
      if (!return_trialwise) result(0, i) = sum;
    }
  }

  return result;
  // NumericVector out(n_particles);
  // for(int i = 0; i < n_particles; i++) out[i] = result(0,i);
  // return(out);
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

