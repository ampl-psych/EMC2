#ifdef _OPENMP
#include <omp.h>
#endif

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
#include "CensorSpec.h"
#include "TruncSpec.h"
using namespace Rcpp;


// =============================================================================
// PipelineCache — pre-computed specs and masks for the parameter pipeline
// =============================================================================

struct PipelineCache {
  std::unordered_set<std::string> postmap_param_set;
  std::vector<TransformSpec>      postmap_specs;
  std::vector<TransformSpec>      premap_specs;       // empty if no premap trend
  std::vector<TransformSpec>      pretransform_specs; // empty if no pretransform trend

  // Masks — std::vector<bool> so safe inside OpenMP regions
  std::vector<bool> mask_premap;            // regular premap designs
  std::vector<bool> mask_premap_reparam;    // reparam targets that are premap
  std::vector<bool> mask_map;               // regular main designs
  std::vector<bool> mask_reparam;           // reparam in main step
};

PipelineCache make_pipeline_cache(
    ParamTable& param_table,
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

  // Initialise design plans
  // param_table.init_design_plan(designs);

  // Figure out which parameters are *targets* for reparameterisations
  std::unordered_set<std::string> reparam_set;
  for (int i = 0; i < n_designs; ++i) {
    Rcpp::RObject dm = designs[i];
    Rcpp::RObject pd_attr = dm.attr("parameter_design");
    if (!Rf_isNull(pd_attr) && Rcpp::as<bool>(pd_attr)) {
      reparam_set.insert(Rcpp::as<std::string>(dnames[i]));
    }
  }

  // Initialise all masks to false
  cache.mask_premap.assign(n_designs, false);
  cache.mask_premap_reparam.assign(n_designs, false);
  cache.mask_map.assign(n_designs, false);
  cache.mask_reparam.assign(n_designs, false);

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

  // --- Constant-column flags ---
  // Start pessimistic
  std::fill(param_table.col_is_constant.begin(),
            param_table.col_is_constant.end(), false);

  // Pass 1: direct design targets are non-constant
  // (skip_self_intercept entries are intercept-only and stay constant)
  for (const DesignEntry& entry : param_table.design_plan) {
    if (!entry.valid) continue;
    if (entry.dm_is_constant) param_table.col_is_constant[entry.out_idx] = true;
  }

  // Trend targets are non-constant
  if (trend_runtime_ptr) {
    for (const std::string& nm : trend_runtime_ptr->all_trend_targets()) {
      auto it = param_table.name_to_base_idx.find(nm);
      if (it != param_table.name_to_base_idx.end())
        param_table.col_is_constant[it->second] = false;
    }
  }

  // Pass 2: reparam targets inherit non-constancy from their inputs
  for (int i = 0; i < n_designs; ++i) {
    if (!cache.mask_reparam[i]) continue;
    const DesignEntry& entry = param_table.design_plan[i];
    if (!entry.valid || !entry.dm_is_constant) continue;

    bool all_inputs_constant = true;
    for (int cidx : entry.coef_idx) {
      if (cidx >= 0 && !param_table.col_is_constant[cidx]) {
        all_inputs_constant = false;
        break;
      }
    }
    if (all_inputs_constant)
      param_table.col_is_constant[entry.out_idx] = true;
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

void run_pars_pipeline(ParamTable&          param_table,
                       TrendRuntime*        trend_runtime,
                       const PipelineCache& cache)
  {
  if (trend_runtime) {
    // 0) Ensure kernels are reset
    trend_runtime->reset_all_kernels();
  }

  // 1) Premap trends: MAP premap trend parameters, TRANSFORM them, RUN kernels+bases
  if (trend_runtime && trend_runtime->has_premap()) {
    param_table.map_from_designs(cache.mask_premap);
    param_table.map_from_designs(cache.mask_premap_reparam);
    if (!cache.premap_specs.empty()) {
      c_do_transform_pt(param_table, cache.premap_specs);
    }
    for (BaseRuntime& base : trend_runtime->premap_bases) {
      trend_runtime->apply_base(base, param_table);
    }
  }

  // 2) Map designs for remaining parameters
  param_table.map_from_designs(cache.mask_map);
  param_table.map_from_designs(cache.mask_reparam);

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
                           const double* rts_ptr,
                           const std::vector<int>& idx_win,
                           const std::vector<int>& idx_los,
                           int n_acc,
                           double* ll_row_ptr,
                           int     ll_row_size,
                           double* ll_buf,
                           RaceScratch& scratch)
{
  const int n_winners = (int)idx_win.size();

  setup.fill_both(rts_ptr, pt, setup.spec,
                  idx_win, idx_los, ll_row_ptr, scratch);
  vec_log(ll_row_ptr, ll_row_size);

  if (n_acc == 1) {
    for (int t = 0; t < n_winners; ++t)
      ll_buf[idx_win[t]] = ll_row_ptr[idx_win[t]];
  } else {
    for (int t = 0; t < n_winners; ++t) {
      // NB: ll_buf is sized n_trials - including the missing (NA) trials. Can't index with t directly...
      const int trial_idx = idx_win[t] / n_acc;
      const int base      = trial_idx * n_acc;
      double ll = 0.0;
      for (int k = 0; k < n_acc; ++k)
        ll += ll_row_ptr[base + k];
      ll_buf[trial_idx] = ll;
    }
  }
}

void c_log_likelihood_DDM(const double* rts,
                          const int* R,
                          const ParamTable& pt,
                          const RaceSpec& spec,
                          const std::vector<int>& idx_all,
                          double* __restrict__ ll_row)
{
  // this is perhaps a bit pointless to have? only for naming convention...
  fill_ddm(rts, R, pt,
           spec,
           idx_all,
           ll_row);
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

void c_log_likelihood_MRI_white(const ParamTable& pt,
                                const MRISpec& spec,
                                const double* __restrict__ y,
                                int n,
                                double* __restrict__ ll_buf)
{
  const double* sigma = pt.base.colptr(spec.col_sigma);

  // Pre-fetch all mean column pointers
  const int nm = spec.n_mean_cols;
  std::vector<const double*> mean_cols(nm);
  for (int j = 0; j < nm; ++j)
    mean_cols[j] = pt.base.colptr(spec.col_means[j]);

  for (int i = 0; i < n; i++) {
    double s = 0.0;
    for (int j = 0; j < nm; j++) s += mean_cols[j][i];
    ll_buf[i] = R::dnorm(y[i], s, sigma[i], true);
  }
}

void c_log_likelihood_MRI_ar1(const ParamTable& pt,
                              const MRISpec& spec,
                              const double* __restrict__ y,
                              int n,
                              double* __restrict__ ll_buf)
{
  const double* sigma = pt.base.colptr(spec.col_sigma);
  const double* rho   = pt.base.colptr(spec.col_rho);

  const int nm = spec.n_mean_cols;
  std::vector<const double*> mean_cols(nm);
  for (int j = 0; j < nm; ++j)
    mean_cols[j] = pt.base.colptr(spec.col_means[j]);

  // t = 0: stationary
  double s_prev = 0.0;
  for (int j = 0; j < nm; j++) s_prev += mean_cols[j][0];
  ll_buf[0] = R::dnorm(y[0], s_prev, sigma[0], true);

  // t >= 1: AR(1) conditional
  for (int i = 1; i < n; i++) {
    double s_curr = 0.0;
    for (int j = 0; j < nm; j++) s_curr += mean_cols[j][i];

    const double rho_i    = rho[i];
    const double cond_sd  = sigma[i] * std::sqrt(1.0 - rho_i * rho_i);
    const double cond_mean = s_curr + rho_i * (y[i - 1] - s_prev);

    ll_buf[i] = R::dnorm(y[i], cond_mean, cond_sd, true);
    s_prev = s_curr;
  }
}

// void c_log_likelihood_MRI_white(NumericMatrix& pars, NumericVector& y,
//                           int n, int m,
//                           double* ll_buf)             // size = n
// {
//   for (int i = 0; i < n; i++) {
//     double s = 0.0;
//     for (int j = 0; j < m - 1; j++) s += pars(i, j);
//     ll_buf[i] = R::dnorm(y[i], s, pars(i, m - 1), true);
//   }
// }
//
//
//
// void c_log_likelihood_MRI_ar1(NumericMatrix& pars, NumericVector& y,
//                               int n, int m,
//                               double* ll_buf)        // size = n
// {
//   // First observation: stationary variance
//   double s = 0.0;
//   for (int j = 0; j < m - 2; j++) s += pars(0, j);
//   ll_buf[0] = R::dnorm(y[0], s, pars(0, m - 1), true);
//
//   // t >= 1: AR(1) conditional likelihood
//   for (int i = 1; i < n; i++) {
//     double s = 0.0;
//     for (int j = 0; j < m - 2; j++) s += pars(i, j);
//     const double sigma_i   = pars(i, m - 1);
//     const double rho_i     = pars(i, m - 2);
//     const double cond_sd   = sigma_i * std::sqrt(1.0 - rho_i * rho_i);
//     const double cond_mean = s + rho_i * (y[i - 1] - ([&]() {
//       double s_prev = 0.0;
//       for (int j = 0; j < m - 2; j++) s_prev += pars(i - 1, j);
//       return s_prev;
//     })());
//     ll_buf[i] = R::dnorm(y[i], cond_mean, cond_sd, true);
//   }
// }


// [[Rcpp::export]]
NumericMatrix calc_ll(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                      List designs, String type, List bounds, List transforms, List pretransforms,
                      CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                      bool return_trialwise = false) {

  const int n_particles = particle_matrix.nrow();
  const int n_rows      = data.nrow();
  const bool has_lR         = (sum(contains(data.names(), "lR")) == 1);
  const int n_lR            = has_lR ? unique(IntegerVector(data["lR"])).length() : 1;
  const int n_choice_trials = n_rows / n_lR;
  const int out_rows        = return_trialwise ? n_choice_trials : 1;

  // Column i = particle i's trialwise buffer (contiguous)
  NumericMatrix result(out_rows, n_particles);  // column i = particle i, contiguous

  std::vector<int>     is_ok(n_rows, 1);
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
    NumericVector rts    = data["rt"];
    const double* rts_ptr = rts.begin();
    IntegerVector R      = data["R"];
    const int*    Rs_ptr = R.begin();

    RaceModelSetup setup = make_race_setup(type, ctx.param_table);
    RaceScratch ddm_scratch;
    CensorSpec censor = make_censor_spec(data, n_rows, /* nacc = */ 1, setup, ctx.param_table, ddm_scratch);
    TruncSpec  trunc  = make_trunc_spec(data, n_rows, /* nacc = */ 1, setup, ctx.param_table, ddm_scratch);
    // test for missingness

    IntegerVector missingness;
    const bool has_missingness = (sum(contains(data.names(), "missingness")) == 1);
    if (has_missingness) missingness = data["missingness"];

    std::vector<int> idx_all;
    idx_all.reserve(n_rows);
    for (int i = 0; i < n_rows; ++i) {
      if (has_missingness && !IntegerVector::is_na(missingness[i])) continue;
      idx_all.push_back(i);
    }

    for (int i = 0; i < n_particles; ++i) {
      // 1) Map p_vector to trialwise parameters
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);

      // 2) calculate raw (compressed) trialwise log-likelihoods (fills ll_buf)
      std::fill(ll_buf.begin(), ll_buf.end(), 1.0); // re-fill row-wise ll -- helps with missngness functionality (all trials skipped get likelihood = 1)
      c_log_likelihood_DDM(rts_ptr, Rs_ptr, ctx.param_table, setup.spec, idx_all, ll_buf.data());

      // 3) Trialwise truncation correction. Needs to be in this order!
      if(trunc.any()) trunc.calculate_normalization_constant();
      if(censor.any()) censor.fill_censored_rows(trunc, ll_buf, min_ll);
      if(trunc.any()) for (int t = 0; t < trunc.n_trials; ++t) ll_buf[t] -= trunc.log_Z[t];

      // 5) Check for bound violations (fills is_ok)
      std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
      for (int t = 0; t < n_choice_trials; ++t) if(!is_ok[t]) ll_buf[t] = min_ll;  // apply min_ll when is_ok = false (overwrite ll_buf)

      // 6) Determine output location (tw is a pointer to the correct address in result) and protect via expand, clamp, sum
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
      // 1) Map p_vector to trialwise parameters
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);

      // 2) calculate raw (compressed) trialwise log-likelihoods (fills ll_buf)
      c_log_likelihood_ordered(pars, data, n_lR, is_probit, ll_buf.data());

      // 3) Handle not-ok parameter values (out of bound)
      std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok [parameter dependent]
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
      lr_all(is_ok, n_lR);                                 // make sure the ok value is shared across accumulators / levels or lR
      for (int t = 0; t < n_choice_trials; ++t) if(!is_ok[t * n_lR]) ll_buf[t] = min_ll;  // apply min_ll when is_ok = false (overwrite ll_buf)

      // 4) Determine output location (tw is a pointer to the correct address in result) and protect via expand, clamp, sum
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
      // 1) Map p_vector to trialwise parameters
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);
      NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);

      // 2) calculate raw (compressed) trialwise log-likelihoods (fills ll_buf)
      c_log_likelihood_multinomial_logit(pars, data, n_lR, ll_buf.data());

      // 3) Handle not-ok parameter values (out of bound)
      std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok [parameter dependent]
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
      lr_all(is_ok, n_lR);                                 // make sure the ok value is shared across accumulators / levels or lR
      for (int t = 0; t < n_choice_trials; ++t) if(!is_ok[t * n_lR]) ll_buf[t] = min_ll;  // apply min_ll when is_ok = false (overwrite ll_buf)

      // 4) Determine output location (tw is a pointer to the correct address in result) and protect via expand, clamp, sum
      double* tw = return_trialwise ? result.column(i).begin() : nullptr;
      const double sum = expand_clamp_sum(ll_buf.data(), expand.begin(), n_exp, min_ll, tw);
      if (!return_trialwise) result(0, i) = sum;
    }

  // -----------------------------------------------------------------------
  // MRI / MRI_AR1
  // -----------------------------------------------------------------------
  } else if (type == "MRI" || type == "MRI_AR1") {
    NumericVector y_rcpp = extract_y(data);
    const double* y      = y_rcpp.begin();
    const bool    is_ar1 = (type == "MRI_AR1");

    // Resolve column indices once — no allocation inside particle loop
    MRISpec spec = make_mri_spec(ctx.param_table, ctx.keep_names, is_ar1);

    for (int i = 0; i < n_particles; ++i) {
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);

      // No get_pars_matrix call — work directly on param_table
      if (is_ar1) c_log_likelihood_MRI_ar1  (ctx.param_table, spec, y, n_choice_trials, ll_buf.data());
      else        c_log_likelihood_MRI_white (ctx.param_table, spec, y, n_choice_trials, ll_buf.data());

      std::fill(is_ok.begin(), is_ok.begin() + n_choice_trials, 1);
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);
      for (int t = 0; t < n_choice_trials; ++t) if (!is_ok[t]) ll_buf[t] = min_ll;

      double* tw = return_trialwise ? result.column(i).begin() : nullptr;
      const double sum = clamp_sum(ll_buf.data(), n_choice_trials, min_ll, tw);
      if (!return_trialwise) result(0, i) = sum;
    }
    // const int     n_pars = p_types.length();
    // NumericVector y      = extract_y(data);
    // const bool    is_ar1 = (type == "MRI_AR1");
    // for (int i = 0; i < n_particles; ++i) {
    //   // 1) Map p_vector to trialwise parameters
    //   if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
    //   run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);
    //   NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);
    //
    //   // 2) Fill log-likelihood buffer
    //   if (is_ar1) c_log_likelihood_MRI_ar1(pars, y, n_rows, n_pars, ll_buf.data());
    //   else        c_log_likelihood_MRI_white(pars, y, n_rows, n_pars, ll_buf.data());
    //
    //   // 3) Check for bound violations (fills is_ok)
    //   std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok
    //   c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
    //   for (int t = 0; t < n_choice_trials; ++t) if(!is_ok[t]) ll_buf[t] = min_ll;  // apply min_ll when is_ok = false (overwrite ll_buf)
    //
    //   // 4) Determine output location (tw is a pointer to the correct address in result) and protect via clamp, sum [[no expanding here, compression doesn't work for MRI]]
    //   double* tw = return_trialwise ? result.column(i).begin() : nullptr;
    //   const double sum = clamp_sum(ll_buf.data(), n_rows, min_ll, tw);
    //   if (!return_trialwise) result(0, i) = sum;
    // }

  // -----------------------------------------------------------------------
  // Race models (RDM, LBA, LNR)
  // -----------------------------------------------------------------------
  } else {
    NumericVector lR     = data["lR"];
    IntegerVector expand = data.attr("expand");
    const int n_exp      = expand.size();
    const int* exp_ptr   = expand.begin();

    NumericVector rts    = data["rt"];
    const double* rts_ptr = rts.begin();

    LogicalVector winner = data["winner"];
    int* win_flag = LOGICAL(winner);

    // test for missingness
    IntegerVector missingness;
    const bool has_missingness = (sum(contains(data.names(), "missingness")) == 1);
    if (has_missingness) missingness = data["missingness"];

    // Precompute winner/loser index lists (once, outside particle loop)
    std::vector<int> idx_win, idx_los;
    idx_win.reserve(n_rows);
    idx_los.reserve(n_rows);

    // int total_n_winners = 0;
    for (int i = 0; i < n_rows; ++i) {
      if (has_missingness && !IntegerVector::is_na(missingness[i])) continue;
      if (win_flag[i]) {
        idx_win.push_back(i);
      } else {
        idx_los.push_back(i);
      }
    }

    // Scratch buffers (reused across particles)
    std::vector<double> ll_row(n_rows);                 // stores (log)likelihood of row in dadm
    // ll_trial sized n_choice_trials; assumes n_winners <= n_choice_trials (one winner per trial at most)
    std::vector<double> ll_trial(n_choice_trials);     // compressed scratch for (log)likelihoods in race (compressed! so needs expanding)

    // Race model setup
    RaceScratch scratch;
    scratch.reserve(n_rows);

    // Race model setup
    RaceModelSetup setup = make_race_setup(type, ctx.param_table);
    // Pre-compute index lists for censored and truncated trials
    CensorSpec censor = make_censor_spec(data, n_choice_trials, n_lR, setup, ctx.param_table, scratch);// make_censor_spec(data, n_rows);
    TruncSpec trunc = make_trunc_spec(data, n_choice_trials, n_lR, setup, ctx.param_table, scratch);

    // Begin particle loop
    for (int i = 0; i < n_particles; ++i) {
      // 1) Map p_vector to trialwise parameters
      if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);

      // 2) calculate raw (compressed) trialwise log-likelihoods (fills ll_buf)
      std::fill(ll_row.begin(), ll_row.end(), 1.0); // re-fill row-wise ll -- helps with RACE functionality (all trials skipped get likelihood = 1)
      c_log_likelihood_race(ctx.param_table, setup, rts_ptr,
                            idx_win, idx_los, n_lR,
                            ll_row.data(), (int)ll_row.size(), ll_trial.data(), scratch);

      // 3) Handle truncation, censoring
      if(trunc.any()) trunc.calculate_normalization_constant();
      if(censor.any()) censor.fill_censored_rows(trunc, ll_trial, min_ll);  // needs to have trunc know the normalization constant, hence this order
      if(trunc.any()) for (int t = 0; t < trunc.n_trials; ++t) ll_trial[t] -= trunc.log_Z[t];  // cannot be before censoring - censoring fills rows, doesn't add

      // 4) Handle not-ok parameter values (out of bound)
      std::fill(is_ok.begin(), is_ok.end(), 1);            // reset is_ok [parameter dependent]
      c_do_bound_pt(ctx.param_table, bound_specs, is_ok);  // fill is_ok by bound
      lr_all(is_ok, n_lR);                             // make sure the ok value is shared across accumulators / levels or lR
      for (int t = 0; t < (int)idx_win.size(); ++t) if(!is_ok[idx_win[t]]) ll_trial[idx_win[t] / n_lR] = min_ll;
      // for(int t=0; t < n_choice_trials; t++) {
      //   for(int k=0; k < n_lR; k++) {
      //     if(!is_ok[t * n_lR + k]) {
      //       // set to min_ll and move to next trial - any subsequent !is_ok checks are unnecessary
      //       ll_trial[t] = min_ll;
      //       break;  // no need to check remaining accumulators for this trial
      //     }
      //   }
      // }

      // 5) Expand, clamp, sum, etc
      double* tw = return_trialwise ? result.column(i).begin() : nullptr;
      const double sum = expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll, tw);
      if (!return_trialwise) result(0, i) = sum;
    }
  }

  return result;
}

// [[Rcpp::export]]
NumericMatrix calc_ll_multithreaded(NumericMatrix particle_matrix, DataFrame data,
                                    NumericVector constants, List designs, String type,
                                    List bounds, List transforms, List pretransforms,
                                    CharacterVector p_types, double min_ll,
                                    Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                                    bool return_trialwise = false,
                                    int n_threads = -1) {

  // Only race models are parallelised for now; all others fall back to serial calc_ll
  const bool is_race = (//type != "DDM"            &&
    type != "ORDERED_PROBIT"  &&
      type != "ORDERED_LOGIT"   &&
      type != "MULTINOMIAL_LOGIT" &&
      type != "MRI"             &&
      type != "MRI_AR1");
  if (!is_race) {
    return calc_ll(particle_matrix, data, constants, designs, type, bounds,
                   transforms, pretransforms, p_types, min_ll, trend,
                   return_trialwise);
  }

  // ---------------------------------------------------------------------------
  // Thread count
  // ---------------------------------------------------------------------------
#ifdef _OPENMP
  const int max_threads    = omp_get_max_threads();
  const int n_threads_used = (n_threads <= 0) ? max_threads
  : std::min(n_threads, max_threads);
#else
  if (n_threads > 1)
    Rcpp::warning("calc_ll_multithreaded: OpenMP not available, running single-threaded.");
  const int n_threads_used = 1;
#endif

  // ---------------------------------------------------------------------------
  // Shared setup (identical to calc_ll)
  // ---------------------------------------------------------------------------
  const int n_particles     = particle_matrix.nrow();
  const int n_rows        = data.nrow();
  const bool has_lR         = (sum(contains(data.names(), "lR")) == 1);
  const int n_lR            = has_lR ? unique(IntegerVector(data["lR"])).length() : 1;
  const int n_choice_trials = n_rows / n_lR;
  const int out_rows        = return_trialwise ? n_choice_trials : 1;

  NumericMatrix result(out_rows, n_particles);

  PipelineContext ctx = make_pipeline_context(particle_matrix, data, constants,
                                              designs, transforms, pretransforms, trend);
  TrendRuntime* trend_runtime_ptr = ctx.trend_runtime ? ctx.trend_runtime.get() : nullptr;

  NumericMatrix   minmax   = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<BoundSpec> bound_specs = make_bound_specs_pt(minmax, mm_names,
                                                           ctx.param_table, bounds);

  PipelineCache cache = make_pipeline_cache(ctx.param_table, designs,
                                            ctx.transform_specs, trend_runtime_ptr);

  // Pre-initialise design_plan to avoid lazy-init race inside map_from_designs
  // ctx.param_table.init_design_plan(designs);

  // Always needed - per-thread clones of paramtable and trend runtimes
  std::vector<ParamTable>                    pt_vec;
  std::vector<std::unique_ptr<TrendRuntime>> tr_vec;
  pt_vec.reserve(n_threads_used);
  tr_vec.reserve(n_threads_used);
  for (int t = 0; t < n_threads_used; ++t) {
    pt_vec.push_back(ctx.param_table.deep_copy());
    if (ctx.trend_runtime) {
      tr_vec.push_back(std::make_unique<TrendRuntime>(
          clone_trend_runtime(*ctx.trend_runtime, pt_vec.back())
      ));
    } else {
      tr_vec.push_back(nullptr);
    }
  }

  // Raw pointer to result data — safe to write column i from thread i
  double* result_ptr = result.begin();

  if(type == "DDM") {
    // DDM shared setup — extracted once, read-only in parallel loop
    NumericVector rts     = data["rt"];
    const double* rts_ptr = rts.begin();
    IntegerVector R_vec   = data["R"];
    const int*    Rs_ptr  = R_vec.begin();

    IntegerVector missingness;
    const bool has_missingness = (sum(contains(data.names(), "missingness")) == 1);
    if (has_missingness) missingness = data["missingness"];

    std::vector<int> idx_all;
    idx_all.reserve(n_rows);
    for (int i = 0; i < n_rows; ++i) {
      if (has_missingness && !IntegerVector::is_na(missingness[i])) continue;
      idx_all.push_back(i);
    }

    RaceModelSetup setup = make_race_setup(type, ctx.param_table);
    RaceScratch    scratch_tmp;
    scratch_tmp.reserve(n_rows);
    CensorSpec censor = make_censor_spec(data, n_choice_trials, n_lR, setup, ctx.param_table, scratch_tmp);
    TruncSpec  trunc  = make_trunc_spec (data, n_choice_trials, n_lR, setup, ctx.param_table, scratch_tmp);

    std::vector<std::vector<double>> ll_trial_vec(n_threads_used, std::vector<double>(n_choice_trials));
    std::vector<std::vector<int>>    is_ok_vec(n_threads_used,    std::vector<int>(n_rows, 1));
    std::vector<std::vector<double>> tw_vec(n_threads_used,
                                            std::vector<double>(return_trialwise ? n_choice_trials : 0));
    std::vector<RaceScratch>  scratch_vec(n_threads_used);  // DDM doesn't use scratch but rebind needs it
    std::vector<TruncSpec>  trunc_vec(n_threads_used, trunc);
    std::vector<CensorSpec> censor_vec(n_threads_used, censor);
    for (int t = 0; t < n_threads_used; ++t) {
      scratch_vec[t].reserve(n_rows);
      censor_vec[t].rebind(pt_vec[t], scratch_vec[t]);
      trunc_vec[t].rebind(pt_vec[t], scratch_vec[t]);
    }

    IntegerVector expand   = data.attr("expand");
    const int     n_exp    = expand.size();
    const int*    exp_ptr  = expand.begin();

#pragma omp parallel for schedule(static) num_threads(n_threads_used)
    for (int i = 0; i < n_particles; ++i) {
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      ParamTable&          pt_local = pt_vec[tid];
      TrendRuntime*        tr_local = tr_vec[tid].get();
      std::vector<double>& ll_trial = ll_trial_vec[tid];
      std::vector<int>&    is_ok    = is_ok_vec[tid];
      CensorSpec&    censor_local   = censor_vec[tid];
      TruncSpec&     trunc_local    = trunc_vec[tid];

      pt_local.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(pt_local, tr_local, cache);

      c_log_likelihood_DDM(rts_ptr, Rs_ptr, pt_local, setup.spec, idx_all, ll_trial.data());

      if(trunc_local.any()) trunc_local.calculate_normalization_constant();
      if(censor_local.any()) censor_local.fill_censored_rows(trunc_local, ll_trial, min_ll);
      if(trunc_local.any()) for (int t = 0; t < trunc_local.n_trials; ++t) ll_trial[t] -= trunc_local.log_Z[t];

      std::fill(is_ok.begin(), is_ok.end(), 1);
      c_do_bound_pt(pt_local, bound_specs, is_ok);
      for (int t = 0; t < n_choice_trials; ++t) if (!is_ok[t]) ll_trial[t] = min_ll;

      if (return_trialwise) {
        std::vector<double>& tw = tw_vec[tid];
        expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll, tw.data());
        double* col = result_ptr + (ptrdiff_t)i * out_rows;
        std::copy(tw.begin(), tw.end(), col);
      } else {
        const double sum = expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll);
        result_ptr[i] = sum;
      }
    }
  } else {
    // ---------------------------------------------------------------------------
    // Race-specific shared setup
    // ---------------------------------------------------------------------------
    // Extract raw pointers from Rcpp vectors BEFORE the parallel region —
    // no Rcpp types are touched inside the parallel loop.
    IntegerVector expand   = data.attr("expand");
    const int     n_exp    = expand.size();
    const int*    exp_ptr  = expand.begin();
    NumericVector rts    = data["rt"];
    const double* rts_ptr = rts.begin();
    LogicalVector winner  = data["winner"];
    const int*    win_flag = LOGICAL(winner);

    NumericVector lR     = data["lR"];

    // missingness handling
    IntegerVector missingness;
    const bool has_missingness = (sum(contains(data.names(), "missingness")) == 1);
    if (has_missingness) missingness = data["missingness"];

    std::vector<int> idx_win, idx_los;
    idx_win.reserve(n_rows);
    idx_los.reserve(n_rows);
    for (int i = 0; i < n_rows; ++i) {
      if(has_missingness && !IntegerVector::is_na(missingness[i])) continue;  // handled by censor
      if(win_flag[i]) {
        idx_win.push_back(i);
      } else {
        idx_los.push_back(i);
      }
    }

    RaceModelSetup setup = make_race_setup(type, ctx.param_table);
    RaceScratch scratch_tmp;
    scratch_tmp.reserve(n_rows);
    CensorSpec censor = make_censor_spec(data, n_choice_trials, n_lR, setup, ctx.param_table, scratch_tmp);
    TruncSpec  trunc  = make_trunc_spec (data, n_choice_trials, n_lR, setup, ctx.param_table, scratch_tmp);

    // Per-thread scratch — all plain std::vector, no Rcpp types
    const int scratch_size = n_rows;
    std::vector<RaceScratch>         scratch_vec(n_threads_used);
    std::vector<std::vector<double>> ll_row_vec(n_threads_used,   std::vector<double>(n_rows, 1.0));
    std::vector<std::vector<double>> ll_trial_vec(n_threads_used, std::vector<double>(n_choice_trials));
    std::vector<std::vector<int>>    is_ok_vec(n_threads_used,    std::vector<int>(n_rows, 1));
    // trialwise output buffer — written then copied to result, avoiding Rcpp inside loop
    std::vector<std::vector<double>> tw_vec(n_threads_used,
                                            std::vector<double>(return_trialwise ? n_choice_trials : 0));
    std::vector<TruncSpec> trunc_vec(n_threads_used, trunc);
    std::vector<CensorSpec> censor_vec(n_threads_used, censor);
    for (int t = 0; t < n_threads_used; ++t) {
      scratch_vec[t].reserve(scratch_size);
      censor_vec[t].rebind(pt_vec[t], scratch_vec[t]);
      trunc_vec[t].rebind(pt_vec[t], scratch_vec[t]);
    }

    // #pragma omp parallel for schedule(dynamic, 4) num_threads(n_threads_used)
#pragma omp parallel for schedule(static) num_threads(n_threads_used)
    for (int i = 0; i < n_particles; ++i) {
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif

      ParamTable&          pt_local  = pt_vec[tid];
      TrendRuntime*        tr_local  = tr_vec[tid].get();
      RaceScratch&         scratch   = scratch_vec[tid];
      std::vector<double>& ll_row    = ll_row_vec[tid];
      std::vector<double>& ll_trial  = ll_trial_vec[tid];
      std::vector<int>&    is_ok     = is_ok_vec[tid];
      CensorSpec& censor_local = censor_vec[tid];
      TruncSpec&  trunc_local  = trunc_vec[tid];

      pt_local.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(pt_local, tr_local, cache);

      std::fill(ll_row.begin(),   ll_row.end(),   1.0);
      std::fill(is_ok.begin(),    is_ok.end(),    1);

      c_log_likelihood_race(pt_local, setup, rts_ptr,
                            idx_win, idx_los, n_lR,
                            ll_row.data(), (int)ll_row.size(), ll_trial.data(), scratch);

      // truncation
      if (trunc_local.any()) trunc_local.calculate_normalization_constant();
      if (censor_local.any()) censor_local.fill_censored_rows(trunc_local, ll_trial, min_ll);
      if(trunc_local.any()) for (int t = 0; t < trunc_local.n_trials; ++t) ll_trial[t] -= trunc_local.log_Z[t];

      c_do_bound_pt(pt_local, bound_specs, is_ok);
      lr_all(is_ok, n_lR);                             // make sure the ok value is shared across accumulators / levels or lR
      for (int t = 0; t < (int)idx_win.size(); ++t) if(!is_ok[idx_win[t]]) ll_trial[idx_win[t] / n_lR] = min_ll;
      // for(int t=0; t < n_choice_trials; t++) {
      //   for(int k=0; k < n_lR; k++) {
      //     if(!is_ok[t * n_lR + k]) {
      //       // set to min_ll and move to next trial - any subsequent !is_ok checks are unnecessary
      //       ll_trial[t] = min_ll;
      //       break;  // no need to check remaining accumulators for this trial
      //     }
      //   }
      // }

      if (return_trialwise) {
        std::vector<double>& tw = tw_vec[tid];
        expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll, tw.data());
        // Write trialwise column i into result — raw pointer, no Rcpp
        double* col = result_ptr + (ptrdiff_t)i * out_rows;
        std::copy(tw.begin(), tw.end(), col);
      } else {
        const double sum = expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll);
        result_ptr[i] = sum;  // out_rows == 1, column i is just element i
      }
    }
  }

  return result;
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
  const int n_lR        = has_lR ? unique(IntegerVector(data["lR"])).length() : 1;

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
    run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);

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


// [[Rcpp::export]]
void omp_diagnostics(int n_threads = -1) {
#ifdef _OPENMP
  Rcpp::Rcout << "OpenMP is compiled in.\n";
  Rcpp::Rcout << "omp_get_max_threads() = " << omp_get_max_threads() << "\n";
  Rcpp::Rcout << "omp_get_num_procs()   = " << omp_get_num_procs()   << "\n";

  if (n_threads > 0) omp_set_num_threads(n_threads);

  std::vector<int> seen(omp_get_max_threads(), 0);

#pragma omp parallel for schedule(static) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
  for (int i = 0; i < 100; ++i) {
    seen[omp_get_thread_num()] = 1;
  }

  int active = 0;
  for (int t = 0; t < (int)seen.size(); ++t) {
    if (seen[t]) {
      Rcpp::Rcout << "  thread " << t << " was active\n";
      ++active;
    }
  }
  Rcpp::Rcout << "Total active threads: " << active << "\n";
#else
  Rcpp::Rcout << "OpenMP is NOT compiled in (_OPENMP not defined).\n";
#endif
}
