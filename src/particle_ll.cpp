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

#include "model_CDM.h"
// for extract_y -- should be moved elsewhere
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
  std::vector<TransformSpec> t_specs = make_transform_specs(particle_matrix, pretransforms);
  ctx.particle_matrix = c_do_transform(particle_matrix, t_specs);

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
  ctx.transform_specs = make_transform_specs(ctx.param_table, transforms);

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
  std::vector<TransformSpec> specs = make_transform_specs(pars, transform);
  // Apply transformation in place and return
  return c_do_transform(pars, specs);
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
      c_do_transform(param_table, cache.premap_specs);
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
      c_do_transform(param_table, cache.pretransform_specs);
    }
    for (BaseRuntime& base : trend_runtime->pretransform_bases) {
      trend_runtime->apply_base(base, param_table);
    }
  }

  // 4) Transforms for all parameters excluding trend pars used so far
  c_do_transform(param_table, cache.postmap_specs);

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

// ---------------------------------------------------------------------------
// c_expand_ordered_cut — raw double arrays, no Rcpp, safe in parallel
// cut_in / cut_out both length n_rows = n_trials * n_lR
// ---------------------------------------------------------------------------
inline void c_expand_ordered_cut(const double* __restrict__ cut_in,
                                 double* __restrict__       cut_out,
                                 int n_rows,
                                 int n_lR)
{
  const int n_trials = n_rows / n_lR;

  if (n_lR == 2) {
    for (int t = 0; t < n_trials; ++t) {
      const int base  = t * n_lR;
      cut_out[base]   = cut_in[base];
      cut_out[base+1] = cut_in[base];
    }
    return;
  }

  for (int t = 0; t < n_trials; ++t) {
    const int base = t * n_lR;
    double current = cut_in[base];
    cut_out[base]  = current;

    for (int r = 1; r < n_lR - 1; ++r) {
      current          += std::exp(cut_in[base + r]);
      cut_out[base + r] = current;
    }
    cut_out[base + n_lR - 1] = cut_out[base + n_lR - 2];
  }
}

double c_ordered_cdf(double x, double location, double scale, bool probit) {
  if (x == R_NegInf) return 0.0;
  if (x == R_PosInf) return 1.0;
  if (probit) return PNORM_STD((x - location) / scale, true, false);
  const double z = (x - location) / scale;
  return 1.0 / (1.0 + std::exp(-z));
}

// =============================================================================
// expand_clamp_sum — final step for all models with expand
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

// is_ok handling - shared by all models
inline void apply_bounds(std::vector<int>&             is_ok,
                         double*                       ll,
                         int                           n_trials,
                         int                           n_lR,
                         double                        min_ll,
                         const std::vector<bool>&      participating)   // length n_rows, all true if no missingness
{
  if (n_lR == 1) {
#pragma omp simd
    for (int t = 0; t < n_trials; ++t)
      ll[t] = (!is_ok[t] && participating[t]) ? min_ll : ll[t];

  } else {
    for (int t = 0; t < n_trials; ++t)
      for (int k = 0; k < n_lR; ++k) {
        const int row = t * n_lR + k;
        if (!participating[row]) continue;
        if (!is_ok[row]) { ll[t] = min_ll; break; }
      }
  }
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
                           double* ll_buf_ptr,
                           RaceScratch& scratch)
{
  const int n_winners = (int)idx_win.size();

  setup.fill_both(rts_ptr, pt, setup.spec,
                  idx_win, idx_los, ll_row_ptr, scratch);
  vec_log(ll_row_ptr, ll_row_size);

  if (n_acc == 1) {
    for (int t = 0; t < n_winners; ++t)
      ll_buf_ptr[idx_win[t]] = ll_row_ptr[idx_win[t]];
  } else {
    for (int t = 0; t < n_winners; ++t) {
      // NB: ll_buf_ptr is sized n_trials - including the missing (NA) trials. Can't index with t directly...
      const int trial_idx = idx_win[t] / n_acc;
      const int base      = trial_idx * n_acc;
      double ll = 0.0;
      for (int k = 0; k < n_acc; ++k)
        ll += ll_row_ptr[base + k];
      ll_buf_ptr[trial_idx] = ll;
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


inline void c_log_likelihood_ordered(const ParamTable&       pt,
                                     const ChoiceOnlySpec&   spec,
                                     const std::vector<int>& idx_win,
                                     const int*              lR_ptr,
                                     int                     n_rows,
                                     int                     n_lR,
                                     double* __restrict__    cut_buf,
                                     bool                    probit,
                                     double* __restrict__    ll_buf)
{
  const double* loc_col = pt.base.colptr(spec.col_location);
  const double* sc_col  = pt.base.colptr(spec.col_scale);
  const double* cut_col = pt.base.colptr(spec.col_cut);

  c_expand_ordered_cut(cut_col, cut_buf, n_rows, n_lR);

  const int n_winners = (int)idx_win.size();
  for (int t = 0; t < n_winners; ++t) {
    const int    i        = idx_win[t];
    const int    level    = lR_ptr[i];
    const double location = loc_col[i];
    const double scale    = sc_col[i];
    const double upper    = (level == n_lR) ? R_PosInf : cut_buf[i];
    const double lower    = (level == 1)    ? R_NegInf : cut_buf[i - 1];
    const double prob     = c_ordered_cdf(upper, location, scale, probit) - c_ordered_cdf(lower, location, scale, probit);
    ll_buf[i/n_lR] = std::log(prob);
  }
  // vec_log(ll_buf, n_rows);
}

inline void c_log_likelihood_multinomial_logit(const ParamTable&       pt,
                                               const ChoiceOnlySpec&   spec,
                                               const std::vector<int>& idx_win,
                                               int                     n_lR,
                                               double* __restrict__    ll_buf)
{
  const double* util_col    = pt.base.colptr(spec.col_utility);
  const int     n_winners   = (int)idx_win.size();

  for (int t = 0; t < n_winners; ++t) {
    const int base = (idx_win[t] / n_lR) * n_lR;  // trial base row

    double max_u = util_col[base];
    for (int r = 1; r < n_lR; ++r)
      if (util_col[base + r] > max_u) max_u = util_col[base + r];

    double denom = 0.0;
    for (int r = 0; r < n_lR; ++r)
      denom += std::exp(util_col[base + r] - max_u);

    const double chosen = std::exp(util_col[idx_win[t]] - max_u);
    ll_buf[idx_win[t] / n_lR] = std::log(chosen / denom);
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
    ll_buf[i] = DNORM(y[i], s, sigma[i]);
  }
  vec_log(ll_buf, n);
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
  ll_buf[0] = DNORM(y[0], s_prev, sigma[0]);

  // t >= 1: AR(1) conditional
  for (int i = 1; i < n; i++) {
    double s_curr = 0.0;
    for (int j = 0; j < nm; j++) s_curr += mean_cols[j][i];

    const double rho_i    = rho[i];
    const double cond_sd  = sigma[i] * std::sqrt(1.0 - rho_i * rho_i);
    const double cond_mean = s_curr + rho_i * (y[i - 1] - s_prev);

    ll_buf[i] = DNORM(y[i], cond_mean, cond_sd);
    s_prev = s_curr;
  }
  vec_log(ll_buf, n);
}

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


  // n_exp needed for out_rows — extract expand early for non-MRI models
  const bool is_mri = (type == "MRI" || type == "MRI_AR1");
  IntegerVector expand;
  const int* exp_ptr = nullptr;
  int n_exp = n_choice_trials;  // MRI has no expand, trialwise = n_choice_trials
  if (!is_mri) {
    expand  = data.attr("expand");
    n_exp   = expand.size();
    exp_ptr = expand.begin();
  }

  const int out_rows = return_trialwise ? n_exp : 1;
  NumericMatrix result(out_rows, n_particles);
  double* result_ptr = result.begin();
  // Column i = particle i's trialwise buffer (contiguous)
  // NumericMatrix result(out_rows, n_particles);  // column i = particle i, contiguous

  std::vector<int>     is_ok(n_rows, 1);
  std::vector<double>  ll_buf(n_choice_trials);  // compressed scratch, reused per particle


  // Shared setup -- context holds the param_table as well as designs, constants, trend etc
  PipelineContext ctx = make_pipeline_context(particle_matrix, data, constants,
                                              designs, transforms, pretransforms, trend);
  TrendRuntime* trend_runtime_ptr = ctx.trend_runtime ? ctx.trend_runtime.get() : nullptr;
  PipelineCache cache = make_pipeline_cache(ctx.param_table, designs,
                                            ctx.transform_specs, trend_runtime_ptr);

  // Bounds — built once from structure
  NumericMatrix   minmax   = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<BoundSpec> bound_specs = make_bound_specs(minmax, mm_names, ctx.param_table, bounds);
  std::vector<bool> participating(n_rows, true);

  // -----------------------------------------------------------------------
  // MRI / MRI_AR1 -- no expand, rt, R attributes
  // -----------------------------------------------------------------------
  if (type == "MRI" || type == "MRI_AR1") {
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
      else        c_log_likelihood_MRI_white(ctx.param_table, spec, y, n_choice_trials, ll_buf.data());

      c_do_bound(ctx.param_table, bound_specs, is_ok);
      apply_bounds(is_ok, ll_buf.data(), n_choice_trials, /* n_lR = */ 1, min_ll, participating);

      double* tw = return_trialwise ? result_ptr + (ptrdiff_t)i * out_rows : nullptr;
      const double sum = clamp_sum(ll_buf.data(), n_choice_trials, min_ll, tw);
      if (!return_trialwise) result(0, i) = sum;
    }
  } else {
    // Missingness — shared across all choice models
    IntegerVector missingness;
    const bool has_missingness = (sum(contains(data.names(), "missingness")) == 1);
    if (has_missingness) missingness = data["missingness"];

    // Winner flag — shared across all choice models except DDM
    // (choice-only models don't use it yet but will)
    const bool has_winner = (sum(contains(data.names(), "winner")) == 1);
    LogicalVector winner;
    int* win_flag = nullptr;
    if (has_winner) {
      winner   = data["winner"];
      win_flag = LOGICAL(winner);
    }

    // Index lists — built once, shared across all choice models
    // idx_all: all non-missing rows (DDM, choice-only)
    // idx_win / idx_los: winner/loser split (race, and future choice-only)
    std::vector<int> idx_all, idx_win, idx_los;
    idx_all.reserve(n_rows);
    if (has_winner) {
      idx_win.reserve(n_rows);
      idx_los.reserve(n_rows);
    }
    for (int i = 0; i < n_rows; ++i) {
      if (has_missingness && !IntegerVector::is_na(missingness[i])) {
        if (missingness[i] == 0) {
          participating[i] = false;  // ignored: bounds checking skipped entirely
        }
        continue;                 // but skip from likelihood - probability calculated by censorspec
      }
      idx_all.push_back(i);
      if (has_winner) {
        if (win_flag[i]) idx_win.push_back(i);
        else             idx_los.push_back(i);
      }
    }

    // -----------------------------------------------------------------------
    // Choice-only models (ORDERED_PROBIT, ORDERED_LOGIT, MULTINOMIAL_LOGIT)
    // -----------------------------------------------------------------------
    if (type == "ORDERED_PROBIT" || type == "ORDERED_LOGIT" || type == "MULTINOMIAL_LOGIT") {
      // Shared choice-only setup
      const bool is_probit = (type == "ORDERED_PROBIT");  // unused for MULTINOMIAL_LOGIT, but fine
      Rcpp::IntegerVector lR = data["lR"];
      const int* lR_ptr = INTEGER(lR);
      std::vector<double> cut_buf(n_rows);
      // ll_buf sized n_choice_trials, filled at trial index
      std::vector<double> ll_buf(n_choice_trials, 0.0);  // default 1.0 = log(1) = 0 for missing
      // to-do -  make log_likelihoods fill out probabilities first and then do vec_log over the buffer.
      // In that case
      ChoiceOnlySpec spec = make_choice_only_spec(ctx.param_table, type);

      for (int i = 0; i < n_particles; ++i) {
        // 1) Map p_vector to trialwise parameters
        if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
        run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);

        // 2) calculate raw (compressed) trialwise log-likelihoods (fills ll_buf)
        if (type == "MULTINOMIAL_LOGIT") {
          c_log_likelihood_multinomial_logit(ctx.param_table, spec, idx_win, n_lR, ll_buf.data());
        } else {
          c_log_likelihood_ordered(ctx.param_table, spec, idx_win, lR_ptr, n_rows, n_lR, cut_buf.data(), is_probit, ll_buf.data());
        }

        // 3) Handle not-ok parameter values (out of bound)
        c_do_bound(ctx.param_table, bound_specs, is_ok);
        apply_bounds(is_ok, ll_buf.data(), n_choice_trials, n_lR, min_ll, participating);

        // 4) Determine output location (tw is a pointer to the correct address in result) and protect via expand, clamp, sum
        double* tw = return_trialwise ? result_ptr + (ptrdiff_t)i * out_rows : nullptr;
        const double sum = expand_clamp_sum(ll_buf.data(), exp_ptr, n_exp, min_ll, tw);
        if (!return_trialwise) result(0, i) = sum;
      }
      // -----------------------------------------------------------------------
      // Continuous-choice-RT models (CDM, PSDM, PHSDM)
      // -----------------------------------------------------------------------
    } else if (type == "CDM" || type == "PSDM" || type == "PHSDM") {
      const bool has_R2 = (sum(contains(data.names(), "R2")) == 1);
      const bool has_R3 = (sum(contains(data.names(), "R3")) == 1);
      NumericVector rts = data["rt"];
      NumericVector Rs  = data["R"];
      NumericVector R2s = has_R2 ? NumericVector(data["R2"]) : NumericVector();
      NumericVector R3s = has_R3 ? NumericVector(data["R3"]) : NumericVector();

      std::vector<double> ll_trial(n_choice_trials, 0.0);     // compressed scratch for (log)likelihoods in race (compressed! so needs expanding)

      for (int i = 0; i < n_particles; ++i) {
        std::fill(ll_trial.begin(), ll_trial.end(), 0.0);
        std::fill(is_ok.begin(),  is_ok.end(),  1);
        if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
        run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);
        NumericMatrix pars = get_pars_matrix(ctx.param_table, ctx.keep_names);

        if (type == "CDM") {
          // Sub-type dispatch: CDM -> HSDM (R,R2,R3) -> SDM (R,R2) -> CDM (R only)
          if(has_R2 && has_R3 && pars.ncol() >= 8)
            c_dHSDM(rts, Rs, R2s, R3s, pars, participating, ll_trial.data(), n_choice_trials);
          else if (has_R2 && pars.ncol() >= 7)
            c_dSDM (rts, Rs, R2s, pars, participating, ll_trial.data(), n_choice_trials);
          else
            c_dCDM (rts, Rs, pars, participating, ll_trial.data(), n_choice_trials);
        } else if (type == "PSDM") {
          c_dPSDM (rts, Rs, pars, participating, ll_trial.data(), n_choice_trials);
        } else {
          c_dPHSDM(rts, Rs, R2s,  pars, participating, ll_trial.data(), n_choice_trials);
        }
        c_do_bound(ctx.param_table, bound_specs, is_ok);
        apply_bounds(is_ok, ll_trial.data(), n_choice_trials, 1, min_ll, participating);

        double* tw = return_trialwise ? result_ptr + (ptrdiff_t)i * out_rows : nullptr;
        const double sum = expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll, tw);
        if (!return_trialwise) result(0, i) = sum;
      }
    // -----------------------------------------------------------------------
    // Discrete-Choice-RT models (DDM, & Race: RDM, LBA, LNR, ...)
    // -----------------------------------------------------------------------
    } else {
      // Shared Choice-RT setup
      NumericVector rts     = data["rt"];
      const double* rts_ptr = rts.begin();

      RaceModelSetup setup = make_race_setup(type, ctx.param_table);
      RaceScratch scratch;
      auto [n_dbl, n_int] = race_scratch_slots(type);
      scratch.reserve(n_dbl, n_int, n_rows);

      CensorSpec censor = make_censor_spec(data, n_choice_trials, n_lR, participating, &setup, ctx.param_table, scratch);  // n_lR equals 1 for DDM
      TruncSpec  trunc  = make_trunc_spec (data, n_choice_trials, n_lR, participating, &setup, ctx.param_table, scratch);

      if (type == "DDM") {
        IntegerVector R      = data["R"];
        const int*    Rs_ptr = R.begin();

        // Loop over particles
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

          // 4) Check for bound violations
          c_do_bound(ctx.param_table, bound_specs, is_ok);
          apply_bounds(is_ok, ll_buf.data(), n_choice_trials, /* n_lR = */ 1, min_ll, participating);

          // 5) Determine output location (tw is a pointer to the correct address in result) and protect via expand, clamp, sum
          double* tw = return_trialwise ? result_ptr + (ptrdiff_t)i * out_rows : nullptr;
          const double sum = expand_clamp_sum(ll_buf.data(), exp_ptr, n_exp, min_ll, tw);
          if (!return_trialwise) result(0, i) = sum;
        }

        // ---- Race models (RDM, LBA, LNR, ...) ----
      } else {
        // Scratch buffers (reused across particles)
        std::vector<double> ll_row(n_rows);                 // stores (log)likelihood of row in dadm
        // ll_trial sized n_choice_trials; assumes n_winners <= n_choice_trials (one winner per trial at most)
        std::vector<double> ll_trial(n_choice_trials);     // compressed scratch for (log)likelihoods in race (compressed! so needs expanding)

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
          c_do_bound(ctx.param_table, bound_specs, is_ok);
          apply_bounds(is_ok, ll_trial.data(), n_choice_trials, n_lR, min_ll, participating);

          // 5) Expand, clamp, sum, etc
          double* tw = return_trialwise ? result_ptr + (ptrdiff_t)i * out_rows : nullptr;
          const double sum = expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll, tw);
          if (!return_trialwise) result(0, i) = sum;
        }
      }
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
  // Shared setup
  // ---------------------------------------------------------------------------
  const int  n_particles     = particle_matrix.nrow();
  const int  n_rows          = data.nrow();
  const bool has_lR          = (sum(contains(data.names(), "lR")) == 1);
  const int  n_lR            = has_lR ? unique(IntegerVector(data["lR"])).length() : 1;
  const int  n_choice_trials = n_rows / n_lR;

  // n_exp needed for out_rows — extract expand early for non-MRI models
  const bool is_mri = (type == "MRI" || type == "MRI_AR1");
  IntegerVector expand;
  const int* exp_ptr = nullptr;
  int n_exp = n_choice_trials;  // MRI has no expand, trialwise = n_choice_trials
  if (!is_mri) {
    expand  = data.attr("expand");
    n_exp   = expand.size();
    exp_ptr = expand.begin();
  }

  const int out_rows = return_trialwise ? n_exp : 1;
  NumericMatrix result(out_rows, n_particles);
  double* result_ptr = result.begin();

  PipelineContext ctx = make_pipeline_context(particle_matrix, data, constants,
                                              designs, transforms, pretransforms, trend);
  TrendRuntime* trend_runtime_ptr = ctx.trend_runtime ? ctx.trend_runtime.get() : nullptr;

  NumericMatrix   minmax   = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<BoundSpec> bound_specs = make_bound_specs(minmax, mm_names,
                                                        ctx.param_table, bounds);
  PipelineCache cache = make_pipeline_cache(ctx.param_table, designs,
                                            ctx.transform_specs, trend_runtime_ptr);

  // ---------------------------------------------------------------------------
  // Per-thread mutable state — always needed
  // ---------------------------------------------------------------------------
  std::vector<ParamTable>                    pt_vec;
  std::vector<std::unique_ptr<TrendRuntime>> tr_vec;
  pt_vec.reserve(n_threads_used);
  tr_vec.reserve(n_threads_used);
  for (int t = 0; t < n_threads_used; ++t)
    pt_vec.push_back(ctx.param_table.deep_copy());
  for (int t = 0; t < n_threads_used; ++t) {
    if (ctx.trend_runtime)
      tr_vec.push_back(std::make_unique<TrendRuntime>(
          clone_trend_runtime(*ctx.trend_runtime, pt_vec[t])));
    else
      tr_vec.push_back(nullptr);
  }

  std::vector<std::vector<double>> ll_buf_vec(n_threads_used, std::vector<double>(n_choice_trials, 1.0));
  std::vector<std::vector<int>>    is_ok_vec (n_threads_used, std::vector<int>   (n_rows,          1));
  // tw_vec sized n_exp — written by expand_clamp_sum then copied into result
  std::vector<std::vector<double>> tw_vec    (n_threads_used,
                                              std::vector<double>(return_trialwise ? n_exp : 0));

  // ---------------------------------------------------------------------------
  // MRI / MRI_AR1
  // ---------------------------------------------------------------------------
  if (is_mri) {
    NumericVector y_rcpp = extract_y(data);
    const double* y      = y_rcpp.begin();
    const bool    is_ar1 = (type == "MRI_AR1");

    MRISpec spec = make_mri_spec(ctx.param_table, ctx.keep_names, is_ar1);
    const std::vector<bool> participating(n_rows, true);

#pragma omp parallel for schedule(static) num_threads(n_threads_used)
    for (int i = 0; i < n_particles; ++i) {
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      ParamTable&          pt_local = pt_vec[tid];
      TrendRuntime*        tr_local = tr_vec[tid].get();
      std::vector<double>& ll_buf   = ll_buf_vec[tid];
      std::vector<int>&    is_ok    = is_ok_vec[tid];

      pt_local.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
      run_pars_pipeline(pt_local, tr_local, cache);

      std::fill(is_ok.begin(), is_ok.end(), 1);
      if (is_ar1) c_log_likelihood_MRI_ar1  (pt_local, spec, y, n_choice_trials, ll_buf.data());
      else        c_log_likelihood_MRI_white(pt_local, spec, y, n_choice_trials, ll_buf.data());

      c_do_bound(pt_local, bound_specs, is_ok);
      apply_bounds(is_ok, ll_buf.data(), n_choice_trials, 1, min_ll, participating);

      if (return_trialwise) {
        std::vector<double>& tw = tw_vec[tid];
        clamp_sum(ll_buf.data(), n_choice_trials, min_ll, tw.data());
        std::copy(tw.begin(), tw.end(), result_ptr + (ptrdiff_t)i * out_rows);
      } else {
        result_ptr[i] = clamp_sum(ll_buf.data(), n_choice_trials, min_ll);
      }
    }

  } else {
    // ---------------------------------------------------------------------------
    // All choice models — shared read-only data
    // ---------------------------------------------------------------------------
    IntegerVector missingness;
    const bool has_missingness = (sum(contains(data.names(), "missingness")) == 1);
    if (has_missingness) missingness = data["missingness"];

    const bool has_winner = (sum(contains(data.names(), "winner")) == 1);
    LogicalVector winner;
    int* win_flag = nullptr;
    if (has_winner) {
      winner   = data["winner"];
      win_flag = LOGICAL(winner);
    }

    std::vector<bool> participating(n_rows, true);

    std::vector<int> idx_all, idx_win, idx_los;
    idx_all.reserve(n_rows);
    if (has_winner) { idx_win.reserve(n_rows); idx_los.reserve(n_rows); }
    for (int i = 0; i < n_rows; ++i) {
      if (has_missingness && !IntegerVector::is_na(missingness[i])) {
        if (missingness[i] == 0) participating[i] = false;
        continue;
      }
      idx_all.push_back(i);
      if (has_winner) {
        if (win_flag[i]) idx_win.push_back(i);
        else             idx_los.push_back(i);
      }
    }

    // -----------------------------------------------------------------------
    // Choice-only branch
    // -----------------------------------------------------------------------
    const bool is_choice_only = (type == "ORDERED_PROBIT" ||
                                 type == "ORDERED_LOGIT"  ||
                                 type == "MULTINOMIAL_LOGIT");
    if (is_choice_only) {
      const bool     is_probit = (type == "ORDERED_PROBIT");
      ChoiceOnlySpec spec      = make_choice_only_spec(ctx.param_table, std::string(type));

      IntegerVector lR_vec = has_lR ? IntegerVector(data["lR"]) : IntegerVector();
      const int* lR_ptr    = has_lR ? INTEGER(lR_vec) : nullptr;

      std::vector<std::vector<double>> cut_buf_vec(n_threads_used, std::vector<double>(n_rows, 0.0));
      for (auto& v : ll_buf_vec) std::fill(v.begin(), v.end(), 0.0);

#pragma omp parallel for schedule(static) num_threads(n_threads_used)
      for (int i = 0; i < n_particles; ++i) {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        ParamTable&          pt_local = pt_vec[tid];
        TrendRuntime*        tr_local = tr_vec[tid].get();
        std::vector<double>& ll_buf   = ll_buf_vec[tid];
        std::vector<int>&    is_ok    = is_ok_vec[tid];
        std::vector<double>& cut_buf  = cut_buf_vec[tid];

        pt_local.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
        run_pars_pipeline(pt_local, tr_local, cache);

        std::fill(ll_buf.begin(), ll_buf.end(), 0.0);
        std::fill(is_ok.begin(),  is_ok.end(),  1);

        if (type == "MULTINOMIAL_LOGIT") {
          c_log_likelihood_multinomial_logit(pt_local, spec, idx_win, n_lR, ll_buf.data());
        } else {
          c_log_likelihood_ordered(pt_local, spec, idx_win, lR_ptr,
                                   n_rows, n_lR, cut_buf.data(), is_probit, ll_buf.data());
        }

        c_do_bound(pt_local, bound_specs, is_ok);
        apply_bounds(is_ok, ll_buf.data(), n_choice_trials, n_lR, min_ll, participating);

        if (return_trialwise) {
          std::vector<double>& tw = tw_vec[tid];
          expand_clamp_sum(ll_buf.data(), exp_ptr, n_exp, min_ll, tw.data());
          std::copy(tw.begin(), tw.end(), result_ptr + (ptrdiff_t)i * out_rows);
        } else {
          result_ptr[i] = expand_clamp_sum(ll_buf.data(), exp_ptr, n_exp, min_ll);
        }
      }

      // -----------------------------------------------------------------------
      // CDM / PSDM / PHSDM branch
      // -----------------------------------------------------------------------
    } else if (type == "CDM" || type == "PSDM" || type == "PHSDM") {
      const bool has_R2 = (sum(contains(data.names(), "R2")) == 1);
      const bool has_R3 = (sum(contains(data.names(), "R3")) == 1);
      NumericVector rts = data["rt"];
      NumericVector Rs  = data["R"];
      NumericVector R2s = has_R2 ? NumericVector(data["R2"]) : NumericVector();
      NumericVector R3s = has_R3 ? NumericVector(data["R3"]) : NumericVector();

      std::vector<std::vector<double>> ll_trial_vec(n_threads_used, std::vector<double>(n_choice_trials, 0.0));

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

        pt_local.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
        run_pars_pipeline(pt_local, tr_local, cache);

        std::fill(ll_trial.begin(), ll_trial.end(), 0.0);
        std::fill(is_ok.begin(),    is_ok.end(),    1);

        NumericMatrix pars = get_pars_matrix(pt_local, ctx.keep_names);

        if (type == "CDM") {
          if      (has_R2 && has_R3 && pars.ncol() >= 8)
            c_dHSDM(rts, Rs, R2s, R3s, pars, participating, ll_trial.data(), n_choice_trials);
          else if (has_R2 && pars.ncol() >= 7)
            c_dSDM (rts, Rs, R2s,       pars, participating, ll_trial.data(), n_choice_trials);
          else
            c_dCDM (rts, Rs,            pars, participating, ll_trial.data(), n_choice_trials);
        } else if (type == "PSDM") {
          c_dPSDM (rts, Rs,      pars, participating, ll_trial.data(), n_choice_trials);
        } else {
          c_dPHSDM(rts, Rs, R2s, pars, participating, ll_trial.data(), n_choice_trials);
        }

        c_do_bound(pt_local, bound_specs, is_ok);
        apply_bounds(is_ok, ll_trial.data(), n_choice_trials, 1, min_ll, participating);

        if (return_trialwise) {
          std::vector<double>& tw = tw_vec[tid];
          expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll, tw.data());
          std::copy(tw.begin(), tw.end(), result_ptr + (ptrdiff_t)i * out_rows);
        } else {
          result_ptr[i] = expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll);
        }
      }

      // -----------------------------------------------------------------------
      // Choice-RT branch (DDM + race)
      // -----------------------------------------------------------------------
    } else {
      NumericVector rts     = data["rt"];
      const double* rts_ptr = rts.begin();

      IntegerVector R_vec;
      const int*    Rs_ptr = nullptr;
      if (type == "DDM") {
        R_vec  = data["R"];
        Rs_ptr = R_vec.begin();
      }

      RaceModelSetup setup = make_race_setup(type, ctx.param_table);
      auto [n_dbl, n_int]  = race_scratch_slots(type);

      RaceScratch scratch_tmp;
      scratch_tmp.reserve(n_dbl, n_int, n_rows);
      CensorSpec censor = make_censor_spec(data, n_choice_trials, n_lR, participating,
                                           &setup, ctx.param_table, scratch_tmp);
      TruncSpec  trunc  = make_trunc_spec (data, n_choice_trials, n_lR, participating,
                                           &setup, ctx.param_table, scratch_tmp);

      std::vector<RaceScratch>         scratch_vec (n_threads_used);
      std::vector<std::vector<double>> ll_row_vec  (n_threads_used, std::vector<double>(n_rows,          1.0));
      std::vector<std::vector<double>> ll_trial_vec(n_threads_used, std::vector<double>(n_choice_trials));
      std::vector<CensorSpec>          censor_vec  (n_threads_used, censor);
      std::vector<TruncSpec>           trunc_vec   (n_threads_used, trunc);

      for (int t = 0; t < n_threads_used; ++t) {
        scratch_vec[t].reserve(n_dbl, n_int, n_rows);
        censor_vec[t].rebind(pt_vec[t], scratch_vec[t]);
        trunc_vec[t].rebind(pt_vec[t], scratch_vec[t]);
      }

#pragma omp parallel for schedule(static) num_threads(n_threads_used)
      for (int i = 0; i < n_particles; ++i) {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        ParamTable&          pt_local     = pt_vec[tid];
        TrendRuntime*        tr_local     = tr_vec[tid].get();
        RaceScratch&         scratch      = scratch_vec[tid];
        std::vector<double>& ll_buf       = ll_buf_vec[tid];
        std::vector<double>& ll_row       = ll_row_vec[tid];
        std::vector<double>& ll_trial     = ll_trial_vec[tid];
        std::vector<int>&    is_ok        = is_ok_vec[tid];
        CensorSpec&          censor_local = censor_vec[tid];
        TruncSpec&           trunc_local  = trunc_vec[tid];

        pt_local.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
        run_pars_pipeline(pt_local, tr_local, cache);

        if (type == "DDM") {
          std::fill(ll_buf.begin(), ll_buf.end(), 1.0);
          std::fill(is_ok.begin(),  is_ok.end(),  1);
          c_log_likelihood_DDM(rts_ptr, Rs_ptr, pt_local, setup.spec, idx_all, ll_buf.data());

          if (trunc_local.any())  trunc_local.calculate_normalization_constant();
          if (censor_local.any()) censor_local.fill_censored_rows(trunc_local, ll_buf, min_ll);
          if (trunc_local.any())  for (int t = 0; t < trunc_local.n_trials; ++t) ll_buf[t] -= trunc_local.log_Z[t];

          c_do_bound(pt_local, bound_specs, is_ok);
          apply_bounds(is_ok, ll_buf.data(), n_choice_trials, 1, min_ll, participating);

          if (return_trialwise) {
            std::vector<double>& tw = tw_vec[tid];
            expand_clamp_sum(ll_buf.data(), exp_ptr, n_exp, min_ll, tw.data());
            std::copy(tw.begin(), tw.end(), result_ptr + (ptrdiff_t)i * out_rows);
          } else {
            result_ptr[i] = expand_clamp_sum(ll_buf.data(), exp_ptr, n_exp, min_ll);
          }

        } else {
          std::fill(ll_row.begin(),  ll_row.end(),  1.0);
          std::fill(is_ok.begin(),   is_ok.end(),   1);
          c_log_likelihood_race(pt_local, setup, rts_ptr,
                                idx_win, idx_los, n_lR,
                                ll_row.data(), (int)ll_row.size(), ll_trial.data(), scratch);

          if (trunc_local.any())  trunc_local.calculate_normalization_constant();
          if (censor_local.any()) censor_local.fill_censored_rows(trunc_local, ll_trial, min_ll);
          if (trunc_local.any())  for (int t = 0; t < trunc_local.n_trials; ++t) ll_trial[t] -= trunc_local.log_Z[t];

          c_do_bound(pt_local, bound_specs, is_ok);
          apply_bounds(is_ok, ll_trial.data(), n_choice_trials, n_lR, min_ll, participating);

          if (return_trialwise) {
            std::vector<double>& tw = tw_vec[tid];
            expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll, tw.data());
            std::copy(tw.begin(), tw.end(), result_ptr + (ptrdiff_t)i * out_rows);
          } else {
            result_ptr[i] = expand_clamp_sum(ll_trial.data(), exp_ptr, n_exp, min_ll);
          }
        }
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
  const int n_lR        = has_lR ? unique(Rcpp::as<IntegerVector>(data["lR"])).length() : 1;

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
  std::vector<BoundSpec> bound_specs = make_bound_specs(minmax, mm_names, ctx.param_table, bounds);
  std::vector<int> is_ok(n_trials, 1);

  List result(n_particles);

  for (int i = 0; i < n_particles; ++i) {
    // Update param_table for this particle (skip refill on first particle)
    if (i > 0) ctx.param_table.fill_from_particle_row(ctx.particle_matrix, i, ctx.pm_col_to_base_idx);
    run_pars_pipeline(ctx.param_table, trend_runtime_ptr, cache);

    std::fill(is_ok.begin(), is_ok.end(), 1);
    c_do_bound(ctx.param_table, bound_specs, is_ok);
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
