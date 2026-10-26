#ifdef _OPENMP
#include <omp.h>
#endif

#include <Rcpp.h>
#include <unordered_map>
#include "LikelihoodContext.h"
#include "openmp.h"

// Utilities first — no dependencies on model types
#include "utility_functions.h"
#include "transform_utils.h"
#include "ParamTable.h"
#include "TrendEngine.h"
#include "math_utils.h"

#include "model_CDM.h"
// for extract_y -- should be moved elsewhere
#include "model_MRI.h"
#include "model_multinomial.h"

// RaceSetup last — references functions defined in model headers above
#include "RaceSetup.h"
#include "CensorSpec.h"
#include "TruncSpec.h"

// Stop-signal models (after RaceSetup.h: they build on model_RDM.h and
// model_exgaussian.h). Header-only; include from this translation unit only.
#include "model_SS_EXG.h"
#include "model_SS_RDEX.h"
#include "ss_fast.h"         // stop-signal: data-only SSSpec + thread-safe per-particle likelihood
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

// =============================================================================
// Likelihood functions — called by LikelihoodContext::Impl::evaluate_particle
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
                           RaceScratch& scratch,
                           const std::vector<int>& is_ok, bool log_p)
{
  const int n_winners = (int)idx_win.size();

  if (log_p) std::fill(ll_row_ptr,ll_row_ptr+ll_row_size,0.0);
  setup.fill_both(rts_ptr, pt, setup.spec,
                  idx_win, idx_los, ll_row_ptr, scratch, is_ok, log_p);
  if (!log_p) vec_log(ll_row_ptr, ll_row_size);

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
                          double* __restrict__ ll_row,
                          const std::vector<int>& is_ok)
{
  // this is perhaps a bit pointless to have? only for naming convention...
  fill_ddm(rts, R, pt,
           spec,
           idx_all,
           ll_row, is_ok);
}


inline void c_log_likelihood_ordered(const ParamTable&       pt,
                                     const ChoiceOnlySpec&   spec,
                                     const std::vector<int>& idx_win,
                                     const int*              lR_ptr,
                                     int                     n_rows,
                                     int                     n_lR,
                                     double* __restrict__    cut_buf,
                                     bool                    probit,
                                     double* __restrict__    ll_buf,
                                     const std::vector<int>& is_ok)
{
  const double* loc_col = pt.base.colptr(spec.col_location);
  const double* sc_col  = pt.base.colptr(spec.col_scale);
  const double* cut_col = pt.base.colptr(spec.col_cut);

  c_expand_ordered_cut(cut_col, cut_buf, n_rows, n_lR);

  const int n_winners = (int)idx_win.size();
  for (int t = 0; t < n_winners; ++t) {
    const int    i        = idx_win[t];
    if (!is_ok[i]) continue;
    const int    level    = lR_ptr[i];
    const double location = loc_col[i];
    const double scale    = sc_col[i];
    const double upper    = (level == n_lR) ? R_PosInf : cut_buf[i];
    const double lower    = (level == 1)    ? R_NegInf : cut_buf[i - 1];
    if (probit) {
      ll_buf[i/n_lR]=pnorm_detail::normal_interval((lower-location)/scale,(upper-location)/scale,true);
    } else {
      const double prob=c_ordered_cdf(upper,location,scale,false)-c_ordered_cdf(lower,location,scale,false);
      ll_buf[i/n_lR]=std::log(prob);
    }
  }
  // vec_log(ll_buf, n_rows);
}

inline void c_log_likelihood_multinomial_logit(const ParamTable&       pt,
                                               const ChoiceOnlySpec&   spec,
                                               const std::vector<int>& idx_win,
                                               int                     n_lR,
                                               double* __restrict__    ll_buf,
                                               const std::vector<int>& is_ok)
{
  const double* util_col    = pt.base.colptr(spec.col_utility);
  const int     n_winners   = (int)idx_win.size();

  for (int t = 0; t < n_winners; ++t) {
    const int base = (idx_win[t] / n_lR) * n_lR;  // trial base row
    if (!is_ok[base]) continue;

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

inline void c_log_likelihood_multinomial_probit(const ParamTable& pt,
                                               const ChoiceOnlySpec& spec,
                                               const std::vector<int>& idx_win,
                                               int n_lR,
                                               double* ll_buf,
                                               const std::vector<int>& is_ok)
{
  const double* utility = pt.base.colptr(spec.col_utility);
  for (int winner : idx_win) {
    if (!is_ok[winner]) continue;
    const int trial = winner / n_lR;
    const int base = trial * n_lR;
    ll_buf[trial] = std::log(multinomial_probit_probability(utility + base, n_lR, winner - base));
  }
}

void c_log_likelihood_MRI_white(const ParamTable& pt,
                                const MRISpec& spec,
                                const std::vector<const double*>& mean_cols,
                                const double* __restrict__ y,
                                int n,
                                double* __restrict__ ll_buf,
                                const std::vector<int>& is_ok)
{
  const double* sigma = pt.base.colptr(spec.col_sigma);

  const int nm = spec.n_mean_cols;

  for (int i = 0; i < n; i++) {
    if (!is_ok[i]) continue;
    double s = 0.0;
    for (int j = 0; j < nm; j++) s += mean_cols[j][i];
    ll_buf[i] = DNORM(y[i], s, sigma[i]);
  }
  vec_log(ll_buf, n);
}

void c_log_likelihood_MRI_ar1(const ParamTable& pt,
                              const MRISpec& spec,
                              const std::vector<const double*>& mean_cols,
                              const double* __restrict__ y,
                              int n,
                              double* __restrict__ ll_buf,
                              const std::vector<int>& is_ok)
{
  const double* sigma = pt.base.colptr(spec.col_sigma);
  const double* rho   = pt.base.colptr(spec.col_rho);

  const int nm = spec.n_mean_cols;

  // t = 0: stationary
  double s_prev = 0.0;
  for (int j = 0; j < nm; j++) s_prev += mean_cols[j][0];
  if (is_ok[0]) ll_buf[0] = DNORM(y[0], s_prev, sigma[0]);

  // t >= 1: AR(1) conditional
  for (int i = 1; i < n; i++) {
    double s_curr = 0.0;
    for (int j = 0; j < nm; j++) s_curr += mean_cols[j][i];

    if (is_ok[i]) {
      const double rho_i = rho[i];
      const double cond_sd = sigma[i] * std::sqrt(1.0 - rho_i * rho_i);
      const double cond_mean = s_curr + rho_i * (y[i - 1] - s_prev);
      ll_buf[i] = DNORM(y[i], cond_mean, cond_sd);
    }
    // The next conditional still needs this row's mean when this row is invalid.
    s_prev = s_curr;
  }
  vec_log(ll_buf, n);
}

namespace emc2 {

// =============================================================================
// LikelihoodContext — shared setup from calc_ll / calc_ll_multithreaded
// Kept between sampler calls; only the particle values change.
// =============================================================================

struct LikelihoodContext::Impl {
  enum Model { MRI, MRI_AR1, ORDERED_PROBIT, ORDERED_LOGIT, MULTINOMIAL_LOGIT, MULTINOMIAL_PROBIT,
               SSEXG, SSRDEX, CDM, PSDM, PHSDM, DDM, RACE };

  // Per-thread mutable state
  struct State {
    ParamTable pt_local;
    std::unique_ptr<TrendRuntime> tr_local;
    std::vector<int> is_ok;
    // ll_buf holds compressed trial log-likelihoods; ll_row holds accumulator rows.
    std::vector<double> ll_buf, ll_row, cut_buf;
    std::vector<const double*> mean_cols;
    RaceScratch scratch;
    CensorSpec censor;
    TruncSpec trunc;
    SSSpec ss;
    NumericMatrix pars; // CDM's native density functions use the Rcpp matrix ABI.
    double* pars_ptr = nullptr;
  };

  List fixed_inputs; // retains designs and custom compiled-kernel pointers
  PipelineContext ctx;
  PipelineCache cache;
  std::vector<TransformSpec> pretransforms;
  std::vector<BoundSpec> bound_specs;
  std::vector<std::pair<int, double>> constants;

  // Index lists — built once, shared across all choice models
  // idx_all: all non-missing rows (DDM)
  // idx_win / idx_los: winner/loser split (race and choice-only)
  std::vector<int> idx_all, idx_win, idx_los, keep_indices;
  std::vector<bool> participating;
  IntegerVector expand, R_vec, lR_vec;
  NumericVector rts, Rs, R2s, R3s, y_rcpp;
  const double* rts_ptr = nullptr;
  const double* y = nullptr;
  const int* Rs_ptr = nullptr;
  const int* lR_ptr = nullptr;
  const int* exp_ptr = nullptr;
  RaceModelSetup setup{};
  ChoiceOnlySpec choice;
  MRISpec mri;
  Model model;
  std::string model_name;
  int n_parameters, n_rows, n_lR, n_choice_trials, n_exp, n_threads_used;
  double min_ll;
  bool safe = true;

  Impl(CharacterVector parameter_names, DataFrame data, NumericVector constant_values,
       List designs, String type, List bounds, List transforms,
       List pretransform_values, double min_ll_, Nullable<List> trend, int n_threads)
    : fixed_inputs(List::create(data, designs, trend)), model_name(type),
      n_parameters(parameter_names.size()), n_rows(data.nrow()), min_ll(min_ll_)
  {
#ifdef _OPENMP
    n_threads_used = n_threads <= 0 ? omp_get_max_threads() : std::min(n_threads, omp_get_max_threads());
#else
    n_threads_used = 1;
#endif
    const bool has_lR = data.containsElementNamed("lR");
    n_lR = has_lR ? unique(as<IntegerVector>(data["lR"])).length() : 1;
    n_choice_trials = n_rows / n_lR;
    n_exp = n_choice_trials;
    participating.assign(n_rows, true);

    // Shared setup -- context holds the param_table, designs, constants and trends.
    NumericMatrix initial(1, n_parameters);
    colnames(initial) = parameter_names;
    pretransforms = make_transform_specs(initial, pretransform_values);
    ctx = make_pipeline_context(initial, data, constant_values, designs,
                                transforms, pretransform_values, trend);
    cache = make_pipeline_cache(ctx.param_table, designs, ctx.transform_specs,
                                ctx.trend_runtime.get());

    // Bounds — built once from structure
    NumericMatrix minmax = bounds["minmax"];
    bound_specs = make_bound_specs(minmax, colnames(minmax), ctx.param_table, bounds);
    for (int j = n_parameters; j < ctx.particle_matrix.ncol(); ++j)
      constants.emplace_back(ctx.pm_col_to_base_idx[j], ctx.particle_matrix(0, j));

    if (ctx.trend_plan) {
      for (const auto& kernel : ctx.trend_plan->kernels)
        if (kernel.second.kernel_type == KernelType::Custom) safe = false;
    }

    // -----------------------------------------------------------------------
    // MRI / MRI_AR1 -- no expand, rt, R attributes
    // -----------------------------------------------------------------------
    if (type == "MRI" || type == "MRI_AR1") {
      model = type == "MRI" ? MRI : MRI_AR1;
      y_rcpp = extract_y(data);
      y = y_rcpp.begin();
      // Resolve column indices once — no allocation inside particle loop
      mri = make_mri_spec(ctx.param_table, ctx.keep_names, model == MRI_AR1);
    } else {
      // All choice models — shared read-only data
      expand = data.attr("expand");
      n_exp = expand.size();
      exp_ptr = expand.begin();
      const bool has_missingness = data.containsElementNamed("missingness");
      const bool has_winner = data.containsElementNamed("winner");
      IntegerVector missingness = has_missingness ? as<IntegerVector>(data["missingness"]) : IntegerVector();
      LogicalVector winner = has_winner ? as<LogicalVector>(data["winner"]) : LogicalVector();
      idx_all.reserve(n_rows);
      if (has_winner) { idx_win.reserve(n_rows); idx_los.reserve(n_rows); }
      for (int i = 0; i < n_rows; ++i) {
        if (has_missingness && !IntegerVector::is_na(missingness[i])) {
          // Ignored rows also skip bounds; censored probabilities are filled separately.
          if (missingness[i] == 0) participating[i] = false;
          continue;
        }
        idx_all.push_back(i);
        if (has_winner) (winner[i] ? idx_win : idx_los).push_back(i);
      }

      // Choice-only models
      if (type == "ORDERED_PROBIT" || type == "ORDERED_LOGIT" ||
          type == "MULTINOMIAL_LOGIT" || type == "MULTINOMIAL_PROBIT") {
        model = type == "ORDERED_PROBIT" ? ORDERED_PROBIT :
                type == "ORDERED_LOGIT" ? ORDERED_LOGIT :
                type == "MULTINOMIAL_LOGIT" ? MULTINOMIAL_LOGIT : MULTINOMIAL_PROBIT;
        choice = make_choice_only_spec(ctx.param_table, std::string(type));
        lR_vec = data["lR"];
        lR_ptr = lR_vec.begin();
      } else if (type == "SSEXG" || type == "SSRDEX") {
        model = type == "SSEXG" ? SSEXG : SSRDEX;
      } else {
        rts = data["rt"];
        rts_ptr = rts.begin();
        // Continuous-choice-RT models (CDM, PSDM, PHSDM)
        if (type == "CDM" || type == "PSDM" || type == "PHSDM") {
          model = type == "CDM" ? CDM : type == "PSDM" ? PSDM : PHSDM;
          safe = false;
          Rs = data["R"];
          if (data.containsElementNamed("R2")) R2s = data["R2"];
          if (data.containsElementNamed("R3")) R3s = data["R3"];
          for (const auto& name : ctx.keep_names)
            keep_indices.push_back(ctx.param_table.base_index_for(as<std::string>(name)));
        } else {
          // Discrete-Choice-RT models (DDM, & Race: RDM, LBA, LNR, ...)
          model = type == "DDM" ? DDM : RACE;
          setup = make_race_setup(type, ctx.param_table);
          if (model == DDM) {
            R_vec = data["R"];
            Rs_ptr = R_vec.begin();
          }
        }
      }
    }
    // These existing native APIs allocate R objects and must run on the R thread.
    if (!safe) n_threads_used = 1;

  }

  // ---------------------------------------------------------------------------
  // Per-thread mutable state — built on the R thread, bound to its own ParamTable
  // ---------------------------------------------------------------------------
  std::vector<std::unique_ptr<State>> make_states() const {
    DataFrame data = as<DataFrame>(fixed_inputs[0]);
    String type(model_name);
    std::vector<std::unique_ptr<State>> states;
    states.reserve(n_threads_used);
    for (int t = 0; t < n_threads_used; ++t) {
      states.emplace_back(new State);
      State& s = *states.back();
      // TrendPlan is shared; kernel buffers and parameter pointers are local.
      s.pt_local = ctx.param_table.deep_copy();
      if (ctx.trend_runtime)
        s.tr_local.reset(new TrendRuntime(clone_trend_runtime(*ctx.trend_runtime, s.pt_local)));
      s.is_ok.assign(n_rows, 1);
      s.ll_buf.resize(n_choice_trials);
      if (model == MRI || model == MRI_AR1) {
        // Pre-fetch all mean column pointers
        for (int col : mri.col_means) s.mean_cols.push_back(s.pt_local.base.colptr(col));
      } else if (model == ORDERED_PROBIT || model == ORDERED_LOGIT) {
        s.cut_buf.resize(n_rows);
      } else if (model == SSEXG || model == SSRDEX) {
        // The SS likelihood handles missingness codes itself (deadline / lower
        // censoring / intrinsic no-response), so no CensorSpec/TruncSpec here.
        s.ss = make_ss_spec(std::string(type), data, s.pt_local);
      } else if (model == CDM || model == PSDM || model == PHSDM) {
        s.pars = NumericMatrix(n_rows, keep_indices.size());
        s.pars_ptr = s.pars.begin();
      } else if (model == DDM || model == RACE) {
        auto slots = race_scratch_slots(type);
        s.scratch.reserve(slots.first, slots.second, n_rows);
        s.censor = make_censor_spec(data, n_choice_trials, n_lR, participating, &setup, s.pt_local, s.scratch);
        s.trunc = make_trunc_spec(data, n_choice_trials, n_lR, participating, &setup, s.pt_local, s.scratch);
        if (model == RACE) s.ll_row.resize(n_rows);
      }
    }
    return states;
  }

  // Pretransform + constants, formerly applied to the whole particle_matrix in
  // make_pipeline_context. Fill one particle directly into the cached ParamTable.
  void fill_particle(State& s, const double* particles, int row, int n_particles) const {
    s.pt_local.reset_base_to_zero();
    for (const TransformSpec& sp : pretransforms) {
      const int col = ctx.pm_col_to_base_idx[sp.col_idx];
      if (col < 0) continue;
      double value = particles[sp.col_idx * n_particles + row];
      if (sp.code == EXP) value = sp.lower + std::exp(value);
      else if (sp.code == PNORM)
        value = sp.lower + (sp.upper - sp.lower) * PNORM_STD(value, true, false);
      std::fill_n(s.pt_local.base.colptr(col), n_rows, value);
    }
    for (const auto& constant : constants)
      if (constant.first >= 0) std::fill_n(s.pt_local.base.colptr(constant.first), n_rows, constant.second);
    run_pars_pipeline(s.pt_local, s.tr_local.get(), cache);
  }

  // ---------------------------------------------------------------------------
  // One particle — common pipeline for serial and parallel evaluation
  // ---------------------------------------------------------------------------
  double evaluate_particle(State& s, const double* particles, int row,
                           int n_particles, double* trialwise) const {
    // 1) Map p_vector to trialwise parameters
    fill_particle(s, particles, row, n_particles);
    // Bounds apply after all mappings and trends. A participating accumulator
    // outside its bounds invalidates the whole trial, including its integrals.
    c_do_bound(s.pt_local, bound_specs, s.is_ok);
    const bool continuous = model == MRI || model == MRI_AR1 || model == CDM || model == PSDM || model == PHSDM;
    const int stride = continuous ? 1 : n_lR;
    for (int t = 0; t < n_choice_trials; ++t) {
      const int base = t * stride;
      bool ok = true;
      for (int k = 0; k < stride; ++k)
        if (participating[base + k] && !s.is_ok[base + k]) { ok = false; break; }
      std::fill_n(s.is_ok.begin() + base, stride, ok);
    }

    // 2) Calculate raw (compressed) trialwise log-likelihoods
    // Preserve the existing initial values for skipped trials.
    std::fill(s.ll_buf.begin(), s.ll_buf.end(), model == DDM ? 1.0 : 0.0);

    switch (model) {
    case MRI:
      c_log_likelihood_MRI_white(s.pt_local, mri, s.mean_cols, y, n_choice_trials, s.ll_buf.data(), s.is_ok);
      break;
    case MRI_AR1:
      c_log_likelihood_MRI_ar1(s.pt_local, mri, s.mean_cols, y, n_choice_trials, s.ll_buf.data(), s.is_ok);
      break;
    case ORDERED_PROBIT:
    case ORDERED_LOGIT:
      c_log_likelihood_ordered(s.pt_local, choice, idx_win, lR_ptr, n_rows, n_lR,
                               s.cut_buf.data(), model == ORDERED_PROBIT, s.ll_buf.data(), s.is_ok);
      break;
    case MULTINOMIAL_LOGIT:
      c_log_likelihood_multinomial_logit(s.pt_local, choice, idx_win, n_lR, s.ll_buf.data(), s.is_ok);
      break;
    case MULTINOMIAL_PROBIT:
      c_log_likelihood_multinomial_probit(s.pt_local, choice, idx_win, n_lR, s.ll_buf.data(), s.is_ok);
      break;
    case SSEXG:
    case SSRDEX:
      ss_log_likelihood(s.ss, min_ll, s.ll_buf.data(), s.is_ok);
      break;
    case CDM:
    case PSDM:
    case PHSDM:
      for (size_t j = 0; j < keep_indices.size(); ++j)
        std::copy_n(s.pt_local.base.colptr(keep_indices[j]), n_rows, s.pars_ptr + j * n_rows);
      if (model == CDM) {
        // Sub-type dispatch: CDM -> HSDM (R,R2,R3) -> SDM (R,R2) -> CDM (R only)
        if (R2s.size() && R3s.size() && keep_indices.size() >= 8)
          c_dHSDM(rts, Rs, R2s, R3s, s.pars, participating, s.ll_buf.data(), n_choice_trials, s.is_ok);
        else if (R2s.size() && keep_indices.size() >= 7)
          c_dSDM(rts, Rs, R2s, s.pars, participating, s.ll_buf.data(), n_choice_trials, s.is_ok);
        else c_dCDM(rts, Rs, s.pars, participating, s.ll_buf.data(), n_choice_trials, s.is_ok);
      } else if (model == PSDM)
        c_dPSDM(rts, Rs, s.pars, participating, s.ll_buf.data(), n_choice_trials, s.is_ok);
      else c_dPHSDM(rts, Rs, R2s, s.pars, participating, s.ll_buf.data(), n_choice_trials, s.is_ok);
      break;
    case DDM:
      c_log_likelihood_DDM(rts_ptr, Rs_ptr, s.pt_local, setup.spec, idx_all, s.ll_buf.data(), s.is_ok);
      break;
    case RACE:
      // Re-fill row-wise likelihoods — skipped accumulators contribute 1.
      std::fill(s.ll_row.begin(), s.ll_row.end(), 1.0);
      c_log_likelihood_race(s.pt_local, setup, rts_ptr, idx_win, idx_los, n_lR,
                            s.ll_row.data(), n_rows, s.ll_buf.data(), s.scratch, s.is_ok, s.trunc.any());
      break;
    }

    // 3) Handle truncation, censoring. Needs to be in this order!
    if (model == DDM || model == RACE) {
      if (s.trunc.any()) s.trunc.calculate_normalization_constant(s.is_ok);
      // Censoring needs the normalization constant and fills rows, doesn't add.
      if (s.censor.any()) s.censor.fill_censored_rows(s.trunc, s.is_ok, s.ll_buf);
      if (s.trunc.any())
        for (int t = 0; t < n_choice_trials; ++t)
          s.ll_buf[t] = is_finite(s.trunc.log_Z[t]) ? s.ll_buf[t]-s.trunc.log_Z[t] : min_ll;
    }

    // 4) Handle not-ok parameter values (out of bound)
    for (int t = 0; t < n_choice_trials; ++t)
      if (!s.is_ok[t * stride]) s.ll_buf[t] = min_ll;

    // 5) Expand, clamp, sum. MRI has no expand.
    return (model == MRI || model == MRI_AR1) ?
      clamp_sum(s.ll_buf.data(), n_choice_trials, min_ll, trialwise) :
      expand_clamp_sum(s.ll_buf.data(), exp_ptr, n_exp, min_ll, trialwise);
  }

  // Column i = particle i's trialwise buffer (contiguous).
  void evaluate(std::vector<std::unique_ptr<State>>& states, const double* particles, int n_particles, double* result, bool trialwise) const {
    if (n_threads_used == 1 || n_particles == 1) {
      for (int i = 0; i < n_particles; ++i) {
        const double sum = evaluate_particle(*states[0], particles, i, n_particles,
                                             trialwise ? result + (ptrdiff_t)i * n_exp : nullptr);
        if (!trialwise) result[i] = sum;
      }
    } else {
#pragma omp parallel for schedule(static) num_threads(n_threads_used)
      for (int i = 0; i < n_particles; ++i) {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        const double sum = evaluate_particle(*states[tid], particles, i, n_particles,
                                             trialwise ? result + (ptrdiff_t)i * n_exp : nullptr);
        if (!trialwise) result[i] = sum;
      }
    }
  }
};

// Copies share Impl's fixed inputs, but need their own buffers and bindings.
struct LikelihoodContext::Workspace {
  std::vector<std::unique_ptr<Impl::State>> states;
  explicit Workspace(const Impl& fixed) : states(fixed.make_states()) {}
};

LikelihoodContext::LikelihoodContext(CharacterVector parameter_names, DataFrame data,
                                     NumericVector constants, List designs, String type,
                                     List bounds, List transforms, List pretransforms,
                                     double min_ll, Nullable<List> trend, int n_threads)
  : impl_(new Impl(parameter_names, data, constants, designs, type, bounds,
                   transforms, pretransforms, min_ll, trend, n_threads)),
    work_(new Workspace(*impl_)) {}
LikelihoodContext::LikelihoodContext(const LikelihoodContext& other)
  : impl_(other.impl_), work_(new Workspace(*impl_)) {}
LikelihoodContext::~LikelihoodContext() = default;
LikelihoodContext::LikelihoodContext(LikelihoodContext&&) noexcept = default;
LikelihoodContext& LikelihoodContext::operator=(LikelihoodContext&&) noexcept = default;
void LikelihoodContext::evaluate(const double* particles, int n_particles,
                                 double* result, bool trialwise) {
  impl_->evaluate(work_->states, particles, n_particles, result, trialwise);
}
int LikelihoodContext::trial_count() const { return impl_->n_exp; }
bool LikelihoodContext::thread_safe() const { return impl_->safe; }

} // namespace emc2

// [[Rcpp::export]]
NumericMatrix calc_ll(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                      List designs, String type, List bounds, List transforms, List pretransforms,
                      CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                      bool return_trialwise = false) {
  emc2::LikelihoodContext context(colnames(particle_matrix), data, constants, designs, type,
                                  bounds, transforms, pretransforms, min_ll, trend);
  NumericMatrix result(return_trialwise ? context.trial_count() : 1, particle_matrix.nrow());
  context.evaluate(particle_matrix.begin(), particle_matrix.nrow(), result.begin(), return_trialwise);
  return result;
}

// [[Rcpp::export]]
NumericMatrix calc_ll_multithreaded(NumericMatrix particle_matrix, DataFrame data,
                                    NumericVector constants, List designs, String type,
                                    List bounds, List transforms, List pretransforms,
                                    CharacterVector p_types, double min_ll,
                                    Rcpp::Nullable<Rcpp::List> trend = R_NilValue,
                                    bool return_trialwise = false, int n_threads = -1) {
  emc2::OpenMPPoolRelease release;
#ifndef _OPENMP
  if (n_threads > 1)
    Rcpp::warning("calc_ll_multithreaded: OpenMP not available, running single-threaded.");
#endif
  emc2::LikelihoodContext context(colnames(particle_matrix), data, constants, designs, type,
                                  bounds, transforms, pretransforms, min_ll, trend, n_threads);
  NumericMatrix result(return_trialwise ? context.trial_count() : 1, particle_matrix.nrow());
  context.evaluate(particle_matrix.begin(), particle_matrix.nrow(), result.begin(), return_trialwise);
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
