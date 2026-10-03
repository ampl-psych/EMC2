#ifndef SUBJECT_PIPELINE_H
#define SUBJECT_PIPELINE_H

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
    TrendRuntime* trend_runtime_ptr,
    bool compute_col_is_constant = true);


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
  int                            n_active_trials;
};

PipelineContext make_pipeline_context(
    Rcpp::NumericMatrix particle_matrix,
    const Rcpp::DataFrame& data,
    const Rcpp::NumericVector& constants,
    const Rcpp::List& designs,
    const Rcpp::List& transforms,
    const Rcpp::List& pretransforms,
    const Rcpp::Nullable<Rcpp::List>& trend,
    const int n_active_trials = -1);


void run_pars_pipeline(ParamTable&          param_table,
                       TrendRuntime*        trend_runtime,
                       const PipelineCache& cache,
                       int row_start = 0,
                       int row_end   = -1);


// For persistent parameter mapping pipelines
struct SubjectPipeline {
  ParamTable                  param_table;
  std::vector<TransformSpec>  transform_specs;
  std::vector<double> particle_values;
  std::vector<int>    pm_col_to_base_idx;

  std::unique_ptr<TrendPlan>    trend_plan;
  std::unique_ptr<TrendRuntime> trend_runtime;

  PipelineCache cache;

  int n_trials  = 0;
  int n_acc     = 0;

  bool has_trend() const { return trend_runtime != nullptr; }
};


SEXP create_subject_pipeline(
    Rcpp::NumericMatrix                pars,
    const Rcpp::List&                  designs,
    const Rcpp::List&                  transform,
    const Rcpp::DataFrame&             data,
    const Rcpp::NumericVector&         constants,
    const Rcpp::List&                  pretransform,
    const Rcpp::Nullable<Rcpp::List>& trend = R_NilValue);

void step_subject_pipeline(
    SEXP                   xptr,
    const Rcpp::List&      new_designs,
    const Rcpp::DataFrame& new_data,
    int                    row_start,
    int                    row_end);

Rcpp::NumericMatrix get_subject_pipeline_result(SEXP xptr);

#endif
