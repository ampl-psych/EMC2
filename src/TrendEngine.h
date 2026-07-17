#ifndef TREND_ENGINE_H
#define TREND_ENGINE_H

#include <Rcpp.h>          // boundary only — not used in hot paths
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <memory>
#include "kernels.h"
#include "ParamTable.h"
#include "Mat.h"

// =============================================================================
// Phase
// =============================================================================

enum class TrendPhase { Premap, Pretransform, Posttransform };

inline TrendPhase parse_phase(const std::string& ph) {
  if (ph == "premap")        return TrendPhase::Premap;
  if (ph == "pretransform")  return TrendPhase::Pretransform;
  if (ph == "posttransform") return TrendPhase::Posttransform;
  Rf_error("Unknown trend phase: '%s'", ph.c_str());
  return TrendPhase::Premap; // unreachable
}

// =============================================================================
// KernelSpec  —  mirrors emc2_kernel; pure C++ after construction
// =============================================================================

struct KernelSpec {
  // identity
  std::string kernel_id;
  std::string kernel;
  KernelType  kernel_type;

  // phase: earliest phase among all bases referencing this kernel (set by R)
  TrendPhase phase;

  // inputs
  std::vector<std::string> cov_names;
  std::vector<std::string> par_input;

  // 'at' filter
  bool        has_at = false;
  std::string at;

  // finalised prefixed parameter names
  std::vector<std::string> pnames;

  // sequential flag
  bool sequential = false;

  // kernel_args
  std::vector<int> q_reset_col;   // length n_trials, or empty
  std::vector<int> belief_reset_col;
  KernelArgs       kernel_args;   // raw-pointer view; rebuilt via build_kernel_args()

  // custom kernel pointer (R_NilValue if not custom)
  SEXP custom_fun = R_NilValue;

  // data-derived fields — plain C++ after construction
  mutable Mat              kernel_input;   // n_trials x (n_cov + n_par)
  std::vector<bool>        first_level;    // length n_trials
  std::vector<int>         expand_idx;     // length n_trials
  std::vector<int>         comp_index;     // first-level row indices

  std::vector<int>         covariate_indices;
  std::vector<int>         par_input_indices;

  void build_kernel_args() {
    kernel_args = KernelArgs{};
    if (!q_reset_col.empty())
      kernel_args.q_reset = q_reset_col.data();
    if (!belief_reset_col.empty())
      kernel_args.belief_reset = belief_reset_col.data();
  }
};

// =============================================================================
// BaseSpec  —  mirrors emc2_base; pure C++ after construction
// =============================================================================

struct BaseSpec {
  std::string base_type;
  std::string target_parameter;
  std::string kernel_id;
  int         kernel_output = 1;
  TrendPhase  phase;

  std::vector<std::string> pnames;

  // covariate coding: always matrices (column-major, n_trials x n_covs)
  // one Mat per named map entry; empty if no coding schemes
  bool             has_covariate_coding = false;
  std::vector<Mat> covariate_coding;
};

// =============================================================================
// TrendPlan  —  owns all specs; built once per subject
// =============================================================================

struct TrendPlan {
  std::unordered_map<std::string, KernelSpec> kernels;

  std::vector<BaseSpec> premap_bases;
  std::vector<BaseSpec> pretransform_bases;
  std::vector<BaseSpec> posttransform_bases;

  std::unordered_set<std::string> premap_params;
  std::unordered_set<std::string> pretransform_params;
  std::unordered_set<std::string> posttransform_params;
  std::unordered_set<std::string> all_trend_params;

  TrendPlan() = default;

  // Rcpp boundary: only place Rcpp types are used
  TrendPlan(const Rcpp::List& trend, const Rcpp::DataFrame& data);

  bool has_premap()        const { return !premap_bases.empty(); }
  bool has_pretransform()  const { return !pretransform_bases.empty(); }
  bool has_posttransform() const { return !posttransform_bases.empty(); }

  // Returns a LogicalVector (Rcpp boundary — used by make_pipeline_cache)
  Rcpp::LogicalVector premap_design_mask(const Rcpp::List& designs) const;

  const std::unordered_set<std::string>& premap_trend_params()       const { return premap_params; }
  const std::unordered_set<std::string>& pretransform_trend_params()  const { return pretransform_params; }
  const std::unordered_set<std::string>& posttransform_trend_params() const { return posttransform_params; }
};

// =============================================================================
// Runtime structs — pure C++
// =============================================================================

struct KernelRuntime {
  const KernelSpec* spec = nullptr;
  std::vector<std::unique_ptr<BaseKernel>> kernel_ptrs;  // one per cov slot
  std::vector<int>                         kernel_par_indices;

  // one pre-allocated 1-column Mat per slot (arity-1 only)
  // variadic: empty, kernel uses ks.kernel_input directly
  std::vector<Mat> slot_inputs;

  // which slots are par_input (need refilling per particle)
  // index into slot_inputs / kernel_ptrs
  std::vector<int> par_input_slot_indices;  // slot index -> par_input index
  // parallel: which par_input string to pull from ParamTable
  std::vector<int> par_input_param_col;     // slot index -> par_input[i]

  // convenience: number of slots
  int n_slots() const { return (int)kernel_ptrs.size(); }

  // convenience: Is this a variatic kernel? Ie., does it require *multiple covariates/par_input* to run?
  bool is_variadic() const {
    if (!spec) return false;
    KernelMeta m = kernel_meta(spec->kernel_type);
    return (m.input_arity == -1 || spec->kernel_type == KernelType::Custom);
  }
};

struct BaseRuntime {
  const BaseSpec*    spec      = nullptr;
  KernelRuntime*     kernel_rt = nullptr;  // non-owning, into TrendRuntime::kernels
  std::vector<int>   base_par_indices;
};

struct TrendRuntime {
  const TrendPlan* plan = nullptr;

  std::unordered_map<std::string, KernelRuntime> kernels;

  std::vector<BaseRuntime> premap_bases;
  std::vector<BaseRuntime> pretransform_bases;
  std::vector<BaseRuntime> posttransform_bases;

  explicit TrendRuntime(const TrendPlan& plan);

  bool has_premap()        const { return plan->has_premap(); }
  bool has_pretransform()  const { return plan->has_pretransform(); }
  bool has_posttransform() const { return plan->has_posttransform(); }

  const std::unordered_set<std::string>& premap_trend_params()       const { return plan->premap_params; }
  const std::unordered_set<std::string>& pretransform_trend_params()  const { return plan->pretransform_params; }
  const std::unordered_set<std::string>& posttransform_trend_params() const { return plan->posttransform_params; }
  const std::unordered_set<std::string>& all_trend_params()           const { return plan->all_trend_params; }

  Rcpp::LogicalVector premap_design_mask(const Rcpp::List& designs) const {
    return plan->premap_design_mask(designs);
  }

  void bind_all_to_paramtable(const ParamTable& pt);
  void apply_base(BaseRuntime& base_rt, ParamTable& pt);
  void reset_all_kernels();

  // Rcpp boundary — diagnostic only, not in hot path
  Rcpp::NumericMatrix all_kernel_outputs(ParamTable& pt);
  Rcpp::NumericMatrix all_kernel_outputs(ParamTable& pt, const std::vector<int>& codes);

private:
  void run_kernel(KernelRuntime& k_rt, ParamTable& pt);
};

// =============================================================================
// Inline helper
// =============================================================================

inline KernelParsView make_kernel_pars_view(const ParamTable& pt,
                                            const std::vector<int>& col_indices)
{
  KernelParsView view;
  view.n_rows = pt.n_trials;
  view.cols.resize(col_indices.size());
  for (size_t i = 0; i < col_indices.size(); ++i)
    view.cols[i] = &pt.base(0, col_indices[i]);
  return view;
}

#endif // TREND_ENGINE_H
