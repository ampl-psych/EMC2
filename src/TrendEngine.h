#include <Rcpp.h>
#include <unordered_set>
#include "kernels.h"
#include "ParamTable.h"

// Phase enum (reuse existing)
enum class TrendPhase { Premap, Pretransform, Posttransform };

inline TrendPhase parse_phase(const std::string& ph) {
  if (ph == "premap")        return TrendPhase::Premap;
  if (ph == "pretransform")  return TrendPhase::Pretransform;
  if (ph == "posttransform") return TrendPhase::Posttransform;
  Rcpp::stop("Unknown trend phase: '%s'", ph.c_str());
}

// ---- Specification side ----

// One trend operation spec, built from one element of the R 'trend' list
struct TrendOpSpec {
  TrendPhase phase;
  std::string kernel;                 // e.g. "lin_incr"
  std::string base_type;              // e.g. "lin"
  std::string target_param;

  bool has_at = false;
  std::string at;
  Rcpp::LogicalVector first_level;
  std::vector<int> expand_idx;        // length = n_trials
  std::vector<int> comp_index;        // compressed row indices

  Rcpp::CharacterVector trend_pnames;
  Rcpp::CharacterVector par_input;
  Rcpp::List spec;                    // original R spec

  KernelType kernel_type;

  // Cached covariate. Since it's now in the spec, it can be re-used across particles
  Rcpp::NumericVector covariate;

  TrendOpSpec(const Rcpp::List& cur, const std::string& name_from_list);

  void make_first_level(const Rcpp::DataFrame& data);
};

// Plan: all spec operations + parameter sets
struct TrendPlan {
  Rcpp::List trend;       // original trend list or R_NilValue
  Rcpp::DataFrame data;

  std::vector<TrendOpSpec> premap_ops;
  std::vector<TrendOpSpec> pretransform_ops;
  std::vector<TrendOpSpec> posttransform_ops;

  std::unordered_set<std::string> premap_trend_params;
  std::unordered_set<std::string> all_trend_params;

  TrendPlan(Rcpp::Nullable<Rcpp::List> trend_, const Rcpp::DataFrame& data_);

  bool has_premap()        const { return !premap_ops.empty(); }
  bool has_pretransform()  const { return !pretransform_ops.empty(); }
  bool has_posttransform() const { return !posttransform_ops.empty(); }

  Rcpp::LogicalVector premap_design_mask(const Rcpp::List& designs) const;
};


inline KernelParsView make_kernel_pars_view(const ParamTable& pt,
                                            const std::vector<int>& col_indices)
{
  KernelParsView view;
  view.n_rows = pt.n_trials;
  view.cols.resize(col_indices.size());
  for (size_t k = 0; k < col_indices.size(); ++k) {
    int base_idx = col_indices[k];
    view.cols[k] = &pt.base(0, base_idx);  // pointer to first element of that column
  }
  return view;
}


inline std::vector<std::string> covariate_names_from_spec(const Rcpp::List& tr) {
  std::vector<std::string> out;

  if (!tr.containsElementNamed("covariate")) {
    return out; // no covariate; will be handled as error later
  }

  SEXP cov_spec = tr["covariate"];

  // Case 1: single string: "cov1"
  if (Rf_isString(cov_spec) && Rf_length(cov_spec) == 1) {
    out.push_back(Rcpp::as<std::string>(cov_spec));
    return out;
  }

  // Case 2: character vector: c("cov1", "cov2", ...)
  if (TYPEOF(cov_spec) == STRSXP) {
    Rcpp::CharacterVector cv(cov_spec);
    out.reserve(cv.size());
    for (int i = 0; i < cv.size(); ++i) {
      out.push_back(Rcpp::as<std::string>(cv[i]));
    }
    return out;
  }

  // Case 3: non‑string covariate (e.g. numeric vector supplied inline):
  // treat as "inline covariate", not multiple columns
  return out;  // empty => we will keep the original spec as a single TrendOp
}



// ---- Runtime side ----

// Runtime info for a single trend op: binds spec to ParamTable and kernel
struct TrendOpRuntime {
  const TrendOpSpec* spec;              // non-owning pointer into TrendPlan
  std::unique_ptr<BaseKernel> kernel_ptr;
  int base_par_idx = -1;                // column in ParamTable::base, or -1 if none
  std::vector<int> kernel_par_indices;  // columns in ParamTable::base

  TrendOpRuntime(const TrendOpSpec* s)
    : spec(s),
      kernel_ptr(make_kernel(s->kernel_type)) {}
};

struct TrendRuntime {
  const TrendPlan* plan;  // non-owning pointer

  std::vector<TrendOpRuntime> premap_ops;
  std::vector<TrendOpRuntime> pretransform_ops;
  std::vector<TrendOpRuntime> posttransform_ops;

  TrendRuntime(const TrendPlan& plan_);

  bool has_premap()        const { return plan->has_premap(); }
  bool has_pretransform()  const { return plan->has_pretransform(); }
  bool has_posttransform() const { return plan->has_posttransform(); }

  const std::unordered_set<std::string>& premap_trend_params() const {
    return plan->premap_trend_params;
  }
  const std::unordered_set<std::string>& all_trend_params() const {
    return plan->all_trend_params;
  }

  Rcpp::LogicalVector premap_design_mask(const Rcpp::List& designs) const {
    return plan->premap_design_mask(designs);
  }

  void bind_all_ops_to_paramtable(const ParamTable& pt);
  void run_kernel_for_op(TrendOpRuntime& op, ParamTable& pt);
  void apply_base_for_op(TrendOpRuntime& op, ParamTable& pt);
};

// helper to bind one op at runtime
inline void bind_trend_op_to_paramtable(TrendOpRuntime& op,
                                        const ParamTable& pt);

