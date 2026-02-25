#include <Rcpp.h>
#include <unordered_set>
#include "kernels.h"
#include "ParamTable.h"

// Phase enum
enum class TrendPhase { Premap, Pretransform, Posttransform };

inline TrendPhase parse_phase(const std::string& ph) {
  if (ph == "premap")        return TrendPhase::Premap;
  if (ph == "pretransform")  return TrendPhase::Pretransform;
  if (ph == "posttransform") return TrendPhase::Posttransform;
  Rcpp::stop("Unknown trend phase: '%s'", ph.c_str());
}

// ---- Specification side ----
enum class InputKind { None, Covariate, ParInput, MultiCovariate, MultiParInput, Combined};


// One TrendOp = one base + one group of kernels. Groups of kernels are combined
struct KernelSlotSpec {
  std::string kernel;
  KernelType kernel_type;

  InputKind input_kind = InputKind::None;

  // For InputKind::Covariate
  Rcpp::NumericMatrix kernel_input;  // data column
  std::string kernel_input_name;

  // For InputKind::ParInput
  std::string par_input_name;

  // For InputKind::Combined
  Rcpp::CharacterVector covariate_names;
  std::vector<int> covariate_indices;
  Rcpp::CharacterVector par_input_names;
  std::vector<int> par_input_indices;

  // For custom kernels
  SEXP custom_fun = R_NilValue;

  // Covariate maps specific to this kernel
  bool has_covariate_maps = false;
  std::vector<Rcpp::NumericVector> covariate_map_cols;
};

// One trend operation spec, built from one element of the R 'trend' list
struct TrendOpSpec {
  TrendPhase   phase;
  std::string  target_param;
  std::string  base_type;

  // Shared across all kernels in this op:
  bool has_at = false;
  std::string at;
  Rcpp::LogicalVector first_level;
  std::vector<int> expand_idx; // length = n_trials
  std::vector<int> comp_index; // compressed row indices

  // Union of all trend parameter names for this op
  Rcpp::CharacterVector trend_pnames;

  // Kernels belonging to this TrendOp
  std::vector<KernelSlotSpec> kernels;

  // Do the kernels in this Op share some latent trace?
  // Only applicable to Pearce-Hall and Delta-Risk, so far
  bool shared_latent = false;

  // maps at Op level
  bool has_covariate_maps = false;
  // one column per map, pre-filtered to the relevant covariate
  std::vector<Rcpp::NumericVector> covariate_map_cols;


  TrendOpSpec() = default;

  void make_first_level(const Rcpp::DataFrame& data);
};

// Plan: all spec operations + parameter sets
struct TrendPlan {
  Rcpp::List trend;       // original trend list or R_NilValue
  Rcpp::DataFrame data;

  Rcpp::List data_covariate_maps;
  bool has_data_covariate_maps = false;

  std::vector<TrendOpSpec> premap_ops;
  std::vector<TrendOpSpec> pretransform_ops;
  std::vector<TrendOpSpec> posttransform_ops;

  std::unordered_set<std::string> premap_trend_params;
  std::unordered_set<std::string> pretransform_trend_params;
  std::unordered_set<std::string> posttransform_trend_params;
  std::unordered_set<std::string> all_trend_params;

  TrendPlan(Rcpp::Nullable<Rcpp::List> trend_, const Rcpp::DataFrame& data_);

  bool has_premap()        const { return !premap_ops.empty(); }
  bool has_pretransform()  const { return !pretransform_ops.empty(); }
  bool has_posttransform() const { return !posttransform_ops.empty(); }

  Rcpp::LogicalVector premap_design_mask(const Rcpp::List& designs) const;

  const std::unordered_set<std::string>& premap_trend_params_set() const {
    return premap_trend_params;
  }
  const std::unordered_set<std::string>& pretransform_trend_params_set() const {
    return pretransform_trend_params;
  }

  const std::unordered_set<std::string>& posttransform_trend_params_set() const {
    return posttransform_trend_params;
  }
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
    return out; // no field at all
  }

  SEXP cov_spec = tr["covariate"];

  // NULL => no covariate
  if (Rf_isNull(cov_spec)) {
    return out;
  }

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

  Rcpp::stop("covariate_names_from_spec: type of covariate not understood. "
              "It should be either a string or a vector of strings");
}

inline std::vector<std::string> par_input_names_from_spec(const Rcpp::List& tr) {
  std::vector<std::string> out;

  if (!tr.containsElementNamed("par_input")) {
    return out;
  }

  SEXP par_input_spec = tr["par_input"];

  // NULL => no par_input
  if (Rf_isNull(par_input_spec)) {
    return out;
  }

  // Case 1: single string: "p1"
  if (Rf_isString(par_input_spec) && Rf_length(par_input_spec) == 1) {
    out.push_back(Rcpp::as<std::string>(par_input_spec));
    return out;
  }

  // Case 2: character vector: c("p1", "p2", ...)
  if (TYPEOF(par_input_spec) == STRSXP) {
    Rcpp::CharacterVector par_input(par_input_spec);
    out.reserve(par_input.size());
    for (int i = 0; i < par_input.size(); ++i) {
      out.push_back(Rcpp::as<std::string>(par_input[i]));
    }
    return out;
  }

  Rcpp::stop("par_input_names_from_spec: type of par_input not understood. "
               "It should be either a string or a vector of strings");

  // Case 3: non‑string input (e.g., numeric vector supplied inline):
  // return out;  // empty => we will keep the original spec as a single TrendOp
}


// ---- Runtime side ----

struct KernelSlotRuntime {
  const KernelSlotSpec* spec;               // non-owning
  std::unique_ptr<BaseKernel> kernel_ptr;   // Delta, Poly*, etc.
  std::vector<int> kernel_par_indices;      // columns in ParamTable::base
};


// Runtime info for a single trend op: binds spec to ParamTable and kernel
// One TrendOp (one base + many kernels) at runtime
struct TrendOpRuntime {
  const TrendOpSpec* spec;  // one base + K kernels

  // base parameters for this op (from spec->trend_pnames)
  std::vector<int> base_par_indices;
  int base_par_idx = -1;    // convenience when only one base param

  std::vector<KernelSlotRuntime> kernels; // 1..K

  // Latent trace of alpha (only for Pearce-Hall atm)
  std::vector<double> shared_alpha;
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
  const std::unordered_set<std::string>& pretransform_trend_params() const {
    return plan->pretransform_trend_params_set();
  }
  const std::unordered_set<std::string>& posttransform_trend_params() const {
    return plan->posttransform_trend_params_set();
  }
  const std::unordered_set<std::string>& all_trend_params() const {
    return plan->all_trend_params;
  }

  Rcpp::LogicalVector premap_design_mask(const Rcpp::List& designs) const {
    return plan->premap_design_mask(designs);
  }

  void bind_all_ops_to_paramtable(const ParamTable& pt);
  void run_kernels_for_op(TrendOpRuntime& op, ParamTable& pt);
  void apply_base_for_op(TrendOpRuntime& op, ParamTable& pt);
  void reset_all_kernels();

  Rcpp::NumericMatrix all_kernel_outputs(ParamTable& pt);
  Rcpp::NumericMatrix all_kernel_outputs(ParamTable& pt, const std::vector<int>& codes);
};

// helper to bind one op at runtime
inline void bind_trend_op_to_paramtable(TrendOpRuntime& op,
                                        const ParamTable& pt);


// // helpers for TrendPlan construction
bool uses_cov_input(const Rcpp::List& tr);

bool uses_par_input(const Rcpp::List& tr);

void init_covariate_for_slot(KernelSlotSpec& slot,
                             const Rcpp::List& spec,
                             const Rcpp::DataFrame& data);

void init_combined_for_slot(KernelSlotSpec& slot,
                             const Rcpp::List& spec,
                             const Rcpp::DataFrame& data);

void init_covariate_maps_for_op(TrendOpSpec& op,
                                const Rcpp::List& tr_i,
                                const Rcpp::DataFrame& data,
                                const Rcpp::List& data_covmaps);
