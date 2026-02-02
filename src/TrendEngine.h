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

// One trend operation, built from one element of the R 'trend' list
struct TrendOp {
  TrendPhase phase;
  std::string kernel;                 // e.g. "lin_incr"
  std::string base_type;              // e.g. "lin"
  std::string target_param;
  Rcpp::CharacterVector trend_pnames;
  Rcpp::CharacterVector par_input;
  Rcpp::List spec;

  // Kernel & base - ptr and parameters
  KernelType kernel_type;
  std::unique_ptr<BaseKernel> kernel_ptr;
  int base_par_idx = -1;                 // column in ParamTable::base, or -1 if none
  std::vector<int> kernel_par_indices;   // columns in ParamTable::base

  TrendOp(const Rcpp::List& cur, const std::string& name_from_list)
    : spec(cur) {
    using Rcpp::as;

    phase        = parse_phase(as<std::string>(cur["phase"]));
    kernel       = as<std::string>(cur["kernel"]);
    base_type    = as<std::string>(cur["base"]);
    target_param = name_from_list;

    if (cur.containsElementNamed("trend_pnames")) {
      trend_pnames = cur["trend_pnames"];
    } else {
      trend_pnames = Rcpp::CharacterVector(0);
    }

    if (cur.containsElementNamed("par_input")) {
      par_input = cur["par_input"];
    } else {
      par_input = Rcpp::CharacterVector(0);
    }

    // init kernel
    kernel_type = to_kernel_type(Rcpp::String(kernel));
    kernel_ptr  = make_kernel(kernel_type);
  }

  bool has_base_par() const { return base_par_idx >= 0; }
};


inline void bind_trend_op_to_paramtable(TrendOp& op,
                                        const ParamTable& pt)
{
  using Rcpp::as;
  using std::string;

  // Determine base type (if any)
  string base_type;
  if (op.spec.containsElementNamed("base")) {
    base_type = as<string>(op.spec["base"]);
  } else {
    base_type = "";   // no special base
  }

  const int n_tp = op.trend_pnames.size();

  op.base_par_idx = -1;
  op.kernel_par_indices.clear();

  int kernel_start = 0;

  if (op.base_type == "lin" || op.base_type == "lin_exp") {
    if (n_tp < 1) {
      Rcpp::stop("TrendOp '%s': base = '%s' but trend_pnames is empty",
                 op.target_param.c_str(), op.base_type.c_str());
    }
    string base_nm = as<string>(op.trend_pnames[0]);
    op.base_par_idx = pt.base_index_for(base_nm);
    kernel_start = 1;
  }

  for (int k = kernel_start; k < n_tp; ++k) {
    string kn = as<string>(op.trend_pnames[k]);
    int idx   = pt.base_index_for(kn);
    op.kernel_par_indices.push_back(idx);
  }
}

struct TrendEngine {
  Rcpp::List trend;
  Rcpp::DataFrame data;

  std::vector<TrendOp> premap_ops;
  std::vector<TrendOp> pretransform_ops;
  std::vector<TrendOp> posttransform_ops;

  std::unordered_set<std::string> premap_trend_params;

  TrendEngine(const Rcpp::List& trend_, const Rcpp::DataFrame& data_)
    : trend(trend_), data(data_) {
    initialize_from_trend();
  }

  void initialize_from_trend();

  bool has_premap() const;
  bool has_pretransform() const;
  bool has_posttransform() const;

  Rcpp::LogicalVector premap_design_mask(const Rcpp::List& designs) const;

  void bind_all_ops_to_paramtable(const ParamTable& pt);

  void run_kernel_for_op(TrendOp& op, ParamTable& pt);

  void apply_base_for_op(TrendOp& op,
                         ParamTable& pt);
  // e.g. helper that runs all premap kernels
  void run_premap_kernels(ParamTable& pt);
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

