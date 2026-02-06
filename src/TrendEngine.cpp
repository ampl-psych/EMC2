#include "TrendEngine.h"
#include <Rcpp.h>
using namespace Rcpp;

// helper function to initialize a covariate
static void init_covariate(TrendOpSpec& op, const Rcpp::DataFrame& data) {
  using namespace Rcpp;
  using std::string;

  if (!op.spec.containsElementNamed("covariate")) {
    stop("TrendOpSpec for target '%s' is missing 'covariate' field",
         op.target_param.c_str());
  }

  SEXP cov_spec = op.spec["covariate"];

  if (Rf_isNull(cov_spec)) {
    stop("TrendOpSpec for target '%s' has NULL 'covariate' but a non-empty covariate is required",
         op.target_param.c_str());
  }

  if (TYPEOF(cov_spec) != STRSXP) {
    stop("TrendOpSpec for target '%s': 'covariate' must be a character vector (names), not numeric",
         op.target_param.c_str());
  }

  CharacterVector cv(cov_spec);
  if (cv.size() != 1) {
    stop("TrendOpSpec for target '%s': covariate must be a single column name after expansion",
         op.target_param.c_str());
  }

  string cov_name = as<string>(cv[0]);
  if (!data.containsElementNamed(cov_name.c_str())) {
    stop("TrendPlan: data has no column named '%s' for covariate",
         cov_name.c_str());
  }

  op.kernel_input = data[cov_name];  // view/copy from data

  if (op.kernel_input.size() != data.nrows()) {
    stop("TrendPlan: covariate for target '%s' has length %d but data has %d rows",
         op.target_param.c_str(),
         op.kernel_input.size(), data.nrows());
  }
}

// helper function to initialize a par_input
// SM: dropped this. The idea was to bind a view to pt to the TrendOp,
// but that was a bit silly since it needs to be a runtime thing (parameters change every particle)
// static void init_par_input(TrendOpSpec& op, const ParamTable& pt) {
//   using namespace Rcpp;
//   using std::string;
//
//   if (!op.spec.containsElementNamed("par_input")) {
//     stop("TrendOpSpec for target '%s' is missing 'par_input' field",
//          op.target_param.c_str());
//   }
//
//   SEXP par_input_spec = op.spec["par_input"];
//
//   if (Rf_isNull(par_input_spec)) {
//     stop("TrendOpSpec for target '%s' has NULL 'par_input' but a non-empty par_input is required",
//          op.target_param.c_str());
//   }
//
//   if (TYPEOF(par_input_spec) != STRSXP) {
//     stop("TrendOpSpec for target '%s': 'par_input' must be a character vector (parameter names), not numeric",
//          op.target_param.c_str());
//   }
//
//   CharacterVector cv(par_input_spec);
//   if (cv.size() != 1) {
//     stop("TrendOpSpec for target '%s': par_input must be a single parameter name after expansion",
//          op.target_param.c_str());
//   }
//
//   string par_input_name = as<string>(cv[0]);
//
//   // Just validate; don't store a view any more
//   (void) pt.base_index_for(par_input_name); // throws if unknown
//
//   op.kernel_input = Rcpp::NumericVector(); // leave empty; we’ll get column at runtime
// }

static void init_covariate_maps_for_op(TrendOpSpec& op,
                                       const Rcpp::DataFrame& data) {
  using namespace Rcpp;
  using std::string;

  if (!op.spec.containsElementNamed("covariate_maps") ||
      Rf_isNull(op.spec["covariate_maps"])) {
    op.has_covariate_maps = false;
    op.covariate_map_cols.clear();
    return;
  }

  if (op.input_kind != InputKind::Covariate) {
    stop("TrendOpSpec for target '%s' has 'covariate_maps' but is not covariate-based",
         op.target_param.c_str());
  }

  SEXP cov_spec = op.spec["covariate"];
  if (!(Rf_isString(cov_spec) && Rf_length(cov_spec) == 1)) {
    stop("TrendOpSpec for target '%s': 'covariate' must be a single name when using covariate_maps",
         op.target_param.c_str());
  }
  string cov_name = as<string>(cov_spec);

  List maps(op.spec["covariate_maps"]);
  const int M = maps.size();
  if (M == 0) {
    op.has_covariate_maps = false;
    op.covariate_map_cols.clear();
    return;
  }

  const int T = data.nrows();
  op.covariate_map_cols.resize(M);

  for (int m = 0; m < M; ++m) {
    NumericMatrix mat(maps[m]);  // shallow wrapper around that element
    if (mat.nrow() != T) {
      stop("TrendOpSpec '%s': covariate_maps[[%d]] has %d rows, expected %d",
           op.target_param.c_str(), m + 1, mat.nrow(), T);
    }

    CharacterVector cn = colnames(mat);
    int col_idx = -1;
    for (int j = 0; j < cn.size(); ++j) {
      if (as<string>(cn[j]) == cov_name) {
        col_idx = j;
        break;
      }
    }
    if (col_idx < 0) {
      stop("TrendOpSpec '%s': covariate_maps[[%d]] has no column '%s'",
           op.target_param.c_str(), m + 1, cov_name.c_str());
    }

    // Store a column view for this map
    op.covariate_map_cols[m] = mat(_, col_idx);
  }

  op.has_covariate_maps = true;
}

// --- TrendOpSpec ---

TrendOpSpec::TrendOpSpec(const Rcpp::List& cur, const std::string& name_from_list)
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

  // if (cur.containsElementNamed("par_input")) {
  //   par_input = cur["par_input"];
  // } else {
  //   par_input = Rcpp::CharacterVector(0);
  // }

  // Handle "at" being possibly NULL
  if (cur.containsElementNamed("at") && !Rf_isNull(cur["at"])) {
    has_at = true;
    at     = as<std::string>(cur["at"]);
  } else {
    has_at = false;
  }

  // init kernel type
  kernel_type = to_kernel_type(Rcpp::String(kernel));
}

void TrendOpSpec::make_first_level(const Rcpp::DataFrame& data) {
  using namespace Rcpp;

  const int n = data.nrows();
  if (n <= 0) {
    stop("make_first_level: data has zero rows");
  }

  // Default: no 'at' => use all rows
  first_level = LogicalVector(n, true);

  if (has_at) {
    if (!data.containsElementNamed(at.c_str())) {
      stop("make_first_level: data has no column named '%s' for 'at'",
           at.c_str());
    }

    SEXP at_col = data[at];
    if (!Rf_inherits(at_col, "factor")) {
      stop("'at' column '%s' must be a factor", at.c_str());
    }

    IntegerVector f(at_col);
    if (f.size() != n) {
      stop("'at' column '%s' has wrong length (%d != %d)",
           at.c_str(), f.size(), n);
    }

    first_level = (f == 1);
  }

  // build expand_idx as cumulative count of TRUEs
  expand_idx.assign(n, 0);
  int count = 0;
  for (int i = 0; i < n; ++i) {
    if (first_level[i]) ++count;
    expand_idx[i] = count;
  }

  if (count == 0) {
    stop("make_first_level: no rows with first 'at' level (== 1) found");
  }

  // Optional consistency check from your original code
  for (int i = 0; i < n; ++i) {
    if (expand_idx[i] == 0) {
      stop("Found rows before first 'at' level within subject. "
             "Cannot anchor expansion.");
    }
  }

  // comp_index = positions where first_level is TRUE
  comp_index.clear();
  comp_index.reserve(count);
  for (int i = 0; i < n; ++i) {
    if (first_level[i]) comp_index.push_back(i);
  }

  if ((int)comp_index.size() != count) {
    stop("make_first_level: comp_index size (%d) != last expand_idx (%d)",
         (int)comp_index.size(), count);
  }
}

// helper
extern std::vector<std::string> covariate_names_from_spec(const Rcpp::List& tr);
extern std::vector<std::string> par_input_names_from_spec(const Rcpp::List& tr);

TrendPlan::TrendPlan(Rcpp::Nullable<Rcpp::List> trend_,
                     const Rcpp::DataFrame& data_)
  : data(data_) {

  if (trend_.isNull()) {
    trend = R_NilValue;
    return;
  }

  trend = Rcpp::List(trend_.get());

  using std::string;
  using namespace Rcpp;

  // helper functions to figure out whether there's covariate_inputs and/or par_inputs
  auto uses_cov_input = [](const Rcpp::List& tr) -> bool {
    if (!tr.containsElementNamed("covariate")) return false;

    SEXP x = tr["covariate"];
    if (Rf_isNull(x)) return false;

    if (TYPEOF(x) == STRSXP) {
      // character vector: length == 0 => no covariate
      return Rf_length(x) > 0;
    }

    // non-string (e.g., numeric vector) => inline covariate
    return true;
  };

  auto uses_par_input = [](const Rcpp::List& tr) -> bool {
    if (!tr.containsElementNamed("par_input")) return false;

    SEXP x = tr["par_input"];
    if (Rf_isNull(x)) return false;

    if (TYPEOF(x) == STRSXP) {
      // empty character vector => no par_input
      return Rf_length(x) > 0;
    }

    // non-string (e.g., numeric vector) => inline par_input
    return true;
  };

  CharacterVector tnames = trend.names();
  if (tnames.size() != trend.size()) {
    stop("TrendPlan: 'trend' must be a named list");
  }

  premap_ops.clear();
  pretransform_ops.clear();
  posttransform_ops.clear();
  premap_trend_params.clear();
  pretransform_trend_params.clear();
  posttransform_trend_params.clear();
  all_trend_params.clear();

  // Loopy over trends
  for (int i = 0; i < trend.size(); ++i) {
    if (trend[i] == R_NilValue) continue;

    Rcpp::List tr_i = trend[i];
    std::string name_i = Rcpp::as<std::string>(tnames[i]);

    // Build once to inspect phase, trend_pnames, etc.
    TrendOpSpec tmp(tr_i, name_i);

    // --- collect trend_pnames into sets -- useful for transformations later on ---
    for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
      std::string pn = Rcpp::as<std::string>(tmp.trend_pnames[k]);
      all_trend_params.insert(pn);
    }

    if (tmp.phase == TrendPhase::Premap) {
      for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
        premap_trend_params.insert(Rcpp::as<std::string>(tmp.trend_pnames[k]));
      }
    }
    if (tmp.phase == TrendPhase::Pretransform) {
      for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
        pretransform_trend_params.insert(Rcpp::as<std::string>(tmp.trend_pnames[k]));
      }
    }
    if (tmp.phase == TrendPhase::Posttransform) {
      for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
        posttransform_trend_params.insert(Rcpp::as<std::string>(tmp.trend_pnames[k]));
      }
    }

    // ---- decide which inputs are actually used in the current trend ----
    bool has_cov_input = uses_cov_input(tr_i);
    bool has_par_input = uses_par_input(tr_i);

    if (!has_cov_input && !has_par_input) {
      Rcpp::stop("TrendPlan: spec for '%s' must have non-empty 'covariate' or 'par_input'",
                 name_i.c_str());
    }

    // Extract names only if that input type is used
    std::vector<std::string> cov_names;
    std::vector<std::string> par_input_names;
    if (has_cov_input) {
      cov_names = covariate_names_from_spec(tr_i);
    }
    if (has_par_input) {
      par_input_names = par_input_names_from_spec(tr_i);
    }

    // Helper: add an op to the correct phase vector and init input
    auto add_op_for_phase = [&](TrendPhase phase,
                                const Rcpp::List& spec_list,
                                bool use_cov_input,
                                bool use_par_input) {
      switch (phase) {
      case TrendPhase::Premap:
        premap_ops.emplace_back(spec_list, name_i);
        if (use_cov_input) {
          init_covariate(premap_ops.back(), data);
          init_covariate_maps_for_op(premap_ops.back(), data);
          premap_ops.back().input_kind = InputKind::Covariate;
        } else if (use_par_input) {
          // par_input: record name
          Rcpp::CharacterVector cv(spec_list["par_input"]);
          if (cv.size() != 1) {
            Rcpp::stop("TrendOpSpec '%s': par_input must be single name after expansion",
                       name_i.c_str());
          }
          premap_ops.back().input_kind = InputKind::ParInput;
          premap_ops.back().par_input_name = Rcpp::as<std::string>(cv[0]);
        } else {
          premap_ops.back().input_kind = InputKind::None;
        }
        break;

      case TrendPhase::Pretransform:
        pretransform_ops.emplace_back(spec_list, name_i);
        if (use_cov_input) {
          init_covariate(pretransform_ops.back(), data);
          init_covariate_maps_for_op(pretransform_ops.back(), data);
          pretransform_ops.back().input_kind = InputKind::Covariate;
        } else if (use_par_input) {
          // par_input: record name
          Rcpp::CharacterVector cv(spec_list["par_input"]);
          if (cv.size() != 1) {
            Rcpp::stop("TrendOpSpec '%s': par_input must be single name after expansion",
                       name_i.c_str());
          }
          pretransform_ops.back().input_kind = InputKind::ParInput;
          pretransform_ops.back().par_input_name = Rcpp::as<std::string>(cv[0]);
        } else {
          pretransform_ops.back().input_kind = InputKind::None;
        }
        break;

      case TrendPhase::Posttransform:
        posttransform_ops.emplace_back(spec_list, name_i);
        if (use_cov_input) {
          init_covariate(posttransform_ops.back(), data);
          init_covariate_maps_for_op(posttransform_ops.back(), data);
          posttransform_ops.back().input_kind = InputKind::Covariate;
        } else if (use_par_input) {
          // par_input: record name
          Rcpp::CharacterVector cv(spec_list["par_input"]);
          if (cv.size() != 1) {
            Rcpp::stop("TrendOpSpec '%s': par_input must be single name after expansion",
                       name_i.c_str());
          }
          posttransform_ops.back().input_kind = InputKind::ParInput;
          posttransform_ops.back().par_input_name = Rcpp::as<std::string>(cv[0]);
        } else {
          posttransform_ops.back().input_kind = InputKind::None;
        }
        break;
      }
    };

    // ---- 1) Covariate-based ops (0, 1, or many) ----
    if (has_cov_input) {
      if (cov_names.empty()) {
        // inline numeric covariate
        add_op_for_phase(tmp.phase, tr_i, /*use_cov_input=*/true, false);
      } else if (cov_names.size() == 1) {
        // single named column
        add_op_for_phase(tmp.phase, tr_i, /*use_cov_input=*/true, false);
      } else {
        // multiple covariate columns: clone spec per column
        for (const std::string& cov_name : cov_names) {
          Rcpp::List tr_copy = Rcpp::clone(tr_i);
          tr_copy["covariate"] = cov_name;
          add_op_for_phase(tmp.phase, tr_copy, /*use_cov_input=*/true, false);
        }
      }
    }

    // ---- 2) par_input-based ops (0, 1, or many) ----
    if (has_par_input) {
      if (par_input_names.empty()) {
        // inline numeric par_input
        add_op_for_phase(tmp.phase, tr_i, /*use_cov_input=*/false, true);
      } else if (par_input_names.size() == 1) {
        // single parameter name
        add_op_for_phase(tmp.phase, tr_i, /*use_cov_input=*/false, true);
      } else {
        // multiple par_input names: clone spec per parameter
        for (const std::string& par_name : par_input_names) {
          Rcpp::List tr_copy = Rcpp::clone(tr_i);
          tr_copy["par_input"] = par_name;
          add_op_for_phase(tmp.phase, tr_copy, /*use_cov_input=*/false, true);
        }
      }
    }
  }

  // Build first-level masks & comp_index/expand_idx
  for (auto& op : premap_ops)        op.make_first_level(data);
  for (auto& op : pretransform_ops)  op.make_first_level(data);
  for (auto& op : posttransform_ops) op.make_first_level(data);

}


Rcpp::LogicalVector TrendPlan::premap_design_mask(const Rcpp::List& designs) const {
  using namespace Rcpp;
  using std::string;

  const int n = designs.size();
  LogicalVector mask(n, false);

  CharacterVector dnames = designs.names();
  if (dnames.size() != n) {
    stop("TrendPlan::premap_design_mask: 'designs' must be a named list");
  }

  for (int i = 0; i < n; ++i) {
    string nm = as<string>(dnames[i]);
    mask[i] = (premap_trend_params.find(nm) != premap_trend_params.end());
  }

  return mask;
}


inline void bind_trend_op_to_paramtable(TrendOpRuntime& op,
                                        const ParamTable& pt)
{
  const TrendOpSpec& spec = *op.spec;

  const int n_tp = spec.trend_pnames.size();
  op.base_par_idx = -1;
  op.base_par_indices.clear();
  op.kernel_par_indices.clear();

  int kernel_start = 0;

  bool needs_base_par =
    (spec.base_type == "lin" || spec.base_type == "exp_lin" ||
    spec.base_type == "lin_exp" || spec.base_type == "centered");

  if (needs_base_par) {
    int n_base_pars = 0;
    if (spec.has_covariate_maps) {
      n_base_pars = static_cast<int>(spec.covariate_map_cols.size());  // one per map
    } else {
      n_base_pars = 1;
    }

    if (n_tp < n_base_pars) {
      Rcpp::stop("TrendOp '%s': base_type '%s' expects at least %d trend_pnames (got %d)",
                 spec.target_param.c_str(), spec.base_type.c_str(), n_base_pars, n_tp);
    }

    // First n_base_pars trend_pnames are base parameters
    for (int b = 0; b < n_base_pars; ++b) {
      std::string base_nm = Rcpp::as<std::string>(spec.trend_pnames[b]);
      int idx = pt.base_index_for(base_nm);
      op.base_par_indices.push_back(idx);
    }
    op.base_par_idx = op.base_par_indices[0];
    kernel_start = n_base_pars;
  }

  // Remaining trend_pnames are kernel parameters
  for (int k = kernel_start; k < n_tp; ++k) {
    std::string kn = Rcpp::as<std::string>(spec.trend_pnames[k]);
    int idx = pt.base_index_for(kn);
    op.kernel_par_indices.push_back(idx);
  }
  // Rprintf("Binding...");
}



TrendRuntime::TrendRuntime(const TrendPlan& plan_)
  : plan(&plan_) {

  premap_ops.reserve(plan->premap_ops.size());
  for (const auto& spec : plan->premap_ops) {
    premap_ops.emplace_back(&spec);
  }

  pretransform_ops.reserve(plan->pretransform_ops.size());
  for (const auto& spec : plan->pretransform_ops) {
    pretransform_ops.emplace_back(&spec);
  }

  posttransform_ops.reserve(plan->posttransform_ops.size());
  for (const auto& spec : plan->posttransform_ops) {
    posttransform_ops.emplace_back(&spec);
  }
}

void TrendRuntime::bind_all_ops_to_paramtable(const ParamTable& pt) {
  for (auto& op : premap_ops)        bind_trend_op_to_paramtable(op, pt);
  for (auto& op : pretransform_ops)  bind_trend_op_to_paramtable(op, pt);
  for (auto& op : posttransform_ops) bind_trend_op_to_paramtable(op, pt);

  // Validate par_inputs --
  auto validate_par_inputs = [&](const TrendOpRuntime& op) {
    const TrendOpSpec& spec = *op.spec;
    if (spec.spec.containsElementNamed("par_input") &&
        !Rf_isNull(spec.spec["par_input"])) {
          SEXP par_spec = spec.spec["par_input"];
          if (TYPEOF(par_spec) == STRSXP) {
            Rcpp::CharacterVector cv(par_spec);
            if (cv.size() == 1) {
              std::string par_name = Rcpp::as<std::string>(cv[0]);
              (void) pt.base_index_for(par_name); // throws if unknown
            }
        }
    }
  };

  for (auto& op : premap_ops)        validate_par_inputs(op);
  for (auto& op : pretransform_ops)  validate_par_inputs(op);
  for (auto& op : posttransform_ops) validate_par_inputs(op);
}

void TrendRuntime::run_kernel_for_op(TrendOpRuntime& op,
                                     ParamTable& pt)
{
  using namespace Rcpp;
  using std::string;

  const TrendOpSpec& spec = *op.spec;

  NumericVector input;

  // waht type of input?
  switch (spec.input_kind) {
  case InputKind::ParInput:
    // Resolve from ParamTable every time (current particle)
    input = pt.column_by_name(spec.par_input_name);
    break;

  case InputKind::Covariate:
    // Use cached data column
    input = spec.kernel_input;
    break;

  case InputKind::None:
  default:
    Rcpp::stop("TrendOp for target '%s' has neither covariate nor par_input input",
               spec.target_param.c_str());
  }

  const int n = input.size();
  if (n != pt.n_trials) {
    Rcpp::stop("TrendRuntime::run_kernel_for_op: input length (%d) != n_trials (%d)",
               n, pt.n_trials);
  }

  const int p = pt.base.ncol();
  for (int idx : op.kernel_par_indices) {
    if (idx < 0 || idx >= p) {
      stop("TrendRuntime::run_kernel_for_op('%s'): kernel_par index %d out of [0,%d)",
           spec.target_param.c_str(), idx, p);
    }
  }

  KernelParsView kp_view = make_kernel_pars_view(pt, op.kernel_par_indices);

  op.kernel_ptr->reset();
  op.kernel_ptr->run(kp_view, input, spec.comp_index);

  if (spec.has_at) {
    op.kernel_ptr->do_expand(spec.expand_idx);
  }
}

void TrendRuntime::apply_base_for_op(TrendOpRuntime& op,
                                     ParamTable& pt)
{
  using namespace Rcpp;

  const TrendOpSpec& spec = *op.spec;

  // Ensure kernel has been run
  if (!op.kernel_ptr->has_run()) {
    run_kernel_for_op(op, pt);
  }

  const std::vector<double>& traj = op.kernel_ptr->get_output();
  const int n = pt.n_trials;
  if ((int)traj.size() != n) {
    stop("TrendRuntime::apply_base_for_op('%s'): trajectory length (%d) != n_trials (%d)",
         spec.target_param.c_str(), (int)traj.size(), n);
  }

  // // SM TEMP
  // // DEBUG: print a few trajectory values for pretransform ops
  // if (spec.phase == TrendPhase::Pretransform) {
  //   Rcout << "apply_base_for_op[pretransform] target=" << spec.target_param
  //         << " traj[0..4]=";
  //   for (int k = 0; k < std::min(n, 5); ++k) {
  //     Rcout << traj[k] << " ";
  //   }
  //   Rcout << "\n";
  // }

  const std::string& base = spec.base_type;

  // Target parameter column
  int target_idx = pt.base_index_for(spec.target_param);
  double* target_col = &pt.base(0, target_idx);

  // Base/trend parameter if needed
  // const double* trend_col = nullptr;
  bool needs_trend_par =
    (base == "lin" || base == "exp_lin" || base == "lin_exp" || base == "centered");

  const double* single_trend_col = nullptr;
  if (needs_trend_par && !spec.has_covariate_maps) {
    if (op.base_par_idx < 0) {
      stop("TrendOp '%s': base_type = '%s' but base_par_idx < 0",
           spec.target_param.c_str(), base.c_str());
    }
    single_trend_col = &pt.base(0, op.base_par_idx);
  }

  for (int r = 0; r < n; ++r) {
    double p = target_col[r];
    if (NumericVector::is_na(p) || std::isnan(p)) {
      continue;
    }

    double q = traj[r];
    if (NumericVector::is_na(q) || std::isnan(q)) {
      q = 0.0;
    }

    double contrib = 0.0;

    if (spec.has_covariate_maps) {
      const int K = static_cast<int>(spec.covariate_map_cols.size());
      if (K != (int)op.base_par_indices.size()) {
        stop("TrendOp '%s': mismatch between covariate_maps (%d) and base_pars (%d)",
             spec.target_param.c_str(), K, (int)op.base_par_indices.size());
      }

      for (int k = 0; k < K; ++k) {
        double base_val = pt.base(r, op.base_par_indices[k]);
        double map_val  = spec.covariate_map_cols[k][r];
        contrib += q * base_val * map_val;
      }
    } else {
      // Single base parameter behaviour
      double tmp = q;
      if (base == "lin" || base == "exp_lin" || base == "lin_exp") {
        double tp = single_trend_col ? single_trend_col[r] : 1.0;
        tmp *= tp;
      } else if (base == "centered") {
        double tp = single_trend_col ? single_trend_col[r] : 1.0;
        tmp = (tmp - 0.5) * tp;
      }
      contrib = tmp;
    }

    double base_term;
    if (base == "exp_lin" || base == "lin_exp") {
      base_term = std::exp(p);
    } else if (base == "lin" || base == "centered" || base == "add" || base.empty()) {
      base_term = p;
    } else {
      base_term = p;
    }

    target_col[r] = contrib + base_term;
  }
}

void TrendRuntime::reset_all_kernels() {
  for (auto& op : premap_ops)        op.kernel_ptr->reset();
  for (auto& op : pretransform_ops)  op.kernel_ptr->reset();
  for (auto& op : posttransform_ops) op.kernel_ptr->reset();
}


Rcpp::NumericMatrix TrendRuntime::all_kernel_outputs(ParamTable& pt) {
  using namespace Rcpp;

  const int n = pt.n_trials;
  const int n_ops =
    static_cast<int>(premap_ops.size() +
    pretransform_ops.size() +
    posttransform_ops.size());

  NumericMatrix out(n, n_ops);
  CharacterVector cn(n_ops);

  int col = 0;

  auto fill_for_ops = [&](std::vector<TrendOpRuntime>& ops, const char* phase_label) {
    for (size_t i = 0; i < ops.size(); ++i) {
      TrendOpRuntime& op = ops[i];
      const TrendOpSpec& spec = *op.spec;

      // Ensure kernel has been run for this ParamTable
      // run_kernel_for_op(op, pt);

      const std::vector<double>& traj = op.kernel_ptr->get_output();
      if (static_cast<int>(traj.size()) != n) {
        stop("TrendRuntime::all_kernel_outputs('%s'): trajectory length (%d) "
               "!= n_trials (%d)",
               op.spec->target_param.c_str(), (int)traj.size(), n);
      }

      // Copy trajectory into column 'col'
      for (int r = 0; r < n; ++r) {
        out(r, col) = traj[r];
      }

      // --- build name: target_param.covariate_name/par_input_name ---
      std::string cov_name;

      if (spec.spec.containsElementNamed("covariate") &&
          !Rf_isNull(spec.spec["covariate"])) {
          SEXP cov_spec = spec.spec["covariate"];
        if (Rf_isString(cov_spec) && Rf_length(cov_spec) == 1) {
          cov_name = Rcpp::as<std::string>(cov_spec);
        } else {
          cov_name = "cov_inline";
        }
      } else if (spec.spec.containsElementNamed("par_input") &&
        !Rf_isNull(spec.spec["par_input"])) {
        SEXP par_spec = spec.spec["par_input"];
        if (Rf_isString(par_spec) && Rf_length(par_spec) == 1) {
          cov_name = Rcpp::as<std::string>(par_spec);
        } else {
          cov_name = "par_inline";
        }
      } else {
        cov_name = "noinput";
      }

      std::string cname = spec.target_param + "." + cov_name;
      cn[col] = cname;

      ++col;
    }
  };

  fill_for_ops(premap_ops,        "premap");
  fill_for_ops(pretransform_ops,  "pretransform");
  fill_for_ops(posttransform_ops, "posttransform");

  colnames(out) = cn;
  return out;
}
