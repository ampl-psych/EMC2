#include "TrendEngine.h"
#include <Rcpp.h>
using namespace Rcpp;


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

  // Optional consistency check (as in your original code)
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


// helpers
extern std::vector<std::string> covariate_names_from_spec(const Rcpp::List& tr);
extern std::vector<std::string> par_input_names_from_spec(const Rcpp::List& tr);

bool uses_cov_input(const Rcpp::List& tr) {
  if (!tr.containsElementNamed("covariate")) return false;
  SEXP x = tr["covariate"];
  if (Rf_isNull(x)) return false;
  if (TYPEOF(x) == STRSXP) {
    // character vector: length > 0 => use covariate
    return Rf_length(x) > 0;
  }
  // non-string (e.g., numeric vector) => inline covariate
  return true;
}

bool uses_par_input(const Rcpp::List& tr) {
  if (!tr.containsElementNamed("par_input")) return false;
  SEXP x = tr["par_input"];
  if (Rf_isNull(x)) return false;
  if (TYPEOF(x) == STRSXP) {
    return Rf_length(x) > 0;
  }
  // non-string (e.g., numeric vector) => inline par_input
  return true;
}

// Initialize covariate for a single kernel slot when covariate is a *single* name
void init_covariate_for_slot(KernelSlotSpec& slot,
                                    const Rcpp::List& spec,
                                    const Rcpp::DataFrame& data) {
  using namespace Rcpp;
  if (!spec.containsElementNamed("covariate")) {
    stop("init_covariate_for_slot: missing 'covariate' field");
  }

  SEXP cov_spec = spec["covariate"];
  if (!Rf_isString(cov_spec) || Rf_length(cov_spec) != 1) {
    stop("init_covariate_for_slot: 'covariate' must be a single column name");
  }

  std::string cov_name = as<std::string>(STRING_ELT(cov_spec, 0));
  if (!data.containsElementNamed(cov_name.c_str())) {
    stop("init_covariate_for_slot: data has no column named '%s'", cov_name.c_str());
  }

  slot.kernel_input = data[cov_name];  // shallow
  if (slot.kernel_input.size() != data.nrows()) {
    stop("init_covariate_for_slot: covariate '%s' length %d != data rows %d",
         cov_name.c_str(), slot.kernel_input.size(), data.nrows());
  }
}


void init_covariate_maps_for_slot(KernelSlotSpec& slot,
                                  const Rcpp::List& spec,    // spec with single 'covariate'
                                  const Rcpp::DataFrame& data,
                                  const Rcpp::List& data_covmaps) {
  using namespace Rcpp;
  using std::string;

  // No map field or NULL: no maps for this slot
  if (!spec.containsElementNamed("map") || Rf_isNull(spec["map"])) {
    slot.has_covariate_maps = false;
    slot.covariate_map_cols.clear();
    return;
  }

  if (data_covmaps.size() == 0) {
    stop("Trend: 'map' specified but data has no 'covariate_maps' attribute");
  }

  if (!spec.containsElementNamed("covariate") ||
      TYPEOF(spec["covariate"]) != STRSXP ||
      Rf_length(spec["covariate"]) != 1) {
    stop("Trend: 'map' requires a single covariate name per kernel slot");
  }

  string cov_name = as<string>(STRING_ELT(spec["covariate"], 0));

  List maps_spec(spec["map"]);
  CharacterVector map_names = maps_spec.names();
  if (map_names.size() != maps_spec.size()) {
    stop("Trend: 'map' list must be named");
  }

  List maps_data(data_covmaps);
  CharacterVector data_map_names = maps_data.names();

  const int M = maps_spec.size();
  slot.covariate_map_cols.resize(M);

  const int T = data.nrows();

  for (int m = 0; m < M; ++m) {
    string map_nm = as<string>(map_names[m]);

    // find numeric map matrix in data_covmaps
    int idx_data = -1;
    for (int j = 0; j < data_map_names.size(); ++j) {
      if (as<string>(data_map_names[j]) == map_nm) {
        idx_data = j;
        break;
      }
    }
    if (idx_data < 0) {
      stop("Trend: covariate_map '%s' not found in data's 'covariate_maps' attribute",
           map_nm.c_str());
    }

    NumericMatrix mat(maps_data[idx_data]);
    if (mat.nrow() != T) {
      stop("Trend: covariate_maps[['%s']] has %d rows, expected %d",
           map_nm.c_str(), mat.nrow(), T);
    }

    // Find column for this covariate
    CharacterVector cn = colnames(mat);
    int col_idx = -1;
    for (int j = 0; j < cn.size(); ++j) {
      if (as<string>(cn[j]) == cov_name) {
        col_idx = j;
        break;
      }
    }
    if (col_idx < 0) {
      stop("Trend: covariate_maps[['%s']] has no column '%s'",
           map_nm.c_str(), cov_name.c_str());
    }

    slot.covariate_map_cols[m] = mat(_, col_idx); // view column
  }

  slot.has_covariate_maps = true;
}


TrendPlan::TrendPlan(Rcpp::Nullable<Rcpp::List> trend_,
                     const Rcpp::DataFrame& data_)
  : data(data_) {

  using namespace Rcpp;
  using std::string;

  if (trend_.isNull()) {
    trend = R_NilValue;
    return;
  }

  trend = Rcpp::List(trend_.get());

  // covariate_maps from data attribute (if any)
  if (data_.hasAttribute("covariate_maps")) {
    SEXP cm = data_.attr("covariate_maps");
    if (!Rf_isNull(cm)) {
      data_covariate_maps = Rcpp::List(cm);
      has_data_covariate_maps = data_covariate_maps.size() > 0;
      // Rprintf("Found a covariate map?!");
    }
  }

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

  const int n_specs = trend.size();

  // loopy loop
  for (int i = 0; i < n_specs; ++i) {
    if (trend[i] == R_NilValue) continue;

    Rcpp::List tr_i = trend[i];
    std::string name_i = Rcpp::as<std::string>(tnames[i]);

    if (!tr_i.containsElementNamed("phase") ||
        !tr_i.containsElementNamed("kernel") ||
        !tr_i.containsElementNamed("base")) {
        stop("TrendPlan: spec '%s' must have 'phase', 'kernel', and 'base'", name_i.c_str());
    }

    TrendPhase phase = parse_phase(Rcpp::as<std::string>(tr_i["phase"]));
    std::string kernel_name = Rcpp::as<std::string>(tr_i["kernel"]);
    std::string base_name   = Rcpp::as<std::string>(tr_i["base"]);

    // Build one TrendOpSpec for this entry. The TrendOpSpec contains one base and a group of kernels
    TrendOpSpec op;
    op.phase        = phase;
    op.target_param = name_i;
    op.base_type    = base_name;

    // 'at' is shared across all kernels in this op
    op.has_at = tr_i.containsElementNamed("at") && !Rf_isNull(tr_i["at"]);
    if (op.has_at) {
      op.at = Rcpp::as<std::string>(tr_i["at"]);
    }

    // trend_pnames from this spec
    // std::unordered_set<std::string> tp_set;
    op.trend_pnames = Rcpp::CharacterVector(0);
    if (tr_i.containsElementNamed("trend_pnames")) {
      Rcpp::CharacterVector tp(tr_i["trend_pnames"]);
      op.trend_pnames = tp;
      for (int k = 0; k < tp.size(); ++k) {
        std::string pn = Rcpp::as<std::string>(tp[k]);
        all_trend_params.insert(pn);
        if (phase == TrendPhase::Premap)        premap_trend_params.insert(pn);
        if (phase == TrendPhase::Pretransform)  pretransform_trend_params.insert(pn);
        if (phase == TrendPhase::Posttransform) posttransform_trend_params.insert(pn);
      }
    }

    // Which kinds of input does this op use?
    const bool has_cov = uses_cov_input(tr_i);
    const bool has_par = uses_par_input(tr_i);

    if (!has_cov && !has_par) {
      stop("TrendPlan: spec '%s' must have non-empty 'covariate' or 'par_input'",
           name_i.c_str());
    }

    KernelType ktype   = to_kernel_type(Rcpp::String(kernel_name));
    KernelMeta kmeta   = kernel_meta(ktype);

    // ---- Covariate-based kernel slots ----
    if (has_cov) {
      SEXP cov_spec = tr_i["covariate"];

      if (TYPEOF(cov_spec) == STRSXP) {
        // character vector: one KernelSlot per covariate
        Rcpp::CharacterVector cv(cov_spec);
        for (int j = 0; j < cv.size(); ++j) {
          std::string cov_name = Rcpp::as<std::string>(cv[j]);

          if (!(kmeta.supports_grouping && kmeta.input_arity == 1)) {
            Rcpp::stop("Kernel '%s' does not support grouped covariates", kernel_name.c_str());
          }

          KernelSlotSpec slot;
          slot.kernel      = kernel_name;
          slot.kernel_type = ktype;
          slot.input_kind  = InputKind::Covariate;

          // spec copy with *single* covariate name
          Rcpp::List tr_copy = Rcpp::clone(tr_i);
          tr_copy["covariate"] = cov_name;

          init_covariate_for_slot(slot, tr_copy, data);

          // Attach maps if present on this spec and in data
          if (has_data_covariate_maps) {
            init_covariate_maps_for_slot(slot, tr_copy, data, data_covariate_maps);
          } else {
            // If the spec has a non-NULL 'map', but data has no maps, error
            if (tr_copy.containsElementNamed("map") && !Rf_isNull(tr_copy["map"])) {
              Rcpp::stop("TrendPlan: spec '%s' has 'map' but data has no 'covariate_maps' attribute",
                         name_i.c_str());
            }
          }

          op.kernels.push_back(std::move(slot));
        }
      }
    }


    // ---- par_input-based kernel slots ----
    if (has_par) {
      SEXP par_spec = tr_i["par_input"];

      if (TYPEOF(par_spec) == STRSXP) {
        CharacterVector pv(par_spec);
        for (int j = 0; j < pv.size(); ++j) {
          std::string par_nm = Rcpp::as<std::string>(pv[j]);

          if (!(kmeta.supports_grouping && kmeta.input_arity == 1)) {
            stop("Kernel '%s' does not support grouped par_input", kernel_name.c_str());
          }

          KernelSlotSpec slot;
          slot.kernel        = kernel_name;
          slot.kernel_type   = ktype;
          slot.input_kind    = InputKind::ParInput;
          slot.par_input_name = par_nm;

          op.kernels.push_back(std::move(slot));
        }
      } else {
        // inline par_input (numeric) not implemented
        stop("TrendPlan: inline numeric 'par_input' not yet supported for '%s'",
             name_i.c_str());
      }
    }

    // Store op in phase-specific vector
    if (phase == TrendPhase::Premap)        premap_ops.push_back(std::move(op));
    else if (phase == TrendPhase::Pretransform)  pretransform_ops.push_back(std::move(op));
    else if (phase == TrendPhase::Posttransform) posttransform_ops.push_back(std::move(op));
  }

  // Build first_level / comp_index / expand_idx once per op
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




TrendRuntime::TrendRuntime(const TrendPlan& plan_)
  : plan(&plan_) {

  premap_ops.reserve(plan->premap_ops.size());
  for (const auto& spec : plan->premap_ops) {
    TrendOpRuntime rt;
    rt.spec = &spec;
    premap_ops.push_back(std::move(rt));
  }

  pretransform_ops.reserve(plan->pretransform_ops.size());
  for (const auto& spec : plan->pretransform_ops) {
    TrendOpRuntime rt;
    rt.spec = &spec;
    pretransform_ops.push_back(std::move(rt));
  }

  posttransform_ops.reserve(plan->posttransform_ops.size());
  for (const auto& spec : plan->posttransform_ops) {
    TrendOpRuntime rt;
    rt.spec = &spec;
    posttransform_ops.push_back(std::move(rt));
  }
}

// Updated to allow for kernelgroups
void TrendRuntime::bind_all_ops_to_paramtable(const ParamTable& pt) {
  auto bind_vec = [&](std::vector<TrendOpRuntime>& rt_ops,
                      const std::vector<TrendOpSpec>& sp_ops) {
    rt_ops.clear();
    rt_ops.reserve(sp_ops.size());

    for (const auto& spec : sp_ops) {
      TrendOpRuntime op_rt;
      op_rt.spec = &spec;

      // ---- 1) Bind base parameters ----
      const int n_tp = spec.trend_pnames.size();
      const std::string& base = spec.base_type;
      // Rprintf("Identified base = %s\n", base.c_str());

      bool needs_base_par = (base == "lin" || base == "exp_lin" || base == "lin_exp" || base == "centered");

      int n_base_pars = 0;
      bool slot_has_maps = !spec.kernels.empty() && spec.kernels[0].has_covariate_maps;
      if (needs_base_par) {
        if (slot_has_maps) {
          n_base_pars = static_cast<int>(spec.kernels[0].covariate_map_cols.size());
        } else {
          n_base_pars = 1;
        }
      }

      if (n_tp < n_base_pars) {
        Rcpp::stop("TrendOp '%s': base_type '%s' expects at least %d trend_pnames (got %d)",
                   spec.target_param.c_str(), base.c_str(), n_base_pars, n_tp);
      }

      op_rt.base_par_indices.clear();
      for (int b = 0; b < n_base_pars; ++b) {
        std::string base_nm = Rcpp::as<std::string>(spec.trend_pnames[b]);
        // Rprintf("base_nm: %s\n", base_nm.c_str());
        int idx = pt.base_index_for(base_nm);
        op_rt.base_par_indices.push_back(idx);
      }
      op_rt.base_par_idx = op_rt.base_par_indices.empty() ? -1 : op_rt.base_par_indices[0];

      // Rprintf("n_base_pars = %d, base_par_indices[0] = %d\n", n_base_pars, op_rt.base_par_idx);


      // ---- 2) Kernel parameter indices (shared across all kernels in this TrendOp) ----
      std::vector<int> kernel_indices;
      kernel_indices.reserve(std::max(0, n_tp - n_base_pars));
      for (int k = n_base_pars; k < n_tp; ++k) {
        std::string kn = Rcpp::as<std::string>(spec.trend_pnames[k]);
        // Rprintf("Identified kernel parameter name: %s\n", kn.c_str());
        int idx = pt.base_index_for(kn);
        kernel_indices.push_back(idx);
      }

      // ---- 3) Build KernelSlotRuntime for each kernel spec ----
      op_rt.kernels.clear();
      op_rt.kernels.reserve(spec.kernels.size());
      for (const auto& kspec : spec.kernels) {
        KernelSlotRuntime k_rt;
        k_rt.spec         = &kspec;
        k_rt.kernel_ptr   = make_kernel(kspec.kernel_type);
        k_rt.kernel_par_indices = kernel_indices;
        op_rt.kernels.push_back(std::move(k_rt));
      }

      rt_ops.push_back(std::move(op_rt));
    }
  };

  bind_vec(premap_ops,        plan->premap_ops);
  bind_vec(pretransform_ops,  plan->pretransform_ops);
  bind_vec(posttransform_ops, plan->posttransform_ops);

  // Optionally: validate par_input names for all kernels
  auto validate_par_inputs = [&](const TrendOpRuntime& op) {
    for (const auto& k_rt : op.kernels) {
      const KernelSlotSpec& kspec = *k_rt.spec;
      if (kspec.input_kind == InputKind::ParInput) {
        (void) pt.base_index_for(kspec.par_input_name); // throws if unknown
      }
    }
  };

  for (const auto& op : premap_ops)        validate_par_inputs(op);
  for (const auto& op : pretransform_ops)  validate_par_inputs(op);
  for (const auto& op : posttransform_ops) validate_par_inputs(op);
}


void TrendRuntime::run_kernels_for_op(TrendOpRuntime& op,
                                      ParamTable& pt)
{
  using namespace Rcpp;

  const TrendOpSpec& spec = *op.spec;
  const int n = pt.n_trials;

  // since kernels are now kernelgroups, loop over the kernels
  for (auto& k_rt : op.kernels) {
    const KernelSlotSpec& kspec = *k_rt.spec;

    // Select input
    NumericVector input;
    switch (kspec.input_kind) {
    case InputKind::ParInput:
      input = pt.column_by_name(kspec.par_input_name);
      break;
    case InputKind::Covariate:
      input = kspec.kernel_input;
      break;
    case InputKind::None:
    default:
      stop("TrendOp '%s': kernel '%s' has neither covariate nor par_input",
           spec.target_param.c_str(), kspec.kernel.c_str());
    }

    if (input.size() != n) {
      stop("TrendRuntime::run_kernels_for_op('%s'): input length (%d) != n_trials (%d)",
           spec.target_param.c_str(), input.size(), n);
    }

    // Parameter columns for this kernel
    KernelParsView kp_view = make_kernel_pars_view(pt, k_rt.kernel_par_indices);

    k_rt.kernel_ptr->reset();
    k_rt.kernel_ptr->run(kp_view, input, spec.comp_index);

    if (spec.has_at) {
      k_rt.kernel_ptr->do_expand(spec.expand_idx);
    }
  }
}

void TrendRuntime::apply_base_for_op(TrendOpRuntime& op,
                                     ParamTable& pt) {
  using namespace Rcpp;

  const TrendOpSpec& spec = *op.spec;
  const int n = pt.n_trials;

  // Ensure kernels have run
  for (auto& krt : op.kernels) {
    if (!krt.kernel_ptr->has_run()) {
      run_kernels_for_op(op, pt);
      break;
    }
  }

  const int K = static_cast<int>(op.kernels.size());
  if (K == 0) return;

  // Gather trajectories
  std::vector<const std::vector<double>*> trajs(K);
  for (int k = 0; k < K; ++k) {
    const auto& v = op.kernels[k].kernel_ptr->get_output();
    if ((int)v.size() != n) {
      stop("apply_base_for_op('%s'): trajectory length mismatch",
           spec.target_param.c_str());
    }
    trajs[k] = &v;
  }

  int target_idx = pt.base_index_for(spec.target_param);
  double* target_col = &pt.base(0, target_idx);

  const std::string& base = spec.base_type;
  bool needs_base_par =
    (base == "lin" || base == "exp_lin" ||
    base == "lin_exp" || base == "centered");

  // Check if any slot has maps; assume all slots share same #maps if they do.
  const KernelSlotSpec& first_slot = *op.kernels[0].spec;
  bool slots_have_maps = first_slot.has_covariate_maps;
  int n_maps = slots_have_maps
  ? static_cast<int>(first_slot.covariate_map_cols.size())
    : 0;

  // Sanity check: all slots should agree on map count
  if (slots_have_maps) {
    for (int k = 1; k < K; ++k) {
      const KernelSlotSpec& s = *op.kernels[k].spec;
      if (!s.has_covariate_maps ||
          (int)s.covariate_map_cols.size() != n_maps) {
        stop("TrendOp '%s': inconsistent map usage across kernel slots",
             spec.target_param.c_str());
      }
    }
  }

  // Base params
  const int n_base_pars = (int)op.base_par_indices.size();
  const bool has_base_pars = n_base_pars > 0;

  const double* base_col0 = nullptr;
  if (has_base_pars) {
    base_col0 = &pt.base(0, op.base_par_indices[0]);
  }

  for (int r = 0; r < n; ++r) {
    double p = target_col[r];
    if (NumericVector::is_na(p) || std::isnan(p)) continue;

    double contrib = 0.0;

    if (slots_have_maps && needs_base_par) {
      // --- maps + linear-type base: use one base param per map, shared across slots ---
      if (n_base_pars != n_maps) {
        stop("TrendOp '%s': #base_pars (%d) != #maps (%d)",
             spec.target_param.c_str(), n_base_pars, n_maps);
      }

      for (int k = 0; k < K; ++k) {
        const KernelSlotSpec& slot = *op.kernels[k].spec;
        double q = (*trajs[k])[r];
        if (NumericVector::is_na(q) || std::isnan(q)) q = 0.0;

        for (int m = 0; m < n_maps; ++m) {
          double base_val = pt.base(r, op.base_par_indices[m]);
          double map_val  = slot.covariate_map_cols[m][r];
          contrib += q * base_val * map_val;
        }
      }

    } else {
      // --- no maps: combine kernel outputs, then apply base type ---
      double q_combined = 0.0;
      for (int k = 0; k < K; ++k) {
        double q = (*trajs[k])[r];
        if (NumericVector::is_na(q) || std::isnan(q)) q = 0.0;
        q_combined += q;
      }

      double tmp = q_combined;
      if (needs_base_par) {
        double tp = base_col0 ? base_col0[r] : 1.0;
        if (base == "lin" || base == "exp_lin" || base == "lin_exp") {
          tmp *= tp;
        } else if (base == "centered") {
          tmp = (tmp - 0.5) * tp;
        }
      }
      contrib = tmp;
    }

    // Base term from target value
    double base_term;
    if (base == "exp_lin" || base == "lin_exp") {
      base_term = std::exp(p);
    } else {
      base_term = p;
    }

    target_col[r] = base_term + contrib;
  }
}

void TrendRuntime::reset_all_kernels() {
  // lambda to loop over a vector of kernels
  auto reset_vec = [](std::vector<TrendOpRuntime>& ops) {
    for (auto& op : ops) {
      for (auto& k_rt : op.kernels) {
        k_rt.kernel_ptr->reset();
      }
    }
  };
  reset_vec(premap_ops);
  reset_vec(pretransform_ops);
  reset_vec(posttransform_ops);
}


Rcpp::NumericMatrix TrendRuntime::all_kernel_outputs(ParamTable& pt) {
  using namespace Rcpp;

  const int n = pt.n_trials;

  // Count total number of kernel slots
  int n_slots = 0;
  for (const auto& op : premap_ops)        n_slots += op.kernels.size();
  for (const auto& op : pretransform_ops)  n_slots += op.kernels.size();
  for (const auto& op : posttransform_ops) n_slots += op.kernels.size();

  NumericMatrix out(n, n_slots);
  CharacterVector cn(n_slots);

  int col = 0;

  auto fill_for_ops = [&](std::vector<TrendOpRuntime>& ops) {
    for (auto& op : ops) {
      const TrendOpSpec& spec = *op.spec;

      for (auto& k_rt : op.kernels) {
        const KernelSlotSpec& kspec = *k_rt.spec;
        const std::vector<double>& traj = k_rt.kernel_ptr->get_output();

        if ((int)traj.size() != n) {
          stop("TrendRuntime::all_kernel_outputs('%s'): trajectory length (%d) != n_trials (%d)",
               spec.target_param.c_str(), (int)traj.size(), n);
        }

        for (int r = 0; r < n; ++r) {
          out(r, col) = traj[r];
        }

        // Build a column name: target_param + "." + input_name
        std::string input_name;
        if (kspec.input_kind == InputKind::Covariate) {
          // try to get column name from attributes if available
          if (kspec.kernel_input.hasAttribute("names")) {
            // optional; often covariate is directly from data[cov_name],
            // so we don't have the name here. You can store the cov_name in KernelSlotSpec if needed.
            input_name = "cov";
          } else {
            input_name = "cov";
          }
        } else if (kspec.input_kind == InputKind::ParInput) {
          input_name = kspec.par_input_name;
        } else {
          input_name = "noinput";
        }

        std::string cname = spec.target_param + "." + input_name;
        cn[col] = cname;

        ++col;
      }
    }
  };

  fill_for_ops(premap_ops);
  fill_for_ops(pretransform_ops);
  fill_for_ops(posttransform_ops);

  colnames(out) = cn;
  return out;
}

// Rcpp::NumericMatrix TrendRuntime::all_kernel_outputs(ParamTable& pt) {
//   using namespace Rcpp;
//
//   const int n = pt.n_trials;
//   const int n_ops =
//     static_cast<int>(premap_ops.size() +
//     pretransform_ops.size() +
//     posttransform_ops.size());
//
//   NumericMatrix out(n, n_ops);
//   CharacterVector cn(n_ops);
//
//   int col = 0;
//
//   auto fill_for_ops = [&](std::vector<TrendOpRuntime>& ops, const char* phase_label) {
//     for (size_t i = 0; i < ops.size(); ++i) {
//       TrendOpRuntime& op = ops[i];
//       const TrendOpSpec& spec = *op.spec;
//
//       // Ensure kernel has been run for this ParamTable
//       // run_kernel_for_op(op, pt);
//
//       const std::vector<double>& traj = op.kernel_ptr->get_output();
//       if (static_cast<int>(traj.size()) != n) {
//         stop("TrendRuntime::all_kernel_outputs('%s'): trajectory length (%d) "
//                "!= n_trials (%d)",
//                op.spec->target_param.c_str(), (int)traj.size(), n);
//       }
//
//       // Copy trajectory into column 'col'
//       for (int r = 0; r < n; ++r) {
//         out(r, col) = traj[r];
//       }
//
//       // --- build name: target_param.covariate_name/par_input_name ---
//       std::string cov_name;
//
//       if (spec.spec.containsElementNamed("covariate") &&
//           !Rf_isNull(spec.spec["covariate"])) {
//           SEXP cov_spec = spec.spec["covariate"];
//         if (Rf_isString(cov_spec) && Rf_length(cov_spec) == 1) {
//           cov_name = Rcpp::as<std::string>(cov_spec);
//         } else {
//           cov_name = "cov_inline";
//         }
//       } else if (spec.spec.containsElementNamed("par_input") &&
//         !Rf_isNull(spec.spec["par_input"])) {
//         SEXP par_spec = spec.spec["par_input"];
//         if (Rf_isString(par_spec) && Rf_length(par_spec) == 1) {
//           cov_name = Rcpp::as<std::string>(par_spec);
//         } else {
//           cov_name = "par_inline";
//         }
//       } else {
//         cov_name = "noinput";
//       }
//
//       std::string cname = spec.target_param + "." + cov_name;
//       cn[col] = cname;
//
//       ++col;
//     }
//   };
//
//   fill_for_ops(premap_ops,        "premap");
//   fill_for_ops(pretransform_ops,  "pretransform");
//   fill_for_ops(posttransform_ops, "posttransform");
//
//   colnames(out) = cn;
//   return out;
// }



// void TrendRuntime::apply_base_for_op(TrendOpRuntime& op,
//                                      ParamTable& pt)
// {
//   using namespace Rcpp;
//
//   const TrendOpSpec& spec = *op.spec;
//
//   // Ensure kernel has been run
//   if (!op.kernel_ptr->has_run()) {
//     run_kernel_for_op(op, pt);
//   }
//
//   const std::vector<double>& traj = op.kernel_ptr->get_output();
//   const int n = pt.n_trials;
//   if ((int)traj.size() != n) {
//     stop("TrendRuntime::apply_base_for_op('%s'): trajectory length (%d) != n_trials (%d)",
//          spec.target_param.c_str(), (int)traj.size(), n);
//   }
//
//   // // SM TEMP
//   // // DEBUG: print a few trajectory values for pretransform ops
//   // if (spec.phase == TrendPhase::Pretransform) {
//   //   Rcout << "apply_base_for_op[pretransform] target=" << spec.target_param
//   //         << " traj[0..4]=";
//   //   for (int k = 0; k < std::min(n, 5); ++k) {
//   //     Rcout << traj[k] << " ";
//   //   }
//   //   Rcout << "\n";
//   // }
//
//   const std::string& base = spec.base_type;
//   if(base == "exp_lin" || base == "lin_exp") {
//     Rcpp::stop("exp_lin and lin_exp are deprecated. Please use pretransform instead.");
//   }
//
//   // Target parameter column
//   int target_idx = pt.base_index_for(spec.target_param);
//   double* target_col = &pt.base(0, target_idx);
//
//   // Base/trend parameter if needed
//   // const double* trend_col = nullptr;
//   bool needs_trend_par =
//     (base == "lin" || base == "exp_lin" || base == "lin_exp" || base == "centered");
//
//   const double* single_trend_col = nullptr;
//   if (needs_trend_par && !spec.has_covariate_maps) {
//     if (op.base_par_idx < 0) {
//       stop("TrendOp '%s': base_type = '%s' but base_par_idx < 0",
//            spec.target_param.c_str(), base.c_str());
//     }
//     single_trend_col = &pt.base(0, op.base_par_idx);
//   }
//
//   for (int r = 0; r < n; ++r) {
//     double p = target_col[r];
//     if (NumericVector::is_na(p) || std::isnan(p)) {
//       continue;
//     }
//
//     double q = traj[r];
//     if (NumericVector::is_na(q) || std::isnan(q)) {
//       q = 0.0;
//     }
//
//     double contrib = 0.0;
//
//     if (spec.has_covariate_maps) {
//       const int K = static_cast<int>(spec.covariate_map_cols.size());
//       if (K != (int)op.base_par_indices.size()) {
//         stop("TrendOp '%s': mismatch between covariate_maps (%d) and base_pars (%d)",
//              spec.target_param.c_str(), K, (int)op.base_par_indices.size());
//       }
//
//       for (int k = 0; k < K; ++k) {
//         double base_val = pt.base(r, op.base_par_indices[k]);
//         double map_val  = spec.covariate_map_cols[k][r]; // column view
//         contrib += q * base_val * map_val;
//       }
//     } else {
//       // Single base parameter behaviour
//       double tmp = q;
//       if (base == "lin" || base == "exp_lin" || base == "lin_exp") {
//         double tp = single_trend_col ? single_trend_col[r] : 1.0;
//         tmp *= tp;
//       } else if (base == "centered") {
//         double tp = single_trend_col ? single_trend_col[r] : 1.0;
//         tmp = (tmp - 0.5) * tp;
//       }
//       contrib = tmp;
//     }
//
//     double base_term;
//     if (base == "exp_lin" || base == "lin_exp") {
//       base_term = std::exp(p);
//     } else if (base == "lin" || base == "centered" || base == "add" || base.empty()) {
//       base_term = p;
//     } else {
//       base_term = p;
//     }
//
//     target_col[r] = contrib + base_term;
//   }
// }



// void TrendRuntime::bind_all_ops_to_paramtable(const ParamTable& pt) {
//   for (auto& op : premap_ops)        bind_trend_op_to_paramtable(op, pt);
//   for (auto& op : pretransform_ops)  bind_trend_op_to_paramtable(op, pt);
//   for (auto& op : posttransform_ops) bind_trend_op_to_paramtable(op, pt);
//
//   // Validate par_inputs --
//   auto validate_par_inputs = [&](const TrendOpRuntime& op) {
//     const TrendOpSpec& spec = *op.spec;
//     if (spec.spec.containsElementNamed("par_input") &&
//         !Rf_isNull(spec.spec["par_input"])) {
//           SEXP par_spec = spec.spec["par_input"];
//           if (TYPEOF(par_spec) == STRSXP) {
//             Rcpp::CharacterVector cv(par_spec);
//             if (cv.size() == 1) {
//               std::string par_name = Rcpp::as<std::string>(cv[0]);
//               (void) pt.base_index_for(par_name); // throws if unknown
//             }
//         }
//     }
//   };
//
//   for (auto& op : premap_ops)        validate_par_inputs(op);
//   for (auto& op : pretransform_ops)  validate_par_inputs(op);
//   for (auto& op : posttransform_ops) validate_par_inputs(op);
// }

// void TrendRuntime::run_kernel_for_op(TrendOpRuntime& op,
//                                      ParamTable& pt)
// {
//   using namespace Rcpp;
//   using std::string;
//
//   const TrendOpSpec& spec = *op.spec;
//
//   NumericVector input;
//
//   // waht type of input?
//   switch (spec.input_kind) {
//   case InputKind::ParInput:
//     // Resolve from ParamTable every time (current particle)
//     input = pt.column_by_name(spec.par_input_name);
//     break;
//
//   case InputKind::Covariate:
//     // Use cached data column
//     input = spec.kernel_input;
//     break;
//
//   case InputKind::None:
//   default:
//     Rcpp::stop("TrendOp for target '%s' has neither covariate nor par_input input",
//                spec.target_param.c_str());
//   }
//
//   const int n = input.size();
//   if (n != pt.n_trials) {
//     Rcpp::stop("TrendRuntime::run_kernel_for_op: input length (%d) != n_trials (%d)",
//                n, pt.n_trials);
//   }
//
//   const int p = pt.base.ncol();
//   for (int idx : op.kernel_par_indices) {
//     if (idx < 0 || idx >= p) {
//       stop("TrendRuntime::run_kernel_for_op('%s'): kernel_par index %d out of [0,%d)",
//            spec.target_param.c_str(), idx, p);
//     }
//   }
//
//   KernelParsView kp_view = make_kernel_pars_view(pt, op.kernel_par_indices);
//
//   op.kernel_ptr->reset();
//   op.kernel_ptr->run(kp_view, input, spec.comp_index);
//
//   if (spec.has_at) {
//     op.kernel_ptr->do_expand(spec.expand_idx);
//   }
// }



// inline void bind_trend_op_to_paramtable(TrendOpRuntime& op,
//                                         const ParamTable& pt)
// {
//   const TrendOpSpec& spec = *op.spec;
//
//   const int n_tp = spec.trend_pnames.size();
//   op.base_par_idx = -1;
//   op.base_par_indices.clear();
//   op.kernel_par_indices.clear();
//
//   int kernel_start = 0;
//
//   bool needs_base_par =
//     (spec.base_type == "lin" || spec.base_type == "exp_lin" ||
//     spec.base_type == "lin_exp" || spec.base_type == "centered");
//
//   if (needs_base_par) {
//     int n_base_pars = 0;
//     if (spec.has_covariate_maps) {
//       n_base_pars = static_cast<int>(spec.covariate_map_cols.size());  // one per map
//     } else {
//       n_base_pars = 1;
//     }
//
//     if (n_tp < n_base_pars) {
//       Rcpp::stop("TrendOp '%s': base_type '%s' expects at least %d trend_pnames (got %d)",
//                  spec.target_param.c_str(), spec.base_type.c_str(), n_base_pars, n_tp);
//     }
//
//     // First n_base_pars trend_pnames are base parameters
//     for (int b = 0; b < n_base_pars; ++b) {
//       std::string base_nm = Rcpp::as<std::string>(spec.trend_pnames[b]);
//       int idx = pt.base_index_for(base_nm);
//       op.base_par_indices.push_back(idx);
//     }
//     op.base_par_idx = op.base_par_indices[0];
//     kernel_start = n_base_pars;
//   }
//
//   // Remaining trend_pnames are kernel parameters
//   for (int k = kernel_start; k < n_tp; ++k) {
//     std::string kn = Rcpp::as<std::string>(spec.trend_pnames[k]);
//     int idx = pt.base_index_for(kn);
//     op.kernel_par_indices.push_back(idx);
//   }
//   // Rprintf("Binding...");
//
//   // After building base_par_indices and kernel_par_indices:
// //#if 0
//   // Rcpp::Rcout << "TrendOp target=" << spec.target_param
//   //             << " kernel=" << spec.kernel
//   //             << " base_type=" << spec.base_type
//   //             << " has_maps=" << spec.has_covariate_maps
//   //             << " n_trend_pnames=" << n_tp
//   //             << " n_base_pars=" << op.base_par_indices.size()
//   //             << " n_kernel_pars=" << op.kernel_par_indices.size()
//   //             << "\n";
// //#endif
// }



// TrendPlan::TrendPlan(Rcpp::Nullable<Rcpp::List> trend_,
//                      const Rcpp::DataFrame& data_)
//   : data(data_) {
//
//   if (trend_.isNull()) {
//     trend = R_NilValue;
//     return;
//   }
//
//   trend = Rcpp::List(trend_.get());
//
//   using std::string;
//   using namespace Rcpp;
//
//   // Get numeric covariate maps from data attribute
//   if (data_.hasAttribute("covariate_maps")) {
//     SEXP cm = data_.attr("covariate_maps");
//     if (!Rf_isNull(cm)) {
//       data_covariate_maps = Rcpp::List(cm);
//       has_data_covariate_maps = data_covariate_maps.size() > 0;
//     }
//   }
//
//   // helper functions to figure out whether there's covariate_inputs and/or par_inputs
//   auto uses_cov_input = [](const Rcpp::List& tr) -> bool {
//     if (!tr.containsElementNamed("covariate")) return false;
//
//     SEXP x = tr["covariate"];
//     if (Rf_isNull(x)) return false;
//
//     if (TYPEOF(x) == STRSXP) {
//       // character vector: length == 0 => no covariate
//       return Rf_length(x) > 0;
//     }
//
//     // non-string (e.g., numeric vector) => inline covariate
//     return true;
//   };
//
//   auto uses_par_input = [](const Rcpp::List& tr) -> bool {
//     if (!tr.containsElementNamed("par_input")) return false;
//
//     SEXP x = tr["par_input"];
//     if (Rf_isNull(x)) return false;
//
//     if (TYPEOF(x) == STRSXP) {
//       // empty character vector => no par_input
//       return Rf_length(x) > 0;
//     }
//
//     // non-string (e.g., numeric vector) => inline par_input
//     return true;
//   };
//
//   CharacterVector tnames = trend.names();
//   if (tnames.size() != trend.size()) {
//     stop("TrendPlan: 'trend' must be a named list");
//   }
//
//   premap_ops.clear();
//   pretransform_ops.clear();
//   posttransform_ops.clear();
//   premap_trend_params.clear();
//   pretransform_trend_params.clear();
//   posttransform_trend_params.clear();
//   all_trend_params.clear();
//
//   // Loopy over trends
//   for (int i = 0; i < trend.size(); ++i) {
//     if (trend[i] == R_NilValue) continue;
//
//     Rcpp::List tr_i = trend[i];
//     std::string name_i = Rcpp::as<std::string>(tnames[i]);
//
//     // Build once to inspect phase, trend_pnames, etc.
//     TrendOpSpec tmp(tr_i, name_i);
//
//     // --- collect trend_pnames into sets -- useful for transformations later on ---
//     for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
//       std::string pn = Rcpp::as<std::string>(tmp.trend_pnames[k]);
//       all_trend_params.insert(pn);
//     }
//
//     if (tmp.phase == TrendPhase::Premap) {
//       for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
//         premap_trend_params.insert(Rcpp::as<std::string>(tmp.trend_pnames[k]));
//       }
//     }
//     if (tmp.phase == TrendPhase::Pretransform) {
//       for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
//         pretransform_trend_params.insert(Rcpp::as<std::string>(tmp.trend_pnames[k]));
//       }
//     }
//     if (tmp.phase == TrendPhase::Posttransform) {
//       for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
//         posttransform_trend_params.insert(Rcpp::as<std::string>(tmp.trend_pnames[k]));
//       }
//     }
//
//     // ---- decide which inputs are actually used in the current trend ----
//     bool has_cov_input = uses_cov_input(tr_i);
//     bool has_par_input = uses_par_input(tr_i);
//
//     if (!has_cov_input && !has_par_input) {
//       Rcpp::stop("TrendPlan: spec for '%s' must have non-empty 'covariate' or 'par_input'",
//                  name_i.c_str());
//     }
//
//     // Extract names only if that input type is used
//     std::vector<std::string> cov_names;
//     std::vector<std::string> par_input_names;
//     if (has_cov_input) {
//       cov_names = covariate_names_from_spec(tr_i);
//     }
//     if (has_par_input) {
//       par_input_names = par_input_names_from_spec(tr_i);
//     }
//
//     // Helper: add an op to the correct phase vector and init input
//     auto add_op_for_phase = [&](TrendPhase phase,
//                                 const Rcpp::List& spec_list,
//                                 bool use_cov_input,
//                                 bool use_par_input) {
//       switch (phase) {
//       case TrendPhase::Premap:
//         premap_ops.emplace_back(spec_list, name_i);
//         if (use_cov_input) {
//           premap_ops.back().input_kind = InputKind::Covariate;
//           init_covariate(premap_ops.back(), data);
//           init_covariate_maps_for_op(premap_ops.back(), data, data_covariate_maps);
//         } else if (use_par_input) {
//           // par_input: record name
//           Rcpp::CharacterVector cv(spec_list["par_input"]);
//           if (cv.size() != 1) {
//             Rcpp::stop("TrendOpSpec '%s': par_input must be single name after expansion",
//                        name_i.c_str());
//           }
//           premap_ops.back().input_kind = InputKind::ParInput;
//           premap_ops.back().par_input_name = Rcpp::as<std::string>(cv[0]);
//         } else {
//           premap_ops.back().input_kind = InputKind::None;
//         }
//         break;
//
//       case TrendPhase::Pretransform:
//         pretransform_ops.emplace_back(spec_list, name_i);
//         if (use_cov_input) {
//           pretransform_ops.back().input_kind = InputKind::Covariate;
//           init_covariate(pretransform_ops.back(), data);
//           init_covariate_maps_for_op(pretransform_ops.back(), data, data_covariate_maps);
//         } else if (use_par_input) {
//           // par_input: record name
//           Rcpp::CharacterVector cv(spec_list["par_input"]);
//           if (cv.size() != 1) {
//             Rcpp::stop("TrendOpSpec '%s': par_input must be single name after expansion",
//                        name_i.c_str());
//           }
//           pretransform_ops.back().input_kind = InputKind::ParInput;
//           pretransform_ops.back().par_input_name = Rcpp::as<std::string>(cv[0]);
//         } else {
//           pretransform_ops.back().input_kind = InputKind::None;
//         }
//         break;
//
//       case TrendPhase::Posttransform:
//         posttransform_ops.emplace_back(spec_list, name_i);
//         if (use_cov_input) {
//           posttransform_ops.back().input_kind = InputKind::Covariate;
//           init_covariate(posttransform_ops.back(), data);
//           init_covariate_maps_for_op(posttransform_ops.back(), data, data_covariate_maps);
//         } else if (use_par_input) {
//           // par_input: record name
//           Rcpp::CharacterVector cv(spec_list["par_input"]);
//           if (cv.size() != 1) {
//             Rcpp::stop("TrendOpSpec '%s': par_input must be single name after expansion",
//                        name_i.c_str());
//           }
//           posttransform_ops.back().input_kind = InputKind::ParInput;
//           posttransform_ops.back().par_input_name = Rcpp::as<std::string>(cv[0]);
//         } else {
//           posttransform_ops.back().input_kind = InputKind::None;
//         }
//         break;
//       }
//     };
//
//     // ---- 1) Covariate-based ops (0, 1, or many) ----
//     if (has_cov_input) {
//       if (cov_names.empty()) {
//         // inline numeric covariate
//         add_op_for_phase(tmp.phase, tr_i, /*use_cov_input=*/true, false);
//       } else if (cov_names.size() == 1) {
//         // single named column
//         add_op_for_phase(tmp.phase, tr_i, /*use_cov_input=*/true, false);
//       } else {
//         // multiple covariate columns: clone spec per column
//         for (const std::string& cov_name : cov_names) {
//           Rcpp::List tr_copy = Rcpp::clone(tr_i);
//           tr_copy["covariate"] = cov_name;
//           add_op_for_phase(tmp.phase, tr_copy, /*use_cov_input=*/true, false);
//         }
//       }
//     }
//
//     // ---- 2) par_input-based ops (0, 1, or many) ----
//     if (has_par_input) {
//       if (par_input_names.empty()) {
//         // inline numeric par_input
//         add_op_for_phase(tmp.phase, tr_i, /*use_cov_input=*/false, true);
//       } else if (par_input_names.size() == 1) {
//         // single parameter name
//         add_op_for_phase(tmp.phase, tr_i, /*use_cov_input=*/false, true);
//       } else {
//         // multiple par_input names: clone spec per parameter
//         for (const std::string& par_name : par_input_names) {
//           Rcpp::List tr_copy = Rcpp::clone(tr_i);
//           tr_copy["par_input"] = par_name;
//           add_op_for_phase(tmp.phase, tr_copy, /*use_cov_input=*/false, true);
//         }
//       }
//     }
//   }
//
//   // Build first-level masks & comp_index/expand_idx
//   for (auto& op : premap_ops)        op.make_first_level(data);
//   for (auto& op : pretransform_ops)  op.make_first_level(data);
//   for (auto& op : posttransform_ops) op.make_first_level(data);
// }

