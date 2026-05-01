#include "TrendEngine.h"
#include <Rcpp.h>
using namespace Rcpp;


// ============================================================
// New helper: init_kernel_args_for_slot()
// Call this during TrendPlan construction, after the slot's
// kernel_type is known, alongside init_covariate_for_slot() etc.
// ============================================================
void init_kernel_args_for_slot(KernelSlotSpec& slot,
                               const Rcpp::List& tr,       // one trend list entry
                               const Rcpp::DataFrame& data)
{
  // Nothing to do if no kernel_args field in the trend spec
  if (!tr.containsElementNamed("kernel_args")) {
    slot.build_kernel_args();
    return;
  }

  SEXP ka_sexp = tr["kernel_args"];
  if (Rf_isNull(ka_sexp)) {
    slot.build_kernel_args();
    return;
  }

  Rcpp::List ka = Rcpp::as<Rcpp::List>(ka_sexp);

  // ---- q_reset_column ----
  if (ka.containsElementNamed("q_reset_column")) {
    SEXP col_name_sexp = ka["q_reset_column"];
    if (!Rf_isNull(col_name_sexp)) {
      std::string col_name = Rcpp::as<std::string>(col_name_sexp);

      if (!data.containsElementNamed(col_name.c_str())) {
        Rcpp::stop("kernel_args$q_reset_column: column '%s' not found in data",
                   col_name.c_str());
      }

      SEXP col = data[col_name.c_str()];

      // Accept logical or integer columns; coerce to integer for raw-pointer access
      if (TYPEOF(col) == LGLSXP) {
        slot.q_reset_col = Rcpp::as<Rcpp::IntegerVector>(col);  // TRUE->1, FALSE->0, NA->NA_INTEGER
      } else if (TYPEOF(col) == INTSXP) {
        slot.q_reset_col = Rcpp::IntegerVector(col);
      } else {
        Rcpp::stop("kernel_args$q_reset_column: column '%s' must be logical or integer",
                   col_name.c_str());
      }

      // Rprintf("q_reset_size=%d\n", (int)slot.q_reset_col.size());
      // Validate length matches data
      if (slot.q_reset_col.size() != data.nrows()) {
        Rcpp::stop("kernel_args$q_reset_column: column '%s' has length %d, expected %d",
                   col_name.c_str(), (int)slot.q_reset_col.size(), (int)data.nrows());
      }
    }
  }

  // Future kernel_args fields parsed here in the same pattern.

  // Sync raw-pointer view from all populated fields
  slot.build_kernel_args();
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

  // Wrap in 1-column matrix
  NumericVector v = data[cov_name];
  int n_trials = data.nrows();
  slot.kernel_input = NumericMatrix(n_trials, 1);
  slot.kernel_input_name = cov_name;
  slot.kernel_input(_, 0) = v;
}


// Initialize covariate for a single kernel slot when covariate is a *single* name
void init_combined_for_slot(KernelSlotSpec& slot,
                            const Rcpp::List& spec,
                            const Rcpp::DataFrame& data) {

  Rcpp::CharacterVector cov_spec = spec["covariate"];
  Rcpp::CharacterVector par_input = spec["par_input"];

  int n_covariates = cov_spec.size();
  int n_par_input  = par_input.size();
  int n_trials     = data.nrow();

  slot.kernel_input = Rcpp::NumericMatrix(n_trials, n_covariates + n_par_input);
  slot.par_input_names = par_input;
  slot.covariate_names = cov_spec;

  Rcpp::CharacterVector col_names(n_covariates + n_par_input);

  // copy covariate columns
  for (int i = 0; i < n_covariates; i++) {
    Rcpp::String cov_name = cov_spec[i];
    slot.kernel_input.column(i) = Rcpp::as<Rcpp::NumericVector>(data[cov_name]);
    col_names[i] = cov_name;
    slot.covariate_indices.push_back(i);
  }

  // set parameter-input column names (values filled later)
  for (int i = 0; i < n_par_input; i++) {
    col_names[n_covariates + i] = par_input[i];
    slot.par_input_indices.push_back(n_covariates+i);
  }

  Rcpp::colnames(slot.kernel_input) = col_names;
}

// Initialize covariate for a kernel slot when covariate is a *vector* of column names.
// Builds an n_trials x N matrix from all named columns.
// Used by variadic kernels (e.g. RescorlaWagner) that consume all covariates at once.
void init_multicovariate_for_slot(KernelSlotSpec& slot,
                                  const Rcpp::List& spec,
                                  const Rcpp::DataFrame& data) {
  using namespace Rcpp;

  if (!spec.containsElementNamed("covariate")) {
    stop("init_multicovariate_for_slot: missing 'covariate' field");
  }

  CharacterVector cov_spec(spec["covariate"]);
  const int n_covs   = cov_spec.size();
  const int n_trials = data.nrow();

  if (n_covs == 0) {
    stop("init_multicovariate_for_slot: 'covariate' must name at least one column");
  }

  slot.kernel_input  = NumericMatrix(n_trials, n_covs);
  slot.covariate_names = cov_spec;

  CharacterVector col_names(n_covs);

  for (int i = 0; i < n_covs; ++i) {
    std::string cov_name = as<std::string>(cov_spec[i]);
    if (!data.containsElementNamed(cov_name.c_str())) {
      stop("init_multicovariate_for_slot: data has no column named '%s'",
           cov_name.c_str());
    }
    slot.kernel_input(_, i) = as<NumericVector>(data[cov_name]);
    col_names[i]             = cov_name;
    slot.covariate_indices.push_back(i);
  }

  colnames(slot.kernel_input) = col_names;

  // kernel_input_name: join all covariate names for display purposes
  // e.g. "cov_A+cov_B+cov_C"
  std::string combined_name;
  for (int i = 0; i < n_covs; ++i) {
    if (i > 0) combined_name += "+";
    combined_name += as<std::string>(cov_spec[i]);
  }
  slot.kernel_input_name = combined_name;
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


// Initialize covariate maps for a variadic kernel slot (e.g. RescorlaWagner).
// Each map in data_covmaps is an [n_rows x n_covs] matrix — stored as-is,
// no column extraction needed. apply_base_for_op indexes them directly.
void init_multicovariate_maps_for_slot(KernelSlotSpec& slot,
                                       const Rcpp::List& spec,
                                       const Rcpp::DataFrame& data,
                                       const Rcpp::List& data_covmaps) {
  using namespace Rcpp;
  using std::string;

  if (!spec.containsElementNamed("map") || Rf_isNull(spec["map"])) {
    slot.has_covariate_maps = false;
    slot.covariate_map_mats.clear();
    return;
  }

  if (data_covmaps.size() == 0) {
    stop("Trend: 'map' specified but data has no 'covariate_maps' attribute");
  }

  List maps_spec(spec["map"]);
  CharacterVector map_names = maps_spec.names();
  if (map_names.size() != maps_spec.size()) {
    stop("Trend: 'map' list must be named");
  }

  List maps_data(data_covmaps);
  CharacterVector data_map_names = maps_data.names();

  const int M = maps_spec.size();
  const int T = data.nrows();

  slot.covariate_map_mats.resize(M);

  for (int m = 0; m < M; ++m) {
    string map_nm = as<string>(map_names[m]);

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

    slot.covariate_map_mats[m] = mat;
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

    // Custom kernels takes *all* covariates and par_input as input
    if(ktype == KernelType::Custom) {
      KernelSlotSpec slot;
      slot.kernel      = kernel_name;
      slot.kernel_type = ktype;
      slot.input_kind  = InputKind::Combined;

      // custom_ptr is stored as an attribute on the trend spec entry
      if (!tr_i.hasAttribute("custom_ptr")) {
        Rcpp::stop("Custom kernel '%s' requires 'custom_ptr' attribute", name_i.c_str());
      }

      SEXP custom_attr = tr_i.attr("custom_ptr");
      if (Rf_isNull(custom_attr)) {
        Rcpp::stop("Custom kernel '%s': 'custom_ptr' attribute is NULL", name_i.c_str());
      }
      slot.custom_fun = custom_attr;  // store SEXP (externalptr)

      // spec copy with *single* covariate name
      Rcpp::List tr_copy = Rcpp::clone(tr_i);

      init_combined_for_slot(slot, tr_copy, data);
      op.kernels.push_back(std::move(slot));
    } else {

      // ---- Covariate-based kernel slots ----
      if (has_cov) {
        SEXP cov_spec = tr_i["covariate"];

        if (TYPEOF(cov_spec) == STRSXP) {
          Rcpp::CharacterVector cv(cov_spec);

          if (kmeta.input_arity == -1) {
            // Variadic kernel (e.g. RescorlaWagner): one slot receives all
            // covariates together as a matrix. No expansion into separate slots.
            KernelSlotSpec slot;
            slot.kernel      = kernel_name;
            slot.kernel_type = ktype;
            slot.input_kind  = InputKind::Covariate;

            // Pass the full covariate vector as-is; init_multicovariate_for_slot
            init_multicovariate_for_slot(slot, tr_i, data);
            init_kernel_args_for_slot(slot, tr_i, data);

            if (has_data_covariate_maps) {
              // no column split of covariate maps - just the entire matrix again
              init_multicovariate_maps_for_slot(slot, tr_i, data, data_covariate_maps);
            } else {
              if (tr_i.containsElementNamed("map") && !Rf_isNull(tr_i["map"])) {
                Rcpp::stop("TrendPlan: spec '%s' has 'map' but data has no "
                             "'covariate_maps' attribute", name_i.c_str());
              }
            }

            op.kernels.push_back(std::move(slot));

          } else {
            // Standard kernels: expand one spec into N single-covariate slots
            for (int j = 0; j < cv.size(); ++j) {
              std::string cov_name = Rcpp::as<std::string>(cv[j]);

              if (!(kmeta.supports_grouping && kmeta.input_arity == 1)) {
                Rcpp::stop("Kernel '%s' does not support grouped covariates",
                           kernel_name.c_str());
              }

              KernelSlotSpec slot;
              slot.kernel      = kernel_name;
              slot.kernel_type = ktype;
              slot.input_kind  = InputKind::Covariate;

              Rcpp::List tr_copy = Rcpp::clone(tr_i);
              tr_copy["covariate"] = cov_name;

              init_covariate_for_slot(slot, tr_copy, data);
              init_kernel_args_for_slot(slot, tr_copy, data);

              if (has_data_covariate_maps) {
                init_covariate_maps_for_slot(slot, tr_copy, data, data_covariate_maps);
              } else {
                if (tr_copy.containsElementNamed("map") && !Rf_isNull(tr_copy["map"])) {
                  Rcpp::stop("TrendPlan: spec '%s' has 'map' but data has no "
                               "'covariate_maps' attribute", name_i.c_str());
                }
              }

              op.kernels.push_back(std::move(slot));
            }
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
            slot.kernel_input  = NumericMatrix(data_.nrow(), 1); // initialize empty matrix

            op.kernels.push_back(std::move(slot));
          }
        } else {
          // inline par_input (numeric) not implemented
          stop("TrendPlan: inline numeric 'par_input' not yet supported for '%s'",
               name_i.c_str());
        }
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

      bool needs_base_par = (base == "lin" || base == "centered");

      int n_base_pars = 0;
      bool slot_has_maps = !spec.kernels.empty() && spec.kernels[0].has_covariate_maps;
      if (needs_base_par) {
        if (slot_has_maps) {
          const KernelSlotSpec& slot0 = spec.kernels[0];

          if (!slot0.covariate_map_mats.empty()) {
            // Variadic kernel (e.g. RescorlaWagner): one base param per map matrix
            n_base_pars = static_cast<int>(slot0.covariate_map_mats.size());
          } else {
            // Single-covariate kernels: one base param per extracted map column
            n_base_pars = static_cast<int>(slot0.covariate_map_cols.size());
          }
        } else {
          // No maps: a single base parameter
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
        k_rt.kernel_ptr   = make_kernel(kspec.kernel_type, kspec.custom_fun);
        k_rt.kernel_ptr->set_kernel_args(kspec.kernel_args);
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
  const int n_trials = pt.n_trials;

  // since kernels are now kernelgroups, loop over the kernels
  for (auto& k_rt : op.kernels) {
    const KernelSlotSpec& kspec = *k_rt.spec;

    NumericMatrix input = kspec.kernel_input;  // shallow handle

    switch (kspec.input_kind) {

    case InputKind::Combined: {
      // Combined: covariate columns already in kernel_input,
      //           par_input columns are filled from ParamTable at runtime.
      for (int i = 0; i < (int)kspec.par_input_indices.size(); ++i) {
      int col_idx = kspec.par_input_indices[i];
      std::string name = as<std::string>(kspec.par_input_names[i]);

      NumericVector v = pt.column_by_name(name);
      if (v.size() != n_trials) {
        stop("TrendRuntime::run_kernels_for_op('%s'): par_input '%s' length (%d) != n_trials (%d)",
             spec.target_param.c_str(), name.c_str(), v.size(), n_trials);
      }
      input(_, col_idx) = v;  // overwrite column (no reallocation)
    }
      break;
    }

    case InputKind::ParInput: {
      // 1-column shape allocated at spec construction; fill from ParamTable.
      NumericVector v = pt.column_by_name(kspec.par_input_name);
      if (v.size() != n_trials) {
        stop("TrendRuntime::run_kernels_for_op('%s'): par_input '%s' length (%d) != n_trials (%d)",
             spec.target_param.c_str(), kspec.par_input_name.c_str(), v.size(), n_trials);
      }
      if (input.ncol() != 1) {
        stop("ParInput kernel '%s' expected 1 input column, got %d",
             kspec.kernel.c_str(), input.ncol());
      }
      input(_, 0) = v;
      break;
    }

    case InputKind::Covariate:
      // Covariate column already stored in kernel_input (1-column matrix)
      // Nothing to do here.
      break;

    case InputKind::None:
    default:
      stop("TrendOp '%s': kernel '%s' has neither covariate nor par_input",
           spec.target_param.c_str(), kspec.kernel.c_str());
    }

    if (input.nrow() != n_trials) {
      stop("TrendRuntime::run_kernels_for_op('%s'): input nrow (%d) != n_trials (%d)",
           spec.target_param.c_str(), input.nrow(), n_trials);
    }


    // Parameter columns for this kernel
    KernelParsView kp_view = make_kernel_pars_view(pt, k_rt.kernel_par_indices);

    k_rt.kernel_ptr->reset();
    k_rt.kernel_ptr->run(kp_view, input, spec.comp_index);

    if (spec.has_at) {
      // record and expand
      k_rt.kernel_ptr->set_expand_idx(spec.expand_idx);
      k_rt.kernel_ptr->do_expand(spec.expand_idx);  // expands out_
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

  // // Gather trajectories
  // std::vector<const std::vector<double>*> trajs(K);
  // for (int k = 0; k < K; ++k) {
  //   const auto& v = op.kernels[k].kernel_ptr->get_output();
  //   if ((int)v.size() != n) {
  //     stop("apply_base_for_op('%s'): trajectory length mismatch",
  //          spec.target_param.c_str());
  //   }
  //   trajs[k] = &v;
  // }

  int target_idx = pt.base_index_for(spec.target_param);
  double* target_col = &pt.base(0, target_idx);

  const std::string& base = spec.base_type;
  bool needs_base_par = (base == "lin" || base == "exp_lin" ||  base == "lin_exp" || base == "centered");
  const int  n_base_pars = (int)op.base_par_indices.size();
  const bool has_base_pars = n_base_pars > 0;
  const double* base_col0 = has_base_pars ? &pt.base(0, op.base_par_indices[0]) : nullptr;

  // Check if any slot has maps; assume all slots share same #maps if they do.
  const KernelSlotSpec& first_slot = *op.kernels[0].spec;
  bool slots_have_maps = first_slot.has_covariate_maps;
  int n_maps = slots_have_maps ? static_cast<int>(first_slot.covariate_map_cols.size()) : 0;

  for (int r = 0; r < n; ++r) {
    double p = target_col[r];
    if (NumericVector::is_na(p) || std::isnan(p)) continue;

    double contrib = 0.0;

    // ---- Accumulate kernel contributions ----
    // Each kernel slot k exposes a KernelOutput with n_cols columns.
    // For all current single-covariate kernels: n_cols == 1.
    // For RescorlaWagner: n_cols == N covariates.
    // covariate_map_cols[col] aligns with KernelOutput column col.

    for (int k = 0; k < K; ++k) {
      KernelOutput traj = op.kernels[k].kernel_ptr->get_output_stream(1);
      const KernelSlotSpec& slot = *op.kernels[k].spec;
      const int n_cols = traj.n_cols;

      if (slots_have_maps && needs_base_par) {

        if (!slot.covariate_map_mats.empty()) {
          // ---- Variadic kernel (e.g. RescorlaWagner) ----
          // covariate_map_mats: one [n_rows x n_covs] matrix per map.
          // One base param per map, applied across all covariate columns.
          const int n_mat_maps = static_cast<int>(slot.covariate_map_mats.size());
          if (n_base_pars != n_mat_maps) {
            stop("TrendOp '%s': #base_pars (%d) != #covariate_map_mats (%d)",
                 spec.target_param.c_str(), n_base_pars, n_mat_maps);
          }

          for (int m = 0; m < n_mat_maps; ++m) {
            double base_val = pt.base(r, op.base_par_indices[m]);
            const Rcpp::NumericMatrix& mat = slot.covariate_map_mats[m];

            for (int col = 0; col < n_cols; ++col) {
              double q = traj(r, col);
              if (NumericVector::is_na(q) || std::isnan(q)) q = 0.0;
              double map_val = mat(r, col);
              contrib += q * base_val * map_val;
            }
          }

        } else {
          // ---- Single-covariate kernel ----
          // covariate_map_cols: one extracted NumericVector per map/covariate.
          // n_cols == 1 always here; one base param per map.
          if (n_maps != n_base_pars) {
            stop("TrendOp '%s': #base_pars (%d) != #maps (%d)",
                 spec.target_param.c_str(), n_base_pars, n_maps);
          }

          double q = traj(r, 0);  // n_cols == 1 always
          if (NumericVector::is_na(q) || std::isnan(q)) q = 0.0;

          for (int m = 0; m < n_maps; ++m) {
            double base_val = pt.base(r, op.base_par_indices[m]);
            double map_val  = slot.covariate_map_cols[m][r];
            contrib += q * base_val * map_val;
          }
        }

      } else {
        // No maps: sum all columns into q_combined, then apply base type.
        double q_combined = 0.0;
        for (int col = 0; col < n_cols; ++col) {
          double q = traj(r, col);
          if (!NumericVector::is_na(q) && !std::isnan(q)) q_combined += q;
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
        contrib += tmp;
      }
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

Rcpp::NumericMatrix TrendRuntime::all_kernel_outputs(ParamTable& pt,
                                                     const std::vector<int>& codes) {
  using namespace Rcpp;

  const int n = pt.n_trials;

  // ---- Pass 1: ensure all kernels have run, then count total output columns ----
  // get_output_stream() must be called after run(), so we run first, then count.

  auto ensure_run = [&](std::vector<TrendOpRuntime>& ops) {
    for (auto& op : ops) {
      for (auto& k_rt : op.kernels) {
        if (!k_rt.kernel_ptr->has_run()) {
          run_kernels_for_op(op, pt);
          break;  // run_kernels_for_op runs all slots in the op at once
        }
      }
    }
  };

  ensure_run(premap_ops);
  ensure_run(pretransform_ops);
  ensure_run(posttransform_ops);

  // Count total columns: RescorlaWagner streams contribute n_cols > 1
  int n_cols_total = 0;
  auto count_cols = [&](const std::vector<TrendOpRuntime>& ops) {
    for (const auto& op : ops) {
      for (const auto& k_rt : op.kernels) {
        for (int code : codes) {
          if (k_rt.kernel_ptr->has_output_stream(code)) {
            KernelOutput ko = k_rt.kernel_ptr->get_output_stream(code);
            n_cols_total += ko.n_cols;
          }
        }
      }
    }
  };

  count_cols(premap_ops);
  count_cols(pretransform_ops);
  count_cols(posttransform_ops);

  NumericMatrix out(n, n_cols_total);
  CharacterVector cn(n_cols_total);

  int col = 0;

  auto fill_for_ops = [&](std::vector<TrendOpRuntime>& ops) {
    for (auto& op : ops) {
      const TrendOpSpec& spec = *op.spec;

      for (auto& k_rt : op.kernels) {
        const KernelSlotSpec& kspec = *k_rt.spec;

        std::string input_name;
        if (kspec.input_kind == InputKind::Covariate) {
          input_name = kspec.kernel_input_name;
        } else if (kspec.input_kind == InputKind::ParInput) {
          input_name = kspec.par_input_name;
        } else if (kspec.input_kind == InputKind::Combined) {
          input_name = "combined";
        } else {
          input_name = "noinput";
        }

        for (int code : codes) {
          if (!k_rt.kernel_ptr->has_output_stream(code)) continue;

          KernelOutput ko = k_rt.kernel_ptr->get_output_stream(code);

          if (ko.n_rows != n) {
            stop("TrendRuntime::all_kernel_outputs('%s'): "
                   "stream rows (%d) != n_trials (%d)",
                   spec.target_param.c_str(), ko.n_rows, n);
          }

          std::string suffix = k_rt.kernel_ptr->output_stream_name(code);

          //
          for (int c = 0; c < ko.n_cols; ++c) {
            const double* src = ko.data + c * ko.n_rows;
            for (int r = 0; r < n; ++r) {
              out(r, col) = src[r];
            }

            std::string cname;
            if (ko.n_cols > 1) {
              // Variadic kernel: name is target_param.covariate_name.stream
              // e.g. "m.covariate1.Qmatrix", "m.covariate2.Qmatrix"
              // input_name (the joined "+"-string) is dropped — redundant and verbose.
              std::string cov_name = (c < (int)kspec.covariate_names.size())
              ? Rcpp::as<std::string>(kspec.covariate_names[c])
                : std::to_string(c + 1);
              cname = spec.target_param + "." + cov_name + "." + suffix;
            } else {
              // Single-covariate kernel: target_param.input_name.stream
              // e.g. "m.covariate1.Q"
              cname = spec.target_param + "." + input_name + "." + suffix;
            }

            cn[col] = cname;
            ++col;
          }
        }
      }
    }
  };

  fill_for_ops(premap_ops);
  fill_for_ops(pretransform_ops);
  fill_for_ops(posttransform_ops);

  colnames(out) = cn;
  return out;
}


Rcpp::NumericMatrix TrendRuntime::all_kernel_outputs(ParamTable& pt) {
  return all_kernel_outputs(pt, std::vector<int>{1}); // only main trajectories
}

// Rcpp::NumericMatrix TrendRuntime::all_kernel_outputs(ParamTable& pt) {
//   using namespace Rcpp;
//
//   const int n = pt.n_trials;
//
//   // Count total number of kernel slots
//   int n_slots = 0;
//   for (const auto& op : premap_ops)        n_slots += op.kernels.size();
//   for (const auto& op : pretransform_ops)  n_slots += op.kernels.size();
//   for (const auto& op : posttransform_ops) n_slots += op.kernels.size();
//
//   NumericMatrix out(n, n_slots);
//   CharacterVector cn(n_slots);
//
//   int col = 0;
//
//   auto fill_for_ops = [&](std::vector<TrendOpRuntime>& ops) {
//     for (auto& op : ops) {
//       const TrendOpSpec& spec = *op.spec;
//
//       for (auto& k_rt : op.kernels) {
//         const KernelSlotSpec& kspec = *k_rt.spec;
//         const std::vector<double>& traj = k_rt.kernel_ptr->get_output();
//
//         if ((int)traj.size() != n) {
//           stop("TrendRuntime::all_kernel_outputs('%s'): trajectory length (%d) != n_trials (%d)",
//                spec.target_param.c_str(), (int)traj.size(), n);
//         }
//
//         for (int r = 0; r < n; ++r) {
//           out(r, col) = traj[r];
//         }
//
//         // Build a column name: target_param + "." + input_name
//         std::string input_name;
//         if (kspec.input_kind == InputKind::Covariate) {
//           // try to get column name from attributes if available
//           if (kspec.kernel_input.hasAttribute("names")) {
//             // optional; often covariate is directly from data[cov_name],
//             // so we don't have the name here. You can store the cov_name in KernelSlotSpec if needed.
//             input_name = "cov";
//           } else {
//             input_name = "cov";
//           }
//         } else if (kspec.input_kind == InputKind::ParInput) {
//           input_name = kspec.par_input_name;
//         } else {
//           input_name = "noinput";
//         }
//
//         std::string cname = spec.target_param + "." + input_name;
//         cn[col] = cname;
//
//         ++col;
//       }
//     }
//   };
//
//   fill_for_ops(premap_ops);
//   fill_for_ops(pretransform_ops);
//   fill_for_ops(posttransform_ops);
//
//   colnames(out) = cn;
//   return out;
// }
