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

  if (Rf_isString(cov_spec)) {
    // covariate specified by column name
    string cov_name = as<string>(cov_spec);
    if (!data.containsElementNamed(cov_name.c_str())) {
      stop("TrendPlan: data has no column named '%s' for covariate",
           cov_name.c_str());
    }
    op.covariate = data[cov_name];
  } else {
    // inline numeric vector
    op.covariate = as<NumericVector>(cov_spec);
  }

  if (op.covariate.size() != data.nrows()) {
    stop("TrendPlan: covariate for target '%s' has length %d but data has %d rows",
         op.target_param.c_str(),
         op.covariate.size(), data.nrows());
  }
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

  if (cur.containsElementNamed("par_input")) {
    par_input = cur["par_input"];
  } else {
    par_input = Rcpp::CharacterVector(0);
  }

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

TrendPlan::TrendPlan(Rcpp::Nullable<Rcpp::List> trend_, const Rcpp::DataFrame& data_)
  : data(data_) {

  if (trend_.isNull()) {
    trend = R_NilValue;
    return;
  }

  trend = Rcpp::List(trend_.get());

  using std::string;
  using namespace Rcpp;

  CharacterVector tnames = trend.names();
  if (tnames.size() != trend.size()) {
    stop("TrendPlan: 'trend' must be a named list");
  }

  premap_ops.clear();
  pretransform_ops.clear();
  posttransform_ops.clear();
  premap_trend_params.clear();
  all_trend_params.clear();

  for (int i = 0; i < trend.size(); ++i) {
    if (trend[i] == R_NilValue) continue;

    Rcpp::List tr_i = trend[i];
    std::string name_i = Rcpp::as<std::string>(tnames[i]);

    // Build once to inspect phase, trend_pnames, etc.
    TrendOpSpec tmp(tr_i, name_i);

    // Extract covariate names (if the spec uses multiple columns)
    std::vector<std::string> cov_names = covariate_names_from_spec(tr_i);

    // Collect trend pnames
    for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
      std::string pn = Rcpp::as<std::string>(tmp.trend_pnames[k]);
      all_trend_params.insert(pn);
    }

    // Register trend_pnames for premap
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

    // 0 or 1 covariate name: keep spec as-is
    if (cov_names.size() <= 1) {
      switch (tmp.phase) {
      case TrendPhase::Premap:
        premap_ops.emplace_back(tr_i, name_i);
        init_covariate(premap_ops.back(), data);
        break;
      case TrendPhase::Pretransform:
        pretransform_ops.emplace_back(tr_i, name_i);
        init_covariate(pretransform_ops.back(), data);
        break;
      case TrendPhase::Posttransform:
        posttransform_ops.emplace_back(tr_i, name_i);
        init_covariate(posttransform_ops.back(), data);
        break;
      }
      continue;
    }

    // Multiple covariate columns: one TrendOpSpec per column
    for (const std::string& cov_name : cov_names) {
      Rcpp::List tr_copy = Rcpp::clone(tr_i);
      tr_copy["covariate"] = cov_name;

      switch (tmp.phase) {
      case TrendPhase::Premap:
        premap_ops.emplace_back(tr_copy, name_i);
        init_covariate(premap_ops.back(), data);
        break;
      case TrendPhase::Pretransform:
        pretransform_ops.emplace_back(tr_copy, name_i);
        init_covariate(pretransform_ops.back(), data);
        break;
      case TrendPhase::Posttransform:
        posttransform_ops.emplace_back(tr_copy, name_i);
        init_covariate(posttransform_ops.back(), data);
        break;
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
  op.kernel_par_indices.clear();

  int kernel_start = 0;

  if (spec.base_type == "lin" || spec.base_type == "exp_lin") {
    if (n_tp < 1) {
      Rcpp::stop("TrendOp '%s': base = '%s' but trend_pnames is empty",
                 spec.target_param.c_str(), spec.base_type.c_str());
    }
    std::string base_nm = Rcpp::as<std::string>(spec.trend_pnames[0]);
    op.base_par_idx = pt.base_index_for(base_nm);
    kernel_start = 1;
  }

  for (int k = kernel_start; k < n_tp; ++k) {
    std::string kn = Rcpp::as<std::string>(spec.trend_pnames[k]);
    int idx   = pt.base_index_for(kn);
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
}

void TrendRuntime::run_kernel_for_op(TrendOpRuntime& op,
                                     ParamTable& pt)
{
  using namespace Rcpp;
  using std::string;

  const TrendOpSpec& spec = *op.spec;

  // 1) Covariate from spec + plan->data
  if (!spec.spec.containsElementNamed("covariate")) {
    stop("TrendOp for target '%s' is missing 'covariate' field",
         spec.target_param.c_str());
  }

  // Look-up covariate
  NumericVector cov = spec.covariate;
  const int n = cov.size();
  if (n != pt.n_trials) {
    stop("TrendRuntime::run_kernel_for_op: covariate length (%d) != n_trials (%d)",
         n, pt.n_trials);
  }

  // kernel_par_indices must be initialized
  // actually not needed -- linear kernels have no params...
  // if (op.kernel_par_indices.empty()) {
  //   stop("TrendRuntime::run_kernel_for_op('%s'): kernel_par_indices is empty; "
  //          "did you call bind_all_ops_to_paramtable()?",
  //          spec.target_param.c_str());
  // }

  const int p = pt.base.ncol();
  for (int idx : op.kernel_par_indices) {
    if (idx < 0 || idx >= p) {
      stop("TrendRuntime::run_kernel_for_op('%s'): kernel_par index %d out of [0,%d)",
           spec.target_param.c_str(), idx, p);
    }
  }

  KernelParsView kp_view = make_kernel_pars_view(pt, op.kernel_par_indices);

  op.kernel_ptr->reset();
  op.kernel_ptr->run(kp_view, cov, spec.comp_index);

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
  const double* trend_col = nullptr;
  bool needs_trend_par =
    (base == "lin" || base == "exp_lin" || base == "lin_exp" || base == "centered");

  if (needs_trend_par) {
    if (op.base_par_idx < 0) {
      stop("TrendOp '%s': base_type = '%s' but base_par_idx < 0 "
             "(did you bind_all_ops_to_paramtable() and set trend_pnames?)",
             spec.target_param.c_str(), base.c_str());
    }
    trend_col = &pt.base(0, op.base_par_idx);
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

    double contrib = q;

    if (base == "lin" || base == "exp_lin" || base == "lin_exp") {
      double tp = trend_col ? trend_col[r] : 1.0;
      contrib *= tp;
    } else if (base == "centered") {
      double tp = trend_col ? trend_col[r] : 1.0;
      contrib = (contrib - 0.5) * tp;
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

      // --- build name: target_param.covariate_name ---
      std::string cov_name;
      if (spec.spec.containsElementNamed("covariate")) {
        SEXP cov_spec = spec.spec["covariate"];

        // covariate given as column name
        if (Rf_isString(cov_spec) && Rf_length(cov_spec) == 1) {
          cov_name = Rcpp::as<std::string>(cov_spec);
        } else {
          // inline numeric vector or something else
          cov_name = "inline";
        }
      } else {
        cov_name = "nocov";
      }

      std::string cname = spec.target_param + "." + cov_name;      cn[col] = cname;

      ++col;
    }
  };

  fill_for_ops(premap_ops,        "premap");
  fill_for_ops(pretransform_ops,  "pretransform");
  fill_for_ops(posttransform_ops, "posttransform");

  colnames(out) = cn;
  return out;
}
