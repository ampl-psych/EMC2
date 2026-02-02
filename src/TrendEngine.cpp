#include "TrendEngine.h"
#include <Rcpp.h>
using namespace Rcpp;


void TrendEngine::run_kernel_for_op(TrendOp& op,
                                    ParamTable& pt)
{
  using namespace Rcpp;
  using std::string;

  // 1) Covariate from op.spec + this->data
  if (!op.spec.containsElementNamed("covariate")) {
    stop("TrendOp for target '%s' is missing 'covariate' field",
         op.target_param.c_str());
  }

  NumericVector cov;
  {
    SEXP cov_spec = op.spec["covariate"];
    if (Rf_isString(cov_spec)) {
      string cov_name = as<string>(cov_spec);
      if (!data.containsElementNamed(cov_name.c_str())) {
        stop("TrendEngine: data has no column named '%s' for covariate",
             cov_name.c_str());
      }
      cov = data[cov_name];
    } else {
      cov = as<NumericVector>(cov_spec);
    }
  }

  const int n = cov.size();
  if (n != pt.n_trials) {
    stop("TrendEngine::run_kernel_for_op: covariate length (%d) != n_trials (%d)",
         n, pt.n_trials);
  }

  // kernel_par_indices must be initialized
  if (op.kernel_par_indices.empty()) {
    stop("TrendEngine::run_kernel_for_op('%s'): kernel_par_indices is empty; "
           "did you call bind_all_ops_to_paramtable()?",
           op.target_param.c_str());
  }

  // Optional: bounds check
  const int p = pt.base.ncol();
  for (int idx : op.kernel_par_indices) {
    if (idx < 0 || idx >= p) {
      stop("TrendEngine::run_kernel_for_op('%s'): kernel_par index %d out of [0,%d)",
           op.target_param.c_str(), idx, p);
    }
  }

  // 2) kernel_pars view directly on ParamTable::base
  KernelParsView kp_view = make_kernel_pars_view(pt, op.kernel_par_indices);

  // 3) comp_idx (no compression yet)
  IntegerVector comp_idx(n);
  for (int i = 0; i < n; ++i) comp_idx[i] = i;

  // 4) Run the kernel
  op.kernel_ptr->reset();
  op.kernel_ptr->run(kp_view, cov, comp_idx);
}

void TrendEngine::apply_base_for_op(TrendOp& op,
                                    ParamTable& pt)
{
  using namespace Rcpp;

  // Ensure kernel has been run
  if (!op.kernel_ptr->has_run()) {
    run_kernel_for_op(op, pt);
  }

  const std::vector<double>& traj = op.kernel_ptr->get_output();
  const int n = pt.n_trials;
  if ((int)traj.size() != n) {
    stop("TrendEngine::apply_base_for_op('%s'): trajectory length (%d) != n_trials (%d)",
         op.target_param.c_str(), (int)traj.size(), n);
  }

  const std::string& base = op.base_type;

  // Target parameter column
  int target_idx = pt.base_index_for(op.target_param);
  double* target_col = &pt.base(0, target_idx);  // will be updated in-place

  // Base/trend parameter if needed (trend_pars in original code)
  const double* trend_col = nullptr;
  bool needs_trend_par =
    (base == "lin" || base == "exp_lin" || base == "lin_exp" || base == "centered");

  if (needs_trend_par) {
    if (op.base_par_idx < 0) {
      stop("TrendOp '%s': base_type = '%s' but base_par_idx < 0 "
             "(did you bind_all_ops_to_paramtable() and set trend_pnames?)",
             op.target_param.c_str(), base.c_str());
    }
    trend_col = &pt.base(0, op.base_par_idx);
  }

  for (int r = 0; r < n; ++r) {
    double p = target_col[r];  // old param[r]
    if (NumericVector::is_na(p) || std::isnan(p)) {
      // If you prefer to treat NA as 0, change this
      continue;
    }

    double q = traj[r];        // kernel_out(r, 0)
    if (NumericVector::is_na(q) || std::isnan(q)) {
      // matches kernel behavior: NA -> no contribution
      q = 0.0;
    }

    // --- Compute contrib as in original inner loop ---
    double contrib = q;

    if (base == "lin" || base == "exp_lin" || base == "lin_exp") {
      double tp = trend_col ? trend_col[r] : 1.0;
      contrib *= tp;  // contrib *= trend_pars(r, map_n)
    } else if (base == "centered") {
      double tp = trend_col ? trend_col[r] : 1.0;
      contrib = (contrib - 0.5) * tp;
    }
    // base == "add": contrib unchanged

    // --- Add parameter term as in the final block of original code ---
    double base_term;
    if (base == "exp_lin" || base == "lin_exp") {
      base_term = std::exp(p);
    } else if (base == "lin" || base == "centered" || base == "add" || base.empty()) {
      base_term = p;
    } else {
      // default: treat unknown base types as linear
      base_term = p;
    }

    // out[r] = contrib + base_term; write back to target in-place
    target_col[r] = contrib + base_term;
  }
}

void TrendEngine::initialize_from_trend() {
  using namespace Rcpp;
  using std::string;

  CharacterVector tnames = trend.names();
  if (tnames.size() != trend.size()) {
    stop("TrendEngine: 'trend' must be a named list");
  }

  premap_ops.clear();
  pretransform_ops.clear();
  posttransform_ops.clear();
  premap_trend_params.clear();

  for (int i = 0; i < trend.size(); ++i) {
    if (trend[i] == R_NilValue) continue;

    Rcpp::List tr = trend[i];
    std::string name_i = Rcpp::as<std::string>(tnames[i]);

    // build once to inspect phase and trend_pnames
    TrendOp tmp(tr, name_i);

    switch (tmp.phase) {
    case TrendPhase::Premap: {
      // emplace a fresh TrendOp (no copy, direct construction)
      premap_ops.emplace_back(tr, name_i);
      // collect pnames from tmp
      for (int k = 0; k < tmp.trend_pnames.size(); ++k) {
        premap_trend_params.insert(Rcpp::as<std::string>(tmp.trend_pnames[k]));
      }
      break;
    }
    case TrendPhase::Pretransform:
      pretransform_ops.emplace_back(tr, name_i);
      break;
    case TrendPhase::Posttransform:
      posttransform_ops.emplace_back(tr, name_i);
      break;
    }
  }
}


bool TrendEngine::has_premap() const       { return !premap_ops.empty(); }
bool TrendEngine::has_pretransform() const { return !pretransform_ops.empty(); }
bool TrendEngine::has_posttransform() const{ return !posttransform_ops.empty(); }

// LogicalVector over designs: TRUE if names(designs)[i] is a trend_pname
// of any *premap* trend.
Rcpp::LogicalVector TrendEngine::premap_design_mask(const Rcpp::List& designs) const {
  using namespace Rcpp;
  using std::string;

  const int n = designs.size();
  LogicalVector mask(n, false);

  CharacterVector dnames = designs.names();
  if (dnames.size() != n) {
    stop("TrendEngine::premap_design_mask: 'designs' must be a named list");
  }

  for (int i = 0; i < n; ++i) {
    string nm = as<string>(dnames[i]);
    mask[i] = (premap_trend_params.find(nm) != premap_trend_params.end());
  }

  return mask;
}

// Call this once you have a ParamTable for this run
void TrendEngine::bind_all_ops_to_paramtable(const ParamTable& pt) {
  for (auto& op : premap_ops)        bind_trend_op_to_paramtable(op, pt);
  for (auto& op : pretransform_ops)  bind_trend_op_to_paramtable(op, pt);
  for (auto& op : posttransform_ops) bind_trend_op_to_paramtable(op, pt);
}

