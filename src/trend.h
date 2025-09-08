#ifndef dynamic_h
#define dynamic_h

#include "utility_functions.h"
#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;


NumericVector run_delta_rcpp(NumericVector q0, NumericVector alpha, NumericVector covariate) {
  int n = covariate.length();
  NumericVector q(n);
  NumericVector pe(n);
  q[0] = q0[0];
  for(int i = 1; i < n; i++) {
    pe[i-1] = covariate[i-1] - q[i-1];
    q[i] = q[i-1] + alpha[i-1] * pe[i-1];
  }

  return q;
}

NumericVector run_delta2_rcpp(NumericVector q0, NumericVector alphaFast,
                              NumericVector propSlow, NumericVector dSwitch,
                              NumericVector covariate) {
  int n = covariate.length();
  NumericVector q(n);
  NumericVector qFast(n);
  NumericVector qSlow(n);
  NumericVector peFast(n);
  NumericVector peSlow(n);

  q[0] = qFast[0] = qSlow[0] = q0[0];
  NumericVector alphaSlow = propSlow * alphaFast;

  for(int i = 1; i < n; i++) {
    peFast[i-1] = covariate[i-1] - qFast[i-1];
    peSlow[i-1] = covariate[i-1] - qSlow[i-1];

    qFast[i] = qFast[i-1] + alphaFast[i-1] * peFast[i-1];
    qSlow[i] = qSlow[i-1] + alphaSlow[i-1] * peSlow[i-1];

    if(std::abs(qFast[i] - qSlow[i]) > dSwitch[i]) {
      q[i] = qFast[i];
    } else {
      q[i] = qSlow[i];
    }
  }

  return q;
}


NumericVector run_kernel_rcpp(NumericMatrix trend_pars, String kernel, NumericMatrix input, int n_base_pars) {
  // Kernels accept any number of input columns; apply per column and sum.
  const int n = input.nrow();
  const int p = input.ncol();
  NumericVector out(n, 0.0);

  for (int c = 0; c < p; ++c) {
    NumericVector covariate = input(_, c);
    NumericVector contrib(n);

    if(kernel == "lin_decr") {
      contrib = -1 * covariate;
    }
    else if(kernel == "lin_incr") {
      contrib = covariate;
    }
    else if(kernel == "exp_decr") {
      contrib = exp(-trend_pars(_, 0 + n_base_pars) * covariate);
    }
    else if(kernel == "exp_incr") {
      contrib = 1 - exp(-trend_pars(_, 0 + n_base_pars) * covariate);
    }
    else if(kernel == "pow_decr") {
      contrib = vector_pow(1 + covariate, -trend_pars(_, 0 + n_base_pars));
    }
    else if(kernel == "pow_incr") {
      contrib = 1 - vector_pow(1 + covariate, -trend_pars(_, 0 + n_base_pars));
    }
    else if(kernel == "poly2") {
      contrib = trend_pars(_, 0 + n_base_pars) * covariate + trend_pars(_, 1 + n_base_pars) * pow(covariate, 2);
    }
    else if(kernel == "poly3") {
      contrib = trend_pars(_, 0 + n_base_pars) * covariate + trend_pars(_, 1 + n_base_pars) * pow(covariate, 2) +
        trend_pars(_, 2 + n_base_pars) * pow(covariate, 3);
    }
    else if(kernel == "poly4") {
      contrib = trend_pars(_, 0 + n_base_pars) * covariate + trend_pars(_, 1 + n_base_pars) * pow(covariate, 2) +
        trend_pars(_, 2 + n_base_pars) * pow(covariate, 3) + trend_pars(_, 3 + n_base_pars) * pow(covariate, 4);
    }
    else if(kernel == "delta") {
      // For sequential kernels, ignore rows where the covariate is NA by operating on compressed series
      LogicalVector good = !is_na(covariate);
      NumericVector cov_tmp = covariate[good];
      NumericMatrix pars_tmp = submat_rcpp(trend_pars, good);
      NumericVector tmp_out = run_delta_rcpp(pars_tmp(_, 0 + n_base_pars), pars_tmp(_, 1 + n_base_pars), cov_tmp);
      // expand back
      int pos = 0;
      for (int i = 0; i < n; ++i) {
        if (good[i]) {
          contrib[i] = tmp_out[pos++];
        } else {
          contrib[i] = 0.0; // ignore NAs
        }
      }
    }
    else if(kernel == "delta2") {
      LogicalVector good = !is_na(covariate);
      NumericVector cov_tmp = covariate[good];
      NumericMatrix pars_tmp = submat_rcpp(trend_pars, good);
      NumericVector tmp_out = run_delta2_rcpp(pars_tmp(_, 0 + n_base_pars), pars_tmp(_, 1 + n_base_pars), pars_tmp(_, 2 + n_base_pars),
                            pars_tmp(_, 3 + n_base_pars), cov_tmp);
      int pos = 0;
      for (int i = 0; i < n; ++i) {
        if (good[i]) {
          contrib[i] = tmp_out[pos++];
        } else {
          contrib[i] = 0.0;
        }
      }
    } else {
      // Unknown kernel
      stop("Unknown kernel type");
    }

    // Replace NA by 0 in contrib (ignore NA values)
    for (int i = 0; i < n; ++i) {
      if (NumericVector::is_na(contrib[i])) contrib[i] = 0.0;
      out[i] += contrib[i];
    }
  }

  return out;
}

// Now accepts the full parameter matrix `pars_full` so we can use par_input columns as inputs too.
// Passes all inputs (covariates + par_input) to kernel in one call; kernel sums across columns.
NumericVector run_trend_rcpp(DataFrame data, List trend, NumericVector param, NumericMatrix trend_pars, NumericMatrix pars_full) {
  String kernel = as<String>(trend["kernel"]);
  String base = as<String>(trend["base"]);
  CharacterVector covnames;
  if (trend.containsElementNamed("covariate") && !Rf_isNull(trend["covariate"])) {
    covnames = trend["covariate"];
  } else {
    covnames = CharacterVector(0);
  }
  CharacterVector par_input;
  if (trend.containsElementNamed("par_input") && !Rf_isNull(trend["par_input"])) {
    par_input = trend["par_input"];
  } else {
    par_input = CharacterVector(0);
  }
  // Initialize output vector with zeros
  int n_trials = param.length();
  NumericVector out(n_trials, 0.0);
  int n_base_pars = 0;
  if(base == "lin" || base == "exp_lin" || base == "centered") {
    n_base_pars = 1;
  }
  // Build input matrix from covariates and par_input columns
  int n_cov = covnames.size();
  // Keep only par_input columns that actually exist in pars_full (if provided)
  CharacterVector pars_full_names = colnames(pars_full);
  std::vector<std::string> par_in_keep;
  for (int i = 0; i < par_input.size(); i++) {
    std::string nm = Rcpp::as<std::string>(par_input[i]);
    bool found = false;
    for (int j = 0; j < pars_full_names.size(); j++) {
      if (nm == Rcpp::as<std::string>(pars_full_names[j])) { found = true; break; }
    }
    if (found) par_in_keep.push_back(nm);
  }
  int n_par_in = (int)par_in_keep.size();

  int n_inputs = n_cov + n_par_in;
  // If no inputs provided, keep out as zeros (base will handle accordingly)
  if (n_inputs > 0) {
    NumericMatrix input_all(n_trials, n_inputs);
    // Fill covariate columns first
    for (int i = 0; i < n_cov; i++) {
      String cur_cov = covnames[i];
      NumericVector covariate = as<NumericVector>(data[cur_cov]);
      input_all(_, i) = covariate;
    }
    // Then par_input columns
    if (n_par_in > 0) {
      CharacterVector par_in_keep_cv(n_par_in);
      for (int i = 0; i < n_par_in; i++) par_in_keep_cv[i] = par_in_keep[i];
      NumericMatrix pin = submat_rcpp_col_by_names(pars_full, par_in_keep_cv);
      for (int j = 0; j < n_par_in; j++) {
        input_all(_, n_cov + j) = pin(_, j);
      }
    }

    // Call kernel once with all inputs; kernel sums across columns and ignores NA
    NumericVector kernel_out = run_kernel_rcpp(trend_pars, kernel, input_all, n_base_pars);

    out = out + kernel_out;
  }
  // Apply base transformation to final summed output
  if(base == "lin") {
    out = param + trend_pars(_, 0) * out;
  }
  else if(base == "exp_lin") {
    out = exp(param) + trend_pars(_, 0) * out;
  }
  else if(base == "centered") {
    out = param + trend_pars(_, 0) * (out - 0.5);
  }
  else if(base == "add") {
    out = param + out;
  }
  return out;
}

// A few unneccessary loops in here, but seems reasonably efficient
inline NumericMatrix prep_trend_phase(DataFrame data, List trend, NumericMatrix pars, String phase) {
  CharacterVector trend_names = trend.names();
  CharacterVector par_names = colnames(pars);
  CharacterVector all_remove;
  for (int i = 0; i < trend.size(); ++i) {
    List cur_trend = trend[i];
    std::string cur_ph = Rcpp::as<std::string>(cur_trend["phase"]);
    std::string phase_s = Rcpp::as<std::string>(Rcpp::wrap(phase));
    if (cur_ph != phase_s) continue;
    String par = trend_names[i];
    CharacterVector trend_pnames = cur_trend["trend_pnames"];
    all_remove = c_add_charvectors(all_remove, trend_pnames);
    LogicalVector par_idx = contains(par_names, par);
    NumericVector param = as<NumericVector>(submat_rcpp_col(pars, par_idx));
    NumericMatrix trend_pars = submat_rcpp_col_by_names(pars, trend_pnames);
    pars(_, as<int>(which_rcpp(par_idx))) = run_trend_rcpp(data, cur_trend, param, trend_pars, pars);
  }
  all_remove = unique(all_remove);
  if (all_remove.size() > 0) {
    CharacterVector pnames = colnames(pars);
    LogicalVector keep = !contains_multiple(pnames, all_remove);
    pars = submat_rcpp_col(pars, keep);
  }
  return(pars);
}

inline NumericMatrix prep_trend_phase_with_pars(DataFrame data, List trend, NumericMatrix pars, String phase, const NumericMatrix& trend_pars_override) {
  CharacterVector trend_names = trend.names();
  CharacterVector par_names = colnames(pars);
  CharacterVector all_remove;
  for (int i = 0; i < trend.size(); ++i) {
    List cur_trend = trend[i];
    // Simple cast to std::string for comparison
    std::string cur_ph = Rcpp::as<std::string>(cur_trend["phase"]);
    std::string phase_s = Rcpp::as<std::string>(Rcpp::wrap(phase));
    if (cur_ph != phase_s) continue;
    String par = trend_names[i];
    CharacterVector trend_pnames = cur_trend["trend_pnames"];
    all_remove = c_add_charvectors(all_remove, trend_pnames);
    LogicalVector par_idx = contains(par_names, par);
    NumericVector param = as<NumericVector>(submat_rcpp_col(pars, par_idx));
    NumericMatrix cur_tp = submat_rcpp_col_by_names(trend_pars_override, trend_pnames);
    pars(_, as<int>(which_rcpp(par_idx))) = run_trend_rcpp(data, cur_trend, param, cur_tp, pars);
  }
  all_remove = unique(all_remove);
  if (all_remove.size() > 0) {
    CharacterVector pnames = colnames(pars);
    LogicalVector keep = !contains_multiple(pnames, all_remove);
    pars = submat_rcpp_col(pars, keep);
  }
  return(pars);
}

// ---- Trend helpers for mapping pipeline ----

// Collect all unique trend parameter names across trend list entries
inline CharacterVector collect_trend_param_names(const List& trend) {
  CharacterVector trend_pnames;
  for (int i = 0; i < trend.size(); ++i) {
    List cur_trend = trend[i];
    CharacterVector cur_names = cur_trend["trend_pnames"];
    trend_pnames = c_add_charvectors(trend_pnames, cur_names);
  }
  return unique(trend_pnames);
}

inline CharacterVector collect_trend_param_names_phase(const List& trend, const std::string& phase) {
  CharacterVector trend_pnames;
  for (int i = 0; i < trend.size(); ++i) {
    List cur_trend = trend[i];
    std::string ph = Rcpp::as<std::string>(cur_trend["phase"]);
    if (ph != phase) continue;
    CharacterVector cur_names = cur_trend["trend_pnames"];
    trend_pnames = c_add_charvectors(trend_pnames, cur_names);
  }
  return unique(trend_pnames);
}

// Build per-trial columns for trend parameters from designs and p_vector, and apply transforms
inline NumericMatrix build_trend_columns_from_design(NumericVector p_vector,
                                                     CharacterVector p_types,
                                                     List designs,
                                                     int n_trials,
                                                     const List& trend,
                                                     const std::vector<TransformSpec>& full_specs) {
  CharacterVector trend_pnames = collect_trend_param_names(trend);
  if (trend_pnames.size() == 0) {
    return NumericMatrix(n_trials, 0); // empty
  }

  // Map parameter name -> index in p_types/designs
  std::unordered_map<std::string,int> name_to_idx;
  for (int i = 0; i < p_types.size(); ++i) {
    name_to_idx[ Rcpp::as<std::string>(p_types[i]) ] = i;
  }

  // Allocate output matrix with columns strictly in trend_pnames order
  NumericMatrix trend_pars(n_trials, trend_pnames.size());
  colnames(trend_pars) = trend_pnames;

  for (int c = 0; c < trend_pnames.size(); ++c) {
    std::string pname = Rcpp::as<std::string>(trend_pnames[c]);
    auto it = name_to_idx.find(pname);
    int idx = it->second;
    NumericMatrix cur_design = designs[idx];
    CharacterVector cur_names = colnames(cur_design);

    // Accumulate p_vector * design columns
    NumericVector acc(n_trials, 0.0);
    for (int j = 0; j < cur_design.ncol(); ++j) {
      String cur_name(cur_names[j]);
      NumericVector tmp = p_vector[cur_name] * cur_design(_, j);
      LogicalVector bad = is_na(tmp) | is_nan(tmp);
      tmp[bad] = 0;
      acc = acc + tmp;
    }
    trend_pars(_, c) = acc;
  }

  // Transform trend parameter columns using precomputed specs
  std::vector<TransformSpec> t_specs = make_transform_specs_from_full(trend_pars, p_types, full_specs);
  trend_pars = c_do_transform(trend_pars, t_specs);
  return trend_pars;
}

// Apply premap trend for a single parameter vector if a trend is defined
inline NumericVector apply_premap_trends(const DataFrame& data,
                                         const List& trend,
                                         const CharacterVector& trend_names,
                                         const String& param_name,
                                         NumericVector param_values,
                                         const NumericMatrix& trend_pars) {
  // Apply all premap trend entries for this param_name sequentially, in order
  NumericVector result = clone(param_values);
  for (int i = 0; i < trend.size(); ++i) {
    if (trend_names[i] == param_name) {
      List cur_trend = trend[i];
      std::string ph = Rcpp::as<std::string>(cur_trend["phase"]);
      if (ph != "premap") continue;
      CharacterVector cur_trend_pnames = cur_trend["trend_pnames"];
      NumericMatrix cur_trend_pars = submat_rcpp_col_by_names(trend_pars, cur_trend_pnames);
      NumericMatrix empty_full(result.size(), 0);
      result = run_trend_rcpp(data, cur_trend, result, cur_trend_pars, empty_full);
    }
  }
  return result;
}

// Fill columns for trend parameters into the output matrix by name (for pretransform case)
inline void fill_trend_columns_for_pretransform(NumericMatrix& pars,
                                                const CharacterVector& p_types,
                                                const NumericMatrix& trend_pars) {
  if (trend_pars.ncol() == 0) return;
  CharacterVector tnames = colnames(trend_pars);
  // Build a name->col map for trend_pars
  std::unordered_map<std::string,int> tmap;
  for (int c = 0; c < tnames.size(); ++c) {
    tmap[Rcpp::as<std::string>(tnames[c])] = c;
  }
  for (int i = 0; i < p_types.size(); ++i) {
    std::string pname = Rcpp::as<std::string>(p_types[i]);
    auto it = tmap.find(pname);
    if (it != tmap.end()) {
      pars(_, i) = trend_pars(_, it->second);
    }
  }
}

#endif
