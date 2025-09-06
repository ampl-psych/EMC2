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


NumericVector run_kernel_rcpp(NumericMatrix trend_pars, String kernel, NumericVector covariate, int n_base_pars) {
  NumericVector out(covariate.length());
  if(kernel == "lin_decr") {
    out = -1 * covariate;
  }
  else if(kernel == "lin_incr") {
    out = covariate;
  }
  else if(kernel == "exp_decr") {
    out = exp(-trend_pars(_, 0 + n_base_pars) * covariate);
  }
  else if(kernel == "exp_incr") {
    out = 1 - exp(-trend_pars(_, 0 + n_base_pars) * covariate);
  }
  else if(kernel == "pow_decr") {
    out = vector_pow(1 + covariate, -trend_pars(_, 0 + n_base_pars));
  }
  else if(kernel == "pow_incr") {
    out = 1 - vector_pow(1 + covariate, -trend_pars(_, 0 + n_base_pars));
  }
  else if(kernel == "poly2") {
    out = trend_pars(_, 0 + n_base_pars) * covariate + trend_pars(_, 1 + n_base_pars) * pow(covariate, 2);
  }
  else if(kernel == "poly3") {
    out = trend_pars(_, 0 + n_base_pars) * covariate + trend_pars(_, 1 + n_base_pars) * pow(covariate, 2) +
      trend_pars(_, 2 + n_base_pars) * pow(covariate, 3);
  }
  else if(kernel == "poly4") {
    out = trend_pars(_, 0 + n_base_pars) * covariate + trend_pars(_, 1 + n_base_pars) * pow(covariate, 2) +
      trend_pars(_, 2 + n_base_pars) * pow(covariate, 3) + trend_pars(_, 3 + n_base_pars) * pow(covariate, 4);
  }
  else if(kernel == "delta") {
    out = run_delta_rcpp(trend_pars(_, 0 + n_base_pars), trend_pars(_, 1 + n_base_pars), covariate);
  }
  else if(kernel == "delta2") {
    out = run_delta2_rcpp(trend_pars(_, 0 + n_base_pars), trend_pars(_, 1 + n_base_pars), trend_pars(_, 2 + n_base_pars),
                          trend_pars(_, 3 + n_base_pars), covariate);
  }

  return out;
}

NumericVector run_trend_rcpp(DataFrame data, List trend, NumericVector param, NumericMatrix trend_pars) {
  String kernel = as<String>(trend["kernel"]);
  String base = as<String>(trend["base"]);
  CharacterVector covnames = trend["covariate"];
  // Initialize output vector with zeros
  int n_trials = param.length();
  NumericVector out(n_trials);
  int n_base_pars = 0;
  if(base == "lin" || base == "exp_lin" || base == "centered") {
    n_base_pars = 1;
  }
  // Loop through covariates
  for(int i = 0; i < covnames.length(); i++) {
    String cur_cov = covnames[i];
    // Get covariate data and handle NAs
    NumericVector covariate = as<NumericVector>(data[cur_cov]);
    LogicalVector NA_idx = is_na(covariate);
    // Create temporary vectors excluding NAs
    NumericVector cov_tmp = covariate[!NA_idx];
    NumericMatrix trend_pars_tmp = submat_rcpp(trend_pars, !NA_idx);
    // For non-delta kernels, filter duplicates
    LogicalVector filter;
    if(kernel != "delta" && kernel != "delta2") {
      // Create matrix of covariate and trend parameters for duplicate checking
      NumericMatrix together(cov_tmp.length(), trend_pars_tmp.ncol() + 1);
      together(_, 0) = cov_tmp;
      for(int j = 0; j < trend_pars_tmp.ncol(); j++) {
        together(_, j + 1) = trend_pars_tmp(_, j);
      }
      filter = !duplicated_matrix(together);
    } else {
      filter = LogicalVector(cov_tmp.length(), true);
    }
    // Run kernel on unique entries
    NumericVector output = run_kernel_rcpp(submat_rcpp(trend_pars_tmp, filter), kernel, cov_tmp[filter], n_base_pars);
    // // Create index for expanding back to full size
    IntegerVector unq_idx = cumsum_logical(filter); // This function is 1-based
    NumericVector expanded_output = c_expand(output, unq_idx); //This function assumes 1-based as well
    // Add to cumulative output
    for(int k = 0, l = 0; k < n_trials; k ++){
      // put back values, expanded output could be shorter, since NAs are excluded
      // hence the loop.
      if(NA_idx[k] == FALSE){
        out[k] = out[k] + expanded_output[l];
        l++;
      }
    }
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
NumericMatrix prep_trend(DataFrame data, List trend, NumericMatrix pars) {
  // Get parameter names
  CharacterVector trend_names = trend.names();
  CharacterVector all_trend_names;
  CharacterVector par_names = colnames(pars);
  // Apply trends
  for(unsigned int i = 0; i < trend_names.length(); i ++) {
    String par = trend_names[i];
    List cur_trend = trend[par];
    // Get trend parameter names
    CharacterVector trend_pnames = cur_trend["trend_pnames"];
    all_trend_names = c_add_charvectors(all_trend_names, trend_pnames);
    LogicalVector idx = contains_multiple(par_names, trend_pnames);
    LogicalVector par_idx = contains(par_names, par);
    // Extract parameter and trend parameter columns
    NumericVector param = as<NumericVector>(submat_rcpp_col(pars, par_idx));
    NumericMatrix trend_pars = submat_rcpp_col(pars, idx);
    // Run trend
    pars(_, as<int>(which_rcpp(par_idx))) = run_trend_rcpp(data, cur_trend, param, trend_pars);
  }
  all_trend_names = unique(all_trend_names);
  pars = submat_rcpp_col(pars, !contains_multiple(par_names, all_trend_names));
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

// Build per-trial columns for trend parameters from designs and p_vector, and apply transforms
inline NumericMatrix build_trend_columns_from_design(NumericVector p_vector,
                                                     CharacterVector p_types,
                                                     List designs,
                                                     int n_trials,
                                                     const List& trend,
                                                     const List& transforms) {
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

  // Transform trend parameter columns
  std::vector<TransformSpec> t_specs = make_transform_specs(trend_pars, transforms);
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
  // Quick check if a trend exists for param_name
  if (!is_true(any(contains(trend_names, param_name)))) {
    return param_values;
  }
  // Fetch specific trend and its parameter subset
  List cur_trend = trend[param_name];
  CharacterVector cur_trend_pnames = cur_trend["trend_pnames"];
  NumericMatrix cur_trend_pars = submat_rcpp_col_by_names(trend_pars, cur_trend_pnames);
  return run_trend_rcpp(data, cur_trend, param_values, cur_trend_pars);
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
