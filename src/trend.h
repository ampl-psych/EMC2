#ifndef dynamic_h
#define dynamic_h

#include "utility_functions.h"
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector run_deltab_rcpp(NumericVector weight, NumericVector q0, NumericVector alpha, NumericVector b, NumericVector covariate) {
  int n = covariate.length();
  NumericVector q(n);
  NumericVector pe(n);
  q[0] = q0[0];
  b[0] = b[0] + weight[0]*q[0];
  for(int i = 1; i < n; i++) {
    pe[i-1] = b[i-1]/covariate[i-1] - q[i-1];
    q[i] = q[i-1] + alpha[i-1] * pe[i-1];
    b[i] = b[i] + weight[i]*q[i];
  }

  return q;
}


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


NumericVector run_kernel_rcpp(NumericMatrix trend_pars, String kernel, NumericVector covariate, int n_base_pars, NumericVector param) {
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
  else if(kernel == "deltab") {
    // weight is always base_pars 0
    out = run_deltab_rcpp(trend_pars(_, 0), trend_pars(_, 1), trend_pars(_, 2), param, covariate);
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

// [[Rcpp::export]]
NumericVector run_trend_rcpp(DataFrame data, List trend, NumericVector param, NumericMatrix trend_pars) {
  String kernel = as<String>(trend["kernel"]);
  String base = as<String>(trend["base"]);
  CharacterVector covnames = trend["covariate"];
  bool filter_lR = false;
  if(trend.containsElementNamed("filter_lR")) {
    filter_lR = as<bool>(trend["filter_lR"]);
  }
  // Initialize output vector with zeros
  int n_trials = param.length();
  NumericVector out(n_trials);
  int n_base_pars = 0;
  if(base == "lin" || base == "exp_lin" || base == "centered") {
    n_base_pars = 1;
    if(kernel == "deltab") {
      n_base_pars = 0;
    }
  }
  // Loop through covariates
  for(int i = 0; i < covnames.length(); i++) {
    String cur_cov = covnames[i];
    // Get covariate data and handle NAs
    NumericVector covariate = as<NumericVector>(data[cur_cov]);
    LogicalVector NA_idx = is_na(covariate);

    // In delta rules, the first trial should always be 'updated' with the q0-value, even when it is set to NA. Otherwise
    // the first NA-trials will have an updated_covariate value of 0, which is incorrect.
    if(kernel != "delta" && kernel != "delta2" && kernel != "deltab") {
      if(NA_idx[0]) {
        covariate[0] = trend_pars(0,1);
        NA_idx[0] = false;
      }
    }

    // Create temporary vectors excluding NAs
    NumericVector cov_tmp = covariate[!NA_idx];
    NumericMatrix trend_pars_tmp = submat_rcpp(trend_pars, !NA_idx);
    NumericVector param_tmp = param[!NA_idx];
    // For non-delta kernels, filter duplicates
    LogicalVector filter;
    if(kernel != "delta" && kernel != "delta2" && kernel != "deltab") {
      // Create matrix of covariate and trend parameters for duplicate checking
      NumericMatrix together(cov_tmp.length(), trend_pars_tmp.ncol() + 1);
      together(_, 0) = cov_tmp;
      for(int j = 0; j < trend_pars_tmp.ncol(); j++) {
        together(_, j + 1) = trend_pars_tmp(_, j);
      }
      filter = !duplicated_matrix(together);
    } else {
      filter = LogicalVector(cov_tmp.length(), true);

      // // SM: filter only first lR levels in dadm
      // bool filter_lR = as<bool>(trend["filter_lR"]);
      // if(filter_lR) {
      //   IntegerVector lRcol = data["lR"];
      //   for(int tr = 0; tr < lRcol.length(); tr++) {
      //     if(lRcol[tr] != 1) {
      //       filter[tr] = false;
      //     }
      //   }
      // }
    }
    // Run kernel on unique entries
    NumericVector output = run_kernel_rcpp(submat_rcpp(trend_pars_tmp, filter), kernel, cov_tmp[filter], n_base_pars, param_tmp[filter]);
    // // Create index for expanding back to full size
    IntegerVector unq_idx = cumsum_logical(filter); // Is 1-based
    NumericVector expanded_output = c_expand(output, unq_idx); //This assumes 1-based as well
    // Add to cumulative output
    for(int k = 0, l = 0; k < n_trials; k ++){
      // put back values, expanded output could be shorter, since NAs are excluded
      if(NA_idx[k] == FALSE){
        out[k] = out[k] + expanded_output[l];
        l++;
      } else {
        if(filter_lR) {
          // If covariates were filtered on lR, forward fill the covariates
          if(k > 0) {
            out[k] = out[k-1];
          }
        }
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

#endif
