#ifndef dynamic_h
#define dynamic_h

#include "utility_functions.h"
#include <Rcpp.h>
using namespace Rcpp;

NumericVector run_delta_i_dyn(double q0, double alpha, NumericVector target, NumericVector winner){

  NumericVector q(target.length());
  q[0] = q0;
  double alpha_use = alpha; //R::pnorm(alpha, 0, 1, TRUE, FALSE);
  double pe;
  for(int t = 1; t < q.length(); t ++){
    if(NumericVector::is_na(target[t-1]) || winner[t-1] != 1){
      q[t] = q[t-1];
    } else{
      pe = target[t-1] - q[t-1];
      q[t] = q[t-1] + alpha_use*pe;
    }
  }
  return(q);
}

NumericVector run_delta_dyn(double q0, double p, NumericMatrix target, bool isRL){
  NumericVector winner(target.nrow());
  NumericVector s_col;
  LogicalVector tmp = contains(colnames(target), "winner");
  if(isRL){
    if(tmp[0] == TRUE){
      winner = target(_, 0);
      s_col = target(_, 1);
    } else{
      winner = target(_, 1);
      s_col = target(_, 0);
    }
    target = target(_, Range(2, target.ncol() - 1));
  } else{
    winner.fill(1);
  }
  NumericMatrix out(target.nrow(), target.ncol());
  NumericVector tmp2;
  for(int r = 0; r < target.ncol(); r ++){
    out(_, r) = run_delta_i_dyn(q0, p, target(_, r), winner);
  }
  NumericVector out_v(target.nrow());
  if(isRL){
    for(int v = 0; v < target.nrow(); v ++){
      out_v[v] = out(v, s_col[v]-1);
    }
  } else{
    out_v = as<NumericVector>(out);
  }
  return(out_v);
}

// Here param is the parameter, and dp are the additional parameters. They are all doubles.
// kernel and base inform us what transformation to be done.
// Data and covnames, tell us what columns in the data are informing our transformations
NumericVector dynfuns_c(double param, NumericVector dp, String kernel, String base, DataFrame data, CharacterVector covnames,
                        bool do_lR) {
  bool isRL = is_true(any(contains(covnames, "winner")));
  NumericVector out(data.nrow());
  NumericVector lR = data["lR"];
  NumericMatrix covariates(data.nrow(), covnames.length());
  CharacterVector data_names = data.names();
  LogicalVector cnames_idx = contains_multiple(data_names, covnames);
  int k = 0;
  CharacterVector col_names(sum(cnames_idx));
  for(int j = 0; j < data.ncol(); j ++){
    if(cnames_idx[j] == TRUE){
      covariates(_, k) = as<NumericVector>(data[j]);
      col_names[k] = data_names[j];
      k++;
    }
  }
  if(do_lR){
    covariates = submat_rcpp(covariates, lR == 1);
  } else{
    covariates = submat_rcpp(covariates, lR < 9999); // I'm doing something wrong here but this fixes it????
  }
  colnames(covariates) = col_names;
  NumericVector res(covariates.nrow());
  //
  // dyntypes
  if(kernel == "lin_decr"){
    out = -1* as<NumericVector>(covariates);
  }
  if(kernel == "lin_incr"){
    out = as<NumericVector>(covariates);
  }
  if(kernel == "exp_decr"){
    out = exp(dp[1]*as<NumericVector>(covariates));
  }
  if(kernel == "exp_incr"){
    out = 1-exp(-dp[1]*as<NumericVector>(covariates));
  }
  if(kernel == "pow_decr"){
    out = pow((1+as<NumericVector>(covariates)), -dp[1]);
  }
  if(kernel == "pow_incr"){
    out = 1-pow((1+as<NumericVector>(covariates)), -dp[1]);
  }
  if(kernel == "poly2"){
    out = dp[0]*as<NumericVector>(covariates) + dp[1]* pow(as<NumericVector>(covariates), 2);
  }
  if(kernel == "poly3"){
    out = dp[0]*as<NumericVector>(covariates) + dp[1]* pow(as<NumericVector>(covariates), 2) + dp[2]* pow(as<NumericVector>(covariates), 3);
  }
  if(kernel == "poly4"){
    out = dp[0]*as<NumericVector>(covariates) + dp[1]* pow(as<NumericVector>(covariates), 2) + dp[2]* pow(as<NumericVector>(covariates), 3) + dp[3]* pow(as<NumericVector>(covariates), 4);
  }
  if(kernel == "delta"){
    out = run_delta_dyn(dp[1], dp[2], covariates, isRL);
  }
  if(kernel == "delta2"){
    // Note: delta2 implementation would go here
  }
  // Map types
  if(base == "lin"){
    out = param + dp[0]*out;
  }
  if(base == "exp_lin"){
    out = exp(param) + dp[0]*out;
  }
  if(base == "add"){
    out = param + out;
  }
  if(do_lR){
    for(int i = 0; i < max(lR); i ++){
      out[lR == (i + 1)] = out;
    }
  }
  return(out);
}

NumericMatrix map_dyn(List dynamic, DataFrame data, NumericVector p, CharacterVector curr_names,
                      LogicalVector isin) {
  NumericVector dp;
  List cur_dynamic;
  CharacterVector dyn_names = dynamic.names();
  NumericVector p_curr = p[curr_names];
  NumericMatrix out(data.nrow(), curr_names.length());
  NumericVector input(data.nrow());
  for(int q = 0; q < isin.length(); q ++){
    if(isin[q] == TRUE){
      String curr_name = curr_names[q];
      cur_dynamic = dynamic[curr_name];
      String kernel = cur_dynamic["kernel"];
      String base = cur_dynamic["base"];
      CharacterVector dpnames = cur_dynamic["dpnames"];
      CharacterVector covnames = cur_dynamic["covnames"];
      bool do_lR = cur_dynamic["lR1"];
      dp = p[dpnames];
      out(_, q) = dynfuns_c(p_curr[q],dp, kernel, base, data, covnames, do_lR);
    } else{
      input.fill(p_curr[q]);
      out(_, q) = input;
    }
  }
  return(out);
}


NumericVector vector_pow(NumericVector x1, NumericVector x2){
  NumericVector out(x1.length());
  for(unsigned int i = 0; i < out.length(); i ++){
    out[i] = pow(x1[i], x2[i]);
  }
  return(out);
}


NumericVector run_delta_i_map(NumericVector q0, NumericVector alpha, NumericVector target, NumericVector winner){
  NumericVector q(target.length());
  q[0] = q0[0];
  double pe;
  for(int t = 1; t < q.length(); t ++){
    if(NumericVector::is_na(target[t-1]) || winner[t-1] != 1){
      q[t] = q[t-1];
    } else{
      pe = target[t-1] - q[t-1];
      q[t] = q[t-1] + alpha[t-1]*pe;
    }
  }
  return(q);
}


NumericVector run_delta_adapt(NumericVector q0, NumericVector p, NumericMatrix target, bool isRL){
  NumericVector winner(target.nrow());
  NumericVector s_col;
  LogicalVector tmp = contains(colnames(target), "winner");
  if(isRL){
    if(tmp[0] == TRUE){
      winner = target(_, 0);
      s_col = target(_, 1);
    } else{
      winner = target(_, 1);
      s_col = target(_, 0);
    }
    target = target(_, Range(2, target.ncol() - 1));
  } else{
    winner.fill(1);
  }

  NumericMatrix out(target.nrow(), target.ncol());
  for(int r = 0; r < target.ncol(); r ++){
    out(_, r) = run_delta_i_map(q0, p, target(_, r), winner);
  }
  NumericVector out_v(target.nrow());
  if(isRL){
    for(int v = 0; v < target.nrow(); v ++){
      out_v[v] = out(v, s_col[v]-1);
    }
  } else{
    out_v = as<NumericVector>(out);
  }
  return(out_v);
}


// Here param is the param parameter, and dp are the additional parameters. They are all vectors.
// Because here they have already been mapped back to the design (so one parameter per observation * number of accumulators)
// kernel and base inform us what transformation to be done.
// Data and covnames, tell us what columns in the data are informing our transformations
NumericVector adaptfuns_c(NumericVector param, NumericMatrix dp, String kernel, String base, DataFrame data, CharacterVector covnames, bool do_lR) {
  NumericVector out(data.nrow());
  bool isRL = is_true(any(contains(covnames, "winner")));
  NumericVector lR = data["lR"];
  NumericMatrix covariates(data.nrow(), covnames.length());
  CharacterVector data_names = data.names();
  LogicalVector cnames_idx = contains_multiple(data_names, covnames);
  int k = 0;
  CharacterVector col_names(sum(cnames_idx));
  for(int j = 0; j < data.ncol(); j ++){
    if(cnames_idx[j] == TRUE){
      covariates(_, k) = as<NumericVector>(data[j]);
      col_names[k] = data_names[j];
      k++;
    }
  }

  if(do_lR){
    covariates = submat_rcpp(covariates, lR == 1);
    dp = submat_rcpp(dp, lR == 1);
  } else{
    covariates = submat_rcpp(covariates, lR < 9999); // I'm doing something wrong here but this fixes it????
  }
  colnames(covariates) = col_names;
  NumericVector res(covariates.nrow());

  // dyntypes
  if(kernel == "lin_decr"){
    res = -1*as<NumericVector>(covariates);
  }
  if(kernel == "lin_incr"){
    res = as<NumericVector>(covariates);
  }
  if(kernel == "exp_decr"){
    res = exp(dp(_, 1)*as<NumericVector>(covariates));
  }
  if(kernel == "exp_incr"){
    res = 1-exp(-dp(_, 1)*as<NumericVector>(covariates));
  }
  if(kernel == "pow_decr"){
    res = vector_pow((1+as<NumericVector>(covariates)), -dp(_, 1));
  }
  if(kernel == "pow_incr"){
    res = 1-vector_pow((1+as<NumericVector>(covariates)), -dp(_, 1));
  }
  if(kernel == "poly2"){
    res = dp(_, 0)*as<NumericVector>(covariates) + dp(_, 1)*pow(as<NumericVector>(covariates), 2);
  }
  if(kernel == "poly3"){
    res = dp(_, 0)*as<NumericVector>(covariates) + dp(_, 1)*pow(as<NumericVector>(covariates), 2) + dp(_, 2)*pow(as<NumericVector>(covariates), 3);
  }
  if(kernel == "poly4"){
    res = dp(_, 0)*as<NumericVector>(covariates) + dp(_, 1)*pow(as<NumericVector>(covariates), 2) + dp(_, 2)*pow(as<NumericVector>(covariates), 3) + dp(_, 3)*pow(as<NumericVector>(covariates), 4);
  }
  if(kernel == "delta"){
    res = run_delta_adapt(dp(_, 1), dp(_, 2), covariates, isRL);
  }
  if(kernel == "delta2"){
    // Note: delta2 implementation would go here
  }
  // Map types
  if(base == "lin"){
    res = param + dp(_, 0)*res;
  }
  if(base == "exp_lin"){
    res = exp(param) + exp(dp(_, 0))*res;
  }
  if(base == "add"){
    res = param + res;
  }
  if(do_lR){
    for(int i = 0; i < max(lR); i ++){
      out[lR == (i + 1)] = res;
    }
  } else{
    out = res;
  }

  return(out);
}



NumericMatrix map_adaptive(List adaptive, NumericMatrix p, DataFrame data) {
  List cur_dynamic;
  CharacterVector adapt_names = adaptive.names();
  CharacterVector curr_names = colnames(p);
  LogicalVector isin = contains_multiple(curr_names, adapt_names);
  for(int q = 0; q < isin.length(); q ++){
    if(isin[q] == TRUE){
      String curr_name = curr_names[q];
      cur_dynamic = adaptive[curr_name];
      String kernel = cur_dynamic["kernel"];
      String base = cur_dynamic["base"];
      CharacterVector covnames = cur_dynamic["covnames"];
      CharacterVector aptypes = cur_dynamic["aptypes"];
      LogicalVector is_adapt = contains_multiple(curr_names, aptypes);
      NumericMatrix dp(p.nrow(), sum(is_adapt));
      int z = 0;
      for(int k = 0; k < isin.length();k++){
        if(is_adapt[k] == TRUE){
          dp(_, z) = p(_,k);
          z++;
        }
      }
      bool do_lR;
      if(is_true(any(contains(covnames, "winner")))){
        do_lR = FALSE;
      } else{
        do_lR = cur_dynamic["lR1"];
      }
      p(_, q) = adaptfuns_c(p(_, q),dp, kernel, base, data, covnames, do_lR);
    }
  }
  return(p);
}


#endif
