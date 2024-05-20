#ifndef fmri_h
#define fmri_h

#include <Rcpp.h>
#include "utility_functions.h"

using namespace Rcpp;


NumericVector transform_fMRI(NumericVector x){
  return(x);
}

NumericMatrix Ntransform_fMRI(NumericMatrix x) {
  NumericMatrix out(clone(x));
  LogicalVector col_idx = contains(colnames(x), "sd");
  for(int i = 0; i < x.ncol(); i ++){
    if(col_idx[i] == FALSE){
      out (_, i) = exp(out(_, i)) + 0.001;
    };
  };
  return(out);
}

double c_log_likelihood_fMRI(NumericMatrix pars, DataFrame data, NumericMatrix designMatrix, double min_ll){
  int n = data.nrow();
  NumericVector y_hat(n);
  NumericVector out(n);
  int n_regr = designMatrix.ncol();
  for(int i = 0; i < n_regr; i ++){
    y_hat = y_hat + designMatrix(_, i) * pars[i];
  }
  double mean_y_hat = mean(y_hat);
  for(int j = 0; j < n; j ++){
    out[j] = R::dnorm4(data[j], y_hat[j] - mean_y_hat, pars(j, n_regr + 1), true);
  }
  out[out < min_ll] = min_ll;
  return(sum(out));
}

#endif