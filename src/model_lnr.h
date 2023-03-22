#ifndef lnr_h
#define lnr_h

#include <Rcpp.h>
#include "utility_functions.h"

using namespace Rcpp;


NumericVector transform_lnr(NumericVector x){
  return(x);
}

NumericMatrix Ntransform_lnr(NumericMatrix x) {
  NumericMatrix out(clone(x));
  LogicalVector col_idx = contains(colnames(x), "m");
  for(int i = 0; i < x.ncol(); i ++){
    if(col_idx[i] == FALSE){
      out (_, i) = exp(out(_, i));
    };
  };
  return(out);
}

NumericVector plnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx){
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(rts[i] - pars(i,2) > 0){
        out[k] = R::plnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), TRUE, FALSE);
      }
      k++;
    }
  }

  return(out);
}

NumericVector dlnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx){
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(rts[i] - pars(i,2) > 0){
        out[k] = R::dlnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), FALSE);
      }
      k++;
    }

  }

  return(out);
}


#endif
