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
  NumericVector m = pars(_, 0);
  NumericVector s = pars(_, 1);
  NumericVector t0 = pars(_, 2);
  m = m[idx];
  s = s[idx];
  t0 = t0[idx];
  rts = rts[idx];
  rts = rts - t0;
  NumericVector out(rts.length());
  for(int i = 0; i < rts.length(); i++){
    if(rts[i] > 0){
      out[i] = R::plnorm(rts[i], m[i], s[i], TRUE, FALSE);
    }
  }

  return(out);
}

NumericVector dlnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx){
  NumericVector m = pars(_, 0);
  NumericVector s = pars(_, 1);
  NumericVector t0 = pars(_, 2);
  m = m[idx];
  s = s[idx];
  t0 = t0[idx];
  rts = rts[idx];
  rts = rts - t0;
  NumericVector out(rts.length());
  for(int i = 0; i < rts.length(); i++){
    if(rts[i] > 0){
      out[i] = R::dlnorm(rts[i], m[i], s[i], FALSE);
    }
  }

  return(out);
}



#endif
