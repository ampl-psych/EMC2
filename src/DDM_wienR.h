#ifndef DDM_wien
#define DDM_wien

#include <Rcpp.h>
using namespace Rcpp;
#include "utility_functions.h"
#include "pdf_fncs.h"
#include "fncs_seven.h"
#include "tools.h"


NumericMatrix Ntransform_DDM(NumericMatrix x) {
  NumericMatrix out(clone(x));
  CharacterVector is_log = {"a","sv","t0","st0","s"};
  CharacterVector is_probit = {"Z","SZ","DP"};
  LogicalVector col_idx_log = contains_multiple(colnames(x), is_log);
  LogicalVector col_idx_probit = contains_multiple(colnames(x), is_probit);
  for(int i = 0; i < x.ncol(); i ++){
    if(col_idx_log[i] == TRUE){
      out (_, i) = exp(out(_, i));
    };
    if(col_idx_probit[i] == TRUE){
      out (_, i) = pnorm_multiple(out(_, i));
    };
  };
  return(out);
}

NumericVector transform_DDM(NumericVector x){
  return(x);
}

NumericVector d_DDM_Wien(NumericVector rts, IntegerVector Rs, NumericMatrix pars){
  int Epsflag = 1;
  double eps = 1e-3;
  int K = 0;
  int Neval = 6000;
  int choice = 0; //the type of integration method to choose.
  //0 = "v", 1 = "a", 2= "sv", 3 = "t0", 4 = "st0", 5 = "s", 6 = "Z", 7 = "SZ",
  int N = rts.length();
  NumericVector out(N);
  for(int i = 0; i < N; i++){
    double pm = (Rs[i]==1) ? -1 : 1;
    if(pars(i,7) == 0 && pars(i, 4) == 0){
      double new_rt = rts[i] - pars(i,3);
      if(new_rt > 0){
        out[i] = dwiener(new_rt*pm, pars(i, 1), pars(i, 0), pars(i, 6), pars(i, 2), eps, K, Epsflag);
      } else{
        out[i] = 	R_NegInf;
      }
    } else{
      double Rval;
      double Rerr;
      double sz = (pars(i,6) < (1 - pars(i,6))) ? 2*pars(i,7)*pars(i,6) : 2*pars(i,7)*(1-pars(i,6));
      ddiff(choice, rts[i], pm, pars(i, 1), pars(i, 0), pars(i, 3), pars(i, 6), sz, pars(i, 2), pars(i,4), eps, K, Epsflag, Neval, &Rval, &Rerr);
      out[i] = log(Rval);
    }
  }
  return(out);
}


#endif

