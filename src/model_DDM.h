#ifndef DDM_wien
#define DDM_wien

#include <Rcpp.h>
using namespace Rcpp;
#include "utility_functions.h"
#include "pdf_fncs.h"
#include "cdf_fncs.h"
#include "ddm_functions_inline.h"
#include "fncs_seven.h"
#include "tools.h"

NumericVector d_DDM_Wien(NumericVector rts, IntegerVector Rs, NumericMatrix pars, std::vector<int> is_ok){
  int Epsflag = 1;
  double eps = 5e-3;
  int K = 0;
  int Neval = 6000;
  int choice = 0; //the type of integration method to choose.
  //0 = "v", 1 = "a", 2= "sv", 3 = "t0", 4 = "st0", 5 = "s", 6 = "Z", 7 = "SZ",
  int N = rts.length();
  NumericVector out(N);
  for(int i = 0; i < N; i++){
    if (!is_ok[i]) {
      out[i] = R_NegInf;
    } else{
      // we divide v, a and sv by s to introduce the scaling parameter s
      double pm = (Rs[i]==1) ? -1 : 1;
      // if sz and st0 are zero we can use simple and fast dwiener function
      if(pars(i,7) == 0 && pars(i, 4) == 0){
        double new_rt = rts[i] - pars(i,3);
        if(new_rt > 0){
          out[i] = dwiener(new_rt*pm, pars(i, 1)/pars(i,5), pars(i, 0)/pars(i,5), pars(i, 6), pars(i, 2)/pars(i,5), eps, K, Epsflag);
        } else{
          out[i] = 	R_NegInf;
        }
      } else{ // otherwise use ddiff function with integration
        double Rval;
        double Rerr;
        double sz = (pars(i,6) < (1 - pars(i,6))) ? 2*pars(i,7)*pars(i,6) : 2*pars(i,7)*(1-pars(i,6));
        ddiff(choice, rts[i], pm, pars(i, 1)/pars(i,5), pars(i, 0)/pars(i,5), pars(i, 3), pars(i, 6), sz, pars(i, 2)/pars(i,5), pars(i,4), eps, K, Epsflag, Neval, &Rval, &Rerr);
        out[i] = log(Rval);
      }
    }
  }
  return(out);
}


// Defective DDM CDF: log P(RT < t, response R) per row. Companion to d_DDM_Wien.
// Simple case (sv=SZ=st0=0) uses pwiener_inline; variability uses pdiff.
// pwiener() has no signed-time argument, so the upper boundary is obtained by
// reflecting the diffusion (v -> -v, w -> 1-w) and evaluating the lower CDF.
NumericVector p_DDM_Wien(NumericVector rts, IntegerVector Rs, NumericMatrix pars, LogicalVector is_ok){
  int Epsflag = 1;
  double eps = 5e-3;
  int K = 0;
  int Neval = 6000;
  int choice = 0; //the type of integration method to choose.
  //0 = "v", 1 = "a", 2= "sv", 3 = "t0", 4 = "st0", 5 = "s", 6 = "Z", 7 = "SZ",
  int N = rts.length();
  NumericVector out(N);
  for(int i = 0; i < N; i++){
    if(is_ok[i] == FALSE){
      out[i] = R_NegInf;
    } else{
      // we divide v, a and sv by s to introduce the scaling parameter s
      double new_rt = rts[i] - pars(i,3);
      if(new_rt <= 0){
        out[i] = R_NegInf;
      } else if(pars(i,2) == 0 && pars(i,7) == 0 && pars(i, 4) == 0){
        // if sv, sz and st0 are zero we can use simple and fast pwiener function
        double v = pars(i, 0)/pars(i,5);
        double w = pars(i, 6);
        if (Rs[i] != 1) { // upper boundary
          v = -v;
          w = 1.0 - w;
        }
        out[i] = pwiener_inline(new_rt, pars(i, 1)/pars(i,5), v, w, eps, K, Epsflag);
      } else{ // otherwise use pdiff function with integration
        double Rval;
        double Rerr;
        double pm = (Rs[i]==1) ? -1 : 1;
        double sz = (pars(i,6) < (1 - pars(i,6))) ? 2*pars(i,7)*pars(i,6) : 2*pars(i,7)*(1-pars(i,6));
        pdiff(choice, rts[i], pm, pars(i, 1)/pars(i,5), pars(i, 0)/pars(i,5), pars(i, 3), pars(i, 6), sz, pars(i, 2)/pars(i,5), pars(i,4), eps, K, Epsflag, Neval, &Rval, &Rerr);
        out[i] = (Rval > 0.0 && R_FINITE(Rval)) ? log(Rval) : R_NegInf;
      }
      // log CDF must be finite and <= 0; guard against numerical spillover.
      if (!R_FINITE(out[i])) {
        out[i] = R_NegInf;
      } else if (out[i] > 0.0) {
        out[i] = 0.0;
      }
    }
  }
  return(out);
}


#endif
