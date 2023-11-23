#ifndef lba_h
#define lba_h

#include <Rcpp.h>
#include "utility_functions.h"

using namespace Rcpp;

double pnormP(double q, double mean = 0.0, double sd = 1.0,
              bool lower = true, bool log = false, bool robust = false){
  if (robust == true){
    if (q < -7.) {
      return 0.;
    } else if (q > 7.){
      return 1.;
    }
  }
  return R::pnorm(q, mean, sd, lower, log);
}

double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log = false, bool robust = false){
  if (robust == true){
    if (x < -7.) {
      return 0.;
    } else if (x > 7.){
      return 1.;
    }
  }
  return R::dnorm(x, mean, sd, log);
}

double plba_norm(double t, double A, double b, double v, double sv,
                 bool posdrift = true, bool robust = false){
  double denom = 1.;
  if (posdrift) {
    denom = pnormP(v / sv, 0., 1., true, false, robust);
    if (denom < 1e-10)
      denom = 1e-10;
  }

  double cdf;

  if (A > 1e-10){
    double zs = t * sv;
    double cmz = b - t * v;
    double xx = cmz - A;
    double cz = cmz / zs;
    double cz_max = xx / zs;
    cdf = (1. + (zs * (dnormP(cz_max, 0., 1., false, robust) - dnormP(cz, 0., 1., false, robust))
                   + xx * pnormP(cz_max, 0., 1., true, false, robust) - cmz * pnormP(cz, 0., 1., true, false, robust))/A) / denom;
  } else {
    cdf = pnormP(b / t, v, sv, false, false, robust) / denom;
  }

  if (cdf < 0.) {
    return 0.;
  } else if (cdf > 1.){
    return 1.;
  }
  return cdf;
}

double dlba_norm(double t, double A,double b, double v, double sv,
                 bool posdrift = true, bool robust = false){
  double denom = 1.;
  if (posdrift) {
    denom = pnormP(v / sv, 0., 1., true, false, robust);
    if (denom < 1e-10)
      denom = 1e-10;
  }

  double pdf;

  if (A > 1e-10){
    double zs = t * sv;
    double cmz = b - t * v;;
    double cz = cmz / zs;
    double cz_max = (cmz - A) / zs;
    pdf = (v * (pnormP(cz, 0., 1., true, false, robust) - pnormP(cz_max, 0., 1., true, false, robust)) +
      sv * (dnormP(cz_max, 0., 1., false, robust) - dnormP(cz, 0., 1., false, robust))) / (A * denom);
  } else {
    pdf = dnormP(b / t, v, sv, false, robust) * b / (t * t * denom);
  }

  if (pdf < 0.) {
    return 0.;
  }
  return pdf;
}

NumericVector dlba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(!NumericVector::is_na(pars(i,0)) & (rts[i] - pars(i,4) > 0) & (pars(i,2) >= 0) & (pars(i,4) > 0.05) &
         ((pars(i,3) > 1e-6) | (pars(i,3) == 0)) &
         ((pars(i,0) > -100) & (pars(i,0) < 100)) & ((pars(i,1) > 1e-3) | (pars(i,1) == 0)) ){
        out[k] = dlba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true, false);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}

NumericVector plba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(!NumericVector::is_na(pars(i,0)) & (rts[i] > 0) & (pars(i,2) >= 0) & (pars(i,4) > 0.05) &
         ((pars(i,3) > 1e-6) | (pars(i,3) == 0)) &
         ((pars(i,0) > -100) & (pars(i,0) < 100)) & ((pars(i,1) > 1e-3) | (pars(i,1) == 0)) ){
        out[k] = plba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true, false);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}

NumericMatrix Ntransform_lba(NumericMatrix x) {
  NumericMatrix out(clone(x));
  LogicalVector col_idx = contains(colnames(x), "v");
  for(int i = 0; i < x.ncol(); i ++){
    if(col_idx[i] == FALSE){
      out (_, i) = exp(out(_, i));
    };
  };
  return(out);
}

NumericVector transform_lba(NumericVector x){
  return(x);
}

// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true, bool robust = false)

{
  int n = t.size();
  NumericVector pdf(n);

  for (int i = 0; i < n; i++){
    pdf[i] = dlba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift, robust);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector plba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true, bool robust = false)

{
  int n = t.size();
  NumericVector cdf(n);

  for (int i = 0; i < n; i++){
    cdf[i] = plba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift, robust);
  }
  return cdf;
}

#endif


