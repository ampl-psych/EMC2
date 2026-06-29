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
#include "RaceSpec.h"   // RaceSpec / RaceScratch for the ddm_truncate backend

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


// Per-row defective DDM CDF: log P(RT < rt, response R). R: 1 = lower, 2 = upper.
// pars order is v, a, sv, t0, st0, s, Z, SZ. Simple case (sv=SZ=st0=0) uses
// pwiener_inline; variability uses pdiff. pwiener() has no signed-time argument,
// so the upper boundary is obtained by reflecting the diffusion (v -> -v,
// w -> 1-w) and evaluating the lower CDF.
inline double ddm_logcdf_one(double rt, int R, double v, double a, double sv,
                             double t0, double st0, double s, double Z, double SZ) {
  const int Epsflag = 1; const double eps = 5e-3; const int K = 0;
  const int Neval = 6000; const int choice = 0;
  double out;
  const double new_rt = rt - t0;            // divide v, a, sv by s for scaling
  if (new_rt <= 0) {
    out = R_NegInf;
  } else if (sv == 0 && SZ == 0 && st0 == 0) {
    double vv = v / s;
    double w  = Z;
    if (R != 1) { vv = -vv; w = 1.0 - w; }  // upper boundary
    out = pwiener_inline(new_rt, a / s, vv, w, eps, K, Epsflag);
  } else {
    double Rval, Rerr;
    const double pm = (R == 1) ? -1 : 1;
    const double sz = (Z < (1 - Z)) ? 2 * SZ * Z : 2 * SZ * (1 - Z);
    pdiff(choice, rt, pm, a / s, v / s, t0, Z, sz, sv / s, st0, eps, K, Epsflag, Neval, &Rval, &Rerr);
    out = (Rval > 0.0 && R_FINITE(Rval)) ? std::log(Rval) : R_NegInf;
  }
  // log CDF must be finite and <= 0; guard against numerical spillover.
  if (!R_FINITE(out)) out = R_NegInf;
  else if (out > 0.0)  out = 0.0;
  return out;
}

// Defective DDM CDF over a matrix of rows. Companion to d_DDM_Wien.
NumericVector p_DDM_Wien(NumericVector rts, IntegerVector Rs, NumericMatrix pars, LogicalVector is_ok){
  int N = rts.length();
  NumericVector out(N);
  for(int i = 0; i < N; i++){
    if(is_ok[i] == FALSE){ out[i] = R_NegInf; continue; }
    out[i] = ddm_logcdf_one(rts[i], Rs[i], pars(i,0), pars(i,1), pars(i,2),
                            pars(i,3), pars(i,4), pars(i,5), pars(i,6), pars(i,7));
  }
  return(out);
}

// DDM survivor backend for TruncSpec/CensorSpec (fill_truncate slot). For each
// active row writes S(bound) = P(RT > bound) = 1 - F_lo(bound) - F_hi(bound),
// the DDM analogue of the race product-of-survivors. n_acc is 1 for DDM, so the
// per-trial reduction in TruncSpec/CensorSpec is a pass-through.
inline void ddm_truncate(const std::vector<int>& idx,
                         const std::vector<double>& bound,
                         const ParamTable& pt,
                         const RaceSpec& spec,
                         double* __restrict__ out,
                         RaceScratch& scratch)
{
  (void)scratch;
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ a   = &pt.base(0, spec.col_a);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);
  const double* __restrict__ st0 = &pt.base(0, spec.col_st0);
  const double* __restrict__ s   = &pt.base(0, spec.col_s);
  const double* __restrict__ Z   = &pt.base(0, spec.col_Z);
  const double* __restrict__ SZ  = &pt.base(0, spec.col_SZ);

  const int n = (int)idx.size();
  for (int j = 0; j < n; ++j) {
    const int i = idx[j];
    const double lo = ddm_logcdf_one(bound[i], 1, v[i], a[i], sv[i], t0[i], st0[i], s[i], Z[i], SZ[i]);
    const double hi = ddm_logcdf_one(bound[i], 2, v[i], a[i], sv[i], t0[i], st0[i], s[i], Z[i], SZ[i]);
    double S = 1.0 - std::exp(lo) - std::exp(hi);   // survivor P(RT > bound)
    if (!std::isfinite(S) || S < 0.0) S = 0.0;
    else if (S > 1.0) S = 1.0;
    out[i] = S;
  }
}


#endif
