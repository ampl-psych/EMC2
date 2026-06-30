#ifndef DDM_wien
#define DDM_wien

#include <Rcpp.h>
using namespace Rcpp;
#include "utility_functions.h"
#include "pdf_fncs.h"
#include "fncs_seven.h"
#include "tools.h"
#include "RaceSpec.h"


void fill_ddm(const NumericVector& rts,
              const IntegerVector& R,
              const ParamTable& pt,
              const RaceSpec& spec,
              const std::vector<int>& idx_all,
              double* __restrict__ ll_row)
{
    const double* __restrict__ rt_ptr = rts.begin();
    const int*    __restrict__ R_ptr  = R.begin();

    const double* __restrict__ v   = &pt.base(0, spec.col_v);
    const double* __restrict__ a   = &pt.base(0, spec.col_a);
    const double* __restrict__ t0  = &pt.base(0, spec.col_t0);
    const double* __restrict__ z   = &pt.base(0, spec.col_Z);
    const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
    const double* __restrict__ s   = &pt.base(0, spec.col_s);
    const double* __restrict__ sz  = &pt.base(0, spec.col_SZ);
    const double* __restrict__ st0 = &pt.base(0, spec.col_st0);
    const int n = (int)idx_all.size();

    int    Epsflag = 1;
    double eps     = 5e-3;
    int    K       = 0;

    for (int j = 0; j < n; ++j) {
        const int i       = idx_all[j];
        const double pm   = (R_ptr[i] == 1) ? -1.0 : 1.0;
        const double teff = rt_ptr[i] - t0[i];

        if (sz[i] == 0.0 && st0[i] == 0.0) {
          if (teff <= 0.0) { ll_row[i] = R_NegInf; continue; }

            double val = dwiener(
                teff * pm,
                a[i]  / s[i],
                v[i]  / s[i],
                z[i],
                sv[i] / s[i],
                eps, K, Epsflag);
          ll_row[i] = std::isfinite(val) ? val : R_NegInf;
        } else {
            int    Neval = 6000;
            int    choice = 0;
            double Rval, Rerr;

            double sz_i = (z[i] < (1.0 - z[i]))
                ? 2.0 * sz[i] *        z[i]
                : 2.0 * sz[i] * (1.0 - z[i]);

            ddiff(choice,
                  rt_ptr[i], pm,
                  a[i]  / s[i],
                  v[i]  / s[i],
                  t0[i], z[i], sz_i,
                  sv[i] / s[i],
                  st0[i],
                  eps, K, Epsflag, Neval, &Rval, &Rerr);

            ll_row[i] = (std::isfinite(Rval) && Rval > 0.0) ? std::log(Rval) : R_NegInf;
        }
    }
}


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


#endif
