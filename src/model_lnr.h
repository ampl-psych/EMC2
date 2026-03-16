#ifndef lnr_h
#define lnr_h

#include <Rcpp.h>
#include "utility_functions.h"
#include "ParamTable.h"

using namespace Rcpp;

struct LNRSpec {
  int col_m;
  int col_s;
  int col_t0;
};

LNRSpec make_lnr_spec(const ParamTable& pt) {
  LNRSpec s;
  s.col_m   = pt.base_index_for("m");
  s.col_s   = pt.base_index_for("s");
  s.col_t0  = pt.base_index_for("t0");
  return s;
}

NumericVector plnr_c_pt(NumericVector rts,
                        const ParamTable& pt,
                        const LNRSpec& spec,
                        LogicalVector idx,
                        double min_ll,
                        LogicalVector is_ok)
{
  // 0 = m, 1 = s, 2 = t0
  const int N       = rts.size();
  const int out_len = sum(idx);

  NumericVector out(out_len);
  double* out_ptr = out.begin();

  const double* rt  = rts.begin();
  const double* m   = &pt.base(0, spec.col_m);
  const double* s   = &pt.base(0, spec.col_s);
  const double* t0  = &pt.base(0, spec.col_t0);

  int* idx_ptr   = LOGICAL(idx);
  int* ok_ptr    = LOGICAL(is_ok);

  int k = 0;
  for (int i = 0; i < N; ++i) {
    if (!idx_ptr[i]) continue;

    if (std::isnan(m[i])) {
      // same behaviour: missing accumulator gets 0, not min_ll
      out_ptr[k++] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff > 0.0 && ok_ptr[i]) {
      out_ptr[k++] = R::plnorm(t_eff, m[i], s[i], /*lower_tail*/ true, /*log_p*/ false);
    } else {
      out_ptr[k++] = min_ll;
    }
  }

  return out;
}

NumericVector dlnr_c_pt(NumericVector rts,
                        const ParamTable& pt,
                        const LNRSpec& spec,
                        LogicalVector idx,
                        double min_ll,
                        LogicalVector is_ok)
{
  const int N       = rts.size();
  const int out_len = sum(idx);

  NumericVector out(out_len);
  double* out_ptr = out.begin();

  const double* rt  = rts.begin();
  const double* m   = &pt.base(0, spec.col_m);
  const double* s   = &pt.base(0, spec.col_s);
  const double* t0  = &pt.base(0, spec.col_t0);

  int* idx_ptr   = LOGICAL(idx);
  int* ok_ptr    = LOGICAL(is_ok);

  int k = 0;
  for (int i = 0; i < N; ++i) {
    if (!idx_ptr[i]) continue;

    if (std::isnan(m[i])) {
      out_ptr[k++] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff > 0.0 && ok_ptr[i]) {
      out_ptr[k++] = R::dlnorm(t_eff, m[i], s[i], /*log_p*/ false);
    } else {
      out_ptr[k++] = min_ll;
    }
  }

  return out;
}

NumericVector plnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  // 0 = m, 1 = s, 2 = t0
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rts[i] - pars(i,2) > 0) && (is_ok[i] == TRUE)){
        out[k] = R::plnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), TRUE, FALSE);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

NumericVector dlnr_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0; // This is a bit tricky, but helps with assigning missing values a zero (instead of min_ll value)
        // which is important for RACE
      } else if((rts[i] - pars(i,2) > 0) && (is_ok[i] == TRUE)){
        out[k] = R::dlnorm(rts[i] - pars(i,2), pars(i, 0), pars(i, 1), FALSE);
      } else{
        out[k] = min_ll;
      }
      k++;
    }

  }

  return(out);
}


#endif
