#ifndef lba_h
#define lba_h

#include <Rcpp.h>
#include "utility_functions.h"
#include "ParamTable.h"

using namespace Rcpp;

struct LBASpec {
  int col_v;
  int col_sv;
  int col_B;
  int col_A;
  int col_t0;
};

LBASpec make_lba_spec(const ParamTable& pt) {
  LBASpec s;
  s.col_v   = pt.base_index_for("v");
  s.col_sv  = pt.base_index_for("sv");
  s.col_B   = pt.base_index_for("B");
  s.col_A   = pt.base_index_for("A");
  s.col_t0  = pt.base_index_for("t0");
  return s;
}


double pnormP(double q, double mean = 0.0, double sd = 1.0,
              bool lower = true, bool log = false){
  return R::pnorm(q, mean, sd, lower, log);
}

double dnormP(double x, double mean = 0.0, double sd = 1.0,
              bool log = false){
  return R::dnorm(x, mean, sd, log);
}

double plba_norm(double t, double A, double b, double v, double sv,
                 bool posdrift = true){
  double denom = 1.;
  if (posdrift) {
    denom = pnormP(v / sv, 0., 1., true, false);
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
    cdf = (1. + (zs * (dnormP(cz_max, 0., 1., false) - dnormP(cz, 0., 1., false))
                   + xx * pnormP(cz_max, 0., 1., true, false) - cmz * pnormP(cz, 0., 1., true, false))/A) / denom;
  } else {
    cdf = pnormP(b / t, v, sv, false, false) / denom;
  }

  if (cdf < 0.) {
    return 0.;
  } else if (cdf > 1.){
    return 1.;
  }
  return cdf;
}

double dlba_norm(double t, double A,double b, double v, double sv,
                 bool posdrift = true){
  double denom = 1.;
  if (posdrift) {
    denom = pnormP(v / sv, 0., 1., true, false);
    if (denom < 1e-10)
      denom = 1e-10;
  }

  double pdf;

  if (A > 1e-10){
    double zs = t * sv;
    double cmz = b - t * v;;
    double cz = cmz / zs;
    double cz_max = (cmz - A) / zs;
    pdf = (v * (pnormP(cz, 0., 1., true, false) - pnormP(cz_max, 0., 1., true, false)) +
      sv * (dnormP(cz_max, 0., 1., false) - dnormP(cz, 0., 1., false))) / (A * denom);
  } else {
    pdf = dnormP(b / t, v, sv, false) * b / (t * t * denom);
  }

  if (pdf < 0.) {
    return 0.;
  }
  return pdf;
}

// LBA
NumericVector dlba_c_pt(NumericVector rts,
                        const ParamTable& pt,
                        const LBASpec& spec,
                        LogicalVector idx,
                        double min_ll,
                        LogicalVector is_ok)
{
  const int N = rts.size();
  const int out_len = sum(idx);
  NumericVector out(out_len);
  double* out_ptr = out.begin();

  const double* rt  = rts.begin();
  const double* v   = &pt.base(0, spec.col_v);
  const double* sv  = &pt.base(0, spec.col_sv);
  const double* B   = &pt.base(0, spec.col_B);
  const double* A   = &pt.base(0, spec.col_A);
  const double* t0  = &pt.base(0, spec.col_t0);

  int* idx_ptr   = LOGICAL(idx);
  int* ok_ptr    = LOGICAL(is_ok);

  int k = 0;
  for (int i = 0; i < N; ++i) {
    if (!idx_ptr[i]) continue;

    if (std::isnan(v[i])) {       // matches NumericVector::is_na(pars(i,0))
      out_ptr[k++] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff > 0.0 && ok_ptr[i]) {
      // b parameter in your code is B + A
      const double A_i = A[i];
      const double B_i = B[i] + A_i;
      out_ptr[k++] = dlba_norm(t_eff, A_i, B_i, v[i], sv[i], true);
    } else {
      out_ptr[k++] = min_ll;
    }
  }
  return out;
}

NumericVector plba_c_pt(NumericVector rts,
                        const ParamTable& pt,
                        const LBASpec& spec,
                        LogicalVector idx,
                        double min_ll,
                        LogicalVector is_ok)
{
  const int N = rts.size();
  const int out_len = sum(idx);
  NumericVector out(out_len);
  double* out_ptr = out.begin();

  const double* rt  = rts.begin();
  const double* v   = &pt.base(0, spec.col_v);
  const double* sv  = &pt.base(0, spec.col_sv);
  const double* B   = &pt.base(0, spec.col_B);
  const double* A   = &pt.base(0, spec.col_A);
  const double* t0  = &pt.base(0, spec.col_t0);

  int* idx_ptr   = LOGICAL(idx);
  int* ok_ptr    = LOGICAL(is_ok);

  int k = 0;
  for (int i = 0; i < N; ++i) {
    if (!idx_ptr[i]) continue;

    if (std::isnan(v[i])) {
      out_ptr[k++] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff > 0.0 && ok_ptr[i]) {
      const double A_i = A[i];
      const double B_i = B[i] + A_i;
      out_ptr[k++] = plba_norm(t_eff, A_i, B_i, v[i], sv[i], true);
    } else {
      out_ptr[k++] = min_ll;
    }
  }
  return out;
}

NumericVector dlba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = dlba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}

NumericVector plba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, sv = 1, B = 2, A = 3, t0 = 4
  int n = sum(idx);
  NumericVector out(n);
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
        out[k] = plba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true);
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true)

{
  int n = t.size();
  NumericVector pdf(n);

  for (int i = 0; i < n; i++){
    pdf[i] = dlba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector plba(NumericVector t,
                   NumericVector A, NumericVector b, NumericVector v, NumericVector sv,
                   bool posdrift = true)

{
  int n = t.size();
  NumericVector cdf(n);

  for (int i = 0; i < n; i++){
    cdf[i] = plba_norm(t[i], A[i], b[i], v[i], sv[i], posdrift);
  }
  return cdf;
}


void dlba_pt_fill(const NumericVector& rts,
                  const ParamTable& pt,
                  const LBASpec& spec,
                  const LogicalVector& winner,
                  double* raw)
{
  const int N = rts.size();

  const double* rt  = rts.begin();
  const double* v   = &pt.base(0, spec.col_v);
  const double* sv  = &pt.base(0, spec.col_sv);
  const double* B   = &pt.base(0, spec.col_B);
  const double* A   = &pt.base(0, spec.col_A);
  const double* t0  = &pt.base(0, spec.col_t0);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (!win_ptr[i])
      continue;  // only winners get pdf

    if (std::isnan(v[i])) {
      raw[i] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0) {
      raw[i] = 0.0;
      continue;
    }

    const double A_i = A[i];
    const double B_i = B[i] + A_i;  // b = B + A, as in dlba_c_pt

    double pdf = dlba_norm(t_eff, A_i, B_i, v[i], sv[i], /*posdrift=*/true);
    if (!std::isfinite(pdf) || pdf < 0.0)
      pdf = 0.0;

    raw[i] = pdf;
  }
}

void plba_pt_fill(const NumericVector& rts,
                  const ParamTable& pt,
                  const LBASpec& spec,
                  const LogicalVector& winner,
                  double* raw)
{
  const int N = rts.size();

  const double* rt  = rts.begin();
  const double* v   = &pt.base(0, spec.col_v);
  const double* sv  = &pt.base(0, spec.col_sv);
  const double* B   = &pt.base(0, spec.col_B);
  const double* A   = &pt.base(0, spec.col_A);
  const double* t0  = &pt.base(0, spec.col_t0);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (win_ptr[i])
      continue;  // only losers get cdf

    if (std::isnan(v[i])) {
      raw[i] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0) {
      raw[i] = 0.0;
      continue;
    }

    const double A_i = A[i];
    const double B_i = B[i] + A_i;

    double cdf = plba_norm(t_eff, A_i, B_i, v[i], sv[i], /*posdrift=*/true);
    if (!std::isfinite(cdf) || cdf < 0.0)
      cdf = 0.0;
    else if (cdf > 1.0)
      cdf = 1.0;

    raw[i] = cdf;
  }
}
#endif


