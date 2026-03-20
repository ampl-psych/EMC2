#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "ParamTable.h"

// Compile-time switch
// #define USE_FAST_PNORM
#ifdef USE_FAST_PNORM
#define PPNORM_STD(x, lower, logp) fast_pnorm_std((x), (lower), (logp))
#else
// exact R implementation
#define PPNORM_STD(x, lower, logp) R::pnorm((x), 0.0, 1.0, (lower), (logp))
#endif

using namespace Rcpp;

const double L_PI = 1.1447298858494001741434;  // std::log(M_PI)

struct RDMSpec {
  int col_v;
  int col_B;
  int col_A;
  int col_t0;
  int col_s;
};

RDMSpec make_rdm_spec(const ParamTable& pt) {
  RDMSpec s;
  s.col_v  = pt.base_index_for("v");
  s.col_B  = pt.base_index_for("B");
  s.col_A  = pt.base_index_for("A");
  s.col_t0 = pt.base_index_for("t0");
  s.col_s  = pt.base_index_for("s");
  return s;
}

// Fast pnorm
// Approximate standard normal CDF Φ(x): taken from
// https://stackoverflow.com/a/23119456 (CC BY-SA 3.0)
constexpr double RT2PI = 2.506628274631000502415765284811;
constexpr double SPLIT = 7.07106781186547;

constexpr double N0 = 220.206867912376;
constexpr double N1 = 221.213596169931;
constexpr double N2 = 112.079291497871;
constexpr double N3 = 33.912866078383;
constexpr double N4 = 6.37396220353165;
constexpr double N5 = 0.700383064443688;
constexpr double N6 = 3.52624965998911e-02;

constexpr double M0 = 440.413735824752;
constexpr double M1 = 793.826512519948;
constexpr double M2 = 637.333633378831;
constexpr double M3 = 296.564248779674;
constexpr double M4 = 86.7807322029461;
constexpr double M5 = 16.064177579207;
constexpr double M6 = 1.75566716318264;
constexpr double M7 = 8.83883476483184e-02;

inline double phi(double x)
{

  const double z = std::fabs(x);
  double c = 0.0;

  if (z <= 37.0) {
    const double e = std::exp(-z * z / 2.0);
    if (z < SPLIT) {
      const double n = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
      const double d = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
      c = e * n / d;
    } else {
      const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
      c = e / (RT2PI * f);
    }
  }
  return x <= 0.0 ? c : 1.0 - c;
}

// Convenience wrapper approximating R::pnorm(q, 0, 1, lower, logp)
inline double fast_pnorm_std(double x, bool lower = true, bool logp = false)
{
  double cdf = phi(x);          // Φ(x)
  if (!lower) {
    cdf = 1.0 - cdf;            // upper tail
  }
  if (logp) {
    // Guard against log(0) (underflow). For extremely small cdf, this will
    // be -inf, which your min_ll handling will clamp anyway.
    if (cdf <= 0.0) return R_NegInf;
    return std::log(cdf);
  }
  return cdf;
}


// RDM
inline double pigt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.){
  //  return 0.;
  //}
  double mu = k / l;
  double lambda = k * k;

  double z1 = std::sqrt(lambda / t) * (1.0 + t / mu);
  double z2 = std::sqrt(lambda / t) * (1.0 - t / mu);

  double p1 = fast_pnorm_std(z1, /*lower=*/false, /*logp=*/false);
  double p2 = fast_pnorm_std(z2, /*lower=*/false, /*logp=*/false);

  return std::exp(std::exp(std::log(2. * lambda) - std::log(mu)) + std::log(p1)) + p2;
}

inline double digt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.) {
  //  return 0.;
  //}
  double lambda = k * k;
  double e;
  if (l == 0.) {
    e = -.5 * lambda / t;
  } else {
    double mu = k / l;
    e = - (lambda / (2. * t)) * ((t * t) / (mu * mu) - 2. * t / mu + 1.);
  }
  return std::exp(e + .5 * std::log(lambda) - .5 * std::log(2. * t * t * t * M_PI));
}

double pigt(double t, double k = 1., double l = 1., double a = .1, double threshold = 1e-10)
{
  if (t <= 0.0) {
    return 0.0;
  }
  if (a < threshold) {
    return pigt0(t, k, l);   // can be dropped if a is bounded away from 0
  }

  const double sqt     = std::sqrt(t);
  const double inv_sqt = 1.0 / sqt;
  const double lgt     = std::log(t);
  double cdf;

  if (std::fabs(l) < threshold) {
    // keep your existing l≈0 branch for now
    double t5a = 2.0 * R::pnorm((k + a) / sqt, 0.0, 1.0, true, false) - 1.0;
    double t5b = 2.0 * R::pnorm((-k - a) / sqt, 0.0, 1.0, true, false) - 1.0;

    double t6a = -0.5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
    double t6b = -0.5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);

    cdf = 1.0 + std::exp(t6a) - std::exp(t6b) + ((-k + a) * t5a - (k - a) * t5b) / (2.0 * a);
  } else {
    const double inv_t = 1.0 / t;

    // t1: simplified as sqrt(t)/(sqrt(2π)) * (exp(...) - exp(...))
    const double tmp1 = k - a - t * l;
    const double tmp2 = a + k - t * l;
    const double t1a  = std::exp(-0.5 * tmp1 * tmp1 * inv_t);
    const double t1b  = std::exp(-0.5 * tmp2 * tmp2 * inv_t);
    const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
    const double t1 = sqt * inv_sqrt_2pi * (t1a - t1b);

    // t2: unchanged algebraically
    const double argA = -(k - a + t * l) * inv_sqt;
    const double argB = -(k + a + t * l) * inv_sqt;
    const double t2a = std::exp(2.0 * l * (k - a) +
                                R::pnorm(argA, 0.0, 1.0, true, true));
    const double t2b = std::exp(2.0 * l * (k + a) +
                                R::pnorm(argB, 0.0, 1.0, true, true));
    const double t2 = a + (t2b - t2a) / (2.0 * l);

    // t4: unchanged
    const double t4a = 2.0 * R::pnorm((k + a) * inv_sqt - sqt * l,
                                      0.0, 1.0, true, false) - 1.0;
    const double t4b = 2.0 * R::pnorm((k - a) * inv_sqt - sqt * l,
                                      0.0, 1.0, true, false) - 1.0;
    const double t4 = 0.5 * (t * l - a - k + 0.5 / l) * t4a
    + 0.5 * (k - a - t * l - 0.5 / l) * t4b;

    cdf = 0.5 * (t4 + t2 + t1) / a;
  }

  if (cdf < 0.0 || std::isnan(cdf)) {
    return 0.0;
  }
  if (cdf > 1.0) {
    return 1.0;
  }
  return cdf;
}


double digt(double t, double k = 1., double l = 1., double a = .1, double threshold = 1e-10)
{
  if (t <= 0.) {
    return 0.0;  // you already ensure t>0 in the caller; this is just extra safety
  }

  // Small-a asymptotic: if your bounds guarantee a >> threshold, you can drop this.
  if (a < threshold) {
    return digt0(t, k, l);
  }

  double pdf;

  if (std::fabs(l) < threshold) {
    // small-l branch; again, can be dropped if |l| is bounded away from 0
    const double inv_two_t = 0.5 / t;
    const double km_a = k - a;
    const double kp_a = k + a;
    const double term =
      std::exp(-km_a * km_a * inv_two_t) -
      std::exp(-kp_a * kp_a * inv_two_t);

    // same expression as before, just with minor micro-optimisations
    pdf = std::exp(
      -0.5 * (M_LN2 + L_PI + std::log(t))  // normalising bits
      + std::log(term)                     // log(term)
      - M_LN2
      - std::log(a)
    );
  } else {
    // main branch
    const double sqt    = std::sqrt(t);
    const double inv_sqt = 1.0 / sqt;
    const double inv_t   = 1.0 / t;

    // t1: simplified and avoiding pow
    const double temp1 = a - k + t * l;
    const double temp2 = a + k - t * l;
    const double t1a   = -0.5 * temp1 * temp1 * inv_t;
    const double t1b   = -0.5 * temp2 * temp2 * inv_t;

    // 1/sqrt(2π)
    const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
    const double t1 = inv_sqrt_2pi * (std::exp(t1a) - std::exp(t1b)) * inv_sqt;

    // t2: no log/exp, just 0.5 * l
    const double arg1 = (-k + a) * inv_sqt + sqt * l;
    const double arg2 = ( k + a) * inv_sqt - sqt * l;
    const double t2a = 2.0 * R::pnorm(arg1, 0.0, 1.0, true, false) - 1.0;
    const double t2b = 2.0 * R::pnorm(arg2, 0.0, 1.0, true, false) - 1.0;
    const double t2  = 0.5 * l * (t2a + t2b);

    const double sum = t1 + t2;

    // pdf = exp(log(t1 + t2) - M_LN2 - log(a)) == (t1 + t2) / (2a)
    if (sum <= 0.0 || !std::isfinite(sum)) {
      pdf = 0.0;
    } else {
      pdf = sum / (2.0 * a);
    }
  }

  if (pdf < 0.0 || std::isnan(pdf)) {
    return 0.0;
  }
  return pdf;
}

NumericVector drdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = digt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}

NumericVector prdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
  //v = 0, B = 1, A = 2, t0 = 3, s = 4
  NumericVector out(sum(idx));
  int k = 0;
  for(int i = 0; i < rts.length(); i++){
    if(idx[i] == TRUE){
      if(NumericVector::is_na(pars(i,0))){
        out[k] = 0;
      } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
        out[k] = pigt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4));
      } else{
        out[k] = min_ll;
      }
      k++;
    }
  }

  return(out);
}


// [[Rcpp::export]]
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt(t[i], B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0[i];
    if (t[i] <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt(t[i], B[i] + .5 * A[i], v[i], .5 * A[i]);
    }
  }
  return cdf;
}

#endif



// RDM model
// Minimal values of A and v/s to avoid numerical instability.
constexpr double A_EPS = 1e-4;   // or 1e-3
constexpr double L_EPS = 1e-4;   // or 1e-3

inline void clamp_a_l(double& a, double& l) {
  if (a < A_EPS) a = A_EPS;
  if (l > -L_EPS && l < L_EPS) l = (l >= 0.0 ? L_EPS : -L_EPS);
}

// Core digt and pigt
inline double digt_core(double t, double k, double l, double a)
{
  if (t <= 0.0) {
    return 0.0;
  }

  const double sqt      = std::sqrt(t);
  const double inv_sqt  = 1.0 / sqt;
  const double inv_t    = 1.0 / t;
  const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);

  // t1 part
  const double temp1 = a - k + t * l;
  const double temp2 = a + k - t * l;
  const double t1a   = -0.5 * temp1 * temp1 * inv_t;
  const double t1b   = -0.5 * temp2 * temp2 * inv_t;
  const double t1    = inv_sqrt_2pi * (std::exp(t1a) - std::exp(t1b)) * inv_sqt;

  // t2 part
  const double arg1 = (-k + a) * inv_sqt + sqt * l;
  const double arg2 = ( k + a) * inv_sqt - sqt * l;

  const double t2a = 2.0 * PPNORM_STD(arg1, /*lower=*/true, /*logp=*/false) - 1.0;
  const double t2b = 2.0 * PPNORM_STD(arg2, /*lower=*/true, /*logp=*/false) - 1.0;
  const double t2  = 0.5 * l * (t2a + t2b);

  const double sum = t1 + t2;

  if (sum <= 0.0 || !std::isfinite(sum)) {
    return 0.0;
  }

  double pdf = sum / (2.0 * a);

  if (!std::isfinite(pdf) || pdf < 0.0) {
    return 0.0;
  }
  return pdf;
}

// pigt core
inline double pigt_core(double t, double k, double l, double a)
{
  if (t <= 0.0) {
    return 0.0;
  }

  const double sqt      = std::sqrt(t);
  const double inv_sqt  = 1.0 / sqt;
  const double inv_t    = 1.0 / t;
  const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);

  // t1 term: sqt / sqrt(2π) * (exp(...) - exp(...))
  const double tmp1 = k - a - t * l;
  const double tmp2 = a + k - t * l;
  const double t1a  = std::exp(-0.5 * tmp1 * tmp1 * inv_t);
  const double t1b  = std::exp(-0.5 * tmp2 * tmp2 * inv_t);
  const double t1   = sqt * inv_sqrt_2pi * (t1a - t1b);

  // t2 term
  const double argA = -(k - a + t * l) * inv_sqt;
  const double argB = -(k + a + t * l) * inv_sqt;

  const double t2a = std::exp(2.0 * l * (k - a) +
                              PPNORM_STD(argA, /*lower=*/true, /*logp=*/true));
  const double t2b = std::exp(2.0 * l * (k + a) +
                              PPNORM_STD(argB, /*lower=*/true, /*logp=*/true));
  const double t2  = a + (t2b - t2a) / (2.0 * l);

  // t4 term
  const double t4a = 2.0 * PPNORM_STD((k + a) * inv_sqt - sqt * l,
                                      /*lower=*/true, /*logp=*/false) - 1.0;
  const double t4b = 2.0 * PPNORM_STD((k - a) * inv_sqt - sqt * l,
                                      /*lower=*/true, /*logp=*/false) - 1.0;
  const double t4  = 0.5 * (t * l - a - k + 0.5 / l) * t4a + 0.5 * (k - a - t * l - 0.5 / l) * t4b;

  double cdf = 0.5 * (t4 + t2 + t1) / a;

  if (!std::isfinite(cdf) || cdf < 0.0) {
    return 0.0;
  }
  if (cdf > 1.0) {
    return 1.0;
  }
  return cdf;
}

// pdf for RDM, winners only
// pdf for RDM, winners only (contiguous + mask)
void drdm_c_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RDMSpec& spec,
                    const LogicalVector& winner,
                    double* raw /* length n_trials */)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* v  = &pt.base(0, spec.col_v);
  const double* B  = &pt.base(0, spec.col_B);
  const double* A  = &pt.base(0, spec.col_A);
  const double* t0 = &pt.base(0, spec.col_t0);
  const double* s  = &pt.base(0, spec.col_s);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (!win_ptr[i])
      continue;  // only winners get pdf; others left unchanged

    double t_eff = rt[i] - t0[i];

    double inv_s = 1.0 / s[i];
    double a     = 0.5 * A[i] * inv_s;
    double l     = v[i]   * inv_s;
    double k     = B[i]   * inv_s + a;

    clamp_a_l(a, l);

    double pdf = digt_core(t_eff, k, l, a);   // handles t<=0 internally

    raw[i] = pdf;  // may be 0, NaN, etc.; interpreted later
  }
}

// cdf for RDM, losers only (contiguous + mask)
void prdm_c_fast(const NumericVector& rts,
                    const ParamTable& pt,
                    const RDMSpec& spec,
                    const LogicalVector& winner,
                    double* raw /* length n_trials */)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* v  = &pt.base(0, spec.col_v);
  const double* B  = &pt.base(0, spec.col_B);
  const double* A  = &pt.base(0, spec.col_A);
  const double* t0 = &pt.base(0, spec.col_t0);
  const double* s  = &pt.base(0, spec.col_s);

  int* win_ptr = LOGICAL(winner);

#pragma omp simd
  for (int i = 0; i < N; ++i) {
    if (win_ptr[i])
      continue;  // only losers get cdf

    double t_eff = rt[i] - t0[i];

    double inv_s = 1.0 / s[i];
    double a     = 0.5 * A[i] * inv_s;
    double l     = v[i]   * inv_s;
    double k     = B[i]   * inv_s + a;

    clamp_a_l(a, l);

    double cdf = pigt_core(t_eff, k, l, a);   // handles t<=0 internally

    raw[i] = cdf;  // may be NaN, 0, 1, etc.; interpreted later
  }
}



inline void prdm_wald_fast_one(const int i,
                                  const double* rt,
                                  const double* v,
                                  const double* B,
                                  const double* t0,
                                  const double* s,
                                  int* win_ptr,
                                  double* raw)
{
  // losers only
  if (win_ptr[i])
    return;

  double t_eff = rt[i] - t0[i];
  double k     = B[i] / s[i];
  double l     = v[i] / s[i];

  double cdf = pigt0(t_eff, k, l);
  raw[i] = cdf;
}

// Wald cdf: losers only via winner mask, 2× unrolled
void prdm_wald_fast(const NumericVector& rts,
                       const ParamTable& pt,
                       const RDMSpec& spec,
                       const LogicalVector& winner,
                       double* raw)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* v  = &pt.base(0, spec.col_v);
  const double* B  = &pt.base(0, spec.col_B);
  const double* t0 = &pt.base(0, spec.col_t0);
  const double* s  = &pt.base(0, spec.col_s);

  int* win_ptr = LOGICAL(winner);

  int i = 0;

  // Main 2× unrolled loop
  for (; i + 1 < N; i += 2) {
    prdm_wald_fast_one(i,   rt, v, B, t0, s, win_ptr, raw);
    prdm_wald_fast_one(i+1, rt, v, B, t0, s, win_ptr, raw);
  }

  // Tail (if N is odd)
  if (i < N) {
    prdm_wald_fast_one(i, rt, v, B, t0, s, win_ptr, raw);
  }
}

inline void drdm_wald_fast_one(const int i,
                                  const double* rt,
                                  const double* v,
                                  const double* B,
                                  const double* t0,
                                  const double* s,
                                  int* win_ptr,
                                  double* raw)
{
  if (!win_ptr[i]) return;  // only winners get pdf

  double t_eff = rt[i] - t0[i];
  double k     = B[i] / s[i];   // A==0 → a=0 → k=B/s
  double l     = v[i] / s[i];

  double pdf = digt0(t_eff, k, l);
  raw[i] = pdf;
}

// Wald pdf: winners only via winner mask, 2× unrolled
void drdm_wald_fast(const NumericVector& rts,
                       const ParamTable& pt,
                       const RDMSpec& spec,
                       const LogicalVector& winner,
                       double* raw)
{
  const int N = rts.size();

  const double* rt = rts.begin();
  const double* v  = &pt.base(0, spec.col_v);
  const double* B  = &pt.base(0, spec.col_B);
  const double* t0 = &pt.base(0, spec.col_t0);
  const double* s  = &pt.base(0, spec.col_s);

  int* win_ptr = LOGICAL(winner);

  int i = 0;

  // Main 2× unrolled loop
  for (; i + 1 < N; i += 2) {
    drdm_wald_fast_one(i,   rt, v, B, t0, s, win_ptr, raw);
    drdm_wald_fast_one(i+1, rt, v, B, t0, s, win_ptr, raw);
  }

  // Tail (if N is odd)
  if (i < N) {
    drdm_wald_fast_one(i, rt, v, B, t0, s, win_ptr, raw);
  }
}
