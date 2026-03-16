#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "ParamTable.h"

// Compile-time switch
#define USE_FAST_PNORM
#ifdef USE_FAST_PNORM
#define PPNORM_STD(x, lower, logp) fast_pnorm_std((x), (lower), (logp))
#else
// exact R implementation
#define PPNORM_STD(x, lower, logp) R::pnorm((x), 0.0, 1.0, (lower), (logp))
#endif

// // [[Rcpp::export]]
// bool rdm_using_fast_pnorm() {
// #ifdef USE_FAST_PNORM
//   return true;
// #else
//   return false;
// #endif
// }

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
double phi(double x)
{
  static const double RT2PI = std::sqrt(4.0 * std::acos(0.0));
  static const double SPLIT = 7.07106781186547;

  static const double N0 = 220.206867912376;
  static const double N1 = 221.213596169931;
  static const double N2 = 112.079291497871;
  static const double N3 = 33.912866078383;
  static const double N4 = 6.37396220353165;
  static const double N5 = 0.700383064443688;
  static const double N6 = 3.52624965998911e-02;

  static const double M0 = 440.413735824752;
  static const double M1 = 793.826512519948;
  static const double M2 = 637.333633378831;
  static const double M3 = 296.564248779674;
  static const double M4 = 86.7807322029461;
  static const double M5 = 16.064177579207;
  static const double M6 = 1.75566716318264;
  static const double M7 = 8.83883476483184e-02;

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
    if (cdf <= 0.0) return -INFINITY;
    return std::log(cdf);
  }
  return cdf;
}

// RDM
double pigt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.){
  //  return 0.;
  //}
  double mu = k / l;
  double lambda = k * k;

  double z1 = std::sqrt(lambda / t) * (1.0 + t / mu);
  double z2 = std::sqrt(lambda / t) * (1.0 - t / mu);

  double p1 = fast_pnorm_std(z1, /*lower=*/false, /*logp=*/false);
  double p2 = fast_pnorm_std(z2, /*lower=*/false, /*logp=*/false);
  // double p1 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. + t/mu), 0., 1., true, false);
  // double p2 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. - t/mu), 0., 1., true, false);

  return std::exp(std::exp(std::log(2. * lambda) - std::log(mu)) + std::log(p1)) + p2;
}

double digt0(double t, double k = 1., double l = 1.){
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

// double pigt(double t, double k = 1, double l = 1, double a = .1, double threshold = 1e-10){
//   if (t <= 0.){
//     return 0.;
//   }
//   if (a < threshold){
//     return pigt0(t, k, l);
//   }
//
//   double sqt = std::sqrt(t);
//   double lgt = std::log(t);
//   double cdf;
//
//   if (l < threshold){
//     double t5a = 2. * R::pnorm((k + a) / sqt, 0., 1., true, false) - 1;
//     double t5b = 2. * R::pnorm((- k - a) / sqt, 0., 1., true, false) - 1;
//
//     double t6a = - .5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
//     double t6b = - .5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);
//
//     cdf = 1. + std::exp(t6a) - std::exp(t6b) + ((- k + a) * t5a - (k - a) * t5b) / (2. * a);
//   } else {
//     double t1a = std::exp(- .5 * std::pow(k - a - t * l, 2) / t);
//     double t1b = std::exp(- .5 * std::pow(a + k - t * l, 2) / t);
//     double t1 = std::exp(.5* (lgt - M_LN2 - L_PI)) * (t1a - t1b);
//
//     double t2a = std::exp(2. * l * (k - a) + R::pnorm(- (k - a + t * l) / sqt, 0., 1., true, true));
//     double t2b = std::exp(2. * l * (k + a) + R::pnorm(- (k + a + t * l) / sqt, 0., 1., true, true));
//     double t2 = a + (t2b - t2a) / (2. * l);
//
//     double t4a = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
//     double t4b = 2. * R::pnorm((k - a) / sqt - sqt * l, 0., 1., true, false) - 1.;
//     double t4 = .5 * (t * l - a - k + .5 / l) * t4a + .5 * (k - a - t * l - .5 / l) * t4b;
//
//     cdf = .5 * (t4 + t2 + t1) / a;
//   }
//   if (cdf < 0. || std::isnan(cdf)) {
//     return 0.;
//   }
//   return cdf;
// }

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
    return 1.0;  // optional clamp, like your LBA code
  }
  return cdf;
}

// double digt(double t, double k = 1., double l = 1., double a = .1, double threshold= 1e-10){
//   if (t <= 0.){
//     return 0.;
//   }
//   if (a < threshold){
//     return digt0(t, k, l);
//   }
//   double pdf;
//   if (l < threshold){
//     double term = std::exp(- (k - a) * (k - a) / (2. * t)) - std::exp(- (k + a) * (k + a) / (2. * t));
//     pdf = std::exp(-.5 * (M_LN2 + L_PI + std::log(t)) + std::log(term) - M_LN2 - std::log(a));
//   } else {
//     double sqt = std::sqrt(t);
//
//     double t1a = - std::pow(a - k + t * l, 2) / (2. * t);
//     double t1b = - std::pow(a + k - t * l, 2) / (2. * t);
//     double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);
//
//     double t2a = 2. * R::pnorm((- k + a) / sqt + sqt * l, 0., 1., true, false) - 1.;
//     double t2b = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
//     double t2 = std::exp(std::log(.5) + std::log(l)) * (t2a + t2b);
//
//     pdf = std::exp(std::log(t1 + t2) - M_LN2 - std::log(a));
//   }
//   if (pdf < 0. || std::isnan(pdf)) {
//     return 0.;
//   }
//   return pdf;
// }
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

// // [[Rcpp::export]]
// NumericVector bench_pnorm_vs_phi(int n) {
//   double sum_exact = 0.0;
//   double sum_approx = 0.0;
//   double x = -5.0;
//   double step = 10.0 / n;
//
//   // warm-up
//   for (int i = 0; i < 100; ++i) {
//     sum_exact  += R::pnorm(x, 0.0, 1.0, true, false);
//     sum_approx += phi(x);
//     x += step;
//   }
//
//   // timing
//   x = -5.0;
//   clock_t t1 = clock();
//   for (int i = 0; i < n; ++i) {
//     sum_exact += R::pnorm(x, 0.0, 1.0, true, false);
//     x += step;
//   }
//   clock_t t2 = clock();
//   x = -5.0;
//   for (int i = 0; i < n; ++i) {
//     sum_approx += phi(x);
//     x += step;
//   }
//   clock_t t3 = clock();
//
//   NumericVector out(4);
//   out[0] = sum_exact;
//   out[1] = sum_approx;
//   out[2] = (double)(t2 - t1) / CLOCKS_PER_SEC;
//   out[3] = (double)(t3 - t2) / CLOCKS_PER_SEC;
//   return out;
// }

// inline, pt
NumericVector drdm_c_pt(NumericVector rts,
                        const ParamTable& pt,
                        const RDMSpec& spec,
                        LogicalVector idx,
                        double min_ll,
                        LogicalVector is_ok)
{
  const int N       = rts.size();
  const int out_len = sum(idx);

  NumericVector out(out_len);
  double* out_ptr = out.begin();

  const double* rt  = rts.begin();
  const double* v   = &pt.base(0, spec.col_v);
  const double* B   = &pt.base(0, spec.col_B);
  const double* A   = &pt.base(0, spec.col_A);
  const double* t0  = &pt.base(0, spec.col_t0);
  const double* s   = &pt.base(0, spec.col_s);

  int* idx_ptr   = LOGICAL(idx);
  int* ok_ptr    = LOGICAL(is_ok);

  const double threshold      = 1e-10;
  const double inv_sqrt_2pi   = 1.0 / std::sqrt(2.0 * M_PI);

  int k = 0;
  for (int i = 0; i < N; ++i) {
    if (!idx_ptr[i]) continue;

    // unused accumulator?
    if (std::isnan(v[i])) {
      out_ptr[k++] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0 || !ok_ptr[i]) {
      out_ptr[k++] = min_ll;
      continue;
    }

    // Map to inverse Gaussian parameters as before
    const double inv_s = 1.0 / s[i];
    const double Ai    = A[i];
    const double Bi    = B[i];
    const double vi    = v[i];

    const double a = 0.5 * Ai * inv_s;
    const double l = vi   * inv_s;
    const double k_param = Bi * inv_s + a;  // B/s + 0.5*A/s

    double pdf = 0.0;

    // Small-a -> use original digt0 as safeguard if desired
    if (a < threshold) {
      pdf = digt0(t_eff, k_param, l);
    } else if (std::fabs(l) < threshold) {
      // small-l branch from your original digt
      const double inv_two_t = 0.5 / t_eff;
      const double km_a = k_param - a;
      const double kp_a = k_param + a;
      const double term =
        std::exp(-km_a * km_a * inv_two_t) -
        std::exp(-kp_a * kp_a * inv_two_t);

      pdf = std::exp(
        -0.5 * (M_LN2 + L_PI + std::log(t_eff)) +
          std::log(term) - M_LN2 - std::log(a)
      );
    } else {
      // main digt branch, inlined and optimised
      const double sqt     = std::sqrt(t_eff);
      const double inv_sqt = 1.0 / sqt;
      const double inv_t   = 1.0 / t_eff;

      const double temp1 = a - k_param + t_eff * l;
      const double temp2 = a + k_param - t_eff * l;
      const double t1a   = -0.5 * temp1 * temp1 * inv_t;
      const double t1b   = -0.5 * temp2 * temp2 * inv_t;

      const double exp1  = std::exp(t1a);
      const double exp2  = std::exp(t1b);

      const double t1 = inv_sqrt_2pi * (exp1 - exp2) * inv_sqt;

      const double arg1 = (-k_param + a) * inv_sqt + sqt * l;
      const double arg2 = ( k_param + a) * inv_sqt - sqt * l;
      const double t2a  = 2.0 * PPNORM_STD(arg1, /*lower=*/true,  /*logp=*/false) - 1.0;
      const double t2b  = 2.0 * PPNORM_STD(arg2, /*lower=*/true,  /*logp=*/false) - 1.0;
      const double t2   = 0.5 * l * (t2a + t2b);

      const double sum = t1 + t2;
      if (sum <= 0.0 || !std::isfinite(sum)) {
        pdf = 0.0;
      } else {
        pdf = sum / (2.0 * a);
      }
    }

    if (pdf < 0.0 || std::isnan(pdf)) {
      pdf = 0.0;
    }
    out_ptr[k++] = pdf;
  }

  return out;
}

NumericVector drdm_c(NumericVector rts,
                     NumericMatrix pars,
                     LogicalVector idx,
                     double min_ll,
                     LogicalVector is_ok)
{
  const int N = rts.size();

  // out has length sum(idx) as before
  const int out_len = sum(idx);
  NumericVector out(out_len);
  double* out_ptr = out.begin();

  const double* rt   = rts.begin();
  const double* p0   = &pars(0, 0); // v
  const double* p1   = &pars(0, 1); // B
  const double* p2   = &pars(0, 2); // A
  const double* p3   = &pars(0, 3); // t0
  const double* p4   = &pars(0, 4); // s

  int* idx_ptr    = LOGICAL(idx);
  int* is_ok_ptr  = LOGICAL(is_ok);

  int k = 0; // index into out
  for (int i = 0; i < N; ++i) {
    if (!idx_ptr[i]) continue;  // skip non‑selected trials

    // If current accumulator is "unused" for this trial
    if (std::isnan(p0[i])) {
      out_ptr[k++] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - p3[i];
    if (t_eff > 0.0 && is_ok_ptr[i]) {
      const double s     = p4[i];
      const double inv_s = 1.0 / s;

      const double k_param = p1[i] * inv_s + 0.5 * p2[i] * inv_s;
      const double l_param = p0[i] * inv_s;
      const double a_param = 0.5 * p2[i] * inv_s;

      out_ptr[k++] = digt(t_eff, k_param, l_param, a_param);
    } else {
      // For invalid times or !is_ok: use minimal density (min_ll is already on the density scale)
      out_ptr[k++] = min_ll;
    }
  }

  return out;
}

// NumericVector drdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
//   //v = 0, B = 1, A = 2, t0 = 3, s = 4
//   NumericVector out(sum(idx));
//   int k = 0;
//   for(int i = 0; i < rts.length(); i++){
//     if(idx[i] == TRUE){
//       if(NumericVector::is_na(pars(i,0))){
//         out[k] = 0;
//       } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
//         out[k] = digt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4));
//       } else{
//         out[k] = min_ll;
//       }
//       k++;
//     }
//   }
//
//   return(out);
// }

NumericVector prdm_c(NumericVector rts,
                     NumericMatrix pars,
                     LogicalVector idx,
                     double min_ll,
                     LogicalVector is_ok)
{
  const int N = rts.size();

  const int out_len = sum(idx);
  NumericVector out(out_len);
  double* out_ptr = out.begin();

  const double* rt   = rts.begin();
  const double* p0   = &pars(0, 0);
  const double* p1   = &pars(0, 1);
  const double* p2   = &pars(0, 2);
  const double* p3   = &pars(0, 3);
  const double* p4   = &pars(0, 4);

  int* idx_ptr    = LOGICAL(idx);
  int* is_ok_ptr  = LOGICAL(is_ok);

  int k = 0;
  for (int i = 0; i < N; ++i) {
    if (!idx_ptr[i]) continue;

    if (std::isnan(p0[i])) {
      out_ptr[k++] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - p3[i];
    if (t_eff > 0.0 && is_ok_ptr[i]) {
      const double s     = p4[i];
      const double inv_s = 1.0 / s;

      const double k_param = p1[i] * inv_s + 0.5 * p2[i] * inv_s;
      const double l_param = p0[i] * inv_s;
      const double a_param = 0.5 * p2[i] * inv_s;

      out_ptr[k++] = pigt(t_eff, k_param, l_param, a_param);
    } else {
      // For invalid times or !is_ok: use minimal CDF (min_ll here is the prob threshold)
      out_ptr[k++] = min_ll;
    }
  }

  return out;
}


NumericVector prdm_c_pt(NumericVector rts,
                        const ParamTable& pt,
                        const RDMSpec& spec,
                        LogicalVector idx,
                        double min_ll,
                        LogicalVector is_ok)
{
  const int N       = rts.size();
  const int out_len = sum(idx);

  NumericVector out(out_len);
  double* out_ptr = out.begin();

  const double* rt  = rts.begin();
  const double* v   = &pt.base(0, spec.col_v);
  const double* B   = &pt.base(0, spec.col_B);
  const double* A   = &pt.base(0, spec.col_A);
  const double* t0  = &pt.base(0, spec.col_t0);
  const double* s   = &pt.base(0, spec.col_s);

  int* idx_ptr   = LOGICAL(idx);
  int* ok_ptr    = LOGICAL(is_ok);

  const double threshold      = 1e-10;
  const double inv_sqrt_2pi   = 1.0 / std::sqrt(2.0 * M_PI);

  int k = 0;
  for (int i = 0; i < N; ++i) {
    if (!idx_ptr[i]) continue;

    if (std::isnan(v[i])) {
      out_ptr[k++] = 0.0;
      continue;
    }

    const double t_eff = rt[i] - t0[i];
    if (t_eff <= 0.0 || !ok_ptr[i]) {
      out_ptr[k++] = min_ll;
      continue;
    }

    const double inv_s = 1.0 / s[i];
    const double Ai    = A[i];
    const double Bi    = B[i];
    const double vi    = v[i];

    const double a = 0.5 * Ai * inv_s;
    const double l = vi   * inv_s;
    const double k_param = Bi * inv_s + a;

    double cdf = 0.0;

    if (a < threshold) {
      cdf = pigt0(t_eff, k_param, l);
    } else if (std::fabs(l) < threshold) {
      // small-l branch from original pigt
      const double sqt = std::sqrt(t_eff);
      const double lgt = std::log(t_eff);

      double t5a = 2.0 * PPNORM_STD((k + a) / sqt,  /*lower=*/true,  /*logp=*/false) - 1.0;
      double t5b = 2.0 * PPNORM_STD((-k - a) / sqt, /*lower=*/true,  /*logp=*/false) - 1.0;

      double t6a = -0.5 * ((k_param + a) * (k_param + a) / t_eff - M_LN2 - L_PI + lgt)
        - std::log(a);
      double t6b = -0.5 * ((k_param - a) * (k_param - a) / t_eff - M_LN2 - L_PI + lgt)
        - std::log(a);

      cdf = 1.0 + std::exp(t6a) - std::exp(t6b)
        + ((-k_param + a) * t5a - (k_param - a) * t5b) / (2.0 * a);
    } else {
      // main pigt branch, inlined and simplified
      const double sqt     = std::sqrt(t_eff);
      const double inv_sqt = 1.0 / sqt;
      const double inv_t   = 1.0 / t_eff;

      const double tmp1 = k_param - a - t_eff * l;
      const double tmp2 = a + k_param - t_eff * l;
      const double t1a  = std::exp(-0.5 * tmp1 * tmp1 * inv_t);
      const double t1b  = std::exp(-0.5 * tmp2 * tmp2 * inv_t);
      const double t1   = sqt * inv_sqrt_2pi * (t1a - t1b);

      const double argA = -(k_param - a + t_eff * l) * inv_sqt;
      const double argB = -(k_param + a + t_eff * l) * inv_sqt;

      const double t2a = std::exp(2.0 * l * (k_param - a) +
                                  PPNORM_STD(argA, /*lower=*/true, /*logp=*/true));
      const double t2b = std::exp(2.0 * l * (k_param + a) +
                                  PPNORM_STD(argB, /*lower=*/true, /*logp=*/true));
      const double t2 = a + (t2b - t2a) / (2.0 * l);

      const double t4a = 2.0 * R::pnorm((k_param + a) * inv_sqt - sqt * l,
                                        0.0, 1.0, true, false) - 1.0;
      const double t4b = 2.0 * R::pnorm((k_param - a) * inv_sqt - sqt * l,
                                        0.0, 1.0, true, false) - 1.0;
      const double t4 = 0.5 * (t_eff * l - a - k_param + 0.5 / l) * t4a
      + 0.5 * (k_param - a - t_eff * l - 0.5 / l) * t4b;

      cdf = 0.5 * (t4 + t2 + t1) / a;
    }

    if (cdf < 0.0 || std::isnan(cdf)) {
      cdf = 0.0;
    } else if (cdf > 1.0) {
      cdf = 1.0;
    }

    out_ptr[k++] = cdf;
  }

  return out;
}

// NumericVector prdm_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
//   //v = 0, B = 1, A = 2, t0 = 3, s = 4
//   NumericVector out(sum(idx));
//   int k = 0;
//   for(int i = 0; i < rts.length(); i++){
//     if(idx[i] == TRUE){
//       if(NumericVector::is_na(pars(i,0))){
//         out[k] = 0;
//       } else if((rts[i] - pars(i,3) > 0) && (is_ok[i] == TRUE)){
//         out[k] = pigt(rts[i] - pars(i,3), pars(i,1)/pars(i,4) + .5 * pars(i,2)/pars(i,4), pars(i,0)/pars(i,4), .5*pars(i,2)/pars(i,4));
//       } else{
//         out[k] = min_ll;
//       }
//       k++;
//     }
//   }
//
//   return(out);
// }


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

