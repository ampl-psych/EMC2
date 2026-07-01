#include <Rcpp.h>
using namespace Rcpp;
#include "model_DDM.h"

// =============================================================================
// Constants
// =============================================================================
static constexpr int    DDM_NEVAL   = 6000;
static constexpr double DDM_EPS     = 5e-3;
static constexpr double DDM_LOG_EPS = -5.298317366548036311258; //std::log(DDM_EPS);

// lnnorm helpers
static constexpr double LNNORM_MAX_X = 38.0;
static constexpr double LNNORM_MIN_X = -1.00e9;
static constexpr double SQR2PI       = 2.506628274631000502415765e0;

// Mathematical constants (only defined if not already provided by <cmath> / Rmath.h)
#ifndef M_LNPI
static constexpr double M_LNPI       = 1.14472988584940017414342735135;
#endif
#ifndef M_PI
static constexpr double M_PI         = 3.14159265358979323846264338328;
#endif
#ifndef M_LN2
static constexpr double M_LN2        = 0.693147180559945309417232121458;
#endif
#ifndef M_LN_SQRT_PI
static constexpr double M_LN_SQRT_PI = 0.572364942924700087071713675677;
#endif


// general helpers
static double logsum(double xa, double xb) {
  double temp;
  if (xa <= R_NegInf) return xb;
  if (xb <= R_NegInf) return xa;
  if (xa > xb) temp = xa + log1p(exp(xb - xa));
  else temp = xb + log1p(exp(xa - xb));
  return temp;
}

static double logdiff(double xa, double xb) {
  double result;
  if (xb <= R_NegInf) return(xa);
  if (xa <= R_NegInf) return(xb);
  if (xa > xb) result = (xa + log1p(-exp(xb - xa)));
  else
    if (xb > xa) result = (xb + log1p(-exp(xa - xb)));
  else result = R_NegInf;
  return result;
}

static double lognormal(double x) {
  return	-0.5*x*x - M_LN_SQRT_PI - 0.5*M_LN2;
}

/*
 The original function lnnorm was written by
 Jean Marie Linhart
 StataCorp LP
 jlinhart@stata.com
 January 4, 2008
 and later modified
 */
static double lnnorm(double z)
{
  int lower ;

  double z2, y, s, p1, q1, p2, q2, t, a1, a2 ;
  double n, m ;

  if (z==0.0e0) return(std::log(0.50e0)) ;

  if (z > LNNORM_MAX_X) return(0.0e0);
  if (z <= LNNORM_MIN_X) return(-0.5e0*z*z);

  if (z<0.0e0) {
    z= -z ;
    lower=1 ;
  }
  else    lower=0 ;
  //upper = !lower ;

  z2 = z*z ;

  y = exp(-0.5*z2) / SQR2PI ;
  n = y/z ;

  if (!( ( z>2.00e0))) {
    z *= y ;
    s=z ;
    t=0.0e0 ;

    for (n=3.0e0;s!=t;n+=2.0e0) {
      t=s ;
      z *= ((z2)/ (n)) ;
      s+=z ;
    }
    if (lower) return(std::log(0.50e0-s)) ;
    return(std::log(0.50e0+s)) ;
  }

  a1=2.0e0 ;
  a2=0.0e0 ;
  n=z2+3.0e0 ;
  p1=1.0e0 ;
  q1=z ;
  p2=(n-1.0e0) ;
  q2=n*z ;
  m = p1/q1 ;
  t = p2/q2 ;
  s = m ;

  for (n+=4.0e0; m!=t && s!=t; n+=4.0e0) {
    a1 -= 8.0 ;
    a2 += a1 ;
    s = a2*p1 + n*p2 ;
    p1=p2 ;
    p2=s ;
    s = a2*q1 + n*q2 ;
    q1=q2 ;
    q2=s ;
    s=m ;
    m=t ;
    if (q2>1.0e30) {
      p1 /= 1.0e30 ;
      p2 /= 1.0e30 ;
      q1 /= 1.0e30 ;
      q2 /= 1.0e30 ;
    }
    t = p2/q2;
  }
  t = lower ? std::log(t) - 0.5*z2 - std::log(SQR2PI) : log1p(-y*t);
  return(t) ;
}

static double logMill(double x) {
  double m;
  if (x > 1.0e5) return -std::log(x);
  m = lnnorm(-x) - lognormal(x);
  return m;
}


// =============================================================================
// PDF series helpers
// =============================================================================

static double ks(double t, double w, double eps) {
  const double K1    = (std::sqrt(2.0 * t) + w) / 2.0;
  const double u_eps = std::fmin(-1.0, M_LN2 + M_LNPI + 2.0 * std::log(t) + 2.0 * eps);
  const double arg   = -t * (u_eps - std::sqrt(-2.0 * u_eps - 2.0));
  const double K2    = (arg > 0.0) ? 0.5 * (std::sqrt(arg) - w) : K1;
  return std::ceil(std::fmax(K1, K2));
}

static double kl(double q, double w, double err) {
  const double K1   = 1.0 / (M_PI * std::sqrt(q));
  const double temp = -2.0 * (std::log(M_PI * q) + err);
  const double K2   = (temp >= 0.0) ? std::sqrt(temp / (M_PI * M_PI * q)) : 0.0;
  return std::ceil(std::fmax(K1, K2));
}

// Short-t series: linear-domain accumulation, log returned
static double logfs_linear(double t, double w, int K) {
  if (w == 0.0) return R_NegInf;
  if (K < 0)    K = 0;

  const double inv_twot = 1.0 / (2.0 * t);
  double sum_plus  = w * std::exp(-w * w * inv_twot);
  double sum_minus = 0.0;

#pragma omp simd reduction(+:sum_plus,sum_minus)
  for (int k = 1; k <= K; ++k) {
    const double t1 = w + 2.0 * k;
    const double t2 = w - 2.0 * k;
    sum_plus  += t1    * std::exp(-t1 * t1 * inv_twot);
    sum_minus += (-t2) * std::exp(-t2 * t2 * inv_twot);
  }

  const double S = sum_plus - sum_minus;
  if (!(S > 0.0)) return R_NegInf;
  return -0.5 * M_LN2 - M_LN_SQRT_PI - 1.5 * std::log(t) + std::log(S);
}

// Large-t series: linear-domain accumulation, log returned
static double logfl_linear(double q, double w, int K) {
  if (w == 0.0 || K <= 0) return R_NegInf;

  const double halfq = q / 2.0;
  double sum_plus = 0.0, sum_minus = 0.0;

#pragma omp simd reduction(+:sum_plus,sum_minus)
  for (int k = 1; k <= K; ++k) {
    const double kpi  = k * M_PI;
    const double base = k * std::sin(kpi * w) * std::exp(-kpi * kpi * halfq);
    sum_plus  += (base > 0.0) ?  base : 0.0;
    sum_minus += (base < 0.0) ? -base : 0.0;
  }

  const double S = sum_plus - sum_minus;
  if (!(S > 0.0)) return R_NegInf;
  return std::log(S) + M_LNPI;
}

// log-density of the Wiener first-passage time
static double dwiener(double q, double a, double vn, double wn, double sv)
{
  if (q == 0.0) return R_NegInf;

  double v, w;
  if (q >= 0.0) { w = 1.0 - wn; v = -vn; }
  else          { q = std::fabs(q); w = wn; v = vn; }

  const double q_asq   = q / (a * a);
  const double eta_sqr = sv * sv;
  const double denom   = 1.0 + eta_sqr * q;
  const double lg1     = (eta_sqr * (a * w) * (a * w) - 2.0 * a * v * w - v * v * q)
    / (2.0 * denom)
    - 2.0 * std::log(a) - 0.5 * std::log(denom);

    const double es  = DDM_LOG_EPS - lg1;
    const double kss = ks(q_asq, w, es);
    const double kll = kl(q_asq, w, es);

    if (2.0 * kss <= kll)
      return lg1 + logfs_linear(q_asq, w, static_cast<int>(kss));
    else
      return lg1 + logfl_linear(q_asq, w, static_cast<int>(kll));
}

// =============================================================================
// CDF series helpers
// =============================================================================

// Asymptotic log-probability of hitting lower boundary (t -> inf)
static double logP(double a, double v, double w) {
  constexpr double em1 = 1.0 - 1.0e-6;
  if (std::fabs(v) == 0.0) return std::log1p(-w);

  const double e = -2.0 * v * a * (1.0 - w);
  double tt;
  if (e < 0.0) {
    tt = std::exp(e);
    if (tt >= em1) return std::log1p(-w);
    tt = std::log1p(-tt) - logdiff(2.0 * v * a * w, e);
  } else {
    tt = std::exp(-e);
    if (tt >= em1) return std::log1p(-w);
    tt = std::log1p(-tt) - std::log1p(-std::exp(2.0 * v * a));
  }
  return tt;
}

static double Ks(double t, double v, double a, double w, double eps) {
  const double K1  = 0.5 * (std::fabs(v) / a * t - w);
  const double arg = std::fmax(0.0, std::fmin(1.0,
                                              std::exp(v * a * w + v * v * t / 2.0 + eps) / 2.0));
  const double K2  = (arg == 0.0) ?  INFINITY
  : (arg == 1.0) ? -INFINITY
  : -std::sqrt(t) / (2.0 * a) * QNORM_STD(arg); //gsl_cdf_ugaussian_Pinv(arg);
  return std::ceil(std::fmax(K1, K1 + K2));
}

static double Kl(double t, double v, double a, double w, double err) {
  const double api   = a / M_PI;
  const double vsq   = v * v;
  const double L1    = std::sqrt(1.0 / t) * api;
  const double inner = err + std::log(M_PI * t / 2.0 * (vsq + (M_PI / a) * (M_PI / a)))
    + v * a * w + vsq * t / 2.0;
  const double L2    = std::sqrt(std::fmax(1.0, -2.0 / t * api * api * inner));
  return std::ceil(std::fmax(L1, L2));
}

// Short-t CDF series: linear-domain accumulation, log returned
static double logFs(double t, double v, double a, double w, int K)
{
  double fplus = -INFINITY, fminus = -INFINITY;
  double sqt = std::sqrt(t), temp = -v * a*w - (v * v)*t / 2;
  double vt = v * t;


  // SM: alternative implementation that skips the lognormal calls in logMill  -- not 100% if correct, so commented out
//   const double c = 0.5 * v * v * t;   // constant across iterations
//
// for (int k = K; k >= 0; k--) {
//     double rj_p = a * (2*k + w);
//     double rj_m = a * (2*k + 2.0 - w);
//
//     double pos1 = -rj_p * v + c + lnnorm(-(rj_p - vt) / sqt);
//     double pos2 =  rj_p * v + c + lnnorm(-(rj_p + vt) / sqt);  // sign flips for +vt
//     double neg1 = -rj_m * v + c + lnnorm(-(rj_m - vt) / sqt);
//     double neg2 =  rj_m * v + c + lnnorm(-(rj_m + vt) / sqt);
//
//     fplus  = logsum(logsum(pos1, pos2), fplus);
//     fminus = logsum(logsum(neg1, neg2), fminus);
// }
  for (int k = K; k >= 0; k--)
  {
    double rj = a*(2 * k + w);
    double dj = lognormal(rj / sqt);
    double pos1 = dj + logMill((rj - vt) / sqt);
    double pos2 = dj + logMill((rj + vt) / sqt);
    fplus = logsum(logsum(pos1, pos2), fplus);
    rj = a*(2.0 * k + 2.0 - w);
    dj =  lognormal(rj / sqt);
    double neg1 = dj + logMill((rj - vt) / sqt);
    double neg2 = dj + logMill((rj + vt) / sqt);
    fminus = logsum(logsum(neg1, neg2), fminus);
  }

  return logdiff(fplus, fminus)+temp;
}

// Large-t CDF series: linear-domain accumulation, log returned
static double logFl(double q, double v, double a, double w, int K)
{
  double fplus = -INFINITY, fminus = -INFINITY;
  double la = std::log(a), lv = std::log(fabs(v));
  double F = -INFINITY;
  for (int k = K; k >= 1; k--) {
    double temp0 = std::log(k * 1.0), temp1 = k * M_PI, temp2 = temp1 * w;
    double check = sin(temp2);
    if (check > 0) {
      double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * ((temp1 / a)*(temp1 / a)) * q + std::log(check);
      fplus = logsum(temp, fplus);
    }
    else if (check < 0)
    {
      double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * ((temp1 / a)*(temp1 / a)) * q + std::log(-check);
      fminus = logsum(temp, fminus);
    }
  }
  F = logdiff(fplus, fminus);
  return (F - v * a * w - 0.5 * (v * v) * q);
}

// log-CDF of the Wiener first-passage time (lower boundary)
static double pwiener(double q, double a, double v, double w)
{
  if (std::isinf(q)) return logP(a, v, w);   // t -> inf: return total mass

  const double Kss = Ks(q, v, a, w, DDM_LOG_EPS);
  const double Kll = Kl(q, v, a, w, DDM_LOG_EPS);

  if (3.0 * Kss < Kll)
    return logFs(q, v, a, w, static_cast<int>(Kss));
  else {
    const double lg = M_LN2 + M_LNPI - 2.0 * std::log(a);
    return logdiff(logP(a, v, w), lg + logFl(q, v, a, w, static_cast<int>(Kll)));
  }
}

// =============================================================================
// Integration helpers
// =============================================================================

static void run_hcubature(int (*integrand)(unsigned, const double*, void*, unsigned, double*),
                          my_params& params,
                          int dim, double t, double t0, double st,
                          bool sv_first_dim,
                          double* val, double* err)
{
  double xmin[3] = {}, xmax[3] = {1.0, 1.0, 1.0};

  if (sv_first_dim) {
    xmin[0] = -1.0; xmax[0] = 1.0;
    for (int i = 1; i < dim; ++i) { xmin[i] = 0.0; xmax[i] = 1.0; }
  } else {
    for (int i = 0; i < dim; ++i) { xmin[i] = 0.0; xmax[i] = 1.0; }
  }
  if (st > 0.0) xmax[dim - 1] = std::fmin(1.0, (t - t0) / st);

  hcubature(integrand, &params, dim, xmin, xmax, DDM_NEVAL,
            DDM_EPS * 0.9, 0.0, val, err);
}

static int int_ddiff(unsigned /*dim*/, const double* x, void* p,
                     unsigned /*fdim*/, double* retval)
{
  const auto* params = static_cast<const my_params*>(p);
  const double omega = params->sw ? params->w + params->sw * (x[0] - 0.5) : params->w;
  const double tau   = params->sw ? (params->st ? params->t0 + params->st * x[1] : params->t0)
    : (params->st ? params->t0 + params->st * x[0] : params->t0);

  if (params->t - tau <= 0.0) { retval[0] = 0.0; return 0; }
  // dwiener returns log-density; exponentiate for the integrand
  retval[0] = std::exp(dwiener(params->low_or_up * (params->t - tau),
                               params->a, params->v, omega, params->sv));
  return 0;
}

static void ddiff(double t, int low_or_up, double a, double v, double t0,
                  double w, double sw, double sv, double st,
                  double* derivF, double* Rerr)
{
  my_params params = {t, low_or_up, a, v, t0, w, sw, sv, st};
  const int dim = (sw != 0.0) + (st != 0.0);
  double val, err;
  run_hcubature(int_ddiff, params, dim, t, t0, st, false, &val, &err);
  *derivF = val;
  if (*Rerr < err) *Rerr = err;
}

static int int_pdiff(unsigned /*dim*/, const double* x, void* p,
                     unsigned /*fdim*/, double* retval)
{
  const auto* params = static_cast<const my_params*>(p);
  const double x0sq  = params->sv ? x[0] * x[0] : 0.0;
  const double y     = params->sv ? x[0] / (1.0 - x0sq) : 0.0;
  const double nu    = params->sv ? params->v + params->sv * y : params->v;
  const double omega = params->sv ? (params->sw ? params->w + params->sw * (x[1] - 0.5) : params->w)
    : (params->sw ? params->w + params->sw * (x[0] - 0.5) : params->w);
  const double tau   = params->sv
  ? (params->sw ? (params->st ? params->t0 + params->st * x[2] : params->t0)
       : (params->st ? params->t0 + params->st * x[1] : params->t0))
    : (params->sw ? (params->st ? params->t0 + params->st * x[1] : params->t0)
         : (params->st ? params->t0 + params->st * x[0] : params->t0));

  if (params->t - tau <= 0.0) { retval[0] = 0.0; return 0; }

  const double wn  = (params->low_or_up == 1) ? 1.0 - omega : omega;
  const double lpW = pwiener(params->t - tau, params->a, -params->low_or_up * nu, wn);
  const double jac = params->sv
  ? -0.5 * y * y - M_LN_SQRT_PI - 0.5 * M_LN2
  + std::log1p(x0sq) - 2.0 * std::log1p(-x0sq)
    : 0.0;

  retval[0] = std::exp(lpW + jac);
  return 0;
}

static void pdiff(double t, int low_or_up, double a, double v, double t0,
                  double w, double sw, double sv, double st,
                  double* derivF, double* Rerr)
{
  my_params params = {t, low_or_up, a, v, t0, w, sw, sv, st};
  const int dim = (sv != 0.0) + (sw != 0.0) + (st != 0.0);
  double val, err;
  run_hcubature(int_pdiff, params, dim, t, t0, st, sv != 0.0, &val, &err);
  *derivF = val;
  if (*Rerr < err) *Rerr = err;
}

static inline double compute_sz(double Z, double SZ) {
  return (Z < 1.0 - Z) ? 2.0 * SZ * Z : 2.0 * SZ * (1.0 - Z);
}

// CDF
static inline double ddm_logcdf_scalar(double rt, int R,
                                       double v, double a, double sv,
                                       double t0, double st0, double s,
                                       double Z, double SZ)
{
  const double new_rt = rt - t0;
  if (new_rt <= 0.0) return R_NegInf;   // semantic guard

  if (sv == 0.0 && SZ == 0.0 && st0 == 0.0) {
    double vv = v / s;
    double w  = Z;
    if (R != 1) { vv = -vv; w = 1.0 - w; }
    return pwiener(new_rt, a / s, vv, w);
  }

  // Slow path: same NaN/Inf note as ddm_logpdf_scalar above.
  const double pm = (R == 1) ? -1.0 : 1.0;
  const double sz = compute_sz(Z, SZ);
  double Rval = 0.0, Rerr = 0.0;
  pdiff(rt, static_cast<int>(pm), a / s, v / s, t0,
        Z, sz, sv / s, st0, &Rval, &Rerr);
  return std::log(Rval);
}



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
  const double* __restrict__ Z   = &pt.base(0, spec.col_Z);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ s   = &pt.base(0, spec.col_s);
  const double* __restrict__ SZ  = &pt.base(0, spec.col_SZ);
  const double* __restrict__ st0 = &pt.base(0, spec.col_st0);
  const int n = (int)idx_all.size();

  // int    Epsflag = 1;
  // double eps     = 5e-3;
  // int    K       = 0;

  for (int j = 0; j < n; ++j) {
    const int i       = idx_all[j];
    const double pm   = (R_ptr[i] == 1) ? -1.0 : 1.0;
    const double teff = rt_ptr[i] - t0[i];

    if (SZ[i] == 0.0 && st0[i] == 0.0) {
      if (teff <= 0.0) { ll_row[i] = R_NegInf; continue; }

      double val = dwiener(
        teff * pm,
        a[i]  / s[i],
        v[i]  / s[i],
        Z[i],
        sv[i] / s[i]);
      ll_row[i] = std::isfinite(val) ? val : R_NegInf;
    } else {
      // int    Neval = 6000;
      // int    choice = 0;
      double Rval, Rerr;

      double sz_i = compute_sz(Z[i], SZ[i]);

      ddiff(rt_ptr[i], pm,
            a[i] / s[i],
            v[i] / s[i],
            t0[i], Z[i], sz_i,
            sv[i] / s[i],
            st0[i], &Rval, &Rerr);

      ll_row[i] = (std::isfinite(Rval) && Rval > 0.0) ? std::log(Rval) : R_NegInf;
    }
  }
}


NumericVector d_DDM_Wien(NumericVector rts, IntegerVector Rs, NumericMatrix pars, std::vector<int> is_ok){
  // int Epsflag = 1;
  // double eps = 5e-3;
  // int K = 0;
  // int Neval = 6000;
  // int choice = 0; //the type of integration method to choose.
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
          out[i] = dwiener(new_rt*pm, pars(i, 1)/pars(i,5), pars(i, 0)/pars(i,5), pars(i, 6), pars(i, 2)/pars(i,5));
        } else{
          out[i] = 	R_NegInf;
        }
      } else{ // otherwise use ddiff function with integration
        double Rval;
        double Rerr;
        double sz = (pars(i,6) < (1 - pars(i,6))) ? 2*pars(i,7)*pars(i,6) : 2*pars(i,7)*(1-pars(i,6));
        ddiff(rts[i], pm, pars(i, 1)/pars(i,5), pars(i, 0)/pars(i,5), pars(i, 3), pars(i, 6), sz, pars(i, 2)/pars(i,5), pars(i,4), &Rval, &Rerr);
        out[i] = log(Rval);
      }
    }
  }
  return(out);
}


void ddm_survivor(const std::vector<int>&     idx,
                  const std::vector<double>&  bound,
                  const ParamTable&           pt,
                  const RaceSpec&             spec,
                  double* __restrict__        out,
                  RaceScratch&                scratch)
{
  const double* __restrict__ bd  = bound.data();
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ a   = &pt.base(0, spec.col_a);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);
  const double* __restrict__ st0 = &pt.base(0, spec.col_st0);
  const double* __restrict__ s   = &pt.base(0, spec.col_s);
  const double* __restrict__ Z   = &pt.base(0, spec.col_Z);
  const double* __restrict__ SZ  = &pt.base(0, spec.col_SZ);

  for (int j = 0; j < (int)idx.size(); ++j) {
    const int i = idx[j];

    if (bd[i] - t0[i] <= 0.0) continue;  // bound before t0 — no mass yet, log(1) = 0

    const double logcdf_upper = ddm_logcdf_scalar(bd[i], 1, v[i], a[i], sv[i], t0[i], st0[i], s[i], Z[i], SZ[i]);
    const double logcdf_lower = ddm_logcdf_scalar(bd[i], 2, v[i], a[i], sv[i], t0[i], st0[i], s[i], Z[i], SZ[i]);
    double S = 1 - (std::exp(logcdf_upper) + std::exp(logcdf_lower));

    // clamp - probability so between 0 and 1, anything else is integration error or numerical precision error
    if (!std::isfinite(S) || S < 0.0) S = 0.0;
    else if (S > 1.0) S = 1.0;
    out[i] = S;
    }
}


void ddm_survivor_with_response(const std::vector<int>&    idx,
                                const std::vector<int>&    winner,
                                const std::vector<double>& lower,
                                const std::vector<double>& upper,
                                int                        n_acc,
                                const ParamTable&          pt,
                                const RaceSpec&            spec,
                                double* __restrict__       out)
{
  // DDM has exactly 2 boundaries: winner=0 → upper (R=1), winner=1 → lower (R=2).
  // P(winner hits first in [lo, hi]) = CDF_winner(hi) - CDF_winner(lo).
  // No race integral needed — the two accumulators are the two boundaries of a
  // single diffusion process, so the winning boundary is determined by R directly.

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
    const int    i  = idx[j];           // base row = trial * n_acc, n_acc is always 1 for DDM
    const int    R  = winner[j] + 1;    // 0-based winner → 1-based R (1=upper, 2=lower)
    const double lo = lower[j];
    const double hi = upper[j];

    // CDF at upper bound
    double cdf_hi = 0.0;
    if (hi - t0[i] > 0.0)
      cdf_hi = std::exp(ddm_logcdf_scalar(hi, R, v[i], a[i], sv[i],
                                          t0[i], st0[i], s[i], Z[i], SZ[i]));

    // CDF at lower bound
    double cdf_lo = 0.0;
    if (lo - t0[i] > 0.0)
      cdf_lo = std::exp(ddm_logcdf_scalar(lo, R, v[i], a[i], sv[i],
                                          t0[i], st0[i], s[i], Z[i], SZ[i]));

    double p = cdf_hi - cdf_lo;
    if (!std::isfinite(p) || p < 0.0) p = 0.0;
    else if (p > 1.0)                  p = 1.0;
    out[j] = p;
  }
}
