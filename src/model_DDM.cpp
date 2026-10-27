#include <Rcpp.h>
using namespace Rcpp;
#include "model_DDM.h"
#include "composite_functions.h"

// =============================================================================
// Constants
// =============================================================================
static constexpr int    DDM_NEVAL   = 6000;
static constexpr double DDM_EPS     = 5e-3;
static constexpr double DDM_LOG_EPS = -5.298317366548036311258; //std::log(DDM_EPS);
static constexpr double DDM_CDF_LOG_EPS = -27.631021115928547651; // log(1e-12)
static constexpr double DDM_LOG_INNER_EPS = -7.600902459542082361471; //std::log(DDM_EPS * 0.1);

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

// The diffusion series contains signed terms; retain their absolute difference.
static double log_absdiff_exp(double a,double b) {
  return log_diff_exp(std::max(a,b),std::min(a,b));
}

static double lnnorm(double z) { return PNORM_STD(z,true,true); }

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
  double sum_plus  = w;
  double sum_minus = 0.0;

#pragma omp simd reduction(+:sum_plus,sum_minus)
  for (int k = 1; k <= K; ++k) {
    const double t1 = w + 2.0 * k;
    const double t2 = w - 2.0 * k;
    sum_plus  += t1    * std::exp(-(t1*t1-w*w)*inv_twot);
    sum_minus += (-t2) * std::exp(-(t2*t2-w*w)*inv_twot);
  }

  const double S = sum_plus - sum_minus;
  if (!(S > 0.0)) return R_NegInf;
  return -0.5 * M_LN2 - M_LN_SQRT_PI - 1.5 * std::log(t) + std::log(S) - w*w*inv_twot;
}

// Large-t series: linear-domain accumulation, log returned
static double logfl_linear(double q, double w, int K) {
  if (w == 0.0 || K <= 0) return R_NegInf;

  const double halfq = q / 2.0;
  double sum_plus = 0.0, sum_minus = 0.0;

#pragma omp simd reduction(+:sum_plus,sum_minus)
  for (int k = 1; k <= K; ++k) {
    const double kpi  = k * M_PI;
    const double base = k * std::sin(kpi * w) * std::exp(-(kpi*kpi-M_PI*M_PI)*halfq);
    sum_plus  += (base > 0.0) ?  base : 0.0;
    sum_minus += (base < 0.0) ? -base : 0.0;
  }

  const double S = sum_plus - sum_minus;
  if (!(S > 0.0)) return R_NegInf;
  return std::log(S) + M_LNPI - M_PI*M_PI*halfq;
}

// log-density of the Wiener first-passage time
static double dwiener(double q, double a, double vn, double wn, double sv,
                      double log_eps = DDM_LOG_EPS)
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

    const double es  = log_eps - lg1;
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
  const double z=2*v*a;
  if (z==0.0) return std::log1p(-w);
  return log1m_exp(-std::fabs(z)*(1-w))-log1m_exp(-std::fabs(z))
         -std::max(z,0.0)*w;
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

// Log CDF from the short-time series.
static double logFs(double t, double v, double a, double w, int K)
{
  // Cancel the Gaussian exponents analytically before evaluating normal tails.
  double fplus = -INFINITY, fminus = -INFINITY;
  const double sqt  = std::sqrt(t);
  const double vt   = v * t;
  const double c    = 0.5 * v * v * t;
  const double temp = -v * a * w - c;

  for (int k = K; k >= 0; k--) {
    // Positive terms: r = a*(2k + w)
    const double rp   = a * (2.0 * k + w);
    const double pos1 = (-rp * v + c) + lnnorm(-(rp - vt) / sqt);
    const double pos2 = (+rp * v + c) + lnnorm(-(rp + vt) / sqt);
    fplus = log_sum_exp(log_sum_exp(pos1, pos2), fplus);

    // Negative terms: r = a*(2k + 2 - w)
    const double rm   = a * (2.0 * k + 2.0 - w);
    const double neg1 = (-rm * v + c) + lnnorm(-(rm - vt) / sqt);
    const double neg2 = (+rm * v + c) + lnnorm(-(rm + vt) / sqt);
    fminus = log_sum_exp(log_sum_exp(neg1, neg2), fminus);
  }

  return log_absdiff_exp(fplus, fminus)+temp;
}

// Log remaining boundary mass from the large-time series.
static double logFl(double q, double v, double a, double w, int K)
{
  double fplus = -INFINITY, fminus = -INFINITY;
  double la = std::log(a), lv = std::log(fabs(v));
  const double half_pia_sq_q  = 0.5 * (M_PI / a) * (M_PI / a) * q;  // (k*pi/a)^2 * q/2, scaled by k^2 in loop

  for (int k = K; k >= 1; k--) {
    const double temp0  = std::log((double)k);
    const double kpi    = k * M_PI;
    const double check  = std::sin(kpi * w);
    if (check == 0.0) continue;

    const double temp = temp0 - log_sum_exp(2.0 * lv, 2.0 * (temp0 + M_LNPI - la))
                              - (double)(k * k) * half_pia_sq_q
                              + std::log(std::fabs(check));
    if (check > 0.0) fplus  = log_sum_exp(temp, fplus);
    else             fminus = log_sum_exp(temp, fminus);
    }
  return log_absdiff_exp(fplus, fminus) - v * a * w - 0.5 * v * v * q;
}

// log-CDF of the Wiener first-passage time (lower boundary)
static double pwiener(double q, double a, double v, double w)
{
  if (std::isinf(q)) return logP(a, v, w);   // t -> inf: return total mass

  const double Kss = Ks(q, v, a, w, DDM_CDF_LOG_EPS);
  const double Kll = Kl(q, v, a, w, DDM_CDF_LOG_EPS);

  if (3.0 * Kss < Kll)
    return logFs(q, v, a, w, static_cast<int>(Kss));
  else {
    const double lg = M_LN2 + M_LNPI - 2.0 * std::log(a);
    return log_absdiff_exp(logP(a, v, w), lg + logFl(q, v, a, w, static_cast<int>(Kll)));
  }
}

static double wiener_logremaining(double q,double a,double v,double w) {
  if (q<=0.0) return logP(a,v,w);
  if (is_pos_inf(q)) return neg_inf();
  const double scaled_time=(q/a)/a;
  if (scaled_time>=0.1) {
    const int K=static_cast<int>(std::ceil(std::sqrt(1+80/(M_PI*M_PI*scaled_time))));
    return M_LN2+M_LNPI-2*std::log(a)+logFl(q,v,a,w,K);
  }
  return log_diff_exp(logP(a,v,w),pwiener(q,a,v,w));
}

// Positive interval mass for one boundary. Use CDFs at early times and
// remaining mass at late times, before averaging over parameter variability.
static double wiener_log_interval(double lower,double upper,double a,double v,double w) {
  lower=std::max(0.0,lower);
  if (!(upper>lower)) return neg_inf();
  if (lower==0.0) return pwiener(upper,a,v,w);
  if (is_pos_inf(upper)) return wiener_logremaining(lower,a,v,w);
  const bool tail=(lower/a)/a>=0.1;
  const double high=tail ? wiener_logremaining(lower,a,v,w) : pwiener(upper,a,v,w);
  const double low=tail ? wiener_logremaining(upper,a,v,w) : pwiener(lower,a,v,w);
  if (is_neg_inf(low)) return high;
  if (high>low && !probability_roundoff_exceeds(-std::expm1(low-high),1.0))
    return log_diff_exp(high,low);
  // A narrow interval can lose accuracy in either difference. Integrate the
  // existing density only in that case, using the same scaling helper.
  const double width=upper-lower;
  return std::log(width)+log_integral_unit([&](double u) {
    return dwiener(-(lower+u*width),a,v,w,0.0,DDM_CDF_LOG_EPS);
  },-1.0);
}

// =============================================================================
// Integration helpers
// =============================================================================

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
                               params->a, params->v, omega, params->sv, DDM_LOG_INNER_EPS));
  return 0;
}

static void ddiff(double t, int low_or_up, double a, double v, double t0,
                  double w, double sw, double sv, double st,
                  double* derivF, double* Rerr)
{
  my_params params = {t, low_or_up, a, v, t0, w, sw, sv, st};
  const int dim = (sw != 0.0) + (st != 0.0);
  double val, err;
  double lo[2]={0,0},hi[2]={1,1};
  if (st>0.0) hi[dim-1]=std::min(1.0,(t-t0)/st);
  hcubature(int_ddiff,&params,dim,lo,hi,DDM_NEVAL,DDM_EPS*0.9,0,&val,&err);
  *derivF = val;
  *Rerr = err;
}

static inline double compute_sz(double Z, double SZ) {
  return (Z < 1.0 - Z) ? 2.0 * SZ * Z : 2.0 * SZ * (1.0 - Z);
}

static inline double ddm_logpdf_scalar(double rt,int response,double v,double a,
                                       double sv,double t0,double st0,double s,
                                       double Z,double SZ) {
  const double q=rt-t0;
  if (q<=0.0) return neg_inf();
  const int side=response==1 ? -1 : 1;
  if (SZ==0.0 && st0==0.0) return dwiener(side*q,a/s,v/s,Z,sv/s);
  double value=0.0,error=0.0;
  ddiff(rt,side,a/s,v/s,t0,Z,compute_sz(Z,SZ),sv/s,st0,&value,&error);
  return std::log(value);
}

void fill_ddm(const double* rt_ptr,
              const int* R_ptr,
              const ParamTable& pt,
              const RaceSpec& spec,
              const std::vector<int>& idx_all,
              double* __restrict__ ll_row,
              const std::vector<int>& is_ok)
{
  const double* __restrict__ v   = &pt.base(0, spec.col_v);
  const double* __restrict__ a   = &pt.base(0, spec.col_a);
  const double* __restrict__ t0  = &pt.base(0, spec.col_t0);
  const double* __restrict__ Z   = &pt.base(0, spec.col_Z);
  const double* __restrict__ sv  = &pt.base(0, spec.col_sv);
  const double* __restrict__ s   = &pt.base(0, spec.col_s);
  const double* __restrict__ SZ  = &pt.base(0, spec.col_SZ);
  const double* __restrict__ st0 = &pt.base(0, spec.col_st0);
  const int n = (int)idx_all.size();

  for (int j=0; j<n; ++j) {
    const int i=idx_all[j];
    if (!is_ok[i]) continue;
    const double value=ddm_logpdf_scalar(rt_ptr[i],R_ptr[i],v[i],a[i],sv[i],
                                         t0[i],st0[i],s[i],Z[i],SZ[i]);
    ll_row[i]=is_finite(value) ? value : neg_inf();
  }
}

NumericVector d_DDM_Wien(NumericVector rts,IntegerVector Rs,NumericMatrix pars) {
  NumericVector out(rts.size());
  for (int i=0; i<rts.size(); ++i)
    out[i]=ddm_logpdf_scalar(rts[i],Rs[i],pars(i,0),pars(i,1),pars(i,2),
                             pars(i,3),pars(i,4),pars(i,5),pars(i,6),pars(i,7));
  return out;
}

// response 0 integrates both boundaries; 1 and 2 select one. This is shared
// by truncation, known-response censoring and unknown-response censoring.
double ddm_log_response_interval(const ParamTable& pt,const RaceSpec& spec,int row,
                                 int response,double lower,double upper) {
  const double s=pt.base(row,spec.col_s),a=pt.base(row,spec.col_a)/s;
  const double v=pt.base(row,spec.col_v)/s,sv=pt.base(row,spec.col_sv)/s;
  const double t0=pt.base(row,spec.col_t0),st=pt.base(row,spec.col_st0);
  const double z=pt.base(row,spec.col_Z),sz=compute_sz(z,pt.base(row,spec.col_SZ));
  lower=std::max(lower,t0);
  if (!(upper>lower)) return neg_inf();
  if (response==0 && lower==t0 && is_pos_inf(upper)) return 0.0;
  auto point=[&](double drift,double start,double onset) {
    const double lo=lower-onset,hi=upper-onset;
    if (response==1) return wiener_log_interval(lo,hi,a,drift,start);
    if (response==2) return wiener_log_interval(lo,hi,a,-drift,1-start);
    return log_sum_exp(wiener_log_interval(lo,hi,a,drift,start),
                       wiener_log_interval(lo,hi,a,-drift,1-start));
  };
  if (sv==0.0 && sz==0.0 && st==0.0) return point(v,z,t0);
  // The density averages normal drift analytically. Integrating decision
  // time against the overlap with the uniform nondecision interval also
  // removes that variability dimension exactly. Only start and time remain.
  lower-=t0;
  upper-=t0;
  std::vector<double> breaks{std::max(0.0,lower-st),upper};
  if (lower>breaks.front() && lower<upper) breaks.push_back(lower);
  if (upper-st>breaks.front() && upper-st<upper) breaks.push_back(upper-st);
  std::sort(breaks.begin(),breaks.end());
  // Relative error is controlled in the requested mass, including tiny CDFs.
  constexpr double tolerance=1e-7;
  double result=neg_inf();
  for (std::size_t k=1;k<breaks.size();++k) {
    const double left=breaks[k-1],right=breaks[k];
    const bool infinite=is_pos_inf(right);
    // Diffusion tails decay on the a^2 time scale. A scale proportional to
    // the absolute bound would hide the mass near a very late lower limit.
    const bool transformed=infinite || right-left>a*a;
    const double width=transformed ? a*a : right-left;
    const double end=infinite || !transformed ? 1.0 : (right-left)/(width+right-left);
    double scale=0.0,largest=neg_inf();
    auto probability=[&](const double* x) {
      const double den=transformed ? 1-end*x[0] : 1.0;
      if (den==0.0) return 0.0;
      const double time=left+width*end*x[0]/den;
      const double overlap=st>0.0
        ? std::max(0.0,std::min(st,upper-time)-std::max(0.0,lower-time))/st : 1.0;
      if (overlap==0.0) return 0.0;
      const double log_weight=std::log(overlap*width*end)-2*std::log(den);
      auto density=[&](const double* u) {
        const double start=sz>0.0 ? z+sz*(u[0]-0.5) : z;
        double lp=neg_inf();
        if (response!=2) lp=dwiener(-time,a,v,start,sv,DDM_CDF_LOG_EPS);
        if (response!=1) lp=log_sum_exp(lp,dwiener(time,a,v,start,sv,DDM_CDF_LOG_EPS));
        lp+=log_weight;
        largest=std::max(largest,lp);
        return std::exp(lp-scale);
      };
      if (sz==0.0) return density(x);
      using Density=decltype(density);
      auto average=[](unsigned,const double* u,void* p,unsigned,double* out) -> int {
        out[0]=(*static_cast<Density*>(p))(u);return 0;
      };
      double a0=0,b0=1,value=0,error=0;
      hcubature(average,&density,1,&a0,&b0,6000,0,tolerance/10,&value,&error);
      return value;
    };
    using Function=decltype(probability);
    auto integrand=[](unsigned,const double* x,void* p,unsigned,double* out) -> int {
      out[0]=(*static_cast<Function*>(p))(x);return 0;
    };
    double lo=0,hi=1,value=0.0,error=0.0;
    hcubature(integrand,&probability,1,&lo,&hi,6000,0,tolerance,&value,&error);
    if (value<1e-280 && is_finite(largest)) {
      scale=largest;
      hcubature(integrand,&probability,1,&lo,&hi,6000,0,tolerance,&value,&error);
    }
    if (error>tolerance*value)
      hcubature(integrand,&probability,1,&lo,&hi,100000,0,tolerance,&value,&error);
    result=log_sum_exp(result,scale+std::log(value));
  }
  return result;
}
