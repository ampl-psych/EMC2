#ifndef model_LBA_h
#define model_LBA_h

#include <Rcpp.h>
#include "RaceSpec.h"
#include "utility_functions.h"
#include "math_utils.h"
#include "pnorm_utils.h"
#include "composite_functions.h"
#include "ParamTable.h"
#include "hcubature.h"

constexpr double LBA_A_ASYMPTOTIC = 1e-10;

// R exports — declared here, defined in model_LBA.cpp
Rcpp::NumericVector dlba(Rcpp::NumericVector t,
                         Rcpp::NumericVector A, Rcpp::NumericVector b,
                         Rcpp::NumericVector v, Rcpp::NumericVector sv,
                         bool posdrift);

Rcpp::NumericVector plba(Rcpp::NumericVector t,
                         Rcpp::NumericVector A, Rcpp::NumericVector b,
                         Rcpp::NumericVector v, Rcpp::NumericVector sv,
                         bool posdrift);

// ---------------------------------------------------------------------------
// Split compute functions — small, always inlined, no ODR issue in header
// ---------------------------------------------------------------------------

namespace lba_detail {
// Uniform average of the point-start distribution, used when the closed
// form loses precision. Scaling retains tails below the probability range.
struct Average {
  double t,v,sv,denom;
  int quantity;
  bool posdrift;
  double log_value(double d) const {
    const double z=(d/t-v)/sv;
    if (quantity==0)
      return (-0.5*z)*z-pnorm_detail::LOG_SQRT_2PI+std::log(d)
             -2*std::log(t)-std::log(sv)-std::log(denom);
    if (quantity==1) return PNORM_STD(z,false,true)-std::log(denom);
    const double residual=posdrift ? denom-PNORM_STD(v/sv,true,false) : PNORM_STD(v/sv,false,false);
    return log_sum_exp(pnorm_detail::normal_interval_width(-v/sv,d/(t*sv),true),
                       residual>0 ? std::log(residual) : neg_inf())-std::log(denom);
  }
};
inline double average(double t,double A,double b,double v,double sv,double denom,
                      int quantity,bool log_p,bool posdrift=true) {
  Average w{t,v,sv,denom,quantity,posdrift};
  double d=quantity==1 ? b-A : b;
  if (quantity==0) {
    const double root=std::hypot(v,2*sv);
    d=std::clamp(t*(v>=0 ? 0.5*(v+root) : 2*sv*sv/(root-v)),b-A,b);
  }
  const double lp=log_integral_unit([&](double u) {
    return w.log_value(b-A+u*A);
  },std::clamp((d-(b-A))/A,0.0,1.0));
  return log_p ? lp : std::exp(lp);
}
}

inline double dlba_core(double t,double A,double b,double v,double sv,double denom,bool log_p=false) {
  const double hi=(b/t-v)/sv, lo=((b-A)/t-v)/sv;
  const double ga=pnorm_detail::gaussian_exponential(lo),gb=pnorm_detail::gaussian_exponential(hi);
  // Bound the density before skipping a raw-probability underflow. The
  // log path still evaluates these tails for conditional likelihoods.
  if (!log_p && ga==0.0 && gb==0.0 && lo*hi>0.0 &&
      -0.5*std::min(lo*lo,hi*hi)+std::log(b)-2*std::log(t)
        -std::log(sv*denom)-pnorm_detail::LOG_SQRT_2PI < -746.0) return 0.0;
  const double delta=-0.5*(A/(t*sv))*(hi+lo);
  const double dg=pnorm_detail::INV_SQRT_2PI*(delta<=0 ? ga*(-std::expm1(delta)) : -gb*(-std::expm1(-delta)));
  const bool lower=lo+hi<0;
  const double pa=pnorm_detail::normal_probability(lo,ga,lower);
  const double pb=pnorm_detail::normal_probability(hi,gb,lower);
  const double dn=v*(lower ? pb-pa : pa-pb);
  const double sum=dn+sv*dg;
  if (probability_roundoff_exceeds(sum,std::fabs(v)*(pa+pb)+std::fabs(sv*dg)))
    return lba_detail::average(t,A,b,v,sv,denom,0,log_p);
  return log_p ? std::log(sum)-std::log(A)-std::log(denom) : sum/(A*denom);
}
inline double dlba_noA(double t,double b,double v,double sv,double denom,bool log_p=false) {
  if (log_p) return LOG_DNORM(b/t,v,sv)+std::log(b)-2*std::log(t)-std::log(denom);
  return DNORM(b/t,v,sv)*b/(t*t*denom);
}

inline double plba_noA(double t,double b,double v,double sv,double denom,
                       bool lower_tail=true,bool log_p=false,bool posdrift=true) {
  const double z=(b/t-v)/sv;
  if (lower_tail) {
    return log_p ? PNORM_STD(z,false,true)-std::log(denom) : PNORM_STD(z,false,false)/denom;
  }
  const double negative=posdrift ? (denom>1e-10 ? PNORM_STD(v/sv,false,false) : 1-denom) : 0.0;
  const double p=(PNORM_STD(z,true,false)-negative)/denom;
  if (p>0.01 && denom>0.1) return log_p ? std::log1p(-PNORM_STD(z,false,false)/denom) : p;
  if (!log_p)
    return (pnorm_detail::normal_interval_width(-v/sv,b/(t*sv))
            +(posdrift ? denom-PNORM_STD(v/sv,true,false) : PNORM_STD(v/sv,false,false)))/denom;
  lba_detail::Average w{t,v,sv,denom,2,posdrift};
  return w.log_value(b);
}
// J(x)=x*Phi(x)+phi(x). In the negative tail use the derivative
// of the scaled normal tail instead of subtracting phi(x) and |x|*Phi(x).
inline double lba_normal_primitive(double x) {
  const double y=std::fabs(x);
  if (y<=0.67448975) return x*PNORM_STD(x,true,false)+DNORM_STD(x);
  double derivative;
  pnorm_detail::normal_tail_factor(y,&derivative);
  const double tail=-pnorm_detail::gaussian_exponential(y)*derivative;
  return x>0 ? x+tail : tail;
}
inline double plba_core(double t,double A,double b,double v,double sv,double denom,
                        bool lower_tail=true,bool log_p=false,bool posdrift=true) {
  const double width=A/(t*sv), hi=(b/t-v)/sv, lo=((b-A)/t-v)/sv;
  // Evaluate the smaller unconditioned normal average first.
  const bool small_cdf=lo+hi>=0.0;
  double small,magnitude;
  if (std::fabs(lo)<4 && std::fabs(hi)<4 && width>0.1) {
    const double ga=pnorm_detail::gaussian_exponential(lo),gb=pnorm_detail::gaussian_exponential(hi);
    const double a=hi*pnorm_detail::normal_probability(hi,gb,!small_cdf);
    const double b=lo*pnorm_detail::normal_probability(lo,ga,!small_cdf);
    const double dg=pnorm_detail::INV_SQRT_2PI*(ga-gb);
    small=(a-b+(small_cdf ? dg : -dg))/width;
    magnitude=(std::fabs(a)+std::fabs(b)+std::fabs(dg))/width;
  } else {
    const double first=lba_normal_primitive(small_cdf ? -lo : hi);
    const double second=lba_normal_primitive(small_cdf ? -hi : lo);
    small=(first-second)/width;
    magnitude=(first+second)/width;
  }
  const double cdf=(small_cdf ? small : 1-small)/denom;
  const double negative=posdrift && !small_cdf ? (denom>1e-10 ? PNORM_STD(v/sv,false,false) : 1-denom) : 0.0;
  const double surv=small_cdf ? 1-cdf : (small-negative)/denom;
  const double p=lower_tail ? cdf : surv;
  const double cdf_error=(magnitude+(small_cdf ? 0.0 : 1.0))/denom;
  const double surv_error=(magnitude+(small_cdf ? 0.0 : negative))/denom;
  const double roundoff_scale=lower_tail ? cdf_error : surv_error;
  if (p>1 || probability_roundoff_exceeds(p,roundoff_scale)) {
    if (p>0.5) {
      const double other=lba_detail::average(t,A,b,v,sv,denom,lower_tail ? 2 : 1,true,posdrift);
      return log_p ? log1m_exp(other) : -std::expm1(other);
    }
    return lba_detail::average(t,A,b,v,sv,denom,lower_tail ? 1 : 2,log_p,posdrift);
  }
  if (!log_p) return p;
  const double other=lower_tail ? surv : cdf;
  if (p>0.5 && probability_roundoff_exceeds(other,lower_tail ? surv_error : cdf_error))
    return log1m_exp(lba_detail::average(t,A,b,v,sv,denom,lower_tail ? 2 : 1,true,posdrift));
  return p>0.5 ? std::log1p(-other) : std::log(p);
}

// ---------------------------------------------------------------------------
// Function declarations — defined in model_LBA.cpp
// ---------------------------------------------------------------------------
void dlba_plba_fast(const double*           rt,
                    const ParamTable&       pt,
                    const RaceSpec&         spec,
                    const std::vector<int>& idx_win,
                    const std::vector<int>& idx_los,
                    double* __restrict__    ll_row,
                    RaceScratch&            scratch,
                    const std::vector<int>& is_ok, bool log_p = false);

void lba_survivor(const std::vector<int>& idx,
                  const std::vector<double>& bound,
                  const ParamTable& pt,
                  const RaceSpec& spec,
                  double* __restrict__ out,
                  RaceScratch& scratch,
                  const std::vector<int>& is_ok, bool log_p = false);

#endif // model_LBA_h
