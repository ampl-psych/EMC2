#ifndef WALD_FUNCTIONS_H
#define WALD_FUNCTIONS_H

#include <cmath>
#include <algorithm>
#include <Rcpp.h>
#include "pnorm_utils.h"
#include "composite_functions.h"
#include "hcubature.h"

constexpr double A_EPS = 1e-4;
constexpr double L_EPS = 1e-4;
constexpr double K_MAX = 1e6;
inline void clamp_l(double& l) {
  l=(l>-L_EPS && l<L_EPS) ? (l>=0 ? L_EPS : -L_EPS) : l;
}
inline void clamp_a_l(double& a,double& l) { a=std::max(a,A_EPS); clamp_l(l); }

inline double wald_limit(double k,double l,double a,bool lower_tail,bool log_p) {
  double lc=0.0;
  if (l<0.0) {
    lc=2*l*k;
    if (a>=A_EPS) lc=2*l*(k-a)+std::log(-std::expm1(4*l*a))-std::log(-4*l*a);
  }
  const double lp=lower_tail ? lc : log1m_exp(lc);
  return log_p ? lp : std::exp(lp);
}

inline double digt0(double t,double k,double l,bool log_p=false) {
  if (is_pos_inf(t)) return log_p ? neg_inf() : 0.0;
  const double z=(k-t*l)/std::sqrt(t);
  if (log_p) return (-0.5*z)*z+std::log(k)-1.5*std::log(t)-pnorm_detail::LOG_SQRT_2PI;
  return pnorm_detail::INV_SQRT_2PI*pnorm_detail::gaussian_exponential(z)*k/(t*std::sqrt(t));
}

// Exp(2*l*d)*Phi(-(d+t*l)/sqrt(t)), without overflow or cancellation
// between large logarithms. z=(d-t*l)/sqrt(t), w=(d+t*l)/sqrt(t).
inline double wald_weighted_tail(double d,double l,double z,double w,bool log_p) {
  if (w>=0.0) {
    const double h=pnorm_detail::normal_tail_factor(w);
    return log_p ? (-0.5*z)*z+std::log(h) : pnorm_detail::gaussian_exponential(z)*h;
  }
  return log_p ? 2*l*d+PNORM_STD(-w,true,true)
               : std::exp(2*l*d)*PNORM_STD(-w,true,false);
}

inline double pigt0(double t,double k,double l,bool lower_tail=true,bool log_p=false) {
  if (is_pos_inf(t)) return wald_limit(k,l,0,lower_tail,log_p);
  const double rt=std::sqrt(t), z=(k-t*l)/rt, w=(k+t*l)/rt;
  if (lower_tail) {
    if (log_p && z<=0) return log1m_exp(pigt0(t,k,l,false,true));
    if (log_p) return log_sum_exp(PNORM_STD(-z,true,true),wald_weighted_tail(k,l,z,w,true));
    return std::min(1.0,PNORM_STD(-z,true,false)+wald_weighted_tail(k,l,z,w,false));
  }
  if (log_p && z>0.67448975) {
    const double lc=log_sum_exp(PNORM_STD(-z,true,true),wald_weighted_tail(k,l,z,w,true));
    if (lc < -composite_detail::LOG2) return log1m_exp(lc);
  }
  if (l>=0.0 && z<=0.67448975) {
    const double h=pnorm_detail::normal_tail_difference(-z,2*k/rt);
    return log_p ? (-0.5*z)*z+std::log(h) : pnorm_detail::gaussian_exponential(z)*h;
  }
  if (l<0.0) {
    const double lp=log_sum_exp(pnorm_detail::normal_interval_width(-w,2*k/rt,true),
                               log1m_exp(2*l*k)+PNORM_STD(-w,true,true));
    return log_p ? lp : std::exp(lp);
  }
  const double p=PNORM_STD(z,true,false)-wald_weighted_tail(k,l,z,w,false);
  return log_p ? std::log(p) : p;
}

namespace wald_detail {
// The closed forms can cancel for a narrow start range or extreme tails.
// Integrate the same point-start distribution only in those ill-conditioned
// cases, scaling out its largest log value before using existing cubature.
inline double average(double t,double k,double l,double a,int quantity,bool log_p) {
  double d=quantity==1 ? k-a : k+a;
  if (quantity==0) {
    const double drift=t*l;
    const double root=std::hypot(drift,2*std::sqrt(t));
    d=std::clamp(drift>=0 ? 0.5*(drift+root) : 2*t/(root-drift),k-a,k+a);
  }
  const double lp=log_integral_unit([&](double u) {
    const double distance=k+a*(2*u-1);
    return quantity==0 ? digt0(t,distance,l,true) : pigt0(t,distance,l,quantity==1,true);
  },std::clamp((d-(k-a))/(2*a),0.0,1.0));
  return log_p ? lp : std::exp(lp);
}
}

inline double digt_core(double t,double k,double l,double a,bool log_p=false) {
  if (is_pos_inf(t)) return log_p ? neg_inf() : 0.0;
  const double rt=std::sqrt(t), za=(k-a-t*l)/rt, zb=(k+a-t*l)/rt;
  const double delta=-2*a*(k-t*l)/t;
  const double ga=pnorm_detail::gaussian_exponential(za), gb=pnorm_detail::gaussian_exponential(zb);
  if (!log_p && ga==0.0 && gb==0.0 && za*zb>0.0 &&
      -0.5*std::min(za*za,zb*zb)+std::log(k+a)-1.5*std::log(t)
        -pnorm_detail::LOG_SQRT_2PI < -746.0) return 0.0;
  const double dg=pnorm_detail::INV_SQRT_2PI*(delta<=0 ? ga*(-std::expm1(delta)) : -gb*(-std::expm1(-delta)))/rt;
  const bool lower=za+zb<0;
  const double pa=pnorm_detail::normal_probability(za,ga,lower);
  const double pb=pnorm_detail::normal_probability(zb,gb,lower);
  const double dn=l*(lower ? pb-pa : pa-pb);
  const double sum=dg+dn;
  if (probability_roundoff_exceeds(sum,std::fabs(dg)+std::fabs(l)*(pa+pb)))
    return wald_detail::average(t,k,l,a,0,log_p);
  return log_p ? std::log(sum)-std::log(2*a) : sum/(2*a);
}

inline double pigt_core(double t,double k,double l,double a,bool lower_tail=true,bool log_p=false) {
  if (is_pos_inf(t)) return wald_limit(k,l,a,lower_tail,log_p);
  const double rt=std::sqrt(t), za=(k-a-t*l)/rt, zb=(k+a-t*l)/rt;
  const double delta=-2*a*(k-t*l)/t;
  const double ga=pnorm_detail::gaussian_exponential(za), gb=pnorm_detail::gaussian_exponential(zb);
  const double dg=pnorm_detail::INV_SQRT_2PI*(delta<=0 ? ga*(-std::expm1(delta)) : -gb*(-std::expm1(-delta)));
  const bool lower=za+zb<0;
  const double pa=pnorm_detail::normal_probability(za,ga,lower);
  const double pb=pnorm_detail::normal_probability(zb,gb,lower);
  const double dn=lower ? pb-pa : pa-pb;
  const double wa=(k-a+t*l)/rt, wb=(k+a+t*l)/rt;
  const double ta=wald_weighted_tail(k-a,l,za,wa,false), tb=wald_weighted_tail(k+a,l,zb,wb,false);
  if (l==0.0) return wald_detail::average(t,k,l,a,lower_tail ? 1 : 2,log_p);
  const double t1=rt*dg, t2=(tb-ta)/(2*l), t4=(t*l-k+0.5/l)*dn;
  const double mean=0.5*(pa+pb);
  const double base=lower==!lower_tail ? mean : 1.0-mean;
  const double correction=(t1+t2+t4)/(2*a);
  const double p=base+(lower_tail ? correction : -correction);
  const double magnitude=std::fabs(base)+(std::fabs(t1)+(std::fabs(ta)+std::fabs(tb))/(2*std::fabs(l))
      +std::fabs(t*l-k+0.5/l)*(pa+pb))/(2*a);
  if (p>0.5 && p<=1.0 && log_p)
    return log1m_exp(pigt_core(t,k,l,a,!lower_tail,true));
  if (p>1.0 || probability_roundoff_exceeds(p,magnitude))
    return wald_detail::average(t,k,l,a,lower_tail ? 1 : 2,log_p);
  return log_p ? std::log(p) : p;
}

inline double digt(double t,double k,double l,double a,bool log_p=false) {
  if (t<=0.0 || is_pos_inf(t)) return log_p ? neg_inf() : 0.0;
  if (a<A_EPS) return digt0(t,k,l,log_p);
  clamp_a_l(a,l);
  return digt_core(t,k,l,a,log_p);
}
inline double pigt(double t,double k,double l,double a,bool lower_tail=true,bool log_p=false) {
  if (t<=0.0) {
    const bool one=!lower_tail;
    return log_p ? (one ? 0.0 : neg_inf()) : double(one);
  }
  if (a<A_EPS) return pigt0(t,k,l,lower_tail,log_p);
  clamp_a_l(a,l);
  return pigt_core(t,k,l,a,lower_tail,log_p);
}
Rcpp::NumericVector dWald(Rcpp::NumericVector t,Rcpp::NumericVector v,Rcpp::NumericVector B,Rcpp::NumericVector A,Rcpp::NumericVector t0);
Rcpp::NumericVector pWald(Rcpp::NumericVector t,Rcpp::NumericVector v,Rcpp::NumericVector B,Rcpp::NumericVector A,Rcpp::NumericVector t0);
#endif
