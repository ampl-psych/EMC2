#ifndef exgaussian_functions_h
#define exgaussian_functions_h

#include <cmath>
#include <algorithm>
#include "composite_functions.h"
#include "pnorm_utils.h"

// Keep the existing lower scale bound. All routes use the same distribution,
// including the parameter-only preparation used by stop-signal quadrature.
static constexpr double SIG_TAU_EPS = 1e-12;

struct ExGaussian {
  double mu, sigma, tau, inv_sigma, ratio, log_tau;
  double lower, upper, log_normalizer;
  enum Approximation { Exact, Normal, ShiftedExponential };
  Approximation approximation;

  ExGaussian(double m, double s, double t,
             double lo=neg_inf(), double hi=pos_inf(),
             double normal_ratio=0.0, double shifted_exp_sigma=0.0)
    : mu(m), sigma(std::max(s,SIG_TAU_EPS)), tau(std::max(t,SIG_TAU_EPS)),
      inv_sigma(1.0/sigma), ratio(sigma/tau), log_tau(std::log(tau)),
      lower(lo), upper(hi), log_normalizer(0.0),
      approximation(shifted_exp_sigma>0 && s<shifted_exp_sigma ? ShiftedExponential :
                    normal_ratio>0 && t<normal_ratio*s ? Normal : Exact) {
    // Approximation thresholds apply to the original, unclamped scales.
    if (approximation!=Exact) {
      sigma=s; tau=t; inv_sigma=1.0/s; log_tau=std::log(t);
    }
    if (!is_neg_inf(lo) || !is_pos_inf(hi))
      log_normalizer=raw_log_interval(lo,hi);
  }

  // exp(a*a/2-a*z)*Phi(z-a), sharing the Gaussian with the density.
  double weighted_tail(double z, bool log_p) const {
    const double w=ratio-z;
    if (w >= 0.0) {
      const double h=pnorm_detail::normal_tail_factor(w);
      return log_p ? (-0.5*z)*z+std::log(h)
                   : pnorm_detail::gaussian_exponential(z)*h;
    }
    const double exponent=ratio*(0.5*ratio-z);
    return log_p ? exponent+PNORM_STD(-w,true,true)
                 : std::exp(exponent)*PNORM_STD(-w,true,false);
  }

  double raw_pdf(double x, bool log_p) const {
    if (is_inf(x)) return log_p ? neg_inf() : 0.0;
    const double z=(x-mu)*inv_sigma;
    if (approximation==Normal) {
      const double lp=(-0.5*z)*z-pnorm_detail::LOG_SQRT_2PI-std::log(sigma);
      return log_p ? lp : std::exp(lp);
    }
    if (approximation==ShiftedExponential) {
      if (x<mu) return log_p ? neg_inf() : 0.0;
      const double lp=-(x-mu)/tau-log_tau;
      return log_p ? lp : std::exp(lp);
    }
    if (log_p) return weighted_tail(z,true)-log_tau;
    // Divide before the final exponential can underflow when tau is tiny.
    if (std::fabs(z)>37.0) return std::exp(weighted_tail(z,true)-log_tau);
    return weighted_tail(z,false)/tau;
  }

  double raw_prob(double x, bool lower_tail, bool log_p) const {
    if (is_nan(x)) return x;
    if (is_inf(x)) {
      const bool one=(x>0)==lower_tail;
      return log_p ? (one ? 0.0 : neg_inf()) : double(one);
    }
    const double z=(x-mu)*inv_sigma;
    if (is_nan(z)) return z;
    if (approximation==Normal) return PNORM_STD(z,lower_tail,log_p);
    if (approximation==ShiftedExponential) {
      const double ls=-std::max(0.0,x-mu)/tau;
      return lower_tail ? (log_p ? log1m_exp(ls) : -std::expm1(ls))
                        : (log_p ? ls : std::exp(ls));
    }
    if (!lower_tail) {
      if (log_p) {
        const double lp=log_sum_exp(PNORM_STD(-z,true,true),weighted_tail(z,true));
        return lp > -composite_detail::LOG2 ? log1m_exp(raw_prob(x,true,true)) : lp;
      }
      return std::min(1.0,PNORM_STD(-z,true,false)+weighted_tail(z,false));
    }
    if (z<=0.0) {
      const double h=pnorm_detail::normal_tail_difference(-z,ratio);
      return log_p ? (-0.5*z)*z+std::log(h)
                   : pnorm_detail::gaussian_exponential(z)*h;
    }
    const double phi=PNORM_STD(z,true,false);
    double cdf=phi-weighted_tail(z,false);
    if (cdf < 0.05*phi) {
      if (ratio<z) {
        const double exponent=ratio*(0.5*ratio-z);
        cdf=pnorm_detail::normal_interval_width(z-ratio,ratio)
            -std::expm1(exponent)*PNORM_STD(z-ratio,true,false);
      } else {
        cdf=pnorm_detail::gaussian_exponential(z)*
            pnorm_detail::normal_tail_difference(-z,ratio);
      }
    }
    if (log_p && cdf>0.5)
      return log1m_exp(raw_prob(x,false,true));
    return log_p ? std::log(cdf) : cdf;
  }

  double raw_log_interval(double lo,double hi) const {
    if (is_nan(lo) || is_nan(hi) || lo>hi) return na_real();
    if (lo==hi) return neg_inf();
    if (approximation==Normal)
      return pnorm_detail::normal_interval((lo-mu)*inv_sigma,(hi-mu)*inv_sigma,true);
    if (approximation==ShiftedExponential) {
      lo=std::max(lo,mu);
      if (hi<=lo) return neg_inf();
      return -(lo-mu)/tau+log1m_exp(-(hi-lo)/tau);
    }
    if (is_neg_inf(lo)) return raw_prob(hi,true,true);
    if (is_pos_inf(hi)) return raw_prob(lo,false,true);
    const double width=hi-lo, centre=lo+0.5*width;
    const double z=(centre-mu)*inv_sigma;
    if (width*inv_sigma*(1.0+std::fabs(z))<0.1) {
      const double ld=raw_pdf(centre,true);
      const double mass=pnorm_detail::integrate_unit([&](double u) {
        return std::exp(raw_pdf(lo+width*u,true)-ld);
      });
      return std::log(width)+ld+std::log(mass);
    }
    const bool use_lower=raw_prob(lo,true,false)<0.5;
    return log_diff_exp(raw_prob(use_lower ? hi : lo,use_lower,true),
                        raw_prob(use_lower ? lo : hi,use_lower,true));
  }

  double pdf(double x,bool log_p) const {
    if (x<=lower || x>=upper) return log_p ? neg_inf() : 0.0;
    if (is_neg_inf(lower) && is_pos_inf(upper)) return raw_pdf(x,log_p);
    const double lp=raw_pdf(x,true)-log_normalizer;
    return log_p ? lp : std::exp(lp);
  }

  double prob(double x,bool lower_tail,bool log_p) const {
    if (is_nan(x)) return x;
    if (is_neg_inf(lower) && is_pos_inf(upper)) return raw_prob(x,lower_tail,log_p);
    if (x<=lower || x>=upper) {
      const bool one=(x>=upper)==lower_tail;
      return log_p ? (one ? 0.0 : neg_inf()) : double(one);
    }
    const double lp=std::min((lower_tail ? raw_log_interval(lower,x)
                                      : raw_log_interval(x,upper))-log_normalizer,0.0);
    return log_p ? lp : std::exp(lp);
  }
};

inline double dexg(double x,double mu=5.,double sigma=1.,double tau=1.,bool log_d=false,
                   double normal_ratio=0.,double shifted_exp_sigma=0.) {
  if (sigma<=0. || tau<=0.) return na_real();
  return ExGaussian(mu,sigma,tau,neg_inf(),pos_inf(),normal_ratio,shifted_exp_sigma).raw_pdf(x,log_d);
}
inline double pexg(double q,double mu=5.,double sigma=1.,double tau=1.,
                   bool lower_tail=true,bool log_p=false,
                   double normal_ratio=0.,double shifted_exp_sigma=0.) {
  if (sigma<=0. || tau<=0.) return na_real();
  return ExGaussian(mu,sigma,tau,neg_inf(),pos_inf(),normal_ratio,shifted_exp_sigma).raw_prob(q,lower_tail,log_p);
}
inline double exg_log_interval(double lower,double upper,double mu,double sigma,double tau,
                               double normal_ratio=0.,double shifted_exp_sigma=0.) {
  if (sigma<=0. || tau<=0.) return na_real();
  return ExGaussian(mu,sigma,tau,neg_inf(),pos_inf(),normal_ratio,shifted_exp_sigma).raw_log_interval(lower,upper);
}
inline double dtexg(double x,double mu=5.,double sigma=1.,double tau=1.,
                    double lower=neg_inf(),double upper=pos_inf(),bool log_d=false,
                    double normal_ratio=0.,double shifted_exp_sigma=0.) {
  if (sigma<=0. || tau<=0. || lower>=upper) return na_real();
  return ExGaussian(mu,sigma,tau,lower,upper,normal_ratio,shifted_exp_sigma).pdf(x,log_d);
}
inline double ptexg(double q,double mu=5.,double sigma=1.,double tau=1.,
                    double lower=neg_inf(),double upper=pos_inf(),
                    bool lower_tail=true,bool log_p=false,
                    double normal_ratio=0.,double shifted_exp_sigma=0.) {
  if (sigma<=0. || tau<=0. || lower>=upper) return na_real();
  return ExGaussian(mu,sigma,tau,lower,upper,normal_ratio,shifted_exp_sigma).prob(q,lower_tail,log_p);
}
#endif
