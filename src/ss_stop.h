#ifndef SS_STOP_H
#define SS_STOP_H

#include <cmath>
#include <algorithm>
#include "pnorm_utils.h"

enum class SSStopFamily { ExGaussian, Normal, Lognormal, Weibull };

// Prepared non-exGaussian stop distributions. Integration uses standardized
// normal/lognormal time or Weibull cumulative hazard, including the Jacobian.
struct SSStop {
  SSStopFamily family;
  double location, scale, lower, inv_scale, log_scale, alpha, log_tail_factor;

  SSStop(SSStopFamily f, double a, double b, double lb=0.0)
    : family(f), location(a), scale(b), lower(lb), inv_scale(1.0/b),
      log_scale(std::log(b)), alpha((lb-a)/b), log_tail_factor(0.0) {
    if (family==SSStopFamily::Normal)
      log_tail_factor=alpha>0 ? std::log(pnorm_detail::normal_tail_factor(alpha))
                             : PNORM_STD(alpha,false,true);
  }

  double log_survivor(double t) const {
    if (family==SSStopFamily::Normal) {
      if (t<=lower) return 0.0;
      if (is_pos_inf(t)) return neg_inf();
      if (alpha>0) {
        const double d=(t-lower)*inv_scale;
        return -d*(alpha+0.5*d) +
          std::log(pnorm_detail::normal_tail_factor(alpha+d))-log_tail_factor;
      }
      return PNORM_STD((t-location)*inv_scale,false,true)-log_tail_factor;
    }
    if (t<=0) return 0.0;
    if (family==SSStopFamily::Lognormal)
      return PNORM_STD((std::log(t)-location)*inv_scale,false,true);
    return -std::exp(location*(std::log(t)-log_scale));
  }

  double time(double x) const {
    if (family==SSStopFamily::Weibull)
      return x<=0 ? 0.0 : std::exp(log_scale+std::log(x)/location);
    if (family==SSStopFamily::Lognormal) return std::exp(location+scale*x);
    return alpha>0 ? lower+scale*(x/(1.0+alpha)) : location+scale*x;
  }

  double coordinate(double t) const {
    if (family==SSStopFamily::Weibull)
      return t<=0 ? 0.0 : std::exp(location*(std::log(t)-log_scale));
    if (family==SSStopFamily::Lognormal)
      return t<=0 ? neg_inf() : (std::log(t)-location)*inv_scale;
    return alpha>0 ? ((t-lower)*inv_scale)*(1.0+alpha) : (t-location)*inv_scale;
  }

  double log_weight(double x) const {
    if (family==SSStopFamily::Weibull) return -x;
    if (family==SSStopFamily::Normal && alpha>0) {
      const double d=x/(1.0+alpha);
      return -d*(alpha+0.5*d)-pnorm_detail::LOG_SQRT_2PI-
        log_tail_factor-std::log1p(alpha);
    }
    return -0.5*x*x-pnorm_detail::LOG_SQRT_2PI-
      (family==SSStopFamily::Normal ? log_tail_factor : 0.0);
  }

  void bounds(double upper, double& lo, double& hi) const {
    // Omitted stop mass is below 3e-20 (10 normal SDs or exp(-45)).
    // Positive normal truncation is scaled by its tail width, sigma/(1+alpha).
    lo=family==SSStopFamily::Weibull ? 0.0 :
       family==SSStopFamily::Normal ? (alpha>0 ? 0.0 : std::max(alpha,-10.0)) : -10.0;
    hi=family==SSStopFamily::Weibull ? 45.0 :
       family==SSStopFamily::Normal && alpha>0 ? 48.0 : 10.0;
    hi=std::min(hi,coordinate(upper));
  }
};

#endif
