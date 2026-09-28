// The Wald (inverse Gaussian) finishing-time distribution of an unpulsed
// racing-diffusion accumulator without start-point variability: drift v,
// threshold b, within-trial noise s, start 0. Used for the rows of a hybrid
// neural likelihood that the network never sees (register_nn_model(hybrid = ),
// R/nn_register.R) and by CRDM() where amp is exactly 0.
//
// The same distribution as EMC2's digt0 / pigt0 (wald_functions.h, k = b / s,
// l = v / s), which it equals to rounding, written in log form: pigt0's
// exp(2 k l) overflows for k l > 354, which lies inside the training box of the
// flows this is paired with (v <= 8, b <= 3, s >= 0.25). Rmath's pnorm, as in
// the flow evaluators, so results compare with the R references; reentrant.

#ifndef EMC2_NLE_WALD_H
#define EMC2_NLE_WALD_H

#include <Rmath.h>
#include <cmath>

namespace nle {

// log density at decision time t > 0
inline double wald_log_pdf(double t, double v, double b, double s) {
  const double k = b / s, l = v / s, d = k - l * t;
  return std::log(k) - 0.5 * std::log(2.0 * M_PI * t * t * t) - 0.5 * d * d / t;
}

// CDF at decision time t > 0
inline double wald_cdf(double t, double v, double b, double s) {
  const double k = b / s, l = v / s, sq = std::sqrt(t);
  const double p = R::pnorm((l * t - k) / sq, 0.0, 1.0, 1, 0) +
    std::exp(2.0 * k * l + R::pnorm(-(l * t + k) / sq, 0.0, 1.0, 1, 1));
  return p > 1.0 ? 1.0 : p;
}

}  // namespace nle

#endif
