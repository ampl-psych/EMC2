#ifndef wald_functions_h
#define wald_functions_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "composite_functions.h"
using namespace Rcpp;

const double L_PI = 1.1447298858494001741434;  // std::log(M_PI)

// helper functions to check for bad parameters / edge cases of Wald distribution
static inline bool wald_bad_params(
    const double x,
    const double mu,
    const double lambda
) {
  return(
    std::isnan(x) || std::isnan(mu) || std::isnan(lambda) ||
    (mu < 0.0) || (lambda < 0.0)
  );
}

static inline bool wald_lower_limit(
    const double x,
    const double mu,
    const double lambda
) {
  return(
    (x < 0.0) || (x == 0.0 && mu > 0.0 && lambda > 0.0) ||
    (x < mu && std::isinf(lambda))
  );
}

static inline bool wald_upper_limit(
    const double x,
    const double mu,
    const double lambda
) {
  return(
    (std::isinf(x)) || (x > 0.0 && lambda == 0.0) ||
    (x > mu && (mu == 0.0 || std::isinf(lambda)))
  );
}

static inline bool wald_spike(
    const double x,
    const double mu,
    const double lambda
) {
  return(
    (x == 0.0 && lambda == 0.0) || (x == mu && (mu == 0.0 || std::isinf(lambda)))
  );
}

// PDF and (C)CDF of Wald distribution, parameterised in terms of drift
// coefficient v and threshold b. Assumes v and b have already been scaled by
// diffusive noise s. Parameters v and b are then used to determine the mean
// (mu = b / v) and shape (lambda = b^2) of the Wald distribution.

// Wald PDF
static double dwald(
    const double x,
    const double v = 1.0,
    const double b = 1.0,
    const bool log_d = false
) {
  double mu = b / v;
  double lambda = b * b;

  if (wald_bad_params(x, mu, lambda)) {
    return NA_REAL;
  }
  if (wald_lower_limit(x, mu, lambda) || wald_upper_limit(x, mu, lambda)) {
    return R_NegInf;
  }
  if (wald_spike(x, mu, lambda)) {
    return R_PosInf;
  }

  double lprob;
  if (std::isinf(mu)) {
    // accurate approximation of Wald LPDF when mu == Inf:
    // Levy distribution with parameters location = 0 and scale = lambda
    lprob = -0.5 * (
      (lambda / x) - std::log(lambda) + std::log(2.0 * M_PI) + 3.0 * std::log(x)
    );
  } else {
    // otherwise, density is given by Wald distribution function with parameters
    // mean = mu, shape = lambda
    double term1 = 0.5 * std::log(lambda / (2.0 * M_PI));
    double term2 = ((x - mu) * (x - mu)) / (x * mu * mu);
    lprob = term1 -1.5 * std::log(x) - 0.5 * lambda * term2;
  }

  return log_d ? lprob : std::exp(lprob);
}

// Wald (C)CDF
static double pwald(
    const double q,
    const double v = 1.0,
    const double b = 1.0,
    const bool lower_tail = true,
    const bool log_p = false
) {
  double mu = b / v;
  double lambda = b * b;
  double out;

  if (wald_bad_params(q, mu, lambda)) {
    return NA_REAL;
  }
  if (wald_lower_limit(q, mu, lambda)) {
    out = lower_tail? R_NegInf : 0.0;
    return log_p? out : std::exp(out);
  }
  if (wald_upper_limit(q, mu, lambda)) {
    out = lower_tail? 0.0 : R_NegInf;
    return log_p? out : std::exp(out);
  }
  if (wald_spike(q, mu, lambda)) {
    return log_p? 0.0 : 1.0;
  }

  if (std::isinf(mu)) {
    // accurate approximation of Wald L(C)CDF when mu == Inf:
    // Chi-Squared distribution with df = 1
    out = R::pchisq((lambda / q), 1.0, !lower_tail, true);
  } else if ((mu / lambda) < 1e-14) {
    // accurate approximation of Wald L(C)CDF when ratio mu / lambda is tiny:
    // Gamma distribution with shape = lambda / mu and scale (inverse rate) =
    // mu^2 / lambda
    out = R::pgamma(q, (lambda / mu), ((mu * mu) / lambda), lower_tail, true);
  } else {
    // otherwise, Wald distribution with parameters mean = mu, shape = lambda
    double q_mu = q / mu;
    double phi_mu = mu / lambda;
    double r = std::sqrt(q_mu * phi_mu);
    double a = R::pnorm(((q_mu - 1.0) / r), 0.0, 1.0, lower_tail, true);
    double b = (2.0 / phi_mu) +
      R::pnorm((-(q_mu + 1.0) / r), 0.0, 1.0, true, true);
    if (lower_tail) {
      out = a + log1p_exp(b - a);
    } else {
      double q_phi_mu = q_mu / (2.0 * phi_mu);
      // if input value q is extremely large, use asymptotic expression given by
      // Giner & Smyth (2016)
      if (q_mu > 1e6 || q_phi_mu > 5e5) {
        out = (1 / phi_mu) - 0.5 * L_PI - std::log(2 * phi_mu) - q_phi_mu -
          1.5 * std::log1p(q_phi_mu);
      } else {
        out = a + log1m_exp(b - a);
      }
    }
  }

  return log_p? out : std::exp(out);
}

// original Wald CDF, parameterised in terms of drift coefficient l and threshold k
double pigt0(double t, double k = 1., double l = 1.){
  //if (t <= 0.){
  //  return 0.;
  //}
  double mu = k / l;
  double lambda = k * k;

  double p1 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. + t/mu), 0., 1., true, false);
  double p2 = 1 - R::pnorm(std::sqrt(lambda/t) * (1. - t/mu), 0., 1., true, false);

  return std::exp(std::exp(std::log(2. * lambda) - std::log(mu)) + std::log(p1)) + p2;
}

// original Wald PDF, parameterised in terms of drift coefficient l and threshold k
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

// CDF for single-boundary diffusion process with random across-trial variability
// in threshold, based on Logan et al. (2014). Parameterised in terms of drift
// coefficient l, threshold k, and uniform variability in threshold +- a.
double pigt(
    double t,
    double k = 1,
    double l = 1,
    double a = .1,
    double threshold = 1e-10,
    bool lower_tail = true,
    bool log_p = false
){
  if (t <= 0.){
    return 0.;
  }
  double cdf;

  if (a < threshold){
    // return pwald(t, l, k, lower_tail, log_p);
    cdf = pigt0(t, k, l);
    if (!lower_tail) {
      return log_p? log1m(cdf) : (1.0 - cdf);
    }
    return log_p? std::log(cdf) : cdf;
  }

  double sqt = std::sqrt(t);
  double lgt = std::log(t);

  if (l < threshold){
    double t5a = 2. * R::pnorm((k + a) / sqt, 0., 1., true, false) - 1;
    double t5b = 2. * R::pnorm((- k - a) / sqt, 0., 1., true, false) - 1;

    double t6a = - .5 * ((k + a) * (k + a) / t - M_LN2 - L_PI + lgt) - std::log(a);
    double t6b = - .5 * ((k - a) * (k - a) / t - M_LN2 - L_PI + lgt) - std::log(a);

    cdf = 1. + std::exp(t6a) - std::exp(t6b) + ((- k + a) * t5a - (k - a) * t5b) / (2. * a);
  } else {
    double t1a = std::exp(- .5 * std::pow(k - a - t * l, 2) / t);
    double t1b = std::exp(- .5 * std::pow(a + k - t * l, 2) / t);
    double t1 = std::exp(.5* (lgt - M_LN2 - L_PI)) * (t1a - t1b);

    double t2a = std::exp(2. * l * (k - a) + R::pnorm(- (k - a + t * l) / sqt, 0., 1., true, true));
    double t2b = std::exp(2. * l * (k + a) + R::pnorm(- (k + a + t * l) / sqt, 0., 1., true, true));
    double t2 = a + (t2b - t2a) / (2. * l);

    double t4a = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
    double t4b = 2. * R::pnorm((k - a) / sqt - sqt * l, 0., 1., true, false) - 1.;
    double t4 = .5 * (t * l - a - k + .5 / l) * t4a + .5 * (k - a - t * l - .5 / l) * t4b;

    cdf = .5 * (t4 + t2 + t1) / a;
  }

  if (cdf < 0. || std::isnan(cdf)) {
    cdf = 0.;
  }

  if (!lower_tail) {
    return log_p? log1m(cdf) : (1.0 - cdf);
  }
  return log_p? std::log(cdf) : cdf;
}

// PDF for single-boundary diffusion process with random across-trial variability
// in threshold, based on Logan et al. (2014). Parameterised in terms of drift
// coefficient l, threshold k, and uniform variability in threshold +- a.
double digt(
    double t,
    double k = 1.,
    double l = 1.,
    double a = .1,
    double threshold= 1e-10,
    bool log_d = false
){
  if (t <= 0.){
    return 0.;
  }
  double pdf;

  if (a < threshold){
    // return dwald(t, l, k, log_d);
    pdf = digt0(t, k, l);
    return log_d? std::log(pdf) : pdf;
  }

  if (l < threshold){
    double term = std::exp(- (k - a) * (k - a) / (2. * t)) - std::exp(- (k + a) * (k + a) / (2. * t));
    pdf = std::exp(-.5 * (M_LN2 + L_PI + std::log(t)) + std::log(term) - M_LN2 - std::log(a));
  } else {
    double sqt = std::sqrt(t);

    double t1a = - std::pow(a - k + t * l, 2) / (2. * t);
    double t1b = - std::pow(a + k - t * l, 2) / (2. * t);
    double t1 = M_SQRT1_2 * (std::exp(t1a) - std::exp(t1b)) / (std::sqrt(M_PI) * sqt);

    double t2a = 2. * R::pnorm((- k + a) / sqt + sqt * l, 0., 1., true, false) - 1.;
    double t2b = 2. * R::pnorm((k + a) / sqt - sqt * l, 0., 1., true, false) - 1.;
    double t2 = std::exp(std::log(.5) + std::log(l)) * (t2a + t2b);

    pdf = std::exp(std::log(t1 + t2) - M_LN2 - std::log(a));
  }

  if (pdf < 0. || std::isnan(pdf)) {
    pdf = 0.;
  }
  return log_d? std::log(pdf) : pdf;
}

#endif
