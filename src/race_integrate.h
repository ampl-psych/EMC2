#ifndef RACE_INTEGRATE_H
#define RACE_INTEGRATE_H

#include <Rcpp.h>
#include <cmath>
#include "gauss.h"  // hcubature + Gauss-Kronrod prototypes
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

using namespace Rcpp;

// Lightweight functor interface for 1D integrands
class Func {
public:
  virtual ~Func() {}
  virtual double operator()(const double& x) const = 0;
};

// Wrapper to use hcubature/Gauss-Kronrod for 1D integration
// Supports finite and infinite bounds via change-of-variable
namespace detail {

struct IntegrandCtx {
  const Func* f;
};

inline int integrand_1d(unsigned /*dim*/, const double* x, void* p,
                        unsigned /*fdim*/, double* retval) {
  IntegrandCtx* ctx = static_cast<IntegrandCtx*>(p);
  const Func* f = ctx->f;
  double val = 0.0;
  if (f) {
    val = (*f)(x[0]);
  }
  if (!R_FINITE(val)) {
    // Return 0 to avoid poisoning the integrator
    val = 0.0;
  }
  retval[0] = val;
  return 0;
}

inline double integrate_finite(const Func& f, double a, double b,
                               double& err_est, int& err_code,
                               size_t maxEval, double abs_tol, double rel_tol) {
  IntegrandCtx ctx{&f};
  double val = 0.0;
  double err = 0.0;
  const unsigned n = 1;
  double lo = a, hi = b;
  int ret = hcubature(integrand_1d, &ctx, n, &lo, &hi,
                      maxEval, abs_tol, rel_tol, &val, &err);
  err_est = err;
  err_code = ret;
  return val;
}

// Change-of-variable helpers for infinite bounds
class UpperInfiniteTransform : public Func {
  const Func& f;
  const double a; // lower finite bound
public:
  UpperInfiniteTransform(const Func& f_, double a_) : f(f_), a(a_) {}
  double operator()(const double& t) const override {
    // Map t in (0,1) to x in (a, +inf): x = a + t/(1-t)
    // Jacobian dx/dt = 1/(1-t)^2
    if (t <= 0.0 || t >= 1.0) return 0.0;
    const double one_mt = 1.0 - t;
    const double x = a + t / one_mt;
    const double jac = 1.0 / (one_mt * one_mt);
    double val = f(x) * jac;
    return R_FINITE(val) ? val : 0.0;
  }
};

class LowerInfiniteTransform : public Func {
  const Func& f;
  const double b; // upper finite bound
public:
  LowerInfiniteTransform(const Func& f_, double b_) : f(f_), b(b_) {}
  double operator()(const double& t) const override {
    // Map t in (0,1) to x in (-inf, b): x = b - t/(1-t)
    // Jacobian dx/dt = 1/(1-t)^2
    if (t <= 0.0 || t >= 1.0) return 0.0;
    const double one_mt = 1.0 - t;
    const double x = b - t / one_mt;
    const double jac = 1.0 / (one_mt * one_mt);
    double val = f(x) * jac;
    return R_FINITE(val) ? val : 0.0;
  }
};

class BothInfiniteTransform : public Func {
  const Func& f;
public:
  explicit BothInfiniteTransform(const Func& f_) : f(f_) {}
  double operator()(const double& u) const override {
    // Map u in (0,1) to x in (-inf, +inf): x = tan(pi*(u-0.5))
    // Jacobian dx/du = pi * sec(pi*(u-0.5))^2
    if (u <= 0.0 || u >= 1.0) return 0.0;
    const double shift = M_PI * (u - 0.5);
    const double x = std::tan(shift);
    const double sec2 = 1.0 / (std::cos(shift) * std::cos(shift));
    const double jac = M_PI * sec2;
    double val = f(x) * jac;
    return R_FINITE(val) ? val : 0.0;
  }
};

} // namespace detail

// Public integrate API to match usage in SS model headers
inline double integrate(const Func& f,
                        double a, double b,
                        double& err_est, int& err_code,
                        int max_subdivisions,
                        double abs_tol, double rel_tol) {
  // Map max_subdivisions to an evaluation budget; allow unlimited when 0
  const size_t maxEval = (max_subdivisions <= 0) ? 0 : static_cast<size_t>(max_subdivisions) * 64;

  // Finite bounds
  if (std::isfinite(a) && std::isfinite(b)) {
    return detail::integrate_finite(f, a, b, err_est, err_code, maxEval, abs_tol, rel_tol);
  }

  // One-sided infinite
  if (std::isfinite(a) && !std::isfinite(b)) {
    detail::UpperInfiniteTransform g(f, a);
    return detail::integrate_finite(g, 0.0, 1.0, err_est, err_code, maxEval, abs_tol, rel_tol);
  }
  if (!std::isfinite(a) && std::isfinite(b)) {
    detail::LowerInfiniteTransform g(f, b);
    return detail::integrate_finite(g, 0.0, 1.0, err_est, err_code, maxEval, abs_tol, rel_tol);
  }

  // Two-sided infinite
  detail::BothInfiniteTransform g(f);
  return detail::integrate_finite(g, 0.0, 1.0, err_est, err_code, maxEval, abs_tol, rel_tol);
}


// ----------------------------------------------------------------------------
// Race helper used by SS models
// ----------------------------------------------------------------------------

// Function pointer types for per-accumulator log-density/log-survivor
typedef NumericVector (*race_logpdf_fn)(NumericVector, NumericMatrix, Rcpp::LogicalVector, double);
typedef NumericVector (*race_logsurv_fn)(NumericVector, NumericMatrix, Rcpp::LogicalVector, double);

class race_f : public Func {
private:
  race_logpdf_fn logpdf_f;   // returns log f for selected indices
  race_logsurv_fn logsurv_f; // returns log S for selected indices
  NumericMatrix pars;
  LogicalVector winner;
  const double min_ll;

public:
  race_f(race_logpdf_fn lpdf,
         race_logsurv_fn lsurv,
         const NumericMatrix& pars_,
         const LogicalVector& winner_,
         double min_ll_,
         bool inputs_are_log_scale = true)
    : logpdf_f(lpdf), logsurv_f(lsurv), pars(pars_), winner(winner_),
      min_ll(min_ll_) {
    (void)inputs_are_log_scale; // parameter kept for API compatibility
  }

  double operator()(const double& t) const override {
    const int n_acc = pars.nrow();
    NumericVector rt(n_acc, t);

    // log density for winner(s)
    NumericVector lw = logpdf_f(rt, pars, winner, min_ll);
    double log_pdf = (lw.size() > 0) ? sum(lw) : R_NegInf; // should be one
    if (!R_FINITE(log_pdf)) log_pdf = min_ll;

    // log survivor for losers
    LogicalVector loser = !winner;
    NumericVector ls = logsurv_f(rt, pars, loser, min_ll);
    double log_surv_sum = 0.0;
    for (int i = 0; i < ls.size(); i++) {
      log_surv_sum += R_FINITE(ls[i]) ? ls[i] : min_ll;
    }

    double out = log_pdf + log_surv_sum;
    if (!R_FINITE(out)) out = min_ll;
    return out;
  }
};

#endif // RACE_INTEGRATE_H
