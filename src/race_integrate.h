#ifndef race_integrate_h
#define race_integrate_h

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>

using namespace Rcpp;
using namespace Numer;

struct IntegrationResult {
  double value;
  double error_estimate;
  int error_code;
};

class race_f : public Func {
private:
  NumericMatrix pars;
  LogicalVector winner;
  NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double);
  NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double);
  double min_ll;
  double min_p;

public:
  race_f(
    NumericMatrix pars_,
    LogicalVector winner_,
    NumericVector (*dfun_)(NumericVector, NumericMatrix, LogicalVector, double),
    NumericVector (*pfun_)(NumericVector, NumericMatrix, LogicalVector, double),
    double min_ll_
  ) :
  pars(pars_),
  winner(winner_),
  dfun(dfun_),
  pfun(pfun_),
  min_ll(min_ll_),
  min_p(std::exp(min_ll_)) {}

  // integrand definition
  double operator()(const double& x) const {
    int n_acc = pars.nrow();
    NumericVector t(n_acc, x);
    // compute probability of winner accumulator finishing at time t
    NumericVector d = dfun(t, pars, winner, min_p);
    if (d.size() != 1) {
      stop("dfun must return a scalar NumericVector of length 1");
    }
    double out = d[0];
    // multiply by probability of loser accumulator(s) not having finished by
    // time t
    if (n_acc > 1) {
      NumericVector p = 1. - pfun(t, pars, !winner, min_p);
      double prod_p = std::accumulate(p.begin(), p.end(), 1., std::multiplies<double>());
      out *= prod_p;
    }
    // output is instantaneous likelihood that winner accumulator finishes at
    // time t, before the loser accumulator(s).
    return out;
  }
};


IntegrationResult my_integrate(
    NumericMatrix pars,
    LogicalVector winner,
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
    double min_ll,
    double lower,
    double upper,
    bool check_err = true,
    double upper_retry = 10.,
    double fail_out = 0.
) {
  // set up the integrand (race model function)
  race_f f(pars, winner, dfun, pfun, min_ll);
  double err_est;
  int err_code;

  // initial attempt at numerical integration
  double res = integrate(f, lower, upper, err_est, err_code);

  if (!check_err) {
    return {res, err_est, err_code};
  }

  // if integration terminated abnormally due to max. number of subdivisions
  // reached, and upper bound was positive infinity, try again with large finite
  // upper bound
  if (err_code == 1 && upper == R_PosInf) {
    double err_est_retry;
    int err_code_retry;
    double res_retry = integrate(f, lower, upper_retry, err_est_retry, err_code_retry);
    // overwrite output elements
    res = res_retry;
    err_est = err_est_retry;
    err_code = err_code_retry;
  }

  // if integration (still) terminated abnormally, use fallback output
  if (err_code > 0) {
    res = fail_out;
  }

  return {res, err_est, err_code};
}


double pr_pt(
    NumericMatrix pars,
    LogicalVector winner,
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
    double min_ll,
    double LT,
    double UT
) {
  // probability in full response window (zero to positive infinity)
  IntegrationResult pr = my_integrate(pars, winner, dfun, pfun, min_ll, 0., R_PosInf);
  if (pr.error_code != 0 || traits::is_nan<REALSXP>(pr.value)) return NA_REAL;

  // probability in truncated response window [LT, UT]
  IntegrationResult pt = my_integrate(pars, winner, dfun, pfun, min_ll, LT, UT);
  if (pt.error_code != 0 || traits::is_nan<REALSXP>(pt.value)) return NA_REAL;

  // handle edge cases before computing and returning ratio of probabilities
  if (pt.value == 0.) return NA_REAL;
  if (pr.value == 0.) return 0.;

  double out = pr.value / pt.value;
  if (traits::is_infinite<REALSXP>(out)) return NA_REAL;

  return std::clamp(out, 0., 1.);
}


double pLU(
    NumericMatrix pars,
    LogicalVector winner,
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double),
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double),
    double min_ll,
    double LT,
    double LC,
    double UC,
    double UT
) {
  // probability in [LT, LC]
  IntegrationResult pL = my_integrate(pars, winner, dfun, pfun, min_ll, LT, LC);
  if (pL.error_code != 0 || traits::is_nan<REALSXP>(pL.value)) return NA_REAL;

  // probability in [UC, UT]
  IntegrationResult pU = my_integrate(pars, winner, dfun, pfun, min_ll, UC, UT);
  if (pU.error_code != 0 || traits::is_nan<REALSXP>(pU.value)) return NA_REAL;

  // clamp probabilities to [0, 1] and return summed probability
  double clampedL = std::clamp(pL.value, 0., 1.);
  double clampedU = std::clamp(pU.value, 0., 1.);

  return clampedL + clampedU;
}

#endif
