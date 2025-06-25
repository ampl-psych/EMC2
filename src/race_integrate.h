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
  NumericVector (*lpdf)(NumericVector, NumericMatrix, LogicalVector, double);
  NumericVector (*lccdf)(NumericVector, NumericMatrix, LogicalVector, double);
  const NumericMatrix pars;
  const LogicalVector winner;
  const double min_ll;
  const bool log_d;

public:
  race_f(
    NumericVector (*lpdf_)(NumericVector, NumericMatrix, LogicalVector, double),
    NumericVector (*lccdf_)(NumericVector, NumericMatrix, LogicalVector, double),
    NumericMatrix pars_,
    LogicalVector winner_,
    double min_ll_,
    bool log_d_ = false
  ) :
  lpdf(lpdf_),
  lccdf(lccdf_),
  pars(pars_),
  winner(winner_),
  min_ll(min_ll_),
  log_d(log_d_) {}

  double operator()(const double& x) const {
    int n_acc = pars.nrow();
    int n_winner = sum(winner);
    int n_loser = sum(!winner);
    // broadcast input response time x to repped vector, one per accumulator
    NumericVector t(n_acc, x);
    // initialise log likelihood
    double log_out = 0.;
    if (n_winner == 1) {
      // compute log density of winner accumulator finishing at time t
      NumericVector log_d = lpdf(t, pars, winner, min_ll);
      log_out += Rcpp::as<double>(log_d);
    }
    if (n_loser >= 1) {
      // (summed) log survival probability of losing accumulator(s)
      NumericVector log_s = lccdf(t, pars, !winner, min_ll);
      log_out += std::accumulate(log_s.begin(), log_s.end(), 0.);
    }
    // output is instantaneous (log) likelihood that winner accumulator finishes
    // at time t and the loser accumulator(s) have not yet finished by time t.
    return log_d ? log_out : std::exp(log_out);
  }
};


IntegrationResult my_integrate(
    const Func& race_pdf,
    double lower,
    double upper,
    bool check_err = true,
    double upper_retry = 10.,
    double fail_out = 0.,
    int max_subdivs = 100,
    double abs_tol = 1e-8,
    double rel_tol = 1e-6
) {

  // initial attempt at numerical integration
  double err_est;
  int err_code;
  double res = integrate(
    race_pdf, lower, upper,
    err_est, err_code,
    max_subdivs, abs_tol, rel_tol
  );

  if (!check_err) {
    return {res, err_est, err_code};
  }

  // integration could have failed for several reasons (7 different values for
  // err_code possible). If due to max. subdivisions reached (err_code == 1)
  // *AND* upper bound was Inf, try again with large finite upper bound.
  if (err_code == 1 && upper == R_PosInf) {
    double err_est_retry;
    int err_code_retry;
    double res_retry = integrate(
      race_pdf, lower, upper_retry,
      err_est_retry, err_code_retry,
      max_subdivs, abs_tol, rel_tol
    );
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
    const Func& race_pdf,
    double LT,
    double UT
) {
  // probability in full response window (zero to positive infinity)
  IntegrationResult pr = my_integrate(race_pdf, 0., R_PosInf);
  if (pr.error_code != 0 || traits::is_nan<REALSXP>(pr.value)) return NA_REAL;

  // probability in truncated response window [LT, UT]
  IntegrationResult pt = my_integrate(race_pdf, LT, UT);
  if (pt.error_code != 0 || traits::is_nan<REALSXP>(pt.value)) return NA_REAL;

  // handle edge cases before computing and returning ratio of probabilities
  if (pt.value == 0.) return NA_REAL;
  if (pr.value == 0.) return 0.;

  double out = pr.value / pt.value;
  if (traits::is_infinite<REALSXP>(out)) return NA_REAL;

  return std::max(0., std::min(1., out));
}


double pLU(
    const Func& race_pdf,
    double LT,
    double LC,
    double UC,
    double UT
) {
  // probability in [LT, LC]
  IntegrationResult pL = my_integrate(race_pdf, LT, LC);
  if (pL.error_code != 0 || traits::is_nan<REALSXP>(pL.value)) return NA_REAL;

  // probability in [UC, UT]
  IntegrationResult pU = my_integrate(race_pdf, UC, UT);
  if (pU.error_code != 0 || traits::is_nan<REALSXP>(pU.value)) return NA_REAL;

  // clamp probabilities to [0, 1] and return summed probability
  // out = std::max(0., std::min(1., out));
  double clampedL = std::max(0., std::min(1., pL.value));
  double clampedU = std::max(0., std::min(1., pU.value));

  return clampedL + clampedU;
}

#endif
