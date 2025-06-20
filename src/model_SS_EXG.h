#ifndef ss_exg_h
#define ss_exg_h

#include <cmath>
#include <Rcpp.h>
#include "utility_functions.h"
#include "exgaussian_functions.h"
#include "race_integrate.h"
using namespace Rcpp;

// wrapper around ex-Gaussian log PDF for race function
NumericVector exg_go_lpdf(
    // single RT, broadcast to repped vector (one per accumulator)
    NumericVector rt,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // accumulator index: for which accumulator(s) should log density be computed?
    LogicalVector idx,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {

  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;

  NumericVector out(n_acc_selected);
  int k = 0;

  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;

    // input args: x, mu, sigma, tau, log_d = TRUE
    double log_d = dexg(rt[i], pars(i, 0), pars(i, 1), pars(i, 2), true);
    out[k] = std::isfinite(log_d) ? log_d : min_ll;

    k++;
  }

  return(out);
}

// wrapper around ex-Gaussian log complementary CDF for race function
NumericVector exg_go_lccdf(
    // single RT, broadcast to repped vector (one per accumulator)
    NumericVector rt,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // accumulator index: for which accumulator(s) should survivor probability be computed?
    LogicalVector idx,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {

  const int n_acc = rt.size();
  const int n_acc_selected = sum(idx);
  if (n_acc_selected == 0) return NA_REAL;

  NumericVector out(n_acc_selected);
  int k = 0;

  for (int i = 0; i < n_acc; i++) {
    if (!idx[i]) continue;

    // input args: q, mu, sigma, tau, lower_tail = FALSE, log_p = TRUE
    double log_s = pexg(rt[i], pars(i, 0), pars(i, 1), pars(i, 2), false, true);
    out[k] = std::isfinite(log_s) ? log_s : min_ll;

    k++;
  }

  return(out);
}

// go race log likelihood function, not accounting for go failure
double ss_exg_go_lpdf(
    // single RT
    double RT,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // index of accumulator winner and loser(s)
    LogicalVector winner,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {
  // set up an instance of a race function (class definition in `race_integrate.h`)
  race_f lpdf(exg_go_lpdf, exg_go_lccdf, pars, winner, min_ll, true);
  // evaluate the race function for the observed RT
  return lpdf(RT);
}

// go vs stop race log likelihood function for the case of stop trials with
// failed inhibition (i.e., stop process losing), not accounting for trigger
// failure and go failure.
double ss_exg_stop_fail_lpdf(
    // single RT
    double RT,
    // single stop signal delay
    double SSD,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // index of accumulator winner and loser(s)
    LogicalVector winner,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {
  // obtain the go race process log likelihood
  double go_lprob = ss_exg_go_lpdf(RT, pars, winner, min_ll);
  // obtain the survivor log probability of the stop process
  // NB SSD subtracted from observed RT to get stop finish time
  // input args: q, muS, sigmaS, tauS, lower_tail = FALSE, log_p = TRUE
  double stop_survivor_lprob = pexg(
    RT - SSD, pars(0, 3), pars(0, 4), pars(0, 5), false, true
  );
  // final output of race model is summed log likelihood
  return go_lprob + stop_survivor_lprob;
}

// go vs stop race likelihood function, for the case of stop trials with
// successful inhibition (i.e., stop process winning). not accounting for
// trigger failure and go failure.
// this function is an integrand (quantity to be integrated)
double exg_stop_success_integrand(
    // stop finish time
    double x,
    // stop signal delay
    double SSD,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // index of accumulator winner and loser(s)
    LogicalVector winner,
    // minimal log likelihood, to protect against numerical issues
    double min_ll
) {
  // obtain the go race process log likelihood (since the go accumulators are
  // losers, this will just be the summed log survivor probability)
  // NB SSD added to stop finish time to get go RT
  double go_lprob = ss_exg_go_lpdf(x + SSD, pars, winner, min_ll);
  // obtain the winner log probability of the stop process
  // input args: x, muS, sigmaS, tauS, log_d = TRUE
  double stop_winner_lprob = dexg(x, pars(0, 3), pars(0, 4), pars(0, 5), true);
  // final output of race model is summed log likelihood, exponentiated to the
  // likelihood to enable numerical integration
  return std::exp(go_lprob + stop_winner_lprob);
}

// wrapper class to make exg_stop_success_integrand compatible with Func
// interface for integration
class exg_ss_integrand : public Func {
private:
  const double SSD;
  const NumericMatrix pars;
  const LogicalVector winner;
  const double min_ll;

public:
  exg_ss_integrand(
    double SSD_,
    NumericMatrix pars_,
    LogicalVector winner_,
    double min_ll_
  ) :
  SSD(SSD_),
  pars(pars_),
  winner(winner_),
  min_ll(min_ll_) {}

  double operator()(const double& x) const {
    return exg_stop_success_integrand(x, SSD, pars, winner, min_ll);
  }
};

// function to compute the stop success integral, not accounting for trigger
// failure and go failure
double ss_exg_stop_success_lpdf(
    // single stop signal delay
    double SSD,
    // parameter values: rows = accumulators, columns = parameters
    NumericMatrix pars,
    // index of accumulator winner and loser(s)
    LogicalVector winner,
    // minimal log likelihood, to protect against numerical issues
    double min_ll,
    // lower limit for integration
    double lower = R_NegInf,
    // upper limit for integration
    double upper = R_PosInf
) {
  // set up an instance of a stop success integrand
  exg_ss_integrand race_integrand(SSD, pars, winner, min_ll);
  // perform integration: likelihood of stop process winning
  IntegrationResult out = my_integrate(race_integrand, lower, upper);
  // check for numerical issues
  bool bad_out = out.error_code != 0 || !traits::is_finite<REALSXP>(out.value);
  // return *log* likelihood
  // NB in the original R code from the DMC toolbox (`my.integrate`), -Inf was
  // returned in case of failed integration
  return bad_out ? min_ll : std::log(out.value);
}

// top-level log-likelihood function for the stop signal task
NumericVector ss_exg_lpdf(
    NumericVector RT,
    IntegerVector R,
    NumericVector SSD,
    NumericVector lR,
    LogicalVector winner,
    NumericMatrix pars,
    LogicalVector is_ok,
    double min_ll
) {
  // pars columns: mu=0, sigma=1, tau=2, muS=3, sigmaS=4, tauS=5, tf=6, gf=7

  NumericVector unique_lR = unique(lR);
  const int n_acc = unique_lR.length();         // number of go accumulators
  const int n_trials = lR.length() / n_acc;
  NumericVector out(n_trials);                  // trial-wise log-likelihood

  // extract go failure and trigger failure parameters as trial-wise vectors
  // from pars matrix; these parameters cannot differ between accumulators
  // within a given trial
  LogicalVector is_first_acc = lR == unique_lR[0]; // arbitrary index for one accumulator
  NumericVector tf = pars(_, 6);
  NumericVector gf = pars(_, 7);
  tf = tf[is_first_acc];
  gf = gf[is_first_acc];

  // local variables for intermediate computation; to be overwritten trial by trial
  double go_lprob;
  double stop_fail_lprob;
  double stop_success_integral;
  double stop_success_lprob;

  // loop over trials
  for (int trial = 0; trial < n_trials; trial++) {

    // get indices for all accumulators corresponding to the current trial
    int start_row = trial * n_acc;
    int end_row = (trial + 1) * n_acc - 1;

    // if any accumulator was flagged as not ok for the current trial, assign
    // minimal log likelihood and skip to next trial
    if (is_true(any(!is_ok[Range(start_row, end_row)]))) {
      out[trial] = min_ll;
      continue;
    }

    // branching logic:
    // first, was a response observed? second, was a stop signal presented?
    bool response_observed = !NumericVector::is_na(R[start_row]);
    bool stop_signal_presented = std::isfinite(SSD[start_row]);

    if (response_observed) {
      if (stop_signal_presented) {
        // stop trial with a response:
        // explained as a probabilistic mixture of
        // (i)  go process "won" because stop process wasn't triggered
        //      (trigger failure); OR
        // (ii) go process beat stop process in fair race (i.e., without
        //      trigger failure)
        // Lastly, this mixture has to be adjusted for lack of go failure
        go_lprob = ss_exg_go_lpdf(
          RT[start_row],
          pars(Range(start_row, end_row), _),
          winner[Range(start_row, end_row)],
          min_ll
        );
        stop_fail_lprob = ss_exg_stop_fail_lpdf(
          RT[start_row],
          SSD[start_row],
          pars(Range(start_row, end_row), _),
          winner[Range(start_row, end_row)],
          min_ll
        );
        // likelihood = (1-gf) x [tf x go_prob + (1-tf) x stop_fail_prob]
        out[trial] = log1m(gf[trial]) + log_mix(tf[trial], go_lprob, stop_fail_lprob);
      } else {
        // go trial with a response:
        // explained by go process without go failure
        go_lprob = ss_exg_go_lpdf(
          RT[start_row],
          pars(Range(start_row, end_row), _),
          winner[Range(start_row, end_row)],
          min_ll
        );
        // likelihood = (1-gf) x go_prob
        out[trial] = log1m(gf[trial]) + go_lprob;
      }
    } else {
      if (stop_signal_presented) {
        // stop trial with no response:
        // explained as a probabilistic mixture of
        // (i)  go failure; OR
        // (ii) stop process beat go process in fair race (i.e., without go
        //      failure and without trigger failure)
        stop_success_integral = ss_exg_stop_success_lpdf(
          SSD[start_row],
          pars(Range(start_row, end_row), _),
          winner[Range(start_row, end_row)],
          min_ll
        );
        stop_success_lprob = log1m(gf[trial]) + log1m(tf[trial]) + stop_success_integral;
        // likelihood = gf + [(1-gf) x (1-tf) x stop_success_integral]
        out[trial] = log_sum_exp(std::log(gf[trial]), stop_success_lprob);
      } else {
        // go trial with no response:
        // explained by go failure
        // likelihood = gf
        out[trial] = std::log(gf[trial]);
      }
    }
  }

  return(out);
}


// FOLLOWING ONLY KEPT IN FOR POSTERITY / TESTING R CODE
NumericVector dexg_c(
    const NumericVector x,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool log_d = false
) {
  int n = x.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = dexg(x[i], mu, sigma, tau, log_d);
  }
  return(out);
}

NumericVector pexg_c(
    const NumericVector q,
    const double mu = 5.,
    const double sigma = 1.,
    const double tau = 1.,
    const bool lower_tail = true,
    const bool log_p = false
) {
  int n = q.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = pexg(q[i], mu, sigma, tau, lower_tail, log_p);
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector dEXGrace(
    NumericMatrix dt, NumericVector mu, NumericVector sigma, NumericVector tau
){
  int n = mu.size();
  NumericVector log_out(dt.nrow());
  log_out = dexg_c(dt(0, _), mu[0], sigma[0], tau[0], true);
  for (int i = 1; i < n; i++){
    log_out += pexg_c(dt(i, _), mu[i], sigma[i], tau[i], false, true);
  }
  return exp(log_out);
}

// [[Rcpp::export]]
NumericVector stopfn_exg(
    NumericVector t, NumericVector mu, NumericVector sigma, NumericVector tau, double SSD
){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dEXGrace(dt, mu, sigma, tau);
}

#endif
