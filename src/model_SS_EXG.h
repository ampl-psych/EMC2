#ifndef ss_exg_h
#define ss_exg_h

#include <cmath>
#include <Rcpp.h>
#include "utility_functions.h"
#include "exgaussian_functions.h"
#include "race_integrate.h"
using namespace Rcpp;

// ----------------------------------------------------------------------------
// TRUNCATED EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

// wrapper around truncated ex-Gaussian log PDF for race function
NumericVector texg_go_lpdf(
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

    // input args: x, mu, sigma, tau, exg_lb, upper = Inf, log_d = TRUE
    double log_d = dtexg(
      rt[i], pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, true
    );
    out[k] = std::isfinite(log_d) ? log_d : min_ll;

    k++;
  }

  return(out);
}

// wrapper around truncated ex-Gaussian log complementary CDF for race function
NumericVector texg_go_lccdf(
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

    // input args: q, mu, sigma, tau, exg_lb, upper = Inf, lower_tail = FALSE, log_p = TRUE
    double log_s = ptexg(
      rt[i], pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, false, true
    );
    out[k] = std::isfinite(log_s) ? log_s : min_ll;

    k++;
  }

  return(out);
}

// go race log likelihood function, not accounting for go failure
double ss_texg_go_lpdf(
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
  race_f lpdf(texg_go_lpdf, texg_go_lccdf, pars, winner, min_ll, true);
  // evaluate the race function for the observed RT
  return lpdf(RT);
}

// go vs stop race log likelihood function for the case of stop trials with
// failed inhibition (i.e., stop process losing), not accounting for trigger
// failure and go failure.
double ss_texg_stop_fail_lpdf(
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
  double go_lprob = ss_texg_go_lpdf(RT, pars, winner, min_ll);
  // obtain the survivor log probability of the stop process
  // NB SSD subtracted from observed RT to get stop finish time
  // input args: q, muS, sigmaS, tauS, exgS_lb, upper = Inf, lower_tail = FALSE, log_p = TRUE
  double stop_survivor_lprob = ptexg(
    RT - SSD, pars(0, 3), pars(0, 4), pars(0, 5), pars(0, 9), R_PosInf, false, true
  );
  if (!traits::is_finite<REALSXP>(stop_survivor_lprob)) {
    stop_survivor_lprob = min_ll;
  }
  // final output of race model is summed log likelihood
  return go_lprob + stop_survivor_lprob;
}

// go vs stop race likelihood function, for the case of stop trials with
// successful inhibition (i.e., stop process winning). not accounting for
// trigger failure and go failure.
// this function is an integrand (quantity to be integrated)
class texg_stop_success_integrand : public Func {
private:
  const double SSD;
  const NumericMatrix pars;
  const int n_go;
  const double min_ll;

public:
  texg_stop_success_integrand(
    double SSD_,
    NumericMatrix pars_,
    double min_ll_
  ) :
  SSD(SSD_),
  pars(pars_),
  n_go(pars_.nrow()),
  min_ll(min_ll_) {}

  double operator()(const double& x) const {
    // log density of stop process winning at time x
    // input args: x, muS, sigmaS, tauS, exgS_lb, upper = Inf, log_d = TRUE
    double log_d = dtexg(
      x, pars(0, 3), pars(0, 4), pars(0, 5), pars(0, 9), R_PosInf, true
    );
    if (!R_FINITE(log_d)) log_d = min_ll;
    // summed log density of go accumulators not yet having finished by time x
    double summed_log_s = 0.0;
    for (int i = 0; i < n_go; ++i) {
      // NB stop signal delay added to finish time x, to account for earlier start of go process
      // input args: q, mu, sigma, tau, exg_lb, upper = Inf, lower_tail = FALSE, log_p = TRUE
      double log_s_i = ptexg(
        x + SSD, pars(i, 0), pars(i, 1), pars(i, 2), pars(i, 8), R_PosInf, false, true
      );
      if (!R_FINITE(log_s_i)) log_s_i = min_ll;
      summed_log_s += log_s_i;
    }
    // return sum of (1) log winner density and (2) sum of log survivor densities,
    // exponentiated to put on likelihood scale so integrator can work with it
    return std::exp(log_d + summed_log_s);
  }
};

// function to compute the stop success integral, not accounting for trigger
// failure and go failure
static inline double ss_texg_stop_success_lpdf(
    double SSD,
    NumericMatrix pars,
    double min_ll
) {
  // set up the integrand
  texg_stop_success_integrand f(SSD, pars, min_ll);
  // declare local variables for the integrator
  double err_est, res;
  int err_code;
  // perform numerical integration
  res = integrate(
    f,          // integrand
    pars(0, 9), // lower limit of integration (= lower bound of stop process)
    R_PosInf,   // upper limit of integration
    err_est,    // placeholder for estimation error
    err_code,   // placeholder for failed integration error code
    100,        // maximum number of subdivisions
    1e-8,       // absolute tolerance
    1e-6        // relative tolerance
  );
  // check for issues with result, and return
  bool bad = (err_code != 0) || !R_FINITE(res) || (res <= 0.0);
  return bad ? min_ll : std::log(res);
}

// top-level log-likelihood function for the stop signal task
NumericVector ss_texg_lpdf(
    NumericVector RT,
    IntegerVector R,
    NumericVector SSD,
    NumericVector lR,
    LogicalVector winner,
    NumericMatrix pars,
    LogicalVector is_ok,
    double min_ll
) {
  // pars columns: mu=0, sigma=1, tau=2, muS=3, sigmaS=4, tauS=5, tf=6, gf=7, exg_lb=8, exgS_lb=9

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
    if (is_ok[trial] != 1) {
      out[trial] = min_ll;
      continue;
    }

    // branching logic:
    // first, was a response observed? second, was a stop signal presented?
    bool response_observed = R[start_row] != NA_INTEGER;
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
        go_lprob = ss_texg_go_lpdf(
          RT[start_row],
          pars(Range(start_row, end_row), _),
          winner[Range(start_row, end_row)],
          min_ll
        );
        stop_fail_lprob = ss_texg_stop_fail_lpdf(
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
        go_lprob = ss_texg_go_lpdf(
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
        stop_success_integral = ss_texg_stop_success_lpdf(
           SSD[start_row],
           pars(Range(start_row, end_row), _),
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


// ----------------------------------------------------------------------------
// REGULAR EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

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
  if (!traits::is_finite<REALSXP>(stop_survivor_lprob)) {
    stop_survivor_lprob = min_ll;
  }
  // final output of race model is summed log likelihood
  return go_lprob + stop_survivor_lprob;
}

// go vs stop race likelihood function, for the case of stop trials with
// successful inhibition (i.e., stop process winning). not accounting for
// trigger failure and go failure.
// this function is an integrand (quantity to be integrated)
class exg_stop_success_integrand : public Func {
private:
  const double SSD;
  const NumericMatrix pars;
  const int n_go;
  const double min_ll;

public:
  exg_stop_success_integrand(
    double SSD_,
    NumericMatrix pars_,
    double min_ll_
  ) :
  SSD(SSD_),
  pars(pars_),
  n_go(pars_.nrow()),
  min_ll(min_ll_) {}

  double operator()(const double &x) const {
    // log density of stop process winning at time x
    // input args: x, muS, sigmaS, tauS, log_d = TRUE
    double log_d = dexg(x, pars(0, 3), pars(0, 4), pars(0, 5), true);
    if (!R_FINITE(log_d)) log_d = min_ll;
    // summed log density of go accumulators not yet having finished by time x
    double summed_log_s = 0.0;
    for (int i = 0; i < n_go; ++i) {
      // NB stop signal delay added to finish time x, to account for earlier start of go process
      // input args: q, mu, sigma, tau, lower_tail = FALSE, log_p = TRUE
      double log_s_i = pexg(x + SSD, pars(i, 0), pars(i, 1), pars(i, 2), false, true);
      if (!R_FINITE(log_s_i)) log_s_i = min_ll;
      summed_log_s += log_s_i;
    }
    // return sum of (1) log winner density and (2) sum of log survivor densities,
    // exponentiated to put on likelihood scale so integrator can work with it
    return std::exp(log_d + summed_log_s);
  }
};

// function to compute the stop success integral, not accounting for trigger
// failure and go failure
static inline double ss_exg_stop_success_lpdf(
    double SSD,
    NumericMatrix pars,
    double min_ll
) {
  // set up the integrand
  exg_stop_success_integrand f(SSD, pars, min_ll);
  // declare local variables for the integrator
  double err_est, res;
  int err_code;
  // perform numerical integration
  res = integrate(
    f,        // integrand
    R_NegInf, // lower limit of integration
    R_PosInf, // upper limit of integration
    err_est,  // placeholder for estimation error
    err_code, // placeholder for failed integration error code
    100,      // maximum number of subdivisions
    1e-8,     // absolute tolerance
    1e-6      // relative tolerance
  );
  // check for issues with result, and return
  bool bad = (err_code != 0) || !R_FINITE(res) || (res <= 0.0);
  return bad ? min_ll : std::log(res);
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
    if (is_ok[trial] != 1) {
      out[trial] = min_ll;
      continue;
    }

    // branching logic:
    // first, was a response observed? second, was a stop signal presented?
    bool response_observed = R[start_row] != NA_INTEGER;
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


// ----------------------------------------------------------------------------
// OLD STUFF BELOW, KEPT FOR TESTING R CODE
// ----------------------------------------------------------------------------

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

NumericVector protect_finite(const NumericVector& x, double min_ll) {
  NumericVector out = clone(x);
  for (int i = 0; i < out.size(); i++) {
    if (!R_FINITE(out[i])) {
      out[i] = min_ll;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector dEXGrace(
    NumericMatrix dt, NumericVector mu, NumericVector sigma, NumericVector tau,
    double min_ll
){
  int n = mu.size();
  NumericVector log_out(dt.nrow());
  log_out = protect_finite(
    dexg_c(dt(0, _), mu[0], sigma[0], tau[0], true),
    min_ll
  );
  for (int i = 1; i < n; i++) {
    log_out += protect_finite(
      pexg_c(dt(i, _), mu[i], sigma[i], tau[i], false, true),
      min_ll
    );
  }
  return exp(log_out);
}

// [[Rcpp::export]]
NumericVector stopfn_exg(
    NumericVector t, NumericVector mu, NumericVector sigma, NumericVector tau,
    double SSD, double min_ll
){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dEXGrace(dt, mu, sigma, tau, min_ll);
}

// old

// [[Rcpp::export]]
NumericVector pEXG_old(NumericVector q,
                   double mu = 5., double sigma = 1., double tau = 1.,
                   bool lower_tail = true, bool log_p = false) {
  int n = q.size();
  if (tau <= 0 || sigma <= 0) {
    NumericVector cdf(n, NA_REAL);
    return cdf;
  }

  NumericVector cdf(n);
  if (sigma < 1e-4){
    for (int i = 0; i < n; i++){
      cdf[i] = R::pexp(q[i] - mu, tau, lower_tail, log_p);
    }
    return cdf;
  }

  for (int i = 0; i < n; i++){
    if (!traits::is_infinite<REALSXP>(q[i])){
//      if (tau > .05 * sigma){
        double z_i = q[i] - mu - (sigma * sigma) / tau;
        cdf[i] = R::pnorm((q[i] - mu) / sigma, 0., 1., true, false) - std::exp(std::log(R::pnorm(z_i / sigma, 0., 1., true, false)) + (std::pow((mu + (sigma * sigma / tau)), 2) - mu * mu - 2. * q[i] * (sigma * sigma / tau)) / (2. * sigma * sigma));
//      } else {
//        cdf[i] = R::pnorm(q[i], mu, sigma, true, false);
//      }
    } else {
      if (q[i] < 0) {
        cdf[i] = 0.;
      } else {
        cdf[i] = 1.;
      }
    }
  }
  if (!lower_tail){
    for(int i = 0; i < n; i++){
      cdf[i] = 1. - cdf[i];
    }
  }
  if (log_p){
    for(int i = 0; i < n; i++){
      cdf[i] = std::log(cdf[i]);
    }
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector dEXG_old(NumericVector x,
                   double mu = 5., double sigma = 1., double tau = 1.,
                   bool log_d = false) {
  int n = x.size();
  if (tau <= 0 || sigma <= 0) {
    NumericVector pdf(n, NA_REAL);
    return pdf;
  }

  NumericVector pdf(n);
  if (sigma < 1e-4){
    for (int i = 0; i < n; i++){
      pdf[i] = R::dexp(x[i] - mu, tau, log_d);
    }
    return pdf;
  }

  for (int i = 0; i < n; i++){
//    if (tau > .05 * sigma){
      double z_i = x[i] - mu - (sigma * sigma) / tau;
      pdf[i] = - std::log(tau) - (z_i + (sigma * sigma)/(2. * tau)) / tau + std::log(R::pnorm(z_i / sigma, 0., 1., true, false));
//    } else {
//      pdf[i] = R::dnorm(x[i], mu, sigma, true);
//    }
  }
  if (!log_d){
    for(int i = 0; i < n; i++){
      pdf[i] = std::exp(pdf[i]);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector dEXGrace_old(NumericMatrix dt,
                       NumericVector mu, NumericVector sigma, NumericVector tau){
  int n = mu.size();
  NumericVector out(dt.nrow());
  out = dEXG_old(dt(0, _), mu[0], sigma[0], tau[0], false);
  for (int i = 1; i < n; i++){
    out = out * pEXG_old(dt(i, _), mu[i], sigma[i], tau[i], false, false);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_exg_old(NumericVector t,
                         NumericVector mu, NumericVector sigma, NumericVector tau,
                         double SSD){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dEXGrace_old(dt, mu, sigma, tau);
}


// [[Rcpp::export]]
NumericVector pTEXG_vec(
    NumericVector q, double mu = 5., double sigma = 1., double tau = 1., double lb = .05,
    bool lower_tail = true, bool log_p = false
) {
  int n = q.size();
  if (tau <= 0. || sigma <= 0.) {
    NumericVector cdf(n, NA_REAL);
    return cdf;
  }
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    cdf[i] = ptexg(q[i], mu, sigma, tau, lb, R_PosInf, lower_tail, log_p);
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector dTEXG_vec(
    NumericVector x, double mu = 5., double sigma = 1., double tau = 1., double lb = .05,
    bool log_d = false
) {
  int n = x.size();
  if (tau <= 0. || sigma <= 0.) {
    NumericVector pdf(n, NA_REAL);
    return pdf;
  }
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    pdf[i] = dtexg(x[i], mu, sigma, tau, lb, R_PosInf, log_d);
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector dTEXGrace(
    NumericMatrix dt,
    NumericVector mu, NumericVector sigma, NumericVector tau, NumericVector lb
){
  int n = mu.size();
  NumericVector out(dt.nrow());
  out = dTEXG_vec(dt(0, _), mu[0], sigma[0], tau[0], lb[0]);
  for (int i = 1; i < n; i++){
    out = out * pTEXG_vec(dt(i, _), mu[i], sigma[i], tau[i], lb[i], false);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector stopfn_texg(
    NumericVector t,
    NumericVector mu, NumericVector sigma, NumericVector tau, NumericVector lb,
    double SSD
){
  NumericVector tmp(mu.size() * t.size());
  tmp = rep_each(t, mu.size()) + SSD;
  NumericMatrix dt(mu.size(), t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dTEXGrace(dt, mu, sigma, tau, lb);
}

#endif
