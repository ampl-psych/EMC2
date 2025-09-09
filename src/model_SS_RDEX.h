#ifndef ss_rdex_h
#define ss_rdex_h

#include <cmath>
#include <Rcpp.h>
#include "wald_functions.h"
#include "exgaussian_functions.h"
#include "race_integrate.h"
using namespace Rcpp;

// ----------------------------------------------------------------------------
// HYBRID WALD / EX-GAUSSIAN FUNCTIONS
// ----------------------------------------------------------------------------

// wrapper around Wald log PDF for race function
NumericVector rdex_go_lpdf(
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

    double dt_i = rt[i] - pars(i, 3);
    double log_d = R_NegInf;
    if (dt_i > 0.) {
      log_d = std::log(
        digt(
          dt_i,
          (pars(i, 1) / pars(i, 4)) + .5 * (pars(i, 2) / pars(i, 4)),
          pars(i, 0) / pars(i, 4),
          .5 * (pars(i, 2) / pars(i, 4))
        )
      );
    }

    out[k] = std::isfinite(log_d) ? log_d : min_ll;

    k++;
  }

  return(out);
}

// wrapper around Wald log complementary CDF for race function
NumericVector rdex_go_lccdf(
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

    double dt_i = rt[i] - pars(i, 3);
    double log_s = 0.; // log(1)
    if (dt_i > 0.) {
      log_s = log1m(
        pigt(
          dt_i,
          (pars(i, 1) / pars(i, 4)) + .5 * (pars(i, 2) / pars(i, 4)),
          pars(i, 0) / pars(i, 4),
          .5 * (pars(i, 2) / pars(i, 4))
        )
      );
    }

    out[k] = std::isfinite(log_s) ? log_s : min_ll;

    k++;
  }

  return(out);
}

// go race log likelihood function, not accounting for go failure
double ss_rdex_go_lpdf(
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
  race_f lpdf(rdex_go_lpdf, rdex_go_lccdf, pars, winner, min_ll, true);
  // evaluate the race function for the observed RT
  return lpdf(RT);
}

// go vs stop race log likelihood function for the case of stop trials with
// failed inhibition (i.e., stop process losing), not accounting for trigger
// failure and go failure.
double ss_rdex_stop_fail_lpdf(
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
  double go_lprob = ss_rdex_go_lpdf(RT, pars, winner, min_ll);
  // obtain the survivor log probability of the stop process
  // NB SSD subtracted from observed RT to get stop finish time
  // input args: q, muS, sigmaS, tauS, exgS_lb, upper = Inf, lower_tail = FALSE, log_p = TRUE
  double stop_survivor_lprob = ptexg(
    RT - SSD, pars(0, 5), pars(0, 6), pars(0, 7), pars(0, 10), R_PosInf, false, true
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
class rdex_stop_success_integrand : public Func {
private:
  const double SSD;
  const NumericMatrix pars;
  const int n_go;
  const double min_ll;

public:
  rdex_stop_success_integrand(
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
      x, pars(0, 5), pars(0, 6), pars(0, 7), pars(0, 10), R_PosInf, true
    );
    if (!R_FINITE(log_d)) log_d = min_ll;
    // summed log density of go accumulators not yet having finished by time x
    double summed_log_s = 0.0;
    for (int i = 0; i < n_go; ++i) {
      // NB stop signal delay added to finish time x, to account for earlier start of go process
      double dt_i = (x + SSD) - pars(i, 3);
      double log_s_i = 0.; // log(1)
      if (dt_i > 0.) {
        log_s_i = log1m(
          pigt(
            dt_i,
            // go pars columns: v=0, B=1, A=2, t0=3, s=4
            (pars(i, 1) / pars(i, 4)) + .5 * (pars(i, 2) / pars(i, 4)),
            pars(i, 0) / pars(i, 4),
            .5 * (pars(i, 2) / pars(i, 4))
          )
        );
      }
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
static inline double ss_rdex_stop_success_lpdf(
    double SSD,
    NumericMatrix pars,
    double min_ll,
    double upper = R_PosInf
) {
  // set up the integrand
  rdex_stop_success_integrand f(SSD, pars, min_ll);
  // declare local variables for the integrator
  double err_est, res;
  int err_code;
  // perform numerical integration
  res = integrate(
    f,            // integrand
    pars(0, 10),  // lower limit of integration (= lower bound of stop process)
    upper,        // upper limit of integration
    err_est,      // placeholder for estimation error
    err_code,     // placeholder for failed integration error code
    100,          // maximum number of subdivisions
    1e-8,         // absolute tolerance
    1e-6          // relative tolerance
  );
  // check for issues with result, and return
  bool bad = (err_code != 0) || !R_FINITE(res) || (res <= 0.0);
  return bad ? min_ll : std::log(res);
}

// top-level log-likelihood function for the stop signal task
NumericVector ss_rdex_lpdf(
    NumericVector RT,
    IntegerVector R,
    NumericVector SSD,
    NumericVector lR,
    LogicalVector winner,
    NumericMatrix pars,
    LogicalVector is_ok,
    double min_ll
) {
  // pars columns: v=0, B=1, A=2, t0=3, s=4, muS=5, sigmaS=6, tauS=7, tf=8, gf=9, exgS_lb=10, b=11

  NumericVector unique_lR = unique(lR);
  const int n_acc = unique_lR.length();         // number of go accumulators
  const int n_trials = lR.length() / n_acc;
  NumericVector out(n_trials);                  // trial-wise log-likelihood

  // extract go failure and trigger failure parameters as trial-wise vectors
  // from pars matrix; these parameters cannot differ between accumulators
  // within a given trial
  LogicalVector is_first_acc = lR == unique_lR[0]; // arbitrary index for one accumulator
  NumericVector tf = pars(_, 8);
  NumericVector gf = pars(_, 9);
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
        go_lprob = ss_rdex_go_lpdf(
          RT[start_row],
            pars(Range(start_row, end_row), _),
            winner[Range(start_row, end_row)],
                  min_ll
        );
        stop_fail_lprob = ss_rdex_stop_fail_lpdf(
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
        go_lprob = ss_rdex_go_lpdf(
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
        stop_success_integral = ss_rdex_stop_success_lpdf(
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

// [[Rcpp::export]]
NumericVector pEXG_RDEX(NumericVector q,
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
      if (tau > .05 * sigma){
        double z_i = q[i] - mu - (sigma * sigma) / tau;
        cdf[i] = R::pnorm((q[i] - mu) / sigma, 0., 1., true, false) - std::exp(std::log(R::pnorm(z_i / sigma, 0., 1., true, false)) + (std::pow((mu + (sigma * sigma / tau)), 2) - mu * mu - 2. * q[i] * (sigma * sigma / tau)) / (2. * sigma * sigma));
      } else {
        cdf[i] = R::pnorm(q[i], mu, sigma, true, false);
      }
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
NumericVector dEXG_RDEX(NumericVector x,
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
    if (tau > .05 * sigma){
      double z_i = x[i] - mu - (sigma * sigma) / tau;
      pdf[i] = - std::log(tau) - (z_i + (sigma * sigma)/(2. * tau)) / tau + std::log(R::pnorm(z_i / sigma, 0., 1., true, false));
    } else {
      pdf[i] = R::dnorm(x[i], mu, sigma, true);
    }
  }
  if (!log_d){
    for(int i = 0; i < n; i++){
      pdf[i] = std::exp(pdf[i]);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector dWald_RDEX_old(NumericVector t, double v,
                             double B, double A, double t0){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0;
    if (t[i] <= 0){
      pdf[i] = 0.;
    } else {
      pdf[i] = digt(t[i], B + .5 * A, v, .5 * A);
    }
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector dWald_RDEX(
    NumericVector t,
    double v, double B, double A, double t0, double s
) {
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++) {
    t[i] = t[i] - t0;
    pdf[i] = 0.;
    if (t[i] > 0.) {
      pdf[i] = digt(t[i], (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald_RDEX_old(NumericVector t, double v,
                             double B, double A, double t0){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    t[i] = t[i] - t0;
    if (t[i] <= 0){
      cdf[i] = 0.;
    } else {
      cdf[i] = pigt(t[i], B + .5 * A, v, .5 * A);
    }
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector pWald_RDEX(
    NumericVector t,
    double v, double B, double A, double t0, double s
) {
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++) {
    t[i] = t[i] - t0;
    cdf[i] = 0.;
    if (t[i] > 0.) {
      cdf[i] = pigt(t[i], (B/s) + .5 * (A/s), (v/s), .5 * (A/s));
    }
  }
  return cdf;
}


// [[Rcpp::export]]
NumericVector pTEXG_RDEX(
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
NumericVector dTEXG_RDEX(
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
NumericVector dRDEXrace_old(NumericMatrix dt,
                            double mu, double sigma, double tau,
                            NumericVector v, NumericVector B, NumericVector A,
                            NumericVector t0, bool exgWinner = true){
  int n = v.size();
  NumericVector out(dt.nrow());
  if (exgWinner){
    out = dEXG_RDEX(dt(0, _), mu, sigma, tau, false);
    out = out * (1. - pWald_RDEX_old(dt(1, _), v[0], B[0], A[0], t0[0]));
  } else {
    out = dWald_RDEX_old(dt(0, _), v[0], B[0], A[0], t0[0]);
    out = out * (1. - pEXG_RDEX(dt(1, _), mu, sigma, tau));
  }
  for (int i = 1; i < n; i++){
    out = out * (1. - pWald_RDEX_old(dt(i + 1, _), v[i], B[i], A[i], t0[i]));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector dRDEXrace(
    NumericMatrix dt,
    double mu, double sigma, double tau, double lb,
    NumericVector v, NumericVector B, NumericVector A, NumericVector t0, NumericVector s,
    bool exgWinner = true
) {
  int n = v.size();
  NumericVector out(dt.nrow());
  if (exgWinner) {
    out = dTEXG_RDEX(dt(0, _), mu, sigma, tau, lb);
    out = out * (1. - pWald_RDEX(dt(1, _), v[0], B[0], A[0], t0[0], s[0]));
  } else {
    out = dWald_RDEX(dt(0, _), v[0], B[0], A[0], t0[0], s[0]);
    out = out * (1. - pTEXG_RDEX(dt(1, _), mu, sigma, tau, lb));
  }
  for (int i = 1; i < n; i++){
    out = out * (1. - pWald_RDEX(dt(i + 1, _), v[i], B[i], A[i], t0[i], s[i]));
  }
  return out;
}


// [[Rcpp::export]]
NumericVector stopfn_rdex_old(NumericVector t, int n_acc,
                              double mu, double sigma, double tau,
                              NumericVector v, NumericVector B, NumericVector A,
                              NumericVector t0, double SSD){
  NumericVector tmp( (n_acc + 1) * t.size());
  tmp = rep_each(t, n_acc + 1) + SSD;
  NumericMatrix dt(n_acc + 1, t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dRDEXrace_old(dt, mu, sigma, tau, v, B, A, t0);
}

// [[Rcpp::export]]
NumericVector stopfn_rdex(
    NumericVector t, int n_acc,
    double mu, double sigma, double tau, double lb,
    NumericVector v, NumericVector B, NumericVector A, NumericVector t0, NumericVector s,
    double SSD
) {
  NumericVector tmp((n_acc + 1) * t.size());
  tmp = rep_each(t, n_acc + 1) + SSD;
  NumericMatrix dt(n_acc + 1, t.size(), tmp.begin());
  dt(0, _) = dt(0, _) - SSD;
  return dRDEXrace(dt, mu, sigma, tau, lb, v, B, A, t0, s);
}

#endif
