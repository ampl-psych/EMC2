#ifndef trend_DBM_h
#define trend_DBM_h

#include <Rcpp.h>
#include <limits>
#include <deque>
using namespace Rcpp;

// ----- Beta distribution helpers --------------------------------------------

// NB assuming inputs a and b are strictly positive
inline double beta_mode(const double a, const double b) {
  // symmetric case: a == b.
  // - if a = b < 1: density diverges at both 0 and 1: no unique mode; return 0.5
  // - if a = b = 1: uniform(0, 1): no unique mode; return 0.5 (i.e., the mean)
  // - if a = b > 1: unimodal symmetric distribution: mode is exactly 0.5
  if (a == b) {
    return 0.5;
  }
  // if either a or b is less than 1, the density diverges at one endpoint.
  // divergence is stronger on the side with the smaller parameter.
  if (a < 1.0 || b < 1.0) {
    return (a < b) ? 0.0 : 1.0;
  }
  // typical case: both a and b > 1
  return (a - 1.0) / (a + b - 2.0);
}

// NB assuming inputs a and b are strictly positive
inline double beta_mean(const double a, const double b) {
  return a / (a + b);
}


// ----- Shannon surprise (bits) ----------------------------------------------

inline void shannon_surprise(NumericVector &x) {
  static const double inv_ln2 = 1.0 / std::log(2.0);
  static const double tiny = std::numeric_limits<double>::min();
  const int n = x.length();
  for (int i = 0; i < n; i++) {
    const double safe_x = (x[i] > 0.0) ? x[i] : tiny;
    x[i] = -std::log(safe_x) * inv_ln2;
  }
}


// ----- Beta-Binomial model --------------------------------------------------

// Input `covariate` assumed to be binary data vector represented as reals;
// Inputs `a` and `b` correspond to Beta shape parameters.
// By default returns the trial-wise posterior mean, but can also return the
// posterior mode if `return_map` is true.
NumericVector run_beta_binomial_basic(
    const NumericVector &covariate,
    const double a,
    const double b,
    const bool return_map
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  // local count variables
  double n_hit = 0.0;
  double n_trial = 0.0;
  // trial-wise shape parameters of Beta distribution
  double a_t;
  double b_t;

  for (int t = 0; t < n_trials; ++t) {
    // prediction before observing trial t
    a_t = a + n_hit;
    b_t = b + (n_trial - n_hit);
    out[t] = return_map ? beta_mode(a_t, b_t) : beta_mean(a_t, b_t);
    // update after observing trial t
    n_hit += covariate[t];
    n_trial += 1.0;
  }

  return out;
}

// Following Meyniel et al. (2016), exponential decay ("leaky integration")
// applied to the count of past events. Each past observation is discounted
// exponentially as:
// weight_after_k_trials = exp(-k/decay)
// where `decay` is the characteristic decay constant (in units of trials).
// This implies a half-life of decay * log(2). That is, after `decay * log(2)`
// trials, the contribution of past observations is reduced by 50%. For example,
// decay = 16 implies a half-life of approximately 11 trials.
// In practice, the decay is implemented by multiplying past accumulated counts
// by `exp(-1/decay)` at each trial.
NumericVector run_beta_binomial_decay(
    const NumericVector &covariate,
    const double a,
    const double b,
    const int decay,
    const bool return_map
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  double n_hit = 0.0;
  double n_trial = 0.0;
  const double decay_factor = std::exp(-1.0 / decay);
  double a_t;
  double b_t;

  for (int t = 0; t < n_trials; ++t) {
    a_t = a + n_hit;
    b_t = b + (n_trial - n_hit);
    out[t] = return_map ? beta_mode(a_t, b_t) : beta_mean(a_t, b_t);
    // decayed update
    n_hit = decay_factor * (n_hit + covariate[t]);
    n_trial = decay_factor * (n_trial + 1.0);
  }

  return out;
}

// Following Meyniel et al. (2016), assume memory of past events is limited to
// a fixed window
NumericVector run_beta_binomial_window(
    const NumericVector &covariate,
    const double a,
    const double b,
    const int window,
    const bool return_map
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  double n_hit = 0.0;
  double n_trial = 0.0;
  // use double-ended queue for fast addition at back & removal at front
  std::deque<double> buf;
  double a_t;
  double b_t;

  for (int t = 0; t < n_trials; ++t) {
    a_t = a + n_hit;
    b_t = b + (n_trial - n_hit);
    out[t] = return_map ? beta_mode(a_t, b_t) : beta_mean(a_t, b_t);
    // update with most recent observation
    n_hit += covariate[t];
    n_trial += 1.0;
    // add most recent observation to end of queue
    buf.push_back(covariate[t]);
    // if queue has become too large, remove oldest (i.e., first) observation
    if ((int)buf.size() > window) {
      n_hit -= buf.front();
      n_trial -= 1.0;
      buf.pop_front();
    }
  }

  return out;
}

// top-level dispatcher of Beta-Binomial model
NumericVector run_beta_binomial(
    const NumericVector covariate,
    const double a,
    const double b,
    const int decay = 0,
    const int window = 0,
    const bool return_map = false,
    const bool return_surprise = false
) {
  if (a <= 0.0 || b <= 0.0) {
    stop("Both shape parameters a and b must be positive.");
  }
  if (decay > 0 && window > 0) {
    stop("Cannot use both decay and window. Choose only one memory constraint.");
  }
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  if (decay > 0) {
    out = run_beta_binomial_decay(covariate, a, b, decay, return_map);
  } else if (window > 0) {
    out = run_beta_binomial_window(covariate, a, b, window, return_map);
  } else {
    out = run_beta_binomial_basic(covariate, a, b, return_map);
  }

  if (return_surprise) {
    shannon_surprise(out);
  }
  return out;
}


// ----- discretised density helpers ------------------------------------------

// normalise a vector in place and return its sum
inline double normalise_inplace(NumericVector &v) {
  double s = 0.0;
  for (double x : v) s += x;
  if (s > 0.0) {
    double inv_s = 1.0 / s;
    for (auto &x : v) x *= inv_s;
  }
  return s;
}

// compute mean for discretised probability distribution
inline double mean_discrete(const NumericVector &x, const NumericVector &w) {
  double s = 0.0;
  for (int i = 0; i < x.size(); i++) {
    s += x[i] * w[i];
  }
  return s;
}

// compute mode for discretised probability distribution
inline double mode_discrete(const NumericVector &x, const NumericVector &w) {
  int max_idx = std::distance(w.begin(), std::max_element(w.begin(), w.end()));
  return x[max_idx];
}


// ------ Dynamic Belief Model ------------------------------------------------

// Based on Yu & Cohen (2008) and Ide et al. (2013).
NumericVector run_dbm(
    const NumericVector covariate,
    const double alpha = 0.8,
    const double pmean = 0.25,
    const double pscale = 10.0,
    const bool return_map = false,
    const bool return_surprise = false,
    const int grid_res = 1e3,
    const double alpha_eps = 1e-10
) {

  if (alpha < 0 || alpha > 1) {
    stop("Mixing coefficient alpha must be in the range [0, 1].");
  }
  if (pmean <= 0 || pmean >= 1) {
    stop("Prior mean must be in the range (0, 1).");
  }
  if (pscale <= 0) {
    stop("Prior scale must be strictly positive.");
  }

  // shape parameters of Beta prior
  // NB in original Matlab code from Jaime Ide, the parameterisation
  // a = pmean * pscale + 1; b = (1 - pmean) * pscale + 1
  // was used. The +1 shift was presumably a pragmatic tweak to avoid the shape
  // parameters ever being <= 1.
  const double a = pmean * pscale;
  const double b = (1.0 - pmean) * pscale;

  // if alpha is (practically) equal to 1, the DBM reduces to a standard
  // Beta-Binomial model, for which exact analytic updates are available, which
  // are cheaper to compute than the discretised grid approach of the DBM.
  if ((1.0 - alpha) < alpha_eps) {
    return run_beta_binomial(covariate, a, b, 0, 0, return_map, return_surprise);
  }

  const int n_trials = covariate.length();
  NumericVector out(n_trials);

  // if alpha is (practically) equal to 0, there is no updating of the fixed
  // prior, hence the output is constant (mean / mode of fixed prior)
  if (alpha < alpha_eps) {
    const double out_val = return_map ? beta_mode(a, b) : beta_mean(a, b);
    std::fill(out.begin(), out.end(), out_val);
    if (return_surprise) {
      shannon_surprise(out);
    }
    return out;
  }

  // declare local variables
  NumericVector prob_grid(grid_res + 1);
  NumericVector DBM_prior(grid_res + 1);
  NumericVector DBM_post(grid_res + 1);
  const int grid_size = prob_grid.size();

  // compute discretised density of fixed Beta prior
  for (int i = 0; i < grid_size; i++) {
    prob_grid[i] = static_cast<double>(i) / (grid_size - 1);
    DBM_prior[i] = R::dbeta(prob_grid[i], a, b, false);
  }
  normalise_inplace(DBM_prior);

  // Bernoulli likelihoods for binary observation X vs. Y
  const NumericVector x_like = clone(prob_grid);
  const NumericVector y_like = 1.0 - prob_grid;

  // initialise predicted probability of observation X with fixed prior
  NumericVector DBM_pred = clone(DBM_prior);

  // trial-wise belief updating
  for (int t = 0; t < n_trials; t++) {
    // update predictive distribution: mixture of previous trial's posterior
    // and fixed prior
    if (t > 0) {
      for (int i = 0; i < grid_size; i++) {
        DBM_pred[i] = alpha * DBM_post[i] + (1.0 - alpha) * DBM_prior[i];
      }
      normalise_inplace(DBM_pred);
    }
    // main trial-wise output: predicted probability of observation X,
    // operationalised as either the mean or mode of the predictive distribution
    out[t] = return_map ? mode_discrete(prob_grid, DBM_pred) : mean_discrete(prob_grid, DBM_pred);
    // update posterior distribution
    if (covariate[t] == 1.0) {
      for (int i = 0; i < grid_size; i++) {
        DBM_post[i] = DBM_pred[i] * x_like[i];
      }
    } else {
      for (int i = 0; i < grid_size; i++) {
        DBM_post[i] = DBM_pred[i] * y_like[i];
      }
    }
    normalise_inplace(DBM_post);
  }

  if (return_surprise) {
    shannon_surprise(out);
  }
  return out;
}

#endif
