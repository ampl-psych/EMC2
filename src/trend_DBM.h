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
// Inputs `a0` and `b0` correspond to Beta shape parameters.
// By default returns the trial-wise posterior mean, but can also return the
// posterior mode if `return_map` is true.
NumericVector run_beta_binomial_basic(
    const NumericVector &covariate,
    const double a0,
    const double b0,
    const bool return_map
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  // local count variables
  double n_hit = 0.0;
  double n_trial = 0.0;
  // trial-wise shape parameters of Beta distribution
  double a_t, b_t;

  for (int t = 0; t < n_trials; ++t) {
    // prediction before observing trial t
    a_t = a0 + n_hit;
    b_t = b0 + (n_trial - n_hit);
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
    const double a0,
    const double b0,
    const int decay,
    const bool return_map
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  double n_hit = 0.0;
  double n_trial = 0.0;
  const double decay_factor = std::exp(-1.0 / decay);
  double a_t, b_t;

  for (int t = 0; t < n_trials; ++t) {
    a_t = a0 + n_hit;
    b_t = b0 + (n_trial - n_hit);
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
    const double a0,
    const double b0,
    const int window,
    const bool return_map
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  double n_hit = 0.0;
  double n_trial = 0.0;
  // use double-ended queue for fast addition at back & removal at front
  std::deque<double> buf;
  double a_t, b_t;

  for (int t = 0; t < n_trials; ++t) {
    a_t = a0 + n_hit;
    b_t = b0 + (n_trial - n_hit);
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
    const double a0,
    const double b0,
    const int decay = 0,
    const int window = 0,
    const bool return_map = false,
    const bool return_surprise = false
) {
  const int n_trials = covariate.length();
  if (n_trials < 1) {
    stop("`covariate` should consist of at least one observation.");
  }
  if (a0 <= 0.0 || b0 <= 0.0) {
    stop("Both shape parameters `a0` and `b0` must be positive.");
  }
  if (decay > 0 && window > 0) {
    stop("Cannot use both `decay` and `window`. Choose only one memory constraint.");
  }
  NumericVector out(n_trials);
  if (decay > 0) {
    out = run_beta_binomial_decay(covariate, a0, b0, decay, return_map);
  } else if (window > 0) {
    out = run_beta_binomial_window(covariate, a0, b0, window, return_map);
  } else {
    out = run_beta_binomial_basic(covariate, a0, b0, return_map);
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
    const double cp = 0.1,
    const double mu0 = 0.25,
    const double s0 = 10.0,
    const bool return_map = false,
    const bool return_surprise = false,
    const int grid_res = 100,
    const double cp_eps = 1e-12
) {
  const int n_trials = covariate.length();
  if (n_trials < 1) {
    stop("`covariate` should consist of at least one observation.");
  }
  for (int i = 0; i < n_trials; i++) {
    if (!(covariate[i] == 0.0 || covariate[i] == 1.0)) {
      stop("All `covariate` entries must be 0.0 or 1.0.");
    }
  }
  if (cp < 0 || cp > 1) {
    stop("Change point probability `cp` must be in the range [0, 1].");
  }
  if (mu0 <= 0 || mu0 >= 1) {
    stop("Prior mean `mu0` must be in the range (0, 1).");
  }
  if (s0 <= 0) {
    stop("Prior scale `s0` must be strictly positive.");
  }
  if (grid_res < 1) {
    stop("`grid_res` must be greater than or equal to 1.");
  }
  if (cp_eps <= 0) {
    stop("`cp_eps` must be strictly positive.");
  }

  // shape parameters of Beta prior
  // NB in original Matlab code from Jaime Ide, the parameterisation
  // a = mu0 * s0 + 1; b = (1 - mu0) * s0 + 1
  // was used. The +1 shift was presumably a pragmatic tweak to avoid the shape
  // parameters ever being <= 1.
  const double a = mu0 * s0;
  const double b = (1.0 - mu0) * s0;

  // if cp is (practically) equal to 0 (i.e., no volatility), the DBM reduces to
  // a standard Beta-Binomial model, for which exact analytic updates are
  // available, which are cheaper to compute than the discretised grid approach
  // of the DBM.
  if (cp < cp_eps) {
    return run_beta_binomial(covariate, a, b, 0, 0, return_map, return_surprise);
  }

  // if cp is (practically) equal to 1 (i.e., pure volatility), there is no
  // actual learning from observations; the beliefs are purely driven by the
  // fixed prior, hence the output is constant (mean / mode of fixed prior)
  if ((1 - cp) < cp_eps) {
    const double out_val = return_map ? beta_mode(a, b) : beta_mean(a, b);
    NumericVector out(n_trials, out_val);
    if (return_surprise) {
      shannon_surprise(out);
    }
    return out;
  }

  // declare local variables
  NumericVector out(n_trials);
  const int grid_size = grid_res + 1;
  NumericVector prob_grid(grid_size), DBM_prior(grid_size), DBM_post(grid_size);

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
        DBM_pred[i] = (1.0 - cp) * DBM_post[i] + cp * DBM_prior[i];
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


// ----- Transition Probability Model -----------------------------------------

NumericVector run_tpm_nocp(
    const NumericVector covariate,
    const double a0 = 1.0,
    const double b0 = 1.0,
    const bool return_surprise = false
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  double a_XX = a0;
  double b_XX = b0;
  double a_YX = a0;
  double b_YX = b0;
  out[0] = beta_mean(a0, b0);

  for (int t = 1; t < n_trials; t++) {
    const int prev = static_cast<int>(covariate[(t - 1)]);
    const int curr = static_cast<int>(covariate[t]);
    // prediction
    if (prev == 1) {
      out[t] = beta_mean(a_XX, b_XX);
    } else {
      out[t] = beta_mean(a_YX, b_YX);
    }
    // update
    if (prev == 1) {
      if (curr == 1) {
        a_XX += 1.0;
      } else {
        b_XX += 1.0;
      }
    } else {
      if (curr == 1) {
        a_YX += 1.0;
      } else {
        b_YX += 1.0;
      }
    }
  }

  if (return_surprise) {
    shannon_surprise(out);
  }
  return out;
}

// pre-computed grid of transition probabilities and all four possible single-
// trial likelihoods
struct tpm_grid {
  NumericVector p_XX;
  NumericVector p_XY;
  NumericVector like_XX;
  NumericVector like_XY;
  NumericVector like_YX;
  NumericVector like_YY;
};

inline tpm_grid build_tpm_grid(const int grid_res) {
  const int resol = grid_res + 1;
  const int n_combi = resol * resol;

  NumericVector grid(resol);
  for (int i = 0; i < resol; i++) {
    grid[i] = static_cast<double>(i) / (resol - 1);
  }

  tpm_grid out{
    NumericVector(n_combi),
    NumericVector(n_combi),
    NumericVector(n_combi),
    NumericVector(n_combi),
    NumericVector(n_combi),
    NumericVector(n_combi)
  };

  int idx = 0;
  for (int i0 = 0; i0 < resol; i0++) {
    const double p_XY_val = grid[i0];
    for (int i1 = 0; i1 < resol; i1++) {
      const double p_XX_val = grid[i1];
      // p(X|X), i.e., prob of current obs = 1 given previous obs = 1
      out.p_XX[idx] = p_XX_val;
      // p(X|Y), i.e., prob of current obs = 1 given previous obs = 0
      out.p_XY[idx] = p_XY_val;
      // likelihood of current obs = 1 given previous obs = 0
      out.like_XY[idx] = p_XY_val;
      // likelihood of current obs = 1 given previous obs = 1
      out.like_XX[idx] = p_XX_val;
      // likelihood of current obs = 0 given previous obs = 0
      out.like_YY[idx] = 1.0 - p_XY_val;
      // likelihood of current obs = 0 given previous obs = 1
      out.like_YX[idx] = 1.0 - p_XX_val;
      ++idx;
    }
  }

  return out;
}

NumericVector run_tpm(
    const NumericVector covariate,
    const double cp,
    const double a0 = 1.0,
    const double b0 = 1.0,
    const bool return_surprise = false,
    const int grid_res = 100,
    const double cp_eps = 1e-12
) {
  const int n_trials = covariate.length();
  if (n_trials < 1) {
    stop("`covariate` should consist of at least one observation.");
  }
  for (int i = 0; i < n_trials; i++) {
    if (!(covariate[i] == 0.0 || covariate[i] == 1.0)) {
      stop("All `covariate` entries must be 0.0 or 1.0.");
    }
  }
  if (cp < 0 || cp > 1) {
    stop("Change point probability `cp` must be in the range [0, 1].");
  }
  if (a0 <= 0.0 || b0 <= 0.0) {
    stop("Both prior shape parameters `a0` and `b0` must be positive.");
  }
  if (grid_res < 1) {
    stop("`grid_res` must be greater than or equal to 1.");
  }
  if (cp_eps <= 0) {
    stop("`cp_eps` must be strictly positive.");
  }

  // if cp is (practically) equal to 0 (i.e., no volatility), fallback to simple
  // Beta-Binomial over transitions.
  if (cp < cp_eps) {
    return run_tpm_nocp(covariate, a0, b0, return_surprise);
  }

  // if cp is (practically) equal to 1 (i.e., pure volatility), beliefs are given
  // by fixed prior
  if ((1 - cp) < cp_eps) {
    NumericVector out(n_trials, beta_mean(a0, b0));
    if (return_surprise) {
      shannon_surprise(out);
    }
    return out;
  }

  // declare local variables
  NumericVector out(n_trials);
  tpm_grid grid = build_tpm_grid(grid_res);
  const int resol = grid_res + 1;
  const int n_combi = resol * resol;
  const double inv_n_minus_1 = 1.0 / (n_combi - 1.0);
  NumericVector TPM_post(n_combi), TPM_pred(n_combi), TPM_update(n_combi), mean_p(n_combi);

  // pre-compute:
  for (int j = 0; j < n_combi; j++) {
    // mean trans probs (i.e., (p(X|X) + p(X|Y)) / 2) for every possible joint
    // combination
    mean_p[j] = 0.5 * (grid.p_XX[j] + grid.p_XY[j]);
    // initialise posterior with Beta prior (for t=0 case)
    TPM_post[j] = R::dbeta(grid.p_XX[j], a0, b0, false) *
      R::dbeta(grid.p_XY[j], a0, b0, false);
  }
  normalise_inplace(TPM_post);

  // loop over trials
  for (int t = 0; t < n_trials; t++) {

    const int curr = static_cast<int>(covariate[t]);
    const int prev = t > 0 ? static_cast<int>(covariate[(t - 1)]) : NA_INTEGER;

    // Build predictive distribution from previous posterior.
    // - With prob 1-cp: environment is stable -> posterior propagates.
    // - With prob cp: change point -> new transition probs independent of old ones,
    //   so posterior mass is redistributed uniformly over all *other* grid points.
    double summed_post = std::accumulate(TPM_post.begin(), TPM_post.end(), 0.0);
    for (int j = 0; j < n_combi; j++) {
      TPM_pred[j] = (1.0 - cp) * TPM_post[j] +
        cp * (summed_post - TPM_post[j]) * inv_n_minus_1;
    }
    normalise_inplace(TPM_pred);

    // Before observing trial t, compute marginal predictive probability of
    // observing X - that is, predictive probability of covariate[t] == 1.0
    if (t == 0) {
      out[t] = mean_discrete(mean_p, TPM_pred);
    } else {
      out[t] = prev == 1 ? mean_discrete(grid.p_XX, TPM_pred) : mean_discrete(grid.p_XY, TPM_pred);
    }

    // After observing trial t, update posterior distribution
    if (t > 0) {
      const NumericVector *like_ptr = nullptr;
      if (prev == 0 && curr == 0) {
        like_ptr = &grid.like_YY;
      } else if (prev == 0 && curr == 1) {
        like_ptr = &grid.like_XY;
      } else if (prev == 1 && curr == 0) {
        like_ptr = &grid.like_YX;
      } else {
        like_ptr = &grid.like_XX;
      }
      for (int j = 0; j < n_combi; j++) {
        const double like_j = (*like_ptr)[j];
        TPM_update[j] = ((1.0 - cp) * like_j * TPM_post[j]) +
          (cp * like_j * (summed_post - TPM_post[j]) * inv_n_minus_1);
      }
      normalise_inplace(TPM_update);
      std::swap(TPM_post, TPM_update);
    }
  }

  if (return_surprise) {
    shannon_surprise(out);
  }
  return out;
}

#endif
