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

inline NumericVector shannon_surprise(
    const NumericVector &pred,
    const NumericVector &obs
) {
  static const double inv_ln2 = 1.0 / std::log(2.0);
  static const double tiny = std::numeric_limits<double>::min();
  const int n = pred.length();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    // clamp predictive probability away from exactly 0 or 1
    const double pred_safe = std::min(std::max(pred[i], tiny), (1.0 - tiny));
    // likelihood of observation
    const double like = (obs[i] == 1.0) ? pred_safe : (1.0 - pred_safe);
    // Shannon surprise
    out[i] = -std::log(like) * inv_ln2;
  }
  return out;
}


// ----- Beta-Binomial model --------------------------------------------------

// Input `covariate` assumed to be binary data vector represented as reals;
// Inputs `a0` and `b0` correspond to Beta shape parameters.
// By default returns the trial-wise posterior mean, but can also return the
// posterior mode if `return_map` is true.
NumericVector run_beta_binomial_basic(
    const NumericVector &covariate,
    const NumericVector &a0,
    const NumericVector &b0,
    const bool return_map
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  // count variables
  double n_hit = 0.0;
  double n_trial = 0.0;

  for (int t = 0; t < n_trials; ++t) {
    // prediction before observing trial t
    const double a_t = a0[t] + n_hit;
    const double b_t = b0[t] + (n_trial - n_hit);
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
    const NumericVector &a0,
    const NumericVector &b0,
    const NumericVector &decay,
    const bool return_map
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  double n_hit = 0.0;
  double n_trial = 0.0;

  for (int t = 0; t < n_trials; ++t) {
    const double a_t = a0[t] + n_hit;
    const double b_t = b0[t] + (n_trial - n_hit);
    out[t] = return_map ? beta_mode(a_t, b_t) : beta_mean(a_t, b_t);
    // decayed update
    const double decay_factor = std::exp(-1.0 / decay[t]);
    n_hit = decay_factor * (n_hit + covariate[t]);
    n_trial = decay_factor * (n_trial + 1.0);
  }

  return out;
}

// Following Meyniel et al. (2016), assume memory of past events is limited to
// a fixed window
struct event {
  double obs;
  int idx;
};

NumericVector run_beta_binomial_window(
    const NumericVector &covariate,
    const NumericVector &a0,
    const NumericVector &b0,
    const NumericVector &window,
    const bool return_map
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  double n_hit = 0.0;
  double n_trial = 0.0;
  // use double-ended queue for fast addition at back & removal at front
  std::deque<event> buf;

  for (int t = 0; t < n_trials; ++t) {
    // prune buffer based on current window[t]
    int w = static_cast<int>(window[t]);
    while (!buf.empty() && (t - buf.front().idx) >= w) {
      n_hit -= buf.front().obs;
      n_trial -= 1.0;
      buf.pop_front();
    }
    // compute prediction
    const double a_t = a0[t] + n_hit;
    const double b_t = b0[t] + (n_trial - n_hit);
    out[t] = return_map ? beta_mode(a_t, b_t) : beta_mean(a_t, b_t);
    // update with most recent observation
    const double obs = covariate[t];
    buf.push_back({obs, t});
    n_hit += obs;
    n_trial += 1.0;
  }

  return out;
}

// top-level dispatcher of Beta-Binomial model
NumericVector run_beta_binomial(
    const NumericVector covariate,
    const NumericVector a0,
    const NumericVector b0,
    const NumericVector decay,
    const NumericVector window,
    const bool return_map = false,
    const bool return_surprise = false
) {
  const int n_trials = covariate.length();
  for (int t = 0; t < n_trials; ++t) {
    if (!(covariate[t] == 0.0 || covariate[t] == 1.0)) {
      stop("All `covariate` entries must be 0.0 or 1.0.");
    }
  }

  bool use_decay = false;
  for (int t = 0; t < n_trials; ++t) {
    if (decay[t] > 0.0) {
      use_decay = true;
      break;
    }
  }
  bool use_window = false;
  for (int t = 0; t < n_trials; ++t) {
    if (window[t] > 0.0) {
      use_window = true;
      break;
    }
  }
  if (use_decay && use_window) {
    stop("Cannot use both `decay` and `window`. Choose only one memory constraint.");
  }

  NumericVector out(n_trials);
  if (use_decay) {
    out = run_beta_binomial_decay(covariate, a0, b0, decay, return_map);
  } else if (use_window) {
    out = run_beta_binomial_window(covariate, a0, b0, window, return_map);
  } else {
    out = run_beta_binomial_basic(covariate, a0, b0, return_map);
  }

  if (return_surprise) {
    out = shannon_surprise(out, covariate);
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
    const NumericVector cp,
    const NumericVector mu0,
    const NumericVector s0,
    const bool return_map = false,
    const bool return_surprise = false
) {
  // resolution of grid for discretising probability distributions.
  // exact setting doesn't seem to matter much for final outcome - mean/mode of
  // the trial-wise predictive distribution - but has huge effect on memory use
  const int grid_res = return_map ? 500 : 100;
  // tolerance for treating cp parameter as practically equal to 0 or 1
  const double cp_eps = 1e-10;

  const int n_trials = covariate.length();
  for (int t = 0; t < n_trials; t++) {
    if (!(covariate[t] == 0.0 || covariate[t] == 1.0)) {
      stop("All `covariate` entries must be 0.0 or 1.0.");
    }
  }

  // shape parameters of Beta prior
  // NB in original Matlab code from Jaime Ide, the parameterisation
  // a = mu0 * s0 + 1; b = (1 - mu0) * s0 + 1
  // was used. The +1 shift was presumably a pragmatic tweak to avoid the shape
  // parameters ever being <= 1.
  NumericVector a(n_trials), b(n_trials);
  for (int t = 0; t < n_trials; t++) {
    a[t] = mu0[t] * s0[t];
    b[t] = (1.0 - mu0[t]) * s0[t];
  }

  // if cp is (practically) equal to 0 (i.e., no volatility), the DBM reduces to
  // a standard Beta-Binomial model, for which exact analytic updates are
  // available, which are cheaper to compute than the discretised grid approach
  // of the DBM.
  bool cp_zero = true;
  for (int t = 0; t < n_trials; ++t) {
    if (cp[t] >= cp_eps) {
      cp_zero = false;
      break;
    }
  }
  if (cp_zero) {
    return run_beta_binomial_basic(covariate, a, b, return_map);
  }

  NumericVector out(n_trials);

  // if cp is (practically) equal to 1 (i.e., pure volatility), there is no
  // actual learning from observations; the beliefs are purely driven by the
  // fixed prior, hence the output is constant (mean / mode of fixed prior)
  bool cp_one = true;
  for (int t = 0; t < n_trials; ++t) {
    if ((1 - cp[t]) >= cp_eps) {
      cp_one = false;
      break;
    }
  }
  if (cp_one) {
    for (int t = 0; t < n_trials; ++t) {
      out[t] = return_map ? beta_mode(a[t], b[t]) : beta_mean(a[t], b[t]);
    }
    if (return_surprise) {
      out = shannon_surprise(out, covariate);
    }
    return out;
  }

  // declare local variables
  const int grid_size = grid_res + 1;
  NumericVector prob_grid(grid_size), DBM_prior(grid_size), DBM_pred(grid_size), DBM_post(grid_size);

  for (int i = 0; i < grid_size; i++) {
    prob_grid[i] = static_cast<double>(i) / (grid_size - 1);
  }

  // Bernoulli likelihoods for binary observation X vs. Y
  const NumericVector x_like = clone(prob_grid);
  const NumericVector y_like = 1.0 - prob_grid;

  // trial-wise belief updating
  for (int t = 0; t < n_trials; t++) {
    // compute discretised density of fixed Beta prior
    for (int i = 0; i < grid_size; i++) {
      DBM_prior[i] = R::dbeta(prob_grid[i], a[t], b[t], false);
    }
    normalise_inplace(DBM_prior);
    // compute predictive distribution
    if (t == 0) {
      // initialise predicted probability of observation X with fixed prior
      DBM_pred = clone(DBM_prior);
    } else {
      // update predictive distribution: mixture of previous trial's posterior
      // and fixed prior
      for (int i = 0; i < grid_size; i++) {
        DBM_pred[i] = (1.0 - cp[t]) * DBM_post[i] + cp[t] * DBM_prior[i];
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
    out = shannon_surprise(out, covariate);
  }
  return out;
}


// ----- Transition Probability Model -----------------------------------------

NumericVector run_tpm_nocp(
    const NumericVector &covariate,
    const NumericVector &a0,
    const NumericVector &b0,
    const bool return_surprise = false
) {
  const int n_trials = covariate.length();
  NumericVector out(n_trials);
  double n_hit_XX = 0.0;
  double n_trial_XX = 0.0;
  double n_hit_XY = 0.0;
  double n_trial_XY = 0.0;

  for (int t = 0; t < n_trials; t++) {
    const double a_XX = a0[t] + n_hit_XX;
    const double b_XX = b0[t] + (n_trial_XX - n_hit_XX);
    const double a_XY = a0[t] + n_hit_XY;
    const double b_XY = b0[t] + (n_trial_XY - n_hit_XY);

    const int prev = t > 0 ? static_cast<int>(covariate[(t - 1)]) : 1;
    const int curr = static_cast<int>(covariate[t]);

    out[t] = prev == 1 ? beta_mean(a_XX, b_XX) : beta_mean(a_XY, b_XY);

    if (prev == 1) {
      n_trial_XX += 1.0;
      if (curr == 1) {
        n_hit_XX += 1.0;
      }
    } else {
      n_trial_XY += 1.0;
      if (curr == 1) {
        n_hit_XY += 1.0;
      }
    }
  }

  if (return_surprise) {
    out = shannon_surprise(out, covariate);
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
      // likelihood of current obs = 1 given previous obs = 1
      out.like_XX[idx] = p_XX_val;
      // likelihood of current obs = 1 given previous obs = 0
      out.like_XY[idx] = p_XY_val;
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
    const NumericVector cp,
    const NumericVector a0,
    const NumericVector b0,
    const bool return_surprise = false
) {
  const int grid_res = 100;
  const double cp_eps = 1e-10;

  const int n_trials = covariate.length();
  for (int i = 0; i < n_trials; i++) {
    if (!(covariate[i] == 0.0 || covariate[i] == 1.0)) {
      stop("All `covariate` entries must be 0.0 or 1.0.");
    }
  }

  // if cp is (practically) equal to 0 (i.e., no volatility), fallback to simple
  // Beta-Binomial over transitions.
  bool cp_zero = true;
  for (int t = 0; t < n_trials; ++t) {
    if (cp[t] >= cp_eps) {
      cp_zero = false;
      break;
    }
  }
  if (cp_zero) {
    return run_tpm_nocp(covariate, a0, b0, return_surprise);
  }

  NumericVector out(n_trials);

  // if cp is (practically) equal to 1 (i.e., pure volatility), beliefs are given
  // by fixed prior
  bool cp_one = true;
  for (int t = 0; t < n_trials; ++t) {
    if ((1 - cp[t]) >= cp_eps) {
      cp_one = false;
      break;
    }
  }
  if (cp_one) {
    for (int t = 0; t < n_trials; ++t) {
      out[t] = beta_mean(a0[t], b0[t]);
    }
    if (return_surprise) {
      out = shannon_surprise(out, covariate);
    }
    return out;
  }

  // declare local variables
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
    TPM_post[j] = R::dbeta(grid.p_XX[j], a0[0], b0[0], false) *
      R::dbeta(grid.p_XY[j], a0[0], b0[0], false);
  }
  normalise_inplace(TPM_post);

  // loop over trials
  for (int t = 0; t < n_trials; t++) {

    const int curr = static_cast<int>(covariate[t]);
    const int prev = (t > 0) ? static_cast<int>(covariate[(t - 1)]) : NA_INTEGER;

    // Build predictive distribution from previous posterior.
    // - With prob 1-cp: environment is stable -> posterior propagates.
    // - With prob cp: change point -> new transition probs independent of old ones,
    //   so posterior mass is redistributed uniformly over all *other* grid points.
    double summed_post = std::accumulate(TPM_post.begin(), TPM_post.end(), 0.0);
    for (int j = 0; j < n_combi; j++) {
      TPM_pred[j] = (1.0 - cp[t]) * TPM_post[j] +
        cp[t] * (summed_post - TPM_post[j]) * inv_n_minus_1;
    }
    normalise_inplace(TPM_pred);

    // Before observing trial t, compute marginal predictive probability of
    // observing X - that is, predictive probability of covariate[t] == 1.0
    if (t == 0) {
      out[t] = mean_discrete(mean_p, TPM_pred);
    } else {
      out[t] = (prev == 1) ? mean_discrete(grid.p_XX, TPM_pred) : mean_discrete(grid.p_XY, TPM_pred);
    }

    // After observing trial t, update posterior distribution
    if (t > 0) {
      const NumericVector *like_ptr = nullptr;
      if (prev == 0) {
        if (curr == 0) {
          like_ptr = &grid.like_YY;
        } else {
          like_ptr = &grid.like_XY;
        }
      } else {
        if (curr == 0) {
          like_ptr = &grid.like_YX;
        } else {
          like_ptr = &grid.like_XX;
        }
      }
      for (int j = 0; j < n_combi; j++) {
        const double like_j = (*like_ptr)[j];
        TPM_update[j] = (1.0 - cp[t]) * like_j * TPM_post[j] +
          cp[t] * like_j * (summed_post - TPM_post[j]) * inv_n_minus_1;
      }
      normalise_inplace(TPM_update);
      std::swap(TPM_post, TPM_update);
    }
  }

  if (return_surprise) {
    out = shannon_surprise(out, covariate);
  }
  return out;
}

#endif
