beta_mode <- function(a, b) {
  if (a == b) {
    return(0.5)
  }
  if (a < 1 || b < 1) {
    if (a < b) {
      return(0)
    }
    return(1)
  }
  return((a - 1) / (a + b - 2))
}

beta_mean <- function(a, b) {
  return(a / (a + b))
}

normalise <- function(x) {
  return(x / sum(x))
}

run_beta_binomial <- function(
    covariate, a0 = 1, b0 = 1, decay = 0, window = 0, return_map = FALSE, return_surprise = FALSE
) {
  if (a0 <= 0 || b0 <= 0) {
    stop("Both shape parameters `a0` and `b0` must be positive.")
  }
  if (decay > 0 && window > 0) {
    stop("Cannot use both `decay` and `window`. Choose at most one memory constraint.")
  }

  out <- numeric(length(covariate))
  n_hit <- 0
  n_trial <- 0
  decay_factor <- exp(-1 / decay)
  event_memory <- numeric(0)

  for (t in seq_along(covariate)) {
    # prediction before observing trial t
    a_t <- a0 + n_hit
    b_t <- b0 + (n_trial - n_hit)
    if (return_map) {
      out[t] <- beta_mode(a_t, b_t)
    } else {
      out[t] <- beta_mean(a_t, b_t)
    }
    # update after observing trial t depends on decay / window memory constraints
    if (decay > 0) {
      # exponential decay
      n_hit <- decay_factor * (n_hit + covariate[t])
      n_trial <- decay_factor * (n_trial + 1)
    } else if (window > 0) {
      # limited memory
      n_hit <- n_hit + covariate[t]
      n_trial <- n_trial + 1
      event_memory <- c(event_memory, covariate[t])
      if (length(event_memory) > window) {
        n_hit <- n_hit - event_memory[1]
        n_trial <- n_trial - 1
        event_memory <- event_memory[-1]
      }
    } else {
      # standard beta-binomial
      n_hit <- n_hit + covariate[t]
      n_trial <- n_trial + 1
    }
  }

  if (return_surprise) {
    out <- -log2(out)
  }
  return(out)
}

run_dbm <- function(
    covariate, cp, mu0, s0,
    return_map = FALSE, return_surprise = FALSE, grid_res = 1e3, cp_eps = 1e-10
) {
  if (cp < 0 || cp > 1) {
    stop("Change point probability `cp` must be in the range [0, 1].")
  }
  if (mu0 <= 0 || mu0 >= 1) {
    stop("Prior mean `mu0` must be in the range (0, 1).")
  }
  if (s0 <= 0) {
    stop("Prior scale `s0` must be strictly positive.")
  }
  a <- mu0 * s0
  b <- (1 - mu0) * s0
  # when cp = 1, DBM is actually fixed belief model a.k.a. beta-binomial
  if ((1 - cp) < cp_eps) {
    return(run_beta_binomial(covariate, a, b, 0, 0, return_map, return_surprise))
  }
  # when cp = 0, predictions are constant, determined purely by fixed prior
  if (cp < cp_eps) {
    if (return_map) {
      out <- rep(beta_mode(a, b), length(covariate))
    } else {
      out <- rep(beta_mean(a, b), length(covariate))
    }
    if (return_surprise) {
      out <- -log2(out)
    }
    return(out)
  }
  # actual DBM
  out <- numeric(length(covariate))
  prob_grid <- (0:grid_res) / grid_res
  DBM_prior <- normalise(stats::dbeta(prob_grid, a, b))
  x_like <- prob_grid
  y_like <- 1 - prob_grid
  DBM_pred <- DBM_prior
  DBM_post <- numeric(length(DBM_prior))
  for (t in seq_along(covariate)) {
    if (t > 1) {
      DBM_pred <- normalise(cp * DBM_post + (1 - cp) * DBM_prior)
    }
    if (return_map) {
      out[t] <- prob_grid[which.max(DBM_pred)]
    } else {
      out[t] <- sum(prob_grid * DBM_pred)
    }
    if (covariate[t] == 1) {
      DBM_post <- normalise(DBM_pred * x_like)
    } else {
      DBM_post <- normalise(DBM_pred * y_like)
    }
  }
  if (return_surprise) {
    out <- -log2(out)
  }
  return(out)
}
