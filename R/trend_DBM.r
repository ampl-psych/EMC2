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
  n_total <- length(covariate)
  if (n_trials < 1) {
    stop("`covariate` should consist of at least one observation.")
  }
  if (a0 <= 0 || b0 <= 0) {
    stop("Both prior shape parameters `a0` and `b0` must be positive.")
  }
  if (decay > 0 && window > 0) {
    stop("Cannot use both `decay` and `window`. Choose at most one memory constraint.")
  }

  out <- numeric(n_total)
  n_hit <- 0
  n_trial <- 0
  decay_factor <- exp(-1 / decay)
  event_memory <- numeric(0)

  for (t in seq_len(n_total)) {
    # prediction before observing trial t
    a_t <- a0 + n_hit
    b_t <- b0 + (n_trial - n_hit)
    if (return_map) {
      out[t] <- beta_mode(a_t, b_t)
    } else {
      out[t] <- beta_mean(a_t, b_t)
    }
    # update after observing trial t depends on decay / window memory constraints
    if (decay > 0 && decay < Inf) {
      # exponential decay
      n_hit <- decay_factor * (n_hit + covariate[t])
      n_trial <- decay_factor * (n_trial + 1)
    } else if (window > 0 && window < n_total) {
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
    return_map = FALSE, return_surprise = FALSE, grid_res = 100, cp_eps = 1e-12
) {
  n_total <- length(covariate)
  if (n_total < 1) {
    stop("`covariate` should consist of at least one observation.")
  }
  if (cp < 0 || cp > 1) {
    stop("Change point probability `cp` must be in the range [0, 1].")
  }
  if (mu0 <= 0 || mu0 >= 1) {
    stop("Prior mean `mu0` must be in the range (0, 1).")
  }
  if (s0 <= 0) {
    stop("Prior scale `s0` must be strictly positive.")
  }
  if (grid_res < 1) {
    stop("`grid_res` must be greater than or equal to 1.")
  }
  if (cp_eps <= 0) {
    stop("`cp_eps` must be strictly positive.")
  }

  a <- mu0 * s0
  b <- (1 - mu0) * s0

  # when cp = 0, DBM is actually fixed belief model a.k.a. beta-binomial
  if (cp < cp_eps) {
    return(run_beta_binomial(covariate, a, b, 0, 0, return_map, return_surprise))
  }

  out <- numeric(n_total)

  # when cp = 1, predictions are constant, determined purely by fixed prior
  if ((1 - cp) < cp_eps) {
    if (return_map) {
      out <- rep(beta_mode(a, b), n_total)
    } else {
      out <- rep(beta_mean(a, b), n_total)
    }
    if (return_surprise) {
      out <- -log2(out)
    }
    return(out)
  }

  # actual DBM
  prob_grid <- (0:grid_res) / grid_res
  DBM_prior <- normalise(stats::dbeta(prob_grid, a, b))
  x_like <- prob_grid
  y_like <- 1 - prob_grid
  DBM_pred <- DBM_prior

  for (t in seq_len(n_total)) {
    if (t > 1) {
      DBM_pred <- normalise((1 - cp) * DBM_post + cp * DBM_prior)
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

run_tpm <- function(
    covariate, cp, a0, b0,
    return_surprise = FALSE, grid_res = 100, cp_eps = 1e-12
) {
  n_total <- length(covariate)
  if (n_total < 1) {
    stop("`covariate` should consist of at least one observation.")
  }
  if (cp < 0 || cp > 1) {
    stop("Change point probability `cp` must be in the range [0, 1].")
  }
  if (a0 <= 0 || b0 <= 0) {
    stop("Both prior shape parameters `a0` and `b0` must be positive.")
  }
  if (grid_res < 1) {
    stop("`grid_res` must be greater than or equal to 1.")
  }
  if (cp_eps <= 0) {
    stop("`cp_eps` must be strictly positive.")
  }

  out <- numeric(n_total)

  # when cp = 0, fallback to simple beta-binomial over transitions
  if (cp < cp_eps) {
    a_XX <- a_XY <- a0
    b_XX <- b_XY <- b0
    out[1] <- beta_mean(a0, b0)
    for (t in 2:n_total) {
      prev <- covariate[(t - 1)]
      curr <- covariate[t]
      if (prev == 1) {
        out[t] <- beta_mean(a_XX, b_XX)
        if (curr == 1) {
          a_XX <- a_XX + 1
        } else {
          b_XX <- b_XX + 1
        }
      } else {
        out[t] <- beta_mean(a_XY, b_XY)
        if (curr == 1) {
          a_XY <- a_XY + 1
        } else {
          b_XY <- b_XY + 1
        }
      }
    }
    if (return_surprise) {
      out <- -log2(out)
    }
    return(out)
  }

  # when cp = 1, predictions are constant, determined purely by fixed prior
  if ((1 - cp) < cp_eps) {
    out <- rep(beta_mean(a0, b0), n_total)
    if (return_surprise) {
      out <- -log2(out)
    }
    return(out)
  }

  # actual TPM
  prob_grid <- (0:grid_res) / grid_res
  n_grid <- length(prob_grid)
  p_XX <- p_XY <- like_XX <- like_XY <- like_YY <- like_YX <- numeric(n_grid^2)
  idx <- 1L
  for (i0 in seq_len(n_grid)) {
    p_XY_val <- prob_grid[i0]
    for (i1 in seq_len(n_grid)) {
      p_XX_val <- prob_grid[i1]
      p_XX[idx] <- p_XX_val
      p_XY[idx] <- p_XY_val
      like_XX[idx] <- p_XX_val
      like_XY[idx] <- p_XY_val
      like_YY[idx] <- 1 - p_XY_val
      like_YX[idx] <- 1 - p_XX_val
      idx <- idx + 1L
    }
  }

  mean_p <- 0.5 * (p_XX + p_XY)
  TPM_post <- normalise(stats::dbeta(p_XX, a0, b0) * stats::dbeta(p_XY, a0, b0))
  inv_n_min_1 <- 1 / (length(TPM_post) - 1)

  for (t in seq_len(n_total)) {
    if (t > 1) {
      prev <- covariate[(t - 1)]
    }
    curr <- covariate[t]

    sum_TPM_post <- sum(TPM_post)

    TPM_pred <- normalise(
      (1 - cp) * TPM_post + cp * (sum_TPM_post - TPM_post) * inv_n_min_1
    )
    if (t == 1) {
      out[t] <- sum(mean_p * TPM_pred)
    } else {
      if (prev == 1) {
        out[t] <- sum(p_XX * TPM_pred)
      } else {
        out[t] <- sum(p_XY * TPM_pred)
      }
    }
    if (t > 1) {
      if (prev == 0) {
        if (curr == 0) {
          like_curr <- like_YY
        } else {
          like_curr <- like_XY
        }
      } else {
        if (curr == 0) {
          like_curr <- like_YX
        } else {
          like_curr <- like_XX
        }
      }
      TPM_post <- normalise(
        (1 - cp) * like_curr * TPM_post +
          cp * like_curr * (sum_TPM_post - TPM_post) * inv_n_min_1
      )
    }
  }

  if (return_surprise) {
    out <- -log2(out)
  }
  return(out)
}
