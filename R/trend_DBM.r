beta_mode <- function(a, b) {
  out <- numeric(length(a))
  # case 1: a == b -> mode = 0.5
  idx_equal <- a == b
  out[idx_equal] <- 0.5
  # case 2: a < 1 or b < 1
  idx_boundary <- (a < 1 | b < 1) & !idx_equal
  #   within that: if a < b -> mode = 0, else -> mode = 1
  idx_left <- idx_boundary & (a < b)
  idx_right <- idx_boundary & !idx_left
  out[idx_left] <- 0
  out[idx_right] <- 1
  # case 3: mode = (a - 1) / (a + b - 2)
  idx_regular <- !(idx_equal | idx_boundary)
  out[idx_regular] <- (a[idx_regular] - 1) /
    (a[idx_regular] + b[idx_regular] - 2)
  return(out)
}

beta_mean <- function(a, b) {
  return(a / (a + b))
}

normalise <- function(x) {
  return(x / sum(x))
}

shannon_surprise <- function(pred, obs) {
  pred_safe <- pmin(pmax(pred, .Machine$double.eps), (1 - .Machine$double.neg.eps))
  return(-log2(ifelse(obs == 1, pred_safe, (1 - pred_safe))))
}

run_beta_binomial <- function(
    covariate, a0, b0, decay, window, return_map = FALSE, return_surprise = FALSE
) {
  n_total <- length(covariate)
  if (!all(covariate) %in% c(0, 1)) {
    stop("All `covariate` entries must be 0 or 1.")
  }
  use_decay <- any(decay > 0)
  use_window <- any(window > 0)
  if (use_decay && use_window) {
    stop("Cannot use both `decay` and `window`. Choose only one memory constraint.")
  }

  out <- numeric(n_total)
  n_hit <- 0
  n_trial <- 0
  buf <- list(obs = numeric(0), idx = integer(0))

  for (t in seq_len(n_total)) {
    # if applicable, prune memory based on current memory window
    if (use_window) {
      while (length(buf[["obs"]]) > 0 && (t - buf[["idx"]][1]) >= window[t]) {
        n_hit <- n_hit - buf[["obs"]][1]
        n_trial <- n_trial - 1
        buf[["obs"]] <- buf[["obs"]][-1]
        buf[["idx"]] <- buf[["idx"]][-1]
      }
    }
    # prediction before observing trial t
    a_t <- a0[t] + n_hit
    b_t <- b0[t] + (n_trial - n_hit)
    if (return_map) {
      out[t] <- beta_mode(a_t, b_t)
    } else {
      out[t] <- beta_mean(a_t, b_t)
    }
    # if current trial is NA, there is nothing to learn, so move on
    if (is.na(covariate[t])) {
      next
    }
    # update after observing trial t depends on decay / window memory constraints
    if (use_decay) {
      # exponential decay
      n_hit <- exp(-1 / decay[t]) * (n_hit + covariate[t])
      n_trial <- exp(-1 / decay[t]) * (n_trial + 1)
    } else if (use_window) {
      # add to limited memory
      buf[["obs"]] <- c(buf[["obs"]], covariate[t])
      buf[["idx"]] <- c(buf[["idx"]], t)
      n_hit <- n_hit + covariate[t]
      n_trial <- n_trial + 1
    } else {
      # standard beta-binomial
      n_hit <- n_hit + covariate[t]
      n_trial <- n_trial + 1
    }
  }

  if (return_surprise) {
    out <- shannon_surprise(out, covariate)
  }
  return(out)
}

run_dbm <- function(
    covariate, cp, mu0, s0, return_map = FALSE, return_surprise = FALSE
) {
  if (return_map) {
    grid_res <- 500
  } else {
    grid_res <- 100
  }
  cp_eps <- 1e-10

  n_total <- length(covariate)
  if (!all(covariate) %in% c(0, 1)) {
    stop("All `covariate` entries must be 0 or 1.")
  }

  a <- mu0 * s0
  b <- (1 - mu0) * s0

  # when cp = 0, DBM is actually fixed belief model a.k.a. beta-binomial
  if (all(cp < cp_eps)) {
    return(run_beta_binomial(covariate, a, b, 0, 0, return_map, return_surprise))
  }

  # when cp = 1, predictions are constant, determined purely by fixed prior
  if (all((1 - cp) < cp_eps)) {
    if (return_map) {
      out <- beta_mode(a, b)
    } else {
      out <- beta_mean(a, b)
    }
    if (return_surprise) {
      out <- shannon_surprise(out, covariate)
    }
    return(out)
  }

  # actual DBM
  out <- numeric(n_total)
  prob_grid <- (0:grid_res) / grid_res
  x_like <- prob_grid
  y_like <- 1 - prob_grid

  for (t in seq_len(n_total)) {
    DBM_prior <- normalise(stats::dbeta(prob_grid, a[t], b[t]))
    if (t == 1) {
      DBM_pred <- DBM_prior
    } else {
      DBM_pred <- normalise((1 - cp[t]) * DBM_post + cp[t] * DBM_prior)
    }
    if (return_map) {
      out[t] <- prob_grid[which.max(DBM_pred)]
    } else {
      out[t] <- sum(prob_grid * DBM_pred)
    }
    if (is.na(covariate[t])) {
      # perhaps controversial, needs more thought
      DBM_post <- DBM_pred
    } else if (covariate[t] == 1) {
      DBM_post <- normalise(DBM_pred * x_like)
    } else {
      DBM_post <- normalise(DBM_pred * y_like)
    }
  }

  if (return_surprise) {
    out <- shannon_surprise(out, covariate)
  }
  return(out)
}

run_tpm_nocp <- function(
    covariate, a0, b0, return_surprise
) {
  n_total <- length(covariate)
  n_hit_XX <- n_trial_XX <- n_hit_XY <- n_trial_XY <- 0
  for (t in seq_len(n_total)) {
    a_XX <- a0[t] + n_hit_XX
    b_XX <- b0[t] + (n_trial_XX - n_hit_XX)
    a_XY <- a0[t] + n_hit_XY
    b_XY <- b0[t] + (n_trial_XY - n_hit_XY)
    if (t == 1) {
      out[t] <- beta_mean(a_XX, b_XX)
      next
    }
    prev <- covariate[(t - 1)]
    curr <- covariate[t]
    if (is.na(prev)) {
      out[t] <- out[(t - 1)]
      next
    }
    if (prev == 1) {
      out[t] <- beta_mean(a_XX, b_XX)
    } else {
      out[t] <- beta_mean(a_XY, b_XY)
    }
    if (is.na(curr)) {
      next
    }
    if (prev == 1) {
      n_trial_XX <- n_trial_XX + 1
      if (curr == 1) {
        n_hit_XX <- n_hit_XX + 1
      }
    } else {
      n_trial_XY <- n_trial_XY + 1
      if (curr == 1) {
        n_hit_XY <- n_hit_XY + 1
      }
    }
  }
  if (return_surprise) {
    out <- shannon_surprise(out, covariate)
  }
  return(out)
}

run_tpm <- function(
    covariate, cp, a0, b0, return_surprise = FALSE
) {
  grid_res <- 100
  cp_eps <- 1e-10

  n_total <- length(covariate)
  if (!all(covariate) %in% c(0, 1)) {
    stop("All `covariate` entries must be 0 or 1.")
  }

  out <- numeric(n_total)

  # when cp = 0, fallback to simple beta-binomial over transitions
  if (all(cp < cp_eps)) {
    return(run_tpm_nocp(covariate, a0, b0, return_surprise))
  }

  # when cp = 1, predictions are constant, determined purely by fixed prior
  if (all((1 - cp) < cp_eps)) {
    out <- beta_mean(a0, b0)
    if (return_surprise) {
      out <- shannon_surprise(out, covariate)
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
  TPM_post <- normalise(
    stats::dbeta(p_XX, a0[1], b0[1]) * stats::dbeta(p_XY, a0[1], b0[1])
  )
  inv_n_min_1 <- 1 / (length(TPM_post) - 1)

  for (t in seq_len(n_total)) {
    if (t > 1) {
      prev <- covariate[(t - 1)]
    }
    curr <- covariate[t]

    sum_TPM_post <- sum(TPM_post)

    TPM_pred <- normalise(
      (1 - cp[t]) * TPM_post + cp[t] * (sum_TPM_post - TPM_post) * inv_n_min_1
    )
    if (t == 1) {
      out[t] <- sum(mean_p * TPM_pred)
    } else {
      if (is.na(prev)) {
        out[t] <- out[(t - 1)]
        next
      }
      if (prev == 1) {
        out[t] <- sum(p_XX * TPM_pred)
      } else {
        out[t] <- sum(p_XY * TPM_pred)
      }
    }
    if (t > 1) {
      if (is.na(curr)) {
        # perhaps controversial, needs more discussion
        TPM_post <- TPM_pred
        next
      }
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
        (1 - cp[t]) * like_curr * TPM_post +
          cp[t] * like_curr * (sum_TPM_post - TPM_post) * inv_n_min_1
      )
    }
  }

  if (return_surprise) {
    out <- shannon_surprise(out, covariate)
  }
  return(out)
}
