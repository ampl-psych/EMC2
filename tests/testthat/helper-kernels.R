# Beta distribution helpers ---------------------------------------------------

beta_mean <- function(a, b) {
  return(a / (a + b))
}

beta_mode <- function(a, b) {
  out <- numeric(length(a))
  # symmetric case: a == b.
  # - if a = b < 1: density diverges at both 0 and 1: no unique mode; return 0.5
  # - if a = b = 1: uniform(0, 1): no unique mode; return 0.5 (i.e., the mean)
  # - if a = b > 1: unimodal symmetric distribution: mode is exactly 0.5
  idx_equal <- a == b
  out[idx_equal] <- 0.5
  # if either a or b is less than 1, the density diverges at one endpoint.
  # divergence is stronger on the side with the smaller parameter.
  idx_boundary <- (a < 1 | b < 1) & !idx_equal
  idx_left <- idx_boundary & (a < b)
  idx_right <- idx_boundary & !idx_left
  out[idx_left] <- 0
  out[idx_right] <- 1
  # typical case: both a and b > 1
  idx_regular <- !(idx_equal | idx_boundary)
  out[idx_regular] <- (a[idx_regular] - 1) /
    (a[idx_regular] + b[idx_regular] - 2)
  return(out)
}

beta_logprecision <- function(a, b) {
  return(2 * log(a + b) + log1p(a + b) - log(a) - log(b))
}

# Discretised probability distribution helpers --------------------------------

normalise <- function(x) {
  return(x / sum(x))
}

mean_discrete <- function(grid, prob) {
  return(sum(grid * prob))
}

mode_discrete <- function(grid, prob) {
  return(grid[which.max(prob)])
}

logprecision_discrete <- function(grid, prob) {
  ceiling <- -log(.Machine$double.eps)
  m <- sum(grid * prob)
  var <- max(0, sum((grid - m)^2 * prob))
  out <- -log(var)
  if (!is.finite(out) || out > ceiling) {
    return(ceiling)
  }
  return(out)
}

# Shannon surprise (bits) -----------------------------------------------------

shannon_surprise <- function(pred, obs) {
  pred_safe <- pmin(pmax(pred, .Machine$double.eps), (1 - .Machine$double.neg.eps))
  return(-log2(ifelse(obs == 1, pred_safe, (1 - pred_safe))))
}


# Beta-binomial model ---------------------------------------------------------

# beta-binomial model, optionally with exponential decay ("leaky integration")
# applied to past observations, or with memory of past events limited to a
# fixed window
beta_binomial <- function(
    x, a0, b0,
    decay = 0, window = 0,
    output_type = c("mean", "mode", "surprise", "log-precision")
) {
  if (!all(x %in% c(0, 1) | is.na(x))) {
    stop("All `x` entries that are not NA must be 0 or 1.")
  }
  output_type <- match.arg(output_type)

  n_total <- length(x)
  pars <- matrix(nrow = 4, ncol = n_total)
  rownames(pars) <- c("a0", "b0", "decay", "window")
  pars["a0", ] <- as.numeric(a0)
  pars["b0", ] <- as.numeric(b0)
  pars["decay", ] <- as.numeric(decay)
  pars["window", ] <- as.numeric(window)

  use_decay <- any(pars["decay", ] > 0)
  use_window <- any(pars["window", ] > 0)
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
      while (length(buf[["obs"]]) > 0 && (t - buf[["idx"]][1]) > pars["window", t]) {
        n_hit <- n_hit - buf[["obs"]][1]
        n_trial <- n_trial - 1
        buf[["obs"]] <- buf[["obs"]][-1]
        buf[["idx"]] <- buf[["idx"]][-1]
      }
    }
    # prediction before observing trial t:
    a_t <- pars["a0", t] + n_hit
    b_t <- pars["b0", t] + (n_trial - n_hit)
    if (output_type == "mode") {
      out[t] <- beta_mode(a_t, b_t)
    } else if (output_type == "log-precision") {
      out[t] <- beta_logprecision(a_t, b_t)
    } else {
      out[t] <- beta_mean(a_t, b_t)
    }
    # if current trial is NA, special handling of update:
    if (is.na(x[t])) {
      # since exponential decay can be thought of as time-based forgetting,
      # we still apply decay to n_hit and n_trial to account for passage of time
      if (use_decay) {
        n_hit <- exp(-1 / pars["decay", t]) * n_hit
        n_trial <- exp(-1 / pars["decay", t]) * n_trial
      }
      # since sliding window is modelling memory capacity for actual observed
      # evidence, not time-based decay, we do not update anything in case of
      # missing observation
      # Same for basic Beta-Binomial model
      next
    }
    # update after observing trial t:
    # depends on decay / window memory constraints
    if (use_decay) {
      # exponential decay
      n_hit <- exp(-1 / pars["decay", t]) * (n_hit + x[t])
      n_trial <- exp(-1 / pars["decay", t]) * (n_trial + 1)
    } else if (use_window) {
      # add to limited memory
      buf[["obs"]] <- c(buf[["obs"]], x[t])
      buf[["idx"]] <- c(buf[["idx"]], t)
      n_hit <- n_hit + x[t]
      n_trial <- n_trial + 1
    } else {
      # standard beta-binomial
      n_hit <- n_hit + x[t]
      n_trial <- n_trial + 1
    }
  }

  if (output_type == "surprise") {
    out <- shannon_surprise(out, x)
  }
  return(out)
}


# Dynamic belief model --------------------------------------------------------

# Dynamic Belief Model (DBM), based on Yu & Cohen (2008), NeurIPS; and
# Ide et al. (2013), JoN
# change point probability (cp) determines the extent to which trial-wise updates
# are "reset" from a fixed prior distribution, which is a Beta distribution
# parameterised with mean (mu0) and scale (s0)
# If cp is practically equal to 0 (i.e., no volatility), the DBM is referred to
# the Fixed Belief Model, which is equivalent to a standard Beta-binomial model.
# If cp is practically equal to 1 (i.e., pure volatility), there is no actual
# learning from input x; the beliefs are purely driven by the fixed prior, hence
# the trial-wise output is constant.
dbm <- function(
    x, cp, mu0, s0,
    output_type = c("mean", "mode", "surprise", "log-precision"),
    grid_res = 100L
) {
  if (!all(x %in% c(0, 1) | is.na(x))) {
    stop("All `x` entries that are not NA must be 0 or 1.")
  }
  output_type <- match.arg(output_type)

  n_total <- length(x)
  pars <- matrix(nrow = 3, ncol = n_total)
  rownames(pars) <- c("cp", "a", "b")
  pars["cp", ] <- as.numeric(cp)
  # shape parameters of Beta prior
  # NB in original Matlab code from Jaime Ide, the parameterisation
  # a = mu0 * s0 + 1; b = (1 - mu0) * s0 + 1
  # was used. The +1 shift was presumably a pragmatic tweak to avoid the shape
  # parameters ever being <= 1.
  pars["a", ] <- as.numeric(mu0) * as.numeric(s0)
  pars["b", ] <- (1 - as.numeric(mu0)) * as.numeric(s0)
  out <- numeric(n_total)

  # discretised density grids
  prob_grid <- (0:grid_res) / grid_res
  # pre-compute Bernoulli likelihoods for binary observation X vs. Y
  x_like <- prob_grid
  y_like <- 1 - prob_grid

  for (t in seq_len(n_total)) {
    # compute discretised Beta prior for trial t
    DBM_prior <- normalise(stats::dbeta(prob_grid, pars["a", t], pars["b", t]))
    # compute predictive distribution:
    if (t == 1) {
      # initialise predicted probability of observation X with fixed prior
      DBM_pred <- DBM_prior
    } else {
      # update predictive distribution: mixture of previous trial's posterior
      # and fixed prior
      DBM_pred <- normalise(
        (1 - pars["cp", t]) * DBM_post + pars["cp", t] * DBM_prior
      )
    }
    # main trial-wise output: mean, mode, or log-precision of predictive distribution
    # for the probability of observation X
    if (output_type == "mode") {
      out[t] <- mode_discrete(prob_grid, DBM_pred)
    } else if (output_type == "log-precision") {
      out[t] <- logprecision_discrete(prob_grid, DBM_pred)
    } else {
      out[t] <- mean_discrete(prob_grid, DBM_pred)
    }
    # update posterior distribution
    if (is.na(x[t])) {
      # NA handling: predictive distribution pushed forward as "posterior" for
      # next trial, effectively causing slight forgetting of past observations
      # (see Zhou et al. 2020, Cogsci)
      DBM_post <- DBM_pred
    } else if (x[t] == 1) {
      DBM_post <- normalise(DBM_pred * x_like)
    } else {
      DBM_post <- normalise(DBM_pred * y_like)
    }
  }

  if (output_type == "surprise") {
    out <- shannon_surprise(out, x)
  }
  return(out)
}


# TODO ------------------------------------------------------------------------

# run_tpm_nocp <- function(
#     covariate, a0, b0, return_surprise
# ) {
#   n_total <- length(covariate)
#   out <- numeric(n_total)
#   n_hit_XX <- n_trial_XX <- n_hit_XY <- n_trial_XY <- 0
#   for (t in seq_len(n_total)) {
#     a_XX <- a0[t] + n_hit_XX
#     b_XX <- b0[t] + (n_trial_XX - n_hit_XX)
#     a_XY <- a0[t] + n_hit_XY
#     b_XY <- b0[t] + (n_trial_XY - n_hit_XY)
#     if (t == 1) {
#       out[t] <- beta_mean(a_XX, b_XX)
#       next
#     }
#     prev <- covariate[(t - 1)]
#     curr <- covariate[t]
#     if (is.na(prev)) {
#       out[t] <- 0.5 * (beta_mean(a_XX, b_XX) + beta_mean(a_XY, b_XY))
#     } else {
#       if (prev == 1) {
#         out[t] <- beta_mean(a_XX, b_XX)
#       } else {
#         out[t] <- beta_mean(a_XY, b_XY)
#       }
#     }
#     if (is.na(curr) || is.na(prev)) {
#       next
#     }
#     if (prev == 1) {
#       n_trial_XX <- n_trial_XX + 1
#       if (curr == 1) {
#         n_hit_XX <- n_hit_XX + 1
#       }
#     } else {
#       n_trial_XY <- n_trial_XY + 1
#       if (curr == 1) {
#         n_hit_XY <- n_hit_XY + 1
#       }
#     }
#   }
#   if (return_surprise) {
#     out <- shannon_surprise(out, covariate)
#   }
#   return(out)
# }

# run_tpm <- function(
#     covariate, cp, a0, b0, return_surprise = FALSE
# ) {
#   if (!all(covariate %in% c(0, 1) | is.na(covariate))) {
#     stop("All `covariate` entries that are not NA must be 0 or 1.")
#   }
#
#   cp_eps <- 1e-10
#
#   # when cp = 0, fallback to simple beta-binomial over transitions
#   if (all(cp < cp_eps)) {
#     return(run_tpm_nocp(covariate, a0, b0, return_surprise))
#   }
#
#   # when cp = 1, predictions are constant, determined purely by fixed prior
#   if (all((1 - cp) < cp_eps)) {
#     out <- beta_mean(a0, b0)
#     if (return_surprise) {
#       out <- shannon_surprise(out, covariate)
#     }
#     return(out)
#   }
#
#   # actual TPM
#   grid_res <- 100
#   n_total <- length(covariate)
#   out <- numeric(n_total)
#   prob_grid <- (0:grid_res) / grid_res
#   n_grid <- length(prob_grid)
#   p_XX <- p_XY <- like_XX <- like_XY <- like_YY <- like_YX <- numeric(n_grid^2)
#   idx <- 1L
#   for (i0 in seq_len(n_grid)) {
#     p_XY_val <- prob_grid[i0]
#     for (i1 in seq_len(n_grid)) {
#       p_XX_val <- prob_grid[i1]
#       p_XX[idx] <- p_XX_val
#       p_XY[idx] <- p_XY_val
#       like_XX[idx] <- p_XX_val
#       like_XY[idx] <- p_XY_val
#       like_YY[idx] <- 1 - p_XY_val
#       like_YX[idx] <- 1 - p_XX_val
#       idx <- idx + 1L
#     }
#   }
#
#   mean_p <- 0.5 * (p_XX + p_XY)
#   TPM_post <- normalise(
#     stats::dbeta(p_XX, a0[1], b0[1]) * stats::dbeta(p_XY, a0[1], b0[1])
#   )
#   inv_n_min_1 <- 1 / (length(TPM_post) - 1)
#
#   prev <- NA
#
#   for (t in seq_len(n_total)) {
#     if (t == 1) {
#       prev <- NA
#     } else {
#       prev <- covariate[(t - 1)]
#     }
#     curr <- covariate[t]
#
#     sum_TPM_post <- sum(TPM_post)
#
#     TPM_pred <- normalise(
#       (1 - cp[t]) * TPM_post + cp[t] * (sum_TPM_post - TPM_post) * inv_n_min_1
#     )
#     if (is.na(prev)) {
#       out[t] <- sum(mean_p * TPM_pred)
#     } else {
#       if (prev == 1) {
#         out[t] <- sum(p_XX * TPM_pred)
#       } else {
#         out[t] <- sum(p_XY * TPM_pred)
#       }
#     }
#
#     if (is.na(curr) || is.na(prev)) {
#       TPM_post <- TPM_pred
#       next
#     }
#
#     if (prev == 0) {
#       if (curr == 0) {
#         like_curr <- like_YY
#       } else {
#         like_curr <- like_XY
#       }
#     } else {
#       if (curr == 0) {
#         like_curr <- like_YX
#       } else {
#         like_curr <- like_XX
#       }
#     }
#     TPM_post <- normalise(
#       (1 - cp[t]) * like_curr * TPM_post +
#         cp[t] * like_curr * (sum_TPM_post - TPM_post) * inv_n_min_1
#     )
#
#   }
#
#   if (return_surprise) {
#     out <- shannon_surprise(out, covariate)
#   }
#   return(out)
# }
