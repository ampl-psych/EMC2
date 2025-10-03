# Wrapper functions for internal summary methods ------------------------------

#' R-hat convergence diagnostic
#'
#' Computes the potential scale reduction factor (R-hat) for each parameter
#' from an MCMC run. Supports both the legacy [coda::gelman.diag()] implementation
#' and the rank-normalized split-\eqn{\hat{R}} diagnostic of Vehtari et al. (2021).
#'
#' The `"new"` method is adapted directly from the \pkg{posterior} package's
#' [`rhat()`](https://mc-stan.org/posterior/reference/rhat.html) implementation,
#' but included here to avoid introducing an additional package dependency.
#' The substantive computation is identical, though some internal details
#' have been simplified.
#'
#' @param mcmc_list A `coda::mcmc.list` object containing MCMC draws.
#' @param version Character string, either `"old"` or `"new"`. `"old"` (default)
#'    calls `gelman_diag_robust()` (based on [coda::gelman.diag()]), while
#'    `"new"` uses a vendored implementation of `posterior::rhat()`.
#'
#' @return A named numeric vector of \eqn{\hat{R}} values, with names
#'   corresponding to parameters.
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: An improved \eqn{\hat{R}} for
#' assessing convergence of MCMC. *Bayesian Analysis*, *16*(2), 667–718.
#'
#' @keywords internal
r_hat <- function(mcmc_list, version = c("old", "new")) {
  version <- match.arg(version)
  if (version == "old") {
    result <- gelman_diag_robust(mcmc_list)
  } else {
    mcmc_mats <- prep_mcmc_diagnostics(mcmc_list)
    result <- vapply(
      X = mcmc_mats,
      FUN = function(x) {
        out <- try(rhat_new(x), silent = TRUE)
        if (is(out, "try-error")) {return(Inf)}
        return(out)
      },
      FUN.VALUE = numeric(1)
    )
  }
  return(result)
}

#' Effective sample size (ESS)
#'
#' Computes the effective sample size (ESS) for each parameter from an MCMC run.
#' Supports both the legacy [coda::effectiveSize()] implementation and the
#' improved version based on Vehtari et al. (2021).
#'
#' The `"new"` method is adapted directly from the \pkg{posterior} package's
#' [`ess_basic()`](https://mc-stan.org/posterior/reference/ess_basic.html)
#' implementation, but included here to avoid introducing an additional
#' package dependency. The substantive computation is identical, though some
#' internal details have been simplified.
#'
#' @param mcmc_list A `coda::mcmc.list` object containing MCMC draws.
#' @param version Character string, either `"old"` or `"new"`. `"old"` (default)
#'    calls [coda::effectiveSize()], while `"new"` uses a vendored implementation
#'    of `posterior::ess_basic()`.
#'
#' @return A named numeric vector of effective sample sizes, with names
#'   corresponding to parameters.
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: An improved \eqn{\hat{R}} for
#' assessing convergence of MCMC. *Bayesian Analysis*, *16*(2), 667–718.
#'
#' @keywords internal
n_eff <- function(mcmc_list, version = c("old", "new")) {
  version <- match.arg(version)
  if (version == "old") {
    result <- coda::effectiveSize(mcmc_list)
  } else {
    mcmc_mats <- prep_mcmc_diagnostics(mcmc_list)
    result <- vapply(
      X = mcmc_mats,
      FUN = function(x) {
        out <- try(ess_basic(x), silent = TRUE)
        if (is(out, "try-error")) {return(0)}
        return(out)
      },
      FUN.VALUE = numeric(1)
    )
  }
  return(result)
}

#' @noRd
prep_mcmc_diagnostics <- function(mcmc_list) {
  stopifnot(is(mcmc_list, "mcmc.list"))
  n_chains <- length(mcmc_list)
  n_iter <- unique(vapply(mcmc_list, nrow, integer(1)))
  if (length(n_iter) > 1L) {
    stop("Chains have unequal numbers of iterations; please trim or pad first.")
  }
  n_iter <- n_iter[1]
  param_names <- colnames(mcmc_list[[1]])
  result <- setNames(
    vector("list", length(param_names)),
    param_names
  )
  for (param in param_names) {
    mat <- matrix(NA_real_, nrow = n_iter, ncol = n_chains)
    for (chain in seq_len(n_chains)) {
      mat[ , chain] <- mcmc_list[[chain]][ , param]
    }
    result[[param]] <- mat
  }
  return(result)
}

# Old Rhat convergence diagnostic ---------------------------------------------

gelman_diag_robust <- function(
    mcl,
    autoburnin = FALSE,
    transform = TRUE,
    omit_mpsrf = TRUE
) {
  mcl <- split_mcl(mcl)
  gd <- try(
    coda::gelman.diag(
      mcl,
      autoburnin = autoburnin,
      transform = transform,
      multivariate = !omit_mpsrf
    ),
    silent = TRUE
  )
  if (is(gd, "try-error")) {
    if (omit_mpsrf) {
      return(list(psrf = matrix(Inf)))
    } else {
      return(list(psrf = matrix(Inf), mpsrf = Inf))
    }
  }
  gd_out <- gd[[1]][ , 1] # Remove CI
  if (!omit_mpsrf) {
    gd_out <- c(gd_out, gd[["mpsrf"]])
    names(gd_out)[length(gd_out)] <- "mpsrf"
  }
  return(gd_out)
}

split_mcl <- function(mcl) {
  if (!is.list(mcl)) {
    mcl <- list(mcl)
  }
  mcl2 <- mcl
  half <- floor(unlist(lapply(mcl, nrow))/2)
  for (i in 1:length(half)) {
    mcl[[i]] <- coda::as.mcmc(
      mcl[[i]][1:half[i], ]
    )
    mcl2[[i]] <- coda::as.mcmc(
      mcl2[[i]][c((half[i]+1):(2*half[i])), ]
    )
  }
  result <- coda::as.mcmc.list(c(mcl, mcl2))
  return(result)
}

# -----------------------------------------------------------------------------

# All of the following code was adapted from the `posterior` package
# (specifically `posterior/R/convergence.R`), corresponding to methods described
# in Vehtari et al. (2021), https://doi.org/10.1214/20-BA1221

# The `posterior` package is released under the BSD 3-Clause License. This
# is a permissive license that is compatible with the `EMC2` package license
# (GPL-3), as long as the original BSD license text is included; the original
# authors are attributed, and the "no endorsement" clause of the BSD license is
# respected.

# Original code copyright (C) 2012–2018 Trustees of Columbia University
# Copyright (C) 2018, 2019 Aki Vehtari, Paul Bürkner
#
# BSD 3-Clause License
#
# Copyright (c) 2021, posterior package authors;
# Stan Developers and their Assignees; Trustees of Columbia University
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# -----------------------------------------------------------------------------

# Rhat convergence diagnostic -------------------------------------------------

#' @noRd
rhat_new <- function(x) {
  rhat_bulk <- rhat_basic(rank_normalise(split_chains(x)))
  rhat_tail <- rhat_basic(rank_normalise(split_chains(fold_draws(x))))
  result <- max(rhat_bulk, rhat_tail)
  return(result)
}

#' @noRd
rhat_basic <- function(x) {
  if (has_bad_draws(x)) {
    return(NA_real_)
  }
  n_iter <- NROW(x)
  chain_mean <- colMeans(x)
  chain_var <- apply(x, 2, stats::var)
  var_between <- n_iter * stats::var(chain_mean)
  var_within <- mean(chain_var)
  result <- sqrt((var_between / var_within + n_iter - 1) / n_iter)
  return(result)
}


# Effective Sample Size (ESS) diagnostics -------------------------------------

#' @noRd
ess_basic <- function(x) {
  result <- ess(split_chains(x))
  return(result)
}

#' @noRd
ess_bulk <- function(x) {
  result <- ess(rank_normalise(split_chains(x)))
  return(result)
}

#' @noRd
ess_tail <- function(x) {
  ess_lower <- ess_quantile(x, probs = 0.05)
  ess_upper <- ess_quantile(x, probs = 0.95)
  result <- min(ess_lower, ess_upper)
  return(result)
}

#' @noRd
ess_mean <- function(x) {
  result <- ess(split_chains(x))
  return(result)
}

#' @noRd
ess_median <- function(x) {
  result <- ess_quantile(x, probs = 0.5, names = FALSE)
  return(result)
}

#' @noRd
ess_sd <- function(x) {
  result <- ess(split_chains((x - mean(x))^2))
  return(result)
}

#' @noRd
ess_quantile <- function(x, probs, names = TRUE) {
  probs <- check_quantile_probs(probs)
  result <- unlist(lapply(probs, ess_quantile_engine, x = x))
  if (names) {
    names(result) <- paste0("ess_q", probs * 100)
  }
  return(result)
}

#' @noRd
ess_quantile_engine <- function(x, prob) {
  if (has_bad_draws(x)) {
    return(NA_real_)
  }
  if (prob == 1) {
    prob <- (length(x) - 0.5) / length(x)
  }
  idx <- x <= stats::quantile(x, prob)
  result <- ess(split_chains(idx))
  return(result)
}

#' @noRd
ess <- function(x) {

  n_chain <- NCOL(x)
  n_iter <- NROW(x)
  if (n_iter < 3L || has_bad_draws(x)) {
    return(NA_real_)
  }
  n_samples <- n_chain * n_iter

  acov <- apply(x, 2, autocovariance)
  acov_means <- rowMeans(acov)
  mean_var <- acov_means[1] * n_iter / (n_iter - 1)
  var_plus <- mean_var * (n_iter - 1) / n_iter
  if (n_chain > 1) {
    var_plus <- var_plus + stats::var(colMeans(x))
  }

  # Geyer's initial positive sequence
  rho_hat_t <- rep.int(0, n_iter)
  t <- 0
  rho_hat_even <- 1
  rho_hat_t[t + 1] <- rho_hat_even
  rho_hat_odd <- 1 - (mean_var - acov_means[t + 2]) / var_plus
  rho_hat_t[t + 2] <- rho_hat_odd
  while (
    (t < NROW(acov) - 5) && !is.nan(rho_hat_even + rho_hat_odd) &&
    (rho_hat_even + rho_hat_odd > 0)
  ) {
    t <- t + 2
    rho_hat_even <- 1 - (mean_var - acov_means[t + 1]) / var_plus
    rho_hat_odd <- 1 - (mean_var - acov_means[t + 2]) / var_plus
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_t[t + 1] <- rho_hat_even
      rho_hat_t[t + 2] <- rho_hat_odd
    }
  }
  max_t <- t # used in the improved estimate
  if (rho_hat_even > 0) {
    rho_hat_t[max_t + 1] <- rho_hat_even
  }

  # Geyer's initial monotone sequence
  t <- 0
  while (t <= (max_t - 4)) {
    t <- t + 2
    if (
      (rho_hat_t[t + 1] + rho_hat_t[t + 2]) > (rho_hat_t[t - 1] + rho_hat_t[t])
    ) {
      rho_hat_t[t + 1] <- (rho_hat_t[t - 1] + rho_hat_t[t]) / 2
      rho_hat_t[t + 2] <- rho_hat_t[t + 1]
    }
  }

  # Geyer's truncated estimate
  # tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t])
  # Improved estimate reduces variance in antithetic case
  tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t]) + rho_hat_t[max_t+1]
  # safety check for negative values and with max ess equal to
  # n_samples*log10(n_samples)
  tau_bound <- 1 / log10(n_samples)
  if (tau_hat < tau_bound) {
    warning("The ESS has been capped to avoid unstable estimates.")
    tau_hat <- tau_bound
  }
  result <- n_samples / tau_hat
  return(result)
}


# Monte Carlo Standard Error (MCSE) diagnostics -------------------------------

#' @noRd
mcse_mean <- function(x) {
  result <- stats::sd(x) / sqrt(ess_mean(x))
  return(result)
}

#' @noRd
mcse_median <- function(x) {
  result <- mcse_quantile(x, probs = 0.5, names = FALSE)
  return(result)
}

#' @noRd
mcse_sd <- function(x) {
  x_c <- x - mean(x)
  ess_x <- ess_mean((x_c)^2)
  # Variance of variance estimate by Kenney and Keeping (1951, p. 141),
  # which doesn't assume normality of sims.
  E_var <- mean(x_c^2)
  var_var <- (mean(x_c^4) - E_var^2) / ess_x
  # The first order Taylor series approximation of variance of sd
  var_sd <- var_var / E_var / 4
  result <- sqrt(var_sd)
  return(result)
}

#' @noRd
mcse_quantile <- function(x, probs, names = TRUE) {
  probs <- check_quantile_probs(probs)
  result <- unlist(lapply(probs, mcse_quantile_engine, x = x))
  if (names) {
    names(result) <- paste0("mcse_q", probs * 100)
  }
  return(result)
}

#' @noRd
mcse_quantile_engine <- function(x, prob) {
  ess_q <- ess_quantile(x, prob)
  p <- c(0.1586553, 0.8413447)
  a <- stats::qbeta(p, ess_q * prob + 1, ess_q * (1 - prob) + 1)
  ssims <- sort(x)
  S <- length(ssims)
  th1 <- ssims[max(floor(a[1] * S), 1)]
  th2 <- ssims[min(ceiling(a[2] * S), S)]
  result <- as.vector((th2 - th1) / 2)
  return(result)
}


# Helper functions ------------------------------------------------------------

#' @noRd
rank_normalise <- function(x) {
  # replace values by ranks, using average rank for ties to maintain the number
  # of unique values of discrete quantities
  r <- rank(as.array(x), ties.method = "average")
  # transform ranks into percentiles, using fractional offset of 3/8 as
  # recommended by Blom (1958)
  c <- 3 / 8
  p <- (r - c) / (length(r) - 2 * c + 1)
  # normalise via inverse normal CDF
  z <- stats::qnorm(p)
  # preserve NA's and original shape
  z[is.na(x)] <- NA_real_
  if (!is.null(dim(x))) {
    z <- array(z, dim = dim(x), dimnames = dimnames(x))
  }
  return(z)
}

#' @noRd
split_chains <- function(x) {
  niter <- NROW(x)
  if (niter == 1L) {
    return(x)
  }
  half <- niter / 2
  result <- cbind(x[1:floor(half), ], x[ceiling(half + 1):niter, ])
  return(result)
}

#' @noRd
fold_draws <- function(x) {
  result <- abs(x - median(x))
  return(result)
}

#' @noRd
has_bad_draws <- function(x, tol = .Machine$double.eps) {
  nonfinite_vals <- anyNA(x) || any(is.infinite(x))
  constant_vals <- abs(max(x) - min(x)) < tol
  return(nonfinite_vals || constant_vals)
}

#' @noRd
autocovariance <- function(x) {
  N <- length(x)
  var_x <- stats::var(x)
  if (var_x == 0) {
    return(rep(0, N))
  }
  M <- stats::nextn(N)
  M_double <- 2 * M
  x_c <- x - mean(x)
  x_c <- c(x_c, rep.int(0, M_double - N))
  # FFT-based unnormalised autocovariances
  ac <- Re(
    stats::fft(
      abs(stats::fft(x_c))^2,
      inverse = TRUE
    )[1:N]
  )
  # use "biased" estimate as recommended by Geyer (1992)
  # direct scaling with var(x) avoids need to compute "mask effect"
  result <- ac / ac[1]
  result <- result * var_x * (N - 1) / N
  return(result)
}

#' @noRd
check_quantile_probs <- function(probs) {
  probs <- as.numeric(probs)
  if (any(probs < 0 | probs > 1)) {
    stop("'probs' must contain values between 0 and 1.")
  }
  return(probs)
}
