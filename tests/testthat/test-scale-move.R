# The interweaving (ASIS) scale move of the standard sampler
# (scale_move_standard, R/variant_standard.R): after the group step it
# rescales a group SD and every subject's deviation from the group mean
# together, which is the non-centred step that breaks the hierarchical
# funnel's autocorrelation (rating-work/sampler/hier/stageH2/REPORT.md).

# ---- (1) the acceptance ratio's prior part, against brute force -----------
# log r = sum_s dl_s + log_prior(delta), with log_prior the change in the
# group density, the covariance prior (with its auxiliary a_j), the prior of
# a_j and the Jacobian of the map
# (alpha_.j, Sigma, a_j) -> (mu_j + f (alpha_.j - mu_j), D Sigma D, a_j / f^2),
# for a parameter in a covariance block and for one with its own variance.
test_that("scale_move_log_prior matches the brute-force ratio", {
  set.seed(5)
  v <- 2; A <- c(.3, .5, 1)
  brute <- function(alpha, mu, Sigma, a, j, delta, blocked) {
    n <- nrow(alpha); p <- ncol(alpha); f <- exp(delta)
    D <- rep(1, p); D[j] <- f
    alpha2 <- alpha; alpha2[, j] <- mu[j] + f * (alpha[, j] - mu[j])
    Sigma2 <- Sigma * outer(D, D)
    a2 <- a; a2[j] <- a[j] / f^2
    lprior <- function(S, a) (if (blocked) {
      log(EMC2:::robust_diwish(S, v + p - 1, 2 * v * diag(1 / a, p)))
    } else sum(EMC2:::logdinvGamma(diag(S), v / 2, v / a))) +
      sum(EMC2:::logdinvGamma(a, 1 / 2, 1 / A^2))
    group <- function(al, S) sum(mvtnorm::dmvnorm(al, mu, S, log = TRUE))
    # Jacobian: n deviations, the free elements of the covariance block, a_j
    log_jac <- n * delta + (if (blocked) (p + 1) * delta else 2 * delta) - 2 * delta
    lprior(Sigma2, a2) - lprior(Sigma, a) + group(alpha2, Sigma2) - group(alpha, Sigma) + log_jac
  }
  for (blocked in c(TRUE, FALSE)) {
    p <- 3; n <- 4
    Sigma <- if (blocked) crossprod(matrix(rnorm(9), 3)) + diag(3) * .2 else diag(exp(rnorm(3)))
    a <- 1 / rgamma(p, 1 / 2, rate = 1 / A^2)
    mu <- rnorm(p); alpha <- mvtnorm::rmvnorm(n, mu, Sigma)
    for (j in 1:p) for (delta in c(-.7, .2, 1.1)) {
      got <- EMC2:::scale_move_log_prior(delta, A[j], a[j])
      expect_equal(got, brute(alpha, mu, Sigma, a, j, delta, blocked), tolerance = 1e-8)
    }
  }
})

# ---- (2) detailed balance on a 2-subject toy -------------------------------
# Two subjects, two parameters, a heavy-tailed (non-Gaussian) likelihood on
# sufficient statistics, so the quadratic surrogate the sweep uses is wrong
# and the delayed-acceptance step has to correct it. A sampler made of exact
# draws of mu | alpha, Sigma, of Sigma | alpha, mu, a and of a | Sigma and an
# exact random-walk step for alpha | mu, Sigma leaves the posterior
# invariant; adding the sweep (scale and location moves) must not change the
# stationary distribution.
toy_ll <- function(pars, dadm, ...) -sum(dadm$T * log1p((dadm$ybar - pars)^2))
toy_pars <- c("a", "b")
toy_T <- c(30, 3)
set.seed(21)
toy_dat <- do.call(rbind, lapply(1:2, function(s) data.frame(
  subjects = s, par = toy_pars, T = toy_T, ybar = c(.4, -.3) * (s - 1.5) * 2 + rnorm(2, 0, .2))))
toy_dat$subjects <- factor(toy_dat$subjects)
toy_design <- design(model = toy_ll, custom_p_vector = toy_pars, report_p_vector = FALSE)
toy_emc <- make_emc(toy_dat, toy_design, type = "standard", n_chains = 2, compress = FALSE)

toy_run <- function(sampler, iter, use_move, seed) {
  set.seed(seed)
  p <- 2; n <- 2; v <- sampler$prior$v; A <- sampler$prior$A
  m0 <- sampler$prior$theta_mu_mean; V0inv <- sampler$prior$theta_mu_invar
  mu <- c(.1, -.1); a <- c(.5, 2)
  Sigma <- diag(c(.3, .3)); alpha <- matrix(0, p, n, dimnames = list(toy_pars, NULL))
  ll <- function(s, x) as.numeric(EMC2:::calc_ll_manager(matrix(x, 1, dimnames = list(NULL, toy_pars)),
                                                         sampler$data[[s]], sampler$model))
  cur_ll <- sapply(1:n, function(s) ll(s, alpha[, s]))
  # a deliberately rough surrogate: curvature from the Gaussian part at the
  # likelihood's mode, which overstates the tails of this likelihood
  lik_prec <- lapply(1:n, function(s) {
    d <- sampler$data[[s]]
    list(prec = diag(2 * d$T), lin = 2 * d$T * d$ybar)
  })
  settings <- EMC2:::scale_move_init(NULL, toy_pars)
  out <- matrix(NA_real_, iter, 9, dimnames = list(NULL, c("s11", "s22", "r12", "a11", "a22", "aux1", "aux2", "mu1", "mu2")))
  for (i in seq_len(iter)) {
    # mu | alpha, Sigma; Sigma | alpha, mu, a ~ IW(v + p - 1 + n, 2 v diag(1/a) + S); a | Sigma
    Sinv <- solve(Sigma); Q <- V0inv + n * Sinv
    mu <- drop(mvtnorm::rmvnorm(1, solve(Q, V0inv %*% m0 + Sinv %*% rowSums(alpha)), solve(Q)))
    r <- alpha - mu
    Sigma <- EMC2:::riwish(v + p - 1 + n, 2 * v * diag(1 / a) + r %*% t(r))
    Sinv <- solve(Sigma)
    a <- 1 / rgamma(p, (v + p) / 2, rate = v * diag(Sinv) + 1 / A^2)
    # alpha_s | Sigma: random-walk Metropolis on the exact posterior
    for (s in 1:n) {
      prop <- alpha[, s] + rnorm(p, 0, .3)
      lp <- ll(s, prop) + mvtnorm::dmvnorm(prop, mu, Sigma, log = TRUE)
      lc <- cur_ll[s] + mvtnorm::dmvnorm(alpha[, s], mu, Sigma, log = TRUE)
      if (log(runif(1)) < lp - lc) { alpha[, s] <- prop; cur_ll[s] <- ll(s, prop) }
    }
    if (use_move) {
      pars <- list(tmu = mu, tvar = Sigma, tvinv = Sinv, a_half = a, subj_mu = matrix(mu, p, n), alpha = alpha)
      sm <- EMC2:::scale_move_standard(sampler, pars, alpha, cur_ll, settings, lik_prec, frozen = i > iter / 4)
      Sigma <- sm$pars$tvar; a <- sm$pars$a_half; mu <- sm$pars$tmu; alpha <- sm$alpha; cur_ll <- sm$ll
      settings <- sm$settings
    }
    out[i, ] <- c(Sigma[1, 1], Sigma[2, 2], Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2]), alpha[1, 1], alpha[2, 2], a, mu)
  }
  list(draws = out[-(1:(iter / 4)), ], settings = settings)
}

test_that("the scale move leaves the posterior invariant on a 2-subject toy", {
  skip_on_cran()
  sampler <- toy_emc[[1]]
  iter <- 12000
  ref <- toy_run(sampler, iter, use_move = FALSE, seed = 1)
  mv <- toy_run(sampler, iter, use_move = TRUE, seed = 2)
  # the move did something: sweeps proposed and accepted, step sizes tuned
  expect_gt(mv$settings$n_out, iter / 4)
  expect_gt(mv$settings$acc_out / mv$settings$n_out, .2)
  expect_true(all(mv$settings$step > 0) && all(mv$settings$loc > 0))
  expect_true(all(mv$settings$acc_loc > iter / 20))
  ess <- function(x) coda::effectiveSize(coda::mcmc(x))
  for (k in colnames(ref$draws)) {
    lg <- k %in% c("s11", "s22", "aux1", "aux2")
    x <- if (lg) log(mv$draws[, k]) else mv$draws[, k]
    y <- if (lg) log(ref$draws[, k]) else ref$draws[, k]
    z <- (mean(x) - mean(y)) / sqrt(var(x) / ess(x) + var(y) / ess(y))
    expect_lt(abs(z), 4)
    expect_lt(abs(log(sd(x) / sd(y))), .25)
  }
})

# ---- (3) the conjugate model: exact posterior, funnel autocorrelation gone --
# y_sj ~ N(alpha_sj, 1 / T_j) on sufficient statistics; the exact posterior
# of the group SDs comes from a Gibbs sampler of the four full conditionals
# (conj_exact below, as in rating-work/sampler/hier/stageH1/exact_gibbs.R):
# the reference moments are from two runs of 400000 sweeps (ESS of log SD
# 3300-6200 for the weak parameters; rating-work/sampler/hier/stageH2/
# diag_ref.R), a 20000-sweep reference being too noisy in the tails. With
# T_j = 2 the third parameter is the funnel: the centred sampler's log group
# SD has an autocorrelation time of 40-200 iterations here (about 500 with
# 20 subjects), the sweep brings it under 15.
conj_ll <- function(pars, dadm, ...) -0.5 * sum(dadm$T * (dadm$ybar - pars)^2)
conj_pars <- c("strong", "mid", "weak")
conj_T <- c(200, 20, 2)
conj_n <- 8
set.seed(11)
conj_alpha <- matrix(rnorm(conj_n * 3, rep(c(.2, -.2, 0), each = conj_n), .3), conj_n, 3)
conj_ybar <- conj_alpha + matrix(rnorm(conj_n * 3), conj_n, 3) / sqrt(rep(conj_T, each = conj_n))
conj_dat <- do.call(rbind, lapply(1:conj_n, function(s) data.frame(
  subjects = s, par = conj_pars, T = conj_T, ybar = conj_ybar[s, ])))
conj_dat$subjects <- factor(conj_dat$subjects)
conj_design <- design(model = conj_ll, custom_p_vector = conj_pars, report_p_vector = FALSE)
conj_emc <- make_emc(conj_dat, conj_design, type = "standard", n_chains = 2, compress = FALSE)
conj_stop <- list(preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
                  sample = list(iter = 1500))

conj_exact <- function(iter, seed = 3, thin = 1) {
  set.seed(seed)
  prior <- conj_emc[[1]]$prior; v <- prior$v; A <- prior$A
  p <- 3; n <- conj_n; Tj <- conj_T; ybar <- conj_ybar
  V0inv <- solve(prior$theta_mu_var); m0 <- prior$theta_mu_mean
  mu <- rep(0, p); Sigma <- diag(.1, p); a <- rep(1, p); alpha <- ybar
  out <- matrix(NA_real_, iter %/% thin, p)
  for (i in seq_len(iter)) {
    Sinv <- solve(Sigma)
    Q <- V0inv + n * Sinv
    mu <- drop(mvtnorm::rmvnorm(1, solve(Q, V0inv %*% m0 + Sinv %*% colSums(alpha)), solve(Q)))
    r <- sweep(alpha, 2, mu)
    Sigma <- EMC2:::riwish(v + p - 1 + n, 2 * v * diag(1 / a) + crossprod(r))
    Sinv <- solve(Sigma)
    a <- 1 / rgamma(p, (v + p) / 2, rate = v * diag(Sinv) + 1 / A^2)
    Qa <- Sinv + diag(Tj); Va <- solve(Qa)
    for (s in 1:n) alpha[s, ] <- mvtnorm::rmvnorm(1, Va %*% (Sinv %*% mu + Tj * ybar[s, ]), Va)
    if (i %% thin == 0) out[i %/% thin, ] <- sqrt(diag(Sigma))
  }
  out[-(1:(nrow(out) / 10)), ]
}
# log group SD under the exact posterior: conj_exact(4e5, seed = 3, thin = 5)
# and seed = 7, averaged
conj_ref <- list(mean = c(-1.077, -1.569, -1.820), sd = c(.276, .721, 1.148),
                 q50 = c(-1.10, -1.44, -1.60), q10 = c(-1.41, -2.355, -3.30))

test_that("with the scale move the conjugate posterior is recovered and the funnel mixes", {
  skip_on_cran()
  skip_on_os("windows")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 500)
  sm <- lapply(emc, function(x) attr(x$samples, "scale_move"))
  expect_named(sm[[1]]$step, conj_pars)
  expect_gt(sm[[1]]$n_out, 100)
  # the surrogate is exact for a Gaussian likelihood, so every sweep is accepted
  expect_gt(sm[[1]]$acc_out / sm[[1]]$n_out, .98)
  # more sample iterations: the step sizes stay frozen
  emc2 <- fit(emc, cores_for_chains = 1, iter = 1600, verbose = FALSE, particle_factor = 20,
              step_size = 500, stop_criteria = list(sample = list(iter = 1600)))
  for (ch in 1:2) {
    sm2 <- attr(emc2[[ch]]$samples, "scale_move")
    expect_identical(sm2$step, sm[[ch]]$step)
    expect_gt(sm2$n_in[1], sm[[ch]]$n_in[1])
  }
  # the exact posterior of the group SDs (mean within 4 Monte-Carlo SEs,
  # spread and median within what 300-400 effective draws of a heavy-tailed
  # log SD resolve)
  sd_draws <- lapply(emc, function(x) {
    i <- x$samples$stage == "sample"
    t(sqrt(apply(x$samples$theta_var[, , i, drop = FALSE], 3, diag)))
  })
  iat <- numeric(3)
  for (j in 1:3) {
    x <- log(do.call(rbind, sd_draws)[, j])
    ess_x <- sum(sapply(sd_draws, function(d) coda::effectiveSize(coda::mcmc(log(d[, j])))))
    z <- (mean(x) - conj_ref$mean[j]) / sqrt(var(x) / ess_x + (conj_ref$sd[j] / 60)^2)
    expect_lt(abs(z), 4)
    expect_lt(abs(log(sd(x) / conj_ref$sd[j])), .3)
    expect_lt(abs(median(x) - conj_ref$q50[j]), .25)
    expect_lt(abs(quantile(x, .1) - conj_ref$q10[j]), .4)
    iat[j] <- length(x) / ess_x
  }
  # the funnel parameter's log SD: a few iterations per effective draw
  expect_lt(iat[3], 15)
})

test_that("options(emc.scale_move = FALSE) and the legacy sampler run no move", {
  skip_on_os("windows")
  stop_short <- list(preburn = list(iter = 5), burn = list(mean_gd = 2.5), adapt = list(min_unique = 5),
                     sample = list(iter = 10))
  op <- options(emc.scale_move = FALSE); on.exit(options(op))
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = stop_short, verbose = FALSE,
             particle_factor = 20, step_size = 5)
  expect_null(attr(emc[[1]]$samples, "scale_move"))
  options(emc.scale_move = NULL, emc.sampler = "legacy")
  emc <- suppressWarnings(fit(conj_emc, cores_for_chains = 1, stop_criteria = stop_short, verbose = FALSE,
                              particle_factor = 20, step_size = 5))
  expect_null(attr(emc[[1]]$samples, "scale_move"))
})
