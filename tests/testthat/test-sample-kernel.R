# The sample stage runs one fixed kernel, built and tuned at the end of adapt,
# and the hierarchical local proposal follows the current group level through
# the subject's likelihood precision (R/sampling.R, R/fitting.R).

conj_pars <- c("a", "b", "c", "d")
conj_T <- c(200, 20, 5, 2)
conj_ybar <- withr::with_seed(11, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                              .rng_sample_kind = "Rejection",
                              t(sapply(1:4, function(s) rnorm(4, 0, .3) + rnorm(4) / sqrt(conj_T))))
colnames(conj_ybar) <- conj_pars
conj_emc <- conj_emc_from(conj_ybar, conj_T)
conj_stop <- list(preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
                  sample = list(iter = 20))
# one fit, shared by the tests that need it
conj_fit <- local({
  emc <- NULL
  function() {
    if (is.null(emc)) emc <<- withr::with_seed(123, .rng_kind = "L'Ecuyer-CMRG",
      fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE, particle_factor = 20, step_size = 20))
    emc
  }
})

test_that("conditional_proposal combines the likelihood's quadratic form with the group level", {
  ybar <- unname(conj_ybar[1, ])
  lik <- list(prec = diag(conj_T), lin = conj_T * ybar)
  cond <- conditional_proposal(lik, diag(1 / .09, 4), setNames(rep(0, 4), conj_pars))
  expect_equal(diag(cond$var), 1 / (conj_T + 1 / .09), tolerance = 1e-6)
  expect_equal(unname(cond$mu), conj_T * ybar / (conj_T + 1 / .09), tolerance = 1e-6)
})

test_that("the sample-stage kernel is fixed", {
  skip_on_cran()
  skip_on_os("windows")
  emc <- conj_fit()
  kernel <- emc[[1]]$sample_kernel
  expect_length(kernel$chains, 2)
  expect_named(kernel$chains[[1]], c("chains_var", "chains_mu", "eff_mu", "eff_var", "prop_var"))
  expect_equal(unname(kernel$lik_prec[[1]]$prec), diag(conj_T), tolerance = 1e-6)
  # adapt waited for the draws the kernel is built from, then ran the tail that tunes it (stored as
  # adapt); every adapt draw precedes the kept ones
  expect_gte(chain_n(emc)[1, "adapt"], adapt_converge$min + 100)
  expect_lte(chain_n(emc)[1, "adapt"], adapt_converge$max + 120)
  st <- emc[[1]]$samples$stage[seq_len(emc[[1]]$samples$idx)]
  expect_lt(max(which(st == "adapt")), min(which(st == "sample")))
  expect_equal(sum(st == "sample"), 20)
  pm <- lapply(emc, function(x) attr(x$samples, "pm_settings"))
  expect_length(pm[[1]][[1]][[1]]$mix, 4)
  sm <- lapply(emc, function(x) attr(x$samples, "scale_move"))
  # more sample iterations: same proposals, step size, mixing weights and particles, same sweep steps
  emc2 <- withr::with_seed(124, .rng_kind = "L'Ecuyer-CMRG",
    fit(emc, cores_for_chains = 1, iter = 40, verbose = FALSE, particle_factor = 20,
        step_size = 20, stop_criteria = list(sample = list(iter = 40))))
  expect_identical(emc2[[1]]$sample_kernel, kernel)
  keep <- c("epsilon", "mix", "n_particles")
  for (ch in 1:2) {
    for (s in 1:4) expect_identical(attr(emc2[[ch]]$samples, "pm_settings")[[s]][[1]][keep], pm[[ch]][[s]][[1]][keep])
    sm2 <- attr(emc2[[ch]]$samples, "scale_move")
    expect_identical(sm2$step, sm[[ch]]$step)
    expect_identical(sm2$r_block, sm[[ch]]$r_block)
    expect_gt(sm2$n_in[1], sm[[ch]]$n_in[1])
  }
})

test_that("adapt's convergence rule can be switched off", {
  skip_on_cran()
  skip_on_os("windows")
  withr::local_seed(123, .rng_kind = "L'Ecuyer-CMRG")
  withr::local_options(emc.adapt_converge = FALSE)
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  # min_unique alone: far fewer adapt iterations than adapt_converge$min, plus the tail
  expect_lt(chain_n(emc)[1, "adapt"], adapt_converge$min)
  expect_length(emc[[1]]$sample_kernel$chains, 2)
})

test_that("the legacy sampler builds no fixed kernel and runs no sweep", {
  skip_on_cran()
  skip_on_os("windows")
  withr::local_seed(123, .rng_kind = "L'Ecuyer-CMRG")
  withr::local_options(emc.sampler = "legacy")
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  expect_null(emc[[1]]$sample_kernel)
  expect_null(emc[[1]]$lik_prec)
  expect_null(attr(emc[[1]]$samples, "scale_move"))
})

test_that("adapt_stop: converged, stalled, at the limit, or carry on", {
  expect_equal(adapt_stop(c(1.35, 1.12), 400), "converged")
  expect_equal(adapt_stop(c(1.5, 1.4, 1.3, 1.25), 600), "")              # still improving
  # no new minimum after the second check
  expect_equal(adapt_stop(c(1.658, 1.657, 2.327), 500), "")
  expect_equal(adapt_stop(c(1.658, 1.657, 2.327, 2.260), 600), "")
  expect_equal(adapt_stop(c(1.658, 1.657, 2.327, 2.260, 2.273), 700), "stalled")
  # a new minimum restarts the count; a tie with the best is no improvement
  expect_equal(adapt_stop(c(2, 1.9, 1.95, 1.96, 1.8, 1.85), 800), "")
  expect_equal(adapt_stop(c(1.5, 1.5, 1.6, 1.5), 600), "stalled")
  expect_equal(adapt_stop(c(1.5, 1.4, 1.3, 1.25), adapt_converge$max), "limit")
})

test_that("the draw-based likelihood precision is used only where the group level is stable over the window", {
  # a mock emc: 3 chains x 250 iterations of a 3-parameter group covariance
  mock <- function(sd_log_sd, seed = 1) {
    set.seed(seed)
    lapply(1:3, function(i) {
      tv <- array(0, c(3, 3, 250))
      for (it in 1:250) tv[, , it] <- diag(exp(2 * (log(c(.3, .5, .2)) + sd_log_sd * rnorm(3))))
      list(type = "standard", nuisance = rep(FALSE, 3), samples = list(theta_var = tv))
    })
  }
  # a stable group level (SD of the log group SD .07): CV of the precision about .14; a funnel (SD 1): far above
  cv_stable <- group_precision_cv(mock(.07), 1:250)
  cv_funnel <- group_precision_cv(mock(1), 1:250)
  expect_lt(cv_stable, .25); expect_gt(cv_stable, .1)
  expect_gt(cv_funnel, 2)
  # the funnel-like window gets no draw-based precision, whatever else it has
  expect_null(window_group_level(mock(1), 1:250, 750))
  # not computable (a zero variance): Inf, so no draw-based precision either
  bad <- mock(.07); bad[[2]]$samples$theta_var[1, 1, 5] <- 0
  expect_identical(group_precision_cv(bad, 1:250), Inf)
})

test_that("adapt_converged reads the window the kernel is built from", {
  skip_on_cran()
  skip_on_os("windows")
  emc <- conj_fit()
  class(emc) <- "emc"
  emc <- restore_duplicates(emc)
  # too few adapt iterations: not yet; at the limit: stop whatever the draws say
  expect_false(as.vector(adapt_converged(emc, adapt_converge$min - 1)))
  expect_true(as.vector(adapt_converged(emc, adapt_converge$max)))
  r <- window_rhat(emc, kernel_window)
  expect_equal(dim(r), c(4, 4))
  a <- adapt_converged(emc, 300)
  expect_identical(as.vector(a), max(r) < adapt_converge$rhat)
  expect_equal(attr(a, "rhat"), max(r))
  # the history of earlier checks is carried and extended
  expect_equal(attr(adapt_converged(emc, 300, history = c(2, 1.5)), "rhat"), c(2, 1.5, max(r)))
  # chains that sit apart in one subject's parameter: not converged
  bad <- emc
  bad[[1]]$samples$alpha[2, 3, ] <- bad[[1]]$samples$alpha[2, 3, ] + 10
  expect_gt(window_rhat(bad, kernel_window)[2, 3], 3)
  expect_false(as.vector(adapt_converged(bad, 300)))
  # ... unless the window Rhat has stopped improving: no new minimum in adapt_converge$stall checks
  expect_true(as.vector(adapt_converged(bad, 600, history = c(1.5, 4, 4, 4)[seq_len(adapt_converge$stall)])))
  # the legacy sampler and single-subject fits keep the old rule
  withr::with_options(list(emc.sampler = "legacy"), expect_true(adapt_converged(bad, 10)))
})

# Draws from the exact conditional posterior of one subject of the conjugate
# model: likelihood precision diag(T), prior N(mu, P^-1).
post_draws <- function(T, ybar, P, mu, N, seed = 5) {
  set.seed(seed)
  S <- solve(diag(T) + P); m <- drop(S %*% (T * ybar + P %*% mu))
  X <- t(mvtnorm::rmvnorm(N, m, S)); rownames(X) <- conj_pars
  X
}

test_that("lik_precision_draws: the draws' estimate where the likelihood shows, the finite differences elsewhere", {
  ybar <- c(.1, -.2, .3, 0); mu <- setNames(c(0, .1, -.1, .05), conj_pars)
  fd <- list(prec = diag(conj_T), lin = conj_T * ybar)
  dimnames(fd$prec) <- list(conj_pars, conj_pars); names(fd$lin) <- conj_pars
  # (a) group SD .3: the likelihood dominates in a and b (T / P = 18, 1.8), the prior in c and d (.45, .18)
  P <- diag(1 / .09, 4)
  X <- post_draws(conj_T, ybar, P, mu, 3000)
  out <- lik_precision_draws(X, P, mu, fd, n_chains = 3)
  expect_equal(out$n_draws, 2); expect_equal(out$n_capped, 0)
  expect_equal(unname(diag(out$prec)), conj_T, tolerance = .15)
  expect_lt(max(abs(out$prec[upper.tri(out$prec)]) / sqrt(outer(conj_T, conj_T))[upper.tri(out$prec)]), .15)
  cond <- conditional_proposal(out, P, mu)
  expect_lt(max(abs(cond$mu - rowMeans(X)) / apply(X, 1, sd)), .05)
  expect_equal(unname(diag(cond$var)), unname(apply(X, 1, var)), tolerance = .1)
  # (b) a collapsed group level (group SD .01): the draws show the prior only, every direction is the
  # finite differences', exactly
  P <- diag(1e4, 4)
  X <- post_draws(conj_T, ybar, P, mu, 3000)
  out <- lik_precision_draws(X, P, mu, fd, n_chains = 3)
  expect_equal(out$n_draws, 0); expect_equal(out$n_capped, 0)
  expect_equal(out$prec, fd$prec, tolerance = 1e-8)
  expect_equal(out$lin, fd$lin, tolerance = 1e-8)
  # (c) finite differences that claim far more than the draws allow in a direction the prior dominates
  # (a curvature taken at one point of a skewed likelihood): capped at tau x the group precision, and the
  # conditional mean stays at the draws' mean
  P <- diag(1 / .09, 4)
  X <- post_draws(conj_T, ybar, P, mu, 3000)
  bad <- fd; bad$prec["d", "d"] <- 2000; bad$lin["d"] <- 2000 * 5
  out <- lik_precision_draws(X, P, mu, bad, n_chains = 3)
  expect_equal(out$n_capped, 1)
  expect_equal(unname(out$prec["d", "d"]), .5 / .09, tolerance = .05)
  cond <- conditional_proposal(out, P, mu)
  expect_lt(max(abs(cond$mu - rowMeans(X)) / apply(X, 1, sd)), .05)
  # the finite-difference quadratic's cross terms with the directions taken from the draws go into its
  # linear term: with a strong a-c cross term the conditional mean is still at the draws' mean
  Lc <- diag(conj_T); Lc[1, 3] <- Lc[3, 1] <- 25; dimnames(Lc) <- list(conj_pars, conj_pars)
  set.seed(6); S <- solve(Lc + P); m <- drop(S %*% (Lc %*% ybar + P %*% mu))
  X <- t(mvtnorm::rmvnorm(3000, m, S)); rownames(X) <- conj_pars
  out <- lik_precision_draws(X, P, mu, list(prec = Lc, lin = setNames(drop(Lc %*% ybar), conj_pars)), n_chains = 3)
  expect_lt(out$n_draws, 4)
  cond <- conditional_proposal(out, P, mu)
  expect_lt(max(abs(cond$mu - rowMeans(X)) / apply(X, 1, sd)), .1)
  # (d) too few distinct draws: no estimate
  expect_null(lik_precision_draws(X[, rep(1:10, 30)], P, mu, fd, n_chains = 3))
})
