# The sample stage runs one fixed kernel, built and tuned at the end of adapt,
# and the hierarchical local proposal follows the current group level through
# the subject's likelihood precision (R/sampling.R, R/fitting.R).

# Conjugate normal model on sufficient statistics: the likelihood precision of
# parameter j is T_j and its linear term T_j * ybar_j, exactly.
conj_ll <- function(pars, dadm, ...) -0.5 * sum(dadm$T * (dadm$ybar - pars)^2)
conj_pars <- c("a", "b", "c", "d")
conj_T <- c(200, 20, 5, 2)
set.seed(11)
conj_dat <- do.call(rbind, lapply(1:4, function(s) data.frame(
  subjects = s, par = conj_pars, T = conj_T, ybar = rnorm(4, 0, .3) + rnorm(4) / sqrt(conj_T))))
conj_dat$subjects <- factor(conj_dat$subjects)
conj_design <- design(model = conj_ll, custom_p_vector = conj_pars, report_p_vector = FALSE)
conj_emc <- make_emc(conj_dat, conj_design, type = "standard", n_chains = 2, compress = FALSE)
conj_stop <- list(preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
                  sample = list(iter = 20))

test_that("lik_precision recovers the likelihood's quadratic form", {
  dadm <- conj_emc[[1]]$data[[1]]
  centre <- setNames(c(.1, -.2, .3, 0), conj_pars)
  # starting steps far too small and far too large for the likelihood's scale
  for (h in list(rep(1e-4, 4), rep(5, 4))) {
    lik <- EMC2:::lik_precision(centre, h, dadm, conj_emc[[1]]$model)
    expect_equal(unname(lik$prec), diag(conj_T), tolerance = 1e-6)
    expect_equal(unname(lik$lin), conj_T * dadm$ybar, tolerance = 1e-6)
  }
  cond <- EMC2:::conditional_proposal(lik, diag(1 / .09, 4), setNames(rep(0, 4), conj_pars))
  expect_equal(diag(cond$var), 1 / (conj_T + 1 / .09), tolerance = 1e-6)
  expect_equal(unname(cond$mu), conj_T * dadm$ybar / (conj_T + 1 / .09), tolerance = 1e-6)
})

test_that("the sample-stage kernel is fixed", {
  skip_on_os("windows")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  kernel <- emc[[1]]$sample_kernel
  expect_length(kernel$chains, 2)
  expect_named(kernel$chains[[1]], c("chains_var", "chains_mu", "eff_mu", "eff_var", "prop_var"))
  expect_equal(unname(kernel$lik_prec[[1]]$prec), diag(conj_T), tolerance = 1e-6)
  # the tail of adapt that tunes the kernel is stored as adapt
  expect_gte(chain_n(emc)[1, "adapt"], 100)
  # Stage H6b: adapt waited for the draws the kernel is built from (at least adapt_converge$min
  # iterations, then the tail), every adapt draw precedes the kept ones, and the kernel has the
  # draw-based likelihood precision
  expect_gte(chain_n(emc)[1, "adapt"], EMC2:::adapt_converge$min + 100)
  expect_lte(chain_n(emc)[1, "adapt"], EMC2:::adapt_converge$max + 120)
  st <- emc[[1]]$samples$stage[seq_len(emc[[1]]$samples$idx)]
  expect_lt(max(which(st == "adapt")), min(which(st == "sample")))
  expect_equal(sum(st == "sample"), 20)
  # ... where the group level is stable over the window it is built from (lik_prec_draws_cv)
  idx <- max(which(st == "adapt")) - 100           # the kernel is built before the 100 tuning iterations
  win <- unique(pmax(1, round(idx - min(250, idx / 1.5)):idx - 1))
  if (EMC2:::group_precision_cv(emc, win) < EMC2:::lik_prec_draws_cv) {
    expect_named(kernel$lik_prec[[1]]$post, c("prec", "lin", "n_draws", "n_capped", "tau", "ess"))
  } else expect_null(kernel$lik_prec[[1]]$post)
  pm <- lapply(emc, function(x) attr(x$samples, "pm_settings"))
  expect_length(pm[[1]][[1]][[1]]$mix, 4)
  # more sample iterations: same proposals, same step size, mixing weights and particles
  emc2 <- fit(emc, cores_for_chains = 1, iter = 40, verbose = FALSE, particle_factor = 20,
              step_size = 20, stop_criteria = list(sample = list(iter = 40)))
  expect_identical(emc2[[1]]$sample_kernel, kernel)
  keep <- c("epsilon", "mix", "n_particles")
  for (ch in 1:2) for (s in 1:4) {
    expect_identical(attr(emc2[[ch]]$samples, "pm_settings")[[s]][[1]][keep], pm[[ch]][[s]][[1]][keep])
  }
})

test_that("adapt's convergence rule can be switched off, and a fit saved before it existed carries on with its kernel", {
  skip_on_os("windows")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  op <- options(emc.adapt_converge = FALSE)
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  options(op)
  # min_unique alone: far fewer adapt iterations than adapt_converge$min, plus the tail
  expect_lt(chain_n(emc)[1, "adapt"], EMC2:::adapt_converge$min)
  expect_length(emc[[1]]$sample_kernel$chains, 2)
  # as a fit saved by a build without the convergence rule, the draw-based precision and the sweep's gate
  for (s in seq_along(emc[[1]]$sample_kernel$lik_prec)) emc[[1]]$sample_kernel$lik_prec[[s]]$post <- NULL
  for (ch in 1:2) {
    sm <- attr(emc[[ch]]$samples, "scale_move")
    attr(emc[[ch]]$samples, "scale_move") <- sm[setdiff(names(sm), c("active_scale", "active_loc", "gate", "mark"))]
  }
  kernel <- emc[[1]]$sample_kernel
  emc2 <- fit(emc, cores_for_chains = 1, iter = 40, verbose = FALSE, particle_factor = 20,
              step_size = 20, stop_criteria = list(sample = list(iter = 40)))
  expect_identical(emc2[[1]]$sample_kernel, kernel)
  expect_equal(unname(chain_n(emc2)[1, "sample"]), 60)
  sm <- attr(emc2[[1]]$samples, "scale_move")
  expect_true(all(sm$active_scale) && all(sm$active_loc))
})

test_that("the legacy sampler keeps adapting in the sample stage", {
  skip_on_os("windows")
  op <- options(emc.sampler = "legacy"); on.exit(options(op))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  expect_null(emc[[1]]$sample_kernel)
  expect_null(emc[[1]]$lik_prec)
})

# ---- Stage H6b: the draw-based likelihood precision ------------------------
test_that("adapt_stop: converged, stalled, at the limit, or carry on", {
  stop_rule <- EMC2:::adapt_stop
  expect_equal(EMC2:::adapt_converge$stall, 3)
  expect_equal(stop_rule(c(1.35, 1.12), 400), "converged")
  expect_equal(stop_rule(c(1.5, 1.4, 1.3, 1.25), 600), "")              # still improving
  # forstmann's adapt (H6b Part B, replicate 1): no new minimum after 400 iterations
  expect_equal(stop_rule(c(1.658, 1.657, 2.327), 500), "")
  expect_equal(stop_rule(c(1.658, 1.657, 2.327, 2.260), 600), "")
  expect_equal(stop_rule(c(1.658, 1.657, 2.327, 2.260, 2.273), 700), "stalled")
  # a new minimum restarts the count; a tie with the best is no improvement
  expect_equal(stop_rule(c(2, 1.9, 1.95, 1.96, 1.8, 1.85), 800), "")
  expect_equal(stop_rule(c(1.5, 1.5, 1.6, 1.5), 600), "stalled")
  expect_equal(stop_rule(c(1.5, 1.4, 1.3, 1.25), EMC2:::adapt_converge$max), "limit")
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
  # Eisenberg-like (SD of log group SD .07): CV of the precision about .14; funnel-like (SD 1): far above
  cv_stable <- EMC2:::group_precision_cv(mock(.07), 1:250)
  cv_funnel <- EMC2:::group_precision_cv(mock(1), 1:250)
  expect_lt(cv_stable, .25); expect_gt(cv_stable, .1)
  expect_gt(cv_funnel, 2)
  expect_equal(EMC2:::lik_prec_draws_cv, .5)
  # the funnel-like window gets no draw-based precision, whatever else it has
  expect_null(EMC2:::window_group_level(mock(1), 1:250, 750))
  # not computable (a zero variance): Inf, so no draw-based precision either
  bad <- mock(.07); bad[[2]]$samples$theta_var[1, 1, 5] <- 0
  expect_identical(EMC2:::group_precision_cv(bad, 1:250), Inf)
})
test_that("adapt_converged reads the window the kernel is built from", {
  skip_on_os("windows")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  class(emc) <- "emc"
  emc <- EMC2:::restore_duplicates(emc)
  # too few adapt iterations: not yet; at the limit: stop whatever the draws say
  expect_false(as.vector(EMC2:::adapt_converged(emc, EMC2:::adapt_converge$min - 1)))
  expect_true(as.vector(EMC2:::adapt_converged(emc, EMC2:::adapt_converge$max)))
  r <- EMC2:::window_rhat(emc, EMC2:::kernel_window)
  expect_equal(dim(r), c(4, 4))
  a <- EMC2:::adapt_converged(emc, 300)
  expect_identical(as.vector(a), max(r) < EMC2:::adapt_converge$rhat)
  expect_equal(attr(a, "rhat"), max(r))
  # the history of earlier checks is carried and extended
  expect_equal(attr(EMC2:::adapt_converged(emc, 300, history = c(2, 1.5)), "rhat"), c(2, 1.5, max(r)))
  # chains that sit apart in one subject's parameter: not converged
  bad <- emc; n <- bad[[1]]$samples$idx
  bad[[1]]$samples$alpha[2, 3, ] <- bad[[1]]$samples$alpha[2, 3, ] + 10
  expect_gt(EMC2:::window_rhat(bad, EMC2:::kernel_window)[2, 3], 3)
  expect_false(as.vector(EMC2:::adapt_converged(bad, 300)))
  # ... unless the window Rhat has stopped improving: no new minimum in adapt_converge$stall checks
  expect_true(as.vector(EMC2:::adapt_converged(bad, 600, history = c(1.5, 4, 4, 4)[seq_len(EMC2:::adapt_converge$stall)])))
  # the legacy sampler and single-subject fits keep the old rule
  op <- options(emc.sampler = "legacy"); expect_true(EMC2:::adapt_converged(bad, 10)); options(op)
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
  out <- EMC2:::lik_precision_draws(X, P, mu, fd, n_chains = 3)
  expect_equal(out$n_draws, 2); expect_equal(out$n_capped, 0); expect_equal(out$tau, .5)
  expect_equal(unname(diag(out$prec)), conj_T, tolerance = .15)
  expect_lt(max(abs(out$prec[upper.tri(out$prec)]) / sqrt(outer(conj_T, conj_T))[upper.tri(out$prec)]), .15)
  cond <- EMC2:::conditional_proposal(out, P, mu)
  expect_lt(max(abs(cond$mu - rowMeans(X)) / apply(X, 1, sd)), .05)
  expect_equal(unname(diag(cond$var)), unname(apply(X, 1, var)), tolerance = .1)
  # (b) a collapsed group level (group SD .01): the draws show the prior only, every direction is the
  # finite differences', exactly
  P <- diag(1e4, 4)
  X <- post_draws(conj_T, ybar, P, mu, 3000)
  out <- EMC2:::lik_precision_draws(X, P, mu, fd, n_chains = 3)
  expect_equal(out$n_draws, 0); expect_equal(out$n_capped, 0)
  expect_equal(out$prec, fd$prec, tolerance = 1e-8)
  expect_equal(out$lin, fd$lin, tolerance = 1e-8)
  # (c) finite differences that claim far more than the draws allow in a direction the prior dominates
  # (a curvature taken at one point of a skewed likelihood): capped at tau x the group precision, and the
  # conditional mean stays at the draws' mean
  P <- diag(1 / .09, 4)
  X <- post_draws(conj_T, ybar, P, mu, 3000)
  bad <- fd; bad$prec["d", "d"] <- 2000; bad$lin["d"] <- 2000 * 5
  out <- EMC2:::lik_precision_draws(X, P, mu, bad, n_chains = 3)
  expect_equal(out$n_capped, 1)
  expect_equal(unname(out$prec["d", "d"]), .5 / .09, tolerance = .05)
  cond <- EMC2:::conditional_proposal(out, P, mu)
  expect_lt(max(abs(cond$mu - rowMeans(X)) / apply(X, 1, sd)), .05)
  # the finite-difference quadratic's cross terms with the directions taken from the draws go into its
  # linear term: with a strong a-c cross term the conditional mean is still at the draws' mean
  Lc <- diag(conj_T); Lc[1, 3] <- Lc[3, 1] <- 25; dimnames(Lc) <- list(conj_pars, conj_pars)
  set.seed(6); S <- solve(Lc + P); m <- drop(S %*% (Lc %*% ybar + P %*% mu))
  X <- t(mvtnorm::rmvnorm(3000, m, S)); rownames(X) <- conj_pars
  out <- EMC2:::lik_precision_draws(X, P, mu, list(prec = Lc, lin = setNames(drop(Lc %*% ybar), conj_pars)), n_chains = 3)
  expect_lt(out$n_draws, 4)
  cond <- EMC2:::conditional_proposal(out, P, mu)
  expect_lt(max(abs(cond$mu - rowMeans(X)) / apply(X, 1, sd)), .1)
  # (d) too few distinct draws: no estimate
  expect_null(EMC2:::lik_precision_draws(X[, rep(1:10, 30)], P, mu, fd, n_chains = 3))
})
