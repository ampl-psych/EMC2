# Chain culling and the drift test: the statistics on a healthy fit, a chain
# held below the others and a subject apart in one chain (simulated in the
# stored draws), the respawn mechanics, persistence, the burn-exit re-seed.
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:3], ]
dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR,
              formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")

test_that("cull settings are validated and adjusted per stage", {
  expect_equal(EMC2:::cull_defaults(list(chain_z = 4))$chain_z, 4)
  expect_null(EMC2:::cull_defaults(FALSE)); expect_null(EMC2:::cull_settings(FALSE, "burn"))
  expect_error(EMC2:::cull_defaults(list(bogus = 1)), "unknown")
  expect_error(EMC2:::cull_defaults("yes"), "TRUE, FALSE or a list")
  b <- EMC2:::cull_settings(TRUE, "burn"); s <- EMC2:::cull_settings(TRUE, "sample")
  expect_lt(b$min_window, s$min_window); expect_gt(b$max_culls, s$max_culls)   # burn: shorter window, more re-seeds
  expect_equal(EMC2:::cull_settings(list(max_culls = 1), "burn")$max_culls, 1)
  expect_equal(EMC2:::cull_settings(list(min_window = 400), "burn")$min_window, 400)
})

test_that("the drift statistic sees a trend, and its Monte Carlo error is summed over the chains' own series", {
  set.seed(1)
  x <- rnorm(600)
  expect_lt(EMC2:::flat_series(x, 1/3, 1/3)[["drift"]], .3)
  expect_gt(EMC2:::flat_series(x + seq(0, 3, length.out = 600), 1/3, 1/3)[["drift"]], 1)
  set.seed(2)
  ar1 <- function(n, rho){ x <- numeric(n); x[1] <- rnorm(1); for(t in 2:n) x[t] <- rho * x[t - 1] + rnorm(1, 0, sqrt(1 - rho^2)); x }
  pool <- function(chains) as.vector(t(sapply(chains, identity)))          # iteration-major: chain varies fastest
  white <- pool(replicate(4, rnorm(1500), simplify = FALSE))
  f <- EMC2:::flat_series(white, 1/3, 1/3, nch = 4)
  expect_gt(f[["ess"]], 5000)
  # mcse of a median difference between thirds of 2000 independent draws, in IQR units: about 1.25*sqrt(2/2000)/1.35
  expect_lt(abs(f[["mcse"]] - 1.25 * sqrt(2 / 2000) / 1.349), .02)
  slow <- pool(replicate(4, ar1(1500, .99), simplify = FALSE))
  g <- EMC2:::flat_series(slow, 1/3, 1/3, nch = 4)
  expect_lt(g[["ess"]], 150)                                               # about 4 x 1500 / 199
  expect_gt(g[["mcse"]], 5 * f[["mcse"]])                                   # the autocorrelation is in the error
})

test_that("stuck chains and subjects are found, respawned and culled", {
  skip_on_cran()
  old <- options(emc.adapt_converge = FALSE); on.exit(options(old))
  args <- list(cores_for_chains = 1, step_size = 50, verbose = FALSE, particle_factor = 30)
  set.seed(11); e <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 3)
  e <- do.call(run_emc, c(list(e, "preburn", sc("preburn", list(iter = 30))), args))
  e <- do.call(run_emc, c(list(e, "burn", sc("burn", list(iter = 100))), args))
  eb <- e                                            # ends in burn: for the burn-exit re-seed below
  e <- do.call(run_emc, c(list(e, "adapt", sc("adapt", list(min_unique = 50))), args))
  e <- do.call(run_emc, c(list(e, "sample", sc("sample", list(iter = 200))), args))
  d <- diagnose(e, window = 150)
  expect_s3_class(d, "emc.diagnose")
  expect_false(any(d$chains$flagged)); expect_equal(nrow(d$subjects), 3)
  expect_true(all(c("mu", "alpha", "LL", "subj_ll") %in% names(d$flat)))
  expect_equal(d$flat$subj_ll$quantile, .99); expect_equal(d$flat$LL$quantile, 1)
  expect_output(print(d), "Class:")
  n <- e[[1]]$samples$idx; w <- (n - 149):n
  # a chain in a different region: far below on the log-likelihood (and so on its own log posterior)
  whole <- e; whole[[1]]$samples$subj_ll[, w] <- whole[[1]]$samples$subj_ll[, w] - 1000
  dw <- diagnose(whole, window = 150)
  expect_true(dw$chains$flagged[1]); expect_false(any(dw$chains$flagged[-1])); expect_true(dw$chains$lp_lower[1])
  expect_equal(dw$class, "region_mode"); expect_equal(dw$stats$ll_worst_chain, 1)
  # tens of units below is a wandering chain, not a different region: not flagged at the defaults
  near <- e; near[[1]]$samples$subj_ll[, w] <- near[[1]]$samples$subj_ll[, w] - 30
  expect_false(any(diagnose(near, window = 150)$chains$flagged))
  # lower on the log-likelihood but higher on its own log posterior: not flagged
  lik <- e; lik[[1]]$samples$subj_ll[, w] <- lik[[1]]$samples$subj_ll[, w] - 40
  p <- nrow(lik[[1]]$samples$theta_mu)
  for(k in w){ lik[[1]]$samples$theta_var[, , k] <- diag(1e-12, p); lik[[1]]$samples$alpha[, , k] <- lik[[1]]$samples$theta_mu[, k] }
  dl <- EMC2:::cull_diagnostics(lik, "sample", EMC2:::cull_defaults(list(window = 150, chain_gap = 20, chain_z = 1)))
  expect_gt(dl$chains$gap[1], 20); expect_false(dl$chains$lp_lower[1]); expect_false(dl$chains$flagged[1])
  # above lp_guard_below the likelihood decides: the same chain is flagged once the gap clears it
  dg <- EMC2:::cull_diagnostics(lik, "sample", EMC2:::cull_defaults(list(window = 150, chain_gap = 20, chain_z = 1, lp_guard_below = 30)))
  expect_false(dg$chains$lp_lower[1]); expect_true(dg$chains$flagged[1])
  # two subjects apart in one chain, the chain's sum unchanged: reported in the subject spread, nothing culled
  cell <- e; cell[[1]]$samples$subj_ll[2, w] <- cell[[1]]$samples$subj_ll[2, w] - 12
  cell[[1]]$samples$subj_ll[3, w] <- cell[[1]]$samples$subj_ll[3, w] + 12
  dc <- diagnose(cell, window = 150)
  expect_false(any(dc$chains$flagged))
  sub <- rownames(e[[1]]$samples$subj_ll)
  expect_setequal(dc$subjects$subject[1:2], sub[2:3]); expect_true(all(dc$subjects$spread[1:2] > 10))
  expect_equal(dc$subjects$worst[dc$subjects$subject == sub[2]], 1); expect_equal(dc$subjects$best[dc$subjects$subject == sub[3]], 1)
  # respawn mechanics: the donor's state and its kernel settings, the recipient's own random number state
  its <- EMC2:::stage_window(whole, "sample", 150)
  marked <- whole; attr(marked[[2]]$samples, "scale_move") <- list(step = c(v = 9)); attr(marked[[1]]$samples, "rng") <- "own"
  r2 <- EMC2:::respawn_chain(marked, 1, 2, its[[2]], 50)
  expect_true(any(apply(whole[[2]]$samples$theta_mu[, its[[2]]], 2, function(v) isTRUE(all.equal(v, r2[[1]]$samples$theta_mu[, n])))))
  expect_identical(attr(r2[[1]]$samples, "scale_move"), list(step = c(v = 9)))
  expect_identical(attr(r2[[1]]$samples, "pm_settings"), attr(marked[[2]]$samples, "pm_settings"))
  expect_identical(attr(r2[[1]]$samples, "rng"), "own")
  # persistence: a flag must repeat before anything is culled; then max_culls holds
  cd <- EMC2:::cull_defaults(list(window = 150))
  s1 <- EMC2:::cull_step(whole, "sample", cd, NULL)
  expect_length(s1$done, 0); expect_equal(s1$state$chain_count, c(1L, 0L, 0L))
  s2 <- EMC2:::cull_step(whole, "sample", cd, s1$state)
  expect_length(s2$done, 1); expect_equal(s2$done[[1]]$chain, 1); expect_match(s2$note, "chain 1 respawned")
  expect_length(s2$emc[[1]]$culls, 1); expect_equal(s2$state$n_culls, 1L); expect_equal(s2$state$chain_count, c(0L, 0L, 0L))
  s3 <- EMC2:::cull_step(whole, "sample", cd, s2$state)
  expect_length(s3$done, 0); expect_identical(s3$emc, whole)
  # burn exit: a chain in a different region at the end of burn is re-seeded, no persistence needed
  nb <- chain_n(eb)[1, "burn"]; wb <- eb[[1]]$samples$idx - seq_len(ceiling(2 * nb / 3)) + 1
  stuck <- eb; stuck[[2]]$samples$subj_ll[, wb] <- stuck[[2]]$samples$subj_ll[, wb] - 1000
  attr(stuck[[3]]$samples, "scale_move") <- list(step = c(m = 7))
  r <- EMC2:::burn_exit_reseed(stuck, EMC2:::cull_defaults(TRUE))
  expect_length(r[[1]]$culls, 1); expect_equal(r[[1]]$culls[[1]]$chain, 2); expect_true(r[[1]]$culls[[1]]$exit)
  expect_false(identical(r[[2]]$samples$theta_mu[, r[[2]]$samples$idx], stuck[[2]]$samples$theta_mu[, stuck[[2]]$samples$idx]))
  expect_length(EMC2:::burn_exit_reseed(eb, EMC2:::cull_defaults(TRUE))[[1]]$culls, 0)   # nothing out of place: untouched
})
