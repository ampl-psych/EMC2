# Chain culling and the flat test: the statistics on a healthy fit, a chain
# held below the others and a subject apart in one chain (simulated in the
# stored draws), the respawn mechanics, persistence, and a cull inside run_emc().
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:3], ]
dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR,
              formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")

test_that("cull settings are validated", {
  d <- EMC2:::cull_defaults(TRUE)
  expect_equal(d$persistence, 2L); expect_equal(d$chain_z, 5); expect_equal(d$chain_gap, 200); expect_equal(d$max_culls, 1L)
  expect_equal(d$lp_guard_below, 1000)
  expect_equal(EMC2:::cull_defaults(list(chain_z = 4))$chain_z, 4)
  expect_null(EMC2:::cull_defaults(FALSE))
  expect_error(EMC2:::cull_defaults(list(bogus = 1)), "unknown")
  expect_error(EMC2:::cull_defaults(list(cell_gap = 5)), "unknown")      # the cell cull is gone
  expect_error(sc("adapt", list(cull = TRUE)), "burn and sample")
  expect_equal(sc("burn", list(cull = TRUE))$cull$min_window, 250L)
  expect_equal(sc("burn", list(cull = TRUE))$cull$max_culls, 2L); expect_equal(sc("sample", list(cull = TRUE))$cull$max_culls, 1L)
  expect_equal(sc("burn", list(cull = list(max_culls = 1)))$cull$max_culls, 1)
  expect_equal(sc("sample", list(cull = TRUE))$cull$min_window, 500L)
  expect_equal(sc("burn", list(cull = list(min_window = 400)))$cull$min_window, 400)
  expect_error(sc("burn", list(max_flat_loc = .3)), "sample stage")
  s <- sc("sample", list(max_flat_loc = .3))
  expect_equal(s$flat_selection, c("mu", "alpha", "LL", "subj_ll")); expect_equal(s$flat_mcse, 0)
  expect_error(sc("sample", list(max_flat_loc = .3, flat_selection = "theta")), "flat_selection")
})

test_that("flat statistics see drift and the flat test passes on noise", {
  set.seed(1)
  x <- rnorm(600)
  expect_lt(EMC2:::flat_series(x, 1/3, 1/3)[["drift"]], .3)
  f <- EMC2:::flat_series(x + seq(0, 3, length.out = 600), 1/3, 1/3)
  expect_gt(f[["drift"]], 1)
  expect_true(EMC2:::flat_done(list(mu = list(drift = .5, mcse = .4)), .3, 2))   # within 2 mcse of zero
  expect_false(EMC2:::flat_done(list(mu = list(drift = .5, mcse = .1)), .3, 2))
  expect_false(EMC2:::flat_done(list(mu = list(drift = .5, mcse = .4)), .3, 0))  # floor off
})

test_that("the Monte Carlo error and the window ESS are summed over the chains' own series", {
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

test_that("thin_to_ess keeps factor * target pooled draws and records what it did", {
  e <- thin_to_ess(samples_LNR, target = 15)       # 3 chains x 50 draws; keep 2 * 15 / 3 = 10 per chain
  expect_equal(unname(chain_n(e)[, "sample"]), rep(10, 3))
  th <- e[[1]]$thinning
  expect_equal(th$n_before, 50); expect_equal(th$n_after, 10); expect_equal(th$interval, 5)
  expect_true(all(is.finite(th$ess_before))); expect_true(all(is.finite(th$ess_after)))
  expect_lte(th$ess_after[["max"]], th$ess_before[["max"]] + 1e-8)
  same <- thin_to_ess(samples_LNR, target = 3000)  # already shorter than 6000 pooled draws: untouched
  expect_equal(unname(chain_n(same)[, "sample"]), rep(50, 3)); expect_equal(same[[1]]$thinning$interval, 1)
  expect_error(thin_to_ess(samples_LNR, target = -1), "target")
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
  d <- chain_diagnostics(e, window = 150)
  expect_s3_class(d, "emc.diagnostics")
  expect_false(any(d$chains$flagged)); expect_equal(nrow(d$subjects), 3)
  expect_true(all(c("mu", "alpha", "LL", "subj_ll") %in% names(d$flat)))
  expect_equal(d$flat$subj_ll$quantile, .99); expect_equal(d$flat$LL$quantile, 1)
  n <- e[[1]]$samples$idx; w <- (n - 149):n
  # a chain in a different region: far below on the log-likelihood (and so on its own log posterior)
  whole <- e; whole[[1]]$samples$subj_ll[, w] <- whole[[1]]$samples$subj_ll[, w] - 1000
  dw <- chain_diagnostics(whole, window = 150)
  expect_true(dw$chains$flagged[1]); expect_false(any(dw$chains$flagged[-1])); expect_true(dw$chains$lp_lower[1])
  # tens of units below is a wandering chain, not a different region: not flagged at the defaults
  near <- e; near[[1]]$samples$subj_ll[, w] <- near[[1]]$samples$subj_ll[, w] - 30
  expect_false(any(chain_diagnostics(near, window = 150)$chains$flagged))
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
  dc <- chain_diagnostics(cell, window = 150)
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
  # persistence: a flag must repeat before anything is culled
  cd <- EMC2:::cull_defaults(list(window = 150))
  dd <- EMC2:::cull_diagnostics(whole, "sample", cd)
  st <- EMC2:::update_cull_state(NULL, dd, 3)
  expect_length(EMC2:::apply_culls(whole, dd, "sample", cd, st)$done, 0)
  st <- EMC2:::update_cull_state(st, dd, 3)
  expect_length(EMC2:::apply_culls(whole, dd, "sample", cd, st)$done, 1)
  # inside run_emc: the chain is culled and the sample stage restarted
  res <- do.call(run_emc, c(list(whole, "sample", sc("sample", list(iter = 100, max_gd = 1.5, cull = list(window = 150, persistence = 1)))),
                            list(cores_for_chains = 1, step_size = 50, verbose = FALSE, max_tries = 3, particle_factor = 30)))
  expect_length(res[[1]]$culls, 1)
  expect_equal(res[[1]]$culls[[1]]$chain, 1); expect_false(res[[1]]$culls[[1]]$exit)
  expect_lte(chain_n(res)[1, "sample"], 200)   # restarted: the old sample draws are gone
  # burn exit: a chain in a different region at the end of burn is re-seeded, no persistence needed
  nb <- chain_n(eb)[1, "burn"]; wb <- eb[[1]]$samples$idx - seq_len(ceiling(2 * nb / 3)) + 1
  stuck <- eb; stuck[[2]]$samples$subj_ll[, wb] <- stuck[[2]]$samples$subj_ll[, wb] - 1000
  attr(stuck[[3]]$samples, "scale_move") <- list(step = c(m = 7))
  r <- EMC2:::burn_exit_reseed(stuck, EMC2:::cull_defaults(TRUE))
  expect_length(r[[1]]$culls, 1); expect_equal(r[[1]]$culls[[1]]$chain, 2); expect_true(r[[1]]$culls[[1]]$exit)
  expect_false(identical(r[[2]]$samples$theta_mu[, r[[2]]$samples$idx], stuck[[2]]$samples$theta_mu[, stuck[[2]]$samples$idx]))
  expect_length(EMC2:::burn_exit_reseed(eb, EMC2:::cull_defaults(TRUE))[[1]]$culls, 0)   # nothing out of place: untouched
})
