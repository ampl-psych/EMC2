# Chain culling and the flat test: the statistics on a healthy fit, a chain
# and a subject held below the others (simulated in the stored draws), the
# respawn mechanics, persistence, and a cull inside run_emc().
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:3], ]
dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR,
              formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")

test_that("cull settings are validated", {
  d <- EMC2:::cull_defaults(TRUE)
  expect_equal(d$persistence, 2L); expect_equal(d$chain_z, 2.5)
  expect_equal(EMC2:::cull_defaults(list(chain_z = 4))$chain_z, 4)
  expect_null(EMC2:::cull_defaults(FALSE))
  expect_error(EMC2:::cull_defaults(list(bogus = 1)), "unknown")
  expect_error(sc("adapt", list(cull = TRUE)), "burn and sample")
  expect_equal(sc("burn", list(cull = TRUE))$cull$min_window, 250L)
  expect_equal(sc("sample", list(cull = TRUE))$cull$min_window, 500L)
  expect_equal(sc("burn", list(cull = list(min_window = 400)))$cull$min_window, 400)
  expect_error(sc("burn", list(max_flat_loc = .3)), "sample stage")
  s <- sc("sample", list(max_flat_loc = .3))
  expect_equal(s$flat_selection, c("mu", "alpha", "LL", "subj_ll")); expect_equal(s$flat_mcse, 2)
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
})

test_that("stuck chains and subjects are found, respawned and culled", {
  skip_on_cran()
  old <- options(emc.adapt_converge = FALSE); on.exit(options(old))
  args <- list(cores_for_chains = 1, step_size = 50, verbose = FALSE, particle_factor = 30)
  set.seed(11); e <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 3)
  e <- do.call(run_emc, c(list(e, "preburn", sc("preburn", list(iter = 30))), args))
  e <- do.call(run_emc, c(list(e, "burn", sc("burn", list(iter = 100))), args))
  e <- do.call(run_emc, c(list(e, "adapt", sc("adapt", list(min_unique = 50))), args))
  e <- do.call(run_emc, c(list(e, "sample", sc("sample", list(iter = 200))), args))
  d <- chain_diagnostics(e, window = 150)
  expect_s3_class(d, "emc.diagnostics")
  expect_false(any(d$chains$flagged)); expect_equal(nrow(d$cells), 0)
  expect_true(all(c("mu", "alpha", "LL", "subj_ll") %in% names(d$flat)))
  expect_equal(d$flat$subj_ll$quantile, .99); expect_equal(d$flat$LL$quantile, 1)
  n <- e[[1]]$samples$idx; w <- (n - 149):n
  # a chain held below the others
  whole <- e; whole[[1]]$samples$subj_ll[, w] <- whole[[1]]$samples$subj_ll[, w] - 8
  dw <- chain_diagnostics(whole, window = 150)
  expect_true(dw$chains$flagged[1]); expect_false(any(dw$chains$flagged[-1]))
  # a subject held below the others, repaired from the best chain for it
  cell <- e; cell[[1]]$samples$subj_ll[2, w] <- cell[[1]]$samples$subj_ll[2, w] - 12
  dc <- chain_diagnostics(cell, window = 150)
  expect_false(any(dc$chains$flagged))
  expect_equal(nrow(dc$cells), 1); expect_equal(dc$cells$chain, 1)
  expect_equal(dc$cells$subject, rownames(e[[1]]$samples$subj_ll)[2])
  # respawn mechanics
  its <- EMC2:::stage_window(cell, "sample", 150)
  r <- EMC2:::respawn_subject(cell, 1, dc$cells$donor, 2, its[[dc$cells$donor]], 50)
  expect_false(identical(r[[1]]$samples$alpha[, 2, n], cell[[1]]$samples$alpha[, 2, n]))
  r2 <- EMC2:::respawn_chain(whole, 1, 2, its[[2]], 50)
  expect_true(any(apply(whole[[2]]$samples$theta_mu[, its[[2]]], 2, function(v) isTRUE(all.equal(v, r2[[1]]$samples$theta_mu[, n])))))
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
  expect_equal(res[[1]]$culls[[1]]$chain, 1)
  expect_lte(chain_n(res)[1, "sample"], 200)   # restarted: the old sample draws are gone
})
