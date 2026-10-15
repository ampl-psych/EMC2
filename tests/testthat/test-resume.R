# Running a stage again on a saved or finished object must work at every
# stage boundary, including more adapt after adapt's tail has tuned the
# (one component longer) sample kernel.
matchfun <- function(d) d$S == d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]
dat$subjects <- droplevels(dat$subjects)

test_that("particle settings are resized to the kernel", {
  mix <- EMC2:::check_mix(c(.1, .3, .3, .3), "adapt")
  expect_length(mix, 3)
  expect_equal(sum(mix), 1)
  expect_length(EMC2:::check_epsilon(c(.5, .6, .7), 10, mix), 2)
  expect_length(EMC2:::check_prop_performance(c(1, 2, 3, 4), "sample"), 4)
  expect_length(EMC2:::check_prop_performance(c(1, 2, 3, 4), "burn"), 3)
  expect_warning(out <- EMC2:::resolve_on_singular(list(max_retries = 2, on_exhausted = "carry_forward")),
                 "on_exhausted")
  expect_equal(out, list(max_retries = 2))
})

test_that("more adapt after the sample-kernel tail", {
  skip_on_cran()
  old <- options(emc.adapt_converge = FALSE); on.exit(options(old))
  RNGkind("L'Ecuyer-CMRG"); set.seed(123)
  des <- design(data = dat, model = LNR, matchfun = matchfun, formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                report_p_vector = FALSE)
  emc <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 2)
  sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
  args <- list(cores_for_chains = 1, particle_factor = 30, step_size = 25, verbose = FALSE)
  emc <- do.call(run_emc, c(list(emc, "preburn", sc("preburn", list(iter = 10))), args))
  emc <- do.call(run_emc, c(list(emc, "burn", sc("burn", list(mean_gd = 2.5))), args))
  emc <- do.call(run_emc, c(list(emc, "adapt", sc("adapt", list(min_unique = 20))), args))
  expect_length(attr(emc[[1]]$samples, "pm_settings")[[1]][[1]]$mix, 4)
  expect_false(is.null(emc[[1]]$sample_kernel))
  n_adapt <- chain_n(emc)[1, "adapt"]
  expect_no_error(more <- do.call(run_emc, c(list(emc, "adapt", sc("adapt", list(iter = 25))), args)))
  expect_gt(chain_n(more)[1, "adapt"], n_adapt)
  expect_no_error(do.call(run_emc, c(list(more, "sample", sc("sample", list(iter = 10))), args)))
  # a sample object sent back to adapt
  samp <- do.call(run_emc, c(list(emc, "sample", sc("sample", list(iter = 10))), args))
  expect_no_error(do.call(run_emc, c(list(samp, "adapt", sc("adapt", list(iter = 25))), args)))
})
