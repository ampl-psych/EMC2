test_that("the sample-stage proposals and tuning stay fixed", {
  skip_on_cran()
  skip_on_os("windows")
  set.seed(11, kind = "Mersenne-Twister")
  dat <- droplevels(forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:4], ])
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR,
                formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                contrasts = list(m = list(lM = ADmat)), report_p_vector = FALSE)
  emc <- make_emc(dat, des, type = "standard", n_chains = 2, rt_resolution = .05)
  criteria <- list(preburn = list(iter = 10), burn = list(mean_gd = 2.5),
                   adapt = list(min_unique = 20), sample = list(iter = 20))
  set.seed(123, kind = "L'Ecuyer-CMRG")
  emc <- fit(emc, cores_for_chains = 1, stop_criteria = criteria, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  set.seed(124)
  continued <- fit(emc, cores_for_chains = 1, verbose = FALSE, particle_factor = 20,
                   step_size = 20, stop_criteria = list(sample = list(iter = 40)))
  expect_length(emc[[1]]$sample_kernel$chains, length(emc))
  expect_identical(continued[[1]]$sample_kernel, emc[[1]]$sample_kernel)
  for (ch in seq_along(emc)) {
    before <- attr(emc[[ch]]$samples, "pm_settings")
    after <- attr(continued[[ch]]$samples, "pm_settings")
    for (s in seq_along(before)) {
      fields <- c("epsilon", "mix", "n_particles")
      expect_identical(after[[s]][[1]][fields], before[[s]][[1]][fields])
    }
    fields <- c("step", "loc", "directions")
    expect_identical(attr(continued[[ch]]$samples, "scale_move")[fields],
                     attr(emc[[ch]]$samples, "scale_move")[fields])
  }
})

# The joint location move's directions are built once per stage: kept across
# the steps of burn, rebuilt when the stage changes, so that loc adapts
# continuously instead of restarting at every step.
test_that("joint translation directions persist within a stage", {
  skip_on_cran()
  old <- options(emc.adapt_converge = FALSE); on.exit(options(old))
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR,
                formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
  sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
  set.seed(7); emc <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 2)
  emc <- run_emc(emc, "preburn", sc("preburn", list(iter = 10)), cores_for_chains = 1, step_size = 25, verbose = FALSE)
  set.seed(1); one <- run_emc(emc, "burn", sc("burn", list(iter = 25)), cores_for_chains = 1, step_size = 25, verbose = FALSE)
  set.seed(1); three <- run_emc(emc, "burn", sc("burn", list(iter = 75)), cores_for_chains = 1, step_size = 25, verbose = FALSE)
  d1 <- attr(one[[1]]$samples, "scale_move")$directions
  d3 <- attr(three[[1]]$samples, "scale_move")$directions
  expect_equal(dim(d1), c(4, 4))
  expect_identical(d3, d1)
  adapt <- run_emc(three, "adapt", sc("adapt", list(iter = 25)), cores_for_chains = 1, step_size = 25, verbose = FALSE)
  expect_false(identical(attr(adapt[[1]]$samples, "scale_move")$directions, d3))
})
