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
    fields <- c("step", "loc")
    expect_identical(attr(continued[[ch]]$samples, "scale_move")[fields],
                     attr(emc[[ch]]$samples, "scale_move")[fields])
  }
})
