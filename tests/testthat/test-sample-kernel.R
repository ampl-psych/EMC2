test_that("fit continuation appends draws and preserves the sample kernel", {
  withr::local_seed(123)
  withr::local_options(emc.adapt_converge = FALSE)
  dat <- data.frame(subjects = factor(rep(1:4, each = 2)), R = factor(1),
                    rt = .2 + exp(c(-.7, -.2, .1, .4, .3, .8, -.4, .6)))
  des <- design(data = dat, model = LNR, formula = list(m ~ 1),
                constants = c(s = 0, t0 = log(.2)), report_p_vector = FALSE)
  emc <- make_emc(dat, des, type = "standard", n_chains = 2)
  criteria <- list(preburn = list(iter = 10), burn = list(iter = 10),
                   adapt = list(iter = 20, min_unique = 5), sample = list(iter = 20))
  emc <- fit(emc, cores_for_chains = 1, stop_criteria = criteria, verbose = FALSE,
             particle_factor = 10, step_size = 20, trim = FALSE)
  path <- tempfile(fileext = ".rds")
  withr::defer(unlink(path))
  saveRDS(emc, path)
  set.seed(124)
  continued <- fit(readRDS(path), cores_for_chains = 1, verbose = FALSE, particle_factor = 10,
                   step_size = 20, stop_criteria = list(sample = list(iter = 20)), trim = FALSE)
  expect_length(emc[[1]]$sample_kernel$chains, length(emc))
  expect_identical(continued[[1]]$sample_kernel, emc[[1]]$sample_kernel)
  for (ch in seq_along(emc)) {
    before_samples <- emc[[ch]]$samples
    after_samples <- continued[[ch]]$samples
    expect_equal(sum(before_samples$stage == "sample"), 20)
    expect_equal(sum(after_samples$stage == "sample"), 40)
    expect_equal(after_samples$idx, before_samples$idx + 20)
    for (field in c("alpha", "subj_ll", "theta_mu", "theta_var", "a_half")) {
      expect_equal(filter_obj(after_samples[[field]], seq_len(before_samples$idx)),
                   before_samples[[field]])
      expect_true(all(is.finite(after_samples[[field]])))
    }
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
