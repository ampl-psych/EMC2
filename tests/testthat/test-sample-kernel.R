test_that("fit continuation appends draws and preserves the sample kernel", {
  withr::local_seed(123)
  withr::local_options(emc.adapt_converge = FALSE)
  emc <- fit_tiny_lnr(tiny_lnr())
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

test_that("run_emc() with an iteration-only criterion draws one block and keeps the kernel", {
  withr::local_seed(123)
  withr::local_options(emc.adapt_converge = FALSE)
  emc <- fit_tiny_lnr(tiny_lnr())
  one <- run_emc(emc, "sample", stop_criteria = list(iter = 20), step_size = 20, max_tries = 1,
                 cores_for_chains = 1, verbose = FALSE, particle_factor = 10, trim = FALSE)
  two <- run_emc(one, "sample", stop_criteria = list(iter = 20), step_size = 20, max_tries = 1,
                 cores_for_chains = 1, verbose = FALSE, particle_factor = 10, trim = FALSE)
  expect_equal(unname(chain_n(one)[, "sample"]), rep(40, 2))
  expect_equal(unname(chain_n(two)[, "sample"]), rep(60, 2))
  expect_identical(two[[1]]$sample_kernel, emc[[1]]$sample_kernel)
  for(ch in seq_along(emc)){
    fields <- c("step", "loc", "directions")
    expect_identical(attr(two[[ch]]$samples, "scale_move")[fields], attr(emc[[ch]]$samples, "scale_move")[fields])
    expect_equal(filter_obj(two[[ch]]$samples$theta_mu, seq_len(emc[[ch]]$samples$idx)), emc[[ch]]$samples$theta_mu)
  }
})
