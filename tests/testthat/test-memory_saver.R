test_that("memory_saver fits predict from supplied data", {
  withr::local_seed(123)
  dat <- data.frame(subjects = factor(rep(1:2, each = 4)), R = factor(1),
                    rt = .2 + exp(c(-.7, -.2, .1, .4, .3, .8, -.4, .6)))
  design_in <- design(data = dat, model = LNR, formula = list(m ~ 1),
                      constants = c(s = 0, t0 = log(.2)), report_p_vector = FALSE)
  emc_mem <- make_emc(dat, design_in, type = "single", n_chains = 2,
                      compress = FALSE, memory_saver = TRUE)
  criteria <- list(preburn = list(iter = 10), burn = list(iter = 10),
                   adapt = list(iter = 20, min_unique = 5), sample = list(iter = 20))
  emc_mem <- fit(emc_mem, cores_for_chains = 1, stop_criteria = criteria,
                 particle_factor = 10, step_size = 20, trim = FALSE, verbose = FALSE)

  expect_error(
    predict(emc_mem, n_post = 1, n_cores = 1),
    "predict\\(\\) requires `data` when the model was created with memory_saver = TRUE"
  )
  post <- predict(emc_mem, data = dat, n_post = 2, n_cores = 1)
  expect_equal(nrow(post), 2 * nrow(dat))
  expect_equal(table(post$postn, post$subjects), matrix(4L, 2, 2), ignore_attr = TRUE)
  expect_true(all(is.finite(post$rt) & post$rt > .2))
  expect_equal(levels(post$R), levels(dat$R))
})
