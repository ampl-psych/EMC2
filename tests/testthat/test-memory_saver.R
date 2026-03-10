test_that("predict requires data for memory_saver fits", {
  dat <- subset(forstmann, subjects == levels(forstmann$subjects)[1])
  dat <- dat[seq_len(100), , drop = FALSE]

  design_in <- design(
    data = dat,
    model = DDM,
    formula = list(v ~ 0 + S, a ~ E, t0 ~ 1, s ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1),
    constants = c(s = log(1))
  )

  emc_mem <- make_emc(
    dat,
    design_in,
    type = "single",
    compress = FALSE,
    memory_saver = TRUE
  )

  expect_error(
    predict(emc_mem, n_post = 1, n_cores = 1),
    "predict\\(\\) requires `data` when the model was created with memory_saver = TRUE"
  )
  expect_error(
    predict(emc_mem, data = dat, n_post = 1, n_cores = 1),
    "Run fit to collect MCMC samples before using this function"
  )
})
