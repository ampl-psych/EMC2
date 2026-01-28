prior <- get_prior(samples_LNR)
design <- get_design(samples_LNR)


RNGkind("L'Ecuyer-CMRG")
set.seed(123)
test_that("priorS3", {
  expect_snapshot(prior)
  expect_snapshot(summary(prior))
  vdiffr::expect_doppelganger("prior_plot", plot(prior, map = FALSE, N = 1e2))
})

test_that("designS3", {
  expect_snapshot(design)
  expect_snapshot(summary(design))
})
