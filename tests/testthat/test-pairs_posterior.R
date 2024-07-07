test_that("pairs_posterior", {
  vdiffr::expect_doppelganger("pairs_plot", pairs_posterior(samples_LNR, selection = "sigma2"))
})

