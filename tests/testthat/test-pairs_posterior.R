test_that("pairs_posterior", {
  vdiffr::expect_doppelganger("pairs_plot", pairs_posterior(samplers_LNR, selection = "variance"))
})
