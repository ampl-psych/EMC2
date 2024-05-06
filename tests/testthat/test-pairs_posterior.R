test_that("pairs_posterior", {
  vdiffr::expect_doppelganger("pairs plot", pairs_posterior(samplers_LNR, selection = "variance"))
})
