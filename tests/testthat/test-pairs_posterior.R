test_that("pairs_posterior", {
  vdiffr::expect_doppelganger(paste0("pairs_plot", Sys.info()[1], Sys.info()[2]), pairs_posterior(samplers_LNR, selection = "variance"))
})
