test_that(paste0("pairs_posterior_", Sys.info()[1], Sys.info()[2]), {
  vdiffr::expect_doppelganger(paste0("pairs_plot", Sys.info()[1], Sys.info()[2]), pairs_posterior(samplers_LNR, selection = "variance"))
})
