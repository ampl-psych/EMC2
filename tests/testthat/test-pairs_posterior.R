test_that(paste0("pairs_posterior_", Sys.info()[1]), {
  vdiffr::expect_doppelganger(paste0("pairs_plot", Sys.info()[1]), pairs_posterior(samplers_LNR, selection = "variance"))
})
