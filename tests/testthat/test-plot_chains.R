test_that("plot_chains", {
  vdiffr::expect_doppelganger(paste0("chain_plots_", Sys.info()[1], Sys.info()[2]), plot_chains(samplers_LNR, selection = "mu"))
})


