test_that("plot_chains", {
  vdiffr::expect_doppelganger("chain_plots", plot_chains(samplers_LNR, selection = "mu"))
})


