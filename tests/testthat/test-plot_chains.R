test_that("plot_chains", {
  vdiffr::expect_doppelganger("chain_plots", plot_chains(samples_LNR, selection = "mu"))
})


