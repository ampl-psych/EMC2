test_that("plot_chains", {
  vdiffr::expect_doppelganger("chain_plots", plot(samples_LNR, selection = "mu"))
})


