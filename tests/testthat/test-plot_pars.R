test_that("plot_pars", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  vdiffr::expect_doppelganger("density_plots", plot_pars(samplers_LNR))
  expect_snapshot(plot_pars(samplers_LNR, selection = "variance", do_plot = FALSE))
})
