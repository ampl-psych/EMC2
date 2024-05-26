test_that(paste0("plot_pars_", Sys.info()[1], Sys.info()[2]), {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  vdiffr::expect_doppelganger(paste0("density_plots_", Sys.info()[1], Sys.info()[2]), plot_pars(samplers_LNR))
  expect_snapshot(plot_pars(samplers_LNR, selection = "variance", do_plot = FALSE))
})
