test_that("plot_relations", {
  vdiffr::expect_doppelganger(paste0("plot_relations_cred_", Sys.info()[1], Sys.info()[2]), plot_relations(samplers_LNR, only_cred = TRUE, plot_cred = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_relations_mean_", Sys.info()[1], Sys.info()[2]), plot_relations(samplers_LNR, plot_means = TRUE, only_cred = FALSE))
})
