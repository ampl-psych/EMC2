test_that("plot_relations", {
  vdiffr::expect_doppelganger("plot_relations_cred", plot_relations(samplers_LNR, only_cred = TRUE, plot_cred = TRUE))
  vdiffr::expect_doppelganger("plot_relations_means", plot_relations(samplers_LNR, plot_means = TRUE, only_cred = FALSE))
})
