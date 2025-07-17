test_that("plot_relations", {
  vdiffr::expect_doppelganger("plot_relations_cred", plot_relations(samples_LNR, only_cred = TRUE, plot_cred = TRUE))
  vdiffr::expect_doppelganger("plot_relations_mean", plot_relations(samples_LNR, plot_means = TRUE, only_cred = FALSE, plot_cred = TRUE))
})
