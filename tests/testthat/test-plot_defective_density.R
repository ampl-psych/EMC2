test_that("plot_defective_density", {
    vdiffr::expect_doppelganger("density_subject",
                      plot_defective_density(forstmann, factors = c("E", "S"), subject = 1,
                                             mfcol = TRUE, layout = c(2,3)))
    vdiffr::expect_doppelganger("density_full_E",
                      plot_defective_density(forstmann, factors = c("E"),
                                             correct_fun = function(d) d$R == d$S,
                      layout = c(3,1)))
    vdiffr::expect_doppelganger("density_full_S",
                                plot_defective_density(forstmann, factors = c("S"),
                                                       correct_fun = function(d) d$R == d$S,
                                                       layout = c(1,3)))
  }
)

