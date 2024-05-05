
test_that("plot_defective_density", {
    expect_doppelganger("density_subject",
                      plot_defective_density(forstmann, factors = c("E", "S"), subject = 1,
                                             mfcol = TRUE, layout = c(2,3)))
    expect_doppelganger("density_full",
                      plot_defective_density(forstmann, factors = c("E", "S"),
                                             correct_fun = function(d) d$R == d$S))
  }
)

