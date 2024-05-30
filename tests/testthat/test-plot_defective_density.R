test_that("density_plot", {
    vdiffr::expect_doppelganger(paste0("density_subject_", Sys.info()[1], Sys.info()[2]),
                      plot_defective_density(forstmann, factors = c("E", "S"), subject = 1,
                                             layout = c(2,3)))
    vdiffr::expect_doppelganger(paste0("density_full_E_", Sys.info()[1], Sys.info()[2]),
                      plot_defective_density(forstmann, factors = c("E"),
                                             correct_fun = function(d) d$R == d$S,
                      layout = c(3,1)))
    vdiffr::expect_doppelganger(paste0("density_full_S_", Sys.info()[1], Sys.info()[2]),
                                plot_defective_density(forstmann, factors = c("S"),
                                                       correct_fun = function(d) d$R == d$S,
                                                       layout = c(1,3)))
  }
)

