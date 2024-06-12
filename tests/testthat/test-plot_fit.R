RNGkind("L'Ecuyer-CMRG")
set.seed(123)
pp <- post_predict(samplers_LNR, n_post = 25)
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:3],]
dat$subjects <- droplevels(dat$subjects)
drt <- function(data) diff(tapply(data$rt,data[,c("S")],mean))

test_that("plot_fit", {
  vdiffr::expect_doppelganger("CDFs_plot_fit", plot_fit(dat, pp, factors = c("S", "E"), layout = c(2,3)))
  vdiffr::expect_doppelganger("dens_plot_fit",
                              plot_fit(forstmann, pp, stat=drt,stat_name="Rt difference",
                                                 main=("Left vs Right")))
})
