RNGkind("L'Ecuyer-CMRG")
set.seed(123)
# Generate some prior predictives
dat <- get_data(samples_LNR)
pp <- predict(get_prior(samples_LNR), n_post = 25, data = dat)
drt <- function(data) diff(tapply(data$rt,data[,c("S")],mean))

test_that("plot_fit", {
  vdiffr::expect_doppelganger("CDF_plot_fit", plot_cdf(dat, prior_predict = pp, factors = c("S", "E"), layout = c(2,3), use_lim = c('real')))
  vdiffr::expect_doppelganger("dens_plot_fit", plot_density(dat, prior_predict = pp, defective_factor = "E"))
  vdiffr::expect_doppelganger("stat_plot_fit", plot_stat(samples_LNR, to_plot = c('real', 'posterior'), stat_fun=drt,stat_name="RT diff Left vs Right",
                                                         n_post = 10))
})
