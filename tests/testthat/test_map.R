dat <- forstmann
dat$covariate <- rnorm(nrow(forstmann))


debug(summary.emc.design)
des <- design(data = dat, formula = list(v ~ covariate, B ~ E, t0 ~ S),
              model = LBA)

test_that("designS3", {
  expect_snapshot(mapped_pars(samples_LNR))
  expect_snapshot(summary(design))
})
