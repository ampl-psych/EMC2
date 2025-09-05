RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR

n_trials <- 10
covariate1 <- rnorm(n_trials*2)
covariate2 <- rnorm(n_trials*2)
# Ensure that NAs are handled correctly in trend
covariate2[1:5] <- NA

trend <- make_trend(par_names = "m", cov_names = list(c("covariate1", "covariate2")), kernels = "exp_incr")

design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)

p_vector <- sampled_pars(design_base, doMap = FALSE)
p_vector[1:6] <- c(-1, 1.5, log(1), log(.2), log(.2), log(.2))

dat <- make_data(p_vector, design_base, n_trials = n_trials, covariates = data.frame(covariate1 = covariate1, covariate2 = covariate2))

LNR2cov <- make_emc(dat, design_base, compress = F, n_chains = 1, type = "single")

test_that("trend", {
  expect_snapshot(init_chains(LNR2cov, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_2types <- make_trend(par_names = c("m", "m_lMd"), cov_names = list(c("covariate1", "covariate2"), "covariate1"), kernels = c("exp_incr", "pow_decr"))

design_base_shared <- design(data = dat,
                             trend = trend_2types,
                             formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                             contrasts = list(lM = ADmat),
                             matchfun = matchfun,
                             model = LNR)

LNR2cov_shared <- make_emc(dat, design_base_shared, compress = FALSE, n_chains = 1, type = "single")

test_that("trend_shared", {
  expect_snapshot(init_chains(LNR2cov_shared, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_premap <- make_trend(
  par_names = c("m", "lMd"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("exp_incr", "poly2"),
  premap = TRUE
)

design_premap <- design(
  data = dat,
  trend = trend_premap,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, lMd.d1 ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_premap <- make_emc(dat, design_premap, compress = FALSE, n_chains = 1, type = "single")
test_that("premap trend works", {
  expect_snapshot(init_chains(LNR_premap, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_pretrans <- make_trend(
  par_names = c("m", "s"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("delta", "exp_decr"),
  premap = FALSE,
  pretransform = TRUE
)

design_pretrans <- design(
  data = dat,
  trend = trend_pretrans,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, s.w ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_pretrans <- make_emc(dat, design_pretrans, compress = FALSE, n_chains = 1, type = "single")

test_that("pretransform trend works", {
  expect_snapshot(init_chains(LNR_pretrans, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_posttrans <- make_trend(
  par_names = c("m", "s"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("pow_decr", "pow_incr"),
  premap = FALSE,
  pretransform = FALSE
)

design_posttrans <- design(
  data = dat,
  trend = trend_posttrans,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, s.w ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)
LNR_posttrans <- make_emc(dat, design_posttrans, compress = FALSE, n_chains = 1, type = "single")

test_that("posttransform trend works", {
  expect_snapshot(init_chains(LNR_posttrans, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_bases <- make_trend(
  par_names = c("m", "s"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("exp_incr", "exp_decr"),
  bases = c("exp_lin", "lin")
)

design_bases <- design(
  data = dat,
  trend = trend_bases,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)
LNR_bases <- make_emc(dat, design_bases, compress = FALSE, n_chains = 1, type = "single")

test_that("different trend base functions work", {
  expect_snapshot(init_chains(LNR_bases, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_poly <- make_trend(
  par_names = c("m", "s"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("poly3", "poly4")
)

design_poly <- design(
  data = dat,
  trend = trend_poly,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_poly <- make_emc(dat, design_poly, compress = FALSE, n_chains = 1, type = "single")

test_that("polynomial trends work", {
  expect_snapshot(init_chains(LNR_poly, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_shared_premap <- make_trend(
  par_names = c("m", "s"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("poly3", "poly4"),
  shared = list(shrd = list("m.d1", "s.d1"))
)

design_shared_premap <- design(
  data = dat,
  trend = trend_shared_premap,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_shared_premap <- make_emc(dat, design_shared_premap, compress = FALSE, n_chains = 1, type = "single")

test_that("share works premap trend", {
  expect_snapshot(init_chains(LNR_shared_premap, particles = 10, cores_per_chain = 1)[[1]]$samples)
})


trend_shared_posttransform <- make_trend(
  par_names = c("m", "s"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("poly3", "poly4"),
  shared = list(shrd = list("m.d1", "s.d1")),
  premap = FALSE,
  pretransform = FALSE
)

design_shared_posttransform <- design(
  data = dat,
  trend = trend_shared_posttransform,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)
LNR_shared_posttransform <- make_emc(dat, design_shared_posttransform, compress = FALSE, n_chains = 1, type = "single")
test_that("share works posttransform trend", {
  expect_snapshot(init_chains(LNR_shared_posttransform, particles = 10, cores_per_chain = 1)[[1]]$samples)
})
