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
covariate2[1:5] <- 0

trend <- make_trend(make_base(type='lin', target_parameter = 'm',
                              kernel=make_kernel(cov_names = c("covariate1", "covariate2"),
                                                 type = "exp_incr")))

design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR, report_p_vector=FALSE)
##mapped_pars(design_base)
p_vector <- sampled_pars(design_base, doMap = FALSE)
p_vector[1:6] <- c(-1, 1.5, log(1), log(.2), log(.2), log(.2))

set.seed(123)
dat <- make_data(p_vector, design_base, n_trials = n_trials,
                 covariates = data.frame(covariate1 = covariate1, covariate2 = covariate2))

LNR2cov <- make_emc(dat, design_base, compress = F, n_chains = 1, type = "single")
# EMC2:::get_pars_oo(p=p_vector, dadm = LNR2cov[[1]]$data[[1]], model = LNR2cov[[1]]$model(), return_all_pars = TRUE)

set.seed(123)
test_that("trend", {
  expect_snapshot(init_chains(LNR2cov, particles = 3, cores_per_chain = 1)[[1]]$samples)
})


# trend_2types <- make_trend(par_names = c("m", "m_lMd"),
#                            cov_names = list(c("covariate1", "covariate2"), "covariate1"),
#                            kernels = c("exp_incr", "pow_decr"))
trend_2types <- make_trend(make_base(type='lin', 'm', make_kernel(type='exp_incr', c('covariate1','covariate2'))),
                           make_base(type='lin', 'm_lMd', make_kernel(type='pow_decr', 'covariate1')))

design_base_shared <- design(data = dat,
                             trend = trend_2types,
                             formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                             contrasts = list(lM = ADmat),
                             matchfun = matchfun,
                             model = LNR, report_p_vector=FALSE)

LNR2cov_shared <- make_emc(dat, design_base_shared, compress = FALSE, n_chains = 1, type = "single")

# set.seed(123)
# p_vector <- init_chains(LNR2cov_shared, particles = 3, cores_per_chain = 1)[[1]]$samples$alpha[,,1]
# EMC2:::get_pars_oo(p=p_vector, dadm = LNR2cov_shared[[1]]$data[[1]], model = LNR2cov_shared[[1]]$model(), return_all_pars = TRUE)

set.seed(123)
test_that("trend_shared", {
  expect_snapshot(init_chains(LNR2cov_shared, particles = 3, cores_per_chain = 1)[[1]]$samples)
})

trend_premap <- make_trend(make_base(type='lin', 'm', make_kernel(type='exp_incr', 'covariate1')),
                           make_base(type='add', 'm_lMd', make_kernel(type='poly2', 'covariate2')))

design_premap <- design(
  data = dat,
  trend = trend_premap,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m_lMd.d1 ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR, report_p_vector=FALSE
)
#mapped_pars(design_premap)

LNR_premap <- make_emc(dat, design_premap, compress = FALSE, n_chains = 1, type = "single")

set.seed(123)
test_that("premap trend works", {
  expect_snapshot(init_chains(LNR_premap, particles = 3, cores_per_chain = 1)[[1]]$samples)
})

trend_pretrans <- make_trend(make_base(type='lin', 'm', make_kernel(type='delta', 'covariate1', at=NULL), phase='pretransform'),
                             make_base(type='lin', 's', make_kernel(type='exp_decr', 'covariate2', at=NULL), phase='pretransform'))

design_pretrans <- design(
  data = dat,
  trend = trend_pretrans,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, s.w ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR, report_p_vector=FALSE
)
#mapped_pars(design_pretrans)

LNR_pretrans <- make_emc(dat, design_pretrans, compress = FALSE, n_chains = 1, type = "single")

set.seed(123)
test_that("pretransform trend works", {
  expect_snapshot(init_chains(LNR_pretrans, particles = 3, cores_per_chain = 1)[[1]]$samples)
})

trend_posttrans <- make_trend(make_base(type='lin','m',make_kernel(type='pow_decr', 'covariate1'), phase='posttransform'),
                              make_base(type='lin','s',make_kernel(type='pow_incr', 'covariate2'), phase='posttransform'))

design_posttrans <- design(
  data = dat,
  trend = trend_posttrans,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, s.w ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR, report_p_vector=FALSE
)
#mapped_pars(design_posttrans)

LNR_posttrans <- make_emc(dat, design_posttrans, compress = FALSE, n_chains = 1, type = "single")
# set.seed(123)
# p_vector <- init_chains(LNR_posttrans, particles = 3, cores_per_chain = 1)[[1]]$samples$alpha[,,1]
# EMC2:::get_pars_oo(p=p_vector, dadm = LNR_posttrans[[1]]$data[[1]], model = LNR_posttrans[[1]]$model(), return_all_pars = TRUE)

set.seed(123)
test_that("posttransform trend works", {
  expect_snapshot(init_chains(LNR_posttrans, particles = 3, cores_per_chain = 1)[[1]]$samples)
})

trend_bases <- make_trend(make_base(type='centered', 'm', make_kernel(type='exp_incr', 'covariate1')),
                          make_base(type='lin', 's', make_kernel(type='exp_decr', 'covariate2')))

design_bases <- design(
  data = dat,
  trend = trend_bases,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR, report_p_vector=FALSE
)
#mapped_pars(design_bases)

LNR_bases <- make_emc(dat, design_bases, compress = FALSE, n_chains = 1, type = "single")

set.seed(123)
test_that("different trend base functions work", {
  expect_snapshot(init_chains(LNR_bases, particles = 3, cores_per_chain = 1)[[1]]$samples)
})


trend_poly <- make_trend(make_base(type='add', 'm', make_kernel(type='poly3', 'covariate1')),
                         make_base(type='add', 's', make_kernel(type='poly4', 'covariate2')))

design_poly <- design(
  data = dat,
  trend = trend_poly,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR, report_p_vector=FALSE
)
#mapped_pars(design_poly)

LNR_poly <- make_emc(dat, design_poly, compress = FALSE, n_chains = 1, type = "single")

set.seed(123)
test_that("polynomial trends work", {
  expect_snapshot(init_chains(LNR_poly, particles = 3, cores_per_chain = 1)[[1]]$samples)
})

## New tests: phase-specific trends and par_input trends

# Three trends on different phases for LNR
trend_phases <- make_trend(make_base(type='lin', 'm', make_kernel(type='lin_incr', 'covariate1'), phase='premap'),
                           make_base(type='lin', 's', make_kernel(type='exp_decr', 'covariate1'), phase='pretransform'),
                           make_base(type='lin', 't0', make_kernel(type='pow_incr', 'covariate2'), phase='posttransform'))



design_phases <- design(
  data = dat,
  trend = trend_phases,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR, report_p_vector=FALSE
)

LNR_phases <- make_emc(dat, design_phases, compress = FALSE, n_chains = 1, type = "single")
# set.seed(123)
# p_vector <- init_chains(LNR_phases, particles = 3, cores_per_chain = 1)[[1]]$samples$alpha[,,1]
# EMC2:::get_pars_oo(p=p_vector, dadm = LNR_phases[[1]]$data[[1]], model = LNR_phases[[1]]$model(), return_all_pars = TRUE)

set.seed(123)
test_that("phase-specific trends work (premap, pretransform, posttransform)", {
  expect_snapshot(init_chains(LNR_phases, particles = 3, cores_per_chain = 1)[[1]]$samples)
})

# Trend where input is another parameter: use t0 as input to a trend on m
trend_par_input <- make_trend(make_base(type='lin', 'm', make_kernel(type='lin_incr', NULL, 't0'), phase="pretransform"))


design_par_input <- design(
  data = dat,
  trend = trend_par_input,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR, report_p_vector=FALSE
)

LNR_par_input <- make_emc(dat, design_par_input, compress = FALSE, n_chains = 1, type = "single")

set.seed(123)
test_that("par_input trend uses t0 as input to m trend", {
  expect_snapshot(init_chains(LNR_par_input, particles = 3, cores_per_chain = 1)[[1]]$samples)
})

trend_shared_premap <- make_trend(make_base(type='add', 'm', kernel=make_kernel(type='poly3', 'covariate1')),
                                  make_base(type='add', 's', kernel=make_kernel(type='poly4', 'covariate2')),
                                  shared = list(shrd=list('m.d1', 's.d1')))


design_shared_premap <- design(
  data = dat,
  trend = trend_shared_premap,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR, report_p_vector=FALSE
)
#mapped_pars(design_shared_premap)

LNR_shared_premap <- make_emc(dat, design_shared_premap, compress = FALSE, n_chains = 1, type = "single")

set.seed(123)
test_that("share works premap trend", {
  expect_snapshot(init_chains(LNR_shared_premap, particles = 3, cores_per_chain = 1)[[1]]$samples)
})


trend_shared_posttransform <- make_trend(make_base(type='add', 'm', kernel=make_kernel(type='poly3', 'covariate1'), phase='posttransform'),
                                  make_base(type='add', 's', kernel=make_kernel(type='poly4', 'covariate2'), phase='posttransform'),
                                  shared = list(shrd=list('m.d1', 's.d1')))


design_shared_posttransform <- design(
  data = dat,
  trend = trend_shared_posttransform,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR, report_p_vector=FALSE
)
LNR_shared_posttransform <- make_emc(dat, design_shared_posttransform, compress = FALSE, n_chains = 1, type = "single")

set.seed(123)
test_that("share works posttransform trend", {
  expect_snapshot(init_chains(LNR_shared_posttransform, particles = 3, cores_per_chain = 1)[[1]]$samples)
})


n_trials <- 10

# Trend uses behavioral covariate rt with delta kernel (forces trial-wise path)
trend_cond <- make_trend(make_base(type='lin', 'm', make_kernel(type='delta', 'trial2', at=NULL), phase='pretransform'))

design_cond <- design(
  factors = list(subjects = 1, S = 1:2),
  Rlevels = 1:2,
  covariates = c("trial2"),
  matchfun = matchfun,
  trend = trend_cond,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  model = LNR, report_p_vector=FALSE
)

p_vec <- sampled_pars(design_cond, doMap = FALSE)
# Set basic LNR params and trend params (base weight w, delta q0/alpha)
p_vec[c("m", "s", "t0")] <- c(-0.5, log(0.3), log(0.2))
p_vec[c("m.w", "m.q0", "m.alpha")] <- c(0.5, 0.0, qnorm(0.2))



set.seed(123)
test_that("trend_conditional", {
  expect_snapshot(attributes(make_data(
    p_vec, design_cond, n_trials = n_trials,
    conditional_on_data = FALSE,                 # force unconditional-on-data stepping
    return_trialwise_parameters = TRUE
  )))
})

# Trend uses multiple same parameters and multiple inputs
trend_mult <- make_trend(make_base(type='lin', 'm', make_kernel(type='exp_incr', 'covariate1'), phase='pretransform'),
                         make_base(type='lin', 'm', make_kernel(type='delta', c('covariate1', 'covariate2'), par_input = "t0"), phase='pretransform'))

design_mult <- design(
  factors = list(subjects = 1, S = 1:2),
  Rlevels = 1:2,
  covariates = c("covariate1", "covariate2"),
  matchfun = matchfun,
  trend = trend_mult,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  model = LNR, report_p_vector=FALSE
)

LNR_multi <- make_emc(dat, design_mult, compress = FALSE, n_chains = 1, type = "single")

set.seed(123)
test_that("trend_multiple", {
  expect_snapshot(init_chains(LNR_multi, particles = 3, cores_per_chain = 1)[[1]]$samples)
})

k1 <- make_kernel(type='delta', c('covariate1', 'covariate2'))
trend <- make_trend(make_base(type='lin', 'm', kernel=k1,
                              coding=list('map1'=function(dadm, covs) {
                                d <- matrix(rnorm(nrow(dadm)*2), ncol=2)
                                colnames(d) <- covs
                                d
                              })),
                    make_base(type='lin', 'm', kernel=k1,
                              coding=list('map2'=function(dadm, covs) {
                                d <- matrix(rnorm(nrow(dadm)*2), ncol=2)
                                colnames(d) <- covs
                                d
                              })))

# covariate coding
design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR, report_p_vector = FALSE)
set.seed(123)
covariate1 <- rnorm(n_trials*2)
covariate2 <- rnorm(n_trials*2)

p_vector <- sampled_pars(design_base, doMap = FALSE)
p_vector[1:6] <- c(-1, 1.5, log(1), log(.2), log(.2), log(.2))

set.seed(123)
dat <- make_data(p_vector, design_base, n_trials = n_trials, covariates = data.frame(covariate1 = covariate1, covariate2 = covariate2))

LNR_covmap <- make_emc(dat, design_base, compress = FALSE, n_chains = 1, type = "single")


set.seed(123)
test_that("trend_covmap", {
  expect_snapshot(init_chains(LNR_covmap, particles = 3, cores_per_chain = 1)[[1]]$samples)
})



# # Manual test of covmaps --------------------------------------------------
# trend <- make_trend(par_names = "m", cov_names = list(c("covariate1", "covariate2")),
#                     kernels = "delta", ffill_na=TRUE,
#                     maps=list('map1'=function(dadm, covs) {
#                       d <- matrix(1, nrow=nrow(dadm), ncol=2)
#                       colnames(d) <- covs
#                       d
#                     },
#                     'map2'=function(dadm, covs) {
#                       d <- matrix(-0.5, nrow=nrow(dadm), ncol=2)
#                       colnames(d) <- covs
#                       d
#                     }))
# design_base <- design(factors = list(subjects = 1, S = 1:2),
#                       Rlevels = 1:2,
#                       covariates = c('covariate1', 'covariate2'),
#                       matchfun = matchfun,
#                       trend = trend,
#                       formula = list(m ~ lM, s ~ 1, t0 ~ 1),
#                       contrasts = list(lM = ADmat),
#                       model = LNR, report_p_vector=FALSE)
# covariate1 <- rep(1, n_trials*2)
# covariate2 <- rep(0.25, n_trials*2)
#
# p_vector <- sampled_pars(design_base, doMap = FALSE)
# p_vector[1:6] <- c(-1, 1.5, log(1), log(.2), log(.2), log(.2))
#
# dat <- make_data(p_vector, design_base, n_trials = n_trials, covariates = data.frame(covariate1 = covariate1, covariate2 = covariate2),
#                  return_trialwise_parameters=TRUE)
#
# LNR_covmap <- make_emc(dat, design_base, compress = FALSE, n_chains = 1, type = "single")
#
# trpars <- attr(dat, 'trialwise_parameters')
# covmaps <- attr(LNR_covmap[[1]]$data[[1]], 'covariate_maps')
# head(trpars)
# head(covmaps[[1]])
# head(covmaps[[2]])
#
# m_without_lMd <- p_vector[1]+ p_vector[['m.w_map1']]*trpars[,4]*covmaps[[1]][,1] + p_vector[['m.w_map1']]*trpars[,5]*covmaps[[1]][,2] +
#                               p_vector[['m.w_map2']]*trpars[,4]*covmaps[[2]][,1] + p_vector[['m.w_map2']]*trpars[,5]*covmaps[[2]][,2]
# m <- m_without_lMd
# m[LNR_covmap[[1]]$data[[1]]$lM==TRUE] <- m[LNR_covmap[[1]]$data[[1]]$lM==TRUE] + 0.5*p_vector[[2]]
# m[LNR_covmap[[1]]$data[[1]]$lM==FALSE] <- m[LNR_covmap[[1]]$data[[1]]$lM==FALSE] - 0.5*p_vector[[2]]
# all(round(m,5) == round(trpars[,1], 5))
#
