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

trend <- make_trend(par_names = "m",
                    cov_names = list(c("covariate1", "covariate2")),
                    kernels = "exp_incr")

design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)
##mapped_pars(design_base)
p_vector <- sampled_pars(design_base, doMap = FALSE)
p_vector[1:6] <- c(-1, 1.5, log(1), log(.2), log(.2), log(.2))

dat <- make_data(p_vector, design_base, n_trials = n_trials,
                 covariates = data.frame(covariate1 = covariate1, covariate2 = covariate2))

LNR2cov <- make_emc(dat, design_base, compress = F, n_chains = 1, type = "single")

test_that("trend", {
  expect_snapshot(init_chains(LNR2cov, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_2types <- make_trend(par_names = c("m", "m_lMd"),
                           cov_names = list(c("covariate1", "covariate2"), "covariate1"),
                           kernels = c("exp_incr", "pow_decr"))

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
  phase = "premap"
)

design_premap <- design(
  data = dat,
  trend = trend_premap,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, lMd.d1 ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)
#mapped_pars(design_premap)

LNR_premap <- make_emc(dat, design_premap, compress = FALSE, n_chains = 1, type = "single")
test_that("premap trend works", {
  expect_snapshot(init_chains(LNR_premap, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_pretrans <- make_trend(
  par_names = c("m", "s"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("delta", "exp_decr"),
  phase = "pretransform", at=NULL  # should be at='lR' in realistic cases but just for this test, set to NULL
)

design_pretrans <- design(
  data = dat,
  trend = trend_pretrans,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, s.w ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)
#mapped_pars(design_pretrans)

LNR_pretrans <- make_emc(dat, design_pretrans, compress = FALSE, n_chains = 1, type = "single")

test_that("pretransform trend works", {
  expect_snapshot(init_chains(LNR_pretrans, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

trend_posttrans <- make_trend(
  par_names = c("m", "s"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("pow_decr", "pow_incr"),
  phase = "posttransform"
)

design_posttrans <- design(
  data = dat,
  trend = trend_posttrans,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, s.w ~ lR),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)
#mapped_pars(design_posttrans)

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
#mapped_pars(design_bases)

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
#mapped_pars(design_poly)

LNR_poly <- make_emc(dat, design_poly, compress = FALSE, n_chains = 1, type = "single")

test_that("polynomial trends work", {
  expect_snapshot(init_chains(LNR_poly, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

## New tests: phase-specific trends and par_input trends

# Three trends on different phases for LNR
trend_phases <- make_trend(
  par_names = c("m", "s", "t0"),
  cov_names = list("covariate1", "covariate1", "covariate2"),
  kernels = c("lin_incr", "exp_decr", "pow_incr"),
  phase = c("premap", "pretransform", "posttransform")
)

design_phases <- design(
  data = dat,
  trend = trend_phases,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_phases <- make_emc(dat, design_phases, compress = FALSE, n_chains = 1, type = "single")

test_that("phase-specific trends work (premap, pretransform, posttransform)", {
  expect_snapshot(init_chains(LNR_phases, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

# Trend where input is another parameter: use t0 as input to a trend on m
trend_par_input <- make_trend(
  par_names = c("m"),
  cov_names = NULL,
  kernels = c("lin_incr"),
  par_input = list(c("t0")),
  phase = "pretransform"
)

design_par_input <- design(
  data = dat,
  trend = trend_par_input,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_par_input <- make_emc(dat, design_par_input, compress = FALSE, n_chains = 1, type = "single")

test_that("par_input trend uses t0 as input to m trend", {
  expect_snapshot(init_chains(LNR_par_input, particles = 10, cores_per_chain = 1)[[1]]$samples)
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
#mapped_pars(design_shared_premap)

LNR_shared_premap <- make_emc(dat, design_shared_premap, compress = FALSE, n_chains = 1, type = "single")

test_that("share works premap trend", {
  expect_snapshot(init_chains(LNR_shared_premap, particles = 10, cores_per_chain = 1)[[1]]$samples)
})


trend_shared_posttransform <- make_trend(
  par_names = c("m", "s"),
  cov_names = list("covariate1", "covariate2"),
  kernels = c("poly3", "poly4"),
  shared = list(shrd = list("m.d1", "s.d1")),
  phase = "posttransform"
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


n_trials <- 10

# Trend uses behavioral covariate rt with delta kernel (forces trial-wise path)
trend_cond <- make_trend(
  par_names = "m",
  cov_names = "trial2",
  kernels = "delta",
  phase = "pretransform", at=NULL  # should be at='lR' for realistic cases
)

design_cond <- design(
  factors = list(subjects = 1, S = 1:2),
  Rlevels = 1:2,
  covariates = c("trial2"),
  matchfun = matchfun,
  trend = trend_cond,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  model = LNR
)

p_vec <- sampled_pars(design_cond, doMap = FALSE)
# Set basic LNR params and trend params (base weight w, delta q0/alpha)
p_vec[c("m", "s", "t0")] <- c(-0.5, log(0.3), log(0.2))
p_vec[c("m.w", "m.q0", "m.alpha")] <- c(0.5, 0.0, qnorm(0.2))



test_that("trend_conditional", {
  expect_snapshot(attributes(make_data(
    p_vec, design_cond, n_trials = n_trials,
    conditional_on_data = FALSE,                 # force unconditional-on-data stepping
    return_trialwise_parameters = TRUE
  )))
})

# Trend uses multiple same parameters and multiple inputs
trend_mult <- make_trend(
  par_names = c("m", "m"),
  cov_names = list("covariate1", c("covariate1", "covariate2")),
  par_input = list(NULL, "t0"),
  kernels = c("exp_incr", "delta"),
  phase = "pretransform",
#  at = "lR",
  ffill_na = FALSE
)

design_mult <- design(
  factors = list(subjects = 1, S = 1:2),
  Rlevels = 1:2,
  covariates = c("trial2", "trial3"),
  matchfun = matchfun,
  trend = trend_mult,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(lM = ADmat),
  model = LNR
)

LNR_multi <- make_emc(dat, design_mult, compress = FALSE, n_chains = 1, type = "single")

test_that("trend_multiple", {
  expect_snapshot(init_chains(LNR_multi, particles = 10, cores_per_chain = 1)[[1]]$samples)
})



# Test handling of NA -----------------------------------------------------
trend <- make_trend(par_names = "m", cov_names = list(c("covariate1", "covariate2")), kernels = "exp_incr", ffill_na=FALSE)
design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)
##mapped_pars(design_base)
p_vector <- sampled_pars(design_base, doMap = FALSE)
p_vector[1:6] <- c(-1, 1.5, log(1), log(.2), log(.2), log(.2))

covariate1 <- rnorm(n_trials*2)
covariate2 <- rnorm(n_trials*2)
# Ensure that NAs are handled correctly in trend
covariate2[c(1:5, 8)] <- NA

dat <- make_data(p_vector, design_base, n_trials = n_trials, covariates = data.frame(covariate1 = covariate1, covariate2 = covariate2), return_trialwise_parameters=TRUE)

test_that("trend_ffillnafalse", {
  expect_snapshot(attr(dat, 'trialwise_parameters'))
})

# with ffill_na
trend <- make_trend(par_names = "m", cov_names = list(c("covariate1", "covariate2")),
                    kernels = "exp_incr", ffill_na=TRUE, at=NULL)
design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)
##mapped_pars(design_base)
p_vector <- sampled_pars(design_base, doMap = FALSE)
p_vector[1:6] <- c(-1, 1.5, log(1), log(.2), log(.2), log(.2))

covariate1 <- rnorm(n_trials*2)
covariate2 <- rnorm(n_trials*2)
# Ensure that NAs are handled correctly in trend
covariate2[c(2:5, 8)] <- NA

dat <- make_data(p_vector, design_base, n_trials = n_trials, covariates = data.frame(covariate1 = covariate1, covariate2 = covariate2), return_trialwise_parameters=TRUE)
test_that("trend_ffillnatrue", {
  expect_snapshot(attr(dat, 'trialwise_parameters'))
})


# Delta rule - always set initial trial to q0
trend <- make_trend(par_names = "m", cov_names = list(c("covariate1", "covariate2")),
                    kernels = "delta", ffill_na=TRUE)
design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)
##mapped_pars(design_base)
p_vector <- sampled_pars(design_base, doMap = FALSE)
p_vector[1:7] <- c(-1, 1.5, log(1), log(.2), 1, .5, qnorm(.2))

covariate1 <- rnorm(n_trials*2)
covariate2 <- rnorm(n_trials*2)
# Ensure that NAs are handled correctly in trend
covariate2[c(1:5, 8)] <- NA

dat <- make_data(p_vector, design_base, n_trials = n_trials, covariates = data.frame(covariate1 = covariate1, covariate2 = covariate2), return_trialwise_parameters=TRUE)
test_that("trend_ffillnatrue_delta", {
  expect_snapshot(attr(dat, 'trialwise_parameters'))
})



# ##
# trend <- make_trend(par_names = "m", cov_names = list(c("covariate1")),
#                     kernels = "delta", ffill_na=TRUE)
# design_base <- design(factors = list(subjects = 1, S = 1:2),
#                       Rlevels = 1:2,
#                       covariates = c('covariate1'),
#                       matchfun = matchfun,
#                       trend = trend,
#                       formula = list(m ~ lM, s ~ 1, t0 ~ 1),
#                       contrasts = list(lM = ADmat),
#                       model = LNR)
# ##mapped_pars(design_base)
# p_vector <- sampled_pars(design_base, doMap = FALSE)
#
# p_vector[1:7] <- c(-1, 1.5, log(1), log(.2), 1, .5, qnorm(.2))
#
# covariate1 <- c(NA, 1, NA, NA, 1, NA, NA, NA, NA, NA)#, 1, NA, 1, 1, NA, rep(NA, 10))
#
# #debug(EMC2:::run_kernel)
# dat <- make_data(p_vector, design_base, n_trials = 5, covariates = data.frame(covariate1 = covariate1), return_trialwise_parameters=TRUE)
# # cbind(attr(dat, 'trialwise_parameters'), rep(covariate1, each=2))
#
#
# ## Test code for comparing the output of run_trend between R and Rcpp
# # signature:
# # NumericVector run_trend_rcpp(DataFrame data, List trend, NumericVector param, NumericMatrix trend_pars, NumericMatrix pars_full) {
# emc <- make_emc(dat, design_base, type='single')
# dadm <- emc[[1]]$data[[1]]
# trend <- emc[[1]]$model()$trend$m
# param <- rep(0,  nrow(dadm))
# trend_pars <- matrix(rep(c(p_vector[5:7]), each=nrow(dadm)), ncol=3, byrow=FALSE)
# trend_pars[,3] <- pnorm(trend_pars[,3])
# pars_full <- matrix(rep(p_vector, each=nrow(dadm)), ncol=length(p_vector), byrow=FALSE)
# pars_full[,ncol(pars_full)] <- pnorm(pars_full[,ncol(pars_full)])
# cv1_updated <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars, pars_full = pars_full, return_kernel = TRUE)
#
#
# old_behavior <- c(0,0,0.5,0.5,0,0,0,0,.6,.6,0,0,0,0,0,0,0,0,0,0)
# ffill_behavior <- c(rep(.5, 8), rep(.6, 12))
# bfill_behavior <- c(rep(.5, 4), rep(.6, 6), rep(0, 10))
# cbind(dadm[,c('trials', 'covariate1')], target=cv1_updated, old_behavior=old_behavior, ffill=ffill_behavior, bfill=bfill_behavior)
#
# # 3 points of failure in the future: 1. hard-coded that the first trial must be q0; 2. "ffill" is actually a "bfill", 3. hard-coded that the last trial must be updated
#
#


# # library(EMC2)
# ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# # We also define a match function for lM
# matchfun=function(d)d$S==d$lR
# trend <- make_trend(par_names = "m", cov_names = list(c("covariate1")),
#                     kernels = "delta", ffill_na=TRUE)
# design_base <- design(factors = list(subjects = 1, S = 1:2),
#                       Rlevels = 1:2,
#                       covariates = c('covariate1'),
#                       matchfun = matchfun,
#                       trend = trend,
#                       formula = list(m ~ lM, s ~ 1, t0 ~ 1),
#                       contrasts = list(lM = ADmat),
#                       model = LNR)
# ##mapped_pars(design_base)
# p_vector <- sampled_pars(design_base, doMap = FALSE)
#
# p_vector[1:7] <- c(-1, 1.5, log(1), log(.2), 1, .5, qnorm(.2))
#
# covariate1 <- c(NA, 1, NA, NA, 1, NA, NA, NA, NA, NA)#, 1, NA, 1, 1, NA, rep(NA, 10))
#
# #debug(EMC2:::run_kernel)
# dat <- make_data(p_vector, design_base, n_trials = 5, covariates = data.frame(covariate1 = covariate1), return_trialwise_parameters=TRUE)
# # cbind(attr(dat, 'trialwise_parameters'), rep(covariate1, each=2))
#
#
# ## Test code for comparing the output of run_trend between R and Rcpp
# # signature:
# # NumericVector run_trend_rcpp(DataFrame data, List trend, NumericVector param, NumericMatrix trend_pars, NumericMatrix pars_full) {
# emc <- make_emc(dat, design_base, type='single')
# dadm <- emc[[1]]$data[[1]]
# trend <- emc[[1]]$model()$trend$m
# param <- rep(0,  nrow(dadm))
# trend_pars <- matrix(rep(c(p_vector[5:7]), each=nrow(dadm)), ncol=3, byrow=FALSE)
# trend_pars[,3] <- pnorm(trend_pars[,3])
# pars_full <- matrix(rep(p_vector, each=nrow(dadm)), ncol=length(p_vector), byrow=FALSE)
# pars_full[,ncol(pars_full)] <- pnorm(pars_full[,ncol(pars_full)])
# cv1_updated <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
#                                      pars_full = pars_full, return_kernel = TRUE, use_ptrs = FALSE)
# cv1_updated
#
# cv1_updated_ptr <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
#                                          pars_full = pars_full, return_kernel = TRUE, use_ptrs = TRUE)
#
# cv1_updated_noptr <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
#                                            pars_full = pars_full, return_kernel = TRUE, use_ptrs = FALSE)
#
#
# n_trials <- 10000
# covariate1 <- rnorm(n_trials)
# covariate2 <- rnorm(n_trials)
# covariate3 <- rnorm(n_trials)
# covariate2[sample(1:n_trials, 50)] <- NA
#
# trend <- make_trend(par_names = "m", cov_names = list(c("covariate1","covariate2","covariate3")),
#                     kernels = "delta", ffill_na=TRUE)
# design_base <- design(factors = list(subjects = 1, S = 1:2),
#                       Rlevels = 1:2,
#                       covariates = c("covariate1","covariate2","covariate3"),
#                       matchfun = matchfun,
#                       trend = trend,
#                       formula = list(m ~ lM, s ~ 1, t0 ~ 1),
#                       contrasts = list(lM = ADmat),
#                       model = LNR)
# dat <- make_data(p_vector, design_base, n_trials = n_trials,
#                  covariates = data.frame(covariate1 = covariate1,
#                                          covariate2 = covariate2,
#                                          covariate3 = covariate3
#                  ), return_trialwise_parameters=TRUE)
# emc <- make_emc(dat, design_base, type='single')
#
#
# emc <- make_emc(dat, design_base, type='single')
# dadm <- emc[[1]]$data[[1]]
# trend <- emc[[1]]$model()$trend$m
# param <- rep(0,  nrow(dadm))
# trend_pars <- matrix(rep(c(p_vector[5:7]), each=nrow(dadm)), ncol=3, byrow=FALSE)
# trend_pars[,3] <- pnorm(trend_pars[,3])
# pars_full <- matrix(rep(p_vector, each=nrow(dadm)), ncol=length(p_vector), byrow=FALSE)
# pars_full[,ncol(pars_full)] <- pnorm(pars_full[,ncol(pars_full)])
#
# cv1_updated_ptr <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
#                                          pars_full = pars_full, return_kernel = TRUE, use_ptrs = TRUE)
#
# cv1_updated_noptr <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
#                                            pars_full = pars_full, return_kernel = TRUE, use_ptrs = FALSE)
#
# microbenchmark::microbenchmark(noptrs=EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
#                                                             pars_full = pars_full, return_kernel = TRUE, use_ptrs = FALSE),
#                                ptrs=EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
#                                                           pars_full = pars_full, return_kernel = TRUE, use_ptrs = TRUE)
# )




