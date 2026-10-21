RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR

n_trials <- 10

covariate1 <- seq(.1, 1, length.out = n_trials * 2)
covariate2 <- rep(c(.2, .5, .8, .4), length.out = n_trials * 2)
covariate2[1:5] <- 0

trend <- make_trend(make_base("m", "lin", make_kernel(c("covariate1", "covariate2"), "exp_incr")))

design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.w ~ 1, m.d_ei ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)
dat <- data.frame(subjects = factor(1), S = factor(rep(1:2, each = n_trials)),
                  R = factor(rep(c(1, 1, 2, 2, 1), length.out = n_trials * 2)),
                  rt = rep(c(.45, .6, .8, 1, 1.2), length.out = n_trials * 2),
                  covariate1 = covariate1, covariate2 = covariate2)

# Small trend effects around a plausible LNR, on the unmapped scale.
trend_parameters <- function(emc) {
  centre <- setNames(rep(.1, length(emc[[1]]$par_names)), emc[[1]]$par_names)
  centre[c("m", "m_lMd", "s", "t0")] <- c(-1, 1.5, log(.5), log(.2))
  likelihood_parameters(centre, sd = .15)
}

LNR2cov <- make_emc(dat, design_base, compress = F, n_chains = 1, type = "single")

test_that("trend", {
  parameters <- trend_parameters(LNR2cov)
  expect_snapshot(round(model_likelihood(LNR2cov, parameters), 6))
})

trend_2types <- make_trend(make_base("m", "lin", make_kernel(c("covariate1", "covariate2"), "exp_incr")),
  make_base("m_lMd", "lin", make_kernel("covariate1", "pow_decr")))

design_base_shared <- design(data = dat,
                             trend = trend_2types,
                             formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.w ~ 1, m.d_ei ~ 1, m_lMd.w ~ 1, m_lMd.d_pd ~ 1),
                             contrasts = list(lM = ADmat),
                             matchfun = matchfun,
                             model = LNR)

LNR2cov_shared <- make_emc(dat, design_base_shared, compress = FALSE, n_chains = 1, type = "single")

test_that("trend_shared", {
  parameters <- trend_parameters(LNR2cov_shared)
  expect_snapshot(round(model_likelihood(LNR2cov_shared, parameters), 6))
})

trend_premap <- make_trend(make_base("m", "lin", make_kernel("covariate1", "exp_incr")),
  make_base("m_lMd", "add", make_kernel("covariate2", "poly2")))

design_premap <- design(
  data = dat,
  trend = trend_premap,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m_lMd.d1 ~ S, m.w ~ 1, m.d_ei ~ 1, m_lMd.d2 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_premap <- make_emc(dat, design_premap, compress = FALSE, n_chains = 1, type = "single")
test_that("premap trend works", {
  parameters <- trend_parameters(LNR_premap)
  expect_snapshot(round(model_likelihood(LNR_premap, parameters), 6))
})

trend_pretrans <- make_trend(make_base("m", "lin", make_kernel("covariate1", "delta", at = NULL), phase = "pretransform"),
  make_base("s", "lin", make_kernel("covariate2", "exp_decr", at = NULL), phase = "pretransform"))

design_pretrans <- design(
  data = dat,
  trend = trend_pretrans,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, s.w ~ lR, m.w ~ 1, m.q0 ~ 1, m.alpha ~ 1, s.d_ed ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_pretrans <- make_emc(dat, design_pretrans, compress = FALSE, n_chains = 1, type = "single")

test_that("pretransform trend works", {
  parameters <- trend_parameters(LNR_pretrans)
  expect_snapshot(round(model_likelihood(LNR_pretrans, parameters), 6))
})

trend_posttrans <- make_trend(make_base("m", "lin", make_kernel("covariate1", "pow_decr"), phase = "posttransform"),
  make_base("s", "lin", make_kernel("covariate2", "pow_incr"), phase = "posttransform"))

design_posttrans <- design(
  data = dat,
  trend = trend_posttrans,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, s.w ~ lR, m.w ~ 1, m.d_pd ~ 1, s.d_pi ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_posttrans <- make_emc(dat, design_posttrans, compress = FALSE, n_chains = 1, type = "single")

test_that("posttransform trend works", {
  parameters <- trend_parameters(LNR_posttrans)
  expect_snapshot(round(model_likelihood(LNR_posttrans, parameters), 6))
})

trend_bases <- make_trend(make_base("m", "centered", make_kernel("covariate1", "exp_incr")),
  make_base("s", "lin", make_kernel("covariate2", "exp_decr")))

design_bases <- design(
  data = dat,
  trend = trend_bases,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.w ~ 1, m.d_ei ~ 1, s.w ~ 1, s.d_ed ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_bases <- make_emc(dat, design_bases, compress = FALSE, n_chains = 1, type = "single")

test_that("different trend base functions work", {
  parameters <- trend_parameters(LNR_bases)
  expect_snapshot(round(model_likelihood(LNR_bases, parameters), 6))
})

trend_poly <- make_trend(make_base("m", "add", make_kernel("covariate1", "poly3")),
  make_base("s", "add", make_kernel("covariate2", "poly4")))

design_poly <- design(
  data = dat,
  trend = trend_poly,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.d1 ~ 1, m.d2 ~ 1, m.d3 ~ 1, s.d1 ~ 1, s.d2 ~ 1, s.d3 ~ 1, s.d4 ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_poly <- make_emc(dat, design_poly, compress = FALSE, n_chains = 1, type = "single")

test_that("polynomial trends work", {
  parameters <- trend_parameters(LNR_poly)
  expect_snapshot(round(model_likelihood(LNR_poly, parameters), 6))
})

## New tests: phase-specific trends and par_input trends

# Three trends on different phases for LNR
trend_phases <- make_trend(make_base("m", "lin", make_kernel("covariate1", "lin_incr")),
  make_base("s", "lin", make_kernel("covariate1", "exp_decr"), phase = "pretransform"),
  make_base("t0", "lin", make_kernel("covariate2", "pow_incr"), phase = "posttransform"))

design_phases <- design(
  data = dat,
  trend = trend_phases,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.w ~ 1, s.w ~ 1, s.d_ed ~ 1, t0.w ~ 1, t0.d_pi ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_phases <- make_emc(dat, design_phases, compress = FALSE, n_chains = 1, type = "single")

test_that("phase-specific trends work (premap, pretransform, posttransform)", {
  parameters <- trend_parameters(LNR_phases)
  expect_snapshot(round(model_likelihood(LNR_phases, parameters), 6))
})

# Trend where input is another parameter: use t0 as input to a trend on m
trend_par_input <- make_trend(make_base("m", "lin", make_kernel(NULL, "lin_incr", par_input = "t0"), phase = "pretransform"))

design_par_input <- design(
  data = dat,
  trend = trend_par_input,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.w ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_par_input <- make_emc(dat, design_par_input, compress = FALSE, n_chains = 1, type = "single")

test_that("par_input trend uses t0 as input to m trend", {
  parameters <- trend_parameters(LNR_par_input)
  expect_snapshot(round(model_likelihood(LNR_par_input, parameters), 6))
})

trend_shared_premap <- make_trend(make_base("m", "add", make_kernel("covariate1", "poly3")),
  make_base("s", "add", make_kernel("covariate2", "poly4")),
  shared = list(shrd = list("m.d1", "s.d1")))

design_shared_premap <- design(
  data = dat,
  trend = trend_shared_premap,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.d2 ~ 1, m.d3 ~ 1, s.d2 ~ 1, s.d3 ~ 1, s.d4 ~ 1, shrd ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)

LNR_shared_premap <- make_emc(dat, design_shared_premap, compress = FALSE, n_chains = 1, type = "single")

test_that("share works premap trend", {
  parameters <- trend_parameters(LNR_shared_premap)
  expect_snapshot(round(model_likelihood(LNR_shared_premap, parameters), 6))
})


trend_shared_posttransform <- make_trend(make_base("m", "add", make_kernel("covariate1", "poly3"), phase = "posttransform"),
  make_base("s", "add", make_kernel("covariate2", "poly4"), phase = "posttransform"),
  shared = list(shrd = list("m.d1", "s.d1")))

design_shared_posttransform <- design(
  data = dat,
  trend = trend_shared_posttransform,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.d2 ~ 1, m.d3 ~ 1, s.d2 ~ 1, s.d3 ~ 1, s.d4 ~ 1, shrd ~ 1),
  contrasts = list(lM = ADmat),
  matchfun = matchfun,
  model = LNR
)
LNR_shared_posttransform <- make_emc(dat, design_shared_posttransform, compress = FALSE, n_chains = 1, type = "single")
test_that("share works posttransform trend", {
  parameters <- trend_parameters(LNR_shared_posttransform)
  expect_snapshot(round(model_likelihood(LNR_shared_posttransform, parameters), 6))
})


# Simulated RTs feed the next trial's delta update.
trend_cond <- make_trend(make_base("m", "lin", make_kernel("rt", "delta"), phase = "pretransform"))
design_cond <- design(
  factors = list(subjects = 1, S = 1:2), Rlevels = 1:2,
  covariates = "rt", trend = trend_cond,
  formula = list(m ~ 1, s ~ 1, t0 ~ 1, m.w ~ 1, m.q0 ~ 1, m.alpha ~ 1), model = LNR
)
p_vec <- c(m = -.5, s = log(.3), t0 = log(.2),
           m.w = .5, m.q0 = 0, m.alpha = qnorm(.2))

test_that("sequential simulation uses RT feedback", {
  withr::local_seed(123, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion")
  simulated <- make_data(
    p_vec, design_cond, n_trials = 10,
    covariates = data.frame(rt = rep(NA_real_, 20)),
    conditional_on_data = FALSE, return_trialwise_parameters = TRUE
  )
  q <- numeric(nrow(simulated))
  for (i in 2:length(q)) q[i] <- q[i - 1] + .2 * (simulated$rt[i - 1] - q[i - 1])
  expect_equal(attr(simulated, "trialwise_parameters")$m,
               rep(-.5 + .5 * q, each = 2), tolerance = 1e-8)
})

# Trend uses multiple same parameters and multiple inputs
trend_mult <- make_trend(make_base("m", "lin", make_kernel("covariate1", "exp_incr"), phase = "pretransform"),
  make_base("m", "lin", make_kernel(c("covariate1", "covariate2"), "delta", par_input = "t0"), phase = "pretransform"))

design_mult <- design(
  factors = list(subjects = 1, S = 1:2),
  Rlevels = 1:2,
  covariates = c("covariate1", "covariate2"),
  matchfun = matchfun,
  trend = trend_mult,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.w_exp_incr ~ 1, m.d_ei ~ 1, m.w_delta ~ 1, m.q0 ~ 1, m.alpha ~ 1),
  contrasts = list(lM = ADmat),
  model = LNR
)

LNR_multi <- make_emc(dat, design_mult, compress = FALSE, n_chains = 1, type = "single")

test_that("trend_multiple", {
  parameters <- trend_parameters(LNR_multi)
  expect_snapshot(round(model_likelihood(LNR_multi, parameters), 6))
})




# covariate maps
kernel <- make_kernel(c("covariate1", "covariate2"), "delta")
trend <- make_trend(
  make_base("m", "lin", kernel, coding = list(map1 = function(dadm, covs) {
    d <- cbind(dadm$covariate1, dadm$covariate2)
    colnames(d) <- covs
    d
  })),
  make_base("m", "lin", kernel, coding = list(map2 = function(dadm, covs) {
    d <- cbind(dadm$covariate2, 1 - dadm$covariate1)
    colnames(d) <- covs
    d
  })))
design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1', 'covariate2'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1, m.w_map1 ~ 1, m.w_map2 ~ 1, m.q0 ~ 1, m.alpha ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)

LNR_covmap <- make_emc(dat, design_base, compress = FALSE, n_chains = 1, type = "single")

test_that("trend_covmap", {
  parameters <- trend_parameters(LNR_covmap)
  expect_snapshot(round(model_likelihood(LNR_covmap, parameters), 6))
})
