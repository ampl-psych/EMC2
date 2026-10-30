ss_trend_design <- function(covariate = "SSD", ...) {
  design(model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
         Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
         formula = list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1,
                        tauS ~ 1, gf ~ 1, tf ~ 1), report_p_vector = FALSE,
         trend = make_trend(make_base("muS", "lin", make_kernel(covariate, "slin_incr"),
                                      phase = "posttransform")),
         transform = list(func = c(muS.w = "exp")), ...)
}
ss_trend_pars <- c(mu = log(.5), mu_lMTRUE = .2, sigma = log(.05), tau = log(.1),
                   muS = log(.12), sigmaS = log(.03), tauS = log(.05),
                   gf = qnorm(.05), tf = qnorm(.1), muS.k_sat = log(5), muS.w = log(.2))

test_that("stop trends match manual likelihoods with caching, compression and threads", {
  # Equal SSDs but different stop means must not share a cached integral.
  dat <- data.frame(subjects = factor(1), S = factor("left", levels = c("left", "right")),
                    R = factor(c("left", NA, NA, "right", NA, "left"), levels = c("left", "right")),
                    rt = c(.5, NA, NA, .65, NA, .55), SSD = c(Inf, .2, .2, .2, .4, Inf),
                    x = c(0, .1, .4, .4, .1, 0))
  dat <- make_missing(rbind(dat, dat), UC = .8, rt_resolution = NULL)
  des <- ss_trend_design("x", covariates = "x")
  plain <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE)
  compressed <- make_emc(dat, des, type = "single", n_chains = 1, compress = TRUE)
  pm <- rbind(ss_trend_pars, ss_trend_pars)
  pm[2, "muS.w"] <- log(.35)
  dadm <- plain[[1]]$data[[1]]
  model <- des$model()
  model$trend <- NULL
  expected <- apply(pm, 1, function(p) {
    pars <- get_pars_matrix_oo(p, dadm, model)
    pars[, "muS"] <- exp(p[["muS"]]) + exp(p[["muS.w"]]) * pmin(1, exp(p[["muS.k_sat"]]) * dadm$x)
    model$log_likelihood(pars, dadm, model)
  })
  expect_equal(as.numeric(model_likelihood(plain, pm)), unname(expected), tolerance = 1e-8)
  expect_equal(as.numeric(model_likelihood(compressed, pm)), unname(expected), tolerance = 1e-8)
  withr::local_options(emc.ll_backend = "multithreaded", emc.n_threads = 2)
  expect_equal(as.numeric(suppressWarnings(model_likelihood(compressed, pm))),
               unname(expected), tolerance = 1e-8)
})

test_that("SSD trends follow separate staircases and deadline outcomes", {
  withr::local_seed(21, .rng_kind = "L'Ecuyer-CMRG")
  des <- ss_trend_design(TC = list(UC = .55), functions = list(SSD =
    make_ssd(SSD0 = .2, stairstep = .05, p_stop = 1, factors = "S", UC = .55)))
  dat <- suppressMessages(make_data(ss_trend_pars, des, n_trials = 12,
                                    return_trialwise_parameters = TRUE))
  for (st in split(dat, dat$S)) {
    expect_equal(st$SSD[1], .2)
    next_ssd <- pmax(0, head(st$SSD, -1) + ifelse(is.na(head(st$R, -1)), .05, -.05))
    expect_equal(tail(st$SSD, -1), next_ssd, tolerance = 1e-12)
  }
  tw <- attr(dat, "trialwise_parameters")
  tw <- tw[!duplicated(tw[, c("subject", "trial")]), , drop = FALSE]
  expect_equal(unname(tw[, "muS"]), .12 + .2 * pmin(1, 5 * dat$SSD), tolerance = 1e-12)
})

test_that("stop-trend prediction reuses observed SSDs or regenerates the staircase", {
  withr::local_seed(24, .rng_kind = "L'Ecuyer-CMRG")
  des <- ss_trend_design(functions = list(SSD = make_ssd(p_stop = 1)))
  dat <- data.frame(subjects = factor(1), S = factor(rep(c("left", "right"), 4)),
                    R = factor("left", levels = c("left", "right")), rt = .6, SSD = .1)
  emc <- make_emc(dat, des, type = "single", n_chains = 1, rt_resolution = 1e-12)
  p <- ss_trend_pars[emc[[1]]$par_names]
  # Fixed posterior draws isolate prediction from sampler convergence.
  emc[[1]]$samples <- list(alpha = array(c(p, p + .01), c(length(p), 1, 2),
                                        dimnames = list(names(p), "1", NULL)),
                           stage = rep("sample", 2), idx = 2L)
  conditional <- predict(emc, n_post = 1, n_cores = 1, conditional_on_data = TRUE)
  unconditional <- suppressMessages(predict(emc, n_post = 1, n_cores = 1))
  expect_equal(conditional$SSD, dat$SSD)
  expect_equal(nrow(unconditional), nrow(dat))
  expect_equal(unconditional$SSD[1], .25)
  expect_false(isTRUE(all.equal(unconditional$SSD, dat$SSD)))
})
