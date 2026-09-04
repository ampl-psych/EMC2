RNGkind("L'Ecuyer-CMRG")
set.seed(321)

make_active_likelihood_design <- function(go, stop) {
  go_formula <- switch(
    go,
    exgaussian = list(mu ~ lM, sigma ~ 1, tau ~ 1),
    racing_diffusion = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1)
  )
  stop_formula <- switch(
    stop,
    exgaussian = list(muS ~ 1, sigmaS ~ 1, tauS ~ 1),
    lognormal = list(meanlogS ~ 1, sdlogS ~ 1),
    weibull = list(shapeS ~ 1, scaleS ~ 1),
    normal = list(muS ~ 1, sigmaS ~ 1)
  )

  design(
    model = function() stop_signal(go = go, stop = stop),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = c(go_formula, stop_formula, list(gf ~ 1, tf ~ 1))
  )
}

active_likelihood_values <- c(
  mu = log(.55),
  mu_lMTRUE = log(.85),
  sigma = log(.07),
  tau = log(.15),
  v = log(1.8),
  v_lMTRUE = log(.35),
  B = log(.8),
  A = log(.10),
  t0 = log(.20),
  s = log(1),
  muS = log(.24),
  sigmaS = log(.06),
  tauS = log(.06),
  meanlogS = log(.24) - .25^2 / 2,
  sdlogS = log(.25),
  shapeS = log(2.4),
  scaleS = log(.27),
  gf = qnorm(.05),
  tf = qnorm(.05)
)

make_active_likelihood_p <- function(design) {
  p <- sampled_pars(design, doMap = FALSE)
  missing <- setdiff(names(p), names(active_likelihood_values))
  if (length(missing)) {
    stop("Missing active likelihood value(s): ", paste(missing, collapse = ", "))
  }
  p[] <- active_likelihood_values[names(p)]
  p
}

expect_stop_signal_likelihood_active <- function(go, stop) {
  design_obj <- make_active_likelihood_design(go, stop)
  p0 <- make_active_likelihood_p(design_obj)
  dat <- make_data(
    p0,
    design_obj,
    n_trials = 40,
    functions = list(SSD = make_ssd())
  )
  emc <- make_emc(
    dat,
    design_obj,
    type = "single",
    n_chains = 1,
    compress = FALSE
  )

  dadm <- emc[[1]]$data[[1]]
  min_total_ll <- length(attr(dadm, "expand")) * log(1e-10)

  ll0 <- calc_ll_manager(t(p0), dadm = dadm, model = emc[[1]]$model)

  p1 <- p0
  if ("mu" %in% names(p1)) {
    p1["mu"] <- p1["mu"] + .10
  } else {
    p1["v"] <- p1["v"] + .10
  }
  ll1 <- calc_ll_manager(t(p1), dadm = dadm, model = emc[[1]]$model)

  expect_true(is.finite(ll0))
  expect_true(is.finite(ll1))
  expect_gt(unname(ll0), min_total_ll)
  expect_gt(unname(ll1), min_total_ll)
  expect_false(isTRUE(all.equal(unname(ll0), unname(ll1))))
}

test_that("default ex-Gaussian go stop-signal likelihoods are active", {
  for (stop in c("exgaussian", "lognormal", "weibull", "normal")) {
    expect_stop_signal_likelihood_active("exgaussian", stop)
  }
})

test_that("default racing-diffusion go stop-signal likelihoods are active", {
  for (stop in c("exgaussian", "lognormal", "weibull", "normal")) {
    expect_stop_signal_likelihood_active("racing_diffusion", stop)
  }
})
