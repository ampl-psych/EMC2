# Integrate over stop probability, not the coordinates used by the C++ kernel.
# R's quantiles also resolve very narrow stops without sharing the C++ integral.
ss_reference_stop_probability <- function(family, p, survivor, upper = Inf) {
  if (family == "normal") {
    log_tail <- pnorm(p["normS_lb"], p["muS"], p["sigmaS"],
                      lower.tail = FALSE, log.p = TRUE)
    quantile <- function(u) qnorm(log_tail + log1p(-u), p["muS"], p["sigmaS"],
                                  lower.tail = FALSE, log.p = TRUE)
    mass <- if (upper <= p["normS_lb"]) 0 else
      -expm1(pnorm(upper, p["muS"], p["sigmaS"],
                    lower.tail = FALSE, log.p = TRUE) - log_tail)
  } else if (family == "lognormal") {
    quantile <- function(u) qlnorm(u, p["meanlogS"], p["sdlogS"])
    mass <- plnorm(upper, p["meanlogS"], p["sdlogS"])
  } else {
    quantile <- function(u) qweibull(u, p["shapeS"], p["scaleS"])
    mass <- pweibull(upper, p["shapeS"], p["scaleS"])
  }
  if (mass == 0) return(0)
  integrate(function(u) survivor(quantile(u) + p["SSD"]), 0, mass,
            abs.tol = 1e-10, rel.tol = 1e-9, subdivisions = 500L)$value
}

test_that("shared stop integrals agree with independent R quadrature", {
  cases <- list(
    normal = list(c(muS = .25, sigmaS = .001, normS_lb = 0),
                  c(muS = .25, sigmaS = .1, normS_lb = 0),
                  c(muS = .2, sigmaS = .001, normS_lb = .25)),
    lognormal = list(c(meanlogS = log(.25), sdlogS = .001),
                     c(meanlogS = log(.25), sdlogS = .3),
                     c(meanlogS = log(.25), sdlogS = 3)),
    weibull = list(c(shapeS = .4, scaleS = .3),
                   c(shapeS = 2.4, scaleS = .3),
                   c(shapeS = 50, scaleS = .3))
  )
  for (go in c("exgaussian", "racing_diffusion")) {
    if (go == "exgaussian") {
      # The enabled approximation is a normal truncated at zero.
      g <- c(mu = .5, sigma = .1, tau = .04, exg_lb = 0,
             exg_normal_ratio = .5, exg_shifted_exp_sigma = 0)
      survivor <- function(t) (pnorm(pmax(t, 0), .5, .1, lower.tail = FALSE) /
                                pnorm(0, .5, .1, lower.tail = FALSE))^2
    } else {
      # A=0 gives two ordinary inverse-Gaussian go runners, with a known CDF.
      g <- c(v = 3, B = 1.2, A = 0, t0 = .15, s = 1)
      survivor <- function(t) {
        dt <- pmax(t - .15, .Machine$double.eps)
        (pnorm((1.2 - 3 * dt) / sqrt(dt)) -
           exp(2 * 3 * 1.2) * pnorm(-(1.2 + 3 * dt) / sqrt(dt)))^2
      }
    }
    for (family in names(cases)) {
      type <- stop_signal(go = go, stop = family)$c_name
      for (s in cases[[family]]) {
        p <- c(g, s, SSD = .2)
        pars <- rbind(p, p)
        for (upper in c(.28, Inf)) {
          expected <- ss_reference_stop_probability(family, p, survivor, upper)
          actual <- ss_stop_probability(pars, 2L, type, upper)
          expect_true(is.finite(actual))
          expect_true(abs(actual - expected) < 1e-7,
                      info = paste(go, family, paste(s, collapse = ","), upper))
        }
      }
    }
  }
})

test_that("non-exGaussian stops retain both go approximation settings in fitting", {
  withr::local_options(emc.ll_backend = "multiprocess")
  formulas <- list(normal = list(muS ~ 1, sigmaS ~ 1),
                   lognormal = list(meanlogS ~ 1, sdlogS ~ 1),
                   weibull = list(shapeS ~ 1, scaleS ~ 1))
  values <- c(mu = log(.5), sigma = log(.1), tau = log(.04),
              muS = log(.24), sigmaS = log(.06), meanlogS = log(.25),
              sdlogS = log(.3), shapeS = log(2.4), scaleS = log(.3),
              gf = qnorm(.05), tf = qnorm(.1))
  dat <- ss_test_data()[c(1, 3, 5, 6, 7), ]
  dat$rt[1:2] <- .56
  settings <- list(default = NULL, normal = list(normal_ratio = .5),
                   exponential = list(shifted_exponential_sigma = .15))
  for (family in names(formulas)) {
    baseline <- NULL
    for (setting in names(settings)) {
      model <- suppressWarnings(stop_signal(go = "exgaussian", stop = family,
                                             exgaussian_approximation = settings[[setting]]))
      des <- design(model = model, factors = list(subjects = 1, S = c("left", "right")),
                    Rlevels = c("left", "right"),
                    formula = c(list(mu ~ 1, sigma ~ 1, tau ~ 1), formulas[[family]],
                                list(gf ~ 1, tf ~ 1)), constants = c(exg_lb = 0),
                    report_p_vector = FALSE)
      p <- make_stop_cpp_p(des, values)
      emc <- make_emc(dat, des, type = "single", n_chains = 1,
                      compress = FALSE, verbose = FALSE)
      proposals <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))
      actual <- as.numeric(calc_ll_manager(proposals, emc[[1]]$data[[1]],
                                            emc[[1]]$model, return_trialwise = TRUE))
      reference <- emc[[1]]$model()
      reference$c_name <- NULL
      expected <- as.numeric(calc_ll_manager(proposals, emc[[1]]$data[[1]],
                                              function() reference))
      expect_true(all(is.finite(actual)))
      expect_equal(sum(actual), expected, tolerance = 1e-7)
      if (setting == "default") baseline <- actual else {
        # Both observed responses and stop nonresponses must use the setting.
        expect_gt(abs(actual[1] - baseline[1]), 1e-4)
        expect_gt(max(abs(actual[3:4] - baseline[3:4])), 1e-4)
        survivor <- if (setting == "normal") {
          function(t) (pnorm(pmax(t, 0), .5, .1, lower.tail = FALSE) /
                         pnorm(0, .5, .1, lower.tail = FALSE))^2
        } else function(t) exp(-2 * pmax(t - .5, 0) / .04)
        natural <- c(muS = .24, sigmaS = .06, normS_lb = 0, meanlogS = log(.25),
                     sdlogS = .3, shapeS = 2.4, scaleS = .3)
        for (i in 3:4) {
          ps <- ss_reference_stop_probability(family, c(natural, SSD = dat$SSD[i]), survivor)
          expect_equal(actual[i], log(.05 + .95 * .9 * ps), tolerance = 1e-7)
        }
      }
    }
  }
})

test_that("stop-triggered responses require stop and go to remain unfinished", {
  withr::local_options(emc.ll_backend = "multiprocess")
  model <- suppressWarnings(stop_signal(stop = "normal",
                                         exgaussian_approximation = list(normal_ratio = .5)))
  des <- design(model = model, factors = list(subjects = 1, S = "stimulus"),
                Rlevels = c("go", "change"),
                functions = list(lI = function(d)
                  factor(ifelse(d$lR == "change", 1, 2), levels = 1:2)),
                formula = list(mu ~ 1, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1,
                               gf ~ 1, tf ~ 1), constants = c(exg_lb = 0),
                report_p_vector = FALSE)
  p <- make_stop_cpp_p(des, c(mu = log(.5), sigma = log(.1), tau = log(.04),
                             muS = log(.24), sigmaS = log(.06),
                             gf = qnorm(.05), tf = qnorm(.1)))
  dat <- data.frame(subjects = factor(1), S = factor("stimulus"), trials = 1:3,
                    R = factor("change", levels = c("go", "change")),
                    rt = c(.35, .55, .8), SSD = .2)
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE, verbose = FALSE)
  actual <- as.numeric(calc_ll_manager(matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
                                        emc[[1]]$data[[1]], emc[[1]]$model,
                                        return_trialwise = TRUE))
  survivor <- function(t) pnorm(pmax(t, 0), .5, .1, lower.tail = FALSE) /
    pnorm(0, .5, .1, lower.tail = FALSE)
  ps <- vapply(dat$rt - .2, function(upper) ss_reference_stop_probability(
    "normal", c(muS = .24, sigmaS = .06, normS_lb = 0, SSD = .2), survivor, upper), 0.0)
  stop_survivor <- pnorm(dat$rt - .2, .24, .06, lower.tail = FALSE) /
    pnorm(0, .24, .06, lower.tail = FALSE)
  st_density <- dnorm(dat$rt - .2, .5, .1) / pnorm(0, .5, .1, lower.tail = FALSE)
  expected <- log(.9 * st_density * (.05 + .95 * (ps + stop_survivor * survivor(dat$rt))))
  expect_equal(actual, expected, tolerance = 1e-7)
})
