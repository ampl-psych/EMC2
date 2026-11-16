# Adapted from the review patch. Fixed draws replace simulation and MCMC.
ssrt_test_fit <- function(model, muS_formula = muS ~ 1, lower = -Inf,
                          type = "standard") {
  withr::local_seed(21)
  des <- suppressMessages(design(
    model = model, factors = list(subjects = 1:3, S = c("left", "right")),
    Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
    formula = list(muS_formula, sigmaS ~ 1, tauS ~ 1),
    constants = c(exgS_lb = lower), report_p_vector = FALSE
  ))
  p <- sampled_pars(des, doMap = FALSE)
  values <- c(muS = log(.2), sigmaS = log(.03), tauS = log(.05),
              muS_Sleft = log(.2),
              muS_Sright = if ("muS" %in% names(p)) log(1.5) else log(.3))
  p[] <- values[names(p)]
  stopifnot(!anyNA(p))
  dat <- ss_test_data()
  n <- nrow(dat)
  dat <- dat[rep(seq_len(n), 3), ]
  dat$subjects <- factor(rep(1:3, each = n))
  emc <- suppressMessages(make_emc(dat, des, type = type, n_chains = 1,
                                  compress = FALSE, verbose = FALSE))
  alpha <- array(rep(p, 12), c(length(p), 3, 4),
                 dimnames = list(names(p), as.character(1:3), NULL))
  mean_rows <- if ("muS" %in% names(p)) "muS" else c("muS_Sleft", "muS_Sright")
  for (s in 1:3) for (i in 1:4) {
    alpha[mean_rows, s, i] <- p[mean_rows] + log(c(1, 1.5, 2)[s] * c(.8, 1, 1.2, 1.4)[i])
  }
  # Zero population variance gives exact expectations for the population path.
  emc[[1]]$samples <- list(
    alpha = alpha, theta_mu = alpha[, 1, ],
    theta_var = array(0, c(length(p), length(p), 4),
                      dimnames = list(names(p), names(p), NULL)),
    stage = rep("sample", 4), idx = 4
  )
  emc
}

expect_ssrt_moments <- function(summary, mean_draws, sd_draws) {
  for (parameter in c("ssrt_mean", "ssrt_sd")) {
    x <- if (parameter == "ssrt_mean") mean_draws else sd_draws
    expected <- c(mean(x), sd(x), quantile(x, c(.025, .5, .975)))
    actual <- summary[summary$parameter == parameter, c("mean", "sd", "2.5%", "50%", "97.5%")]
    expect_equal(nrow(actual), 1L)
    expect_equal(unname(unlist(actual)), unname(expected), tolerance = 1e-10)
  }
}

for (model_name in c("SSEXG", "SSRDEX")) {
  model <- get(model_name)
  test_that(paste(model_name, "individual and population summaries use the fixed draws"), {
    emc <- ssrt_test_fit(model)
    ind <- ssrt_summary(emc, level = "individual")
    expect_equal(nrow(ind), 6)
    expect_setequal(unique(ind$parameter), c("ssrt_mean", "ssrt_sd"))
    for (s in 1:3) {
      expect_ssrt_moments(ind[ind$subject == s, ],
                          .2 * c(1, 1.5, 2)[s] * c(.8, 1, 1.2, 1.4) + .05,
                          rep(sqrt(.03^2 + .05^2), 4))
    }
    withr::local_seed(12)
    rng <- .Random.seed
    expect_silent(pop <- ssrt_summary(emc, n_population = 8, seed = 1))
    expect_identical(.Random.seed, rng)
    expect_equal(nrow(pop), 2)
    expect_ssrt_moments(pop, .2 * c(.8, 1, 1.2, 1.4) + .05,
                        rep(sqrt(.03^2 + .05^2), 4))
    details <- attr(pop, "ssrt_details")
    expect_equal(details$population_evaluations, rep(32, 4))
    expect_true(all(details$invalid_population_evaluations == 0))
    expect_false(details$tail_dominated[details$parameter == "ssrt_mean"])
    expect_error(ssrt_summary(emc, level = "condition"), "not modeled")

    emc[[1]]$samples$theta_var["muS", "muS", ] <- .04
    varying <- ssrt_summary(emc, n_population = 8, seed = 1, details = TRUE)
    expect_gt(varying$mean[varying$parameter == "between_subject_sd"], 0)
    expect_identical(varying, ssrt_summary(emc, n_population = 8, seed = 1, details = TRUE))

    single <- ssrt_test_fit(model, type = "single")
    expect_identical(ssrt_summary(single), ssrt_summary(single, level = "individual"))
    expect_error(ssrt_summary(single, level = "population"), 'type = "single"')
  })

  test_that(paste(model_name, "stop effects are mapped to each subject and condition"), {
    for (formula in list(muS ~ S, muS ~ 0 + S)) {
      emc <- ssrt_test_fit(model, formula)
      ind <- ssrt_summary(emc, level = "condition", condition_level = "individual")
      expect_equal(nrow(ind), 12)
      for (s in 1:3) for (cell in c("left", "right")) {
        mu <- if (cell == "left") .2 else .3
        expect_ssrt_moments(ind[ind$subject == s & ind$S == cell, ],
                            mu * c(1, 1.5, 2)[s] * c(.8, 1, 1.2, 1.4) + .05,
                            rep(sqrt(.03^2 + .05^2), 4))
      }
      pop <- ssrt_summary(emc, level = "condition", n_population = 8, seed = 1)
      expect_setequal(as.character(unique(pop$S)), c("left", "right"))
      for (cell in c("left", "right")) {
        mu <- if (cell == "left") .2 else .3
        expect_ssrt_moments(pop[pop$S == cell, ], mu * c(.8, 1, 1.2, 1.4) + .05,
                            rep(sqrt(.03^2 + .05^2), 4))
      }
      expect_error(ssrt_summary(emc, level = "individual"), "condition effects")
      dat <- get_data(emc)
      selected <- ssrt_summary(emc, level = "condition", condition_level = "individual",
                               subject = "2", data = dat[nrow(dat):1, ])
      expected <- ind[ind$subject == "2", ]
      key <- c("S", "parameter")
      selected <- selected[do.call(order, selected[key]), ]
      expected <- expected[do.call(order, expected[key]), ]
      rownames(selected) <- rownames(expected) <- NULL
      attr(selected, "ssrt_details") <- attr(expected, "ssrt_details") <- NULL
      expect_identical(selected, expected)
    }
  })
}

test_that("SSEXG and SSRDEX lower-truncated moments agree with independent integration", {
  # Direct ex-Gaussian density; do not use the moment helper under test.
  reference <- vapply(.2 * c(.8, 1, 1.2, 1.4), function(mu) {
    density <- function(t) exp((mu - t) / .05 + .5 * (.03 / .05)^2 +
      pnorm((t - mu) / .03 - .03 / .05, log.p = TRUE)) / .05
    mass <- integrate(density, .25, Inf, rel.tol = 1e-10)$value
    mean <- integrate(function(t) t * density(t), .25, Inf, rel.tol = 1e-10)$value / mass
    second <- integrate(function(t) t^2 * density(t), .25, Inf, rel.tol = 1e-10)$value / mass
    c(mean = mean, sd = sqrt(second - mean^2))
  }, numeric(2))
  for (model in list(SSEXG, SSRDEX)) {
    ind <- ssrt_summary(ssrt_test_fit(model, lower = .25), level = "individual")
    expect_ssrt_moments(ind[ind$subject == "1", ], reference["mean", ], reference["sd", ])
  }
})

test_that("population diagnostics distinguish omitted, infinite and partial draws", {
  dat <- data.frame(subjects = factor(rep(1:3, 2)), S = rep(c("A", "B"), each = 3))
  pars <- array(0, c(6, 5, 2), dimnames = list(NULL, NULL, c("meanlogS", "sdlogS")))
  pars[, , "meanlogS"] <- log(.25)
  pars[, , "sdlogS"] <- .2
  pars[1, 2, "meanlogS"] <- NA_real_
  pars[1:3, 3, "meanlogS"] <- NaN
  pars[1, 4, "meanlogS"] <- 1000
  pars[1, 5, "meanlogS"] <- 400
  out <- ssrt_summary_condition_population(dat, pars, list(stop_family = "lognormal"),
                                            c(.025, .5, .975), "S", 3)
  affected <- out[out$S == "A", ]
  expect_equal(affected$invalid_population_evaluations, c(5, 6, 5, 6))
  expect_equal(affected$omitted_draws, c(1, 3, 3, 1))
  expect_equal(affected$infinite_draws, c(1, 1, 1, 2))
  expect_equal(affected$partial_draws, c(1, 0, 0, 1))
  expect_true(all(out$population_evaluations == 15))
  expect_true(all(out$posterior_draws == 5))
  expect_true(all(out$invalid_population_evaluations[out$S == "B"] == 0))
  expect_warning(compact <- ssrt_format_output(out, FALSE), "Non-finite population SSRT")
  expect_identical(attr(compact, "ssrt_details"), out)
  expect_false("omitted_draws" %in% names(compact))
})

test_that("upper-tail messages use finite draws and identify affected conditions", {
  dat <- data.frame(subjects = factor(rep(1:2, 2)), S = rep(c("A", "B"), each = 2))
  pars <- array(0, c(4, 100, 2), dimnames = list(NULL, NULL, c("meanlogS", "sdlogS")))
  pars[, , "meanlogS"] <- log(.25)
  pars[, , "sdlogS"] <- .1
  summarize <- function(p) ssrt_summary_condition_population(
    dat, p, list(stop_family = "lognormal"), c(.1, .5, .9), "S", 2)
  expect_silent(ssrt_format_output(summarize(pars), TRUE))
  pars[1:2, 1, "meanlogS"] <- Inf
  pars[1:2, 2, "meanlogS"] <- NA_real_
  out <- summarize(pars)
  expect_false(any(out$tail_dominated, na.rm = TRUE))
  expect_silent(suppressWarnings(ssrt_format_output(out, TRUE)))

  pars[1:2, 100, "meanlogS"] <- log(1000)
  out <- summarize(pars)
  expect_identical(out$tail_dominated[out$parameter == "ssrt_mean"], c(TRUE, FALSE))
  expect_true(all(is.na(out$tail_dominated[out$parameter != "ssrt_mean"])))
  expect_message(expect_warning(ssrt_format_output(out, TRUE), "Non-finite"),
                 "condition\\(s\\) A.*finite draws")

  pars[1:2, , "meanlogS"] <- NA_real_
  out <- summarize(pars)
  expect_true(is.na(out$tail_dominated[out$S == "A" & out$parameter == "ssrt_mean"]))
  expect_silent(suppressWarnings(ssrt_format_output(out, TRUE)))
})
