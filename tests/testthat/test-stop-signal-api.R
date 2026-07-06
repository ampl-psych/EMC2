test_that("stop_signal accepts the supported model combinations", {
  go_families <- c("exgaussian", "racing_diffusion")
  stop_families <- c("exgaussian", "lognormal", "weibull")

  for (go in go_families) {
    for (stop in stop_families) {
      model <- stop_signal(go = go, stop = stop)

      expect_type(model, "list")
      expect_identical(model$ss_info$go_family, go)
      expect_identical(model$ss_info$stop_family, stop)
      expect_true(all(c("tf", "gf") %in% names(model$p_types)))
      expect_true(all(c("tf", "gf") %in% names(model$transform$func)))
      expect_true(all(c("dfunG", "pfunG", "dfunS", "pfunS", "sfun", "rfun") %in% names(model)))
    }
  }
})

test_that("stop_signal rejects unsupported model families", {
  expect_error(stop_signal(go = "wald"), "'arg' should be one of")
  expect_error(stop_signal(stop = "wald"), "'arg' should be one of")
})

test_that("stop_signal exposes the expected parameter names", {
  expect_named(
    stop_signal(go = "exgaussian", stop = "exgaussian")$p_types,
    c("mu", "sigma", "tau", "muS", "sigmaS", "tauS", "tf", "gf", "exg_lb", "exgS_lb")
  )
  expect_named(
    stop_signal(go = "exgaussian", stop = "lognormal")$p_types,
    c("mu", "sigma", "tau", "meanlogS", "sdlogS", "tf", "gf", "exg_lb")
  )
  expect_named(
    stop_signal(go = "exgaussian", stop = "weibull")$p_types,
    c("mu", "sigma", "tau", "shapeS", "scaleS", "tf", "gf", "exg_lb")
  )
  expect_named(
    stop_signal(go = "racing_diffusion", stop = "exgaussian")$p_types,
    c("v", "B", "A", "t0", "s", "muS", "sigmaS", "tauS", "tf", "gf", "exgS_lb")
  )
  expect_named(
    stop_signal(go = "racing_diffusion", stop = "lognormal")$p_types,
    c("v", "B", "A", "t0", "s", "meanlogS", "sdlogS", "tf", "gf")
  )
  expect_named(
    stop_signal(go = "racing_diffusion", stop = "weibull")$p_types,
    c("v", "B", "A", "t0", "s", "shapeS", "scaleS", "tf", "gf")
  )
})

test_that("stop_signal routes only existing C++ stop-signal backends", {
  expect_identical(stop_signal(go = "exgaussian", stop = "exgaussian")$c_name, "SSEXG")
  expect_identical(stop_signal(go = "exgaussian", stop = "lognormal")$c_name, "SSLNORM")
  expect_identical(stop_signal(go = "exgaussian", stop = "weibull")$c_name, "SSWEIBULL")
  expect_identical(stop_signal(go = "racing_diffusion", stop = "exgaussian")$c_name, "SSRDEX")
  expect_identical(stop_signal(go = "racing_diffusion", stop = "lognormal")$c_name, "SSRDLNORM")
  expect_identical(stop_signal(go = "racing_diffusion", stop = "weibull")$c_name, "SSRDWEIBULL")
})

test_that("SSEXG and SSRDEX are compatibility wrappers for stop_signal", {
  expect_identical(SSEXG()$p_types, stop_signal(go = "exgaussian", stop = "exgaussian")$p_types)
  expect_identical(SSEXG()$c_name, stop_signal(go = "exgaussian", stop = "exgaussian")$c_name)

  expect_identical(SSRDEX()$p_types, stop_signal(go = "racing_diffusion", stop = "exgaussian")$p_types)
  expect_identical(SSRDEX()$c_name, stop_signal(go = "racing_diffusion", stop = "exgaussian")$c_name)
})

test_that("design accepts a stop_signal model object directly", {
  direct <- design(
    model = stop_signal(),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(
      mu ~ lM,
      sigma ~ 1,
      tau ~ 1,
      muS ~ 1,
      sigmaS ~ 1,
      tauS ~ 1,
      gf ~ 1,
      tf ~ 1
    ),
    report_p_vector = FALSE
  )

  wrapped <- design(
    model = function() stop_signal(),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(
      mu ~ lM,
      sigma ~ 1,
      tau ~ 1,
      muS ~ 1,
      sigmaS ~ 1,
      tauS ~ 1,
      gf ~ 1,
      tf ~ 1
    ),
    report_p_vector = FALSE
  )

  expect_true(is.function(direct$model))
  expect_identical(direct$model()$p_types, wrapped$model()$p_types)
  expect_identical(direct$model()$c_name, wrapped$model()$c_name)
})
