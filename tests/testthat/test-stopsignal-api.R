test_that("stopsignal accepts the supported model combinations", {
  go_families <- c("exgaussian", "racing_diffusion")
  stop_families <- c("exgaussian", "lognormal", "weibull")

  for (go in go_families) {
    for (stop in stop_families) {
      model <- stopsignal(go = go, stop = stop)

      expect_type(model, "list")
      expect_identical(model$ss_info$go_family, go)
      expect_identical(model$ss_info$stop_family, stop)
      expect_true(all(c("tf", "gf") %in% names(model$p_types)))
      expect_true(all(c("tf", "gf") %in% names(model$transform$func)))
      expect_true(all(c("dfunG", "pfunG", "dfunS", "pfunS", "sfun", "rfun") %in% names(model)))
    }
  }
})

test_that("stopsignal rejects unsupported model families", {
  expect_error(stopsignal(go = "wald"), "'arg' should be one of")
  expect_error(stopsignal(stop = "wald"), "'arg' should be one of")
})

test_that("stopsignal exposes the expected parameter names", {
  expect_named(
    stopsignal(go = "exgaussian", stop = "exgaussian")$p_types,
    c("mu", "sigma", "tau", "muS", "sigmaS", "tauS", "tf", "gf", "exg_lb", "exgS_lb")
  )
  expect_named(
    stopsignal(go = "exgaussian", stop = "lognormal")$p_types,
    c("mu", "sigma", "tau", "meanlogS", "sdlogS", "tf", "gf", "exg_lb")
  )
  expect_named(
    stopsignal(go = "exgaussian", stop = "weibull")$p_types,
    c("mu", "sigma", "tau", "shapeS", "scaleS", "tf", "gf", "exg_lb")
  )
  expect_named(
    stopsignal(go = "racing_diffusion", stop = "exgaussian")$p_types,
    c("v", "B", "A", "t0", "s", "muS", "sigmaS", "tauS", "tf", "gf", "exgS_lb")
  )
  expect_named(
    stopsignal(go = "racing_diffusion", stop = "lognormal")$p_types,
    c("v", "B", "A", "t0", "s", "meanlogS", "sdlogS", "tf", "gf")
  )
  expect_named(
    stopsignal(go = "racing_diffusion", stop = "weibull")$p_types,
    c("v", "B", "A", "t0", "s", "shapeS", "scaleS", "tf", "gf")
  )
})

test_that("stopsignal routes only existing C++ stop-signal backends", {
  expect_identical(stopsignal(go = "exgaussian", stop = "exgaussian")$c_name, "SSEXG")
  expect_identical(stopsignal(go = "racing_diffusion", stop = "exgaussian")$c_name, "SSRDEX")

  expect_null(stopsignal(go = "exgaussian", stop = "lognormal")$c_name)
  expect_null(stopsignal(go = "exgaussian", stop = "weibull")$c_name)
  expect_null(stopsignal(go = "racing_diffusion", stop = "lognormal")$c_name)
  expect_null(stopsignal(go = "racing_diffusion", stop = "weibull")$c_name)
})

test_that("SSEXG and SSRDEX are compatibility wrappers for stopsignal", {
  expect_identical(SSEXG()$p_types, stopsignal(go = "exgaussian", stop = "exgaussian")$p_types)
  expect_identical(SSEXG()$c_name, stopsignal(go = "exgaussian", stop = "exgaussian")$c_name)

  expect_identical(SSRDEX()$p_types, stopsignal(go = "racing_diffusion", stop = "exgaussian")$p_types)
  expect_identical(SSRDEX()$c_name, stopsignal(go = "racing_diffusion", stop = "exgaussian")$c_name)
})
