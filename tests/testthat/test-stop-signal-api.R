test_that("stop_signal accepts the supported model combinations", {
  go_families <- c("exgaussian", "racing_diffusion")
  stop_families <- c("exgaussian", "lognormal", "weibull", "normal")

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
    c(
      "mu", "sigma", "tau", "muS", "sigmaS", "tauS", "tf", "gf",
      "exg_lb", "exgS_lb", "exg_normal_ratio", "exg_shifted_exp_sigma"
    )
  )
  expect_named(
    stop_signal(go = "exgaussian", stop = "lognormal")$p_types,
    c(
      "mu", "sigma", "tau", "meanlogS", "sdlogS", "tf", "gf",
      "exg_lb", "exg_normal_ratio", "exg_shifted_exp_sigma"
    )
  )
  expect_named(
    stop_signal(go = "exgaussian", stop = "weibull")$p_types,
    c(
      "mu", "sigma", "tau", "shapeS", "scaleS", "tf", "gf",
      "exg_lb", "exg_normal_ratio", "exg_shifted_exp_sigma"
    )
  )
  expect_named(
    stop_signal(go = "exgaussian", stop = "normal")$p_types,
    c(
      "mu", "sigma", "tau", "muS", "sigmaS", "tf", "gf",
      "exg_lb", "normS_lb", "exg_normal_ratio", "exg_shifted_exp_sigma"
    )
  )
  expect_named(
    stop_signal(go = "racing_diffusion", stop = "exgaussian")$p_types,
    c(
      "v", "B", "A", "t0", "s", "muS", "sigmaS", "tauS", "tf", "gf",
      "exgS_lb", "exg_normal_ratio", "exg_shifted_exp_sigma"
    )
  )
  expect_named(
    stop_signal(go = "racing_diffusion", stop = "lognormal")$p_types,
    c("v", "B", "A", "t0", "s", "meanlogS", "sdlogS", "tf", "gf")
  )
  expect_named(
    stop_signal(go = "racing_diffusion", stop = "weibull")$p_types,
    c("v", "B", "A", "t0", "s", "shapeS", "scaleS", "tf", "gf")
  )
  expect_named(
    stop_signal(go = "racing_diffusion", stop = "normal")$p_types,
    c("v", "B", "A", "t0", "s", "muS", "sigmaS", "tf", "gf", "normS_lb")
  )
})

test_that("ex-Gaussian approximations are off by default and opt-in", {
  model <- stop_signal(go = "exgaussian", stop = "exgaussian")
  expect_equal(
    model$p_types[c("exg_normal_ratio", "exg_shifted_exp_sigma")],
    c(exg_normal_ratio = 0, exg_shifted_exp_sigma = 0)
  )
  expect_equal(model$bound$minmax[, c("exg_normal_ratio", "exg_shifted_exp_sigma")],
               cbind(exg_normal_ratio = c(0, Inf), exg_shifted_exp_sigma = c(0, Inf)))
  expect_equal(unname(model$bound$exception["exg_normal_ratio"]), 0)
  expect_equal(unname(model$bound$exception["exg_shifted_exp_sigma"]), 0)

  expect_warning(
    stop_signal(
      go = "exgaussian",
      stop = "exgaussian",
      exgaussian_approximation = list(
        normal_ratio = .05,
        shifted_exponential_sigma = 1e-5
      )
    ),
    "changes the ex-Gaussian likelihood"
  )
  model <- suppressWarnings(stop_signal(
    go = "exgaussian",
    stop = "exgaussian",
    exgaussian_approximation = list(
      normal_ratio = .05,
      shifted_exponential_sigma = 1e-5
    )
  ))
  expect_equal(unname(model$p_types["exg_normal_ratio"]), .05)
  expect_equal(unname(model$p_types["exg_shifted_exp_sigma"]), 1e-5)
  expect_equal(unname(model$bound$exception["exg_normal_ratio"]), 0)
  expect_equal(unname(model$bound$exception["exg_shifted_exp_sigma"]), 0)

  expect_warning(
    stop_signal(
      go = "racing_diffusion",
      stop = "normal",
      exgaussian_approximation = list(normal_ratio = .05)
    ),
    "does not contain an ex-Gaussian runner"
  )

  expect_error(
    stop_signal(exgaussian_approximation = list(normal_ratio = -1)),
    "non-negative numeric scalar"
  )
  expect_error(
    stop_signal(exgaussian_approximation = list(foo = .05)),
    "Unknown `exgaussian_approximation`"
  )
})

test_that("lognormal stop sdlog has a numerical upper guardrail", {
  expect_equal(
    stop_signal(go = "exgaussian", stop = "lognormal")$bound$minmax[, "sdlogS"],
    c(1e-4, 3)
  )
  expect_equal(
    stop_signal(go = "racing_diffusion", stop = "lognormal")$bound$minmax[, "sdlogS"],
    c(1e-4, 3)
  )
})

test_that("stop_signal routes only existing C++ stop-signal backends", {
  expect_identical(stop_signal(go = "exgaussian", stop = "exgaussian")$c_name, "SSEXG")
  expect_identical(stop_signal(go = "exgaussian", stop = "lognormal")$c_name, "SSLNORM")
  expect_identical(stop_signal(go = "exgaussian", stop = "weibull")$c_name, "SSWEIBULL")
  expect_identical(stop_signal(go = "exgaussian", stop = "normal")$c_name, "SSNORM")
  expect_identical(stop_signal(go = "racing_diffusion", stop = "exgaussian")$c_name, "SSRDEX")
  expect_identical(stop_signal(go = "racing_diffusion", stop = "lognormal")$c_name, "SSRDLNORM")
  expect_identical(stop_signal(go = "racing_diffusion", stop = "weibull")$c_name, "SSRDWEIBULL")
  expect_identical(stop_signal(go = "racing_diffusion", stop = "normal")$c_name, "SSRDNORM")
})

test_that("SSEXG and SSRDEX use zero lower bounds by default", {
  expect_identical(SSEXG()$p_types[c("exg_lb", "exgS_lb")],
                   c(exg_lb = 0, exgS_lb = 0))
  expect_identical(SSRDEX()$p_types["exgS_lb"], c(exgS_lb = 0))
})

test_that("legacy, direct and wrapped SS constructors agree through fitting", {
  withr::local_options(emc.ll_backend = "multiprocess")
  settings <- list(default = NULL, normal = list(normal_ratio = 3),
                   exponential = list(shifted_exponential_sigma = .1))
  for (go in c("exgaussian", "racing_diffusion")) {
    constructor <- if (go == "exgaussian") SSEXG else SSRDEX
    for (setting in settings) {
      legacy <- suppressWarnings(constructor(exgaussian_approximation = setting))
      direct <- suppressWarnings(stop_signal(go = go, exgaussian_approximation = setting))
      for (field in c("p_types", "transform", "bound", "ss_info"))
        expect_identical(legacy[[field]], direct[[field]], info = paste(go, field))
      expect_identical(legacy$c_name, direct$c_name)

      designs <- list(
        ss_test_design(go, model = direct),
        ss_test_design(go, model = function() suppressWarnings(
          stop_signal(go = go, exgaussian_approximation = setting))),
        ss_test_design(go, model = function() suppressWarnings(
          constructor(exgaussian_approximation = setting)))
      )
      p <- make_stop_cpp_p(designs[[1]])
      likelihoods <- lapply(designs, function(des) {
        expect_true(is.function(des$model))
        expect_identical(sampled_pars(des, doMap = FALSE),
                         sampled_pars(designs[[1]], doMap = FALSE))
        emc <- make_emc(ss_test_data(), des, type = "single", n_chains = 1,
                        compress = FALSE, verbose = FALSE)
        calc_ll_manager(rbind(p, p + .01), emc[[1]]$data[[1]], emc[[1]]$model,
                         return_trialwise = TRUE)
      })
      expect_true(all(is.finite(likelihoods[[1]])))
      expect_identical(likelihoods[[1]], likelihoods[[2]])
      expect_identical(likelihoods[[1]], likelihoods[[3]])
    }
  }
})
