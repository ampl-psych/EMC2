RNGkind("L'Ecuyer-CMRG")
set.seed(123)

make_sslnorm_cpp_design <- function() {
  design(
    model = function() stop_signal(go = "exgaussian", stop = "lognormal"),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(mu ~ lM, sigma ~ 1, tau ~ 1,
                   meanlogS ~ 1, sdlogS ~ 1, gf ~ 1, tf ~ 1)
  )
}

make_ssrdlnorm_cpp_design <- function() {
  design(
    model = function() stop_signal(go = "racing_diffusion", stop = "lognormal"),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1,
                   meanlogS ~ 1, sdlogS ~ 1, gf ~ 1, tf ~ 1)
  )
}


test_that("SSLNORM C++ likelihood matches R reference path", {
  design_ssl <- make_sslnorm_cpp_design()
  p_ssl <- make_stop_cpp_p(
    design_ssl,
    c(
      mu = log(.6),
      mu_lMTRUE = log(.8),
      sigma = log(.05),
      tau = log(.2),
      meanlogS = log(.22),
      sdlogS = log(.25),
      gf = qnorm(.1),
      tf = qnorm(.1)
    )
  )

  expect_identical(design_ssl$model()$c_name, "SSLNORM")
  expect_stop_cpp_equal_r(design_ssl, p_ssl)
})

test_that("SSRDLNORM C++ likelihood matches R reference path", {
  design_ssrdl <- make_ssrdlnorm_cpp_design()
  p_ssrdl <- make_stop_cpp_p(
    design_ssrdl,
    c(
      v = log(1.8),
      v_lMTRUE = log(.35),
      B = log(.8),
      A = log(.10),
      t0 = log(.20),
      s = log(1),
      meanlogS = log(.22),
      sdlogS = log(.25),
      gf = qnorm(.1),
      tf = qnorm(.1)
    )
  )

  expect_identical(design_ssrdl$model()$c_name, "SSRDLNORM")
  expect_stop_cpp_equal_r(design_ssrdl, p_ssrdl)
})
