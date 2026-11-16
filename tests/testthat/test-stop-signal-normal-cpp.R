RNGkind("L'Ecuyer-CMRG")
set.seed(123)

make_ssnorm_cpp_design <- function() {
  design(
    model = function() stop_signal(go = "exgaussian", stop = "normal"),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(mu ~ lM, sigma ~ 1, tau ~ 1,
                   muS ~ 1, sigmaS ~ 1, gf ~ 1, tf ~ 1)
  )
}

make_ssrdnorm_cpp_design <- function() {
  design(
    model = function() stop_signal(go = "racing_diffusion", stop = "normal"),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1,
                   muS ~ 1, sigmaS ~ 1, gf ~ 1, tf ~ 1)
  )
}


test_that("SSNORM C++ likelihood matches R reference path", {
  design_ssn <- make_ssnorm_cpp_design()
  p_ssn <- make_stop_cpp_p(
    design_ssn,
    c(
      mu = log(.6),
      mu_lMTRUE = log(.8),
      sigma = log(.05),
      tau = log(.2),
      muS = log(.24),
      sigmaS = log(.05),
      gf = qnorm(.1),
      tf = qnorm(.1)
    )
  )

  expect_identical(design_ssn$model()$c_name, "SSNORM")
  expect_stop_cpp_equal_r(design_ssn, p_ssn)
})

test_that("SSRDNORM C++ likelihood matches R reference path", {
  design_ssrdn <- make_ssrdnorm_cpp_design()
  p_ssrdn <- make_stop_cpp_p(
    design_ssrdn,
    c(
      v = log(1.8),
      v_lMTRUE = log(.35),
      B = log(.8),
      A = log(.10),
      t0 = log(.20),
      s = log(1),
      muS = log(.24),
      sigmaS = log(.05),
      gf = qnorm(.1),
      tf = qnorm(.1)
    )
  )

  expect_identical(design_ssrdn$model()$c_name, "SSRDNORM")
  expect_stop_cpp_equal_r(design_ssrdn, p_ssrdn)
})
