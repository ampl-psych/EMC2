RNGkind("L'Ecuyer-CMRG")
set.seed(123)

make_ssweibull_cpp_design <- function() {
  design(
    model = function() stop_signal(go = "exgaussian", stop = "weibull"),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(mu ~ lM, sigma ~ 1, tau ~ 1,
                   shapeS ~ 1, scaleS ~ 1, gf ~ 1, tf ~ 1)
  )
}

make_ssrdweibull_cpp_design <- function() {
  design(
    model = function() stop_signal(go = "racing_diffusion", stop = "weibull"),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1,
                   shapeS ~ 1, scaleS ~ 1, gf ~ 1, tf ~ 1)
  )
}


test_that("SSWEIBULL C++ likelihood matches R reference path", {
  design_ssw <- make_ssweibull_cpp_design()
  p_ssw <- make_stop_cpp_p(
    design_ssw,
    c(
      mu = log(.6),
      mu_lMTRUE = log(.8),
      sigma = log(.05),
      tau = log(.2),
      shapeS = log(2.5),
      scaleS = log(.22),
      gf = qnorm(.1),
      tf = qnorm(.1)
    )
  )

  expect_identical(design_ssw$model()$c_name, "SSWEIBULL")
  expect_stop_cpp_equal_r(design_ssw, p_ssw)
})

test_that("SSRDWEIBULL C++ likelihood matches R reference path", {
  design_ssrdw <- make_ssrdweibull_cpp_design()
  p_ssrdw <- make_stop_cpp_p(
    design_ssrdw,
    c(
      v = log(1.8),
      v_lMTRUE = log(.35),
      B = log(.8),
      A = log(.10),
      t0 = log(.20),
      s = log(1),
      shapeS = log(2.5),
      scaleS = log(.22),
      gf = qnorm(.1),
      tf = qnorm(.1)
    )
  )

  expect_identical(design_ssrdw$model()$c_name, "SSRDWEIBULL")
  expect_stop_cpp_equal_r(design_ssrdw, p_ssrdw)
})
