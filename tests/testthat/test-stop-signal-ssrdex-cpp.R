RNGkind("L'Ecuyer-CMRG")
set.seed(123)

make_ssrdex_cpp_design <- function() {
  design(
    model = function() stop_signal(go = "racing_diffusion", stop = "exgaussian"),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1,
                   muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
  )
}


test_that("SSRDEX C++ likelihood matches R reference path", {
  design_ssrdex <- make_ssrdex_cpp_design()
  p_ssrdex <- make_stop_cpp_p(
    design_ssrdex,
    c(
      v = log(1.8),
      v_lMTRUE = log(.35),
      B = log(.8),
      A = log(.10),
      t0 = log(.20),
      s = log(1),
      muS = log(.17),
      sigmaS = log(.025),
      tauS = log(.05),
      gf = qnorm(.1),
      tf = qnorm(.1)
    )
  )

  expect_identical(design_ssrdex$model()$c_name, "SSRDEX")
  expect_stop_cpp_equal_r(design_ssrdex, p_ssrdex)
})
