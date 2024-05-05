test_that("compare", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  expect_snapshot(
    compare(list(samplers_LNR))
  )
})

test_that("savage-dickey", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  expect_snapshot(
    savage_dickey(samplers_LNR, parameter = "m", do_plot = F))
  expect_snapshot(
    savage_dickey(samplers_LNR, fun = function(d) d["m"] - d["m_lMd"],
                  H0 = -0.5, do_plot = F))
})
