RNGkind("L'Ecuyer-CMRG")
set.seed(123)

test_that("compare", {
  expect_snapshot(
    compare(list(samplers_LNR), cores_for_props = 1)
  )
})

test_that("savage-dickey", {
  expect_snapshot(
    round(savage_dickey(samplers_LNR, parameter = "m", do_plot = F, H0 = -1), 3))
  expect_snapshot(
    round(savage_dickey(samplers_LNR, fun = function(d) d["m"] - d["m_lMd"],
                  H0 = -0.5, do_plot = F), 3))
})
