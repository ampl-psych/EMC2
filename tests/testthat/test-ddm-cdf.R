# Tests for the C++ defective DDM CDF (p_DDM_Wien / p_DDM_Wien_R), validated
# against the independent WienR::pWDM reference (which the R pDDM() wraps).
# pars columns are the DDM p_types order: v, a, sv, t0, st0, s, Z, SZ.

mk_pars <- function(v, a, sv, t0, st0, s, Z, SZ, n) {
  matrix(c(v, a, sv, t0, st0, s, Z, SZ), nrow = n, ncol = 8, byrow = TRUE)
}

# cpp log-CDF -> CDF for response `resp` (1 = lower, 2 = upper) at times `rt`
cpp_cdf <- function(rt, resp, v, a, sv, t0, st0, s, Z, SZ) {
  n <- length(rt)
  exp(p_DDM_Wien_R(rt, rep(as.integer(resp), n),
                   mk_pars(v, a, sv, t0, st0, s, Z, SZ, n), rep(TRUE, n)))
}

wienr_cdf <- function(rt, resp, v, a, sv, t0, st0, s, Z, SZ) {
  rr <- if (resp == 1) "lower" else "upper"
  WienR::pWDM(rt, response = rr, a = a / s, v = v / s, t0 = t0, w = Z,
              sw = SZ, sv = sv / s, st0 = st0, precision = 5e-3)$value
}

rt <- c(0.35, 0.5, 0.7, 1.0, 1.5)

test_that("p_DDM_Wien matches WienR::pWDM with no between-trial variability", {
  grid <- list(
    c(v =  1.0, a = 1.2, t0 = 0.20, Z = 0.50),
    c(v = -0.5, a = 1.5, t0 = 0.15, Z = 0.40),
    c(v =  2.0, a = 0.8, t0 = 0.25, Z = 0.60),
    c(v =  0.0, a = 1.0, t0 = 0.10, Z = 0.50)
  )
  for (p in grid) for (resp in 1:2) {
    cpp <- cpp_cdf(rt, resp, p["v"], p["a"], 0, p["t0"], 0, 1, p["Z"], 0)
    ref <- wienr_cdf(rt, resp, p["v"], p["a"], 0, p["t0"], 0, 1, p["Z"], 0)
    expect_equal(cpp, ref, tolerance = 1e-8)
  }
})

test_that("p_DDM_Wien matches WienR::pWDM with between-trial variability", {
  grid <- list(
    c(v = 1.0, a = 1.2, sv = 0.5, t0 = 0.20, st0 = 0.00, Z = 0.50, SZ = 0.0),
    c(v = 0.8, a = 1.3, sv = 0.4, t0 = 0.20, st0 = 0.05, Z = 0.55, SZ = 0.1),
    c(v = 1.2, a = 1.0, sv = 0.0, t0 = 0.20, st0 = 0.08, Z = 0.50, SZ = 0.0)
  )
  for (p in grid) for (resp in 1:2) {
    cpp <- cpp_cdf(rt, resp, p["v"], p["a"], p["sv"], p["t0"], p["st0"], 1, p["Z"], p["SZ"])
    ref <- wienr_cdf(rt, resp, p["v"], p["a"], p["sv"], p["t0"], p["st0"], 1, p["Z"], p["SZ"])
    # variability integration at precision 5e-3 -> ~1e-3 ABSOLUTE agreement
    # (relative error is larger at small CDF values, so check on the probability scale)
    expect_lt(max(abs(cpp - ref)), 1e-3)
  }
})

test_that("p_DDM_Wien honours the scaling parameter s", {
  # dividing v, a, sv by s must reproduce the s = 1 CDF on the rescaled params
  cpp_s  <- cpp_cdf(rt, 1, v = 2.0, a = 2.4, sv = 0, t0 = 0.2, st0 = 0, s = 2, Z = 0.5, SZ = 0)
  cpp_s1 <- cpp_cdf(rt, 1, v = 1.0, a = 1.2, sv = 0, t0 = 0.2, st0 = 0, s = 1, Z = 0.5, SZ = 0)
  expect_equal(cpp_s, cpp_s1, tolerance = 1e-10)
})

test_that("DDM CDF is a valid defective distribution", {
  # monotone non-decreasing in t, in [0, 1], and the two response CDFs sum to 1
  tgrid <- seq(0.3, 8, length.out = 40)
  lo <- cpp_cdf(tgrid, 1, 1.0, 1.2, 0, 0.2, 0, 1, 0.5, 0)
  hi <- cpp_cdf(tgrid, 2, 1.0, 1.2, 0, 0.2, 0, 1, 0.5, 0)
  expect_true(all(diff(lo) >= -1e-10))
  expect_true(all(diff(hi) >= -1e-10))
  expect_true(all(lo >= -1e-12 & lo <= 1 + 1e-12))
  expect_true(all(hi >= -1e-12 & hi <= 1 + 1e-12))
  # total mass -> 1 as t grows large
  expect_equal(lo[length(tgrid)] + hi[length(tgrid)], 1, tolerance = 1e-4)
})

test_that("p_DDM_Wien returns -Inf for rt <= t0 and respects is_ok", {
  pars <- mk_pars(1.0, 1.2, 0, 0.3, 0, 1, 0.5, 0, 3)
  out <- p_DDM_Wien_R(c(0.1, 0.2, 0.5), c(1L, 1L, 1L), pars, rep(TRUE, 3))
  expect_true(is.infinite(out[1]) && out[1] < 0)  # rt < t0
  expect_true(is.infinite(out[2]) && out[2] < 0)  # rt < t0
  expect_true(is.finite(out[3]))
  out2 <- p_DDM_Wien_R(c(0.5, 0.5), c(1L, 1L), mk_pars(1, 1.2, 0, 0.2, 0, 1, 0.5, 0, 2),
                       c(TRUE, FALSE))
  expect_true(is.finite(out2[1]) && is.infinite(out2[2]))
})
