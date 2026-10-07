# lik_precision() (R/sampling.R): the quadratic likelihood surrogate that the
# subject step's local proposal and the interweaving sweep run on. Its step
# search must settle, and a likelihood floored at a model bound must get no
# precision rather than the cliff read as a huge curvature.

# A custom likelihood on sufficient statistics, parameters (g, h, c): a
# Gaussian in g (precision T), a heavy-tailed curve in h (flat near its mode,
# steep further out, so a multiplicative step search overshoots and
# undershoots), and c, flat below a bound and floored beyond it.
lp_ll <- function(pars, dadm, ...) {
  if (pars[3] >= dadm$bound[1]) return(dadm$floor[1])
  -0.5 * dadm$T[1] * (dadm$ybar[1] - pars[1])^2 - 40 * log1p((pars[2] / .05)^2)
}
lp_make <- function(floor) {
  dat <- data.frame(subjects = factor(1:2), T = 400, ybar = .3, bound = 1, floor = floor)
  des <- design(model = lp_ll, custom_p_vector = c("g", "h", "c"), report_p_vector = FALSE)
  make_emc(dat, des, type = "standard", n_chains = 1, compress = FALSE, verbose = FALSE)[[1]]
}

test_that("lik_precision settles its step and reads a cliff as no precision", {
  s <- lp_make(-1e4)
  centre <- c(g = .3, h = 0, c = 0)
  for (h0 in c(.001, .1, 3)) {        # whatever the posterior scale handed in
    lik <- EMC2:::lik_precision(centre, rep(h0, 3), s$data[[1]], s$model)
    # the Gaussian parameter: exact precision and linear term
    expect_equal(lik$prec["g", "g"], 400, tolerance = 1e-3)
    expect_equal(lik$lin[["g"]], 400 * .3, tolerance = 1e-3)
    # the heavy-tailed one: a finite curvature of the right order (2 * 40 / .05^2 = 32000 at the mode)
    expect_gt(lik$prec["h", "h"], 32000 / 3)
    expect_lt(lik$prec["h", "h"], 32000 * 3)
    # the floored one: nothing, including its cross terms
    expect_equal(unname(lik$prec["c", ]), c(0, 0, 0))
    expect_equal(unname(lik$prec[, "c"]), c(0, 0, 0))
    expect_equal(lik$lin[["c"]], 0)
  }
})

test_that("lik_precision returns NULL off the likelihood's support", {
  s <- lp_make(-Inf)
  expect_null(EMC2:::lik_precision(c(g = .3, h = 0, c = 2), rep(.1, 3), s$data[[1]], s$model))
})
