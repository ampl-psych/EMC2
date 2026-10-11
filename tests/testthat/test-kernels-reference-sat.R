# The reference-point option (make_kernel(reference = )), the saturation-point
# parameterisation (sat_incr / sat_decr), and the non-finite covariate handling
# of the slin_* / sat_* kernels.
#
# Non-finite covariates (Inf as well as NA) mean "this row has no covariate" for
# the slin_* / sat_* kernels, which return 0 there and leave the target
# parameter untrended on that row. Stop-signal designs rely on it: SSD is Inf on
# go trials. Every other kernel expects a finite covariate, and design() refuses
# an NA covariate for them: what a missing covariate means is the user's call.

make_tiny_design <- function(trend, ...) {
  cov_names <- trend$kernels[[1]]$cov_names
  for (nm in names(list(...))) cov_names <- c(cov_names, nm)
  design(factors = list(subjects = 1, S = 1), Rlevels = 1,
         covariates = cov_names, matchfun = function(x) x$S == x$R, trend = trend,
         formula = list(m ~ 1, s ~ 1, t0 ~ 1),
         model = LNR, report_p_vector = FALSE)
}
make_minimal_emc <- function(trend, n_trials = 6, covariate1 = 1:6, ...) {
  des <- make_tiny_design(trend, ...)
  p_vector <- sampled_pars(des, doMap = FALSE)
  dat <- make_data(p_vector, des, n_trials = n_trials,
                   covariates = data.frame(covariate1 = covariate1))
  make_emc(dat, des, type = "single")
}
cov1 <- c(0, .1, .3, 1, Inf, NA)          # two non-finite rows on purpose

test_that("slin_* / sat_* return 0 for a non-finite covariate; the other kernels refuse NA", {
  pars <- list(slin_incr = c(m.k_sat = log(2)), slin_decr = c(m.k_sat = log(2)),
               sat_incr = c(m.s_sat = log(.5)), sat_decr = c(m.s_sat = log(.5)))
  for (ty in names(pars)) {
    emc <- make_minimal_emc(make_trend(make_base("m", "lin", make_kernel("covariate1", ty))),
                            covariate1 = cov1)
    out <- as.numeric(apply_kernel(pars[[ty]], emc))
    expect_equal(out[!is.finite(cov1)], rep(0, sum(!is.finite(cov1))), info = ty)
    expect_true(all(is.finite(out)), info = ty)
  }
  for (ty in c("lin_incr", "exp_decr", "pow_incr", "poly2")) {
    expect_error(
      make_minimal_emc(make_trend(make_base("m", "lin", make_kernel("covariate1", ty))),
                       covariate1 = cov1),
      "NA value found in covariate", info = ty)
  }
})

test_that("sat_incr / sat_decr are slin_* reparameterised by the saturation point", {
  s_sat <- .5
  emc_sat <- make_minimal_emc(make_trend(make_base("m", "lin", make_kernel("covariate1", "sat_incr"))),
                              covariate1 = cov1)
  got <- as.numeric(apply_kernel(c(m.s_sat = log(s_sat)), emc_sat))
  want <- ifelse(is.finite(cov1), pmin(1, cov1 / s_sat), 0)
  expect_equal(got, want)
  # identical to slin_incr with k_sat = 1 / s_sat
  emc_slin <- make_minimal_emc(make_trend(make_base("m", "lin", make_kernel("covariate1", "slin_incr"))),
                               covariate1 = cov1)
  expect_equal(got, as.numeric(apply_kernel(c(m.k_sat = log(1 / s_sat)), emc_slin)))
  emc_dec <- make_minimal_emc(make_trend(make_base("m", "lin", make_kernel("covariate1", "sat_decr"))),
                              covariate1 = cov1)
  expect_equal(as.numeric(apply_kernel(c(m.s_sat = log(s_sat)), emc_dec)), -want)
  # the saturation point is in covariate units: the plateau starts exactly at s_sat
  emc2 <- make_minimal_emc(make_trend(make_base("m", "lin", make_kernel("covariate1", "sat_incr"))),
                           covariate1 = c(.4, .5, .6, .7, .8, .9))
  o2 <- as.numeric(apply_kernel(c(m.s_sat = log(.5)), emc2))
  expect_equal(o2, c(.8, 1, 1, 1, 1, 1))
})

test_that("reference = c0 gives k(c) - k(c0) row by row, exactly 0 at c = c0", {
  cv <- c(.1, .2, .3, .5, .8, 1.2)                         # includes the reference .3
  types <- c("lin_incr", "lin_decr", "exp_incr", "exp_decr", "pow_incr", "pow_decr",
             "poly2", "poly3", "slin_incr", "slin_decr", "sat_incr", "sat_decr")
  pars <- list(lin_incr = NULL, lin_decr = NULL,
               exp_incr = c(m.d_ei = log(2)), exp_decr = c(m.d_ed = log(2)),
               pow_incr = c(m.d_pi = log(2)), pow_decr = c(m.d_pd = log(2)),
               poly2 = c(m.d1 = .5, m.d2 = .1), poly3 = c(m.d1 = .5, m.d2 = .1, m.d3 = .01),
               slin_incr = c(m.k_sat = log(2)), slin_decr = c(m.k_sat = log(2)),
               sat_incr = c(m.s_sat = log(.5)), sat_decr = c(m.s_sat = log(.5)))
  for (ty in types) {
    tr_raw <- make_trend(make_base("m", "lin", make_kernel("covariate1", ty)))
    tr_ref <- make_trend(make_base("m", "lin", make_kernel("covariate1", ty, reference = .3)))
    raw <- as.numeric(apply_kernel(pars[[ty]], make_minimal_emc(tr_raw, covariate1 = cv)))
    at0 <- as.numeric(apply_kernel(pars[[ty]], make_minimal_emc(tr_raw, covariate1 = rep(.3, 6))))
    ref <- as.numeric(apply_kernel(pars[[ty]], make_minimal_emc(tr_ref, covariate1 = cv)))
    expect_equal(ref, raw - at0, info = ty)
    expect_equal(ref[cv == .3], 0, info = ty)
  }
  # slin_* / sat_*: non-finite rows stay at 0 (no trend), finite rows are shifted
  for (ty in c("slin_incr", "sat_decr")) {
    tr_ref <- make_trend(make_base("m", "lin", make_kernel("covariate1", ty, reference = .3)))
    out <- as.numeric(apply_kernel(pars[[ty]], make_minimal_emc(tr_ref, covariate1 = cov1)))
    expect_equal(out[!is.finite(cov1)], rep(0, sum(!is.finite(cov1))), info = ty)
    expect_true(all(is.finite(out)), info = ty)
  }
})

test_that("reference gives the same trend on the one-pass and trial-by-trial simulators", {
  # the shift is a per-row constant, so it cannot depend on which trials exist yet
  tr <- make_trend(make_base("m", "lin", make_kernel("covariate1", "sat_incr", reference = .3)))
  des <- design(factors = list(subjects = 1, S = 1:2), Rlevels = 1:2, covariates = "covariate1",
                matchfun = function(x) x$S == x$lR, trend = tr,
                formula = list(m ~ lM, s ~ 1, t0 ~ 1), model = LNR, report_p_vector = FALSE)
  p <- sampled_pars(des, doMap = FALSE)
  p[c("m", "m_lMTRUE", "s", "t0", "m.w", "m.s_sat")] <- c(-.5, .3, log(.3), log(.2), .5, log(.4))
  cvals <- c(.1, .2, .3, .5, .8, 1.2)
  cv <- data.frame(covariate1 = rep(cvals, 2))              # n_trials is per design cell
  set.seed(1)
  one <- make_data(p, des, n_trials = 6, covariates = cv, return_trialwise_parameters = TRUE)
  set.seed(1)
  tbt <- make_data(p, des, n_trials = 6, covariates = cv, return_trialwise_parameters = TRUE,
                   conditional_on_data = FALSE)
  tw_one <- attr(one, "trialwise_parameters")
  tw_tbt <- attr(tbt, "trialwise_parameters")
  kcol <- "k1.covariate1.covariate"
  expect_equal(unname(tw_tbt[, "m"]), unname(tw_one[, "m"]))
  expect_equal(unname(tw_tbt[, kcol]), unname(tw_one[, kcol]))
  # both equal the closed form min(1, c/s) - min(1, .3/s), a per-row function of c
  c_rows <- one$covariate1[rep(seq_len(nrow(one)), each = nrow(tw_one) / nrow(one))]
  expect_equal(unname(tw_one[, kcol]), pmin(1, c_rows / .4) - min(1, .3 / .4), tolerance = 1e-12)
})

test_that("reference is validated and refused where it has no meaning", {
  expect_error(make_kernel("covariate1", "delta", reference = .3), "not supported for the sequential")
  expect_error(make_kernel("covariate1", "lin_incr", reference = NA), "single finite number")
  expect_error(make_kernel("covariate1", "lin_incr", reference = c(.1, .2)), "single finite number")
  expect_error(make_kernel("covariate1", "lin_incr", reference = "a"), "single finite number")
  expect_error(make_kernel(c("covariate1", "covariate2"), "lin_incr", reference = .3), "exactly one covariate")
  expect_null(make_kernel("covariate1", "lin_incr")$reference)
  expect_equal(make_kernel("covariate1", "lin_incr", reference = .3)$reference, .3)
})

test_that("a saturation point can be bounded with the pnorm transform", {
  # s_sat bounded to the covariate range: the sampled parameter cannot reach the
  # flat regions where the kernel is constant and the weight is unidentified
  tr <- make_trend(make_base("m", "lin", make_kernel("covariate1", "sat_incr")))
  des <- design(factors = list(subjects = 1, S = 1), Rlevels = 1, covariates = "covariate1",
                matchfun = function(x) x$S == x$R, trend = tr,
                formula = list(m ~ 1, s ~ 1, t0 ~ 1), model = LNR, report_p_vector = FALSE,
                transform = list(func = c(m.s_sat = "pnorm"),
                                 lower = c(m.s_sat = .05), upper = c(m.s_sat = .40)))
  expect_equal(unname(des$model()$transform$func["m.s_sat"]), "pnorm")
  expect_equal(unname(des$model()$transform$lower["m.s_sat"]), .05)
  expect_equal(unname(des$model()$transform$upper["m.s_sat"]), .40)
  # the sampled scale maps onto [lower, upper]: extreme draws stay inside
  emc <- make_minimal_emc(tr, covariate1 = c(.02, .05, .1, .2, .3, .4))
  emc_b <- make_emc(get_data(emc), des, type = "single")
  for (z in c(-4, 0, 4)) {
    k <- as.numeric(apply_kernel(c(m.s_sat = z), emc_b))
    # k = min(1, c/s) with s in [.05,.40] => k is bounded by c/.05 and c/.40
    expect_true(all(k <= pmin(1, c(.02, .05, .1, .2, .3, .4) / .05) + 1e-9))
    expect_true(all(k >= pmin(1, c(.02, .05, .1, .2, .3, .4) / .40) - 1e-9))
  }
})
