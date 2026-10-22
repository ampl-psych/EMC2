# The scaled DDM card (kind "flow_scaled", DDMnle): register_nn_model() wires
# it without a scaffold, its R path is the NLE project's reference
# (tests/testthat/port/ddm_scaled.R, ddm_sc_density), the compiled path
# (src/model_NN.h, nn_trial_ll_scaled) equals the R path, and the three SZ
# conventions meet at their boundaries.
#
# Two cards are used: a random-weight fixture (always; its numbers model
# nothing, the agreement is the test) and, when EMC2_NLE_SCALED_CARD names
# the NLE project's ddm_sc_s1106_e245_weights.json, the real DDMnle card
# with its golden vectors.

nle <- function(f) getFromNamespace(f, "EMC2")
port <- new.env()
sys.source(test_path("port", "ddm_scaled.R"), envir = port)

# A scaled card with random weights: inputs kappa, eta, Z, SZ (eta and SZ
# zero-box), the DDMnle box and spline range.
scaled_fixture <- function(seed = 1) {
  set.seed(seed)
  layer <- function(n_in, n_out, sd = 1 / sqrt(n_in))
    list(W = matrix(rnorm(n_in * n_out, sd = sd), n_in), b = rnorm(n_out, sd = .1))
  zc <- list(SZ = 0.05, eta = 0.01)
  nat_lo <- c(kappa = -10, eta = 0, Z = 0.1, SZ = 0); nat_hi <- c(kappa = 10, eta = 8, Z = 0.9, SZ = 0.99)
  tf <- c(kappa = "identity", eta = "zerobox", Z = "probit", SZ = "zerobox")
  to_net <- function(x, p) switch(tf[[p]], identity = x, probit = qnorm(x), zerobox = asinh(x / zc[[p]]))
  list(card_version = 1, model = "DDMsc", budget = "fixture", checkpoint_step_flow = 0,
       context_names = names(tf), classifier_context_names = names(tf),
       context_encoding = "zerobox", zerobox_c = zc,
       scaling = list(kind = "ddm_unit", pars = c("v", "a", "sv", "t0", "st0", "s", "Z", "SZ"), sz = "proportion"),
       context_transforms = as.list(tf), parameter_scale = "sampled",
       bounds_natural = list(lower = unname(nat_lo), upper = unname(nat_hi)),
       bounds_sampled = list(lower = unname(mapply(to_net, nat_lo, names(tf))),
                             upper = unname(mapply(to_net, nat_hi, names(tf)))),
       scaler = list(mean = c(0, 1, 0, 1), scale = c(4, 1.5, 1, 1.2)),
       classifier_scaler = list(mean = c(0, 1, 0, 1), scale = c(4, 1.5, 1, 1.2)),
       flow_mlp = list(activation = "gelu_tanh", use_norm = FALSE,
                       layers = list(layer(5, 32), layer(32, 16), layer(16, 3 * 64 + 1, sd = .05)), norms = NULL),
       classifier_mlp = list(activation = "gelu_tanh", layers = list(layer(4, 16), layer(16, 1))),
       spline = list(num_bins = 64, num_splines = 1, range_min = -12, range_max = 2, min_bin_size = 1e-4,
                     min_knot_slope = 1e-4, boundary_slopes = "identity", output_transform = "exp",
                     base_distribution = "standard_normal"),
       support = list(time = c(0, 1e6)), mass_tolerance = 1e-3,
       beyond_support = "identity tails beyond the spline range in log u carry the remaining mass",
       validation = list(note = "random-weight fixture; nothing validated"),
       ll_floor_log = log(1e-10))
}
# written as JSON, the card format the exporter writes and the NLE reference
# evaluator (port$ddm_sc_card_evaluator) reads
write_card <- function(card, f = tempfile(fileext = ".json")) {
  jsonlite::write_json(card, f, auto_unbox = TRUE, digits = NA)
  f
}
fixture_file <- write_card(scaled_fixture())

real_card <- Sys.getenv("EMC2_NLE_SCALED_CARD", "")
cards <- c(fixture = fixture_file, if (nzchar(real_card) && file.exists(real_card)) c(DDMnle = real_card))

# Natural-scale DDM parameters in the registered model's convention (SZ the
# proportion), with exact zeros in the chosen regime, inside the network's box.
draw_pars <- function(n, regime = character(0), seed = 7) {
  set.seed(seed)
  out <- NULL
  while (is.null(out) || nrow(out) < n) {
    Z <- runif(1, .2, .8)
    p <- c(v = runif(1, -3, 3), a = runif(1, .4, 2.5), sv = runif(1, .05, 2), t0 = runif(1, .1, .5),
           st0 = runif(1, .001, .3), s = runif(1, .3, 2), Z = Z, SZ = runif(1, .05, .9))
    p[regime] <- 0
    if (abs(p[["v"]] * p[["a"]] / p[["s"]]^2) > 9.5 || p[["sv"]] * p[["a"]] / p[["s"]]^2 > 7.5) next
    out <- rbind(out, p)
  }
  rownames(out) <- NULL
  out
}
with_width <- function(pars) { pars[, "SZ"] <- nn_sz_convert(pars[, "SZ"], pars[, "Z"], "proportion", "width"); pars }
regimes <- list(full = character(0), sv0 = "sv", SZ0 = "SZ", st00 = "st0", st0only = c("sv", "SZ"),
                SZonly = c("sv", "st0"), all0 = c("sv", "SZ", "st0"))

test_that("a scaled card registers as flow_scaled with the DDM's parameterisation and no scaffold", {
  for (nm in names(cards)) {
    m <- register_nn_model(cards[[nm]], twin = DDM)
    ml <- m(); reg <- ml$nn
    expect_identical(reg$kind, "flow_scaled")
    expect_true(reg$scaled)
    expect_identical(reg$pars, c("kappa", "eta", "Z", "SZ"))
    expect_identical(reg$sc_pars, c("v", "a", "sv", "t0", "st0", "s", "Z", "SZ"))
    expect_identical(reg$scaling$kind, "ddm_unit")
    expect_identical(ml$p_types, DDM()$p_types)
    expect_identical(ml$transform, DDM()$transform)
    expect_identical(ml$bound$minmax, DDM()$bound$minmax)
    expect_equal(ml$bound$exception, DDM()$bound$exception)
    expect_identical(ml$c_name, "NN")
    expect_identical(reg$native$kind, "flow_joint")
    expect_true(reg$native$scaled)
    expect_identical(reg$native$sc_pars, reg$sc_pars)
    expect_false(reg$cdf)
    expect_length(reg$refuse_default, 0)
    expect_identical(reg$twin, "DDM")
  }
  expect_error(register_nn_model(fixture_file, twin = DDM, cdf = TRUE), "no CDF of rt")
  expect_error(register_nn_model(fixture_file, twin = DDM, pre = "t0"), "pre does not apply")
  expect_error(register_nn_model(fixture_file, twin = DDM, exceptions = c(st0 = 0)), "exceptions do not apply")
  expect_error(register_nn_model(fixture_file, twin = DDM, pars = c("k", "e", "Z", "SZ")), "cannot be renamed")
  bad <- scaled_fixture(); bad$context_names <- c("eta", "kappa", "Z", "SZ")
  bad$context_transforms <- bad$context_transforms[c("eta", "kappa", "Z", "SZ")]
  bad$bounds_natural <- lapply(bad$bounds_natural, function(x) x[c(2, 1, 3, 4)])
  bad$bounds_sampled <- lapply(bad$bounds_sampled, function(x) x[c(2, 1, 3, 4)])
  expect_error(register_nn_model(write_card(bad), twin = DDM), "kappa, eta, Z, SZ in that order")
})

test_that("the R path is the NLE reference (ddm_sc_density, SZ as the width) to 1e-12, in every regime", {
  for (nm in names(cards)) {
    m <- register_nn_model(cards[[nm]], twin = DDM); ml <- m()
    ev <- port$ddm_sc_card_evaluator(cards[[nm]])
    for (rg in names(regimes)) {
      pars <- draw_pars(60, regimes[[rg]], seed = 11)
      rt <- pars[, "t0"] + c(-0.01, 0, 0.002, 0.03, 0.2, 0.6)[(seq_len(60) %% 6) + 1] * (pars[, "a"] / pars[, "s"])^2 +
        pars[, "st0"] * runif(60)
      R <- sample(1:2, 60, TRUE)
      ours <- ml$dfun(rt, R, pars)
      ref <- port$ddm_sc_density(rt, R, with_width(pars), ev)
      expect_lt(max(abs(ours - ref)), 1e-12, label = paste(nm, rg))
      expect_true(all(ours[rt <= pars[, "t0"]] == 0))
    }
  }
})

test_that("with the network replaced by the analytic unit-scale DDM, the likelihood is dDDM with st0 (oracle)", {
  m <- register_nn_model(fixture_file, twin = DDM); reg <- m()$nn
  set.seed(20261013)
  worst0 <- worst <- 0
  for (k in 1:120) {
    regime <- regimes[[(k %% length(regimes)) + 1]]
    p <- draw_pars(1, regime, seed = 1000 + k)[1, ]
    if (k %% 2 == 0) p[["st0"]] <- 0
    sc <- (p[["a"]] / p[["s"]])^2
    rt <- p[["t0"]] + c(-0.01, 0, p[["st0"]] / 3, p[["st0"]], sc * rexp(6, 1 / .3) + p[["st0"]])
    R <- factor(sample(1:2, length(rt), TRUE), levels = 1:2)
    P <- matrix(p, length(rt), length(p), byrow = TRUE, dimnames = list(NULL, names(p)))
    a <- nle("dDDM")(rt, R, with_width(P), precision = 1e-8); a[rt <= p[["t0"]]] <- 0
    o <- nle("nn_eval_scaled")(reg, rt, as.integer(R), P, "pdf", evaluator = port$ddm_sc_oracle)
    big <- a > 1
    d <- abs(log(o[big]) - log(a[big]))
    if (p[["st0"]] == 0) worst0 <- max(worst0, d) else worst <- max(worst, d)
  }
  expect_lt(worst0, 1e-9)      # exact rescaling and shift (the NLE test: 8e-12)
  expect_lt(worst, 1e-4)       # the CDF difference against dDDM's own st0 integration (NLE: <= 5e-5)
})

# The compiled path through a design, as a fit evaluates it, against the R
# path, in the six regimes of the NLE project's eight-case review.
scaled_cell <- function(m, regime) {
  reg <- m()$nn
  PM <- c(v = 1, a = log(.8), t0 = log(.3), Z = qnorm(.5), sv = log(.1), SZ = qnorm(.2), st0 = log(.05))
  PS <- c(v = .15, a = .15, t0 = .1, Z = .05, sv = .15, SZ = .15, st0 = .15)
  pin <- stats::setNames(rep(-Inf, length(regime)), regime)
  # For scaled cards: parameters in reg$fixed (e.g. st0=0) must be constants, not sampled
  fixed_pin <- if (length(reg$fixed)) stats::setNames(rep(-Inf, length(reg$fixed)), names(reg$fixed)) else NULL
  all_pin <- c(pin, fixed_pin[setdiff(names(fixed_pin), names(pin))])
  free <- setdiff(names(PM), names(all_pin))
  suppressWarnings(nn_cell(m, mean = PM[free], sd = PS[free],
                           formula = lapply(free, function(p) stats::as.formula(paste(p, "~ 1"))),
                           constants = if (length(all_pin)) all_pin else NULL,
                           factors = list(S = c("left", "right")), Rlevels = c("left", "right")))
}

test_that("the compiled likelihood equals the R path to 1e-11 over prior draws in each regime", {
  skip_on_cran()
  for (nm in names(cards)) {
    m <- register_nn_model(cards[[nm]], twin = DDM)
    for (rg in c("full", "sv0", "SZ0", "st00", "st0only", "SZonly")) {
      cell <- scaled_cell(m, regimes[[rg]])
      # the R path: the design's prepared model (it carries pre_transform) with
      # the compiled likelihood switched off
      prepared <- cell$design$model()
      r_path <- function() { ml <- prepared; ml$c_name <- NULL; ml }
      set.seed(7)
      dat <- suppressMessages(make_data(cell$mean, cell$control_design, n_trials = 100))
      X <- nle("nn_cell_draws")(cell, if (nm == "DDMnle") 200 else 60)
      dadm <- nle("nn_dadm")(dat, cell$design)
      native <- as.numeric(nle("calc_ll_manager")(X, dadm, cell$design$model))
      rp <- as.numeric(nle("calc_ll_manager")(X, dadm, r_path))
      expect_true(all(is.finite(native)), label = paste(nm, rg))
      expect_lt(max(abs(native - rp)), 1e-11, label = paste(nm, rg))
      # the mapped values of the pinned parameters are exactly zero
      mp <- mapped_pars(cell$design, cell$mean)
      for (p in regimes[[rg]]) expect_true(all(mp[, p] == 0))
    }
  }
})

test_that("the SZ conventions meet at their three boundaries (golden rows)", {
  m <- register_nn_model(cards[[1]], twin = DDM); ml <- m()
  ev <- port$ddm_sc_card_evaluator(cards[[1]])
  grid <- expand.grid(Z = c(.1, .3, .5, .7, .9), SZ = c(0, .05, .2, .5, .99))
  # 1. EMC2's DDM(): proportion -> width through Ttransform; the conversion inverts it exactly
  pars <- cbind(v = 1, a = 1, sv = .1, t0 = .2, st0 = 0, s = 1, Z = grid$Z, SZ = grid$SZ)
  width <- DDM()$Ttransform(pars, NULL)[, "SZ"]
  expect_equal(width, nn_sz_convert(grid$SZ, grid$Z, "proportion", "width"), tolerance = 1e-15)
  expect_equal(nn_sz_convert(width, grid$Z, "width", "proportion"), grid$SZ, tolerance = 1e-15)
  # 2. the registered model: Ttransform leaves the proportion as it is
  expect_identical(ml$Ttransform(pars, NULL)[, "SZ"], pars[, "SZ"])
  # 3. the network: the registered model (proportion) equals the reference
  #    evaluated on the width, and the width-based inputs are the proportion again
  rt <- rep(0.5, nrow(pars)); R <- rep(1:2, length.out = nrow(pars))
  expect_lt(max(abs(ml$dfun(rt, R, pars) - port$ddm_sc_density(rt, R, with_width(pars), ev))), 1e-12)
  U_ref <- port$ddm_sc_unit(with_width(pars))
  expect_equal(U_ref[, "SZ"], nle("nn_scaled_inputs")(pars)[, "SZ"], tolerance = 1e-14)
})

test_that("DDMnle's golden vectors and total mass (real card only)", {
  skip_if(!("DDMnle" %in% names(cards)), "EMC2_NLE_SCALED_CARD not set")
  g <- nn_golden_check(cards[["DDMnle"]], twin = DDM, tol = 1e-9)
  expect_true(g$pass)
  expect_output(print(g), "^PASS")
})

test_that("a card without card_version is version 0; a version this EMC2 does not read is refused by name", {
  c0 <- scaled_fixture(); c0$card_version <- NULL
  expect_equal(register_nn_model(write_card(c0), twin = DDM)()$nn$card_version, 0)
  expect_equal(register_nn_model(fixture_file, twin = DDM)()$nn$card_version, 1)
  c2 <- scaled_fixture(); c2$card_version <- 2
  expect_error(register_nn_model(write_card(c2), twin = DDM), "Card version 2 is not supported")
  # the same through an .rds bundle
  f2 <- tempfile(fileext = ".rds"); saveRDS(c2, f2)
  expect_error(register_nn_model(f2, twin = DDM), "Card version 2 is not supported")
  # a version-1 card without a validation record says so when registered
  c1 <- scaled_fixture(); c1$validation <- NULL
  expect_message(register_nn_model(write_card(c1), twin = DDM), "records no validation")
})
