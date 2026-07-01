RNGkind("L'Ecuyer-CMRG")

# ── Purpose ────────────────────────────────────────────────────────────────────
# Independent verification that the C++ calc_ll path for SSEXG produces values
# that match analytically-derived formulas for every distinct trial type.
#
# THIS FILE TESTS THE CODE, not the math.  For each of the four trial types
# (go response, stop-trial go response, NR-go, NR-stop) we:
#   1. construct a single-trial dadm from scratch (no simulation)
#   2. call EMC2:::calc_ll directly
#   3. compare against a pure-R TEXG reference that does not use any EMC2 code
#
# ──────────────────────────────────────────────────────────────────────────────


# ── Module-level setup ────────────────────────────────────────────────────────

.lIfun_go <- function(d) factor(rep(2, nrow(d)), levels = 1:2)

.make_design_verify <- function() {
  design(
    model    = SSEXG,
    factors  = list(subjects = 1, S = c("left", "right")),
    Rlevels  = c("left", "right"),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    functions = list(lI = .lIfun_go, SSD = function(d) rep(Inf, nrow(d))),
    formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
                    muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1),
    verbose = FALSE
  )
}

# Fixed parameter values (all in sampled / transformed space)
.p_fixed <- list(
  mu_m   = 0.50,   # match accumulator mean
  mu_mm  = 0.65,   # mismatch accumulator mean
  sig    = 0.05,
  tau    = 0.15,
  muS    = 0.20,
  sigS   = 0.04,
  tauS   = 0.08,
  gf     = 0.08,
  tf     = 0.08,
  lb     = 0.05,   # exg_lb = exgS_lb
  lbS    = 0.05
)

.make_p_vector <- function(design_ss) {
  pv <- sampled_pars(design_ss, doMap = FALSE)
  pv["mu_lMFALSE"] <- log(.p_fixed$mu_mm)
  pv["mu_lMTRUE"]  <- log(.p_fixed$mu_m)
  pv["sigma"]      <- log(.p_fixed$sig)
  pv["tau"]        <- log(.p_fixed$tau)
  pv["muS"]        <- log(.p_fixed$muS)
  pv["sigmaS"]     <- log(.p_fixed$sigS)
  pv["tauS"]       <- log(.p_fixed$tauS)
  pv["gf"]         <- qnorm(.p_fixed$gf)
  pv["tf"]         <- qnorm(.p_fixed$tf)
  pv
}

# Build a single-trial data frame for design_model (no simulation).
# SSD is pre-specified so the Ffunction is skipped.
# R must be a factor with levels c("left","right") even when NA.
.make_trial_dat <- function(S_val, R_val, rt_val, SSD_val) {
  data.frame(
    subjects = factor("1"),
    S        = factor(S_val, levels = c("left", "right")),
    R        = factor(R_val, levels = c("left", "right")),
    rt       = rt_val,
    SSD      = SSD_val,
    trials   = 1L,
    stringsAsFactors = FALSE
  )
}

# Wrapper: call C++ log-likelihood on dadm
.cpp_ll_single <- function(p_vector, dadm) {
  model   <- attr(dadm, "model")()
  p_mat   <- matrix(p_vector, nrow = 1); colnames(p_mat) <- names(p_vector)
  p_types <- names(model$p_types)
  consts  <- attr(dadm, "constants"); if (is.null(consts)) consts <- NA
  designs <- setNames(
    lapply(p_types, function(p) {
      d <- attr(dadm, "designs")[[p]]
      d[attr(d, "expand"), , drop = FALSE]
    }), p_types)
  as.numeric(
    EMC2:::calc_ll(p_mat, dadm,
                   constants     = consts,
                   designs       = designs,
                   type          = model$c_name,
                   bound         = model$bound,
                   transform     = model$transform,
                   pretransforms = model$pre_transform,
                   p_types       = p_types,
                   min_ll        = log(1e-10),
                   trend         = model$trend)
  )
}

# ── Pure-R TEXG reference (no EMC2 code used) ─────────────────────────────────

.dexg_r <- function(x, mu, sigma, tau) {
  # ExGaussian PDF: convolution of Normal(mu, sigma^2) and Exp(1/tau)
  z <- (x - mu - sigma^2 / tau) / sigma
  exp(-log(tau) - (x - mu) / tau + sigma^2 / (2 * tau^2) + pnorm(z, log.p = TRUE))
}

.pexg_r <- function(x, mu, sigma, tau) {
  # ExGaussian CDF
  pnorm((x - mu) / sigma) -
    exp(-(x - mu) / tau + sigma^2 / (2 * tau^2)) *
    pnorm((x - mu - sigma^2 / tau) / sigma)
}

.dtexg_r <- function(x, mu, sigma, tau, lb = 0.05) {
  # Lower-truncated ExGaussian PDF  (truncation at lb)
  .dexg_r(x, mu, sigma, tau) / (1 - .pexg_r(lb, mu, sigma, tau))
}

.ptexg_r <- function(x, mu, sigma, tau, lb = 0.05, lower.tail = TRUE) {
  # Lower-truncated ExGaussian CDF / survival
  nc  <- 1 - .pexg_r(lb, mu, sigma, tau)
  cdf <- pmax(0, pmin(1,
               (.pexg_r(x, mu, sigma, tau) - .pexg_r(lb, mu, sigma, tau)) / nc))
  if (lower.tail) cdf else 1 - cdf
}

# Probability that stop process finishes before BOTH go accumulators.
# Integral: ∫ f_stop(s) · S_go1(s+SSD) · S_go2(s+SSD) ds  from lbS to ub
.pstop_r <- function(SSD,
                     muS, sigS, tauS, lbS,
                     mu1, sig1, tau1, lb1,
                     mu2, sig2, tau2, lb2) {
  ub <- muS + 8 * sigS + 16 * tauS    # same heuristic as C++ (k_sigma=8, k_tau=16)
  integrate(function(s) {
    .dtexg_r(s, muS, sigS, tauS, lbS) *
      .ptexg_r(s + SSD, mu1, sig1, tau1, lb1, lower.tail = FALSE) *
      .ptexg_r(s + SSD, mu2, sig2, tau2, lb2, lower.tail = FALSE)
  }, lower = lbS, upper = ub, rel.tol = 1e-7)$value
}


# ── Tests ─────────────────────────────────────────────────────────────────────

# ── Trial type 1: go trial, observed response ─────────────────────────────────
# LL = log(1-gf) + log(dtexg(rt, mu_match)) + log(Stexg(rt, mu_mismatch))
test_that("SSEXG C++ LL matches formula: go trial with observed response", {
  design_ss <- .make_design_verify()
  p_vector  <- .make_p_vector(design_ss)

  RT <- 0.60   # S="left" so lR="left" is the match (mu=0.50)
  dat  <- .make_trial_dat("left", "left", RT, Inf)
  dadm <- EMC2:::design_model(dat, design_ss, verbose = FALSE)

  p <- .p_fixed
  # log(1-gf) + log f_match(rt) + log S_mismatch(rt)
  expected <- log(1 - p$gf) +
    log(.dtexg_r(RT, p$mu_m,  p$sig, p$tau, p$lb)) +
    log(.ptexg_r(RT, p$mu_mm, p$sig, p$tau, p$lb, lower.tail = FALSE))

  expect_equal(.cpp_ll_single(p_vector, dadm), expected, tolerance = 1e-5)
})


# ── Trial type 2: stop trial, go wins ─────────────────────────────────────────
# LL = log(1-gf) + [log f_match(rt) + log S_mismatch(rt)] +
#      log(tf + (1-tf) · S_stop(rt-SSD))
test_that("SSEXG C++ LL matches formula: stop trial, go wins", {
  design_ss <- .make_design_verify()
  p_vector  <- .make_p_vector(design_ss)

  RT <- 0.65; SSD <- 0.20
  dat  <- .make_trial_dat("left", "left", RT, SSD)
  dadm <- EMC2:::design_model(dat, design_ss, verbose = FALSE)

  p <- .p_fixed
  go_lprob <- log(.dtexg_r(RT, p$mu_m,  p$sig, p$tau, p$lb)) +
              log(.ptexg_r(RT, p$mu_mm, p$sig, p$tau, p$lb, lower.tail = FALSE))
  S_stop   <- .ptexg_r(RT - SSD, p$muS, p$sigS, p$tauS, p$lbS, lower.tail = FALSE)
  expected <- log(1 - p$gf) + go_lprob + log(p$tf + (1 - p$tf) * S_stop)

  expect_equal(.cpp_ll_single(p_vector, dadm), expected, tolerance = 1e-5)
})


# ── Trial type 3: go trial, no response (go failure) ─────────────────────────
# LL = log(gf)
test_that("SSEXG C++ LL matches formula: NR go trial (go failure)", {
  design_ss <- .make_design_verify()
  p_vector  <- .make_p_vector(design_ss)

  dat  <- .make_trial_dat("left", NA, NA, Inf)
  dadm <- EMC2:::design_model(dat, design_ss, verbose = FALSE)

  expected <- log(.p_fixed$gf)
  expect_equal(.cpp_ll_single(p_vector, dadm), expected, tolerance = 1e-5)
})


# ── Trial type 4: stop trial, stop wins (NR) ──────────────────────────────────
# LL = log(gf + (1-gf)·(1-tf)·pStop)
# This is the trial type that was systematically broken by the -ffast-math bug:
#   std::isfinite / std::isinf were constant-folded under -ffinite-math-only,
#   causing stop_signal_presented to be always true AND the wrong GSL integrator
#   to be chosen, making pStop → min_ll and LL → log(gf).
test_that("SSEXG C++ LL matches formula: NR stop trial (stop wins)", {
  design_ss <- .make_design_verify()
  p_vector  <- .make_p_vector(design_ss)

  SSD <- 0.20
  dat  <- .make_trial_dat("left", NA, NA, SSD)
  dadm <- EMC2:::design_model(dat, design_ss, verbose = FALSE)

  p <- .p_fixed
  pstop    <- .pstop_r(SSD,
                       p$muS, p$sigS, p$tauS, p$lbS,
                       p$mu_m,  p$sig, p$tau, p$lb,
                       p$mu_mm, p$sig, p$tau, p$lb)
  expected <- log(p$gf + (1 - p$gf) * (1 - p$tf) * pstop)

  # Verify expected is clearly above log(gf): if pstop is non-trivial the
  # two values must differ by at least 0.1.
  expect_gt(expected - log(p$gf), 0.1)

  expect_equal(.cpp_ll_single(p_vector, dadm), expected, tolerance = 1e-4)
})


# ── Trial type 5: stop trial, NR, second SSD value ────────────────────────────
# Repeat for SSD=0.35 to confirm the integration grid is correct for a
# tighter SSD (stop has less time → pStop is lower).
test_that("SSEXG C++ LL matches formula: NR stop trial with SSD=0.35", {
  design_ss <- .make_design_verify()
  p_vector  <- .make_p_vector(design_ss)

  SSD <- 0.35
  dat  <- .make_trial_dat("left", NA, NA, SSD)
  dadm <- EMC2:::design_model(dat, design_ss, verbose = FALSE)

  p <- .p_fixed
  pstop    <- .pstop_r(SSD,
                       p$muS, p$sigS, p$tauS, p$lbS,
                       p$mu_m,  p$sig, p$tau, p$lb,
                       p$mu_mm, p$sig, p$tau, p$lb)
  expected <- log(p$gf + (1 - p$gf) * (1 - p$tf) * pstop)

  expect_equal(.cpp_ll_single(p_vector, dadm), expected, tolerance = 1e-4)
})
