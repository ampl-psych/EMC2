# Tests for the Gauss-Legendre stop-success integral variants (EXG + RDEX).
# Verifies (Q4) GL agrees with the adaptive integrate() route, and (Q5) that
# accuracy is tunable via n_nodes.

# ---- helpers to build a one-accumulator parameter matrix -------------------

make_texg_pars <- function(muG, sigG, tauG, lbG, muS, sigS, tauS, lbS) {
  # column layout expected by ss_texg_stop_success_*: muG=0 sigG=1 tauG=2
  # muS=3 sigS=4 tauS=5 (6,7 unused) lbG=8 lbS=9   (0-based -> R 1-based below)
  matrix(c(muG, sigG, tauG, muS, sigS, tauS, 0, 0, lbG, lbS),
         nrow = 1, ncol = 10)
}

make_rdex_pars <- function(v, B, A, t0, s, muS, sigS, tauS, lbS) {
  # go: v=0 B=1 A=2 t0=3 s=4 ; stop: muS=5 sigS=6 tauS=7 (8,9 unused) lbS=10
  matrix(c(v, B, A, t0, s, muS, sigS, tauS, 0, 0, lbS),
         nrow = 1, ncol = 11)
}

# ---- Q4: GL agrees with adaptive integrate() -------------------------------

test_that("EXG stop-success: GL matches adaptive integrate", {
  p   <- make_texg_pars(muG = 0.50, sigG = 0.05, tauG = 0.08, lbG = 0,
                        muS = 0.20, sigS = 0.03, tauS = 0.05, lbS = 0)
  SSD <- 0.15
  adaptive <- ss_texg_stop_success_value(SSD, p, method = "integrate")
  gl       <- ss_texg_stop_success_value(SSD, p, method = "gl", n_nodes = 128L)
  expect_gt(adaptive, 0)
  expect_equal(gl, adaptive, tolerance = 1e-5)
})

test_that("RDEX stop-success: GL matches adaptive integrate", {
  p   <- make_rdex_pars(v = 2.5, B = 0.6, A = 0.2, t0 = 0.15, s = 1.0,
                        muS = 0.20, sigS = 0.03, tauS = 0.05, lbS = 0)
  SSD <- 0.15
  adaptive <- ss_rdex_stop_success_value(SSD, p, method = "integrate")
  gl       <- ss_rdex_stop_success_value(SSD, p, method = "gl", n_nodes = 128L)
  expect_gt(adaptive, 0)
  expect_equal(gl, adaptive, tolerance = 1e-5)
})

test_that("R route GL matches R route integrate (EXG)", {
  muG<-0.50; sigG<-0.05; tauG<-0.08; lbG<-0
  muS<-0.20; sigS<-0.03; tauS<-0.05; lbS<-0; SSD<-0.15
  mu<-c(muS,muG); sigma<-c(sigS,sigG); tau<-c(tauS,tauG); lb<-c(lbS,lbG)
  a <- stop_success_texg_R(mu,sigma,tau,lb,SSD, method="integrate")
  g <- stop_success_texg_R(mu,sigma,tau,lb,SSD, method="gl", n_nodes=128L)
  expect_gt(a, 0)
  expect_equal(g, a, tolerance = 1e-5)
})

test_that("C++ and R routes agree (EXG)", {
  muG<-0.50; sigG<-0.05; tauG<-0.08; lbG<-0
  muS<-0.20; sigS<-0.03; tauS<-0.05; lbS<-0; SSD<-0.15
  p  <- make_texg_pars(muG,sigG,tauG,lbG,muS,sigS,tauS,lbS)
  cpp <- ss_texg_stop_success_value(SSD, p, method = "gl", n_nodes = 128L)
  r   <- stop_success_texg_R(c(muS,muG), c(sigS,sigG), c(tauS,tauG),
                             c(lbS,lbG), SSD, method = "gl", n_nodes = 128L)
  expect_equal(r, cpp, tolerance = 1e-6)
})

# ---- Q5: accuracy is tunable via n_nodes -----------------------------------

test_that("GL converges to the adaptive answer as n_nodes increases", {
  p   <- make_texg_pars(muG = 0.50, sigG = 0.05, tauG = 0.08, lbG = 0,
                        muS = 0.20, sigS = 0.03, tauS = 0.05, lbS = 0)
  SSD <- 0.15
  truth <- ss_texg_stop_success_value(SSD, p, method = "integrate")
  errs  <- sapply(c(4L, 8L, 16L, 32L, 64L, 128L), function(n)
    abs(ss_texg_stop_success_value(SSD, p, method = "gl", n_nodes = n) - truth))
  # GL error need not fall strictly at every step (it can wiggle near machine
  # precision); the meaningful claims are: tiny at the top end, and the
  # high-node error is far smaller than the low-node error.
  expect_true(errs[length(errs)] < 1e-6)
  expect_lt(errs[length(errs)], errs[1])
  expect_lt(min(errs[4:6]), errs[1] / 100)   # 32+ nodes >= 2 orders better than 4
})

# ---- stop_method config plumbing (the live C++ dispatcher) -----------------

test_that("live dispatcher honours emc2_set_stop_method, incl. n_nodes", {
  on.exit(emc2_set_stop_method("auto", 64L), add = TRUE)   # restore default
  SSD <- 0.15
  p <- make_texg_pars(0.50, 0.05, 0.08, 0, 0.20, 0.02, 0.04, 0)

  emc2_set_stop_method("integrate")
  expect_identical(emc2_get_stop_method()$method, "integrate")
  expect_identical(ss_texg_stop_success_value(SSD, p, method = "live"),
                   ss_texg_stop_success_value(SSD, p, method = "integrate"))

  emc2_set_stop_method("gl", 128L)
  expect_identical(ss_texg_stop_success_value(SSD, p, method = "live"),
                   ss_texg_stop_success_value(SSD, p, method = "gl",
                                              n_nodes = 128L))

  # stop_n_nodes is honoured: a 4-node rule is visibly different, 128 is not
  emc2_set_stop_method("gl", 4L)
  gl4 <- ss_texg_stop_success_value(SSD, p, method = "live")
  emc2_set_stop_method("gl", 128L)
  gl128 <- ss_texg_stop_success_value(SSD, p, method = "live")
  ref <- ss_texg_stop_success_value(SSD, p, method = "integrate")
  expect_gt(abs(gl4 - ref), abs(gl128 - ref))
  expect_equal(gl128, ref, tolerance = 1e-5)

  emc2_set_stop_method("auto", 64L)
  expect_identical(ss_texg_stop_success_value(SSD, p, method = "live"),
                   ss_texg_stop_success_value(SSD, p, method = "auto"))

  # RDEX live path reads the same config
  pr <- make_rdex_pars(v = 2.5, B = 0.6, A = 0.2, t0 = 0.15, s = 1.0,
                       muS = 0.20, sigS = 0.03, tauS = 0.05, lbS = 0)
  emc2_set_stop_method("integrate")
  expect_identical(ss_rdex_stop_success_value(SSD, pr, method = "live"),
                   ss_rdex_stop_success_value(SSD, pr, method = "integrate"))
  emc2_set_stop_method("gl", 96L)
  expect_identical(ss_rdex_stop_success_value(SSD, pr, method = "live"),
                   ss_rdex_stop_success_value(SSD, pr, method = "gl",
                                              n_nodes = 96L))
})

# ---- "auto"/"analytic" dispatch (n_go == 1 closed form + GL fallback) -------

test_that("auto matches integrate across a parameter grid (n_go = 1)", {
  grid <- expand.grid(sigma = c(.05, .03, .02, .01, .005, .002),
                      tau   = c(.04, .08, .12),
                      SSD   = c(0, .15, .3))
  for (i in seq_len(nrow(grid))) {
    p <- make_texg_pars(muG = 0.50, sigG = 0.05, tauG = 0.08, lbG = 0.05,
                        muS = 0.20, sigS = grid$sigma[i], tauS = grid$tau[i],
                        lbS = 0.05)
    ref  <- ss_texg_stop_success_value(grid$SSD[i], p, method = "integrate")
    auto <- ss_texg_stop_success_value(grid$SSD[i], p, method = "auto")
    expect_equal(auto, ref, tolerance = 1e-4,
                 info = sprintf("sigma=%.3f tau=%.3f SSD=%.2f",
                                grid$sigma[i], grid$tau[i], grid$SSD[i]))
  }
})

test_that("n_go = 1 truncated analytic form is exact (vs tight integrate)", {
  p   <- make_texg_pars(muG = 0.50, sigG = 0.05, tauG = 0.08, lbG = 0.05,
                        muS = 0.20, sigS = 0.03, tauS = 0.05, lbS = 0.05)
  for (SSD in c(0, .1, .2, .3, .5)) {
    expect_identical(ss_texg_stop_success_auto_branch(SSD, p), "analytic_trunc")
    ana <- ss_texg_stop_success_value(SSD, p, method = "analytic")
    # tight adaptive reference (small tolerances, wide window)
    ref <- ss_texg_stop_success_value(SSD, p, method = "integrate",
                                      max_subdiv = 200, abs_tol = 1e-12,
                                      rel_tol = 1e-10, k_sigma = 10,
                                      k_tau = 20)
    expect_equal(ana, ref, tolerance = 1e-7, info = sprintf("SSD=%.2f", SSD))
  }
})

test_that("n_go = 1 full-line analytic form matches R reference (lb = -Inf)", {
  muG <- .5; sigG <- .05; tauG <- .08; muS <- .2; sigS <- .03; tauS <- .05
  SSD <- .15
  p <- make_texg_pars(muG, sigG, tauG, -Inf, muS, sigS, tauS, -Inf)
  expect_identical(ss_texg_stop_success_auto_branch(SSD, p),
                   "analytic_fullline")
  ana <- ss_texg_stop_success_value(SSD, p, method = "analytic")
  # exact-exG R reference over the full line
  dexg_ref <- function(x, mu, sigma, tau)
    exp(-log(tau) + (mu - x)/tau + sigma^2/(2*tau^2) +
          pnorm((x - mu)/sigma - sigma/tau, log.p = TRUE))
  sexg_ref <- function(t, mu, sigma, tau)
    pnorm(-(t - mu)/sigma) +
      exp((mu - t)/tau + sigma^2/(2*tau^2) +
            pnorm((t - mu)/sigma - sigma/tau, log.p = TRUE))
  ref <- integrate(function(x) dexg_ref(x, muS, sigS, tauS) *
                     sexg_ref(x + SSD, muG, sigG, tauG),
                   lower = -Inf, upper = Inf,
                   rel.tol = 1e-12, abs.tol = 1e-14)$value
  expect_equal(ana, ref, tolerance = 1e-10)
})

test_that("R and C++ analytic routes agree to near machine precision", {
  for (SSD in c(0, .15, .4)) {
    p <- make_texg_pars(0.50, 0.05, 0.08, 0.05, 0.20, 0.03, 0.05, 0.05)
    cpp <- ss_texg_stop_success_value(SSD, p, method = "analytic")
    r   <- stop_success_texg_analytic1_R(mu = c(.2, .5), sigma = c(.03, .05),
                                         tau = c(.05, .08), lb = c(.05, .05),
                                         SSD = SSD)
    expect_equal(r, cpp, tolerance = 1e-12)
  }
})

test_that("auto guards fall back to GL and still match integrate", {
  SSD <- 0.15
  # guard 1: sigS > 10*tauS (ANALYTIC_MAX_SIGS_OVER_TAUS; small-tauS region
  # where logC0 cancellation/dexg tail erode the closed form -> GL fallback)
  p_guard <- make_texg_pars(0.50, 0.05, 0.08, 0.05, 0.20, 0.10, 0.005, 0.05)
  expect_identical(ss_texg_stop_success_auto_branch(SSD, p_guard), "gl_guard")
  expect_equal(ss_texg_stop_success_value(SSD, p_guard, method = "auto"),
               ss_texg_stop_success_value(SSD, p_guard, method = "integrate"),
               tolerance = 1e-4)
  # guard 2: kinked domain (lb_go > lbS + SSD)
  p_kink <- make_texg_pars(0.50, 0.05, 0.08, 0.30, 0.20, 0.03, 0.05, 0.05)
  expect_identical(ss_texg_stop_success_auto_branch(SSD, p_kink), "gl_guard")
  expect_equal(ss_texg_stop_success_value(SSD, p_kink, method = "auto"),
               ss_texg_stop_success_value(SSD, p_kink, method = "integrate"),
               tolerance = 1e-4)
  # n_go >= 2 and finite upper limits always take GL
  p1 <- make_texg_pars(0.50, 0.05, 0.08, 0.05, 0.20, 0.03, 0.05, 0.05)
  p2 <- rbind(p1, make_texg_pars(0.45, 0.06, 0.12, 0.05,
                                 0.20, 0.03, 0.05, 0.05))
  expect_identical(ss_texg_stop_success_auto_branch(SSD, p2), "gl_ngo2plus")
  expect_identical(ss_texg_stop_success_auto_branch(SSD, p1, upper = 0.4),
                   "gl_finite_upper")
  expect_equal(ss_texg_stop_success_value(SSD, p2, method = "auto"),
               ss_texg_stop_success_value(SSD, p2, method = "integrate"),
               tolerance = 1e-4)
  expect_equal(ss_texg_stop_success_value(SSD, p1, method = "auto", upper = 0.4),
               ss_texg_stop_success_value(SSD, p1, method = "integrate",
                                          upper = 0.4),
               tolerance = 1e-4)
})

test_that("R route auto/analytic matches R route integrate", {
  muG <- 0.50; sigG <- 0.05; tauG <- 0.08; lbG <- 0.05
  muS <- 0.20; sigS <- 0.03; tauS <- 0.05; lbS <- 0.05; SSD <- 0.15
  mu <- c(muS, muG); sigma <- c(sigS, sigG); tau <- c(tauS, tauG)
  lb <- c(lbS, lbG)
  ref  <- stop_success_texg_R(mu, sigma, tau, lb, SSD, method = "integrate")
  auto <- stop_success_texg_R(mu, sigma, tau, lb, SSD, method = "auto")
  ana  <- stop_success_texg_R(mu, sigma, tau, lb, SSD, method = "analytic")
  expect_equal(auto, ref, tolerance = 1e-6)
  expect_identical(auto, ana)   # both take the analytic branch here
  # guard-tripping params (sigS = 0.10 > 10*tauS = 0.05) fall back to
  # (density-bumped) GL in R too
  auto_g <- stop_success_texg_R(c(.2, .5), c(.10, .05), c(.005, .08), lb, SSD,
                                method = "auto")
  ref_g  <- stop_success_texg_R(c(.2, .5), c(.10, .05), c(.005, .08), lb, SSD,
                                method = "integrate")
  expect_equal(auto_g, ref_g, tolerance = 1e-4)
})

test_that("pbvn_R bivariate normal CDF is correct", {
  # independence
  expect_equal(pbvn_R(0.3, -0.7, 0), pnorm(0.3) * pnorm(-0.7),
               tolerance = 1e-14)
  # perfect correlation: P(X<=h, X<=k) = pnorm(min(h, k))
  expect_equal(pbvn_R(0.5, 1.2, 1), pnorm(0.5), tolerance = 1e-14)
  # antithetic: P(X<=h, -X<=k) = max(0, pnorm(h) - pnorm(-k))
  expect_equal(pbvn_R(0.5, 1.2, -1), pnorm(0.5) - pnorm(-1.2),
               tolerance = 1e-14)
  # symmetry in arguments
  expect_equal(pbvn_R(0.4, -1.1, 0.6), pbvn_R(-1.1, 0.4, 0.6),
               tolerance = 1e-14)
  # against a 2D quadrature reference at moderate and high correlation
  for (rho in c(-0.85, -0.4, 0.2, 0.6, 0.95)) {
    h <- 0.3; k <- -0.5
    ref <- integrate(function(x) sapply(x, function(xi)
      dnorm(xi) * pnorm((k - rho * xi) / sqrt(1 - rho^2))),
      lower = -Inf, upper = h, rel.tol = 1e-12)$value
    expect_equal(pbvn_R(h, k, rho), ref, tolerance = 1e-10,
                 info = sprintf("rho=%.2f", rho))
  }
})

# ---- end-to-end: stop_method threads from SSEXG() through the likelihood ----

test_that("SSEXG(stop_method) is honoured end-to-end via calc_ll_manager", {
  skip_if_not_installed("EMC2")
  on.exit(emc2_set_stop_method("auto", 64L), add = TRUE)

  make_dadm <- function(stop_method, stop_n_nodes = 64L) {
    des <- design(
      model    = SSEXG(stop_method = stop_method, stop_n_nodes = stop_n_nodes),
      factors  = list(subjects = 1, S = c("left", "right")),
      Rlevels  = c("left", "right"),
      matchfun = function(d) as.character(d$S) == as.character(d$lR),
      functions = list(
        lI  = function(d) factor(rep(2, nrow(d)), levels = 1:2),
        SSD = function(d) rep(Inf, nrow(d))),
      formula  = list(mu ~ 0 + lM, sigma ~ 1, tau ~ 1,
                      muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1),
      verbose  = FALSE
    )
    # two stop trials (one successful stop, one failed) + one go trial
    dat <- data.frame(
      subjects = factor(rep("1", 3)),
      S        = factor(c("left", "left", "right"),
                        levels = c("left", "right")),
      R        = factor(c(NA, "left", "right"),
                        levels = c("left", "right")),
      rt       = c(NA, 0.55, 0.61),
      SSD      = c(0.2, 0.25, Inf),
      trials   = 1:3
    )
    list(des = des, dadm = design_model(dat, des, verbose = FALSE))
  }

  pv_for <- function(des) {
    pv <- sampled_pars(des, doMap = FALSE)
    pv["mu_lMFALSE"] <- log(0.65); pv["mu_lMTRUE"] <- log(0.50)
    pv["sigma"] <- log(0.05); pv["tau"] <- log(0.15)
    pv["muS"] <- log(0.20); pv["sigmaS"] <- log(0.04); pv["tauS"] <- log(0.08)
    pv["gf"] <- qnorm(0.05); pv["tf"] <- qnorm(0.05)
    pv
  }

  ll_for <- function(stop_method, stop_n_nodes = 64L) {
    x <- make_dadm(stop_method, stop_n_nodes)
    pv <- pv_for(x$des)
    proposals <- matrix(pv, nrow = 2, ncol = length(pv), byrow = TRUE)
    colnames(proposals) <- names(pv)
    model_fun <- attr(x$dadm, "model")
    lls <- calc_ll_manager(proposals, x$dadm, model_fun)
    # the model list carries the user's choice
    expect_identical(model_fun()$stop_method, stop_method)
    lls[1]
  }

  ll_int  <- ll_for("integrate")
  expect_identical(emc2_get_stop_method()$method, "integrate")
  ll_gl   <- ll_for("gl")
  expect_identical(emc2_get_stop_method()$method, "gl")
  ll_auto <- ll_for("auto")
  expect_identical(emc2_get_stop_method()$method, "auto")
  ll_gl4  <- ll_for("gl", stop_n_nodes = 4L)
  expect_identical(emc2_get_stop_method()$n_nodes, 4L)

  expect_equal(ll_gl, ll_int, tolerance = 1e-5)
  expect_equal(ll_auto, ll_int, tolerance = 1e-5)
  # a 4-node rule is measurably worse than the 64-node default
  expect_gt(abs(ll_gl4 - ll_int), abs(ll_gl - ll_int))

  # R likelihood route honours the same choice (sfun closure)
  x <- make_dadm("gl", 128L)
  pv <- pv_for(x$des)
  r_gl  <- calc_ll_R(pv, attr(x$dadm, "model")(), x$dadm)
  x2 <- make_dadm("integrate")
  r_int <- calc_ll_R(pv, attr(x2$dadm, "model")(), x2$dadm)
  expect_equal(r_gl, r_int, tolerance = 1e-5)
})

test_that("tight stop density needs more nodes (low n under-resolves)", {
  # very peaked stop process: small sigS and tauS
  p   <- make_texg_pars(muG = 0.50, sigG = 0.05, tauG = 0.08, lbG = 0,
                        muS = 0.20, sigS = 0.005, tauS = 0.008, lbS = 0)
  SSD <- 0.15
  truth   <- ss_texg_stop_success_value(SSD, p, method = "integrate")
  err_low  <- abs(ss_texg_stop_success_value(SSD, p, method="gl", n_nodes=8L)  - truth)
  err_high <- abs(ss_texg_stop_success_value(SSD, p, method="gl", n_nodes=256L) - truth)
  expect_gt(err_low, err_high)          # more nodes -> closer to truth
  expect_lt(err_high, 1e-4)
})

test_that("gl_auto_nodes_R: density floor, quantised, capped", {
  # narrow window: no bump (returns the floor)
  expect_identical(gl_auto_nodes_R(64L, 0.10, 0.20, 0.05), 64L)
  # wide window vs a sharp peak: bump up, quantised to a multiple of 32
  n <- gl_auto_nodes_R(64L, 0.05, 0.55, 0.025)   # (ub-lo)/sigS = 20 -> want 120
  expect_identical(n, 128L)
  expect_identical(n %% 32L, 0L)
  # never below the requested floor; capped at 256
  expect_gte(gl_auto_nodes_R(96L, 0.10, 0.20, 0.05), 96L)
  expect_identical(gl_auto_nodes_R(64L, 0.0, 10.0, 0.01), 256L)  # want huge
  # shrinking the window past the floor leaves the floor untouched
  expect_identical(gl_auto_nodes_R(64L, 0.10, 0.10, 0.05), 64L)
})

test_that("auto stays accurate for a small-tauS sharp stop peak (n_go=1)", {
  # Colleague's regime: tauS small relative to sigS. With sigS <= 10*tauS the
  # exact closed form is used; push tauS smaller so the guard hands off to the
  # density-bumped GL route, and check both still match adaptive integrate.
  SSD <- 0.20
  # (a) analytic region: sigS/tauS = 1.9 (true colleague values)
  p_ana <- make_texg_pars(muG = 0.436, sigG = 0.044, tauG = 0.070, lbG = 0.05,
                          muS = 0.144, sigS = 0.025, tauS = 0.013, lbS = 0.05)
  expect_identical(ss_texg_stop_success_auto_branch(SSD, p_ana), "analytic_trunc")
  expect_equal(ss_texg_stop_success_value(SSD, p_ana, method = "auto"),
               ss_texg_stop_success_value(SSD, p_ana, method = "integrate"),
               tolerance = 1e-6)
  # (b) guard handoff: sigS/tauS = 12.5 -> GL fallback, density bump keeps it tight
  p_gl <- make_texg_pars(muG = 0.436, sigG = 0.044, tauG = 0.070, lbG = 0.05,
                         muS = 0.144, sigS = 0.025, tauS = 0.002, lbS = 0.05)
  expect_identical(ss_texg_stop_success_auto_branch(SSD, p_gl), "gl_guard")
  expect_equal(ss_texg_stop_success_value(SSD, p_gl, method = "auto"),
               ss_texg_stop_success_value(SSD, p_gl, method = "integrate"),
               tolerance = 1e-4)
})
