# Natural-scale race flows (cards converted from Malte Lueken's flows:
# context_scale = "natural", input_scaling, affine, pre from the card), the
# hybrid race (register_nn_model(hybrid = )) and CRDM(), its twin. The fixtures
# (golden/malte, built by make_fixtures.R there) are the NLE project's three
# random-weight cards: their densities model nothing, but the evaluator must
# reproduce the NLE R port (tests/testthat/port/flow_race.R) on them. The
# golden values guard the evaluator against later changes.

nle <- function(f) getFromNamespace(f, "EMC2")
port <- new.env()
sys.source(test_path("port", "flow_race.R"), envir = port)

fixture <- function(nm) test_path("golden", "malte", paste0(nm, ".rds"))
golden <- function(nm) readRDS(test_path("golden", "malte", paste0(nm, "_golden.rds")))
port_card <- function(nm) {
  fl <- readRDS(fixture(nm))
  fl$mlp$layers <- lapply(fl$mlp$layers, function(l) list(W = as.matrix(l$W), b = as.numeric(l$b)))
  fl$input_scaling <- if (!is.null(fl$input_scaling)) lapply(fl$input_scaling, as.numeric)
  fl$affine <- isTRUE(fl$affine)
  fl
}
rdm_pars <- c("v", "s", "B"); crdm_pars <- c("v", "amp", "tau", "s", "B")
mk <- list(
  rdm_final = register_nn_model(fixture("rdm_final"), pars = rdm_pars, twin = RDM, exceptions = c(A = 0)),
  rdm_plain = register_nn_model(fixture("rdm_plain"), pars = rdm_pars, twin = RDM, exceptions = c(A = 0)),
  hybrid = register_nn_model(fixture("crdm_final"), pars = crdm_pars, twin = CRDM, exceptions = c(A = 0),
                             hybrid = "amp"),
  allflow = register_nn_model(fixture("crdm_final"), pars = crdm_pars, twin = CRDM, exceptions = c(A = 0)))

# natural-scale rows for a card's golden table (decision time dt, t0 added)
nat_rows <- function(g, t0 = .2) {
  nat <- cbind(v = g$v, B = g$b, A = 0, t0 = t0, s = g$s)
  if (!is.null(g$amp)) nat <- cbind(nat, amp = g$amp, tau = g$tau)
  nat
}
wald_pdf <- function(t, v, s, b) b / (s * sqrt(2 * pi * t^3)) * exp(-(b - v * t)^2 / (2 * s^2 * t))
wald_cdf <- function(t, v, s, b) pnorm((v * t - b) / (s * sqrt(t))) +
  exp(2 * v * b / s^2 + pnorm(-(v * t + b) / (s * sqrt(t)), log.p = TRUE))

test_that("registration: what is read from the card", {
  for (nm in names(mk)) {
    ml <- mk[[nm]](); reg <- ml$nn
    expect_identical(reg$pre, "t0")
    expect_true(all(reg$transforms == "identity"))
    expect_identical(ml$c_name, "NN")
    expect_equal(unname(ml$bound$minmax[, "A"]), c(0, 0))
    expect_equal(ml$bound$exception[["A"]], 0)
    expect_equal(ml$p_types[["A"]], -Inf)
    expect_equal(reg$lower_s, unname(reg$lower))               # the sampled box is the natural box
  }
  expect_identical(mk$rdm_final()$nn$label, "rdm_final.rds")
  mh <- mk$hybrid()
  expect_identical(mh$nn$hybrid, "amp")
  expect_null(mk$allflow()$nn$hybrid)
  for (ml in list(mh, mk$allflow())) {                 # closed at amp = 0 with or without hybrid
    expect_equal(ml$p_types[["amp"]], -Inf)
    expect_equal(ml$bound$exception[["amp"]], 0)
    expect_equal(unname(ml$bound$minmax[, "amp"]), c(0, 1))
    expect_identical(ml$nn$closed, "amp")
    expect_length(ml$nn$refuse_default, 0)
  }
  # the other edges stay open: v = 0 is the box's edge and RDM's exception, but not its default
  expect_false("v" %in% names(mk$rdm_final()$bound$exception))
  expect_identical(mk$hybrid()$nn$twin, "CRDM")
})

test_that("a card called card.json is labelled by its model field", {
  d <- file.path(tempfile(), "somewhere"); dir.create(d, recursive = TRUE)
  saveRDS(readRDS(fixture("rdm_plain")), file.path(d, "card.rds"))
  m <- register_nn_model(file.path(d, "card.rds"), pars = rdm_pars, twin = RDM, exceptions = c(A = 0))
  expect_identical(m()$nn$label, "rdm_plain")
})

test_that("refusals", {
  f <- fixture("rdm_final")
  expect_error(register_nn_model(f, pars = rdm_pars, twin = RDM), "would ignore|A = 0")
  expect_error(register_nn_model(f, pars = rdm_pars, twin = RDM, exceptions = c(A = .3)), "A = 0")
  expect_error(register_nn_model(f, pars = rdm_pars, twin = RDM, exceptions = c(A = 0), hybrid = "amp"),
               "not an input")
  expect_error(register_nn_model(f, pars = rdm_pars, twin = RDM, exceptions = c(A = 0), pre = "s"),
               "contradicts the card")
  expect_error(register_nn_model(f, pars = rdm_pars, twin = RDM, exceptions = c(A = 0),
                                 pre = function(rt, pars) rt - pars[, "t0"]), "contradicts the card")
  expect_error(register_nn_model(fixture("crdm_final"), pars = crdm_pars, twin = CRDM, exceptions = c(A = 0),
                                 hybrid = "tau"), "starts at 0")
  # a card with both standardisations; an affine card without its constants; a wrong output width
  bad <- function(change) {
    fl <- change(readRDS(f)); p <- tempfile(fileext = ".rds"); saveRDS(fl, p)
    register_nn_model(p, pars = rdm_pars, twin = RDM, exceptions = c(A = 0))
  }
  expect_error(bad(function(fl) { fl$scaler <- list(mean = rep(0, 3), scale = rep(1, 3)); fl }),
               "both a scaler and input_scaling")
  expect_error(bad(function(fl) { fl$affine_scale <- NULL; fl }), "affine_scale")
  expect_error(bad(function(fl) { fl$affine <- FALSE; fl }), "outputs")
  # a design that samples A, or moves it, is refused
  expect_error(suppressMessages(design(factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
                                       matchfun = function(d) d$S == d$lR, model = mk$rdm_final,
                                       formula = list(v ~ lM, B ~ 1, t0 ~ 1, A ~ 1))), "fixed at 0")
})

test_that("dfun/pfun == golden values (NLE R port and float64 reference)", {
  for (nm in c("rdm_final", "rdm_plain", "crdm_final")) {
    g <- golden(nm); m <- mk[[if (nm == "crdm_final") "allflow" else nm]]()
    nat <- nat_rows(g)
    d <- m$dfun(g$dt + .2, nat); p <- m$pfun(g$dt + .2, nat)
    fin <- g$port_log_pdf > -690                       # representable as a density
    expect_gt(sum(fin), 100)
    expect_lt(max(abs(log(d[fin]) - g$port_log_pdf[fin])), 1e-9)
    expect_lt(max(abs(log(d[fin]) - g$log_pdf[fin])), 1e-9)
    expect_true(all(d[!fin] < 1e-299))
    expect_lt(max(abs(p - g$port_cdf)), 1e-10)
    # the upper tail, where jax's logsf is unreliable: log_sf_exact
    up <- g$log_sf_exact < -1e-3 & g$log_sf_exact > -30
    expect_lt(max(abs(log1p(-p[up]) - g$log_sf_exact[up]) * exp(g$log_sf_exact[up])), 1e-10)
  }
})

test_that("dfun/pfun == the NLE R port; 0 for rt <= t0 and outside the box", {
  set.seed(5); n <- 40
  for (nm in c("rdm_final", "rdm_plain")) {
    fl <- port_card(nm); m <- mk[[nm]]()
    for (k in 1:8) {
      t0 <- runif(1, .15, .4)
      nat <- cbind(v = runif(1, .05, 7.9), B = runif(1, .3, 3.4), A = 0, t0 = t0, s = runif(1, .3, 3.4))[rep(1, n), ]
      if (k == 8) nat[, "v"] <- 8.5
      rt <- c(runif(10, t0 - .05, t0), runif(n - 10, t0 + .01, t0 + 2.5))
      ref <- port$flow_eval(fl, unname(nat[1, c("v", "s", "B")]), rt[rt > t0] - t0)
      d <- m$dfun(rt, nat); p <- m$pfun(rt, nat)
      expect_true(all(d[rt <= t0] == 0) && all(p[rt <= t0] == 0))
      if (k == 8) { expect_true(all(d == 0) && all(p == 0)); next }
      fin <- ref$pdf > 1e-290
      expect_lt(max(abs(log(d[rt > t0][fin]) - ref$log_pdf[fin])), 1e-10)
      expect_lt(max(abs(p[rt > t0] - ref$cdf)), 1e-10)
    }
  }
})

test_that("hybrid: Wald where amp is exactly 0, the network elsewhere; all-flow sends every row to the network", {
  set.seed(6); n <- 40; fl <- port_card("crdm_final")
  ctx <- c("v", "amp", "tau", "s", "B")
  for (k in 1:6) {
    t0 <- runif(1, .15, .4)
    nat <- cbind(v = rep(runif(2, .05, 7.9), n / 2), B = runif(1, .3, 2.9), A = 0, t0 = t0, s = runif(1, .3, 2.9),
                 tau = runif(1, .02, .45), amp = rep(c(runif(1, .02, .98), 0), n / 2))
    rt <- t0 + runif(n, .01, 2.5)
    un <- nat[, "amp"] == 0
    net <- function(rows) port$flow_eval(fl, unname(nat[rows[1], ctx]), rt[rows] - t0)
    for (hy in c(TRUE, FALSE)) {
      m <- mk[[if (hy) "hybrid" else "allflow"]]()
      d <- m$dfun(rt, nat); p <- m$pfun(rt, nat)
      r <- net(which(!un))
      expect_lt(max(abs(log(d[!un]) - r$log_pdf)[r$pdf > 1e-290]), 1e-10)
      expect_lt(max(abs(p[!un] - r$cdf)), 1e-10)
      if (hy) {
        w <- wald_pdf(rt[un] - t0, nat[un, "v"], nat[un, "s"], nat[un, "B"])
        expect_lt(max(abs(log(d[un]) - log(w))[w > 1e-290]), 1e-10)
        expect_lt(max(abs(p[un] - wald_cdf(rt[un] - t0, nat[un, "v"], nat[un, "s"], nat[un, "B"]))), 1e-12)
        # the RDM's own Wald (pWald multiplies exp(2 k l) by a small tail
        # probability, which costs it digits; the hybrid's is in log form)
        rd <- RDM()
        expect_lt(max(abs(d[un] - rd$dfun(rt[un], nat[un, ])) / pmax(d[un], 1e-300)), 1e-9)
        expect_lt(max(abs(p[un] - rd$pfun(rt[un], nat[un, ]))), 1e-8)
      } else {
        r0 <- net(which(un))
        expect_lt(max(abs(log(d[un]) - r0$log_pdf)[r0$pdf > 1e-290]), 1e-10)
      }
    }
  }
  # the box on v, s, B applies to the Wald rows
  m <- mk$hybrid()
  nat <- cbind(v = c(2, 8.5, 2, 2), B = c(1, 1, 3.2, 1), A = 0, t0 = .2, s = c(1, 1, 1, .2), tau = .1, amp = 0)
  expect_equal(m$dfun(rep(.8, 4), nat) > 0, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(m$pfun(rep(.8, 4), nat) > 0, c(TRUE, FALSE, FALSE, FALSE))
  # where pigt0's exp(2 k l) would overflow the Wald stays finite
  nat <- cbind(v = 7.9, B = 2.9, A = 0, t0 = .2, s = .26, tau = .1, amp = 0)
  expect_equal(m$pfun(.2 + 2.9 / 7.9, nat), wald_cdf(2.9 / 7.9, 7.9, .26, 2.9), tolerance = 1e-12)
})

DES <- list(factors = list(subjects = 1, S = c("left", "right"), D = c("left", "right")), Rlevels = c("left", "right"),
            matchfun = function(d) d$S == d$lR,
            functions = list(lD = function(d) factor(d$D == d$lR, levels = c(FALSE, TRUE))))
p_c <- c(v = log(2), v_lMTRUE = .5, B = log(1), t0 = log(.3), tau = log(.1), amp_lDTRUE = log(.4))
f_c <- list(v ~ lM, B ~ 1, t0 ~ 1, tau ~ 1, amp ~ 0 + lD)
crdm_design <- function(model) suppressMessages(do.call(design, c(DES, list(
  formula = f_c, constants = c(amp_lDFALSE = -Inf), model = model, report_p_vector = FALSE))))

test_that("the design gives exact zeros, and the compiled likelihood is the R path's", {
  r_path <- function(model) { ml <- model(); ml$c_name <- NULL; function() ml }
  for (nm in c("hybrid", "allflow")) {
    des <- crdm_design(mk[[nm]])
    mp <- mapped_pars(des, p_c)
    lD <- as.logical(as.character(mp$lD))
    expect_true(all(mp$amp[!lD] == 0) && all(abs(mp$amp[lD] - .4) < 1e-12) && all(mp$A == 0))
    set.seed(2); dat <- make_data(p_c, des, n_trials = 40)
    expect_true(all(is.finite(dat$rt)) && min(dat$rt) > .3)
    dadm <- nle("nn_dadm")(dat, des)
    P <- rbind(p_c, p_c + .05, p_c - .05, p_c + c(3, 0, 0, 0, 0, 0))   # the last leaves the box
    twN <- nle("nn_trialwise")(P, dadm, des$model)
    twR <- nle("nn_trialwise")(P, dadm, r_path(des$model))
    expect_lt(max(abs(twN - twR)), 1e-9)
    expect_true(all(twN[, 4] == log(1e-10)))
    expect_gt(mean(twN[, 1:3] > log(1e-10)), .5)
    old <- options(emc.ll_backend = "multithreaded", emc.n_threads = 2)
    expect_identical(nle("nn_trialwise")(P, dadm, des$model), twN)
    options(old)
  }
})

test_that("CRDM(): one Volterra solve per distinct row, Wald at amp = 0", {
  tw <- CRDM(); set.seed(7); U <- 5; n <- 12
  th <- cbind(v = runif(U, .5, 6), B = runif(U, .4, 2.5), A = 0, t0 = runif(U, .15, .4), s = runif(U, .4, 2.5),
              amp = c(runif(U - 2, .05, .95), 0, 0), tau = runif(U, .02, .45))
  row <- sample(rep(seq_len(U), n))
  pars <- cbind(th[row, ], b = th[row, "B"]); rt <- pars[, "t0"] + runif(U * n, -.03, .9)
  d <- tw$dfun(rt, pars); p <- tw$pfun(rt, pars)
  rd <- rp <- numeric(length(rt))
  for (u in seq_len(U)) {
    i <- which(row == u & rt > th[u, "t0"]); t <- rt[i] - th[u, "t0"]
    if (th[u, "amp"] == 0) {
      rd[i] <- wald_pdf(t, th[u, "v"], th[u, "s"], th[u, "B"]); rp[i] <- wald_cdf(t, th[u, "v"], th[u, "s"], th[u, "B"])
    } else {
      o <- nle("dcrdm_volterra")(t, th[u, "v"], th[u, "amp"], th[u, "tau"], th[u, "s"], th[u, "B"], 5e-4, 1)
      rd[i] <- o$pdf; rp[i] <- o$cdf
    }
  }
  expect_lt(max(abs(d - rd)) / max(rd), 1e-12)
  expect_lt(max(abs(p - rp)), 1e-12)
  expect_true(all(d[rt <= pars[, "t0"]] == 0) && all(p[rt <= pars[, "t0"]] == 0))
  # amp = 0 is the RDM; rows the model cannot represent give 0, not an error
  un <- pars[, "amp"] == 0
  expect_lt(max(abs(d[un] - RDM()$dfun(rt[un], pars[un, ]))), 1e-9)
  bad <- pars[1:3, ]; bad[1, "tau"] <- -1; bad[2, "v"] <- NA; bad[3, "b"] <- 0; bad[, "amp"] <- .3
  expect_identical(tw$dfun(bad[, "t0"] + .5, bad), c(0, 0, 0))
  # a coarser grid is close; the solver is tested at amp > 0 (at amp = 0 it is the Wald at any dt)
  pu <- which(!un & rt > pars[, "t0"])
  expect_lt(max(abs(CRDM(dt = .002)$dfun(rt[pu], pars[pu, ]) - d[pu])) / max(d), .02)
  expect_gt(max(abs(CRDM(dt = .002)$dfun(rt[pu], pars[pu, ]) - d[pu])), 0)
  # the pulse matters
  flat <- pars[pu, ]; flat[, "amp"] <- 0
  expect_gt(max(abs(tw$dfun(rt[pu], flat) - d[pu])), .01)
})

test_that("CRDM(): simulator", {
  des <- crdm_design(CRDM)
  set.seed(4); dat <- make_data(p_c, des, n_trials = 300)
  set.seed(4); dat2 <- make_data(p_c, des, n_trials = 300)
  expect_identical(dat, dat2)
  expect_true(all(is.finite(dat$rt)) && min(dat$rt) > .3)
  # P(pulsed accumulator wins) by quadrature on the Volterra grid
  vm <- exp(p_c[["v"]] + p_c[["v_lMTRUE"]]); vx <- exp(p_c[["v"]])
  for (cong in c(TRUE, FALSE)) {
    vp <- if (cong) vm else vx; vw <- if (cong) vx else vm
    g <- nle("crdm_volterra_grid")(vp, .4, .1, 1, 1, 5e-4, 8)
    pr <- sum(g$pdf * (1 - wald_cdf(g$t, vw, 1, 1))) * 5e-4
    sel <- (dat$S == dat$D) == cong
    won <- as.character(dat$R[sel]) == as.character(dat$D[sel])
    expect_lt(abs(mean(won) - pr), 4 * sqrt(pr * (1 - pr) / sum(sel)))
  }
  # every accumulator pulsed: simulated to t_max; no finisher by then is an error
  lR <- factor(c("a", "b"))
  pars <- cbind(v = c(2, 1), B = 1, A = 0, t0 = .2, s = 1, amp = .3, tau = .1, b = 1)[rep(1:2, 50), ]
  set.seed(1); o <- nle("rCRDM")(factor(rep(c("a", "b"), 50)), pars)
  expect_true(all(is.finite(o$rt)) && mean(o$R == "a") > .5)
  expect_error(nle("rCRDM")(factor(rep(c("a", "b"), 50)), pars, t_max = .05), "no finisher")
  pars[, "A"] <- .1
  expect_error(nle("rCRDM")(factor(rep(c("a", "b"), 50)), pars), "A must be exactly 0")
})

test_that("nn_cell(functions = ): network and control designs; -Inf constants; effects", {
  sdv <- c(v = .15, v_lMTRUE = .1, B = .15, t0 = .1, tau = .2, amp_lDTRUE = .15)
  cell <- do.call(nn_cell, c(list(mk$hybrid, mean = p_c, sd = sdv, formula = f_c, constants = c(amp_lDFALSE = -Inf)),
                             DES[c("Rlevels", "matchfun", "functions")], list(factors = DES$factors[-1])))
  expect_false(is.null(cell$control_design))
  expect_identical(cell$control_label, "CRDM")
  expect_true(all(cell$box$inside))
  amp <- cell$box[cell$box$parameter == "amp", ]
  expect_equal(sort(amp$prior_lower)[1], 0)             # the unpulsed cells: exactly 0, inside
  # an effect is held to +-k sd of the value it makes
  v <- cell$box[cell$box$parameter == "v", ]
  expect_equal(max(v$prior_upper), exp(log(2) + .5 + 4 * sqrt(.15^2 + .1^2)))
  set.seed(1); th <- nle("nn_cell_draws")(cell, 200)
  expect_true(all(abs(th[, "v"] + th[, "v_lMTRUE"] - log(2) - .5) <= 4 * sqrt(.15^2 + .1^2)))
  # a prior that leaves the region is still refused
  sdv[["v"]] <- .3
  expect_error(do.call(nn_cell, c(list(mk$hybrid, mean = p_c, sd = sdv, formula = f_c,
                                       constants = c(amp_lDFALSE = -Inf)),
                                  DES[c("Rlevels", "matchfun", "functions")], list(factors = DES$factors[-1]))),
               "leaves the network's training region")
})
