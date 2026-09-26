# Zero-box DDM flows (register_nn_model(): the "zerobox" transform, exact
# zeros, bound exceptions and the st0 = 0 support rule). The fixture is the
# shipped cap256w_c4 relabelled as zero-box, exactly as the NLE project's
# handover/zerobox/make_zerobox_fixture.R builds ~/R/nle-cards/
# ddm_zerobox_fixture.json: its numbers model nothing, but the evaluator must
# reproduce the NLE R port (tests/testthat/port, copied from the NLE project's
# port/R at the zero-box handover) on it, including at exact zeros. This
# guards exact zeros against later evaluator changes.

port <- new.env()
sys.source(test_path("port", "flow_race.R"), envir = port)
sys.source(test_path("port", "flow_ddm.R"), envir = port)

ZC <- list(sv = 0.01, SZ = 0.05, st0 = 0.002)
zb_relabel <- function(fl) {
  fl$context_encoding <- "zerobox"; fl$zerobox_c <- ZC
  fl$context_transforms[names(ZC)] <- "zerobox"
  cn <- fl$context_names; lo <- fl$bounds_sampled$lower; hi <- fl$bounds_sampled$upper
  for (p in names(ZC)) {
    j <- match(p, cn); lo[j] <- 0
    hi[j] <- asinh(c(sv = 4, SZ = 0.99, st0 = 0.49)[[p]] / ZC[[p]])
  }
  fl$bounds_sampled <- list(lower = lo, upper = hi)
  bn <- fl$bounds_natural
  for (p in names(ZC)) bn$lower[match(p, names(fl$context_transforms))] <- 0
  fl$bounds_natural <- bn
  fl
}
zb_bundle <- function(n_members = 1) {
  b <- readRDS(system.file("extdata", "flownn", "ddm_cap256w_c4.rds", package = "EMC2"))
  b$members <- rep(lapply(b$members, zb_relabel), n_members)
  f <- tempfile(fileext = ".rds"); saveRDS(b, f); f
}
zb_file <- zb_bundle()
m <- register_nn_model(zb_file, twin = DDM)
fl <- zb_relabel(readRDS(system.file("extdata", "flownn", "ddm_cap256w_c4.rds", package = "EMC2"))$members[[1]])

nat_of <- function(th, n) cbind(v = th[["v"]], a = exp(th[["a"]]), t0 = exp(th[["t0"]]), s = exp(th[["s"]]),
                               Z = pnorm(th[["Z"]]), SZ = pnorm(th[["SZ"]]), sv = exp(th[["sv"]]),
                               st0 = exp(th[["st0"]]))[rep(1, n), , drop = FALSE]
regimes <- list(pos = character(0), sv0 = "sv", SZ0 = "SZ", st00 = "st0", all0 = c("sv", "SZ", "st0"))

test_that("registration: zerobox transform, -Inf defaults, bound exceptions at 0", {
  ml <- m(); reg <- ml$nn
  zb <- c("SZ", "sv", "st0")
  expect_identical(unname(reg$transform_codes[match(zb, reg$pars)]), rep(3L, 3))
  expect_equal(reg$zerobox_c[match(zb, reg$pars)], unname(unlist(ZC[zb])))
  expect_true(all(is.infinite(ml$p_types[zb]) & ml$p_types[zb] < 0))
  expect_equal(unname(ml$bound$exception[zb]), c(0, 0, 0))
  # the analytic DDM floor, not a new one; the top is the training region's
  expect_equal(unname(ml$bound$minmax[1, zb]), unname(DDM()$bound$minmax[1, zb]))
  expect_equal(unname(ml$bound$minmax[2, "sv"]), 4)
  expect_length(intersect(reg$refuse_default, zb), 0)   # default 0 is inside the (closed) region
})

test_that("dfun/pfun == NLE R port, all positive and at exact zeros (golden values too)", {
  th0 <- c(v = 1.2, a = log(1.1), t0 = log(.3), s = log(1), Z = qnorm(.45), SZ = qnorm(.2), sv = log(.3), st0 = log(.08))
  rt <- c(.25, .3, .31, .45, .7, 1.2); R <- c(1L, 2L, 1L, 2L, 1L, 2L)
  # recorded from the NLE port (port/R/flow_ddm.R) on ~/R/nle-cards/ddm_zerobox_fixture.json, 2026-09-26
  gold <- list(
    pos  = c(-15.082672076346, -7.531474916538, -15.258188426895, -6.545539080203, -13.857676361102, -7.170337706493),
    sv0  = c(-21.829079548792, -16.329119350432, -22.074081501597, -16.714080216525, -22.884879049757, -16.948714325035),
    SZ0  = c(-15.108373547922, -9.541257233638, -15.296968716826, -7.954731251435, -13.650932347141, -8.933519773459),
    st00 = c(-Inf, -Inf, 0.102988274161, -2.118695790422, -0.330146544808, -2.004858371291),
    all0 = c(-Inf, -Inf, -8.487090847495, -1.143323166980, -1.207261400418, -0.236310911134))
  for (g in names(regimes)) {
    th <- th0; th[regimes[[g]]] <- -Inf
    nat <- nat_of(th, length(rt))
    ref <- port$ddm_eval(fl, th, rt, R)
    lp <- log(m()$dfun(rt, R, nat)); cdf <- m()$pfun(rt, R, nat)
    expect_identical(is.finite(lp), is.finite(ref$log_pdf), label = g)
    expect_identical(is.finite(lp), is.finite(gold[[g]]), label = g)
    fin <- is.finite(ref$log_pdf)
    expect_lt(max(abs(lp[fin] - ref$log_pdf[fin])), 1e-7)
    expect_lt(max(abs(lp[fin] - gold[[g]][fin])), 1e-7)
    expect_lt(max(abs(cdf - ref$cdf)), 1e-7)
  }
})

test_that("random theta in every regime match the port; st0 = 0 support rule on both sides of t0", {
  set.seed(7)
  for (g in names(regimes)) for (k in 1:8) {
    th <- c(v = rnorm(1, 1, .5), a = log(runif(1, .5, 2)), t0 = log(runif(1, .15, .4)), s = log(runif(1, .5, 1.5)),
            Z = qnorm(runif(1, .3, .7)), SZ = qnorm(runif(1, .05, .6)), sv = log(runif(1, .05, 2)),
            st0 = log(runif(1, .01, .3)))
    th[regimes[[g]]] <- -Inf
    t0 <- exp(th[["t0"]])
    rt <- c(pmax(t0 - runif(10, 0, .1), 1e-3), t0, t0 + runif(20, 0, 1.5))
    R <- sample(1:2, length(rt), TRUE)
    ref <- port$ddm_eval(fl, th, rt, R)
    lp <- log(m()$dfun(rt, R, nat_of(th, length(rt))))
    expect_identical(is.finite(lp), is.finite(ref$log_pdf))
    fin <- is.finite(ref$log_pdf)
    expect_lt(max(abs(lp[fin] - ref$log_pdf[fin])), 1e-7)
    if ("st0" %in% regimes[[g]]) expect_true(all(!is.finite(lp[rt <= t0])) && all(is.finite(lp[rt > t0])))
    else expect_true(all(is.finite(lp)))
  }
})

test_that("ensembles apply the support rule too", {
  m2 <- register_nn_model(zb_bundle(2), twin = DDM)
  th <- c(v = 1, a = log(1), t0 = log(.3), s = 0, Z = 0, SZ = -Inf, sv = -Inf, st0 = -Inf)
  rt <- c(.2, .3, .5, .9); R <- c(1L, 2L, 1L, 2L)
  expect_equal(m2()$dfun(rt, R, nat_of(th, 4)), m()$dfun(rt, R, nat_of(th, 4)))
  expect_identical(m2()$dfun(rt, R, nat_of(th, 4))[1:2], c(0, 0))
})

test_that("compiled likelihood == R path with omitted st0 and constant sv = -Inf; mapped zeros", {
  des <- suppressMessages(design(model = m, factors = list(subjects = 1, S = c("left", "right")),
                                 Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
                                 formula = list(v ~ S, a ~ 1, t0 ~ 1, Z ~ 1, SZ ~ 1),
                                 constants = c(sv = -Inf)))
  p <- sampled_pars(des); p[] <- c(-1, 1.5, 0, log(.3), 0, qnorm(.2))
  mp <- mapped_pars(des, p)
  expect_true(all(mp[, "sv"] == 0) && all(mp[, "st0"] == 0))
  set.seed(9)
  dat <- make_data(p, des, n_trials = 50)
  dat$rt[1] <- 0.2                                  # below t0: zero density when st0 = 0
  emc <- suppressMessages(make_emc(dat, des, type = "single", compress = FALSE))
  dadm <- emc[[1]]$data[[1]]
  ll_native <- c(EMC2:::calc_ll_manager(t(p), dadm, emc[[1]]$model))
  pars <- EMC2:::get_pars_matrix_oo(t(p), dadm, emc[[1]]$model)
  ll_R <- sum(pmax(log(m()$dfun(dadm$rt, dadm$R, pars)), log(1e-10)))
  expect_equal(ll_native, ll_R, tolerance = 1e-10)
  # a value between 0 and DDM's floor is rejected by the bound, 0 is not
  des2 <- suppressMessages(design(model = m, factors = list(subjects = 1, S = c("left", "right")),
                                  Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
                                  formula = list(v ~ S, a ~ 1, t0 ~ 1, Z ~ 1, SZ ~ 1, sv ~ 1)))
  e2 <- suppressMessages(make_emc(dat, des2, type = "single", compress = FALSE))
  q <- c(p, sv = log(.005))
  expect_equal(c(EMC2:::calc_ll_manager(t(q), e2[[1]]$data[[1]], e2[[1]]$model)), log(1e-10) * nrow(dat))
  q[["sv"]] <- -Inf
  expect_equal(c(EMC2:::calc_ll_manager(t(q), e2[[1]]$data[[1]], e2[[1]]$model)), ll_native, tolerance = 1e-10)
})

test_that("malformed zerobox cards are refused", {
  b <- readRDS(zb_file)
  no_c <- b; no_c$members[[1]]$zerobox_c <- NULL
  f1 <- tempfile(fileext = ".rds"); saveRDS(no_c, f1)
  expect_error(register_nn_model(f1, twin = DDM), "zerobox_c gives no constant")
  extra <- b; extra$members[[1]]$zerobox_c$v <- 0.1
  f2 <- tempfile(fileext = ".rds"); saveRDS(extra, f2)
  expect_error(register_nn_model(f2, twin = DDM), "does not mark")
  pos <- function() { l <- DDM(); l$p_types[["sv"]] <- log(.05); l }
  expect_error(register_nn_model(zb_file, twin = pos), "needs the defaults")
  expect_error(register_nn_model(zb_file, p_types = c(DDM()$p_types[setdiff(names(DDM()$p_types), "st0")],
                                                      st0 = log(.05))), "needs the defaults")
})
