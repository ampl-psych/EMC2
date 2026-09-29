# MTLNR (correlated multiple-threshold LNR) and the rating layer.
# Ground truth for the density: DMC's n1PDF.3CLNRcorr (golden/mtlnr/dmc_density.rds,
# Stage 0) where its computation is well conditioned, and a 256-bit MPFR evaluation
# of the same formula (golden/mtlnr/mpfr_density.rds) everywhere.

golden_mtlnr <- function(f) readRDS(test_path("golden", "mtlnr", f))

# Two rows per trial (accumulators in level order): m, s, t0, rho, d1 .. d{K-1}.
# d is a 2 x (K - 1) matrix, row = accumulator.
mtlnr_pars <- function(n, m, s, t0, rho, d = NULL) {
  K1 <- if (is.null(d)) 0 else ncol(d)
  p <- cbind(m = rep(m, n), s = rep(s, n), t0 = t0, rho = rho)
  if (K1 > 0) {
    dd <- d[rep(1:2, n), , drop = FALSE]
    colnames(dd) <- paste0("d", seq_len(K1))
    p <- cbind(p, dd)
  }
  p
}

mtlnr_f <- function(r, c, K, m, s, t0, rho, d)
  function(rt) EMC2:::dMTLNR(rt, rep(r, length(rt)), rep(c, length(rt)),
                             mtlnr_pars(length(rt), m, s, t0, rho, d), K)

mtlnr_prob <- function(r, c, K, m, s, t0, rho, d, upper = Inf)
  stats::integrate(mtlnr_f(r, c, K, m, s, t0, rho, d), t0, upper,
                   subdivisions = 2000L, rel.tol = 1e-10)$value

# The golden grid as dMTLNR input
golden_eval <- function(g) {
  n <- nrow(g)
  pars <- matrix(NA, 2 * n, 6, dimnames = list(NULL, c("m", "s", "t0", "rho", "d1", "d2")))
  pars[seq(1, 2 * n, 2), ] <- cbind(g$meanlog_acc1, g$sdlog_acc1, g$t0, g$rho, g$d1_acc1, g$d2_acc1)
  pars[seq(2, 2 * n, 2), ] <- cbind(g$meanlog_acc2, g$sdlog_acc2, g$t0, g$rho, g$d1_acc2, g$d2_acc2)
  EMC2:::dMTLNR(g$rt, g$choice, g$rating, pars, 3)
}

test_that("MTLNR density matches DMC where DMC is well conditioned, MPFR everywhere", {
  g <- golden_mtlnr("dmc_density.rds")$dmc_density
  ref <- golden_mtlnr("mpfr_density.rds")$mpfr_density
  d <- golden_eval(g)
  expect_true(all(is.finite(d)) && all(d > 0))
  # MPFR reference: every cell, including DMC's underflows and cancellations
  expect_lt(max(abs(d - ref) / ref), 1e-12)
  # DMC: condition number of its subtraction of two lower-tail CDFs (ratings 1, 2)
  w1 <- g$choice == 1
  mw <- ifelse(w1, g$meanlog_acc1, g$meanlog_acc2); sw <- ifelse(w1, g$sdlog_acc1, g$sdlog_acc2)
  ml <- ifelse(w1, g$meanlog_acc2, g$meanlog_acc1); sl <- ifelse(w1, g$sdlog_acc2, g$sdlog_acc1)
  d1l <- ifelse(w1, g$d1_acc2, g$d1_acc1); d2l <- ifelse(w1, g$d2_acc2, g$d2_acc1)
  mc <- ml + g$rho * (log(g$dt) - mw) * sl / sw
  sc <- sl * sqrt(1 - g$rho^2)
  P <- function(x) plnorm(x, mc, sc)
  t1 <- ifelse(g$rating == 2, P(g$dt / d1l), ifelse(g$rating == 1, P(g$dt / d2l), 0))
  t2 <- ifelse(g$rating == 2, P(g$dt / d2l), ifelse(g$rating == 1, P(g$dt), 0))
  cond <- ifelse(g$rating == 3, 1, (t1 + t2) / pmax(t1 - t2, 1e-300))
  well <- g$density > 0 & cond < 1e3
  expect_gt(sum(well), 1300)
  expect_lt(max(abs(d[well] - g$density[well]) / g$density[well]), 1e-12)
  # elsewhere the disagreement is DMC's error: it is never closer to MPFR than we are
  expect_true(all(abs(g$density - ref)[!well] >= abs(d - ref)[!well]))
})

test_that("rho = 0 gives the independent product form", {
  dt <- c(.02, .1, .3, .7, 1.5, 4)
  d <- rbind(c(.25, .7), c(.4, .8))
  for (r in 1:2) for (c in 1:3) {
    f <- mtlnr_f(r, c, 3, m = c(-.4, .2), s = c(.8, 1.3), t0 = .2, rho = 0, d = d)(.2 + dt)
    w <- r; l <- 3 - r
    full <- c(0, d[l, ], 1)
    hi <- full[3 - c + 2]; lo <- full[3 - c + 1]
    Sl <- function(x) plnorm(x, c(-.4, .2)[l], c(.8, 1.3)[l], lower.tail = FALSE)
    ind <- dlnorm(dt, c(-.4, .2)[w], c(.8, 1.3)[w]) * (Sl(dt / hi) - if (lo > 0) Sl(dt / lo) else 0)
    expect_equal(f, ind, tolerance = 1e-12)
  }
})

test_that("K = 1 is DMC's binary correlated LNR, and K = 3 sums to it over ratings", {
  # verbatim from DMC dmc/models/LNR-Correlated/dists.R (n1PDFfixedt0.clnr)
  n1PDFfixedt0.clnr=function(dt,meanlog,sdlog,corr)
  {
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
    if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,dim(dt)[2]),nrow=n_acc)
    if (!is.matrix(sdlog))     sdlog <- matrix(rep(sdlog,dim(dt)[2]),nrow=n_acc)
    if(length(corr==1)) corr<-rep(corr[1],dim(dt)[2])
    dt[1,] <- dlnorm(dt[1,],meanlog[1,],sdlog[1,])
    if (dim(meanlog)[1]==2) dt[1,] <- dt[1,]*plnorm(dt[2,],
      (meanlog[2,]+sdlog[2,]/sdlog[1,]*corr*(log(dt[2,])-meanlog[1,])),
      (sdlog[2,]*sqrt(1-corr^2)),lower.tail=FALSE)
    pmax(dt[1,],0)
  }
  dt <- c(.05, .1, .2, .4, .8, 1.5, 3)
  m <- c(-.3, -.1); s <- c(.7, 1.4); d <- rbind(c(.2, .5), c(.4, .8))
  for (rho in c(-.9, -.5, 0, .5, .9, .99)) for (r in 1:2) {
    o <- if (r == 1) 1:2 else 2:1
    ref <- n1PDFfixedt0.clnr(dt, m[o], s[o], rho)
    k1 <- mtlnr_f(r, 1, 1, m, s, .3, rho, NULL)(.3 + dt)
    ok <- ref > 1e-300
    expect_equal(k1[ok], ref[ok], tolerance = 1e-12)
    k3 <- Reduce(`+`, lapply(1:3, function(c) mtlnr_f(r, c, 3, m, s, .3, rho, d)(.3 + dt)))
    expect_equal(k3, k1, tolerance = 1e-12)
  }
})

test_that("the 2K defective densities integrate to 1", {
  cases <- list(
    list(K = 3, m = c(-.3, -.1), s = c(1, 1), rho = .5, d = rbind(c(.3, .6), c(.3, .6))),
    list(K = 3, m = c(-1.2, .4), s = c(.6, 1.1), rho = .99, d = rbind(c(.1, .35), c(.9, .95))),
    list(K = 3, m = c(-.3, -.1), s = c(.7, 1.4), rho = -.9, d = rbind(c(.2, .4), c(.5, .8))),
    list(K = 2, m = c(-.5, 0), s = c(.8, .9), rho = .3, d = rbind(.6, .7)),
    list(K = 4, m = c(-.5, 0), s = c(.8, .9), rho = .7, d = rbind(c(.2, .5, .8), c(.3, .6, .9))),
    list(K = 1, m = c(-.5, 0), s = c(.8, .9), rho = -.4, d = NULL))
  for (cs in cases) {
    tot <- sum(sapply(1:2, function(r) sapply(seq_len(cs$K), function(c)
      mtlnr_prob(r, c, cs$K, cs$m, cs$s, .25, cs$rho, cs$d))))
    expect_equal(tot, 1, tolerance = 1e-6)
  }
})

test_that("log density stays finite for rho near 1 and tiny category probabilities", {
  d <- rbind(c(.5, .99), c(.5, .99))
  for (rho in c(-.999999, .999999)) for (c in 1:3) {
    p <- mtlnr_pars(4, c(-1, 1), c(.5, .5), .2, rho, d)
    ld <- EMC2:::dMTLNR(.2 + c(.01, .1, 1, 10), rep(1, 4), rep(c, 4), p, 3, log = TRUE)
    expect_true(all(is.finite(ld)))
  }
})

# simulated proportions and rt quantiles against the integrated density
check_sim <- function(prop, quants, qp, n, K, m, s, t0, rho, d) {
  for (r in 1:2) for (c in 1:K) {
    p <- mtlnr_prob(r, c, K, m, s, t0, rho, d)
    expect_lt(abs(prop[r, c] - p), 5 * sqrt(p * (1 - p) / n) + 1e-6)
    if (p * n > 5000) {
      Fq <- sapply(quants[[r]][[c]], function(q) mtlnr_prob(r, c, K, m, s, t0, rho, d, upper = q)) / p
      expect_lt(max(abs(Fq - qp)), 0.01)
    }
  }
}

test_that("density agrees with DMC's simulator (1e6 trials per set)", {
  sims <- golden_mtlnr("dmc_density.rds")$sim_summaries
  for (ss in sims) {
    ps <- ss$params
    d <- cbind(ps$d1, ps$d2)
    check_sim(unclass(ss$proportions), ss$quantiles, ss$quantile_probs, ss$n, 3,
              ps$meanlog, ps$sdlog, ps$t0, ps$corr_v, d)
  }
})

test_that("rMTLNR agrees with the density", {
  sims <- golden_mtlnr("dmc_density.rds")$sim_summaries
  set.seed(20260930)
  n <- 2e5
  qp <- c(.1, .3, .5, .7, .9)
  for (ss in sims) {
    ps <- ss$params
    d <- cbind(ps$d1, ps$d2)
    pars <- mtlnr_pars(n, ps$meanlog, ps$sdlog, ps$t0, ps$corr_v, d)
    lR <- factor(rep(c("a", "b"), n))
    sim <- EMC2:::rMTLNR(lR, pars, 3)
    R <- as.integer(sim$R)
    prop <- table(factor(R, levels = 1:2), factor(sim$RR, levels = 1:3)) / n
    quants <- lapply(1:2, function(r) lapply(1:3, function(c)
      stats::quantile(sim$rt[R == r & sim$RR == c], qp)))
    check_sim(unclass(prop), quants, qp, n, 3, ps$meanlog, ps$sdlog, ps$t0, ps$corr_v, d)
  }
})

mtlnr_design <- function(model = MTLNR, rlev = c("left", "right"), ...) {
  S <- rlev
  design(factors = list(subjects = 1, S = S), Rlevels = rlev, model = model,
         matchfun = function(d) d$S == d$lR,
         formula = list(m ~ lM, s ~ 1, t0 ~ 1, rho ~ 1, c1 ~ lR, c2 ~ lR),
         report_p_vector = FALSE, ...)
}
mtlnr_p <- c(m = -.5, m_lMTRUE = -.6, s = log(.8), t0 = log(.3), rho = qnorm(.75),
             c1 = log(.4), c1_lRright = .1, c2 = log(.6), c2_lRright = -.2)

test_that("likelihood: compression on and off agree with the direct density", {
  des <- mtlnr_design()
  set.seed(5)
  dat <- make_data(mtlnr_p, des, n_trials = 400)
  expect_true(all(c("R", "rt", "RR") %in% names(dat)))
  # a coarse rt resolution makes many trials share rt, so only the RR key keeps them apart
  res <- .05
  dc <- EMC2:::design_model(dat, des, compress = TRUE, verbose = FALSE, rt_resolution = res)
  du <- EMC2:::design_model(dat, des, compress = FALSE, verbose = FALSE, rt_resolution = res)
  expect_lt(nrow(dc), nrow(du))
  pm <- t(as.matrix(mtlnr_p))
  llc <- EMC2:::calc_ll_manager(pm, dc, des$model)
  llu <- EMC2:::calc_ll_manager(pm, du, des$model)
  expect_equal(llc, llu, tolerance = 1e-12)
  # direct: natural-scale parameters by hand
  rho <- 2 * pnorm(mtlnr_p[["rho"]]) - 1
  cl <- exp(c(mtlnr_p[["c1"]], mtlnr_p[["c1"]] + mtlnr_p[["c1_lRright"]]))
  c2 <- exp(c(mtlnr_p[["c2"]], mtlnr_p[["c2"]] + mtlnr_p[["c2_lRright"]]))
  d <- cbind(d1 = exp(-(cl + c2)), d2 = exp(-cl))
  n <- nrow(dat)
  match1 <- dat$S == "left"
  m1 <- ifelse(match1, -1.1, -.5); m2 <- ifelse(match1, -.5, -1.1)
  pars <- cbind(m = as.vector(rbind(m1, m2)), s = .8, t0 = .3, rho = rho,
                d[rep(1:2, n), ])
  rt <- floor(dat$rt / res) * res
  direct <- sum(pmax(log(1e-10), EMC2:::dMTLNR(rt, as.integer(dat$R), dat$RR, pars, 3, log = TRUE)))
  expect_equal(llc, direct, tolerance = 1e-10)
})

test_that("mapped parameters report natural-scale thresholds", {
  des <- mtlnr_design()
  mp <- mapped_pars(des, mtlnr_p)
  left <- mp[mp$lR == "left", ][1, ]
  expect_equal(left$d2, exp(-.4), tolerance = 1e-3)
  expect_equal(left$d1, exp(-(.4 + .6)), tolerance = 1e-3)
  expect_equal(left$rho, .5, tolerance = 1e-3)
})

test_that("simulation keeps RR on the conditional, unconditional and RACE paths", {
  des <- mtlnr_design()
  set.seed(6)
  d1 <- make_data(mtlnr_p, des, n_trials = 50)
  d2 <- make_data(mtlnr_p, des, n_trials = 50, conditional_on_data = FALSE)
  for (dd in list(d1, d2)) {
    expect_true("RR" %in% names(dd))
    expect_true(all(dd$RR %in% 1:3))
  }
  d3 <- d1
  d3$RACE <- factor(rep(2, nrow(d3)))
  d3 <- make_data(mtlnr_p, des, data = d3)
  expect_true(all(d3$RR %in% 1:3))
  # K = 2 through a wrapper
  des2 <- design(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"),
                 model = function() MTLNR(n_ratings = 2), matchfun = function(d) d$S == d$lR,
                 formula = list(m ~ lM, s ~ 1, t0 ~ 1, rho ~ 1, c1 ~ lR), report_p_vector = FALSE)
  p2 <- c(m = -.5, m_lMTRUE = -.6, s = log(.8), t0 = log(.3), rho = 0, c1 = log(.4), c1_lRb = 0)
  dk2 <- make_data(p2, des2, n_trials = 50)
  expect_true(all(dk2$RR %in% 1:2))
})

test_that("rating layer: data checks, design refusals, fold/unfold, thresholds", {
  des <- mtlnr_design()
  set.seed(7)
  dat <- make_data(mtlnr_p, des, n_trials = 30)
  bad <- dat; bad$RR <- factor(bad$RR)
  expect_error(EMC2:::design_model(bad, des, verbose = FALSE), "RR must be numeric")
  bad <- dat; bad$RR[1] <- 4
  expect_error(EMC2:::design_model(bad, des, verbose = FALSE), "RR must lie in 1 .. 3")
  bad <- dat; bad$RR[1] <- 1.5
  expect_error(EMC2:::design_model(bad, des, verbose = FALSE), "integer-valued")
  bad <- dat; bad$RR[1] <- NA
  expect_error(EMC2:::design_model(bad, des, verbose = FALSE), "RR is NA")
  bad <- dat; bad$RR <- NULL
  expect_error(EMC2:::design_model(bad, des, verbose = FALSE), "rating column RR")
  bad <- dat; bad$RR[bad$RR == 2] <- 1
  expect_warning(EMC2:::design_model(bad, des, verbose = FALSE), "never occur")
  bad <- dat; bad$rt[1] <- NA
  expect_error(EMC2:::design_model(bad, des, verbose = FALSE), "censored")
  # design: 2 choices only, rho trial level, RR not a predictor, RR not a covariate
  expect_error(mtlnr_design(rlev = c("a", "b", "c")), "exactly 2 response levels")
  expect_error(suppressMessages(
    design(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"),
           model = MTLNR, formula = list(m ~ 1, rho ~ lR), report_p_vector = FALSE)),
    "trial-level")
  expect_error(design(data = dat, model = MTLNR, formula = list(m ~ RR), report_p_vector = FALSE),
               "cannot be used as predictors")
  dd <- design(data = dat, model = MTLNR, matchfun = function(d) d$S == d$lR,
               formula = list(m ~ lM, s ~ 1, t0 ~ 1, rho ~ 1, c1 ~ lR, c2 ~ lR),
               report_p_vector = FALSE)
  expect_false("RR" %in% dd$Fcovariates)
  expect_error(MTLNR(n_ratings = c(3, 2)), "unequal rating counts")
  # fold / unfold
  f <- EMC2:::rating_fold(dat$R, dat$RR, 3)
  expect_equal(levels(f), c("left.3", "left.2", "left.1", "right.1", "right.2", "right.3"))
  u <- EMC2:::rating_unfold(f)
  expect_equal(u$R, dat$R)
  expect_equal(u$RR, as.integer(dat$RR))
  # thresholds: ordered by construction, criterion 1 gives d_{K-1}
  pm <- cbind(c1 = c(.2, 1), c2 = c(.3, .1), c3 = c(.5, 2))
  th <- EMC2:::rating_add_thresholds(pm, 4)
  expect_equal(unname(th[, "d3"]), exp(-c(.2, 1)))
  expect_equal(unname(th[, "d1"]), exp(-c(1, 3.1)))
  expect_true(all(th[, "d1"] < th[, "d2"] & th[, "d2"] < th[, "d3"]))
  # rating from evidence: small proportion = high rating
  dm <- cbind(d1 = .3, d2 = .6)
  expect_equal(EMC2:::rating_from_evidence(c(.1, .4, .9), dm[rep(1, 3), ], 3), c(3, 2, 1))
})

test_that("an MTLNR fit samples and predict() keeps RR", {
  des <- mtlnr_design()
  set.seed(8)
  dat <- make_data(mtlnr_p, des, n_trials = 200)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 2, verbose = FALSE))
  emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 5), cores_for_chains = 1,
                 cores_per_chain = 1, verbose = FALSE)
  expect_true(all(is.finite(emc[[1]]$samples$subj_ll)))
  expect_equal(get_data(emc)$RR, dat$RR)
  pp <- predict(emc, n_post = 2, n_cores = 1)
  expect_true(all(c("R", "rt", "RR") %in% names(pp)))
  expect_equal(nrow(pp), 2 * nrow(dat))
  expect_true(all(pp$RR %in% 1:3))
})
