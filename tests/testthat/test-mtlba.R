# MTLBA (multiple-threshold LBA, uncorrelated: rho = 0) and the rating layer's
# don't-know response. Ground truth: golden/mtlba (Stage M0) -- the don't-know
# paper's n1PDF.DK, BETE's n1PDF_MTR_1 and two independent integrations of the
# kernel, plus 1e6-trial simulator summaries. Closed forms there are accurate
# absolutely, not relatively, in the far tail (density < 1e-8), so parity is
# absolute everywhere and relative only in the body.

golden_mtlba <- function(f) readRDS(test_path("golden", "mtlba", f))

# Two rows per trial: v, sv, A, b, t0 and threshold proportions d1 .. d{K-1}
# (d: 2 x (K - 1) matrix, row = accumulator); rho (one value) if given.
mtlba_pars <- function(n, v, sv, A, b, t0, d = NULL, rho = NULL) {
  r <- function(x) rep_len(x, 2 * n)
  p <- cbind(v = r(v), sv = r(sv), A = r(A), b = r(b), t0 = r(t0))
  if (!is.null(rho)) p <- cbind(p, rho = rho)
  if (!is.null(d) && ncol(d) > 0) {
    dd <- d[rep(1:2, n), , drop = FALSE]
    colnames(dd) <- paste0("d", seq_len(ncol(d)))
    p <- cbind(p, dd)
  }
  p
}

# Independent log-scale references (no G function, no closed form): the stopping
# density as an integral over the drift, and the band probability as the mean over
# the start point of a normal interval probability.
ref_logf <- function(t, A, b, v, sv) {
  lo <- (b - A) / t; hi <- b / t
  g <- function(x) log(x) + stats::dnorm(x, v, sv, log = TRUE)
  M <- max(g(seq(lo, hi, length.out = 201)))
  I <- stats::integrate(function(x) exp(g(x) - M), lo, hi, rel.tol = 1e-13, abs.tol = 0)$value
  log(I) + M - log(A) - stats::pnorm(v / sv, log.p = TRUE)
}
ref_logband <- function(t, A, v, sv, lo, hi) {
  lp <- function(s) {
    zl <- (pmax(0, (lo - s) / t) - v) / sv
    zh <- ((hi - s) / t - v) / sv
    EMC2:::log_pnorm_interval(zl, zh)
  }
  M <- max(lp(seq(0, A, length.out = 201)))
  I <- stats::integrate(function(s) exp(lp(s) - M), 0, A, rel.tol = 1e-12, abs.tol = 0)$value
  log(I / A) + M - stats::pnorm(v / sv, log.p = TRUE)
}

# set list(v, sv, A, b, t0, d[, rho]) -> density function of rt for one response
mtlba_f <- function(R, RR, K, rule, dk, st, quad = EMC2:::mtlba_quad_rule())
  function(rt) EMC2:::dMTLBA(rt, rep(R, length(rt)), rep(RR, length(rt)),
                             mtlba_pars(length(rt), st$v, st$sv, st$A, st$b, st$t0, st$d, st$rho),
                             K, rule, dk, quad = quad)

mtlba_prob <- function(R, RR, K, rule, dk, st, upper = Inf) {
  t0 <- min(st$t0)
  brk <- t0 + c(0, .1, .2, .35, .5, .8, 1.2, 2, 3.5, 6, Inf)
  brk <- unique(c(brk[brk < upper], upper))
  f <- mtlba_f(R, RR, K, rule, dk, st)
  sum(sapply(seq_len(length(brk) - 1), function(i)
    stats::integrate(f, brk[i], brk[i + 1], subdivisions = 2000L, rel.tol = 1e-11,
                     abs.tol = 0)$value))
}

# the responses of a model: (R, RR) in folded order, DK once
mtlba_responses <- function(K, dk) {
  out <- rbind(cbind(R = 1, RR = K:1), cbind(R = 2, RR = 1:K))
  if (dk) out <- out[!(out[, "R"] == 2 & out[, "RR"] == 1), , drop = FALSE]
  out
}

# BETE sets (mtlba_ref_lib.R): thresholds d_k = A + B (k/K)^pow -> proportions (k/K)^pow
bete_set <- function(st, K) {
  d <- if (K > 1) matrix(unlist(lapply(1:2, function(a) ((1:(K - 1)) / K)^st$pow[a])),
                         nrow = 2, byrow = TRUE) else NULL
  list(v = st$v, sv = st$sv, A = st$A, b = st$A + st$B, t0 = st$t0, d = d)
}
# BETE response r = 1 .. 2K (2K - 1 with don't know) -> (R, RR)
bete_resp <- function(r, K, dk = FALSE) {
  if (!dk) return(if (r <= K) c(1, K - r + 1) else c(2, r - K))
  if (r < K) c(1, K - r + 1) else if (r == K) c(1, 1) else c(2, r - K + 1)
}

test_that("response lookup reproduces BETE's kernel indices (BE and TC, K = 1 .. 4)", {
  g <- golden_mtlba("bete_density.rds")
  cl <- unique(g$cells[, c("K", "rule", "R", "sa", "tt", "lt", "ut")])
  for (i in seq_len(nrow(cl))) {
    x <- cl[i, ]
    rr <- bete_resp(x$R, x$K)
    lk <- EMC2:::mtlba_lookup(rr[1], rr[2], x$K, x$rule)
    expect_equal(c(lk$s, lk$j + 1, lk$lo + 1, lk$lo + 2), c(x$sa, x$tt, x$lt, x$ut))
  }
  # don't-know tables (odd number of responses): the middle code holds exactly the two
  # rating-1 kernels
  for (nm in c("BE3", "BE5", "BE7", "TC3", "TC5", "TC7")) {
    tab <- g$dL[[nm]]
    K <- (as.integer(substr(nm, 3, 3)) + 1) / 2
    rule <- substr(nm, 1, 2)
    mid <- tab[tab[, "RR"] == K, , drop = FALSE]
    ours <- t(sapply(1:2, function(R) unlist(EMC2:::mtlba_lookup(R, 1, K, rule))))
    ours <- cbind(ours[, "s"], ours[, "j"] + 1, ours[, "lo"] + 1)
    expect_setequal(paste(mid[, "sa"], mid[, "tt"], mid[, "lt"]),
                    paste(ours[, 1], ours[, 2], ours[, 3]))
  }
})

test_that("density matches the BETE / reference grid at rho = 0", {
  g <- golden_mtlba("ref_density.rds")
  cells <- g$cells[g$cells$rho == 0, ]
  ours <- numeric(nrow(cells))
  for (i in seq_len(nrow(cells))) {
    x <- cells[i, ]
    st <- bete_set(g$sets[[x$set]], x$K)
    rr <- bete_resp(x$R, x$K)
    ours[i] <- mtlba_f(rr[1], rr[2], x$K, x$rule, FALSE, st)(x$rt)
  }
  expect_true(all(is.finite(ours)) && all(ours >= 0))
  expect_lt(max(abs(ours - cells$dens_2d)), 1e-13)
  expect_lt(max(abs(ours - cells$dens_bete)), 1e-13)
  body <- cells$dens_2d > 1e-8
  expect_gt(sum(body), 900)
  expect_lt(max(abs(ours[body] / cells$dens_2d[body] - 1)), 1e-10)
})

test_that("don't-know density matches the DK paper (n1PDF.DK)", {
  g <- golden_mtlba("dk_density.rds")
  cells <- g$cells
  ours <- numeric(nrow(cells))
  for (i in seq_len(nrow(cells))) {
    x <- cells[i, ]
    ps <- g$sets[[x$set]]
    st <- list(v = ps$v, sv = ps$sv, A = ps$A, b = ps$b, t0 = ps$t0,
               d = matrix((ps$d - ps$A) / (ps$b - ps$A), ncol = 1))
    ours[i] <- mtlba_f(x$R, x$RR, 2, "BE", TRUE, st)(x$rt)
  }
  expect_lt(max(abs(ours - cells$dens)), 1e-13)
  body <- cells$kernel_ref > 1e-8
  expect_lt(max(abs(ours[body] / cells$kernel_ref[body] - 1)), 1e-10)
})

test_that("K = 1 and the don't-know response are sums of the finer responses", {
  st <- list(v = c(2, 1.2), sv = c(1, .8), A = c(.5, .7), b = c(1.4, 1.6), t0 = .2,
             d = rbind(c(.3, .7), c(.2, .5)))
  rt <- .2 + c(.05, .15, .3, .6, 1, 2, 4)
  for (rule in c("BE", "TC")) {
    # don't know = sum of the two rating-1 terms
    dk <- mtlba_f(1, 1, 3, rule, TRUE, st)(rt)
    two <- mtlba_f(1, 1, 3, rule, FALSE, st)(rt) + mtlba_f(2, 1, 3, rule, FALSE, st)(rt)
    expect_equal(dk, two, tolerance = 1e-13)
    # R is ignored on a don't-know trial
    expect_equal(mtlba_f(2, 1, 3, rule, TRUE, st)(rt), dk, tolerance = 1e-15)
    # K = 1: BE and TC coincide
    s1 <- st; s1$d <- NULL
    expect_equal(mtlba_f(1, 1, 1, rule, FALSE, s1)(rt), mtlba_f(1, 1, 1, "BE", FALSE, s1)(rt))
  }
  # BE: summing over ratings gives K = 1
  s1 <- st; s1$d <- NULL
  for (R in 1:2) {
    k3 <- Reduce(`+`, lapply(1:3, function(RR) mtlba_f(R, RR, 3, "BE", FALSE, st)(rt)))
    expect_equal(k3, mtlba_f(R, 1, 1, "BE", FALSE, s1)(rt), tolerance = 1e-12)
  }
})

test_that("the 2K (2K - 1 with don't know) defective densities integrate to 1", {
  sets <- list(
    list(v = c(2.5, 1.5), sv = c(1, 1), A = c(.5, .5), b = c(1.5, 1.5), t0 = .2),
    list(v = c(3, 1), sv = c(.7, 1.2), A = c(.4, .8), b = c(1.6, 1.5), t0 = c(.15, .25)),
    list(v = c(.8, -.3), sv = c(1, .6), A = c(1, .3), b = c(1.5, 1.8), t0 = .3),
    list(v = c(1.5, 1), sv = c(1, 1), A = c(0, 0), b = c(1, 1.2), t0 = .2))
  for (st in sets) for (K in 1:4) for (rule in c("BE", "TC")) for (dk in c(FALSE, if (K > 1) TRUE)) {
    set.seed(K)
    s <- st
    if (K > 1) s$d <- matrix(t(apply(matrix(runif(2 * (K - 1), .05, .95), 2), 1, sort)), nrow = 2,
                             byrow = K == 2)
    rsp <- mtlba_responses(K, dk)
    tot <- sum(apply(rsp, 1, function(r) mtlba_prob(r[1], r[2], K, rule, dk, s)))
    expect_equal(tot, 1, tolerance = 1e-7)
  }
})

test_that("log density keeps relative accuracy in the far tail", {
  cases <- expand.grid(t = c(.02, .1, 2, 20), lo = c(0, .9, 1.2), v = c(-1, 1, 3))
  for (i in seq_len(nrow(cases))) {
    x <- cases[i, ]
    hi <- if (x$lo == 0) .9 else if (x$lo == .9) 1.2 else 1.5
    ours <- EMC2:::mtlba_logband(x$t, .5, x$v, .4, x$lo, hi, x$lo == 0)
    ref <- ref_logband(x$t, .5, x$v, .4, x$lo, hi)
    expect_true(is.finite(ours))
    expect_lt(abs(ours - ref), 1e-8 * max(1, abs(ref)))
  }
  # stopping density
  for (t in c(.02, .1, .5, 2, 20)) for (v in c(-1, 1, 3)) {
    ours <- EMC2:::mtlba_logf(t, .5, 1.5, v, .4)
    ref <- ref_logf(t, .5, 1.5, v, .4)
    expect_lt(abs(ours - ref), 1e-8 * max(1, abs(ref)))
  }
  # threshold just above A, small sv: the pnorm difference is ~1 while both dnorms
  # are ~exp(-1000) (overflowed to +Inf before M2)
  ours <- EMC2:::mtlba_logf(.257, .5, .5001, 1, .02)
  expect_true(is.finite(ours))
  expect_lt(abs(ours - ref_logf(.257, .5, .5001, 1, .02)), 1e-8)
})

test_that("non-finite parameters give -Inf, not an error", {
  st <- list(v = c(2, 1), sv = c(1, 1), A = c(.5, .5), b = c(1.5, 1.5), t0 = .2,
             d = rbind(c(.3, .6), c(.3, .6)))
  for (bad in list(c(sv = Inf), c(A = Inf), c(v = NaN), c(b = Inf), c(sv = 0))) {
    pars <- mtlba_pars(3, st$v, st$sv, st$A, st$b, st$t0, st$d)
    pars[2, names(bad)] <- bad
    for (rule in c("BE", "TC")) {
      ld <- EMC2:::dMTLBA(rep(.7, 3), c(1, 2, 1), c(3, 2, 1), pars, 3, rule, log = TRUE)
      expect_true(!anyNA(ld) && all(is.finite(ld[2:3])))
    }
  }
})

test_that("density agrees with the reference simulators (1e6 trials per set)", {
  s <- golden_mtlba("sim_summaries.rds")
  check <- function(p, q, qp, n, rsp, K, rule, dk, st) {
    for (i in seq_len(nrow(rsp))) {
      pe <- mtlba_prob(rsp[i, 1], rsp[i, 2], K, rule, dk, st)
      expect_lt(abs(p[i] - pe), 4.5 * sqrt(pe * (1 - pe) / n) + 1e-7)
      if (pe * n > 20000) {
        Fq <- sapply(q[i, ], function(u) mtlba_prob(rsp[i, 1], rsp[i, 2], K, rule, dk, st, upper = u)) / pe
        expect_lt(max(abs(Fq - qp)), 4.5 * .5 / sqrt(pe * n))
      }
    }
  }
  # don't-know paper simulator: responses 1, 2 (definitive), DK
  for (x in s$dk) {
    ps <- x$pars
    st <- list(v = ps$v, sv = ps$sv, A = ps$A, b = ps$b, t0 = ps$t0,
               d = matrix((ps$d - ps$A) / (ps$b - ps$A), ncol = 1))
    check(x$p, x$q, x$probs, x$n, rbind(c(1, 2), c(2, 2), c(1, 1)), 2, "BE", TRUE, st)
  }
  # BETE simulators at rho = 0, with and without don't know
  for (x in s$bete[grepl("_0$", names(s$bete))]) {
    K <- x$config$K; rule <- as.character(x$config$rule)
    pm <- x$pars
    d <- if (K > 1) (pm[, paste0("DT", 1:(K - 1)), drop = FALSE] - pm[, "A"]) / (pm[, "b"] - pm[, "A"])
    st <- list(v = pm[, "v"], sv = pm[, "sv"], A = pm[, "A"], b = pm[, "b"], t0 = pm[, "t0"], d = d)
    rsp <- t(sapply(seq_len(2 * K), bete_resp, K = K))
    check(x$p, x$q, x$probs, x$n, rsp, K, rule, FALSE, st)
    y <- x$dont_know
    rsp <- t(sapply(seq_len(2 * K - 1), bete_resp, K = K, dk = TRUE))
    check(y$p, y$q, y$probs, y$n, rsp, K, rule, TRUE, st)
  }
})

test_that("rMTLBA agrees with the density (BE, TC, don't know)", {
  st <- list(v = c(2.2, 1.3), sv = c(.8, 1.1), A = c(.5, .7), b = c(1.4, 1.6), t0 = c(.2, .25),
             d = rbind(c(.3, .7), c(.2, .5)))
  n <- 2e5
  qp <- c(.1, .3, .5, .7, .9)
  lR <- factor(rep(c("a", "b"), n))
  for (rule in c("BE", "TC")) for (dk in c(FALSE, TRUE)) {
    set.seed(if (rule == "BE") 31 else 32)
    sim <- EMC2:::rMTLBA(lR, mtlba_pars(n, st$v, st$sv, st$A, st$b, st$t0, st$d), 3, rule, dk)
    expect_true(all(sim$RR %in% 1:3) && !anyNA(sim$R))
    if (dk) expect_true(all(sim$R[sim$RR == 1] == "a"))
    rsp <- mtlba_responses(3, dk)
    for (i in seq_len(nrow(rsp))) {
      sel <- sim$RR == rsp[i, 2] & (if (dk && rsp[i, 2] == 1) TRUE else as.integer(sim$R) == rsp[i, 1])
      pe <- mtlba_prob(rsp[i, 1], rsp[i, 2], 3, rule, dk, st)
      expect_lt(abs(mean(sel) - pe), 4.5 * sqrt(pe * (1 - pe) / n))
      if (pe * n > 5000) {
        Fq <- sapply(stats::quantile(sim$rt[sel], qp), function(u)
          mtlba_prob(rsp[i, 1], rsp[i, 2], 3, rule, dk, st, upper = u)) / pe
        expect_lt(max(abs(Fq - qp)), 4.5 * .5 / sqrt(pe * n))
      }
    }
  }
})

mtlba_design <- function(K = 3, rule = "BE", dk = FALSE, rlev = c("left", "right"), ...) {
  crit <- if (K > 1) paste0("c", seq_len(K - 1)) else character(0)
  form <- c(list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
            lapply(crit, function(cn) stats::as.formula(paste(cn, "~ lR"))))
  mod <- local({ k <- K; r <- rule; d <- dk; function() MTLBA(n_ratings = k, dont_know = d, rule = r) })
  design(factors = list(subjects = 1, S = rlev), Rlevels = rlev, model = mod,
         matchfun = function(d) d$S == d$lR, formula = form, report_p_vector = FALSE, ...)
}
mtlba_p <- function(K) {
  p <- c(v = 1, v_lMTRUE = 1.5, B = log(1), A = log(.5), t0 = log(.25))
  for (k in seq_len(K - 1)) p[paste0("c", k, c("", "_lRright"))] <- c(log(.6), .15)
  p
}

test_that("K = 1 equals the LBA log-likelihood (BE and TC)", {
  set.seed(40)
  lba <- design(factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
                model = LBA, matchfun = function(d) d$S == d$lR,
                formula = list(v ~ lM, sv ~ lM, B ~ lR, A ~ 1, t0 ~ 1), report_p_vector = FALSE)
  p <- c(v = 1, v_lMTRUE = 1.5, sv = 0, sv_lMTRUE = log(.7), B = log(1), B_lRright = .2,
         A = log(.5), t0 = log(.25))
  dat <- make_data(p, lba, n_trials = 300)
  dat$RR <- 1
  props <- rbind(p, mvtnorm::rmvnorm(20, p, diag(length(p)) * .1))
  colnames(props) <- names(p)
  d_lba <- EMC2:::design_model(dat, lba, compress = FALSE, verbose = FALSE)
  ll_lba <- EMC2:::calc_ll_manager(props, d_lba, lba$model)
  tw_lba <- EMC2:::calc_ll_manager(props, d_lba, lba$model, return_trialwise = TRUE)
  for (rule in c("BE", "TC")) {
    mt <- suppressMessages(design(
      factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
      model = local({ r <- rule; function() MTLBA(1, rule = r) }),
      matchfun = function(d) d$S == d$lR,
      formula = list(v ~ lM, sv ~ lM, B ~ lR, A ~ 1, t0 ~ 1), report_p_vector = FALSE))
    d_mt <- EMC2:::design_model(dat, mt, compress = FALSE, verbose = FALSE)
    ll <- EMC2:::calc_ll_manager(props, d_mt, mt$model)
    tw <- sapply(seq_len(nrow(props)), function(k)
      EMC2:::trial_ll_mtlba(EMC2:::get_pars_matrix_oo(props[k, ], d_mt, mt$model()), d_mt, 1,
                            rule, FALSE))
    expect_equal(as.vector(ll), colSums(tw), tolerance = 1e-12)
    # identical where EMC2's LBA is well conditioned ...
    tl <- matrix(tw_lba, ncol = nrow(props))
    body <- tl > -10
    expect_lt(max(abs(tw - tl)[body]), 1e-10)
    # ... and elsewhere (LBA's 1 - CDF loses digits in its far tail) ours is the
    # accurate one
    expect_lt(max(abs(tw - tl)), 1e-6)
    off <- which(abs(tw - tl) > 1e-10, arr.ind = TRUE)
    for (r in utils::head(seq_len(nrow(off)), 10)) {
      i <- off[r, 1]; k <- off[r, 2]
      pm <- EMC2:::get_pars_matrix_oo(props[k, ], d_mt, mt$model())
      rows <- c(2 * i - 1, 2 * i)
      w <- as.integer(d_mt$R[rows[1]])
      pw <- pm[rows[w], ]; pl <- pm[rows[3 - w], ]
      t <- d_mt$rt[rows[1]] - pw[["t0"]]
      ref <- ref_logf(t, pw[["A"]], pw[["b"]], pw[["v"]], pw[["sv"]]) +
        ref_logband(t, pl[["A"]], pl[["v"]], pl[["sv"]], 0, pl[["b"]])
      expect_lt(abs(tw[i, k] - ref), 1e-12)
      expect_lt(abs(tw[i, k] - ref), abs(tl[i, k] - ref))
    }
  }
})

test_that("likelihood: compression on and off agree with the direct density", {
  for (rule in c("BE", "TC")) for (dk in c(FALSE, TRUE)) {
    K <- if (dk) 2 else 3
    des <- mtlba_design(K, rule, dk)
    p <- mtlba_p(K)
    set.seed(50)
    dat <- make_data(p, des, n_trials = 300)
    if (dk) expect_true(all(dat$R[dat$RR == 1] == "left"))
    res <- .05
    dc <- EMC2:::design_model(dat, des, compress = TRUE, verbose = FALSE, rt_resolution = res)
    du <- EMC2:::design_model(dat, des, compress = FALSE, verbose = FALSE, rt_resolution = res)
    expect_lt(nrow(dc), nrow(du))
    pm <- t(as.matrix(p))
    llc <- EMC2:::calc_ll_manager(pm, dc, des$model)
    llu <- EMC2:::calc_ll_manager(pm, du, des$model)
    expect_equal(llc, llu, tolerance = 1e-12)
    # direct, natural-scale parameters by hand
    n <- nrow(dat)
    m1 <- dat$S == "left"
    v <- as.vector(rbind(ifelse(m1, 2.5, 1), ifelse(m1, 1, 2.5)))
    C <- matrix(sapply(seq_len(K - 1), function(k) exp(p[[paste0("c", k)]] + c(0, .15))), nrow = 2)
    if (K > 2) for (k in 2:(K - 1)) C[, k] <- C[, k - 1] + C[, k]
    d <- exp(-C[, rev(seq_len(K - 1)), drop = FALSE])
    pars <- mtlba_pars(n, 0, 1, .5, 1.5, .25, d)
    pars[, "v"] <- v
    rt <- floor(dat$rt / res) * res
    direct <- sum(pmax(log(1e-10), EMC2:::dMTLBA(rt, as.integer(dat$R), dat$RR, pars, K, rule, dk,
                                                  log = TRUE)))
    expect_equal(as.vector(llc), direct, tolerance = 1e-10)
  }
})

test_that("design and data checks: don't-know coding, rho, K = 1 with don't know", {
  expect_error(MTLBA(1, dont_know = TRUE), "n_ratings >= 2")
  expect_error(MTLBA(3, rule = "XX"))
  des <- mtlba_design(2, "BE", TRUE)
  set.seed(60)
  dat <- make_data(mtlba_p(2), des, n_trials = 50)
  bad <- dat
  bad$R[bad$RR == 1][1] <- "right"
  expect_error(EMC2:::design_model(bad, des, verbose = FALSE), "first response level")
  # without don't know the same data are fine
  expect_no_error(EMC2:::design_model(bad, suppressMessages(mtlba_design(2, "BE", FALSE)),
                                      verbose = FALSE))
  bad <- dat; bad$rt[1] <- NA
  expect_error(EMC2:::design_model(bad, des, verbose = FALSE), "censored")
  expect_error(suppressMessages(
    design(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = MTLBA,
           formula = list(v ~ 1, rho ~ lR), report_p_vector = FALSE)), "trial-level")
  expect_error(mtlba_design(3, rlev = c("a", "b", "c")), "exactly 2 response levels")
  mp <- mapped_pars(mtlba_design(3), mtlba_p(3))
  expect_equal(unique(mp$rho), 0)
  expect_equal(mp$b, mp$A + mp$B)
})

test_that("rating layer: don't-know fold/unfold, summaries and zROC", {
  R <- factor(c("n", "o", "n", "o", "n"), levels = c("n", "o"))
  RR <- c(3, 3, 1, 2, 2)
  f <- EMC2:::rating_fold(R, RR, 3, dont_know = TRUE)
  expect_equal(levels(f), c("n.3", "n.2", "DK", "o.2", "o.3"))
  expect_equal(as.character(f), c("n.3", "o.3", "DK", "o.2", "n.2"))
  u <- EMC2:::rating_unfold(f)
  expect_equal(u$R, R)
  expect_equal(u$RR, as.integer(RR))
  des <- mtlba_design(2, "TC", TRUE)
  set.seed(61)
  dat <- make_data(mtlba_p(2), des, n_trials = 200)
  s <- rating_summary(dat, factors = "S", dont_know = TRUE)
  expect_equal(levels(s$resp), c("left.2", "DK", "right.2"))
  expect_equal(as.vector(tapply(s$p, s$cell, sum)), c(1, 1))
  z <- zroc(dat, "S", "left", dont_know = TRUE)
  expect_equal(nrow(z), 2)
  nz <- dat[dat$S == "left", ]
  expect_equal(z$F[1], mean(!(nz$R == "left" & nz$RR == 2)))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_ratings(dat, factors = "S", dont_know = TRUE))
  expect_silent(plot_zroc(dat, "S", "left", dont_know = TRUE))
})

test_that("MTLBA fits sample and predict() keeps the don't-know coding", {
  des <- mtlba_design(2, "BE", TRUE)
  set.seed(70)
  dat <- make_data(mtlba_p(2), des, n_trials = 150)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 2, verbose = FALSE))
  emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 5), cores_for_chains = 1,
                 cores_per_chain = 1, verbose = FALSE)
  expect_true(all(is.finite(emc[[1]]$samples$subj_ll)))
  expect_equal(get_data(emc)$RR, dat$RR)
  pp <- predict(emc, n_post = 2, n_cores = 1)
  expect_true(all(c("R", "rt", "RR") %in% names(pp)))
  expect_equal(nrow(pp), 2 * nrow(dat))
  expect_true(all(pp$RR %in% 1:2))
  expect_true(all(pp$R[pp$RR == 1] == "left"))
  # the emc object tells the rating helpers about don't know
  s <- rating_summary(emc)
  expect_true("DK" %in% levels(s$resp))
  gp <- get_pars(emc, map = TRUE, add_recalculated = TRUE, selection = "alpha",
                 return_mcmc = FALSE, merge_chains = TRUE)
  expect_equal(gp["d1_lRright", 1, ], exp(-gp["c1_lRright", 1, ]), tolerance = 1e-6)
})

# ---- C++ kernel (src/model_mtlba.h, c_name "MTLBA") ----

mtlba_r_model <- function(model) {
  ml <- model()
  ml$c_name <- NULL
  function() ml
}
# error relative to max(1, |reference|): trial log-likelihoods can be near 0
mtlba_err <- function(a, b) max(abs(a - b) / pmax(abs(b), 1))
mtlba_omp <- function(expr, threads = 4) {
  old <- options(emc.ll_backend = "multithreaded", emc.n_threads = threads)
  on.exit(options(old))
  expr
}

test_that("C++ likelihood equals the R reference (K = 1 .. 4, BE/TC, don't know, extremes)", {
  for (K in 1:4) for (rule in c("BE", "TC")) for (dk in c(FALSE, if (K > 1) TRUE)) {
    crit <- if (K > 1) paste0("c", seq_len(K - 1)) else character(0)
    form <- c(list(v ~ lM, sv ~ lM, B ~ lR, A ~ 1, t0 ~ lR),
              lapply(crit, function(cn) stats::as.formula(paste(cn, "~ lR"))))
    mod <- local({ k <- K; r <- rule; d <- dk
      function() MTLBA(n_ratings = k, dont_know = d, rule = r) })
    des <- suppressMessages(design(
      factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
      model = mod, matchfun = function(d) d$S == d$lR, formula = form, report_p_vector = FALSE))
    p <- sampled_pars(des, doMap = FALSE)
    p[] <- 0
    p[c("v", "v_lMTRUE", "sv", "sv_lMTRUE", "B", "B_lRright", "A", "t0", "t0_lRright")] <-
      c(1, 1.5, 0, log(.7), log(1), .2, log(.5), log(.25), .1)
    for (cn in crit) p[c(cn, paste0(cn, "_lRright"))] <- c(log(.6), .15)
    set.seed(100 + 10 * K + dk)
    dat <- make_data(p, des, n_trials = 120)
    props <- rbind(p, mvtnorm::rmvnorm(30, p, diag(length(p)) * .3))
    colnames(props) <- names(p)
    ext <- matrix(p, nrow = 10, ncol = length(p), byrow = TRUE, dimnames = list(NULL, names(p)))
    ext[1, "sv"] <- log(.02)              # sharply peaked drift: far tails everywhere
    ext[2, "v"] <- -4                     # mean drift far below 0 (positive truncation)
    ext[3, "v"] <- 12                     # very fast
    ext[4, "A"] <- -Inf                   # A = 0 exactly (the bound's exception)
    ext[5, "A"] <- log(1e-3)              # tiny A (not the A = 0 limit)
    ext[6, "t0_lRright"] <- 1.5           # t0 above many rts on one accumulator
    ext[7, "B"] <- 800                    # B overflows: non-finite, out of bounds
    ext[8, "v"] <- NaN                    # non-finite drift
    ext[9, "sv_lMTRUE"] <- -800           # sv underflows to 0: out of bounds
    ext[10, "t0"] <- log(.01)             # t0 below its bound: every trial at min_ll
    if (K > 1) {
      ext[1, "c1"] <- log(8)              # highest rating nearly certain
      ext[2, "c1"] <- log(1e-4)           # highest rating nearly impossible
      ext[3, crit[K - 1]] <- 700          # a criterion near overflow: an empty band
    }
    props <- rbind(props, ext)
    for (compress in c(TRUE, FALSE)) {
      dadm <- EMC2:::design_model(dat, des, compress = compress, verbose = FALSE)
      llr <- EMC2:::calc_ll_manager(props, dadm, mtlba_r_model(des$model))
      llc <- EMC2:::calc_ll_manager(props, dadm, des$model)
      expect_lt(mtlba_err(llc, llr), 1e-12)
      twr <- apply(props, 1, function(pv)
        EMC2:::trial_ll_mtlba(EMC2:::get_pars_matrix_oo(pv, dadm, des$model()), dadm, K, rule, dk))
      twc <- EMC2:::calc_ll_manager(props, dadm, des$model, return_trialwise = TRUE)
      expect_lt(mtlba_err(as.vector(twc), as.vector(twr)), 1e-11)
      n_ext <- nrow(props) - 0:3
      expect_equal(as.vector(llc[n_ext]), rep(log(1e-10) * nrow(dat), 4))
      expect_true(all(is.finite(llc)))
      # OpenMP gives exactly the serial numbers
      expect_identical(mtlba_omp(EMC2:::calc_ll_manager(props, dadm, des$model)), llc)
      expect_identical(mtlba_omp(EMC2:::calc_ll_manager(props, dadm, des$model,
                                                        return_trialwise = TRUE)), twc)
    }
  }
})

test_that("C++ likelihood refuses what it does not define", {
  des <- mtlba_design(3, "TC")
  set.seed(80)
  dat <- make_data(mtlba_p(3), des, n_trials = 30)
  dadm <- EMC2:::design_model(dat, des, verbose = FALSE)
  pm <- t(as.matrix(mtlba_p(3)))
  for (omp in c(FALSE, TRUE)) {
    run <- function(d, model = des$model)
      if (omp) mtlba_omp(EMC2:::calc_ll_manager(pm, d, model)) else EMC2:::calc_ll_manager(pm, d, model)
    expect_true(is.finite(run(dadm)))
    bad <- dadm; bad$missingness <- NA_integer_; bad$missingness[1:2] <- 2L; bad$UC <- 2
    expect_error(run(bad), "response retained")
    bad$R[1:2] <- NA; bad$missingness[1:2] <- 4L
    expect_error(run(bad), "missingness code 4")
    bad$missingness[1:2] <- 1L
    expect_error(run(bad), "no LC column")
    bad <- dadm; bad$missingness <- NA_integer_; bad$R[1:2] <- NA
    expect_error(run(bad), "not censored")
    bad <- dadm; bad$RR[1] <- 4
    expect_error(run(bad), "RR must be an integer in 1 .. 3")
    bad <- dadm; bad$RR <- NULL
    expect_error(run(bad), "no rating column")
    bad <- dadm; bad$winner[1:2] <- TRUE
    expect_error(run(bad), "missing response")
    bad <- dadm; attr(bad, "constants")[["rho"]] <- .3
    expect_error(run(bad), "rho != 0")
    bad <- dadm; attr(bad, "c_options") <- NULL
    expect_error(run(bad), "no c_options attribute")
  }
})

test_that("same-seed fits are identical on the serial and OpenMP paths", {
  des <- mtlba_design(2, "TC", TRUE)
  set.seed(81)
  dat <- make_data(mtlba_p(2), des, n_trials = 100)
  fit_once <- function() {
    set.seed(82)
    emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 2, verbose = FALSE))
    emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 5), cores_for_chains = 1,
                   cores_per_chain = 1, verbose = FALSE)
    lapply(emc, function(ch) ch$samples)
  }
  s_serial <- fit_once()
  s_omp <- mtlba_omp(fit_once())
  for (ch in 1:2) {
    expect_identical(s_omp[[ch]]$alpha, s_serial[[ch]]$alpha)
    expect_identical(s_omp[[ch]]$subj_ll, s_serial[[ch]]$subj_ll)
  }
})

# ---- censoring and truncation (Stage M2c) ----

# natural-scale parameters of one trial (two rows) for the direct R functions
mtlba_trial_pars <- function(K, t0 = c(.2, .2)) {
  p <- cbind(v = c(1.2, .6), sv = c(.8, 1), A = c(.5, .4), b = c(1.4, 1.2), t0 = t0)
  if (K > 1) {
    d <- rbind(seq(.3, .8, length.out = K - 1), seq(.2, .7, length.out = K - 1))
    colnames(d) <- paste0("d", seq_len(K - 1))
    p <- cbind(p, d)
  }
  p
}
mtlba_bnd <- function(miss = NA_integer_, LT = 0, UT = Inf, LC = 0, UC = Inf)
  list(miss = miss, LT = LT, UT = UT, LC = LC, UC = UC)

test_that("censored + observed mass is 1 and the truncated density integrates to 1", {
  LC <- .45; UC <- 1.6; LT <- .3; UT <- 2.5
  for (rule in c("BE", "TC")) for (K in 1:4) for (dk in c(FALSE, if (K > 1) TRUE)) {
    p <- mtlba_trial_pars(K)
    rsp <- mtlba_responses(K, dk)
    ll_miss <- function(rt, R, RR, bnd)
      EMC2:::mtlba_logdens_missing(rt, R, RR, p, K, rule, dk, bnd)
    dens <- function(rt, R, RR, bnd = mtlba_bnd())
      vapply(rt, function(x) exp(ll_miss(x, R, RR, bnd)), numeric(1))
    obs <- sum(mapply(function(R, RR)
      stats::integrate(dens, LC, UC, R = R, RR = RR, rel.tol = 1e-11)$value, rsp[, "R"], rsp[, "RR"]))
    cens <- function(m, bnd = mtlba_bnd(m, LC = LC, UC = UC)) exp(ll_miss(NA, NA, NA, bnd))
    expect_equal(obs + cens(1L) + cens(2L), 1, tolerance = 1e-9)
    expect_equal(cens(3L), cens(1L) + cens(2L), tolerance = 1e-12)
    # truncated: observed density plus censored mass inside (LT, UT) is 1
    tb <- mtlba_bnd(LT = LT, UT = UT)
    tr <- sum(mapply(function(R, RR)
      stats::integrate(dens, LT, UT, R = R, RR = RR, bnd = tb, rel.tol = 1e-11)$value,
      rsp[, "R"], rsp[, "RR"]))
    expect_equal(tr, 1, tolerance = 1e-9)
    tcb <- function(m) cens(m, mtlba_bnd(m, LT = LT, UT = UT, LC = LC, UC = UC))
    tr_obs <- sum(mapply(function(R, RR)
      stats::integrate(dens, LC, UC, R = R, RR = RR, bnd = tb, rel.tol = 1e-11)$value,
      rsp[, "R"], rsp[, "RR"]))
    expect_equal(tr_obs + tcb(1L) + tcb(2L), 1, tolerance = 1e-9)
    # the stopping-time survivor is 1 - the integrated response density
    y <- .7
    st <- EMC2:::mtlba_logstop(y, p, EMC2:::mtlba_thresholds(p, K), 1, 2, K, rule)
    below <- sum(mapply(function(R, RR)
      stats::integrate(dens, .2, .2 + y, R = R, RR = RR, rel.tol = 1e-11)$value, rsp[, "R"], rsp[, "RR"]))
    expect_equal(exp(st$P), below, tolerance = 1e-9)
    expect_equal(exp(st$S), 1 - below, tolerance = 1e-9)
  }
})

test_that("censoring window terms keep relative accuracy in the tails", {
  # a lower cut far below the mass and an upper cut far above it: P ~ 1e-20 or
  # less, where 1 - survivor would round to 0
  for (rule in c("BE", "TC")) {
    K <- 3
    p <- mtlba_trial_pars(K)
    lp <- EMC2:::mtlba_logdens_missing(NA, NA, NA, p, K, rule, FALSE, mtlba_bnd(1L, LC = .23))
    expect_true(is.finite(lp) && lp < -30)
    # (the LBA's upper tail is heavy unless a drift is well above 0)
    pf <- p
    pf[1, c("v", "sv")] <- c(3, .3)
    lu <- EMC2:::mtlba_logdens_missing(NA, NA, NA, pf, K, rule, FALSE, mtlba_bnd(2L, UC = 40))
    expect_true(is.finite(lu) && lu < -30)
    # the same through an independent split: P(dt < y) as an integral of the
    # response densities on the log scale
    y <- .03
    dens <- function(rt, R, RR)
      vapply(rt, function(x) exp(EMC2:::mtlba_logdens(x, R, RR, p, K, rule, FALSE) - lp), numeric(1))
    rsp <- mtlba_responses(K, FALSE)
    I <- sum(mapply(function(R, RR) stats::integrate(dens, .2, .2 + y, R = R, RR = RR,
                                                     rel.tol = 1e-11)$value, rsp[, "R"], rsp[, "RR"]))
    expect_equal(I, 1, tolerance = 1e-8)
  }
})

test_that("make_data with a censoring/truncation list round-trips (BE, TC, don't know)", {
  for (rule in c("BE", "TC")) for (dk in c(FALSE, TRUE)) {
    K <- if (dk) 2 else 3
    des <- mtlba_design(K, rule, dk)
    p <- mtlba_p(K)
    set.seed(90 + dk)
    LC <- .5; UC <- 1.3
    dat <- make_data(p, des, n_trials = 2000, TC = list(LC = LC, UC = UC))
    cl <- !is.na(dat$missingness) & dat$missingness == 1L
    cu <- !is.na(dat$missingness) & dat$missingness == 2L
    expect_true(any(cl) && any(cu))
    expect_true(all(is.na(dat$R[cl | cu])) && all(is.na(dat$RR[cl | cu])) &&
                  all(is.na(dat$rt[cl | cu])))
    obs <- is.na(dat$missingness)
    expect_true(all(dat$rt[obs] >= LC & dat$rt[obs] <= UC))
    expect_true(all(dat$RR[obs] %in% seq_len(K)))
    # censored counts agree with the model's window probabilities
    du <- EMC2:::design_model(dat, des, compress = FALSE, verbose = FALSE, rt_resolution = 1e-9)
    pm <- t(as.matrix(p))
    for (m in 1:2) {
      q <- du
      q$missingness <- m
      q$R[] <- NA
      pr <- exp(as.vector(EMC2:::calc_ll_manager(pm, q, des$model, return_trialwise = TRUE)))
      n_obs <- sum(dat$missingness == m, na.rm = TRUE)
      expect_lt(abs(n_obs - sum(pr)) / sqrt(sum(pr * (1 - pr))), 4)
    }
    # the censored data go through make_emc and both likelihoods agree
    emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1, verbose = FALSE))
    dadm <- emc[[1]]$data[[1]]
    expect_true(any(!is.na(dadm$missingness)))
    llc <- EMC2:::calc_ll_manager(pm, dadm, des$model)
    llr <- EMC2:::calc_ll_manager(pm, dadm, mtlba_r_model(des$model))
    expect_true(is.finite(llc))
    expect_lt(mtlba_err(llc, llr), 1e-12)
    # truncation: nothing outside the window survives
    datt <- make_data(p, des, n_trials = 300, TC = list(LT = .45, UT = 1.5))
    expect_true(all(datt$rt >= .45 & datt$rt <= 1.5))
    dt <- EMC2:::design_model(datt, des, verbose = FALSE)
    expect_equal(as.vector(EMC2:::calc_ll_manager(pm, dt, des$model)),
                 EMC2:::calc_ll_manager(pm, dt, mtlba_r_model(des$model)), tolerance = 1e-12)
  }
})

test_that("a don't-know data set with the DK paper's lower cut (.4 s) evaluates", {
  # the DK paper (Experiment 2) dropped responses faster than .4 s: fitted here
  # either as lower censoring or as lower truncation
  des <- mtlba_design(2, "BE", TRUE, rlev = c("new", "old"))
  p <- c(v = 1, v_lMTRUE = 1.2, B = log(1.2), A = log(.6), t0 = log(.2),
         c1 = log(.5), c1_lRold = .1)
  set.seed(95)
  for (TC in list(list(LC = .4), list(LT = .4))) {
    dat <- make_data(p, des, n_trials = 500, TC = TC)
    expect_true(any(dat$RR == 1, na.rm = TRUE))
    expect_true(all(dat$R[!is.na(dat$RR) & dat$RR == 1] == "new"))
    dadm <- EMC2:::design_model(dat, des, verbose = FALSE)
    pm <- t(as.matrix(p))
    llc <- EMC2:::calc_ll_manager(pm, dadm, des$model)
    expect_true(is.finite(llc))
    expect_equal(as.vector(llc), EMC2:::calc_ll_manager(pm, dadm, mtlba_r_model(des$model)),
                 tolerance = 1e-12)
  }
})

test_that("C++ equals the R reference with censoring and truncation (OpenMP = serial)", {
  for (K in 1:3) for (rule in c("BE", "TC")) for (dk in c(FALSE, if (K > 1) TRUE)) {
    crit <- if (K > 1) paste0("c", seq_len(K - 1)) else character(0)
    form <- c(list(v ~ lM, sv ~ lM, B ~ lR, A ~ 1, t0 ~ 1),
              lapply(crit, function(cn) stats::as.formula(paste(cn, "~ lR"))))
    mod <- local({ k <- K; r <- rule; d <- dk
      function() MTLBA(n_ratings = k, dont_know = d, rule = r) })
    des <- suppressMessages(design(
      factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
      model = mod, matchfun = function(d) d$S == d$lR, formula = form, report_p_vector = FALSE))
    p <- sampled_pars(des, doMap = FALSE)
    p[] <- 0
    p[c("v", "v_lMTRUE", "sv", "sv_lMTRUE", "B", "B_lRright", "A", "t0")] <-
      c(1, 1.5, 0, log(.7), log(1), .2, log(.5), log(.25))
    for (cn in crit) p[c(cn, paste0(cn, "_lRright"))] <- c(log(.6), .15)
    set.seed(200 + 10 * K + dk)
    # mixed: lower and upper censoring, an unknown-rt trial (both), truncation
    dat <- make_data(p, des, n_trials = 150, TC = list(LT = .3, LC = .55, UC = 1.2, UT = 3))
    dat$missingness[1:2] <- 3L; dat$R[1:2] <- NA; dat$RR[1:2] <- NA; dat$rt[1:2] <- NA
    expect_true(all(c(1L, 2L, 3L) %in% dat$missingness))
    props <- rbind(p, mvtnorm::rmvnorm(30, p, diag(length(p)) * .3))
    colnames(props) <- names(p)
    ext <- matrix(p, nrow = 7, ncol = length(p), byrow = TRUE, dimnames = list(NULL, names(p)))
    ext[1, "sv"] <- log(.02)
    ext[2, "v"] <- -4
    ext[3, "v"] <- 12
    ext[4, "A"] <- -Inf
    ext[5, "t0"] <- log(.6)               # t0 above LT and LC: empty lower windows
    ext[6, "B"] <- 800
    ext[7, "v"] <- NaN
    if (K > 1) ext[1, "c1"] <- log(8)
    props <- rbind(props, ext)
    for (compress in c(TRUE, FALSE)) {
      dadm <- EMC2:::design_model(dat, des, compress = compress, verbose = FALSE)
      llr <- EMC2:::calc_ll_manager(props, dadm, mtlba_r_model(des$model))
      llc <- EMC2:::calc_ll_manager(props, dadm, des$model)
      expect_lt(mtlba_err(llc, llr), 1e-10)
      twr <- apply(props, 1, function(pv)
        EMC2:::trial_ll_mtlba(EMC2:::get_pars_matrix_oo(pv, dadm, des$model()), dadm, K, rule, dk))
      twc <- EMC2:::calc_ll_manager(props, dadm, des$model, return_trialwise = TRUE)
      expect_lt(mtlba_err(as.vector(twc), as.vector(twr)), 1e-10)
      expect_true(all(is.finite(llc)))
      expect_identical(mtlba_omp(EMC2:::calc_ll_manager(props, dadm, des$model)), llc)
      expect_identical(mtlba_omp(EMC2:::calc_ll_manager(props, dadm, des$model,
                                                        return_trialwise = TRUE)), twc)
    }
  }
})

test_that("censoring refusals: retained response, unequal t0, contamination", {
  des <- mtlba_design(2, "BE")
  p <- mtlba_p(2)
  set.seed(97)
  expect_error(make_emc(make_data(p, des, n_trials = 100, TC = list(LC = .5, LCresponse = TRUE)),
                        des, type = "single", verbose = FALSE), "response retained")
  # a missing response that is not censored
  dc <- make_data(p, des, n_trials = 50)
  dc$R[1] <- NA; dc$RR[1] <- NA
  expect_error(EMC2:::design_model(dc, des, verbose = FALSE), "contamination")
  # t0 by accumulator is fine without censoring but refused with it (R and C++)
  dt0 <- suppressMessages(design(
    factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
    model = function() MTLBA(n_ratings = 2), matchfun = function(d) d$S == d$lR,
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ lR, c1 ~ lR), report_p_vector = FALSE))
  p0 <- c(p, t0_lRright = .1)
  dat <- make_data(p0, dt0, n_trials = 50)
  pm <- t(as.matrix(p0))
  expect_true(is.finite(EMC2:::calc_ll_manager(pm, EMC2:::design_model(dat, dt0, verbose = FALSE),
                                               dt0$model)))
  datc <- make_data(p0, dt0, n_trials = 50, TC = list(UC = 1))
  dadm <- EMC2:::design_model(datc, dt0, verbose = FALSE)
  expect_error(EMC2:::calc_ll_manager(pm, dadm, dt0$model), "t0 equal")
  expect_error(mtlba_omp(EMC2:::calc_ll_manager(pm, dadm, dt0$model)), "t0 equal")
  expect_error(EMC2:::calc_ll_manager(pm, dadm, mtlba_r_model(dt0$model)), "t0 equal")
})

test_that("truncation normaliser cache: many distinct parameter rows (covariate)", {
  # a trial-level covariate gives every trial its own parameters: more keys
  # than the per-call cache holds, so eviction and misses are exercised
  for (rule in c("BE", "TC")) {
    des <- suppressMessages(design(
      factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
      covariates = "x", model = local({ r <- rule; function() MTLBA(n_ratings = 3, rule = r) }),
      matchfun = function(d) d$S == d$lR,
      formula = list(v ~ lM + x, B ~ 1, A ~ 1, t0 ~ 1, c1 ~ lR, c2 ~ lR), report_p_vector = FALSE))
    p <- c(v = 1, v_lMTRUE = 1.5, v_x = .3, B = log(1), A = log(.5), t0 = log(.25),
           c1 = log(.6), c1_lRright = .15, c2 = log(.6), c2_lRright = .15)
    set.seed(98)
    dd <- data.frame(subjects = factor(1), S = factor(rep(c("left", "right"), 150)),
                     x = round(stats::rnorm(300), 2), R = factor("left", levels = c("left", "right")),
                     rt = 1, RR = 1)
    dat <- make_data(p, des, data = dd, TC = list(LT = .35, LC = .5, UC = 1.3, UT = 2.5))
    dadm <- EMC2:::design_model(dat, des, verbose = FALSE)
    pm <- rbind(p, p + .05)
    llc <- EMC2:::calc_ll_manager(pm, dadm, des$model, return_trialwise = TRUE)
    llr <- apply(pm, 1, function(pv)
      EMC2:::trial_ll_mtlba(EMC2:::get_pars_matrix_oo(pv, dadm, des$model()), dadm, 3, rule, FALSE))
    expect_lt(mtlba_err(as.vector(llc), as.vector(llr)), 1e-10)
    expect_identical(mtlba_omp(EMC2:::calc_ll_manager(pm, dadm, des$model, return_trialwise = TRUE)),
                     llc)
  }
})

# ---- correlated drift rates (rho != 0): R reference (Stage M4) ----

test_that("bivariate normal upper orthant (Genz) matches mvtnorm", {
  set.seed(90)
  n <- 400
  h <- rnorm(n, 0, 2.5); k <- rnorm(n, 0, 2.5)
  r <- c(runif(n / 2, -1, 1), sample(c(-1, 1), n / 2, TRUE) * (1 - 10^-runif(n / 2, 1, 6)))
  ours <- EMC2:::mtlba_bvu(h, k, r)
  mv <- sapply(seq_len(n), function(i)
    mvtnorm::pmvnorm(lower = c(h[i], k[i]), upper = c(Inf, Inf),
                     corr = matrix(c(1, r[i], r[i], 1), 2))[1])
  expect_lt(max(abs(ours - mv)), 1e-14)
  expect_equal(EMC2:::mtlba_bvu(c(.3, -1), c(-.2, .4), c(0, 0)),
               pnorm(-c(.3, -1)) * pnorm(-c(-.2, .4)), tolerance = 1e-15)
  expect_true(all(is.nan(EMC2:::mtlba_bvu(c(0, 0, NA), c(0, 0, 0), c(1, -1, .2)))))
  # Gauss-Legendre rule: exact for polynomials of degree 2n - 1
  g <- EMC2:::mtlba_gl(16)
  expect_equal(sum(g$w), 2, tolerance = 1e-15)
  expect_equal(sum(g$w * g$x^30), 2 / 31, tolerance = 1e-14)
})

test_that("correlated density: quadrature = BETE closed form = brute-force 2-D (BE and TC)", {
  g <- golden_mtlba("ref_density.rds")
  cells <- g$cells[g$cells$rho != 0, ]
  cells$id <- seq_len(nrow(cells))
  ref <- def <- numeric(nrow(cells))
  grp <- split(cells, list(cells$set, cells$K, cells$rule, cells$rho, cells$R), drop = TRUE)
  for (x in grp) {
    st <- bete_set(g$sets[[x$set[1]]], x$K[1])
    st$rho <- x$rho[1]
    rr <- bete_resp(x$R[1], x$K[1])
    ref[x$id] <- mtlba_f(rr[1], rr[2], x$K[1], x$rule[1], FALSE, st,
                         EMC2:::mtlba_quad_reference())(x$rt)
    def[x$id] <- mtlba_f(rr[1], rr[2], x$K[1], x$rule[1], FALSE, st)(x$rt)
  }
  expect_true(all(is.finite(ref)) && all(ref >= 0))
  # the three evaluations (reference quadrature rule)
  expect_lt(max(abs(ref - cells$dens_2d)), 1e-13)
  expect_lt(max(abs(ref - cells$dens_1d)), 1e-13)
  expect_lt(max(abs(ref - cells$dens_bete)), 1e-13)
  body <- cells$dens_2d > 1e-8
  expect_gt(sum(body), 3500)
  expect_lt(max(abs(ref[body] / cells$dens_2d[body] - 1)), 1e-10)
  # the fit rule: within the 1e-6 log-likelihood target of the reference
  expect_lt(max(abs(log(def[body]) - log(ref[body]))), 1e-6)
  expect_lt(max(abs(def - ref)), 1e-9)
})

mtlba_rho_sets <- function() list(
  list(v = c(2.2, 1.3), sv = c(.8, 1.1), A = c(.5, .7), b = c(1.4, 1.6), t0 = c(.2, .2),
       d = rbind(c(.3, .7), c(.2, .5)), rho = -.8),
  list(v = c(1.5, 1.8), sv = c(1.2, .9), A = c(.6, .3), b = c(1.2, 1.5), t0 = c(.25, .25),
       d = rbind(c(.25, .6), c(.4, .8)), rho = .6))

test_that("correlated density joins the closed form at rho = 0 and integrates to 1", {
  rt <- .2 + c(.1, .2, .35, .6, 1, 2)
  for (st in mtlba_rho_sets()) for (rule in c("BE", "TC")) {
    s0 <- st; s0$rho <- 0
    for (i in 1:6) {
      rr <- mtlba_responses(3, FALSE)[i, ]
      f0 <- mtlba_f(rr[1], rr[2], 3, rule, FALSE, s0)(rt)
      fe <- sapply(c(-1e-4, 1e-4), function(e) {
        s1 <- st; s1$rho <- e
        mtlba_f(rr[1], rr[2], 3, rule, FALSE, s1)(rt)
      })
      # even part in rho is O(rho^2)
      expect_lt(max(abs(rowMeans(fe) - f0)), 1e-6)
    }
    for (dk in c(FALSE, TRUE)) {
      rsp <- mtlba_responses(3, dk)
      tot <- sum(apply(rsp, 1, function(r) mtlba_prob(r[1], r[2], 3, rule, dk, st)))
      expect_equal(tot, 1, tolerance = 1e-7)
    }
    # don't know = the sum of the two rating-1 kernels
    expect_equal(mtlba_f(1, 1, 3, rule, TRUE, st)(rt),
                 mtlba_f(1, 1, 3, rule, FALSE, st)(rt) + mtlba_f(2, 1, 3, rule, FALSE, st)(rt),
                 tolerance = 1e-12)
  }
  # rho outside (-1, 1) or non-finite: -Inf, not an error
  st <- mtlba_rho_sets()[[1]]
  for (r in c(1, -1, NA, Inf)) {
    st$rho <- r
    expect_identical(suppressWarnings(mtlba_f(1, 2, 3, "BE", FALSE, st)(.6)), 0)
  }
  # A = 0 on the stopping accumulator: no integral left
  st <- mtlba_rho_sets()[[2]]; st$A <- c(0, .3); s2 <- st; s2$A <- c(1e-9, .3)
  expect_equal(mtlba_f(1, 2, 3, "BE", FALSE, st)(rt), mtlba_f(1, 2, 3, "BE", FALSE, s2)(rt),
               tolerance = 1e-6)
})

test_that("correlated censoring and truncation: censored + observed mass is 1", {
  LC <- .45; UC <- 1.6; LT <- .3; UT <- 2.5
  for (rho in c(-.8, .6)) for (rule in c("BE", "TC")) for (K in c(1, 3)) for (dk in c(FALSE, if (K > 1) TRUE)) {
    p <- cbind(mtlba_trial_pars(K), rho = rho)
    rsp <- mtlba_responses(K, dk)
    ll_miss <- function(rt, R, RR, bnd)
      EMC2:::mtlba_logdens_missing(rt, R, RR, p, K, rule, dk, bnd)
    dens <- function(rt, R, RR, bnd = mtlba_bnd()) {
      n <- length(rt)
      exp(EMC2:::mtlba_logdens_missing(rt, rep(R, n), rep(RR, n), p[rep(1:2, n), ], K, rule, dk,
                                       lapply(bnd, rep_len, n)))
    }
    obs <- sum(mapply(function(R, RR)
      stats::integrate(dens, LC, UC, R = R, RR = RR, rel.tol = 1e-10)$value, rsp[, "R"], rsp[, "RR"]))
    cens <- function(m, bnd = mtlba_bnd(m, LC = LC, UC = UC)) exp(ll_miss(NA, NA, NA, bnd))
    expect_equal(obs + cens(1L) + cens(2L), 1, tolerance = 1e-8)
    expect_equal(cens(3L), cens(1L) + cens(2L), tolerance = 1e-12)
    tb <- mtlba_bnd(LT = LT, UT = UT)
    tr <- sum(mapply(function(R, RR)
      stats::integrate(dens, LT, UT, R = R, RR = RR, bnd = tb, rel.tol = 1e-10)$value,
      rsp[, "R"], rsp[, "RR"]))
    expect_equal(tr, 1, tolerance = 1e-8)
    # stop terms: P + S = 1, and S = 1 - the integrated response density
    Tm <- EMC2:::mtlba_thresholds(p, K)
    st <- EMC2:::mtlba_logstop_any(c(.7, Inf, 0), p[rep(1:2, 3), ], Tm[rep(1:2, 3), ],
                                   c(1, 3, 5), c(2, 4, 6), K, rule)
    expect_equal(exp(st$P[1]) + exp(st$S[1]), 1, tolerance = 1e-12)
    expect_identical(c(st$P[2:3], st$S[2:3]), c(0, -Inf, -Inf, 0))
    below <- sum(mapply(function(R, RR)
      stats::integrate(dens, .2, .9, R = R, RR = RR, rel.tol = 1e-10)$value, rsp[, "R"], rsp[, "RR"]))
    expect_equal(exp(st$P[1]), below, tolerance = 1e-8)
  }
  # rho -> 0: the correlated stop terms join the closed form
  for (rule in c("BE", "TC")) {
    p <- mtlba_trial_pars(3)
    Tm <- EMC2:::mtlba_thresholds(p, 3)
    cf <- EMC2:::mtlba_logstop(.7, p, Tm, 1, 2, 3, rule)
    q <- EMC2:::mtlba_stop_rho1(.7, p[1, "A"], p[1, "v"], p[1, "sv"], Tm[1, -1], p[2, "A"], p[2, "v"],
                                p[2, "sv"], Tm[2, -1], 1e-9, 3, rule)
    expect_equal(q, c(cf$P, cf$S), tolerance = 1e-7)
  }
})

test_that("rMTLBA agrees with the correlated density (rho = -.8, -.4, .4, .8)", {
  st <- list(v = c(2.2, 1.3), sv = c(.8, 1.1), A = c(.5, .7), b = c(1.4, 1.6), t0 = c(.2, .25),
             d = rbind(c(.3, .7), c(.2, .5)))
  n <- 1e5
  qp <- c(.1, .3, .5, .7, .9)
  lR <- factor(rep(c("a", "b"), n))
  for (rho in c(-.8, -.4, .4, .8)) for (rule in c("BE", "TC")) {
    dk <- rule == "TC" && rho > 0
    st$rho <- rho
    set.seed(round(100 + 10 * rho) + (rule == "TC"))
    sim <- EMC2:::rMTLBA(lR, mtlba_pars(n, st$v, st$sv, st$A, st$b, st$t0, st$d, rho), 3, rule, dk)
    expect_true(all(sim$RR %in% 1:3) && !anyNA(sim$R) && all(is.finite(sim$rt)))
    rsp <- mtlba_responses(3, dk)
    for (i in seq_len(nrow(rsp))) {
      sel <- sim$RR == rsp[i, 2] & (if (dk && rsp[i, 2] == 1) TRUE else as.integer(sim$R) == rsp[i, 1])
      pe <- mtlba_prob(rsp[i, 1], rsp[i, 2], 3, rule, dk, st)
      expect_lt(abs(mean(sel) - pe), 4.5 * sqrt(pe * (1 - pe) / n))
      if (pe * n > 5000) {
        Fq <- sapply(stats::quantile(sim$rt[sel], qp), function(u)
          mtlba_prob(rsp[i, 1], rsp[i, 2], 3, rule, dk, st, upper = u)) / pe
        expect_lt(max(abs(Fq - qp)), 4.5 * .5 / sqrt(pe * n))
      }
    }
  }
  # the drift sampler: both positive, the truncated bivariate normal's moments
  set.seed(95)
  m <- 2e5
  dr <- EMC2:::mtlba_rdrift_rho(rep(.5, m), rep(1, m), rep(-.2, m), rep(1.5, m), rep(-.7, m))
  expect_true(all(dr$v1 > 0 & dr$v2 > 0))
  S <- matrix(c(1, -.7 * 1.5, -.7 * 1.5, 1.5^2), 2)
  mom <- EMC2:::mtlba_bvu(-.5, .2 / 1.5, -.7)
  e1 <- stats::integrate(function(x) x * dnorm(x, .5, 1) *
                           pnorm((-.2 - .7 * 1.5 * (x - .5)) / (1.5 * sqrt(1 - .49))), 0, Inf)$value / mom
  expect_lt(abs(mean(dr$v1) - e1), 4.5 * sd(dr$v1) / sqrt(m))
})

test_that("correlated model: design, simulation and the R likelihood through EMC2", {
  for (rule in c("BE", "TC")) for (dk in c(FALSE, TRUE)) {
    K <- if (dk) 2 else 3
    crit <- paste0("c", seq_len(K - 1))
    mod <- local({ k <- K; r <- rule; d <- dk
      function() { m <- MTLBA(n_ratings = k, dont_know = d, rule = r); m$c_name <- NULL; m } })
    des <- suppressMessages(design(
      factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
      model = mod, matchfun = function(d) d$S == d$lR,
      formula = c(list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, rho ~ S),
                  lapply(crit, function(cn) stats::as.formula(paste(cn, "~ lR")))),
      report_p_vector = FALSE))
    p <- c(mtlba_p(K), rho = qnorm((1 - .6) / 2), rho_Sright = qnorm((1 + .3) / 2) - qnorm((1 - .6) / 2))
    p <- p[names(sampled_pars(des))]
    mp <- mapped_pars(des, p)
    expect_equal(sort(unique(round(mp$rho, 10))), c(-.6, .3))
    set.seed(96)
    dat <- make_data(p, des, n_trials = 150)
    expect_true(all(is.finite(dat$rt)) && all(dat$RR %in% seq_len(K)))
    res <- .05
    dc <- EMC2:::design_model(dat, des, compress = TRUE, verbose = FALSE, rt_resolution = res)
    du <- EMC2:::design_model(dat, des, compress = FALSE, verbose = FALSE, rt_resolution = res)
    pm <- t(as.matrix(p))
    llc <- EMC2:::calc_ll_manager(pm, dc, des$model)
    expect_equal(llc, EMC2:::calc_ll_manager(pm, du, des$model), tolerance = 1e-12)
    # direct, natural-scale parameters by hand
    n <- nrow(dat)
    m1 <- dat$S == "left"
    C <- matrix(sapply(seq_len(K - 1), function(k) exp(p[[paste0("c", k)]] + c(0, .15))), nrow = 2)
    if (K > 2) for (k in 2:(K - 1)) C[, k] <- C[, k - 1] + C[, k]
    pars <- mtlba_pars(n, 0, 1, .5, 1.5, .25, exp(-C[, rev(seq_len(K - 1)), drop = FALSE]))
    pars[, "v"] <- as.vector(rbind(ifelse(m1, 2.5, 1), ifelse(m1, 1, 2.5)))
    pars <- cbind(pars, rho = rep(ifelse(m1, -.6, .3), each = 2))
    rt <- floor(dat$rt / res) * res
    direct <- sum(pmax(log(1e-10), EMC2:::dMTLBA(rt, as.integer(dat$R), dat$RR, pars, K, rule, dk,
                                                  log = TRUE)))
    expect_equal(as.vector(llc), direct, tolerance = 1e-10)
    # censoring and truncation evaluate, and differ from the uncorrelated value
    cd <- make_data(p, des, n_trials = 150, TC = list(LC = .5, UC = 1.8, LT = .3, UT = 4))
    dm <- EMC2:::design_model(cd, des, verbose = FALSE)
    ll <- EMC2:::calc_ll_manager(pm, dm, des$model)
    p0 <- p; p0[c("rho", "rho_Sright")] <- 0
    expect_true(is.finite(ll) && abs(ll - EMC2:::calc_ll_manager(t(as.matrix(p0)), dm, des$model)) > 1e-3)
  }
  # the compiled likelihood does not have rho yet (Stage M5): it stops
  des <- suppressMessages(design(
    factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
    model = function() MTLBA(2), matchfun = function(d) d$S == d$lR,
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, rho ~ 1, c1 ~ 1), report_p_vector = FALSE))
  p <- c(v = 1, v_lMTRUE = 1.5, B = 0, A = log(.5), t0 = log(.25), rho = .4, c1 = log(.6))
  set.seed(97)
  dat <- make_data(p, des, n_trials = 20)
  dadm <- EMC2:::design_model(dat, des, verbose = FALSE)
  expect_error(EMC2:::calc_ll_manager(t(as.matrix(p)), dadm, des$model), "rho != 0")
})
