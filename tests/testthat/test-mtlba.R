# MTLBA (multiple-threshold LBA, uncorrelated: rho = 0) and the rating layer's
# don't-know response. Ground truth: golden/mtlba (Stage M0) -- the don't-know
# paper's n1PDF.DK, BETE's n1PDF_MTR_1 and two independent integrations of the
# kernel, plus 1e6-trial simulator summaries. Closed forms there are accurate
# absolutely, not relatively, in the far tail (density < 1e-8), so parity is
# absolute everywhere and relative only in the body.

golden_mtlba <- function(f) readRDS(test_path("golden", "mtlba", f))

# Two rows per trial: v, sv, A, b, t0 and threshold proportions d1 .. d{K-1}
# (d: 2 x (K - 1) matrix, row = accumulator).
mtlba_pars <- function(n, v, sv, A, b, t0, d = NULL) {
  r <- function(x) rep_len(x, 2 * n)
  p <- cbind(v = r(v), sv = r(sv), A = r(A), b = r(b), t0 = r(t0))
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

# set list(v, sv, A, b, t0, d) -> density function of rt for one response
mtlba_f <- function(R, RR, K, rule, dk, st)
  function(rt) EMC2:::dMTLBA(rt, rep(R, length(rt)), rep(RR, length(rt)),
                             mtlba_pars(length(rt), st$v, st$sv, st$A, st$b, st$t0, st$d),
                             K, rule, dk)

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
           formula = list(v ~ 1, rho ~ 1), report_p_vector = FALSE)), "not implemented")
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
