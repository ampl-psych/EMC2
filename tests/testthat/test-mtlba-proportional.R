# MTLBA(parameterization = "proportional"): the source paper's nested parameterisation
# (B = log b; c1 = the top confidence threshold as a probit proportion of b; A = the
# start-point range as a proportion of that threshold; each lower threshold a proportion
# of the one above, nested above A or, when the lowest threshold is free, on T_2 alone).
# It is implemented as an exact translation onto the additive model, so the checks are:
# the hand formulas, proportional = additive at the translated parameters (R and C++,
# every option), C++ = R, simulation, fitting and mapping, and the DMC code's density
# reached from DMC's sampled values (golden/mtlba/dmc_p3_density.rds).

prop_golden <- function(f) readRDS(test_path("golden", "mtlba", f))
prop_omp <- function(expr, threads = 2) {
  old <- options(emc.ll_backend = "multithreaded", emc.n_threads = threads)
  on.exit(options(old))
  expr
}
# R reference trial log-likelihoods (trial_ll_mtlba through the model's Ttransform),
# particles in columns as calc_ll_manager(return_trialwise = TRUE)
prop_tw_r <- function(props, dadm, des) {
  m <- des$model()
  apply(props, 1, function(pv)
    EMC2:::trial_ll_mtlba(EMC2:::get_pars_matrix_oo(pv, dadm, m), dadm, m$n_ratings, m$rule, m$dont_know,
                          guess = m$guess, gb_constraint = m$gb_constraint, drifts = m$drifts))
}
# error relative to max(1, |reference|), elementwise
prop_err <- function(a, b) abs(a - b) / pmax(abs(b), 1)

# Design with every threshold quantity by accumulator (treatment-coded lR, two cells), so
# a proportional parameter vector translates into an additive one cell by cell.
prop_design <- function(K, rule = "BE", lowest = "auto", dk = FALSE, rho = FALSE, guess = "none",
                        drifts = "any", param = "proportional") {
  crit <- if (K > 1) paste0("c", seq_len(K - 1)) else character(0)
  pg <- if (guess != "none" && K > 1) paste0("pg", seq_len(K - 1)) else character(0)
  form <- c(list(v ~ lM, sv ~ lM, B ~ lR, A ~ lR, t0 ~ 1), if (rho) list(rho ~ 1),
            if (guess == "accumulator") list(gp ~ lR, gB ~ lR),
            if (guess == "trial") list(gp ~ 1, gB ~ lR),
            lapply(crit, function(cn) stats::as.formula(paste(cn, "~ lR"))),
            lapply(pg, function(cn) stats::as.formula(paste(cn, "~ 1"))))
  mod <- local({
    k <- K; r <- rule; lw <- lowest; d <- dk; dr <- drifts; pa <- param
    g <- if (guess == "none") "accumulator" else guess
    function() MTLBA(n_ratings = k, dont_know = d, rule = r, lowest_threshold = lw, drifts = dr,
                     guess = g, parameterization = pa)
  })
  suppressMessages(design(
    factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
    model = mod, matchfun = function(d) d$S == d$lR, formula = form, report_p_vector = FALSE))
}

# A proportional parameter vector for prop_design (drift and guessing values as the M7c tests).
prop_p <- function(des, K, rho = 0, guess = "none") {
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- 0
  p[c("v", "v_lMTRUE", "sv", "sv_lMTRUE", "B", "B_lRright", "A", "A_lRright", "t0")] <-
    c(1, 1.5, 0, log(.7), log(1.5), .2, qnorm(.6), .2, log(.25))
  if (rho != 0) p["rho"] <- qnorm((1 + rho) / 2)
  if (K > 1) {
    p[c("c1", "c1_lRright")] <- c(qnorm(.7), .15)
    for (k in seq_len(K - 2) + 1) p[paste0("c", k, c("", "_lRright"))] <- c(qnorm(.45), .15)
  }
  if (guess != "none") {
    p[c("gp", "gB", "gB_lRright")] <- c(qnorm(.25), qnorm(.5), .3)
    if (guess == "accumulator") p["gp_lRright"] <- -.4
    for (k in seq_len(K - 1)) p[paste0("pg", k)] <- qnorm(.45)
  }
  p
}

# Natural rows of a proportional model (get_pars_matrix_oo, two rows per trial) -> the
# additive model's sampled A, B and criteria per accumulator, from the natural thresholds:
# A, B = b - A, c1 = -log d_{K-1} and c_k = log(d_{K-k+1} / d_{K-k}) above A, the free
# lowest threshold c_{K-1} = log(T_2 / T_1) (K = 2: log(b / T_1)).
prop_to_additive <- function(pm, K, free) {
  A <- pm[, "A"]; b <- pm[, "b"]
  out <- cbind(A = log(A), B = log(b - A))
  if (K > 1) {
    Tm <- EMC2:::mtlba_thresholds(pm, K)   # columns T_0 .. T_K
    d <- (Tm[, 2:K, drop = FALSE] - A) / (b - A)
    cc <- sapply(seq_len(K - 1), function(k) {
      j <- K - k
      if (k == 1) { if (K == 2 && free) log(b / Tm[, 2]) else -log(d[, j]) } else
        if (j == 1 && free) log(Tm[, 3] / Tm[, 2]) else log(d[, j + 1] / d[, j])
    })
    cc <- matrix(cc, nrow(pm))
    colnames(cc) <- paste0("c", seq_len(K - 1))
    out <- cbind(out, log(cc))
  }
  out
}

# The additive parameter vector at which the additive model is the proportional one at p:
# A, B and the criteria (all ~ lR) from the first trial's two rows, everything else copied.
prop_translate <- function(p, des_p, des_a, dadm_p, K, free) {
  pm <- EMC2:::get_pars_matrix_oo(p, dadm_p, des_p$model())
  lv <- as.integer(dadm_p$lR)
  rows <- c(which(lv == 1)[1], which(lv == 2)[1])
  tr <- prop_to_additive(pm[rows, , drop = FALSE], K, free)
  q <- sampled_pars(des_a, doMap = FALSE)
  q[] <- p[names(q)]
  for (nm in colnames(tr)) q[c(nm, paste0(nm, "_lRright"))] <- c(tr[1, nm], tr[2, nm] - tr[1, nm])
  q
}

test_that("proportional: factory (p_types, transforms, bounds, c_options); additive unchanged", {
  for (K in 1:4) {
    crit <- if (K > 1) paste0("c", seq_len(K - 1)) else character(0)
    a <- MTLBA(K)
    m <- MTLBA(K, parameterization = "proportional")
    expect_identical(names(m$p_types), names(a$p_types))
    expect_identical(m$c_options[["proportional"]], 1)
    expect_identical(a$c_options[["proportional"]], 0)
    expect_identical(a$parameterization, "additive")
    expect_equal(unname(m$p_types[c("A", crit)]), rep(0, K))
    expect_identical(unname(m$transform$func[c("B", "A", crit)]), c("exp", "pnorm", rep("pnorm", K - 1)))
    expect_identical(unname(a$transform$func[c("B", "A", crit)]), rep("exp", K + 1))
    expect_equal(unname(m$bound$minmax[, c("A", crit), drop = FALSE]), matrix(c(0, Inf, rep(c(0, 1), K - 1)), 2))
    expect_equal(unname(a$bound$minmax[, "A"]), c(1e-4, Inf))
    expect_identical(m$recalculated_p_types, "A")
    expect_null(a$recalculated_p_types)
  }
  expect_error(MTLBA(3, parameterization = "nested"))
})

test_that("proportional: thresholds and A follow the hand formulas (K = 1 .. 4, BE free / TC above A)", {
  for (K in 1:4) for (rule in c("BE", "TC")) {
    des <- prop_design(K, rule)
    set.seed(30 + K)
    p <- prop_p(des, K) + stats::rnorm(length(prop_p(des, K)), 0, .3)
    dat <- make_data(p, des, n_trials = 4)
    dadm <- suppressWarnings(EMC2:::design_model(dat, des, compress = FALSE, verbose = FALSE))
    pm <- EMC2:::get_pars_matrix_oo(p, dadm, des$model())
    right <- as.numeric(dadm$lR == "right")
    cell <- function(nm) p[[nm]] + right * p[[paste0(nm, "_lRright")]]
    b <- exp(cell("B"))
    expect_equal(pm[, "b"], b, tolerance = 1e-14)
    expect_equal(pm[, "B"], b, tolerance = 1e-14)   # B reports b itself
    if (K == 1) {
      expect_equal(pm[, "A"], b * pnorm(cell("A")), tolerance = 1e-14)
      next
    }
    Tm <- EMC2:::mtlba_thresholds(pm, K)
    Th <- matrix(NA, nrow(pm), K)
    Th[, K] <- b
    Th[, K - 1] <- b * pnorm(cell("c1"))
    A <- Th[, K - 1] * pnorm(cell("A"))
    expect_equal(pm[, "A"], A, tolerance = 1e-14)
    if (K > 2) for (j in (K - 2):1) {
      pj <- pnorm(cell(paste0("c", K - j)))
      Th[, j] <- if (j == 1 && rule == "BE") Th[, 2] * pj else A + (Th[, j + 1] - A) * pj
    }
    expect_equal(Tm[, 2:(K + 1)], Th, tolerance = 1e-14)
    # every threshold but a free T_1 above A, and the ordering
    expect_true(all(Th[, K - 1] > A))
    if (K > 2) expect_true(all(t(apply(Th, 1, diff)) > 0))
    if (rule == "TC") expect_true(all(Th[, 1] > A))
  }
})

test_that("proportional = additive at the translated parameters, in R and C++; C++ = R; OpenMP identical", {
  # (the full grid, K 1-4 x BE/TC x rho 0/-.5/.8 x guessing none/accumulator/trial x
  # censoring+truncation x don't know with 12 vectors each, is run by the stage's check
  # script; this is a thinned version of it)
  cfg <- expand.grid(K = 1:4, rule = c("BE", "TC", "BEa"), rho = c(0, -.5, .8),
                     guess = c("none", "accumulator", "trial"), cens = c(FALSE, TRUE),
                     dk = c(FALSE, TRUE), stringsAsFactors = FALSE)
  cfg <- cfg[!(cfg$dk & cfg$K == 1), ]
  set.seed(4401)
  cfg <- cfg[sort(sample(nrow(cfg), 16)), ]
  for (i in seq_len(nrow(cfg))) {
    x <- cfg[i, ]
    K <- x$K; rule <- if (x$rule == "TC") "TC" else "BE"
    lowest <- if (x$rule == "BEa") "above_A" else "auto"
    free <- rule == "BE" && lowest == "auto" && K > 1
    drifts <- if (i %% 2) "any" else "positive"
    dp <- prop_design(K, rule, lowest, x$dk, x$rho != 0, x$guess, drifts)
    da <- prop_design(K, rule, lowest, x$dk, x$rho != 0, x$guess, drifts, param = "additive")
    p <- prop_p(dp, K, x$rho, x$guess)
    set.seed(4500 + i)
    dat <- make_data(p, dp, n_trials = 40)
    if (x$cens) {
      dat <- make_missing(dat, LC = .45, UC = 1.6, LT = .3, UT = 4)
      dat$missingness[1:2] <- 3L; dat$R[1:2] <- NA; dat$RR[1:2] <- NA; dat$rt[1:2] <- NA
    }
    props <- rbind(p, mvtnorm::rmvnorm(3, p, diag(length(p)) * .1))
    colnames(props) <- names(p)
    dadm_p <- EMC2:::design_model(dat, dp, compress = !x$cens, verbose = FALSE, rt_resolution = 1e-6)
    dadm_a <- EMC2:::design_model(dat, da, compress = !x$cens, verbose = FALSE, rt_resolution = 1e-6)
    qa <- t(apply(props, 1, prop_translate, des_p = dp, des_a = da, dadm_p = dadm_p, K = K, free = free))
    twc <- EMC2:::calc_ll_manager(props, dadm_p, dp$model, return_trialwise = TRUE)
    twr <- prop_tw_r(props, dadm_p, dp)
    tac <- EMC2:::calc_ll_manager(qa, dadm_a, da$model, return_trialwise = TRUE)
    tar <- prop_tw_r(qa, dadm_a, da)
    expect_true(all(is.finite(twc)) && all(twc > log(1e-10) * 1.0001))
    expect_lt(max(prop_err(twc, tac)), 1e-13)
    expect_lt(max(prop_err(twr, tar)), 1e-13)
    expect_lt(max(prop_err(twc, twr)), 1e-10)
    expect_identical(prop_omp(EMC2:::calc_ll_manager(props, dadm_p, dp$model, return_trialwise = TRUE)), twc)
  }
})

test_that("proportional: C++ = R at extreme parameters (proportions at 0 and 1, A = 0, overflow)", {
  for (K in c(1, 3)) for (rule in c("BE", "TC")) {
    dp <- prop_design(K, rule, rho = TRUE)
    p <- prop_p(dp, K, -.5)
    set.seed(4600 + K)
    dat <- make_missing(make_data(p, dp, n_trials = 40), LC = .45, UC = 1.6)
    dadm <- EMC2:::design_model(dat, dp, verbose = FALSE, rt_resolution = 1e-6)
    ext <- matrix(p, nrow = 7, ncol = length(p), byrow = TRUE, dimnames = list(NULL, names(p)))
    ext[1, "A"] <- -Inf                 # proportion 0: A = 0 (the LBA limit)
    ext[2, "A"] <- 9                    # proportion 1: A = T_{K-1} (K = 1: A = b)
    ext[3, "A"] <- -5                   # tiny A
    ext[4, "B"] <- 800                  # overflow: out of bounds
    ext[5, "v"] <- NaN
    if (K > 1) {
      ext[6, "c1"] <- 9                 # proportion 1: out of bounds
      ext[7, "c2"] <- -9                # T_1 near 0 (BE) / near A (TC)
    }
    llc <- EMC2:::calc_ll_manager(ext, dadm, dp$model, return_trialwise = TRUE)
    llr <- prop_tw_r(ext, dadm, dp)
    expect_lt(max(prop_err(llc, llr)), 1e-10)
    expect_equal(as.vector(llc[, 4:5]), rep(log(1e-10), 2 * nrow(llc)))
    if (K > 1) expect_equal(as.vector(llc[, 6]), rep(log(1e-10), nrow(llc)))
    expect_identical(prop_omp(EMC2:::calc_ll_manager(ext, dadm, dp$model, return_trialwise = TRUE)), llc)
  }
})

test_that("proportional: make_data equals make_data at the translated additive parameters", {
  for (x in list(list(K = 3, rule = "BE", rho = -.5, guess = "accumulator", dk = FALSE),
                 list(K = 3, rule = "TC", rho = 0, guess = "trial", dk = FALSE),
                 list(K = 2, rule = "BE", rho = .8, guess = "none", dk = TRUE),
                 list(K = 4, rule = "BE", rho = 0, guess = "accumulator", dk = FALSE))) {
    dp <- prop_design(x$K, x$rule, dk = x$dk, rho = x$rho != 0, guess = x$guess)
    da <- prop_design(x$K, x$rule, dk = x$dk, rho = x$rho != 0, guess = x$guess, param = "additive")
    p <- prop_p(dp, x$K, x$rho, x$guess)
    set.seed(1); d0 <- make_data(p, dp, n_trials = 2)
    dadm <- suppressWarnings(EMC2:::design_model(d0, dp, compress = FALSE, verbose = FALSE))
    q <- prop_translate(p, dp, da, dadm, x$K, x$rule == "BE" && x$K > 1)
    set.seed(4700); dpd <- make_data(p, dp, n_trials = 300)
    set.seed(4700); dad <- make_data(q, da, n_trials = 300)
    expect_identical(dpd$R, dad$R)
    expect_identical(dpd$RR, dad$RR)
    expect_equal(dpd$rt, dad$rt, tolerance = 1e-12)
    # and with censoring
    set.seed(4701); cp <- make_data(p, dp, n_trials = 100, TC = list(LC = .45, UC = 1.6))
    set.seed(4701); ca <- make_data(q, da, n_trials = 100, TC = list(LC = .45, UC = 1.6))
    expect_identical(cp$missingness, ca$missingness)
    expect_identical(cp$RR, ca$RR)
  }
})

test_that("proportional: fits sample and predict; map = TRUE reports the natural A and b", {
  for (rule in c("BE", "TC")) {
    des <- suppressMessages(design(
      factors = list(subjects = 1, S = c("new", "old")), Rlevels = c("new", "old"),
      model = local({ r <- rule; function() MTLBA(3, rule = r, parameterization = "proportional") }),
      matchfun = function(d) d$S == d$lR, report_p_vector = FALSE,
      formula = list(v ~ lM, B ~ lR, A ~ 1, t0 ~ 1, rho ~ 1, c1 ~ lR, c2 ~ 1)))
    p <- c(v = 1, v_lMTRUE = 1.5, B = log(1.4), B_lRold = .1, A = qnorm(.6), t0 = log(.25),
           rho = qnorm(.2), c1 = qnorm(.7), c1_lRold = .2, c2 = qnorm(.3))
    set.seed(74)
    dat <- make_data(p, des, n_trials = 150)
    emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 2, verbose = FALSE))
    emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 5), cores_for_chains = 1,
                   cores_per_chain = 1, verbose = FALSE)
    expect_true(all(is.finite(emc[[1]]$samples$subj_ll)))
    pp <- predict(emc, n_post = 2, n_cores = 1)
    expect_equal(nrow(pp), 2 * nrow(dat))
    expect_true(all(pp$RR %in% 1:3))
    # mapped: A by accumulator (it scales with b and c1), as the natural start-point range
    al <- get_pars(emc, selection = "alpha", return_mcmc = FALSE, merge_chains = TRUE)
    gp <- get_pars(emc, map = TRUE, add_recalculated = TRUE, selection = "alpha",
                   return_mcmc = FALSE, merge_chains = TRUE)
    expect_true(all(c("A_lRnew", "A_lRold", "b_lRnew", "b_lRold", "B_lRnew") %in% dimnames(gp)[[1]]))
    expect_false("A" %in% dimnames(gp)[[1]])
    for (k in c(1, dim(al)[3])) {
      mp <- mapped_pars(des, al[, 1, k], digits = 12)
      nat <- unique(mp[, c("lR", "A", "b")])
      expect_equal(gp["A_lRold", 1, k], nat$A[nat$lR == "old"], tolerance = 1e-6)
      expect_equal(gp["A_lRnew", 1, k], nat$A[nat$lR == "new"], tolerance = 1e-6)
      expect_equal(gp["b_lRold", 1, k], nat$b[nat$lR == "old"], tolerance = 1e-6)
      expect_equal(gp["B_lRold", 1, k], gp["b_lRold", 1, k])
      expect_true(nat$A[1] < nat$b[1])
    }
    gp0 <- get_pars(emc, map = TRUE, selection = "alpha", return_mcmc = FALSE, merge_chains = TRUE)
    expect_true(all(c("A_lRnew", "A_lRold") %in% dimnames(gp0)[[1]]))
  }
})

test_that("proportional: DMC's sampled values reproduce the DMC code's density (golden dmc_p3_density)", {
  g <- prop_golden("dmc_p3_density.rds")
  # DMC's sampled values (transform.dmc of MTRp3_censor.R inverted): b, zd2 = qnorm(d2 / b),
  # zd1 = qnorm(d1 / d2), zA = qnorm(A / d2), zfb = qnorm((b - fb) / (b - A)), zf = qnorm(pf),
  # z2 = qnorm(p2), z3 = qnorm(p3), zcorr = qnorm((corr_v / tmp + 1) / 2) with
  # tmp = min(sd_v ratio); mapped by the M8p.5 table: B = log b, c1 = zd2, c2 = zd1, A = zA,
  # gB = -zfb, gp = zf, pg by stick-breaking, rho = zcorr (tmp = 1; the m5 set has sd_v by
  # accumulator, tmp != 1, never so in the paper's models: rho from corr_v there)
  emc2_p <- function(s, des) {
    zd2 <- qnorm(s$d2 / s$b); zd1 <- qnorm(s$d1 / s$d2); zA <- qnorm(s$A / s$d2)
    tmp <- min(s$sd_v[1] / s$sd_v[2], s$sd_v[2] / s$sd_v[1])
    zcorr <- qnorm((s$corr_v / tmp + 1) / 2)
    rho <- if (tmp == 1 || s$corr_v == 0) zcorr else qnorm((s$corr_v + 1) / 2)
    cells <- list(B = log(s$b), c1 = zd2, c2 = zd1, A = zA, v = s$mean_v, sv = log(s$sd_v))
    if (any(s$pf > 0)) {
      zfb <- qnorm((s$b - s$fb) / (s$b - s$A)); zf <- qnorm(s$pf); z2 <- qnorm(s$p2); z3 <- qnorm(s$p3)
      pg1 <- qnorm((1 - pnorm(z2)) * (1 - pnorm(z3)))
      cells <- c(cells, list(gB = -zfb, gp = zf, pg1 = pg1, pg2 = qnorm(pnorm(z2) / (1 - pnorm(pg1)))))
    }
    p <- sampled_pars(des, doMap = FALSE)
    p[] <- NA
    for (nm in names(cells)) p[c(nm, paste0(nm, "_lROLD"))] <- c(cells[[nm]][1], diff(cells[[nm]]))
    p[c("t0", "rho")] <- c(log(s$t0), rho)
    expect_false(anyNA(p))
    p
  }
  mk <- function(guess) {
    form <- c(list(v ~ lR, sv ~ lR, B ~ lR, A ~ lR, t0 ~ 1, rho ~ 1, c1 ~ lR, c2 ~ lR),
              if (guess) list(gB ~ lR, gp ~ lR, pg1 ~ lR, pg2 ~ lR))
    suppressMessages(design(factors = list(subjects = 1, S = c("NEW", "OLD")), Rlevels = c("NEW", "OLD"),
                            model = function() MTLBA(3, parameterization = "proportional"),
                            matchfun = function(d) d$S == d$lR, formula = form, report_p_vector = FALSE))
  }
  cells <- g$cells
  for (sn in names(g$sets)) {
    s <- g$sets[[sn]]
    des <- mk(any(s$pf > 0))
    p <- emc2_p(s, des)
    w <- which(cells$set == sn)
    dat <- data.frame(subjects = factor(1), S = factor("NEW", levels = c("NEW", "OLD")),
                      R = factor(ifelse(cells$order[w] == "n1", "NEW", "OLD"), levels = c("NEW", "OLD")),
                      RR = cells$r2[w], rt = cells$rt[w])
    dadm <- EMC2:::design_model(dat, des, compress = FALSE, verbose = FALSE, rt_resolution = 1e-15)
    # the R reference through the factory's Ttransform (no likelihood floor) ...
    pm <- EMC2:::get_pars_matrix_oo(p, dadm, des$model())
    ours <- EMC2:::dMTLBA(dat$rt, as.integer(dat$R), dat$RR, pm, 3, "BE", guess = "accumulator")
    ref <- cells$dens[w]
    body <- ref > 1e-8
    expect_lt(max(abs(ours[body] / ref[body] - 1)), 1e-8)
    expect_lt(max(abs(ours - ref)), 1e-12)
    # ... and the compiled likelihood (floored at 1e-10: the body only)
    cpp <- exp(as.vector(EMC2:::calc_ll_manager(t(p), dadm, des$model, return_trialwise = TRUE)))
    expect_lt(max(abs(cpp[body] / ref[body] - 1)), 1e-8)
    # censored masses P(t0 < rt < CF) and P(rt > CS), through the compiled censoring terms
    cs <- g$cens[g$cens$set == sn, ]
    cd <- data.frame(subjects = factor(1), S = factor("NEW", levels = c("NEW", "OLD")),
                     R = factor(NA, levels = c("NEW", "OLD")), RR = NA_real_, rt = NA_real_,
                     missingness = ifelse(cs$kind == "fast", 1L, 2L),
                     LC = ifelse(cs$kind == "fast", cs$bound, 0), UC = ifelse(cs$kind == "slow", cs$bound, Inf))
    cdadm <- suppressWarnings(EMC2:::design_model(cd, des, compress = FALSE, verbose = FALSE))
    cm <- exp(as.vector(EMC2:::calc_ll_manager(t(p), cdadm, des$model, return_trialwise = TRUE)))
    expect_equal(cm, cs$mass, tolerance = 1e-8)
  }
})
