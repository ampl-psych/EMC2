# MTLBA(guess_ratings = "binomial"): a guessed rating is 1 + Binomial(K - 1, q), with one
# p_type pq (q = pnorm(pq)) in place of pg1 .. pg{K-1}. It is implemented as the exact
# sequential pg_k = P(k) / P(>= k) of the free model, so the checks are: the hand
# probabilities, binomial = free at the mapped pg (R and C++, every option), C++ = R,
# q near 0 and 1, simulation, fitting and mapping.

bin_omp <- function(expr, threads = 2) {
  old <- options(emc.ll_backend = "multithreaded", emc.n_threads = threads)
  on.exit(options(old))
  expr
}
# R reference trial log-likelihoods, particles in columns as calc_ll_manager(return_trialwise = TRUE)
bin_tw_r <- function(props, dadm, des) {
  m <- des$model()
  apply(props, 1, function(pv)
    EMC2:::trial_ll_mtlba(EMC2:::get_pars_matrix_oo(pv, dadm, m), dadm, m$n_ratings, m$rule, m$dont_know,
                          guess = m$guess, gb_constraint = m$gb_constraint, drifts = m$drifts))
}
# error relative to max(1, |reference|), elementwise
bin_err <- function(a, b) abs(a - b) / pmax(abs(b), 1)

# Guessing design with the guess-rating parameters by accumulator: pq ~ lR, or for the free
# comparison model pg1 .. ~ 0 + lR with an identity transform, so that it takes the mapped pg
# exactly (through the probit scale a pg near 1 would come back an ulp off, which log1p(-pg)
# turns into differences of up to ~1e-10 in a trial log-likelihood)
bin_design <- function(K, rule = "BE", dk = FALSE, rho = FALSE, guess = "accumulator",
                       param = "additive", gr = "binomial", drifts = "any") {
  crit <- paste0("c", seq_len(K - 1))
  pg <- paste0("pg", seq_len(K - 1))
  form <- c(list(v ~ lM, sv ~ lM, B ~ lR, A ~ lR, t0 ~ 1), if (rho) list(rho ~ 1),
            if (guess == "accumulator") list(gp ~ lR, gB ~ lR) else list(gp ~ 1, gB ~ lR),
            lapply(c(crit, if (gr == "binomial") "pq"), function(cn) stats::as.formula(paste(cn, "~ lR"))),
            if (gr == "free") lapply(pg, function(cn) stats::as.formula(paste(cn, "~ 0 + lR"))))
  mod <- local({
    k <- K; r <- rule; d <- dk; g <- guess; pa <- param; gg <- gr; dr <- drifts; pgn <- pg
    function() {
      m <- MTLBA(n_ratings = k, dont_know = d, rule = r, guess = g, parameterization = pa,
                 guess_ratings = gg, drifts = dr)
      if (gg == "free") m$transform$func[pgn] <- "identity"
      m
    }
  })
  suppressMessages(design(
    factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
    model = mod, matchfun = function(d) d$S == d$lR, formula = form, report_p_vector = FALSE))
}

# A binomial parameter vector for bin_design, q = (q_left, q_right)
bin_p <- function(des, K, q = c(.3, .7), rho = 0, guess = "accumulator", param = "additive") {
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- 0
  p[c("v", "v_lMTRUE", "sv", "sv_lMTRUE", "B", "B_lRright", "A", "A_lRright", "t0")] <-
    c(1, 1.5, 0, log(.7), log(1.5), .2, if (param == "additive") log(.5) else qnorm(.6), .2, log(.25))
  if (rho != 0) p["rho"] <- qnorm((1 + rho) / 2)
  for (k in seq_len(K - 1))
    p[paste0("c", k, c("", "_lRright"))] <- c(if (param == "additive") log(.6) else qnorm(.6), .15)
  p[c("gp", "gB", "gB_lRright")] <- c(qnorm(.25), qnorm(.5), .3)
  if (guess == "accumulator") p["gp_lRright"] <- -.4
  p[c("pq", "pq_lRright")] <- c(qnorm(q[1]), qnorm(q[2]) - qnorm(q[1]))
  p
}

# The free parameter vector at which the free model (bin_design(gr = "free")) is the binomial
# one at p: the pg the binomial Ttransform gives each accumulator (first trial's two rows).
bin_translate <- function(p, des_b, des_f, dadm_b, K) {
  pm <- EMC2:::get_pars_matrix_oo(p, dadm_b, des_b$model())
  lv <- as.integer(dadm_b$lR)
  rows <- c(which(lv == 1)[1], which(lv == 2)[1])
  f <- sampled_pars(des_f, doMap = FALSE)
  f[] <- p[match(names(f), names(p))]
  for (k in seq_len(K - 1)) f[paste0("pg", k, c("_lRleft", "_lRright"))] <- pm[rows, paste0("pg", k)]
  f
}

test_that("binomial: factory (p_types, transforms, c_options); free default unchanged", {
  for (K in 1:4) {
    pg <- if (K > 1) paste0("pg", seq_len(K - 1)) else character(0)
    f <- MTLBA(K)
    f2 <- MTLBA(K, guess_ratings = "free")
    m <- MTLBA(K, guess_ratings = "binomial")
    expect_identical(f$guess_ratings, "free")
    expect_identical(m$guess_ratings, "binomial")
    for (x in c("p_types", "transform", "bound", "c_options")) expect_identical(f[[x]], f2[[x]])
    expect_identical(f$c_options[["binomial_guess"]], 0)
    expect_identical(names(f$p_types)[seq_along(pg) + length(f$p_types) - length(pg)], pg)
    expect_identical(unname(f$p_types[pg]), qnorm(as.numeric(seq_along(pg) == 1)))
    if (K == 1) {
      expect_identical(m$c_options[["binomial_guess"]], 0)   # one rating: nothing to guess
      expect_false("pq" %in% names(m$p_types))
      next
    }
    expect_identical(m$c_options[["binomial_guess"]], 1)
    expect_identical(names(m$p_types), c(setdiff(names(f$p_types), pg), "pq"))
    expect_identical(m$p_types[["pq"]], -Inf)   # q = 0: guessing off guesses the lowest rating
    expect_identical(m$transform$func[["pq"]], "pnorm")
    expect_false(any(pg %in% names(m$transform$func)))
    expect_identical(MTLBA(K, guess_ratings = "binomial", parameterization = "proportional")$c_options[
      c("proportional", "binomial_guess")], c(proportional = 1, binomial_guess = 1))
  }
  expect_error(MTLBA(3, guess_ratings = "logit"))
})

test_that("binomial: guess-rating probabilities are Binomial(K - 1, q), K = 2 .. 4, q from 0 to 1", {
  qs <- c(0, 1e-300, 1e-10, 1e-4, .1, .3, .5, .7, .9, 1 - 1e-4, 1 - 1e-10, 1)
  # the sequential pg in closed forms well conditioned for every q
  closed <- list(`2` = function(q) cbind(1 - q),
                 `3` = function(q) cbind((1 - q)^2, 2 * (1 - q) / (2 - q)),
                 `4` = function(q) cbind((1 - q)^3, 3 * (1 - q)^2 / (3 - 3 * q + q^2), 3 * (1 - q) / (3 - 2 * q)))
  for (K in 2:4) {
    pg <- EMC2:::mtlba_binomial_pg(qs, K)
    expect_identical(colnames(pg), paste0("pg", seq_len(K - 1)))
    ref <- closed[[as.character(K)]](qs)
    expect_lt(max(abs(pg / ref - 1)[ref > 0]), 1e-14)
    expect_true(all(pg[ref == 0] == 0))
    # the implied distribution through mtlba_guess_pars (the likelihood's and simulator's reader)
    pars <- cbind(A = .5, b = 1.5, gp = .3, gB = .5, pg)
    lpg <- EMC2:::mtlba_guess_pars(pars, K)$lpg
    dens <- t(sapply(qs, function(q) dbinom(0:(K - 1), K - 1, q)))
    expect_lt(max(abs(exp(lpg) - dens)), 4e-16)                     # absolute, every q
    mid <- qs >= .1 & qs <= .9
    expect_lt(max(abs(exp(lpg[mid, ]) / dens[mid, ] - 1)), 1e-13)   # relative, q away from 0 and 1
    expect_equal(rowSums(exp(lpg)), rep(1, length(qs)), tolerance = 1e-15)
    # K = 3 by hand: ((1 - q)^2, 2 q (1 - q), q^2)
    if (K == 3) expect_equal(exp(lpg[qs == .3, ]), c(.49, .42, .09), tolerance = 1e-14)
    # q near 1: every probability (here also the tiny ones) to relative precision
    hi <- qs %in% c(1 - 1e-4, 1 - 1e-10)
    expect_lt(max(abs(exp(lpg[hi, ]) / dens[hi, ] - 1)), 1e-12)
  }
  expect_identical(EMC2:::mtlba_binomial_pg(NaN, 3)[1, ], c(pg1 = NaN, pg2 = NaN))
})

test_that("binomial = free at the mapped pg, in R and C++; C++ = R; OpenMP identical", {
  # (the full grid, K 2-4 x BE/TC x additive/proportional x guess accumulator/trial x
  # censoring x don't know with q at 1e-10, .5 and 1 - 1e-10 and random vectors, is run by
  # the stage's check script; this is a thinned version of it)
  cfg <- expand.grid(K = 2:4, rule = c("BE", "TC"), param = c("additive", "proportional"),
                     guess = c("accumulator", "trial"), cens = c(FALSE, TRUE), dk = c(FALSE, TRUE),
                     stringsAsFactors = FALSE)
  set.seed(5401)
  cfg <- cfg[sort(sample(nrow(cfg), 12)), ]
  for (i in seq_len(nrow(cfg))) {
    x <- cfg[i, ]
    rho <- c(0, -.5, .8)[i %% 3 + 1]
    db <- bin_design(x$K, x$rule, x$dk, rho != 0, x$guess, x$param)
    df <- bin_design(x$K, x$rule, x$dk, rho != 0, x$guess, x$param, gr = "free")
    p <- bin_p(db, x$K, rho = rho, guess = x$guess, param = x$param)
    set.seed(5500 + i)
    dat <- make_data(p, db, n_trials = 40)
    if (x$cens) {
      dat <- make_missing(dat, LC = .45, UC = 1.6, LT = .3, UT = 4)
      dat$missingness[1:2] <- 3L; dat$R[1:2] <- NA; dat$RR[1:2] <- NA; dat$rt[1:2] <- NA
    }
    props <- rbind(p, mvtnorm::rmvnorm(3, p, diag(length(p)) * .1))
    colnames(props) <- names(p)
    # q at the extremes and at .5, on one accumulator or both
    props[2, c("pq", "pq_lRright")] <- c(qnorm(1e-10), qnorm(1 - 1e-10) - qnorm(1e-10))
    props[3, c("pq", "pq_lRright")] <- c(qnorm(.5), qnorm(1e-10) - qnorm(.5))
    props[4, c("pq", "pq_lRright")] <- c(qnorm(1 - 1e-10), 0)
    dadm_b <- EMC2:::design_model(dat, db, compress = !x$cens, verbose = FALSE, rt_resolution = 1e-6)
    dadm_f <- EMC2:::design_model(dat, df, compress = !x$cens, verbose = FALSE, rt_resolution = 1e-6)
    pf <- t(apply(props, 1, bin_translate, des_b = db, des_f = df, dadm_b = dadm_b, K = x$K))
    twc <- EMC2:::calc_ll_manager(props, dadm_b, db$model, return_trialwise = TRUE)
    twr <- bin_tw_r(props, dadm_b, db)
    tfc <- EMC2:::calc_ll_manager(pf, dadm_f, df$model, return_trialwise = TRUE)
    tfr <- bin_tw_r(pf, dadm_f, df)
    expect_true(all(is.finite(twc)))
    expect_lt(max(bin_err(twc, tfc)), 1e-13)
    expect_lt(max(bin_err(twr, tfr)), 1e-13)
    expect_lt(max(bin_err(twc, twr)), 1e-10)
    expect_identical(bin_omp(EMC2:::calc_ll_manager(props, dadm_b, db$model, return_trialwise = TRUE)), twc)
  }
})

test_that("binomial: C++ = R at q = 0 and 1 exactly and at non-finite pq", {
  for (K in c(2, 3)) for (rule in c("BE", "TC")) {
    db <- bin_design(K, rule, rho = TRUE)
    p <- bin_p(db, K, rho = -.5)
    set.seed(5600 + K)
    dat <- make_missing(make_data(p, db, n_trials = 40), LC = .45, UC = 1.6)
    dadm <- EMC2:::design_model(dat, db, verbose = FALSE, rt_resolution = 1e-6)
    ext <- matrix(p, nrow = 4, ncol = length(p), byrow = TRUE, dimnames = list(NULL, names(p)))
    ext[1, c("pq", "pq_lRright")] <- c(-Inf, 0)    # q = 0: every guess is rating 1
    ext[2, c("pq", "pq_lRright")] <- c(Inf, 0)     # q = 1: every guess is rating K
    ext[3, c("pq", "pq_lRright")] <- c(-40, 80)    # q = 0 on one accumulator, 1 on the other
    ext[4, "pq"] <- NaN
    llc <- EMC2:::calc_ll_manager(ext, dadm, db$model, return_trialwise = TRUE)
    llr <- bin_tw_r(ext, dadm, db)
    expect_lt(max(bin_err(llc[, 1:3], llr[, 1:3])), 1e-10)
    expect_true(all(llc[, 1:3] >= log(1e-10)) && all(is.finite(llc)))
    expect_identical(bin_omp(EMC2:::calc_ll_manager(ext, dadm, db$model, return_trialwise = TRUE)), llc)
  }
})

test_that("binomial: make_data equals make_data of the free model at the mapped pg", {
  for (x in list(list(K = 3, rule = "BE", rho = -.5, guess = "accumulator", dk = FALSE, param = "additive"),
                 list(K = 3, rule = "TC", rho = 0, guess = "trial", dk = FALSE, param = "proportional"),
                 list(K = 2, rule = "BE", rho = .8, guess = "accumulator", dk = TRUE, param = "additive"),
                 list(K = 4, rule = "BE", rho = 0, guess = "accumulator", dk = FALSE, param = "proportional"))) {
    db <- bin_design(x$K, x$rule, x$dk, x$rho != 0, x$guess, x$param)
    df <- bin_design(x$K, x$rule, x$dk, x$rho != 0, x$guess, x$param, gr = "free")
    p <- bin_p(db, x$K, q = c(.2, .65), rho = x$rho, guess = x$guess, param = x$param)
    p["gp"] <- qnorm(.4)
    set.seed(1); d0 <- make_data(p, db, n_trials = 2)
    dadm <- suppressWarnings(EMC2:::design_model(d0, db, compress = FALSE, verbose = FALSE))
    f <- bin_translate(p, db, df, dadm, x$K)
    set.seed(5700); dbd <- make_data(p, db, n_trials = 300)
    set.seed(5700); dfd <- make_data(f, df, n_trials = 300)
    expect_identical(dbd$R, dfd$R)
    expect_identical(dbd$RR, dfd$RR)
    expect_identical(dbd$rt, dfd$rt)
    set.seed(5701); cb <- make_data(p, db, n_trials = 100, TC = list(LC = .45, UC = 1.6))
    set.seed(5701); cf <- make_data(f, df, n_trials = 100, TC = list(LC = .45, UC = 1.6))
    expect_identical(cb$missingness, cf$missingness)
    expect_identical(cb$RR, cf$RR)
  }
  # q -> 0 guesses rating 1, q -> 1 rating K (the guessing-only limit, gp = 1 on both accumulators)
  db <- bin_design(3)
  for (q in c(1e-10, 1 - 1e-10)) {
    p <- bin_p(db, 3, q = c(q, q))
    p[c("gp", "gp_lRright")] <- c(Inf, 0)
    set.seed(5702)
    d <- make_data(p, db, n_trials = 50)
    expect_true(all(d$RR == if (q < .5) 1 else 3))
  }
})

test_that("binomial: a fit samples and predicts; mapped_pars and get_pars show pq and the implied pg", {
  des <- suppressMessages(design(
    factors = list(subjects = 1, S = c("new", "old")), Rlevels = c("new", "old"),
    model = function() MTLBA(3, guess_ratings = "binomial"),
    matchfun = function(d) d$S == d$lR, report_p_vector = FALSE,
    formula = list(v ~ lM, B ~ lR, A ~ 1, t0 ~ 1, c1 ~ lR, c2 ~ 1, gp ~ 1, gB ~ 1, pq ~ 1)))
  expect_true("pq" %in% names(sampled_pars(des)))
  expect_false(any(c("pg1", "pg2") %in% names(sampled_pars(des))))
  p <- c(v = 1, v_lMTRUE = 1.5, B = log(1.4), B_lRold = .1, A = log(.5), t0 = log(.25),
         c1 = log(.6), c1_lRold = .2, c2 = log(.8), gp = qnorm(.2), gB = 0, pq = qnorm(.3))
  mp <- mapped_pars(des, p, digits = 12)
  expect_equal(unique(mp$pq), .3)
  expect_equal(unique(mp$pg1), .49)
  expect_equal(unique(mp$pg2), .42 / .51)
  set.seed(75)
  dat <- make_data(p, des, n_trials = 150)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 2, verbose = FALSE))
  emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 5), cores_for_chains = 1,
                 cores_per_chain = 1, verbose = FALSE)
  expect_true(all(is.finite(emc[[1]]$samples$subj_ll)))
  pp <- predict(emc, n_post = 2, n_cores = 1)
  expect_equal(nrow(pp), 2 * nrow(dat))
  expect_true(all(pp$RR %in% 1:3))
  gp0 <- get_pars(emc, map = TRUE, selection = "alpha", return_mcmc = FALSE, merge_chains = TRUE)
  expect_true("pq" %in% dimnames(gp0)[[1]])
  expect_false("pg1" %in% dimnames(gp0)[[1]])
  gp <- get_pars(emc, map = TRUE, add_recalculated = TRUE, selection = "alpha",
                 return_mcmc = FALSE, merge_chains = TRUE)
  expect_true(all(c("pq", "pg1", "pg2") %in% dimnames(gp)[[1]]))
  ref <- EMC2:::mtlba_binomial_pg(gp["pq", 1, ], 3)
  expect_equal(gp["pg1", 1, ], ref[, 1], tolerance = 1e-6)
  expect_equal(gp["pg2", 1, ], ref[, 2], tolerance = 1e-6)
})
