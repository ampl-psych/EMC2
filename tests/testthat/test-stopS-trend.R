## SSD-dependent stop-signal parameters via trends (dEXG3 / make_ssd_trend()):
## C++ likelihood == R reference (two ways), memo-key correctness, compression
## and multithread invariance, staircase simulation on the trial-by-trial path,
## the tf(SSD) variant, and SSRDEX transfer.
RNGkind("L'Ecuyer-CMRG")
set.seed(321)

ss_matchfun <- function(d) d$S == d$lR

ss_design_tr <- function(model, trend = make_ssd_trend(), ...) {
  design(model = model, factors = list(subjects = 1, S = c("left", "right")),
         Rlevels = c("left", "right"), matchfun = ss_matchfun, report_p_vector = FALSE,
         trend = trend, ...)
}

exg_formula <- list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
rdex_formula <- list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
exg_vals <- c(mu = log(.5), mu_lMTRUE = .2, sigma = log(.05), tau = log(.1), muS = log(.12),
              sigmaS = log(.03), tauS = log(.05), gf = qnorm(.05), tf = qnorm(.1))
rdex_vals <- c(v = log(2), v_lMTRUE = .5, B = log(1), A = log(.3), t0 = log(.15), muS = log(.12),
               sigmaS = log(.03), tauS = log(.05), gf = qnorm(.05), tf = qnorm(.1))

# named parameter vector in sampled_pars() order
p_named <- function(des, values) {
  p <- sampled_pars(des, doMap = FALSE)
  stopifnot(setequal(names(p), names(values)))
  p[] <- values[names(p)]
  p
}

# C++ (c_name) and R (c_name nulled) log-likelihoods of one p_vector on a dadm
ss_ll_cr <- function(dat, des, p_vector, compress = FALSE) {
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = compress, verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  pm <- matrix(p_vector, nrow = 1, dimnames = list(NULL, names(p_vector)))
  ll_c <- as.numeric(calc_ll_manager(pm, dadm, des$model))
  model_r <- des$model()
  model_r$c_name <- NULL
  ll_r <- as.numeric(calc_ll_manager(pm, dadm, function() model_r))
  c(cpp = ll_c, r = ll_r)
}

# Independent R reference: untrended design, parameter matrix with muS
# substituted trial-wise by the dEXG3 formula, fed to the R likelihood.
ss_ll_manual_muS <- function(dat, des0, p_tr) {
  emc <- make_emc(dat, des0, type = "single", n_chains = 1, compress = FALSE, verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  model_r <- des0$model()
  model_r$c_name <- NULL
  p0 <- p_tr[names(sampled_pars(des0, doMap = FALSE))]
  pars <- get_pars_matrix_oo(p0, dadm, model_r)
  k <- exp(p_tr[["muS.k_sat"]]); w <- exp(p_tr[["muS.w"]])
  ssd <- dadm$SSD
  pars[, "muS"] <- pars[, "muS"] + ifelse(is.finite(ssd), w * pmin(1, k * ssd), 0)
  as.numeric(model_r$log_likelihood(pars, dadm, model_r))
}

# with k = 5 the ramp clips inside the SSD range: .1 -> .5, .2 -> 1, .3/.4 -> 1
fixed_ssd <- make_ssd(staircase = FALSE, values = c(.1, .2, .3, .4))
tr_vals <- c(muS.k_sat = log(5), muS.w = log(.2))

test_that("make_ssd_trend builds the dEXG3 specification", {
  tr <- make_ssd_trend()
  expect_s3_class(tr, "emc2_trend")
  expect_equal(get_trend_pnames(tr), c("muS.k_sat", "muS.w"))
  expect_equal(tr$bases[[1]]$phase, "posttransform")
  expect_equal(unname(tr$bases[[1]]$transforms["muS.w"]), "exp")
  expect_equal(unname(tr$kernels[[1]]$transforms["muS.k_sat"]), "exp")
  tr2 <- make_ssd_trend(c("muS", "tf"))
  expect_setequal(get_trend_pnames(tr2), c("muS.k_sat", "muS.w", "tf.k_sat", "tf.w"))
  expect_equal(tr2$bases[[2]]$phase, "pretransform")
  expect_equal(unname(tr2$bases[[2]]$transforms["tf.w"]), "identity")
  tr3 <- make_ssd_trend(c("muS", "tf"), shared_k = TRUE)
  expect_setequal(get_trend_pnames(tr3), c("muS.k_sat", "muS.w", "tf.w"))
  expect_error(make_ssd_trend("mu"), "target must be")
  des <- ss_design_tr(SSEXG, formula = exg_formula)
  expect_true(all(c("muS.k_sat", "muS.w") %in% names(sampled_pars(des))))
  expect_equal(unname(des$model()$transform$func[c("muS.k_sat", "muS.w")]), c("exp", "exp"))
})

test_that("SSEXG + muS(SSD): C++ == R reference, with and without a deadline", {
  des <- ss_design_tr(SSEXG, formula = exg_formula)
  p <- p_named(des, c(exg_vals, tr_vals))
  set.seed(11)
  dat0 <- make_data(p, des, n_trials = 60, functions = list(SSD = fixed_ssd))
  expect_true(any(is.na(dat0$rt) & is.finite(dat0$SSD)))
  ll <- ss_ll_cr(dat0, des, p)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-8)
  set.seed(11)
  datU <- make_data(p, des, n_trials = 60, functions = list(SSD = fixed_ssd), TC = list(UC = 0.8))
  expect_true(any(datU$missingness %in% 2L & is.finite(datU$SSD)))
  llU <- ss_ll_cr(datU, des, p)
  expect_equal(llU[["cpp"]], llU[["r"]], tolerance = 1e-8)
  expect_false(isTRUE(all.equal(ll[["cpp"]], llU[["cpp"]])))
  # perturbed proposals (both clipped and unclipped regimes)
  set.seed(12)
  pm <- rbind(p, matrix(rnorm(15 * length(p), rep(p, each = 15), .2), 15))
  colnames(pm) <- names(p)
  pm[, "muS.k_sat"] <- log(c(5, seq(.5, 12, length.out = 15)))
  emc <- make_emc(datU, des, type = "single", n_chains = 1, compress = FALSE, verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  model_r <- des$model(); model_r$c_name <- NULL
  ll_c <- as.numeric(calc_ll_manager(pm, dadm, des$model))
  ll_r <- as.numeric(calc_ll_manager(pm, dadm, function() model_r))
  expect_equal(ll_c, ll_r, tolerance = 1e-7)
})

test_that("SSEXG + muS(SSD): C++ == hand-substituted R likelihood", {
  des  <- ss_design_tr(SSEXG, formula = exg_formula)
  des0 <- design(model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
                 Rlevels = c("left", "right"), matchfun = ss_matchfun, report_p_vector = FALSE,
                 formula = exg_formula)
  p <- p_named(des, c(exg_vals, tr_vals))
  set.seed(13)
  dat0 <- make_data(p, des, n_trials = 60, functions = list(SSD = fixed_ssd))
  expect_equal(ss_ll_cr(dat0, des, p)[["cpp"]], ss_ll_manual_muS(dat0, des0, p), tolerance = 1e-8)
  set.seed(13)
  datU <- make_data(p, des, n_trials = 60, functions = list(SSD = fixed_ssd), TC = list(UC = 0.8))
  expect_equal(ss_ll_cr(datU, des, p)[["cpp"]], ss_ll_manual_muS(datU, des0, p), tolerance = 1e-8)
  # a vanishing weight reproduces the untrended likelihood exactly
  p_off <- p; p_off["muS.w"] <- -30
  p0 <- p[names(sampled_pars(des0, doMap = FALSE))]
  emc0 <- make_emc(datU, des0, type = "single", n_chains = 1, compress = FALSE, verbose = FALSE)
  ll0 <- as.numeric(calc_ll_manager(matrix(p0, 1, dimnames = list(NULL, names(p0))),
                                    emc0[[1]]$data[[1]], des0$model))
  expect_equal(ss_ll_cr(datU, des, p_off)[["cpp"]], ll0, tolerance = 1e-10)
})

test_that("muS(SSD): trialwise muS follows the formula and is untouched on go trials", {
  des <- ss_design_tr(SSEXG, formula = exg_formula)
  p <- p_named(des, c(exg_vals, tr_vals))
  set.seed(14)
  dat <- make_data(p, des, n_trials = 60, functions = list(SSD = fixed_ssd),
                   return_trialwise_parameters = TRUE)
  tw <- attr(dat, "trialwise_parameters")
  tw1 <- tw[!duplicated(tw[, "trials"]), , drop = FALSE]     # first accumulator row per trial
  expect_equal(nrow(tw1), nrow(dat))
  k <- exp(p[["muS.k_sat"]]); w <- exp(p[["muS.w"]])
  expected <- exp(p[["muS"]]) + ifelse(is.finite(dat$SSD), w * pmin(1, k * dat$SSD), 0)
  expect_equal(unname(tw1[, "muS"]), expected, tolerance = 1e-12)
  expect_equal(unname(tw1[!is.finite(dat$SSD), "muS"]), rep(exp(p[["muS"]]), sum(!is.finite(dat$SSD))))
  expect_true(any(abs(tw1[, "muS"] - exp(p[["muS"]]) - w) < 1e-12))   # clipped trials exist
})

test_that("memo key: trials with equal SSD and go parameters but different muS", {
  # trend muS on an extra random covariate, so muS differs across trials that
  # share SSD, deadline and go parameters
  tr_x <- make_trend(make_base("muS", "lin", make_kernel("x", "sat_lin"),
                               phase = "posttransform", transforms = list(w = "exp")))
  des_x <- ss_design_tr(SSEXG, trend = tr_x, formula = exg_formula, covariates = "x")
  des <- ss_design_tr(SSEXG, formula = exg_formula)
  p <- p_named(des, c(exg_vals, tr_vals))
  set.seed(15)
  dat <- make_data(p, des, n_trials = 80, functions = list(SSD = fixed_ssd), TC = list(UC = 0.8))
  dat$x <- runif(nrow(dat), 0, 2)
  p_x <- p_named(des_x, c(exg_vals, muS.k_sat = log(1), muS.w = log(.15)))
  ll <- ss_ll_cr(dat, des_x, p_x)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-8)
  # and the memo must not confuse them with the SSD-only trend either
  emc <- make_emc(dat, des_x, type = "single", n_chains = 1, compress = FALSE, verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  set.seed(16)
  pm <- rbind(p_x, matrix(rnorm(10 * length(p_x), rep(p_x, each = 10), .2), 10))
  colnames(pm) <- names(p_x)
  ll1 <- as.numeric(calc_ll_manager(pm, dadm, des_x$model))
  old <- options(emc.ll_backend = "multithreaded", emc.n_threads = 2)
  on.exit(options(old))
  ll2 <- suppressWarnings(as.numeric(calc_ll_manager(pm, dadm, des_x$model)))
  expect_equal(ll1, ll2, tolerance = 0)
})

test_that("muS(SSD): compression invariance and multithread equality", {
  des <- ss_design_tr(SSEXG, formula = exg_formula)
  p <- p_named(des, c(exg_vals, tr_vals))
  set.seed(17)
  dat <- make_data(p, des, n_trials = 80, functions = list(SSD = fixed_ssd), TC = list(UC = 0.8))
  llF <- ss_ll_cr(dat, des, p, compress = FALSE)
  llT <- ss_ll_cr(dat, des, p, compress = TRUE)
  expect_equal(llF[["cpp"]], llT[["cpp"]], tolerance = 1e-10)
  expect_equal(llT[["cpp"]], llT[["r"]], tolerance = 1e-8)
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE, verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  set.seed(18)
  pm <- rbind(p, matrix(rnorm(10 * length(p), rep(p, each = 10), .2), 10))
  colnames(pm) <- names(p)
  ll1 <- as.numeric(calc_ll_manager(pm, dadm, des$model))
  old <- options(emc.ll_backend = "multithreaded", emc.n_threads = 2)
  on.exit(options(old))
  ll2 <- suppressWarnings(as.numeric(calc_ll_manager(pm, dadm, des$model)))
  expect_equal(ll1, ll2, tolerance = 0)
})

test_that("staircase + trend on SSD simulates trial by trial with a valid ladder", {
  des <- ss_design_tr(SSEXG, formula = exg_formula)
  p <- p_named(des, c(exg_vals, tr_vals))
  stair <- make_ssd(staircase = TRUE, SSD0 = .25, stairstep = .05, p_stop = .3)
  set.seed(19)
  expect_message(
    dat <- make_data(p, des, n_trials = 200, functions = list(SSD = stair),
                     return_trialwise_parameters = TRUE),
    "trial by trial")
  expect_s3_class(attr(dat, "staircase"), "emc_staircase")
  st <- dat[is.finite(dat$SSD), ]
  expect_gt(nrow(st), 30)
  expect_equal(st$SSD[1], .25)
  d <- diff(st$SSD)
  expect_true(all(abs(d) < 1e-9 | abs(abs(d) - .05) < 1e-9))
  # up after a stop success, down after a signal-respond (0 = clamped at stairmin)
  up <- is.na(st$R)[-nrow(st)]
  expect_true(all((d[d != 0] > 0) == up[d != 0]))
  expect_true(all(st$SSD >= 0))
  # per-trial muS uses the realised SSD
  tw <- attr(dat, "trialwise_parameters")
  tw1 <- tw[!duplicated(tw[, c("subject", "trial")]), , drop = FALSE]
  k <- exp(p[["muS.k_sat"]]); w <- exp(p[["muS.w"]])
  expected <- exp(p[["muS"]]) + ifelse(is.finite(dat$SSD), w * pmin(1, k * dat$SSD), 0)
  expect_equal(unname(tw1[, "muS"]), expected, tolerance = 1e-12)
  # explicit conditional_on_data = TRUE is refused
  expect_error(make_data(p, des, n_trials = 50, functions = list(SSD = stair), conditional_on_data = TRUE),
               "trial by trial")
})

test_that("staircase + trend: deadline steps the ladder up after late responses; grouped ladders", {
  des <- ss_design_tr(SSEXG, formula = exg_formula, TC = list(UC = 0.55))
  p <- p_named(des, c(exg_vals, tr_vals))
  stair <- make_ssd(staircase = TRUE, SSD0 = .2, stairstep = .05, p_stop = .3)
  set.seed(20)
  dat <- suppressMessages(make_data(p, des, n_trials = 200, functions = list(SSD = stair)))
  expect_true(all(dat$UC == 0.55))
  st <- dat[is.finite(dat$SSD), ]
  expect_true(any(st$missingness %in% 2L))                   # late/withheld coded 2
  expect_true(all(is.na(st$rt[st$missingness %in% 2L])))
  d <- diff(st$SSD)
  up <- is.na(st$R)[-nrow(st)]                                # late responses are NA after make_missing
  expect_true(all((d[d != 0] > 0) == up[d != 0]))
  # separate ladders per S, each starting at SSD0, one subject
  stairS <- make_ssd(staircase = TRUE, SSD0 = .3, stairstep = .05, p_stop = .3, factors = "S")
  des2 <- design(model = SSEXG, factors = list(subjects = 1:2, S = c("left", "right")),
                 Rlevels = c("left", "right"), matchfun = ss_matchfun, report_p_vector = FALSE,
                 formula = exg_formula, trend = make_ssd_trend())
  set.seed(21)
  dat2 <- suppressMessages(make_data(p, des2, n_trials = 200, functions = list(SSD = stairS)))
  for (s in levels(dat2$subjects)) for (S in c("left", "right")) {
    st <- dat2[is.finite(dat2$SSD) & dat2$subjects == s & dat2$S == S, ]
    expect_equal(st$SSD[1], .3)
    d <- diff(st$SSD)
    expect_true(all(abs(d) < 1e-9 | abs(abs(d) - .05) < 1e-9))
    up <- is.na(st$R)[-nrow(st)]
    expect_true(all((d[d != 0] > 0) == up[d != 0]))
  }
})

test_that("tf(SSD) variant: pretransform on the probit scale, negative weight allowed", {
  des <- ss_design_tr(SSEXG, trend = make_ssd_trend("tf"), formula = exg_formula)
  expect_setequal(setdiff(names(sampled_pars(des)), names(exg_vals)), c("tf.k_sat", "tf.w"))
  p <- p_named(des, c(exg_vals, tf.k_sat = log(4), tf.w = -1))
  set.seed(22)
  dat <- make_data(p, des, n_trials = 60, functions = list(SSD = fixed_ssd),
                   return_trialwise_parameters = TRUE)
  tw <- attr(dat, "trialwise_parameters")
  tw1 <- tw[!duplicated(tw[, "trials"]), , drop = FALSE]
  expected <- pnorm(p[["tf"]] + ifelse(is.finite(dat$SSD), -1 * pmin(1, 4 * dat$SSD), 0))
  expect_equal(unname(tw1[, "tf"]), expected, tolerance = 1e-12)
  ll <- ss_ll_cr(dat, des, p)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-8)
  # combined muS + tf trends
  des2 <- ss_design_tr(SSEXG, trend = make_ssd_trend(c("muS", "tf")), formula = exg_formula)
  p2 <- p_named(des2, c(exg_vals, tr_vals, tf.k_sat = log(4), tf.w = -1))
  ll2 <- ss_ll_cr(dat, des2, p2)
  expect_equal(ll2[["cpp"]], ll2[["r"]], tolerance = 1e-8)
})

test_that("SSRDEX + muS(SSD): C++ == R (the trend transfers unchanged)", {
  des <- ss_design_tr(SSRDEX, formula = rdex_formula)
  p <- p_named(des, c(rdex_vals, tr_vals))
  set.seed(23)
  dat <- make_data(p, des, n_trials = 60, functions = list(SSD = fixed_ssd), TC = list(UC = 0.9))
  ll <- ss_ll_cr(dat, des, p)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-8)
})

test_that("trended SSEXG: init, predict and stop-signal plots run", {
  des <- ss_design_tr(SSEXG, formula = exg_formula)
  p <- p_named(des, c(exg_vals, tr_vals))
  set.seed(24)
  dat <- suppressMessages(make_data(p, des, n_trials = 120, functions = list(SSD = make_ssd(p_stop = .3))))
  emc <- make_emc(dat, des, type = "single", n_chains = 2, verbose = FALSE)
  emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 5), cores_for_chains = 1,
                 cores_per_chain = 1, verbose = FALSE)
  samples <- emc[[1]]$samples
  expect_named(samples$alpha[, 1, 1], names(p))
  expect_true(all(is.finite(samples$subj_ll)))
  pp <- predict(emc, n_post = 3, n_cores = 1)
  expect_true(all(c("SSD", "R", "rt") %in% names(pp)))
  expect_equal(nrow(pp), 3 * nrow(dat))
  expect_equal(pp$SSD[seq_len(nrow(dat))], dat$SSD)       # observed SSDs reused
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_ss_if(dat, post_predict = pp, factors = "S", probs = seq(0, 1, .5)))
  expect_no_error(plot_ss_srrt(dat, post_predict = pp, factors = "S", probs = seq(0, 1, .5)))
})
