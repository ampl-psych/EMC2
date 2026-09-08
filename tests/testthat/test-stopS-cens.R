## Stop-signal models on the censoring framework: C++ likelihood == R reference
## for every missingness code, compression invariance, and simulation coding.
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

ss_matchfun <- function(d) d$S == d$lR

ss_design <- function(model, ...) {
  design(model = model, factors = list(subjects = 1, S = c("left", "right")),
         Rlevels = c("left", "right"), matchfun = ss_matchfun, report_p_vector = FALSE, ...)
}

exg_formula <- list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
exg_p <- c(log(.5), .2, log(.05), log(.1), log(.2), log(.03), log(.05), qnorm(.05), qnorm(.1))
rdex_formula <- list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
rdex_p <- c(log(2), .5, log(1), log(.3), log(.15), log(.2), log(.03), log(.05), qnorm(.05), qnorm(.1))

p_vec <- function(des, values) {
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- values
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

fixed_ssd <- make_ssd(staircase = FALSE, values = c(.2, .3, .4))

test_that("SSEXG: simulated missingness codes follow the deadline convention", {
  des <- ss_design(SSEXG, formula = exg_formula)
  p <- p_vec(des, exg_p)
  set.seed(1)
  dat0 <- make_data(p, des, n_trials = 40, functions = list(SSD = fixed_ssd))
  expect_true("missingness" %in% names(dat0))
  # no deadline: intrinsic no-responses (stop success / go failure) are code 2
  # with UC = Inf; observed responses are NA
  expect_true(all(dat0$missingness[is.na(dat0$rt)] == 2L))
  expect_true(all(is.na(dat0$missingness[!is.na(dat0$rt)])))
  expect_true(all(is.infinite(dat0$UC)))
  # deadline: no-responses and too-slow responses are unified into code 2
  set.seed(1)
  datU <- make_data(p, des, n_trials = 40, functions = list(SSD = fixed_ssd), TC = list(UC = 0.9))
  expect_true(all(is.na(datU$rt[datU$missingness %in% 2L])))
  expect_true(all(datU$rt[is.na(datU$missingness)] <= 0.9))
  expect_true(all(datU$UC == 0.9))
  # design-level TC is honoured by make_data
  desU <- ss_design(SSEXG, formula = exg_formula, TC = list(UC = 0.9))
  set.seed(1)
  datD <- make_data(p, desU, n_trials = 40, functions = list(SSD = fixed_ssd))
  expect_true(any(datD$missingness %in% 2L))
  expect_true(all(datD$UC == 0.9))
})

test_that("SSEXG: C++ == R for uncensored, deadline (UC) and lower (LC) censoring", {
  des <- ss_design(SSEXG, formula = exg_formula)
  p <- p_vec(des, exg_p)
  set.seed(2)
  dat0 <- make_data(p, des, n_trials = 40, functions = list(SSD = fixed_ssd))
  ll <- ss_ll_cr(dat0, des, p)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-8)
  set.seed(2)
  datU <- make_data(p, des, n_trials = 40, functions = list(SSD = fixed_ssd), TC = list(UC = 0.9))
  ll <- ss_ll_cr(datU, des, p)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-8)
  set.seed(2)
  datL <- make_data(p, des, n_trials = 60, functions = list(SSD = fixed_ssd), TC = list(LC = 0.45, UC = 0.9))
  expect_true(any(datL$missingness %in% 1L))
  ll <- ss_ll_cr(datL, des, p)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-6)  # numerical integration
})

test_that("SSRDEX: C++ == R for uncensored, deadline (UC) and lower (LC) censoring", {
  des <- ss_design(SSRDEX, formula = rdex_formula)
  p <- p_vec(des, rdex_p)
  set.seed(3)
  dat0 <- make_data(p, des, n_trials = 40, functions = list(SSD = fixed_ssd))
  ll <- ss_ll_cr(dat0, des, p)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-8)
  set.seed(3)
  datU <- make_data(p, des, n_trials = 40, functions = list(SSD = fixed_ssd), TC = list(UC = 0.9))
  ll <- ss_ll_cr(datU, des, p)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-8)
  set.seed(3)
  datL <- make_data(p, des, n_trials = 40, functions = list(SSD = fixed_ssd), TC = list(LC = 0.35, UC = 0.9))
  expect_true(any(datL$missingness %in% 1L))
  ll <- ss_ll_cr(datL, des, p)
  expect_equal(ll[["cpp"]], ll[["r"]], tolerance = 1e-6)
})

test_that("SSEXG: staircase SSDs survive into the dadm and compression is invariant", {
  des <- ss_design(SSEXG, formula = exg_formula)
  p <- p_vec(des, exg_p)
  set.seed(4)
  dat <- make_data(p, des, n_trials = 60, functions = list(SSD = make_ssd()), TC = list(UC = 0.9))
  ssd <- dat$SSD[is.finite(dat$SSD)]
  expect_gt(length(unique(ssd)), 1)             # the ladder moved
  expect_true(all(ssd >= 0 & ssd < 0.9))
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE, verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  expect_equal(dadm$SSD[dadm$lR == levels(dadm$lR)[1]], dat$SSD)
  llF <- ss_ll_cr(dat, des, p, compress = FALSE)
  llT <- ss_ll_cr(dat, des, p, compress = TRUE)
  expect_equal(llF[["cpp"]], llT[["cpp"]], tolerance = 1e-10)
  expect_equal(llF[["cpp"]], llF[["r"]], tolerance = 1e-8)
})

test_that("censored RDM (non-SS) likelihood is compression invariant and matches R", {
  des <- ss_design(RDM, formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1))
  p <- p_vec(des, c(log(3), log(1), log(1), log(0.5), log(0.2)))
  set.seed(5)
  dat <- make_data(p, des, n_trials = 60)
  mdat <- make_missing(dat, LC = 0.4, UC = 0.85, rt_resolution = 1/60)
  expect_true(any(mdat$missingness %in% 1L) && any(mdat$missingness %in% 2L))
  ll_of <- function(compress) {
    emc <- make_emc(mdat, des, type = "single", rt_resolution = 1/60, n_chains = 1,
                    compress = compress, verbose = FALSE)
    dadm <- emc[[1]]$data[[1]]
    pm <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))
    list(ll = as.numeric(calc_ll_manager(pm, dadm, des$model)), dadm = dadm)
  }
  llF <- ll_of(FALSE); llT <- ll_of(TRUE)
  expect_equal(llF$ll, llT$ll, tolerance = 1e-10)
  # hand computation: observed = pdf(winner) * S(losers); code 1 = 1 - S_race(LC);
  # code 2 = S_race(UC)
  dadm <- llF$dadm
  pars <- get_pars_matrix_oo(p, dadm, des$model())
  nacc <- 2; n_tr <- nrow(dadm) / nacc
  ll_r <- numeric(n_tr)
  for (t in seq_len(n_tr)) {
    rows <- ((t - 1) * nacc + 1):(t * nacc)
    P <- pars[rows, , drop = FALSE]
    miss <- dadm$missingness[rows[1]]
    if (is.na(miss)) {
      w <- dadm$winner[rows]
      ll_r[t] <- log(dRDM(dadm$rt[rows[1]], P[w, , drop = FALSE])) +
        sum(log(1 - pRDM(dadm$rt[rows[1]], P[!w, , drop = FALSE])))
    } else if (miss == 2L) {
      ll_r[t] <- log(prod(1 - pRDM(rep(dadm$UC[rows[1]], nacc), P)))
    } else if (miss == 1L) {
      ll_r[t] <- log(1 - prod(1 - pRDM(rep(dadm$LC[rows[1]], nacc), P)))
    }
  }
  expect_equal(sum(pmax(ll_r, log(1e-10))), llF$ll, tolerance = 1e-10)
})

test_that("real-data censoring routes: make_missing() then design(data=), or design(TC=)", {
  des <- ss_design(SSEXG, formula = exg_formula)
  p <- p_vec(des, exg_p)
  set.seed(6)
  raw <- make_data(p, des, n_trials = 40, functions = list(SSD = fixed_ssd))
  raw <- raw[, c("subjects", "trials", "S", "R", "SSD", "rt")]   # like real data
  # route 1: censor the data first; censoring columns must not become covariates
  dat <- make_missing(raw, UC = 0.9, rt_resolution = NULL)
  des1 <- design(model = SSEXG, data = dat, matchfun = ss_matchfun, formula = exg_formula,
                 report_p_vector = FALSE)
  expect_false(any(c("LT", "UT", "LC", "UC", "missingness") %in% des1$Fcovariates))
  emc1 <- make_emc(dat, des1, type = "single", n_chains = 1, verbose = FALSE)
  d1 <- emc1[[1]]$data[[1]]
  expect_true(all(c("UC", "missingness") %in% names(d1)))
  expect_true(any(d1$missingness %in% 2L))
  # route 2: TC in the design, raw data -> make_emc applies it
  des2 <- design(model = SSEXG, data = raw, matchfun = ss_matchfun, formula = exg_formula,
                 TC = list(UC = 0.9), report_p_vector = FALSE)
  expect_message(emc2 <- make_emc(raw, des2, type = "single", n_chains = 1, verbose = FALSE),
                 "truncation/censoring")
  d2 <- emc2[[1]]$data[[1]]
  expect_true(all(d2$UC == 0.9))
  expect_true(any(d2$missingness %in% 2L))
  # a truncation window is still refused for stop-signal models
  datT <- make_missing(raw, UT = 0.9, rt_resolution = NULL)
  desT <- design(model = SSEXG, data = datT, matchfun = ss_matchfun, formula = exg_formula,
                 report_p_vector = FALSE)
  expect_error(make_emc(datT, desT, type = "single", n_chains = 1, verbose = FALSE),
               "Truncation")
})

test_that("multithreaded backend gives identical stop-signal likelihoods", {
  des <- ss_design(SSEXG, formula = exg_formula)
  p <- p_vec(des, exg_p)
  set.seed(7)
  dat <- make_data(p, des, n_trials = 60, functions = list(SSD = make_ssd()), TC = list(UC = 0.9))
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE, verbose = FALSE)
  dadm <- emc[[1]]$data[[1]]
  pm <- rbind(p, matrix(rnorm(20 * length(p), rep(p, each = 20), .3), 20))
  colnames(pm) <- names(p)
  ll1 <- as.numeric(calc_ll_manager(pm, dadm, des$model))
  old <- options(emc.ll_backend = "multithreaded", emc.n_threads = 2)
  on.exit(options(old))
  ll2 <- suppressWarnings(as.numeric(calc_ll_manager(pm, dadm, des$model)))  # warns if no OpenMP
  expect_equal(ll1, ll2, tolerance = 0)
  # trialwise too
  tw1 <- as.numeric(calc_ll_manager(pm[1, , drop = FALSE], dadm, des$model, return_trialwise = TRUE))
  options(emc.ll_backend = "multiprocess")
  tw0 <- as.numeric(calc_ll_manager(pm[1, , drop = FALSE], dadm, des$model, return_trialwise = TRUE))
  expect_equal(tw1, tw0, tolerance = 0)
  expect_equal(sum(tw0), ll1[1], tolerance = 1e-10)
})
