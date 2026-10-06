## make_stop_data(): latent go-race winner (goR, goRT) and stop finishing time
## (SSRT) for SSEXG / SSRDEX. Checks the latent times against the likelihood's
## distributions, their exact consistency with the simulated R/rt, failure
## rates, censoring/truncation, trends (both simulation paths) and that R/rt
## are unchanged by the latent option.
RNGkind("L'Ecuyer-CMRG")
set.seed(4321)

ss_mf <- function(d) as.numeric(d$S) == as.numeric(d$lR)
fixed_ssd <- function(d) rep(c(Inf, .2, Inf, .5), length.out = nrow(d))

des_exg <- design(model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
                  Rlevels = c("left", "right"), matchfun = ss_mf, report_p_vector = FALSE,
                  formula = list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1),
                  constants = c(exg_lb = .05, exgS_lb = .15))
p_exg <- c(mu = log(.6), mu_lMTRUE = log(.8), sigma = log(.05), tau = log(.2), muS = log(.2),
           sigmaS = log(.03), tauS = log(.05), gf = qnorm(.1), tf = qnorm(.1))

des_rdex <- design(model = SSRDEX, factors = list(subjects = 1, S = c("left", "right")),
                   Rlevels = c("left", "right"), matchfun = ss_mf, report_p_vector = FALSE,
                   formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
p_rdex <- c(v = log(2), v_lMTRUE = .5, B = log(1), A = log(.3), t0 = log(.15), muS = log(.12),
            sigmaS = log(.03), tauS = log(.05), gf = qnorm(.05), tf = qnorm(.1))

# Exact relations between latent times and the simulated response
expect_ss_consistent <- function(d, go_levels = levels(d$goR)) {
  stop_trial <- is.finite(d$SSD)
  resp_go <- !is.na(d$R) & as.character(d$R) %in% go_levels
  go_resp <- !stop_trial & resp_go
  expect_equal(d$rt[go_resp], d$goRT[go_resp])
  expect_equal(as.character(d$R[go_resp]), as.character(d$goR[go_resp]))
  expect_true(all(is.na(d$rt[!stop_trial & !resp_go]) | !as.character(d$R[!stop_trial & !resp_go]) %in% go_levels))
  # signal-respond: the go winner beat the stop racer
  sr <- stop_trial & resp_go
  expect_equal(d$rt[sr], d$goRT[sr])
  expect_true(all(d$goRT[sr] < d$SSD[sr] + d$SSRT[sr]))
  # no go response on a stop trial: go failed or the stop racer won
  ns <- stop_trial & !resp_go
  expect_true(all(is.infinite(d$goRT[ns]) | d$goRT[ns] > d$SSD[ns] + d$SSRT[ns]))
  expect_true(all(is.na(d$SSRT[!stop_trial])))
}

test_that("SSEXG latent times follow the truncated ex-Gaussian race", {
  d <- make_stop_data(p_exg, des_exg, n_trials = 3000, functions = list(SSD = fixed_ssd),
                      censor_go = FALSE)
  expect_true(all(c("goR", "goRT", "SSRT") %in% names(d)))
  expect_ss_consistent(d)
  st <- is.finite(d$SSD)
  ssrt <- d$SSRT[st & is.finite(d$SSRT)]
  expect_gt(min(ssrt), .15)          # stop lower bound exgS_lb, from stop onset
  expect_gt(suppressWarnings(ks.test(ssrt, function(q) pTEXG_RDEX(q, .2, .03, .05, .15)))$p.value, .001)
  # go winner, per stimulus: 1 - prod(1 - F) over the two go racers
  for (s in c("left", "right")) {
    g <- d$goRT[d$S == s & is.finite(d$goRT)]
    expect_gt(min(g), .05)           # go lower bound exg_lb
    Fg <- function(q) 1 - (1 - pTEXG_RDEX(q, .6 * .8, .05, .2, .05)) * (1 - pTEXG_RDEX(q, .6, .05, .2, .05))
    expect_gt(suppressWarnings(ks.test(g, Fg))$p.value, .001)
  }
  # trigger and go failure rates
  expect_equal(mean(is.infinite(d$SSRT[st])), .1, tolerance = .03 / .1)
  expect_equal(mean(is.infinite(d$goRT)), .1, tolerance = .03 / .1)
  expect_true(all(is.na(d$goR[is.infinite(d$goRT)])))
})

test_that("SSRDEX latent times follow the Wald go race and ex-Gaussian stop", {
  d <- make_stop_data(p_rdex, des_rdex, n_trials = 3000, functions = list(SSD = fixed_ssd),
                      censor_go = FALSE)
  expect_ss_consistent(d)
  ok <- is.finite(d$goRT)
  Fg <- function(q) 1 - (1 - pWald_RDEX(q, exp(log(2) + .5), 1, .3, .15, 1)) *
    (1 - pWald_RDEX(q, 2, 1, .3, .15, 1))
  expect_gt(suppressWarnings(ks.test(d$goRT[ok], Fg))$p.value, .001)
  st <- is.finite(d$SSD) & is.finite(d$SSRT)
  expect_gt(min(d$SSRT[st]), .05)
  expect_gt(suppressWarnings(ks.test(d$SSRT[st], function(q) pTEXG_RDEX(q, .12, .03, .05, .05)))$p.value, .001)
  expect_equal(mean(is.infinite(d$goRT)), .05, tolerance = .02 / .05)
})

test_that("stop-triggered accumulators are not part of the go race", {
  des_st <- design(model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
                   Rlevels = c("left", "right", "st"), report_p_vector = FALSE,
                   matchfun = function(d) as.character(d$S) == as.character(d$lR),
                   functions = list(lI = function(d) factor(ifelse(d$lR == "st", 1, 2), levels = 1:2)),
                   formula = list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  d <- make_stop_data(p_exg, des_st, n_trials = 500, functions = list(SSD = fixed_ssd),
                      censor_go = FALSE)
  expect_true(any(d$R == "st", na.rm = TRUE))
  expect_false(any(d$goR == "st", na.rm = TRUE))
  expect_ss_consistent(d, go_levels = c("left", "right"))
})

test_that("latent option leaves R, rt and SSD unchanged (staircase SSDs)", {
  set.seed(11)
  a <- make_data(p_exg, des_exg, n_trials = 200, functions = list(SSD = make_ssd()))
  set.seed(11)
  b <- make_stop_data(p_exg, des_exg, n_trials = 200, functions = list(SSD = make_ssd()),
                      stop_on_go_trials = TRUE)
  expect_identical(a$R, b$R)
  expect_identical(a$rt, b$rt)
  expect_identical(a$SSD, b$SSD)
  expect_identical(a$missingness, b$missingness)
  expect_false(any(c("goR", "goRT", "SSRT", "goMissingness") %in% names(a)))
  expect_ss_consistent(transform(b, SSRT = ifelse(is.finite(SSD), SSRT, NA)))
  # counterfactual stop racers on go trials: bounded below, trigger failures Inf
  go_ssrt <- b$SSRT[!is.finite(b$SSD)]
  expect_false(anyNA(go_ssrt))
  expect_gt(min(go_ssrt), .15)
  expect_true(any(is.infinite(go_ssrt)))
})

test_that("go race is censored like rt; truncation drops the same trials", {
  TC <- list(LC = .4, UC = .7)
  set.seed(5)
  raw <- make_stop_data(p_exg, des_exg, n_trials = 500, functions = list(SSD = fixed_ssd),
                        TC = TC, censor_go = FALSE)
  set.seed(5)
  cen <- make_stop_data(p_exg, des_exg, n_trials = 500, functions = list(SSD = fixed_ssd), TC = TC)
  expect_identical(raw$rt, cen$rt)
  lo <- raw$goRT < .4
  hi <- raw$goRT > .7
  expect_true(any(lo) && any(hi) && any(is.infinite(raw$goRT)))
  expect_true(all(cen$goMissingness[lo] == 1L))
  expect_true(all(cen$goMissingness[hi] == 2L))
  expect_true(all(is.na(cen$goMissingness[!lo & !hi])))
  expect_true(all(is.na(cen$goRT[lo | hi])) && all(is.na(cen$goR[lo | hi])))
  expect_equal(cen$goRT[!lo & !hi], raw$goRT[!lo & !hi])
  # censored go RTs on go trials coincide with censored rt
  go <- !is.finite(cen$SSD)
  expect_identical(cen$missingness[go], cen$goMissingness[go])
  # go failures stay Inf without a finite UC
  set.seed(5)
  lc_only <- make_stop_data(p_exg, des_exg, n_trials = 500, functions = list(SSD = fixed_ssd),
                            TC = list(LC = .4))
  expect_true(any(is.infinite(lc_only$goRT)))
  # LT/UT: same trials removed as for the observed data
  set.seed(6)
  a <- make_data(p_exg, des_exg, n_trials = 500, functions = list(SSD = fixed_ssd),
                 TC = list(LT = .3, UT = .9))
  set.seed(6)
  b <- make_stop_data(p_exg, des_exg, n_trials = 500, functions = list(SSD = fixed_ssd),
                      TC = list(LT = .3, UT = .9))
  expect_identical(a$trials, b$trials)
  expect_identical(a$rt, b$rt)
})

test_that("trends on stop and go parameters carry into the latent times", {
  # stop: muS rises with SSD (saturating linear kernel)
  trS <- make_trend(make_base("muS", "lin", make_kernel("SSD", "slin_incr"), phase = "posttransform"))
  desS <- design(model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
                 Rlevels = c("left", "right"), matchfun = ss_mf, report_p_vector = FALSE,
                 trend = trS, transform = list(func = c(muS.w = "exp")),
                 formula = list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  pS <- sampled_pars(desS, doMap = FALSE)
  pS[] <- c(p_exg, muS.k_sat = log(5), muS.w = log(.2))[names(pS)]
  d <- make_stop_data(pS, desS, n_trials = 2000, functions = list(SSD = fixed_ssd), censor_go = FALSE)
  expect_ss_consistent(d)
  s2 <- d$SSRT[d$SSD == .2 & is.finite(d$SSRT)]
  s5 <- d$SSRT[d$SSD == .5 & is.finite(d$SSRT)]
  # muS(SSD) = .2 + .2 * min(1, 5 * SSD): .4 at both SSDs here, so compare to no trend
  expect_gt(mean(s2), .2 + .05 + .15)
  expect_gt(suppressWarnings(ks.test(s5, function(q) pTEXG_RDEX(q, .4, .03, .05, .05)))$p.value, .001)

  # go: mu rises with a covariate
  trG <- make_trend(make_base("mu", "lin", make_kernel("x", "lin_incr")))
  desG <- design(model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
                 Rlevels = c("left", "right"), matchfun = ss_mf, report_p_vector = FALSE,
                 covariates = "x", trend = trG,
                 formula = list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  pG <- sampled_pars(desG, doMap = FALSE)
  pG[] <- c(p_exg, mu.w = .5)[names(pG)]
  x <- rep(0:1, length.out = 4000)
  dG <- make_stop_data(pG, desG, n_trials = 2000, functions = list(SSD = fixed_ssd),
                       covariates = data.frame(x = x), censor_go = FALSE)
  expect_ss_consistent(dG)
  m <- tapply(dG$goRT[is.finite(dG$goRT)], dG$x[is.finite(dG$goRT)], mean)
  expect_gt(m[["1"]] - m[["0"]], .1)
})

test_that("make_stop_data on an emc object (conditional and trial-by-trial paths)", {
  des <- design(model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
                Rlevels = c("left", "right"), matchfun = ss_mf, report_p_vector = FALSE,
                functions = list(SSD = make_ssd(p_stop = .3)),
                formula = list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  set.seed(24)
  dat <- suppressMessages(make_data(p_exg, des, n_trials = 100))
  emc <- make_emc(dat, des, type = "single", n_chains = 2, verbose = FALSE)
  emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 5), cores_for_chains = 1,
                 cores_per_chain = 1, verbose = FALSE)
  pp <- make_stop_data(emc, n_post = 3, conditional_on_data = TRUE, censor_go = FALSE)
  expect_equal(nrow(pp), 3 * nrow(dat))
  expect_true(all(c("postn", "goR", "goRT", "SSRT") %in% names(pp)))
  expect_equal(pp$SSD[seq_len(nrow(dat))], dat$SSD)
  for (i in unique(pp$postn)) expect_ss_consistent(pp[pp$postn == i, ])
  # staircase re-run trial by trial
  ppu <- suppressMessages(make_stop_data(emc, n_post = 2, conditional_on_data = FALSE,
                                         censor_go = FALSE))
  expect_true(all(c("goR", "goRT", "SSRT") %in% names(ppu)))
  expect_false(isTRUE(all.equal(ppu$SSD[seq_len(nrow(dat))], dat$SSD)))
  expect_true(is.factor(ppu$goR))
  for (i in unique(ppu$postn)) expect_ss_consistent(ppu[ppu$postn == i, ])
})

test_that("make_stop_data rejects non stop-signal models", {
  des <- design(model = LNR, factors = list(subjects = 1, S = 1:2), Rlevels = 1:2,
                matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                report_p_vector = FALSE)
  expect_error(make_stop_data(c(m = 0, m_lMTRUE = 1, s = 0, t0 = log(.2)), des, n_trials = 2),
               "stop-signal model")
})
