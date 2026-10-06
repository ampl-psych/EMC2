# The interweaving sweep (scale_move_standard, R/variant_standard.R). The long statistical checks of its
# stationary distribution run only with EMC2_SLOW_TESTS=true.

# log r = sum_s dl_s + log_prior(delta): the change in the group density, the covariance prior (with its
# auxiliary a_j), the prior of a_j and the Jacobian of the map
# (alpha_.j, Sigma, a_j) -> (mu_j + f (alpha_.j - mu_j), D Sigma D, a_j / f^2),
# for a parameter in a covariance block and for one with its own variance.
test_that("scale_move_log_prior matches the brute-force ratio", {
  set.seed(5)
  v <- 2; A <- c(.3, .5, 1)
  brute <- function(alpha, mu, Sigma, a, j, delta, blocked) {
    n <- nrow(alpha); p <- ncol(alpha); f <- exp(delta)
    D <- rep(1, p); D[j] <- f
    alpha2 <- alpha; alpha2[, j] <- mu[j] + f * (alpha[, j] - mu[j])
    Sigma2 <- Sigma * outer(D, D)
    a2 <- a; a2[j] <- a[j] / f^2
    lprior <- function(S, a) (if (blocked) {
      log(EMC2:::robust_diwish(S, v + p - 1, 2 * v * diag(1 / a, p)))
    } else sum(EMC2:::logdinvGamma(diag(S), v / 2, v / a))) +
      sum(EMC2:::logdinvGamma(a, 1 / 2, 1 / A^2))
    group <- function(al, S) sum(mvtnorm::dmvnorm(al, mu, S, log = TRUE))
    # Jacobian: n deviations, the free elements of the covariance block, a_j
    log_jac <- n * delta + (if (blocked) (p + 1) * delta else 2 * delta) - 2 * delta
    lprior(Sigma2, a2) - lprior(Sigma, a) + group(alpha2, Sigma2) - group(alpha, Sigma) + log_jac
  }
  for (blocked in c(TRUE, FALSE)) {
    p <- 3; n <- 4
    Sigma <- if (blocked) crossprod(matrix(rnorm(9), 3)) + diag(3) * .2 else diag(exp(rnorm(3)))
    a <- 1 / rgamma(p, 1 / 2, rate = 1 / A^2)
    mu <- rnorm(p); alpha <- mvtnorm::rmvnorm(n, mu, Sigma)
    for (j in 1:p) for (delta in c(-.7, .2, 1.1)) {
      got <- EMC2:::scale_move_log_prior(delta, A[j], a[j])
      expect_equal(got, brute(alpha, mu, Sigma, a, j, delta, blocked), tolerance = 1e-8)
    }
  }
})

test_that("the sweep runs in adapt and sample unless switched off, and never under the legacy sampler", {
  expect_identical(EMC2:::scale_move_stages(), c("adapt", "sample"))
  withr::with_options(list(emc.scale_move = FALSE), expect_length(EMC2:::scale_move_stages(), 0))
  withr::with_options(list(emc.sampler = "legacy"), expect_length(EMC2:::scale_move_stages(), 0))
})

test_that("scale_move_gate drops only the useless moves of a sweep that cannot be split further", {
  pn <- c("a", "b", "c")
  mk <- function(blocks, step, loc, acc = 10, n = 100) {
    s <- EMC2:::scale_move_init(NULL, pn); s$blocks <- blocks
    s$n_in[] <- n; s$acc_real[] <- acc; s$acc_loc_real[] <- acc; s$step[] <- step; s$loc[] <- loc
    s
  }
  expect_true(all(mk(3, 1, 1)$active_scale) && all(mk(3, 1, 1)$active_loc))
  # 200 subjects. a, b: the posterior variance of the log group SD is the Gibbs step's own (1 / 400);
  # c: a funnel, a hundred times wider. Group means: posterior variance = group variance / n for all.
  var_lsd <- c(.0025, .0025, .25); var_mu <- rep(4e-4, 3); mean_var <- rep(.08, 3)
  gate <- function(sm, decide = TRUE, vm = var_mu) EMC2:::scale_move_gate(sm, var_lsd, vm, mean_var, 200, decide)
  # realised acceptance .1; u_scale = .1 step^2 / var_lsd = 1.6, 1.6e-4, 4e-3 against kappa x Gibbs = .1, .1, .001
  #                         u_loc   = .1 loc^2 / var_mu   = .4, 6e-5, 6e-5    against .1
  chain <- function(blocks) mk(blocks, step = c(.2, .002, .1), loc = c(.04, .0005, .0005))
  # recorded only
  out <- gate(list(chain(3), chain(3)), decide = FALSE)
  expect_true(all(out[[1]]$active_scale) && all(out[[1]]$active_loc))
  expect_identical(out[[2]]$mark$n_in, out[[2]]$n_in)
  # one chain's sweep can still be split: nothing is dropped
  out <- gate(list(chain(3), chain(2)))
  expect_true(all(out[[1]]$active_scale) && all(out[[1]]$active_loc))
  # as many blocks as parameters in both chains: b loses both moves, c (the funnel) keeps its scale move
  out <- gate(list(chain(3), chain(3)))
  expect_identical(unname(out[[1]]$active_scale), c(TRUE, FALSE, TRUE))
  expect_identical(unname(out[[1]]$active_loc), c(TRUE, FALSE, FALSE))
  expect_identical(out[[1]]$active_scale, out[[2]]$active_scale)
  # a dropped move stays dropped whatever its counters do later
  again <- lapply(out, function(x) { x$n_in[] <- 200; x$acc_real[] <- 60; x$acc_loc_real[] <- 60; x$blocks <- 2; x$step[] <- 1; x$loc[] <- 1; x })
  out2 <- gate(again)
  expect_identical(unname(out2[[1]]$active_scale), c(TRUE, FALSE, TRUE))
  expect_identical(unname(out2[[1]]$active_loc), c(TRUE, FALSE, FALSE))
  # no location moves (a group design): none is active, the scale moves are judged as before
  out <- gate(list(chain(3), chain(3)), vm = NULL)
  expect_false(any(out[[1]]$active_loc))
  expect_identical(unname(out[[1]]$active_scale), c(TRUE, FALSE, TRUE))
})

# ---- a 2-subject toy -------------------------------------------------------
# Two subjects, two parameters, a heavy-tailed likelihood on sufficient statistics, so the quadratic
# surrogate the sweep uses is wrong and the delayed-acceptance step has to correct it.
toy_ll <- function(pars, dadm, ...) -sum(dadm$T * log1p((dadm$ybar - pars)^2))
toy_pars <- c("a", "b")
toy_T <- c(30, 3)
toy_dat <- withr::with_seed(21, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                            .rng_sample_kind = "Rejection", do.call(rbind, lapply(1:2, function(s) data.frame(
  subjects = s, par = toy_pars, T = toy_T, ybar = c(.4, -.3) * (s - 1.5) * 2 + rnorm(2, 0, .2)))))
toy_dat$subjects <- factor(toy_dat$subjects)
toy_design <- design(model = toy_ll, custom_p_vector = toy_pars, report_p_vector = FALSE)
toy_emc <- make_emc(toy_dat, toy_design, type = "standard", n_chains = 2, compress = FALSE)

toy_subject_ll <- function(sampler, s, x)
  as.numeric(EMC2:::calc_ll_manager(matrix(x, 1, dimnames = list(NULL, toy_pars)), sampler$data[[s]], sampler$model))

# A state of the toy and its surrogate: "rough" takes the curvature of the Gaussian part at the
# likelihood's mode (overstating its tails); "flat" has none, so only the exact check knows the likelihood.
toy_state <- function(sampler, surrogate = c("rough", "flat"), alpha = matrix(c(.3, -.2, -.3, .25), 2, 2)) {
  surrogate <- match.arg(surrogate)
  p <- 2; n <- 2; mu <- c(.1, -.1); a <- c(.5, 2); Sigma <- diag(c(.3, .3))
  dimnames(alpha) <- list(toy_pars, NULL)
  lik_prec <- lapply(1:n, function(s) {
    d <- sampler$data[[s]]
    if (surrogate == "flat") list(prec = diag(0, 2), lin = rep(0, 2)) else list(prec = diag(2 * d$T), lin = 2 * d$T * d$ybar)
  })
  list(pars = list(tmu = mu, tvar = Sigma, tvinv = solve(Sigma), a_half = a, subj_mu = matrix(mu, p, n), alpha = alpha),
       alpha = alpha, ll = sapply(1:n, function(s) toy_subject_ll(sampler, s, alpha[, s])), lik_prec = lik_prec)
}

test_that("the sweep keeps its state consistent whether its blocks pass or fail", {
  sampler <- toy_emc[[1]]
  for (surrogate in c("rough", "flat")) {
    set.seed(8)
    x <- toy_state(sampler, surrogate)
    st <- EMC2:::scale_move_init(NULL, toy_pars)
    for (i in 1:50) {
      sm <- EMC2:::scale_move_standard(sampler, x$pars, x$alpha, x$ll, st, x$lik_prec)
      expect_equal(sm$ll, sapply(1:2, function(s) toy_subject_ll(sampler, s, sm$alpha[, s])))
      expect_identical(sm$pars$alpha, sm$alpha)
      expect_equal(sm$pars$subj_mu, matrix(sm$pars$tmu, 2, 2))
      expect_equal(sm$pars$tvinv %*% sm$pars$tvar, diag(2))
      x$pars <- sm$pars; x$alpha <- sm$alpha; x$ll <- sm$ll; st <- sm$settings
    }
    # the surrogate's acceptance counts every move the exact check then let through, and more
    expect_true(all(st$acc_in >= st$acc_real) && all(st$acc_loc >= st$acc_loc_real))
    expect_gt(st$acc_out, 0); expect_lt(st$acc_out, st$n_out)
  }
})

# The reference: exact draws of mu | alpha, Sigma, of Sigma | alpha, mu, a, of a | Sigma, and a random-walk
# step for alpha | mu, Sigma; adding the sweep must not change the stationary distribution.
toy_run <- function(sampler, iter, use_move, seed, surrogate = c("rough", "flat")) {
  set.seed(seed)
  x <- toy_state(sampler, match.arg(surrogate), alpha = matrix(0, 2, 2))
  p <- 2; n <- 2; v <- sampler$prior$v; A <- sampler$prior$A
  m0 <- sampler$prior$theta_mu_mean; V0inv <- sampler$prior$theta_mu_invar
  mu <- x$pars$tmu; a <- x$pars$a_half; Sigma <- x$pars$tvar; alpha <- x$alpha; cur_ll <- x$ll
  settings <- EMC2:::scale_move_init(NULL, toy_pars)
  out <- matrix(NA_real_, iter, 9, dimnames = list(NULL, c("s11", "s22", "r12", "a11", "a22", "aux1", "aux2", "mu1", "mu2")))
  for (i in seq_len(iter)) {
    # mu | alpha, Sigma; Sigma | alpha, mu, a ~ IW(v + p - 1 + n, 2 v diag(1/a) + S); a | Sigma
    Sinv <- solve(Sigma); Q <- V0inv + n * Sinv
    mu <- drop(mvtnorm::rmvnorm(1, solve(Q, V0inv %*% m0 + Sinv %*% rowSums(alpha)), solve(Q)))
    r <- alpha - mu
    Sigma <- EMC2:::riwish(v + p - 1 + n, 2 * v * diag(1 / a) + r %*% t(r))
    Sinv <- solve(Sigma)
    a <- 1 / rgamma(p, (v + p) / 2, rate = v * diag(Sinv) + 1 / A^2)
    for (s in 1:n) {
      prop <- alpha[, s] + rnorm(p, 0, .3)
      lp <- toy_subject_ll(sampler, s, prop) + mvtnorm::dmvnorm(prop, mu, Sigma, log = TRUE)
      lc <- cur_ll[s] + mvtnorm::dmvnorm(alpha[, s], mu, Sigma, log = TRUE)
      if (log(runif(1)) < lp - lc) { alpha[, s] <- prop; cur_ll[s] <- toy_subject_ll(sampler, s, prop) }
    }
    if (use_move) {
      pars <- list(tmu = mu, tvar = Sigma, tvinv = Sinv, a_half = a, subj_mu = matrix(mu, p, n), alpha = alpha)
      sm <- EMC2:::scale_move_standard(sampler, pars, alpha, cur_ll, settings, x$lik_prec, frozen = i > iter / 4)
      Sigma <- sm$pars$tvar; a <- sm$pars$a_half; mu <- sm$pars$tmu; alpha <- sm$alpha; cur_ll <- sm$ll
      settings <- sm$settings
    }
    out[i, ] <- c(Sigma[1, 1], Sigma[2, 2], Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2]), alpha[1, 1], alpha[2, 2], a, mu)
  }
  list(draws = out[-(1:(iter / 4)), ], settings = settings)
}

expect_same_posterior <- function(mv, ref) {
  ess <- function(x) coda::effectiveSize(coda::mcmc(x))
  for (k in colnames(ref$draws)) {
    lg <- k %in% c("s11", "s22", "aux1", "aux2")
    x <- if (lg) log(mv$draws[, k]) else mv$draws[, k]
    y <- if (lg) log(ref$draws[, k]) else ref$draws[, k]
    expect_lt(abs((mean(x) - mean(y)) / sqrt(var(x) / ess(x) + var(y) / ess(y))), 4)
    expect_lt(abs(log(sd(x) / sd(y))), .25)
  }
}

test_that("the sweep leaves the toy's posterior invariant, with a rough and with a flat surrogate", {
  skip_on_cran(); skip_if_not(slow_tests(), "EMC2_SLOW_TESTS is not true")
  sampler <- toy_emc[[1]]
  iter <- 12000
  ref <- toy_run(sampler, iter, use_move = FALSE, seed = 1)
  # rough surrogate: most blocks pass, and the realised acceptance of the scale and location moves sits
  # near target x block acceptance
  mv <- toy_run(sampler, iter, use_move = TRUE, seed = 2)
  st <- mv$settings
  expect_gt(st$n_out, iter / 4); expect_gt(st$acc_out / st$n_out, .2)
  expect_gt(st$r_block, .6); expect_lt(st$r_block, 1)
  for (rate in list(st$acc_real / st$n_in, st$acc_loc_real / st$n_in)) expect_true(all(rate > .2 & rate < .4))
  expect_true(all(st$acc_in / st$n_in > .22 & st$acc_in / st$n_in < .45))
  expect_same_posterior(mv, ref)
  # flat surrogate: only the exact check knows the likelihood; the realised acceptance is held near
  # target x floor and the steps stay bounded
  mv <- toy_run(sampler, iter, use_move = TRUE, seed = 4, surrogate = "flat")
  st <- mv$settings
  expect_gt(sum(st$acc_in), 1.5 * sum(st$acc_real))
  expect_lt(st$acc_out / st$n_out, EMC2:::scale_move_r_floor); expect_lt(st$r_block, .5)
  target <- EMC2:::scale_move_target * EMC2:::scale_move_r_floor
  for (rate in list(st$acc_real / st$n_in, st$acc_loc_real / st$n_in)) expect_true(all(rate > .6 * target & rate < 1.6 * target))
  expect_true(all(st$step < 15) && all(st$loc < 15))
  expect_same_posterior(mv, ref)
})

test_that("the sweep leaves a gated parameter alone", {
  sampler <- toy_emc[[1]]
  set.seed(3)
  x <- toy_state(sampler)
  settings <- EMC2:::scale_move_init(NULL, toy_pars); settings$iter <- 1; settings$blocks <- 2
  settings$step[] <- .3; settings$loc[] <- .1
  # parameter b out of the sweep: its group variance, mean, subjects and step never change
  settings$active_scale[] <- c(TRUE, FALSE); settings$active_loc[] <- c(TRUE, FALSE)
  moved <- FALSE; same <- TRUE; st <- settings; pr <- x$pars; al <- x$alpha; l <- x$ll
  for (i in 1:50) {
    sm <- EMC2:::scale_move_standard(sampler, pr, al, l, st, x$lik_prec)
    pr <- sm$pars; al <- sm$alpha; l <- sm$ll; st <- sm$settings
    same <- same && identical(pr$tvar[2, 2], x$pars$tvar[2, 2]) && identical(pr$tmu[2], x$pars$tmu[2]) &&
      identical(al[2, ], x$alpha[2, ])
    if (pr$tvar[1, 1] != x$pars$tvar[1, 1]) moved <- TRUE
  }
  expect_true(same)
  expect_true(moved)
  expect_identical(st$step[["b"]], .3); expect_identical(st$loc[["b"]], .1)
  expect_equal(unname(st$acc_in[2] + st$acc_loc[2]), 0)
  # nothing left in the sweep: it returns what it was given
  settings$active_scale[] <- FALSE; settings$active_loc[] <- FALSE
  sm <- EMC2:::scale_move_standard(sampler, x$pars, x$alpha, x$ll, settings, x$lik_prec)
  expect_identical(sm$pars, x$pars); expect_identical(sm$alpha, x$alpha); expect_identical(sm$settings, settings)
})

test_that("the sweep's first steps start at the moves' conditional scale when that is below the default", {
  sampler <- toy_emc[[1]]
  set.seed(4)
  x <- toy_state(sampler)
  run1 <- function(curv) {
    lik_prec <- lapply(1:2, function(s) list(prec = diag(curv), lin = rep(0, 2)))
    EMC2:::scale_move_standard(sampler, x$pars, x$alpha, x$ll, EMC2:::scale_move_init(NULL, toy_pars), lik_prec, gain = 0)$settings
  }
  a <- x$pars$a_half
  # a flat surrogate: the default starts (.5, .2), limited only by the prior part of the scale move
  flat <- run1(c(0, 0))
  A <- rep(sampler$prior$A, length.out = 2)
  expect_equal(unname(flat$step), pmin(.5, 2.4 / sqrt(4 / (A^2 * a))))
  expect_equal(unname(flat$loc), pmin(.2, 2.4 / sqrt(diag(sampler$prior$theta_mu_invar))))
  # a sharp one in parameter a: its steps start at 2.4 / sqrt(sum of curvatures), b's as before
  sharp <- run1(c(1e4, 0))
  R <- x$alpha - x$pars$tmu
  expect_equal(sharp$step[["a"]], unname(2.4 / sqrt(sum(R[1, ]^2) * 1e4 + 4 / (A[1]^2 * a[1]))))
  expect_equal(sharp$loc[["a"]], 2.4 / sqrt(2 * 1e4 + sampler$prior$theta_mu_invar[1, 1]))
  expect_equal(sharp$step[["b"]], flat$step[["b"]]); expect_equal(sharp$loc[["b"]], flat$loc[["b"]])
})

test_that("the block's exact check is timed serially and forked, and gives the same values either way", {
  skip_on_cran()
  skip_on_os("windows")
  sampler <- toy_emc[[1]]
  x <- toy_state(sampler)
  run <- function(n_cores, check_serial = NULL) {
    set.seed(9)
    st <- EMC2:::scale_move_init(NULL, toy_pars); st$iter <- 1; st$check_serial <- check_serial
    pr <- x$pars; al <- x$alpha; l <- x$ll
    for (i in 1:8) {
      sm <- EMC2:::scale_move_standard(sampler, pr, al, l, st, x$lik_prec, n_cores = n_cores)
      pr <- sm$pars; al <- sm$alpha; l <- sm$ll; st <- sm$settings
    }
    list(pars = pr, alpha = al, ll = l, settings = st)
  }
  one <- run(1); two <- run(2); ser <- run(2, TRUE); par <- run(2, FALSE)
  for (r in list(two, ser, par)) {
    expect_identical(r$pars, one$pars); expect_identical(r$alpha, one$alpha); expect_identical(r$ll, one$ll)
  }
  expect_null(one$settings$check_serial)                      # one core: nothing to choose
  expect_true(is.logical(two$settings$check_serial))          # timed, and decided
  expect_length(two$settings$check_time$parallel, 3)
})

# ---- the conjugate model: the exact posterior, and the funnel mixes ---------
# With T_j = 2 the third parameter is the funnel. conj_ref: log group SD under the exact posterior, from
# conj_exact(4e5, seed = 3, thin = 5) and seed = 7, averaged.
conj_pars <- c("strong", "mid", "weak")
conj_T <- c(200, 20, 2)
conj_n <- 8
conj_ybar <- withr::with_seed(11, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                              .rng_sample_kind = "Rejection", {
  alpha <- matrix(rnorm(conj_n * 3, rep(c(.2, -.2, 0), each = conj_n), .3), conj_n, 3)
  alpha + matrix(rnorm(conj_n * 3), conj_n, 3) / sqrt(rep(conj_T, each = conj_n))
})
colnames(conj_ybar) <- conj_pars
conj_emc <- conj_emc_from(conj_ybar, conj_T)
conj_ref <- list(mean = c(-1.077, -1.569, -1.820), sd = c(.276, .721, 1.148),
                 q50 = c(-1.10, -1.44, -1.60), q10 = c(-1.41, -2.355, -3.30))

# Gibbs sampler of the four full conditionals: the exact posterior of the group SDs
conj_exact <- function(iter, seed = 3, thin = 1) {
  set.seed(seed)
  prior <- conj_emc[[1]]$prior; v <- prior$v; A <- prior$A
  p <- 3; n <- conj_n; Tj <- conj_T; ybar <- conj_ybar
  V0inv <- solve(prior$theta_mu_var); m0 <- prior$theta_mu_mean
  mu <- rep(0, p); Sigma <- diag(.1, p); a <- rep(1, p); alpha <- ybar
  out <- matrix(NA_real_, iter %/% thin, p)
  for (i in seq_len(iter)) {
    Sinv <- solve(Sigma)
    Q <- V0inv + n * Sinv
    mu <- drop(mvtnorm::rmvnorm(1, solve(Q, V0inv %*% m0 + Sinv %*% colSums(alpha)), solve(Q)))
    r <- sweep(alpha, 2, mu)
    Sigma <- EMC2:::riwish(v + p - 1 + n, 2 * v * diag(1 / a) + crossprod(r))
    Sinv <- solve(Sigma)
    a <- 1 / rgamma(p, (v + p) / 2, rate = v * diag(Sinv) + 1 / A^2)
    Qa <- Sinv + diag(Tj); Va <- solve(Qa)
    for (s in 1:n) alpha[s, ] <- mvtnorm::rmvnorm(1, Va %*% (Sinv %*% mu + Tj * ybar[s, ]), Va)
    if (i %% thin == 0) out[i %/% thin, ] <- sqrt(diag(Sigma))
  }
  out[-(1:(nrow(out) / 10)), ]
}

test_that("with an exact (Gaussian) surrogate every block passes and the realised acceptance is the surrogate's", {
  sampler <- conj_emc[[1]]
  set.seed(6)
  alpha <- t(conj_ybar); n <- conj_n
  mu <- rowMeans(alpha); Sigma <- diag(.1, 3)
  pars <- list(tmu = mu, tvar = Sigma, tvinv = solve(Sigma), a_half = rep(1, 3), subj_mu = matrix(mu, 3, n), alpha = alpha)
  ll <- sapply(1:n, function(s) as.numeric(EMC2:::calc_ll_manager(
    matrix(alpha[, s], 1, dimnames = list(NULL, conj_pars)), sampler$data[[s]], sampler$model)))
  lik_prec <- lapply(sampler$data, function(d) list(prec = diag(d$T), lin = d$T * d$ybar))
  st <- EMC2:::scale_move_init(NULL, conj_pars)
  for (i in 1:100) {
    sm <- EMC2:::scale_move_standard(sampler, pars, alpha, ll, st, lik_prec)
    pars <- sm$pars; alpha <- sm$alpha; ll <- sm$ll; st <- sm$settings
  }
  expect_gt(st$acc_out / st$n_out, .98)
  expect_identical(st$acc_real, st$acc_in)
  expect_identical(st$acc_loc_real, st$acc_loc)
  # ... so the running block acceptance stays at 1 and the target is scale_move_target itself
  expect_identical(st$r_block, 1)
})

test_that("with the sweep fit() recovers the conjugate posterior and the funnel mixes", {
  skip_on_cran(); skip_if_not(slow_tests(), "EMC2_SLOW_TESTS is not true")
  skip_on_os("windows")
  withr::local_seed(123, .rng_kind = "L'Ecuyer-CMRG")
  conj_stop <- list(preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
                    sample = list(iter = 1500))
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 500)
  # the exact posterior of the group SDs (mean within 4 Monte-Carlo SEs, spread and median within what
  # 300-400 effective draws of a heavy-tailed log SD resolve)
  sd_draws <- lapply(emc, function(x) {
    i <- x$samples$stage == "sample"
    t(sqrt(apply(x$samples$theta_var[, , i, drop = FALSE], 3, diag)))
  })
  iat <- numeric(3)
  for (j in 1:3) {
    x <- log(do.call(rbind, sd_draws)[, j])
    ess_x <- sum(sapply(sd_draws, function(d) coda::effectiveSize(coda::mcmc(log(d[, j])))))
    z <- (mean(x) - conj_ref$mean[j]) / sqrt(var(x) / ess_x + (conj_ref$sd[j] / 60)^2)
    expect_lt(abs(z), 4)
    expect_lt(abs(log(sd(x) / conj_ref$sd[j])), .3)
    expect_lt(abs(median(x) - conj_ref$q50[j]), .25)
    expect_lt(abs(quantile(x, .1) - conj_ref$q10[j]), .4)
    iat[j] <- length(x) / ess_x
  }
  # the funnel parameter's log SD: a few iterations per effective draw
  expect_lt(iat[3], 15)
})
