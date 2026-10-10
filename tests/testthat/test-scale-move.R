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

toy_pars <- c("a", "b")

toy_subject_ll <- function(sampler, s, x)
  as.numeric(EMC2:::calc_ll_manager(matrix(x, 1, dimnames = list(NULL, toy_pars)), sampler$data[[s]], sampler$model))

toy_state <- function(sampler, alpha = matrix(c(.3, -.2, -.3, .25), 2, 2)) {
  p <- 2; n <- 2; mu <- c(.1, -.1); a <- c(.5, 2); Sigma <- diag(c(.3, .3))
  dimnames(alpha) <- list(toy_pars, NULL)
  list(pars = list(tmu = mu, tvar = Sigma, tvinv = solve(Sigma), a_half = a, subj_mu = matrix(mu, p, n), alpha = alpha),
       alpha = alpha, ll = sapply(1:n, function(s) toy_subject_ll(sampler, s, alpha[, s])))
}

# Independent Gibbs/Metropolis reference for the posterior-preservation check.
toy_run <- function(sampler, iter, use_move, seed) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  x <- toy_state(sampler, alpha = matrix(0, 2, 2))
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
      sm <- EMC2:::scale_move_standard(sampler, pars, alpha, cur_ll, settings, frozen = i > iter / 4)
      Sigma <- sm$pars$tvar; a <- sm$pars$a_half; mu <- sm$pars$tmu; alpha <- sm$alpha; cur_ll <- sm$ll
      settings <- sm$settings
    }
    out[i, ] <- c(Sigma[1, 1], Sigma[2, 2], Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2]), alpha[1, 1], alpha[2, 2], a, mu)
  }
  out[-seq_len(iter / 4), ]
}

test_that("scale and location moves preserve the posterior", {
  skip_on_cran()
  skip_on_ci()
  set.seed(21, kind = "Mersenne-Twister")
  dat <- do.call(rbind, lapply(1:2, function(s) data.frame(
    subjects = s, par = toy_pars, T = c(30, 3),
    ybar = c(.4, -.3) * (s - 1.5) * 2 + rnorm(2, 0, .2))))
  dat$subjects <- factor(dat$subjects)
  ll <- function(pars, dadm, ...) -sum(dadm$T * log1p((dadm$ybar - pars)^2))
  des <- design(model = ll, custom_p_vector = toy_pars, report_p_vector = FALSE)
  sampler <- make_emc(dat, des, type = "standard", n_chains = 2, compress = FALSE)[[1]]
  ref <- toy_run(sampler, 12000, use_move = FALSE, seed = 1)
  for (seed in c(2, 4)) {
    draws <- toy_run(sampler, 12000, use_move = TRUE, seed = seed)
    for (k in colnames(ref)) {
      lg <- k %in% c("s11", "s22", "aux1", "aux2")
      x <- if (lg) log(draws[, k]) else draws[, k]
      y <- if (lg) log(ref[, k]) else ref[, k]
      se <- sqrt(var(x) / coda::effectiveSize(x) + var(y) / coda::effectiveSize(y))
      expect_lt(abs((mean(x) - mean(y)) / se), 4)
      expect_lt(abs(log(sd(x) / sd(y))), .25)
    }
  }
})
