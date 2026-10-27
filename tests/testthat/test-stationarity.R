test_that("a native particle update preserves a known normal posterior", {
  withr::local_seed(123, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion")
  # One lognormal observation and an N(2, 2) prior give an N(0, 1) posterior.
  # Each subject supplies an independent target draw for one particle update.
  n <- 4000L
  dat <- data.frame(subjects = factor(seq_len(n)), R = factor(1), rt = .2 + exp(-2))
  des <- design(data = dat, model = LNR, formula = list(m ~ 1),
                constants = c(s = log(sqrt(2)), t0 = log(.2)), report_p_vector = FALSE)
  pri <- prior(des, type = "single", theta_mu_mean = c(m = 2), theta_mu_var = matrix(2))
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
                  prior_list = pri, rt_resolution = 1e-12)
  x <- rnorm(n)
  emc[[1]]$samples$alpha[1, , 1] <- x
  emc[[1]]$samples$subj_ll[, 1] <- dlnorm(exp(-2), x, sqrt(2), log = TRUE)
  emc[[1]]$samples$idx <- 1L
  emc[[1]]$init <- TRUE
  emc[[1]]$chains_mu <- rep(list(c(m = 0)), n)
  emc[[1]]$chains_var <- rep(list(matrix(1)), n)
  settings <- list(list(mix = c(.97, .01, .02), epsilon = c(1, 1),
                        n_particles = 10, iter = 0))
  attr(emc[[1]]$samples, "pm_settings") <- rep(list(settings), n)
  updated <- EMC2:::run_native_chains(emc, "adapt", 1, particle_factor = 10, cores_for_chains = 1)
  draws <- updated[[1]]$samples$alpha[1, , 2]
  expect_true(any(draws != x))
  expect_lt(abs(mean(draws)), 6 / sqrt(n))
  expect_lt(abs(var(draws) - 1), 6 * sqrt(2 / (n - 1)))
  probability <- pnorm(2, lower.tail = FALSE)
  expect_lt(abs(mean(draws > 2) - probability), 6 * sqrt(probability * (1 - probability) / n))
})

test_that("hierarchical updates with affine and spread moves preserve a known joint posterior", {
  withr::local_seed(124, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion")
  withr::local_options(emc.scale_move = TRUE)
  n <- 2000L
  y <- rbind(c(-.8, 0, 1.3), c(-.5, .4, 1))
  dat <- data.frame(subjects = factor(rep(1:3, each = 2)),
                    C = factor(rep(1:2, 3)), R = factor(1), rt = .2 + exp(c(y)))
  des <- design(data = dat, model = LNR, formula = list(m ~ 0 + C),
                constants = c(s = 0, t0 = log(.2)), report_p_vector = FALSE)
  mean0 <- c(.2, -.3)
  variance0 <- diag(c(1, .7))
  df <- 8
  A <- c(.7, .9)
  pri <- prior(des, type = "standard", theta_mu_mean = mean0,
                theta_mu_var = variance0, v = df, A = A)
  template <- make_emc(dat, des, type = "standard", prior_list = pri,
                       n_chains = 1, rt_resolution = 1e-12)[[1]]
  template$init <- TRUE
  template$samples$idx <- 1L
  template$chains_mu <- rep(list(c(0, 0)), 3)
  template$chains_var <- rep(list(diag(2)), 3)
  # Fix the direction independently of the target draws below.
  template$lik_prec <- structure(vector("list", 3),
    affine_generators = list(rbind(c(.2, .4, .3), c(-.3, .1, -.2), c(0, 0, 0))))
  # Separate particle blocks must condition on earlier updates.
  attr(template$data, "components") <- 1:2

  # y_s | alpha_s ~ N(alpha_s, I). Integrating alpha and mu gives this
  # Gaussian marginal. Its zero-residual density at Sigma = 0 supplies the
  # rejection bound.
  base <- diag(6) + kronecker(matrix(1, 3, 3), variance0)
  bound <- -.5 * (6 * log(2*pi) + as.numeric(determinant(base, logarithm = TRUE)$modulus))
  precision0 <- solve(variance0)
  emc <- vector("list", n)
  for (i in seq_len(n)) {
    repeat {
      a <- 1 / rgamma(2, .5, rate = 1/A^2)
      Sigma <- solve(rWishart(1, df + 1, diag(a / (2*df)))[, , 1])
      marginal <- base + kronecker(diag(3), Sigma)
      log_accept <- mvtnorm::dmvnorm(c(y), rep(mean0, 3), marginal, log = TRUE) - bound
      if (log(runif(1)) < log_accept) break
    }
    precision <- solve(Sigma + diag(2))
    mu_var <- solve(precision0 + 3 * precision)
    mu <- c(mvtnorm::rmvnorm(1, mu_var %*% (precision0 %*% mean0 + precision %*% rowSums(y)),
                             mu_var))
    alpha_var <- solve(solve(Sigma) + diag(2))
    alpha_mean <- alpha_var %*% (y + c(solve(Sigma) %*% mu))
    alpha <- alpha_mean + t(mvtnorm::rmvnorm(3, sigma = alpha_var))
    chain <- template
    chain$samples$theta_mu[, 1] <- mu
    chain$samples$theta_var[, , 1] <- Sigma
    chain$samples$last_theta_var_inv <- solve(Sigma)
    chain$samples$a_half[, 1] <- a
    chain$samples$alpha[, , 1] <- alpha
    chain$samples$subj_ll[, 1] <- colSums(matrix(dnorm(c(y), c(alpha), 1, log = TRUE) - c(y), 2))
    emc[[i]] <- chain
  }
  updated <- EMC2:::run_native_chains(emc, "adapt", 1, particle_factor = 10, cores_for_chains = 1)
  moments <- function(chain, i) {
    s <- chain$samples
    mu <- s$theta_mu[, i]
    alpha <- s$alpha[, , i]
    Sigma <- s$theta_var[, , i]
    c(mu, mu^2, diag(Sigma), log(diag(Sigma)), Sigma[1, 2] / sqrt(prod(diag(Sigma))),
      log(s$a_half[, i]), alpha, alpha^2, mu * rowMeans(alpha), alpha[1, ] * alpha[2, ])
  }
  # Each column is an independent before/after pair, so its Monte Carlo SE
  # needs no convergence estimate or effective sample size.
  changes <- vapply(seq_len(n), function(i) moments(updated[[i]], 2) - moments(emc[[i]], 1),
                     numeric(length(moments(emc[[1]], 1))))
  expect_lt(max(abs(rowMeans(changes)) / (apply(changes, 1, sd) / sqrt(n))), 6)
  expect_true(any(updated[[1]]$samples$alpha[, , 2] != emc[[1]]$samples$alpha[, , 1]))
  expect_true(any(updated[[1]]$samples$theta_mu[, 2] != emc[[1]]$samples$theta_mu[, 1]))
  expect_true(any(updated[[1]]$samples$theta_var[, , 2] != emc[[1]]$samples$theta_var[, , 1]))
  expect_equal(attr(updated[[1]]$samples, "scale_move")$uses, 1)
  # One adaptation step increases .2 on acceptance and decreases it on rejection.
  for (move in c("affine_step", "spread_step")) {
    steps <- vapply(updated, function(x) attr(x$samples, "scale_move")[[move]], 0)
    expect_true(any(steps > .2) && any(steps < .2))
  }
})
