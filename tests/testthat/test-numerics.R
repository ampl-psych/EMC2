test_that("Gaussian conditioning handles correlated, reordered coefficients", {
  L <- rbind(c(.2, 0, 0, 0), c(1, .3, 0, 0), c(-.4, .8, 1, 0), c(.5, -.2, .7, .4))
  covariance <- tcrossprod(L)
  mu <- c(-1, .5, 2, -.3)
  dependent <- c(3, 1)
  given <- c(4, 2)
  x <- c(.9, -2)
  precision <- solve(covariance)
  expected_var <- solve(precision[dependent, dependent])
  expected_mean <- mu[dependent] - expected_var %*% precision[dependent, given] %*% (x - mu[given])
  conditional <- condMVN(mu, covariance, dependent, given, x)
  expect_equal(conditional$condMean, c(expected_mean), tolerance = 1e-10)
  expect_equal(conditional$condVar, expected_var, tolerance = 1e-10)
})

test_that("covariance transformations include the full Jacobian", {
  L <- rbind(c(.2, 0, 0), c(.7, 1.4, 0), c(-.4, .8, .6))
  covariance <- tcrossprod(L)
  raw <- L
  diag(raw) <- log(diag(L))
  raw <- raw[lower.tri(raw, diag = TRUE)]
  expect_equal(unwind_chol(covariance), raw, tolerance = 1e-10)
  expect_equal(unwind_chol(raw, reverse = TRUE), covariance, tolerance = 1e-10)
  # Differentiate the unique covariance entries, including log diagonals.
  jacobian <- vapply(seq_along(raw), function(i) {
    step <- numeric(length(raw))
    step[i] <- 1e-5
    difference <- (unwind_chol(raw + step, TRUE) - unwind_chol(raw - step, TRUE)) / 2e-5
    difference[lower.tri(difference, diag = TRUE)]
  }, numeric(length(raw)))
  expect_equal(calc_log_jac_chol(raw), as.numeric(determinant(jacobian)$modulus), tolerance = 1e-7)
})

test_that("group designs preserve means and cross-parameter covariance", {
  withr::local_seed(123, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion")
  designs <- list(cbind(1, c(-2, .5, 3)), matrix(c(1, 2, -1), ncol = 1))
  beta <- c(.3, -.7, 1.2)
  expected <- rbind(c(1.7, -.05, -1.8), c(1.2, 2.4, -1.2))
  expect_equal(calculate_subject_means(designs, beta), expected)
  covariance <- matrix(c(.25, -.45, -.45, 1), 2)
  n <- 10000L
  draws <- draw_alpha_from_design(designs, matrix(beta, 3, n), array(covariance, c(2, 2, n)))
  residuals <- matrix(sweep(draws, 1:2, expected), nrow = 2)
  expect_lt(max(abs(rowMeans(residuals)) / sqrt(diag(covariance) / ncol(residuals))), 6)
  se <- sqrt((covariance^2 + outer(diag(covariance), diag(covariance))) / (ncol(residuals) - 1))
  expect_lt(max(abs(cov(t(residuals)) - covariance) / se), 6)
})

test_that("bridge integration recovers normal-model evidence", {
  y <- c(-.2, 1.1, 2.7)
  prior_mean <- -.5
  prior_sd <- 1.4
  noise <- .8
  variance <- 1 / (1 / prior_sd^2 + length(y) / noise^2)
  mean <- variance * (prior_mean / prior_sd^2 + sum(y) / noise^2)
  log_joint <- function(x) dnorm(x, prior_mean, prior_sd, log = TRUE) +
    vapply(x, function(mu) sum(dnorm(y, mu, noise, log = TRUE)), numeric(1))
  # Fixed quadrature points and a deliberately mismatched warp/proposal.
  z <- qnorm((seq_len(4096) - .5) / 4096)
  posterior <- mean + sqrt(variance) * z
  centre <- mean + .4
  scale <- 2 * sqrt(variance)
  log_pair <- function(x) {
    a <- log_joint(x)
    b <- log_joint(2 * centre - x)
    pmax(a, b) + log1p(exp(-abs(a - b)))
  }
  result <- run.iterative.scheme(
    q11 = log_pair(posterior), q12 = dnorm((posterior - centre) / scale, log = TRUE),
    q21 = log_pair(centre + scale * z), q22 = dnorm(z, log = TRUE),
    L = matrix(scale), r0 = 1, tol = 1e-10, criterion = "r", silent = TRUE,
    maxiter = 100, neff = NULL)
  expected <- mvtnorm::dmvnorm(y, rep(prior_mean, length(y)),
    diag(noise^2, length(y)) + prior_sd^2, log = TRUE)
  expect_equal(c(result$logml), expected, tolerance = 1e-4)
})

test_that("sample thinning selects the same iterations in every stored parameter", {
  means <- outer(1:4, 1:4, "+")
  samples <- list(samples = list(
    theta_mu = means, alpha = array(1:16, c(4, 1, 4)),
    theta_var = array(1:64, c(4, 4, 4)), stage = rep("sample", 4), idx = 4))
  thinned <- remove_samples(samples, thin = 2)$samples
  expect_equal(thinned$theta_mu, means[, c(1, 3)])
  expect_equal(thinned$alpha, samples$samples$alpha[, , c(1, 3), drop = FALSE])
  expect_equal(thinned$theta_var, samples$samples$theta_var[, , c(1, 3), drop = FALSE])
})
