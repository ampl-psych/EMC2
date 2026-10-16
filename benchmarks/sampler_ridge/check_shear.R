library(EMC2)
library(posterior)

# Check the shear against full multivariate-normal and inverse-Wishart densities.
set.seed(71)
errors <- replicate(100, {
  p <- 5
  covariance <- crossprod(matrix(rnorm(p * p), p)) + diag(.1, p)
  precision <- solve(covariance)
  a <- exp(rnorm(p))
  mu <- rnorm(p)
  mu_covariance <- crossprod(matrix(rnorm(p * p), p)) + diag(1, p)
  mu_precision <- solve(mu_covariance)
  alpha <- matrix(rnorm(p * 7), p) + mu
  delta <- rnorm(1)
  pair <- sample(p, 2)
  j <- pair[1]; k <- pair[2]
  m <- diag(p); m[j, k] <- delta
  weights <- rexp(7)
  offset <- sum(weights * (mu[k] - alpha[k, ])) / (sum(weights) + mu_precision[j, j])
  new_mu <- mu
  new_mu[j] <- new_mu[j] + delta * offset
  new_alpha <- m %*% (alpha - mu) + new_mu
  proposed <- m %*% covariance %*% t(m)
  stopifnot(abs(precision[j, j] - solve(proposed)[j, j]) < 1e-8)
  normal <- function(x, v) sum(mvtnorm::dmvnorm(t(x), sigma = v, log = TRUE))
  iw <- function(v) -(2 + 2 * p) / 2 * determinant(v, logarithm = TRUE)$modulus -
    sum(diag(solve(v)) * 2 / a)
  direct <- normal(new_alpha - new_mu, proposed) - normal(alpha - mu, covariance) +
    iw(proposed) - iw(covariance) +
    mvtnorm::dmvnorm(new_mu, sigma = mu_covariance, log = TRUE) -
    mvtnorm::dmvnorm(mu, sigma = mu_covariance, log = TRUE)
  change <- delta * offset
  short <- 2 / a[k] * (2 * delta * precision[j, k] - delta^2 * precision[j, j]) -
    change * sum(mu_precision[j, ] * mu) - .5 * change^2 * mu_precision[j, j]
  reverse_offset <- sum(weights * (new_mu[k] - new_alpha[k, ])) /
    (sum(weights) + mu_precision[j, j])
  stopifnot(abs(offset - reverse_offset) < 1e-10)
  abs(direct - short)
})
stopifnot(max(errors) < 1e-8)
cat("Maximum full-density ratio error", max(errors), "\n")

# A constant compiled likelihood leaves the known Huang-Wand prior as target.
x <- matrix(0, 8, 2, dimnames = list(NULL, c("b1", "b2")))
des <- design_fmri(setNames(rep(list(x), 8), 1:8), model = MRI,
  transform = list(func = c(sd = "pnorm"), lower = c(sd = 1), upper = c(sd = 1)))
dat <- data.frame(subjects = factor(rep(1:8, each = 8)), run = 1,
  time = rep(1:8, 8), y = 0)
set.seed(72, kind = "L'Ecuyer-CMRG")
emc <- make_emc(dat, des, n_chains = 4)
options(emc.adapt_converge = FALSE)
for (stage in c("preburn", "burn", "adapt", "sample")) {
  emc <- run_emc(emc, stage, stop_criteria = list(iter = if (stage == "sample") 20000 else 200,
    min_unique = if (stage == "adapt") 0 else NULL), particle_factor = 20,
    step_size = 200, max_tries = 1, trim = FALSE, cores_for_chains = 4)
}
vars <- lapply(emc, function(x) {
  s <- x$samples
  idx <- which(s$stage == "sample")
  v <- s$theta_var[, , idx]
  cbind(mu = s$theta_mu[1, idx], log_sd = log(sqrt(v[1, 1, ])),
        correlation = v[1, 2, ] / sqrt(v[1, 1, ] * v[2, 2, ]))
})
draws <- as_draws_array(aperm(simplify2array(vars), c(1, 3, 2)))
print(summarise_draws(draws, "mean", "sd", "rhat", "ess_bulk"))
observed <- apply(do.call(rbind, vars), 2, quantile, c(.1, .5, .9))
expected <- cbind(mu = qnorm(c(.1, .5, .9)),
  log_sd = log(.3 * qt(c(.55, .75, .95), df = 2)),
  correlation = c(-.8, 0, .8))
print(observed)
print(expected)
stopifnot(max(abs(observed - expected)) < .08)
