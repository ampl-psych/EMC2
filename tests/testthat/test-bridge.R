# Closed-form log densities, independent of the bridge's density helpers.
bridge_normal <- function(x, mu, covariance) {
  z <- x - mu
  -.5 * (length(z) * log(2*pi) + log(det(covariance)) +
           drop(crossprod(z, solve(covariance, z))))
}

bridge_inverse_gamma <- function(x, shape, rate) {
  shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - rate / x
}

bridge_inverse_wishart <- function(covariance, df, scale) {
  p <- nrow(covariance)
  df/2 * log(det(scale)) - df*p/2 * log(2) -
    p*(p - 1)/4 * log(pi) - sum(lgamma((df + 1 - seq_len(p))/2)) -
    (df + p + 1)/2 * log(det(covariance)) -
    sum(diag(scale %*% solve(covariance)))/2
}

bridge_alpha <- list(
  rbind(c(-.3, .7, 1.2), c(.4, -.2, .8)),
  rbind(c(.9, -.4, .2), c(-.6, .3, 1.1)))
bridge_draws <- do.call(cbind, bridge_alpha)
bridge_mean <- c(.1, -.4, .6)
bridge_covariance <- tcrossprod(rbind(c(.8, 0, 0), c(.2, 1.1, 0), c(-.3, .4, .7)))

test_that("single bridge uses the full multivariate normal prior", {
  chain <- list(n_pars = 3, n_subjects = 2,
                prior = list(theta_mu_mean = bridge_mean, theta_mu_var = bridge_covariance))
  info <- make_info(chain, "single")
  packed <- bridge_add_group(bridge_draws, chain, 1:2, "single")
  expected <- vapply(1:2, function(i) sum(vapply(bridge_alpha, function(x)
    bridge_normal(x[i, ], bridge_mean, bridge_covariance), numeric(1))), numeric(1))
  expect_equal(packed, bridge_draws)
  expect_equal(h.unnormalized.posterior(packed, NULL, info, 1, hyper_only = TRUE),
               expected, tolerance = 1e-10)
})

for (variant in c("full", "diagonal", "blocked")) {
  test_that(paste(variant, "standard bridge matches covariance priors and Jacobians"), {
    blocks <- switch(variant, full = list(1:3), diagonal = list(1L, 2L, 3L),
                     blocked = list(1:2, 3L))
    groups <- rep(seq_along(blocks), lengths(blocks))
    prior <- list(theta_mu_mean = bridge_mean, theta_mu_var = bridge_covariance,
                  v = 4.2, A = c(1.2, .7, 1.9))
    chain <- list(n_pars = 3, n_subjects = 2, par_names = c("x", "y", "z"),
                  prior = prior, par_group = groups, is_blocked = lengths(blocks)[groups] > 1,
                  samples = list(theta_mu = cbind(c(.2, -.1, .4), c(-.3, .5, .1)),
                                 a_half = cbind(c(.8, 1.4, 2.1), c(1.3, .6, .9)),
                                 theta_var = array(0, c(3, 3, 2))))
    expected <- numeric(2)
    coordinates <- vector("list", 2)
    for (i in 1:2) {
      L <- (1 + .2 * (i - 1)) * rbind(c(.7, 0, 0), c(.3, 1.1, 0), c(-.2, .4, .9))
      L[outer(groups, groups, "!=")] <- 0
      covariance <- tcrossprod(L)
      chain$samples$theta_var[, , i] <- covariance
      mu <- chain$samples$theta_mu[, i]
      a <- chain$samples$a_half[, i]
      raw <- switch(variant,
        full = c(log(L[1,1]), L[2,1], L[3,1], log(L[2,2]), L[3,2], log(L[3,3])),
        diagonal = log(diag(covariance)),
        blocked = c(log(covariance[3,3]), log(L[1,1]), L[2,1], log(L[2,2])))
      jacobian <- switch(variant,
        full = log(8 * L[1,1]^4 * L[2,2]^3 * L[3,3]^2),
        diagonal = sum(log(diag(covariance))),
        blocked = log(4 * L[1,1]^3 * L[2,2]^2 * covariance[3,3]))
      coordinates[[i]] <- c(mu, log(a), raw)
      lp <- bridge_normal(mu, bridge_mean, bridge_covariance) +
        sum(bridge_inverse_gamma(a, .5, 1/prior$A^2)) + sum(log(a)) + jacobian
      for (block in blocks) {
        d <- length(block)
        lp <- lp + if (d == 1)
          bridge_inverse_gamma(covariance[block, block], prior$v/2, prior$v/a[block]) else
          bridge_inverse_wishart(covariance[block, block], prior$v + d - 1,
                                 diag(2 * prior$v/a[block]))
      }
      expected[i] <- lp + sum(vapply(bridge_alpha, function(x)
        bridge_normal(x[i, ], mu, covariance), numeric(1)))
    }
    idx <- c(2, 1)
    info <- make_info(chain, "standard")
    packed <- bridge_add_group(bridge_draws[idx, ], chain, idx, "standard")
    expect_equal(unname(packed),
                 unname(cbind(bridge_draws, do.call(rbind, coordinates))[idx, ]), tolerance = 1e-10)
    expect_equal(h.unnormalized.posterior(packed, NULL, info, 1, hyper_only = TRUE),
                 expected[idx], tolerance = 1e-10)
  })
}

for (correlated in c(FALSE, TRUE)) {
  test_that(paste("SEM bridge matches Gaussian marginalization with",
                  if (correlated) "correlated factors" else "independent factors"), {
    settings <- list(
      Lambda_mat = rbind(c(1, 0), c(0, 1), c(Inf, Inf)),
      B_mat = rbind(c(0, 0), c(Inf, 0)),
      K_mat = matrix(c(Inf, 0, 0)), G_mat = matrix(c(0, Inf)),
      factor_groups = if (correlated) c(1, 1) else 1:2,
      covariates = matrix(c(-.7, 1.1)))
    prior <- list(theta_mu_mean = bridge_mean, theta_mu_var = c(1.2, .8, 1.5),
                  lambda_mean = rbind(c(.1, .3), c(-.2, .4), c(.2, -.3)),
                  lambda_var = rbind(c(.7, 1.2), c(.8, .6), c(.4, .9)),
                  B_var = .8, K_var = 1.3, G_var = .6,
                  a_e = c(2, 3, 4), b_e = c(.7, 1.1, 1.3),
                  a_d = c(5, 5), b_d = c(1.2, .8))
    chain <- list(n_pars = 3, n_subjects = 2, n_factors = 2, prior = prior,
                  sem_settings = settings,
                  samples = list(theta_mu = cbind(c(.2, -.1, .4), c(-.3, .5, .1)),
                    lambda = array(0, c(3, 2, 2)), B = array(0, c(2, 2, 2)),
                    K = array(0, c(3, 1, 2)), G = array(0, c(2, 1, 2)),
                    epsilon_inv = cbind(c(1.4, .7, 2.2), c(.6, 1.8, 1.1)),
                    delta_inv = array(0, c(2, 2, 2))))
    expected <- numeric(2)
    coordinates <- vector("list", 2)
    for (i in 1:2) {
      mu <- chain$samples$theta_mu[, i]
      eps <- chain$samples$epsilon_inv[, i]
      l <- c(.5, -.4) + .2*(i - 1)
      b <- .3 - .1*(i - 1)
      k <- .15 + .1*(i - 1)
      g <- -.2 - .1*(i - 1)
      L <- (1 + .2*(i - 1)) * rbind(c(.9, 0), c(if (correlated) .25 else 0, 1.3))
      precision <- tcrossprod(L)
      chain$samples$lambda[, , i] <- rbind(c(1, 0), c(0, 1), l)
      chain$samples$B[2, 1, i] <- b
      chain$samples$K[1, 1, i] <- k
      chain$samples$G[2, 1, i] <- g
      chain$samples$delta_inv[, , i] <- precision
      raw <- if (correlated) c(log(L[1,1]), L[2,1], log(L[2,2])) else log(diag(precision))
      coordinates[[i]] <- c(mu, l, b, k, g, log(eps), raw)
      lp <- sum(dnorm(mu, bridge_mean, sqrt(prior$theta_mu_var), log = TRUE)) +
        sum(dnorm(l, c(.2, -.3), sqrt(c(.4, .9)), log = TRUE)) +
        dnorm(b, 0, sqrt(.8), log = TRUE) + dnorm(k, 0, sqrt(1.3), log = TRUE) +
        dnorm(g, 0, sqrt(.6), log = TRUE) +
        sum(dgamma(eps, prior$a_e, rate = prior$b_e, log = TRUE) + log(eps))
      if (correlated) {
        # Covariance ~ IW(df, S) implies precision ~ Wishart(df, inverse(S)).
        df <- prior$a_d[1]
        lp <- lp + (df - 3)/2 * log(det(precision)) -
          sum(prior$b_d * diag(precision))/2 + df/2 * sum(log(prior$b_d)) -
          df * log(2) - .5 * log(pi) - lgamma(df/2) - lgamma((df - 1)/2) +
          log(4 * L[1,1]^3 * L[2,2]^2)
      } else {
        lp <- lp + sum(dgamma(diag(precision), prior$a_d, rate = prior$b_d, log = TRUE) +
                         log(diag(precision)))
      }
      # z2 = b*z1 + g*x + noise; alpha3 = mu3 + l1*z1 + l2*z2 + noise.
      H <- rbind(c(1, 0), c(b, 1), c(l[1] + b*l[2], l[2]))
      covariance <- H %*% solve(precision) %*% t(H) + diag(1/eps)
      expected[i] <- lp + sum(vapply(1:2, function(s)
        bridge_normal(bridge_alpha[[s]][i, ], mu + c(k, g, l[2]*g)*settings$covariates[s, 1],
                      covariance), numeric(1)))
    }
    idx <- c(2, 1)
    info <- make_info(chain, "SEM")
    packed <- bridge_add_group(bridge_draws[idx, ], chain, idx, "SEM")
    expect_equal(unname(packed),
                 unname(cbind(bridge_draws, do.call(rbind, coordinates))[idx, ]), tolerance = 1e-10)
    expect_equal(h.unnormalized.posterior(packed, NULL, info, 1, hyper_only = TRUE),
                 expected[idx], tolerance = 1e-10)
  })
}
