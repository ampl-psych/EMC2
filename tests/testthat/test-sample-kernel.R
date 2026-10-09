# Start at the known target: a single particle sweep must preserve it.
particle_test_kernel <- function(mu, Sigma, log_lik, components = rep(1, length(mu)),
                                 shared = components, mix = c(1, 0, 0), chains_mu = mu) {
  model <- function() list(log_likelihood = function(pars, ...) log_lik(pars))
  settings <- rep(list(list(mix = mix, epsilon = c(1, 1), n_particles = 25)),
                  length(unique(components)))
  tune <- list(components = components, shared_ll_idx = shared, exact = TRUE, frozen = TRUE)
  function(x) EMC2:::new_particle(
    1, data.frame(subjects = 1), settings, chains_mu = chains_mu, chains_var = Sigma,
    prev_ll = log_lik(x), parameters = list(
      alpha = matrix(x, ncol = 1, dimnames = list(names(mu), NULL)),
      subj_mu = matrix(mu, ncol = 1), tvar = Sigma),
    model = model, stage = "adapt", type = "standard", tune = tune)
}

test_that("independent particles preserve an imbalanced mixture", {
  set.seed(701)
  log_lik <- function(x) log(.98 * dnorm(x, -6) + .02 * dnorm(x, 6)) - dnorm(x, -6, log = TRUE)
  step <- particle_test_kernel(c(x = -6), matrix(1), log_lik,
                               mix = c(.98, 0, .02), chains_mu = c(x = 6))
  initial <- rnorm(6000, ifelse(runif(6000) < .02, 6, -6))
  draws <- vapply(initial, function(x) step(x)$proposal, numeric(1))
  expect_lt(abs(mean(draws > 0) - .02), .009)
})

test_that("sequential particle blocks preserve a correlated Gaussian target", {
  set.seed(702)
  Sigma <- matrix(c(1, .8, .8, 1), 2)
  initial <- mvtnorm::rmvnorm(3000, sigma = Sigma)
  for (shared_likelihood in c(FALSE, TRUE)) {
    # The correlation comes either from the prior or from one shared likelihood.
    log_lik <- if (shared_likelihood) function(x) {
      -.5 * (sum(x^2) - 1.6 * prod(x)) / .36 + .5 * sum(x^2)
    } else function(x) 0
    step <- particle_test_kernel(c(x = 0, y = 0),
      if (shared_likelihood) diag(2) else Sigma, log_lik,
      components = 1:2, shared = if (shared_likelihood) c(1, 1) else 1:2)
    out <- lapply(seq_len(nrow(initial)), function(i) step(initial[i, ]))
    draws <- t(vapply(out, `[[`, numeric(2), "proposal"))
    expect_lt(max(abs(cov(draws) - Sigma)), .1)
    expect_lt(max(abs(colMeans(draws))), .06)
    expect_equal(vapply(out, `[[`, numeric(1), "ll"), apply(draws, 1, log_lik),
                 tolerance = 1e-12)
  }
})

test_that("the sample-stage proposals and tuning stay fixed", {
  skip_on_cran()
  skip_on_os("windows")
  set.seed(11, kind = "Mersenne-Twister")
  pars <- c("a", "b", "c", "d")
  T <- c(200, 20, 5, 2)
  dat <- expand.grid(par = pars, subjects = factor(1:4))
  dat$T <- T
  dat$ybar <- unlist(lapply(1:4, function(s) rnorm(4, 0, .3) + rnorm(4) / sqrt(T)))
  ll <- function(pars, dadm, ...) -.5 * sum(dadm$T * (dadm$ybar - pars)^2)
  des <- design(model = ll, custom_p_vector = pars, report_p_vector = FALSE)
  emc <- make_emc(dat, des, type = "standard", n_chains = 2, compress = FALSE)
  criteria <- list(preburn = list(iter = 10), burn = list(mean_gd = 2.5),
                   adapt = list(min_unique = 20), sample = list(iter = 20))
  set.seed(123, kind = "L'Ecuyer-CMRG")
  emc <- fit(emc, cores_for_chains = 1, stop_criteria = criteria, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  set.seed(124)
  continued <- fit(emc, cores_for_chains = 1, verbose = FALSE, particle_factor = 20,
                   step_size = 20, stop_criteria = list(sample = list(iter = 40)))
  expect_length(emc[[1]]$sample_kernel$chains, length(emc))
  expect_identical(continued[[1]]$sample_kernel, emc[[1]]$sample_kernel)
  for (ch in seq_along(emc)) {
    before <- attr(emc[[ch]]$samples, "pm_settings")
    after <- attr(continued[[ch]]$samples, "pm_settings")
    for (s in seq_along(before)) {
      fields <- c("epsilon", "mix", "n_particles")
      expect_identical(after[[s]][[1]][fields], before[[s]][[1]][fields])
    }
    fields <- c("step", "loc")
    expect_identical(attr(continued[[ch]]$samples, "scale_move")[fields],
                     attr(emc[[ch]]$samples, "scale_move")[fields])
  }
})

test_that("likelihood bounds do not create spurious proposal precision", {
  model <- function() list(log_likelihood = function(pars, ...) {
    if (pars[2] >= 1) -1e4 else -200 * (pars[1] - .3)^2
  })
  lik <- EMC2:::lik_precision(c(a = .3, b = 0), c(.1, .1), data.frame(subjects = 1), model)
  expect_equal(unname(lik$prec), diag(c(400, 0)), tolerance = 1e-3)
  expect_equal(unname(lik$lin), c(120, 0), tolerance = 1e-3)
})
