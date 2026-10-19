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
  expect_lt(abs(mean(draws)), 6 / sqrt(n))
  expect_lt(abs(var(draws) - 1), 6 * sqrt(2 / (n - 1)))
  probability <- pnorm(2, lower.tail = FALSE)
  expect_lt(abs(mean(draws > 2) - probability), 6 * sqrt(probability * (1 - probability) / n))
})
