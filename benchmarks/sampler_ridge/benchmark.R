library(EMC2)
library(posterior)

# A Gaussian regression with an unidentified direction, using the compiled MRI
# likelihood and the ordinary full-covariance hierarchical particle sampler.
ridge_model <- function(seed = 41, strength = 100, p = 4, subjects = 20,
                        observations = 24, weak = 1, weak_strength = 0) {
  set.seed(seed, kind = "Mersenne-Twister")
  rotation <- qr.Q(qr(matrix(rnorm(p * p), p)))
  strong <- seq_len(p - weak)
  direction <- rotation[, -strong, drop = FALSE]
  measurements <- qr.Q(qr(matrix(rnorm(observations * length(strong)), observations)))
  x <- strength * measurements %*% t(rotation[, strong, drop = FALSE])
  if (weak_strength > 0) {
    extra <- qr.Q(qr(cbind(measurements,
      matrix(rnorm(observations * weak), observations))))[, -strong, drop = FALSE]
    x <- x + weak_strength * extra %*% t(direction)
  }
  colnames(x) <- paste0("b", seq_len(p))
  dm <- setNames(rep(list(x), subjects), seq_len(subjects))
  des <- design_fmri(dm, model = MRI)
  dat <- data.frame(subjects = factor(rep(seq_len(subjects), each = nrow(x))),
                    run = 1, time = rep(seq_len(nrow(x)), subjects),
                    y = rnorm(subjects * nrow(x), sd = .3))
  pr <- prior(des, mu_mean = c(setNames(rep(0, p), colnames(x)), sd = log(.3)),
              mu_sd = 1)
  emc <- make_emc(dat, des, prior_list = pr, n_chains = 4)
  # Start each chain on the ridge, but in a different place. The same saved
  # initial states are used for both samplers.
  for (ch in seq_along(emc)) {
    start <- c(-2, -1, 1, 2)[(ch + seq_len(weak) - 2) %% 4 + 1]
    mu <- c(direction %*% start, log(.3))
    emc[ch] <- init_chains(emc[ch], start_mu = mu,
      start_var = diag(c(rep(.0001, p), .01)), particles = 1000,
      cores_for_chains = 1)
  }
  list(emc = emc, direction = direction, seed = seed, strength = strength,
       weak_strength = weak_strength)
}

run_ridge <- function(problem, label, draws = 1000, warmup = 600) {
  old_options <- options(emc.adapt_converge = FALSE)
  on.exit(options(old_options))
  set.seed(problem$seed + 100, kind = "L'Ecuyer-CMRG")
  seconds <- system.time({
    fit <- problem$emc
    for (stage in c("preburn", "burn", "adapt", "sample")) {
      count <- switch(stage, preburn = 50, burn = 100, adapt = warmup, sample = draws)
      fit <- run_emc(fit, stage, stop_criteria = list(iter = count, min_unique =
        if (stage == "adapt") 0 else NULL), trim = FALSE, step_size = 100,
        cores_for_chains = 4, cores_per_chain = 1)
    }
  })[["elapsed"]]
  mu <- get_pars(fit, selection = "mu", stage = "sample", return_mcmc = FALSE)
  # Project only for diagnostics. The sampler never receives this direction.
  p <- nrow(problem$direction)
  ridge <- aperm(simplify2array(lapply(mu, function(x)
    t(t(problem$direction) %*% x[seq_len(p), ]))), c(1, 3, 2))
  d <- summarise_draws(as_draws_array(ridge),
                       "mean", "sd", "rhat", "ess_bulk")
  diagnostics <- do.call(rbind, lapply(c("mu", "sigma2", "alpha"), function(group) {
    pars <- get_pars(fit, selection = group, stage = "sample", by_subject = TRUE)
    do.call(rbind, lapply(names(pars), function(name) {
      out <- as.data.frame(summarise_draws(as_draws_array(pars[[name]]), "rhat", "ess_bulk"))
      out$group <- group
      out$variable <- paste(name, out$variable, sep = ".")
      out
    }))
  }))
  result <- list(fit = fit, ridge = ridge, seconds = seconds,
    diagnostics = diagnostics,
    summary = data.frame(sampler = label, seed = problem$seed,
      direction = seq_len(nrow(d)), seconds = seconds, mean = d$mean, sd = d$sd, rhat = d$rhat,
      ess = d$ess_bulk, ess_second = d$ess_bulk / seconds))
  print(result$summary)
  result
}
