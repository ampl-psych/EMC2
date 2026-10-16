library(EMC2)
library(posterior)

# Positive regression effects become indistinguishable from zero in the
# negative log-scale tail. Both effects share an estimated population covariance.
tail_model <- function(seed = 61, subjects = 23, p = 11) {
  set.seed(seed)
  x <- diag(p)[rep(seq_len(p), each = 12), ]
  colnames(x) <- paste0("b", seq_len(p))
  des <- design_fmri(setNames(rep(list(x), subjects), seq_len(subjects)),
    model = MRI, pre_transform = list(func = c(b1 = "exp", b2 = "exp")))
  beta <- matrix(rnorm(subjects * p, sd = .3), subjects, p)
  beta[, 3] <- c(rep(0, 8), runif(subjects - 8, .5, 2))
  beta[, 1:2] <- -2 - 12 * beta[, 3]
  beta[1:8, 1:2] <- cbind(rnorm(8, -2, .3), rnorm(8, -2.5, .3))
  natural <- beta
  natural[, 1:2] <- exp(natural[, 1:2])
  dat <- data.frame(subjects = factor(rep(seq_len(subjects), each = nrow(x))),
    run = 1, time = rep(seq_len(nrow(x)), subjects),
    y = c(t(natural %*% t(x))) + rnorm(subjects * nrow(x), sd = .03))
  pr <- prior(des, mu_mean = c(setNames(rep(0, p), colnames(x)), sd = log(.03)), mu_sd = 1)
  emc <- make_emc(dat, des, prior_list = pr, n_chains = 4)
  for (ch in seq_along(emc)) {
    root <- diag(c(rep(.2, p), .1))
    root[3, 3] <- .7
    root[1:2, 3] <- -c(4, 8, 16, 32)[ch]
    covariance <- tcrossprod(root)
    emc[ch] <- init_chains(emc[ch], start_mu = c(rep(0, p), log(.03)),
      start_var = covariance, particles = 5000, cores_for_chains = 1)
  }
  list(emc = emc, seed = seed)
}

run_tails <- function(problem, label, warmup = 400, draws = 2000,
                      stages = c("preburn", "burn", "adapt", "sample")) {
  old <- options(emc.adapt_converge = FALSE)
  on.exit(options(old))
  set.seed(problem$seed + 100, kind = "L'Ecuyer-CMRG")
  fit <- problem$emc
  seconds <- system.time({
    for (stage in stages) {
      count <- switch(stage, preburn = 50, burn = 100, adapt = warmup, sample = draws)
      fit <- run_emc(fit, stage, stop_criteria = list(iter = count, min_unique =
        if (stage == "adapt") 0 else NULL), trim = FALSE, step_size = 100,
        max_tries = 1, cores_for_chains = 4, cores_per_chain = 1)
    }
  })[["elapsed"]]
  diagnostics <- tail_diagnostics(fit)
  print(aggregate(cbind(rhat, ess_bulk) ~ group, diagnostics,
    function(x) c(min = min(x), max = max(x))))
  cat(label, seconds, "seconds\n")
  list(fit = fit, diagnostics = diagnostics, seconds = seconds, label = label)
}

tail_diagnostics <- function(fit) {
  do.call(rbind, lapply(c("mu", "sigma2", "correlation", "alpha"), function(group) {
    pars <- get_pars(fit, selection = group, stage = "sample", by_subject = TRUE, flatten = TRUE)
    do.call(rbind, lapply(names(pars), function(name) {
      d <- as.data.frame(summarise_draws(as_draws_array(pars[[name]]),
        "mean", "sd", "rhat", "ess_bulk"))
      d$group <- group
      d$variable <- paste(name, d$variable, sep = ".")
      d
    }))
  }))
}
