# Load the installed EMC2 package, then source this file.
# fit is a saved emc object with a frozen sample kernel. Warm-up is excluded.
benchmark_sampler <- function(fit, draws = 1000, cores = 1, seed = 4101) {
  fit <- EMC2:::restore_sample_kernel(EMC2:::restore_duplicates(fit))
  fit <- subset(fit, stage = "sample", filter = EMC2::chain_n(fit)[1, "sample"] - 1)
  set.seed(seed)
  seconds <- system.time({
    out <- EMC2:::run_native_chains(fit, stage = "sample", iter = draws,
      n_cores = cores, cores_for_chains = length(fit))
  })[["elapsed"]]
  class(out) <- "emc"
  out <- subset(out, stage = "sample", filter = 1)
  diagnostics <- do.call(rbind, lapply(c("alpha", "mu", "sigma2"), function(group) {
    pars <- EMC2::get_pars(out, selection = group, stage = "sample", by_subject = TRUE)
    do.call(rbind, lapply(names(pars), function(name) {
      d <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(pars[[name]]),
        "mean", "sd", "rhat", "ess_bulk", "ess_tail"))
      d$variable <- paste(name, d$variable, sep = ".")
      d$group <- group
      d
    }))
  }))
  particles <- draws * sum(unlist(lapply(fit, function(chain)
    lapply(attr(chain$samples, "pm_settings"), function(subject)
      vapply(subject, function(block) block$n_particles + 1, 0)))))
  chain <- fit[[1]]
  moves <- if (chain$type != "single" && "sample" %in% EMC2:::scale_move_stages())
    chain$n_pars * (1 + (chain$type != "standard" ||
      (is.null(chain$group_designs) && length(chain$prior$theta_mu_mean) == chain$n_pars))) else 0
  joint <- draws * length(fit) * chain$n_subjects * moves
  list(fit = out, seconds = seconds, diagnostics = diagnostics,
       likelihoods = c(particle = particles, joint = joint, total = particles + joint))
}
