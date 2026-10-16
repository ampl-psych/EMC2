library(EMC2)
library(posterior)

# Resume the same saved states under three priors, with a fixed short budget.
compare_bisset_priors <- function(input_file, output_dir, cores_per_chain = 8,
                                 labels = c("default", "centred", "regularised")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  old <- options(emc.adapt_converge = FALSE)
  on.exit(options(old))
  summaries <- list()
  for (label in labels) {
    emc <- EMC2:::restore_duplicates(readRDS(input_file))
    p <- emc[[1]]$prior
    if (label != "default") p$theta_mu_mean[c("S1.sigma", "S1.tau")] <- log(.05)
    if (label == "regularised") p$v <- 8
    for (ch in seq_along(emc)) emc[[ch]]$prior <- p
    emc <- EMC2:::map_pm_settings(emc, function(x) { x$n_particles <- 25; x })
    emc[[1]]$sample_kernel <- NULL
    set.seed(415, kind = "L'Ecuyer-CMRG")
    context <- EMC2:::sampler_contexts(emc, cores_per_chain)
    seconds <- system.time({
      emc <- run_emc(emc, "adapt", stop_criteria = list(iter = 100, min_unique = 0),
        step_size = 50, max_tries = 1, particle_factor = 25 / sqrt(emc[[1]]$n_pars),
        trim = FALSE, cores_for_chains = 3,
        cores_per_chain = cores_per_chain, .context = context)
      # run_emc adds 100 discarded iterations to tune the frozen sample kernel.
      saveRDS(emc, file.path(output_dir, paste0(label, "_warmup.rds")))
      emc <- run_emc(emc, "sample", stop_criteria = list(iter = 200),
        step_size = 100, particle_factor = 25 / sqrt(emc[[1]]$n_pars),
        trim = FALSE, cores_for_chains = 3,
        cores_per_chain = cores_per_chain, .context = context)
    })[["elapsed"]]
    diagnostics <- do.call(rbind, lapply(c("mu", "sigma2", "alpha"), function(group) {
      pars <- get_pars(emc, selection = group, stage = "sample", by_subject = TRUE)
      do.call(rbind, lapply(names(pars), function(name) {
        d <- as.data.frame(summarise_draws(as_draws_array(pars[[name]]),
          "mean", "sd", "rhat", "ess_bulk"))
        d$group <- group
        d$variable <- paste(name, d$variable, sep = ".")
        d
      }))
    }))
    summary <- data.frame(prior = label, seconds = seconds,
      rhat_mu = max(diagnostics$rhat[diagnostics$group == "mu"]),
      rhat_sigma2 = max(diagnostics$rhat[diagnostics$group == "sigma2"]),
      rhat_alpha = max(diagnostics$rhat[diagnostics$group == "alpha"]))
    print(summary)
    saveRDS(list(fit = emc, prior = p, diagnostics = diagnostics, summary = summary,
      session = sessionInfo()), file.path(output_dir, paste0(label, ".rds")))
    summaries[[label]] <- summary
  }
  summary <- do.call(rbind, summaries)
  write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
  summary
}
