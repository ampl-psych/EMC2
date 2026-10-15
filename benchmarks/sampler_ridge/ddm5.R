library(EMC2)
library(posterior)

prepare_ddm5 <- function() {
  setup <- readRDS("paper/script/Comparison_Stan/ddm5_delta90_4x2/setup.rds")
  set.seed(124, kind = "L'Ecuyer-CMRG")
  seconds <- system.time({
    emc <- make_emc(setup$data, setup$design, n_chains = 4, compress = FALSE,
                    type = "diagonal", rt_resolution = 1e-12)
    emc <- init_chains(emc, cores_for_chains = 4, cores_per_chain = 2)
  })[["elapsed"]]
  list(emc = emc, setup_seconds = seconds)
}

run_ddm5 <- function(input, label) {
  old_options <- options(emc.adapt_converge = FALSE)
  on.exit(options(old_options))
  set.seed(125, kind = "L'Ecuyer-CMRG")
  emc <- input$emc
  seconds <- c(setup = input$setup_seconds)
  seconds["context"] <- system.time({
    context <- EMC2:::sampler_contexts(emc, 2)
  })[["elapsed"]]
  iterations <- c(preburn = 50, burn = 150, adapt = 300, sample = 1000)
  for (stage in names(iterations)) {
    seconds[stage] <- system.time({
      emc <- run_emc(emc, stage, stop_criteria = list(iter = iterations[[stage]],
        min_unique = if (stage == "adapt") 0 else NULL), step_size = 50,
        max_tries = 1, trim = FALSE, cores_for_chains = 4, cores_per_chain = 2,
        .context = context)
    })[["elapsed"]]
    cat(label, stage, seconds[stage], "seconds\n")
  }
  diagnostics <- do.call(rbind, lapply(c("mu", "sigma2", "alpha"), function(group) {
    pars <- get_pars(emc, selection = group, stage = "sample", by_subject = TRUE)
    do.call(rbind, lapply(names(pars), function(name) {
      out <- as.data.frame(summarise_draws(as_draws_array(pars[[name]]), "rhat", "ess_bulk"))
      out$group <- group
      out$variable <- paste(name, out$variable, sep = ".")
      out
    }))
  }))
  metrics <- aggregate(ess_bulk ~ group, diagnostics, min)
  metrics$ess_second <- metrics$ess_bulk / sum(seconds)
  print(seconds)
  print(metrics)
  cat("Maximum Rhat", max(diagnostics$rhat), "\n")
  list(fit = emc, seconds = seconds, diagnostics = diagnostics, metrics = metrics,
       label = label, library = find.package("EMC2"), session = sessionInfo())
}
