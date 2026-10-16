library(EMC2)
library(posterior)

# Call once with the current installed sampler and the saved Bisset fit.
bisset_likelihood_reference <- function(input_file, output_dir) {
  emc <- EMC2:::restore_duplicates(readRDS(input_file))
  alpha <- lapply(seq_along(emc[[1]]$data), function(s) {
    do.call(rbind, lapply(emc, function(chain) {
      i <- tail(which(chain$samples$stage == "adapt"), 2)
      t(chain$samples$alpha[, s, i])
    }))
  })
  ll <- lapply(seq_along(alpha), function(s) {
    as.numeric(EMC2:::calc_ll_manager(alpha[[s]], emc[[1]]$data[[s]], emc[[1]]$model))
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(alpha = alpha, ll = ll, prior = emc[[1]]$prior),
    file.path(output_dir, "likelihood_reference.rds"))
}

# Run with the installed AccumulatR branch in a separate R library.
run_bisset_legacy <- function(analysis_file, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir)
  old <- setwd(dirname(analysis_file))
  on.exit(setwd(old))
  lines <- readLines(basename(analysis_file))
  setup <- new.env(parent = globalenv())
  set.seed(415, kind = "L'Ecuyer-CMRG")
  eval(parse(text = lines[seq_len(grep("^emc <- fit", lines)[1] - 1L)]), setup)
  emc <- setup$emc
  saveRDS(emc, file.path(output_dir, "input.rds"))

  # Compare the two bridges at identical saved participant coefficients.
  reference <- readRDS(file.path(output_dir, "likelihood_reference.rds"))
  ll <- lapply(seq_along(emc[[1]]$data), function(s) {
    as.numeric(EMC2:::calc_ll_manager(reference$alpha[[s]],
      emc[[1]]$data[[s]], emc[[1]]$model))
  })
  difference <- max(abs(unlist(ll) - unlist(reference$ll)))
  cat("Maximum bridge likelihood difference:", difference, "\n")
  stopifnot(difference < 1e-7,
    identical(emc[[1]]$prior$theta_mu_mean, reference$prior$theta_mu_mean),
    identical(emc[[1]]$prior$theta_mu_var, reference$prior$theta_mu_var),
    identical(emc[[1]]$prior$A, reference$prior$A),
    identical(emc[[1]]$prior$v, reference$prior$v))

  diagnostics <- function(emc, stage, last = 500) {
    draws <- lapply(emc, function(chain) {
      s <- chain$samples
      i <- tail(which(s$stage == stage), last)
      p <- dim(s$alpha)[1]
      n <- dim(s$alpha)[2]
      mu <- t(s$theta_mu[, i, drop = FALSE])
      variance <- t(vapply(i, function(j) diag(s$theta_var[, , j]), numeric(p)))
      alpha <- t(matrix(s$alpha[, , i, drop = FALSE], p * n))
      colnames(mu) <- paste0("mu.", rownames(s$theta_mu))
      colnames(variance) <- paste0("variance.", rownames(s$theta_mu))
      colnames(alpha) <- paste0("alpha.", rep(dimnames(s$alpha)[[2]], each = p),
        ".", rep(dimnames(s$alpha)[[1]], n))
      cbind(mu, variance, alpha)
    })
    a <- aperm(abind::abind(draws, along = 3), c(1, 3, 2))
    as.data.frame(summarise_draws(as_draws_array(a), "mean", "sd", "rhat", "ess_bulk"))
  }

  seconds <- c()
  for (stage in c("preburn", "burn", "adapt", "sample")) {
    iterations <- c(preburn = 50, burn = 2000, adapt = 800, sample = 1000)[[stage]]
    criteria <- list(iter = iterations)
    if (stage == "adapt") criteria$min_unique <- 0
    seconds[stage] <- system.time({
      emc <- run_emc(emc, stage, stop_criteria = criteria,
        step_size = if (stage == "preburn") 50 else 100,
        max_tries = 1, trim = FALSE, verbose = TRUE,
        cores_for_chains = 3, cores_per_chain = 23,
        fileName = file.path(output_dir, "checkpoint.RData"))
    })[["elapsed"]]
    d <- diagnostics(emc, stage)
    saveRDS(list(fit = emc, diagnostics = d, seconds = seconds,
      session = sessionInfo()), file.path(output_dir, paste0(stage, ".rds")))
    print(aggregate(rhat ~ group,
      transform(d, group = sub("\\..*", "", variable)), max))
  }
}
