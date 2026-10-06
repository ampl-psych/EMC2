get_stop_criteria <- function(stage, stop_criteria, type){
  if(is.null(stop_criteria)){
    if(stage == "preburn"){
      stop_criteria$iter <- 50
    }
    if(stage == "burn"){
      stop_criteria$mean_gd <- 1.1
      stop_criteria$omit_mpsrf <- TRUE
      if(type != "single"){
        stop_criteria$selection <- c("alpha", "mu")
      } else{
        stop_criteria$selection <- c("alpha")
      }

    }
    if(stage == "adapt"){
      stop_criteria$min_unique <- 150
    }
    if(stage == "sample"){
      stop_criteria$max_gd <- 1.1
      stop_criteria$omit_mpsrf <- TRUE
      if(type != "single"){
        stop_criteria$selection <- c("alpha", "mu")
      } else{
        stop_criteria$selection <- c("alpha")
      }
    }
  }
  if(!is.null(stop_criteria$max_gd) || !is.null(stop_criteria$mean_gd)){
    if(is.null(stop_criteria$selection)) stop_criteria$selection <- c('alpha', 'mu')
    if(is.null(stop_criteria$omit_mpsrf)) stop_criteria$omit_mpsrf <- TRUE
  }
  if(!is.null(stop_criteria$gd_quantile)){
    q <- stop_criteria$gd_quantile
    if(!is.numeric(q) || length(q) != 1 || is.na(q) || q <= 0 || q > 1) stop("gd_quantile must be a single number in (0, 1]")
    if(is.null(stop_criteria$max_gd)) stop("gd_quantile only applies to max_gd")
  }
  # min_es also needs a selection: without one, check_progress() has nothing to
  # take the effective size of and the criterion is silently always satisfied.
  if(!is.null(stop_criteria$min_es) && is.null(stop_criteria$selection)){
    stop_criteria$selection <- if(type != "single") c('alpha', 'mu') else 'alpha'
  }
  if(stage == "adapt" & is.null(stop_criteria$min_unique)) stop_criteria$min_unique <- 600
  if(stage != "adapt" & !is.null(stop_criteria$min_unique)) stop("min_unique only applicable for adapt stage, try min_es instead.")
  return(stop_criteria)
}

#' Fine-Tuned Model Estimation
#'
#' Although typically users will rely on ``fit``, this function can be used for more fine-tuned specification of estimation needs.
#' The function will throw an error if a stage is skipped,
#' the stages have to be run in order ("preburn", "burn", "adapt", "sample").
#' More details can be found in the ``fit`` help files (``?fit``).
#'
#' @param emc An emc object
#' @param stage A string. Indicates which stage is to be run, either `preburn`, `burn`, `adapt` or `sample`
#' @param search_width A double. Tunes target acceptance probability of the MCMC process.
#' This fine-tunes the width of the search space to obtain the desired acceptance probability.
#' 1 is the default width, increases lead to broader search.
#' @param step_size An integer. After each step, the stopping requirements as
#' specified by `stop_criteria` are checked and, in the stages before `sample`, proposal distributions are updated. Defaults to 100.
#' @param verbose Logical. Whether to print messages between each step with the current status regarding the stop_criteria.
#' @param verboseProgress Logical. Whether to print a progress bar within each step or not. Will print one progress bar for each chain and only if cores_for_chains = 1.
#' @param fileName A string. If specified will autosave emc at this location on every iteration.
#' @param particle_factor An integer. `particle_factor` multiplied by the square root of the number of sampled parameters determines the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain.
#' Parallelizes across participant calculations. Only available on Linux or Mac OS.
#' For Windows, only parallelization across chains (``cores_for_chains``) is available.
#' @param cores_for_chains An integer. How many cores to use across chains.
#' Defaults to the number of chains. the total number of cores used is equal to ``cores_per_chain`` * ``cores_for_chains``.
#' @param max_tries An integer. How many times should it try to meet the finish
#' conditions as specified by stop_criteria? Defaults to 20. max_tries is ignored if the required number of iterations has not been reached yet.
#' @param n_blocks An integer. Number of blocks. Will block the parameter chains such that they are updated in blocks. This can be helpful in extremely tough models with a large number of parameters.
#' @param stop_criteria A list. Defines the stopping criteria and for which types of parameters these should hold. See ``?fit``.
#' @param thin A boolean. If `TRUE` will automatically thin the MCMC samples, closely matched to the ESS.
#' Can also be set to a double, in which case 1/thin of the chain will be removed (does not have to be an integer).
#' @param trim A boolean. If `TRUE` will automatically remove redundant samples (i.e. from preburn, burn, adapt).
#' @param on_singular A list or `NULL` (default). Controls recovery when the group-level
#' covariance becomes computationally singular during sampling. `NULL` errors immediately.
#' A list may set `max_retries` (integer, re-draw the group step on failure), `on_exhausted`
#' (`"error"` or `"carry_forward"` the previous group parameters), `max_carry_forward`
#' (consecutive carried-forward iterations before giving up).
#' @param r_cores An integer for number of cores to use in R-based likelihood calculations, default 1.
#' @export
#' @return An emc object
#' @examples \donttest{
#' # First define a design
#' design_in <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1),
#'                            constants=c(s=log(1)))
#' # Then make the emc, we've omitted a prior here for brevity so default priors will be used.
#' emc <- make_emc(forstmann, design_in, compress = FALSE)
#'
#' # Now for example we can specify that we only want to run the "preburn" phase
#' # for MCMC 10 iterations
#' # emc <- run_emc(emc, stage = "preburn", stop_criteria = list(iter = 10), cores_for_chains = 1)
#'}

run_emc <- function(emc, stage, stop_criteria,
                    search_width = 1, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                    fileName = NULL,particle_factor=50, cores_per_chain = 1,
                    cores_for_chains = length(emc), max_tries = 20, n_blocks = 1,
                    thin = FALSE, trim = TRUE, on_singular = NULL, r_cores=1){
  if(length(emc) == 1) stop("run_emc() currently requires n_chains > 1.")
  emc <- restore_duplicates(emc)
  if(Sys.info()[1] == "Windows" & cores_per_chain > 1) stop("only cores_for_chains can be set on Windows")
  if (verbose) message(paste0("Running ", stage, " stage"))
  total_iters_stage <- chain_n(emc)[,stage][1]
  if(stage != "preburn"){
    iter <- stop_criteria[["iter"]] + total_iters_stage
  } else{
    iter <- stop_criteria[["iter"]]
  }
  progress <- check_progress(emc, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains, verbose, n_blocks = n_blocks)
  emc <- progress$emc
  progress <- progress[!names(progress) == 'emc']
  # We need to multiply step_size by thin to make an accurate guess for good step_size.
  cur_thin <- ifelse(is.numeric(thin), thin, 1)
  while(!progress$done){
    # More adapt draws: a sample-stage kernel built earlier is out of date
    if(stage == "adapt") emc[[1]]$sample_kernel <- NULL
    emc <- reset_pm_settings(emc, stage)
    # Remove redundant samples
    if(trim){
      emc <- fit_remove_samples(emc)
    }
    if(!is.null(progress$n_blocks)) n_blocks <- progress$n_blocks
    emc <- add_proposals(emc, stage, cores_per_chain*cores_for_chains, n_blocks)
    last_stage <- get_last_stage(emc)
    if(stage == "preburn"){
      sub_emc <- emc
    } else if(stage != last_stage){
      sub_emc <- subset(emc, filter = chain_n(emc)[1,last_stage]-1, stage = last_stage)
    } else{
      sub_emc <- subset(emc, filter = chain_n(emc)[1,stage] - 1, stage = stage)
    }
    t0 <- Sys.time()
    # Actual sampling
    sub_emc <- auto_mclapply(sub_emc,run_stages, stage = stage, iter= progress$step_size*max(1,cur_thin),
                             verbose=verbose,  verboseProgress = verboseProgress,
                             particle_factor=particle_factor,search_width=search_width,
                             n_cores=cores_per_chain, mc.cores = cores_for_chains,
                             on_singular = on_singular, r_cores = r_cores)

    # A chain that errors inside the parallel sampler comes back as a try-error
    # (or a value with no $samples) rather than a sampler object. Report which
    # chain failed and why here, instead of failing later with a cryptic error
    # in chain_n() ("$ operator is invalid for atomic vectors").
    check_chain_failures(sub_emc, stage, fileName)

    class(sub_emc) <- "emc"
    if(cores_for_chains > 1) sub_emc <- fix_custom_kernel_pointers(sub_emc, emc)
    if(stage != 'preburn'){
      if(is.numeric(thin)){
        sub_emc <- subset(sub_emc, stage = c("preburn", "burn", "adapt", "sample"), thin = thin)
      } else if(thin){
        sub_emc <- auto_thin(sub_emc, stage = c("preburn", "burn", "adapt", "sample"))
        # Update current rough guess for thinning:
        cur_thin <- progress$step_size/chain_n(sub_emc)[1,stage]
      }
      emc <- concat_emc(emc, sub_emc, progress$step_size, stage)
    } else{
      emc <- sub_emc
    }
    progress <- check_progress(emc, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains,
                               verbose, progress,n_blocks)
    emc <- progress$emc
    progress <- progress[!names(progress) == 'emc']
    # the interweaving sweep's gate, once its steps have had adapt_converge$min sweeps to settle
    if(stage == "adapt" && !legacy_sampler()){
      emc <- gate_scale_move(emc, kernel_window, decide = chain_n(emc)[1, "adapt"] >= adapt_converge$min, verbose = verbose)
    }
    if(!is.null(fileName)){
      emc <- strip_duplicates(emc)
      fileName <- fix_fileName(fileName)
      class(emc) <- "emc"
      save(emc, file = fileName)
      emc <- restore_duplicates(emc)
    }

    elapsed <- Sys.time() - t0
    if (verbose) {
      # get Gelman's diagnostic
      gd <- progress$gd$gd
      gd_message <- NULL
      if(!is.null(stop_criteria$mean_gd)) gd_message <- paste0("Mean Rhat=", round(mean(gd), 3))
      if(!is.null(stop_criteria$max_gd)) gd_message <- paste0("Max Rhat=", round(max(gd), 3))
      # get min effective sample size (ESS)
      ess_message <- NULL
      if(stage == 'sample' & !is.null(progress$curr_min_es)) {
        ess_message <- paste0(" | min ESS=", round(progress$curr_min_es))
      }
      # Get current iteration count for the sample stage
      current_iters <- if (stage == "sample") chain_n(emc)[1, stage] else NULL
      target_iters  <- if (stage == "sample") stop_criteria[["iter"]] else NULL

      # this one is still a little buggy, should fix
      rem <- estimate_remaining_total_time(stage= stage,tries_done = progress$trys, elapsed_dt = elapsed, max_tries = max_tries,
                                           current_iters = current_iters, target_iters = target_iters, step_size = progress$step_size)

      message(sprintf("[%s | try=%d | iters=%d%s%s] Duration: %s - ETA: %s-%s",
                      stage, progress$trys, progress$total_iters,
                      ifelse(!is.null(gd_message), paste0(" | ", gd_message), ""),
                      ifelse(!is.null(ess_message), ess_message, ""),
                      format_duration(elapsed),
                      format_duration(rem$min_time),
                      format_duration(rem$max_time)))
    }
  }

  if(stage == "adapt" && !legacy_sampler() && is.null(emc[[1]]$sample_kernel) &&
     get_last_stage(emc) == "adapt"){
    emc <- tune_sample_kernel(emc, verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName, particle_factor = particle_factor, search_width = search_width,
                              cores_per_chain = cores_per_chain, cores_for_chains = cores_for_chains,
                              n_blocks = n_blocks, on_singular = on_singular, r_cores = r_cores)
  }

  emc <- strip_duplicates(emc)
  class(emc) <- "emc"
  return(emc)
}

# The sample stage runs one fixed kernel: a kernel re-chosen from the chain's
# own recent draws is not a valid MCMC kernel. It is built here, at the end of
# adapt, from the last kernel_window adapt iterations (add_proposals), and
# tuned in one more 100-iteration adapt step (labelled adapt, not kept).
# Because it is never revisited, adapt with a group level also waits
# (adapt_converged) until the largest Rhat of alpha over that window is below
# adapt_converge$rhat, read from adapt_converge$min iterations on and given up
# after adapt_converge$stall checks without improvement or at
# adapt_converge$max. options(emc.adapt_converge = FALSE) turns this off.
kernel_window <- 250
adapt_converge <- list(min = 250, max = 1000, rhat = 1.2, stall = 3)

# Have the adapt draws the sample-stage kernel would be built from converged,
# or stopped improving? history: the window Rhats of the earlier checks of
# this adapt run (check_progress keeps them). Returns TRUE / FALSE with the
# history including this check as attribute "rhat".
adapt_converged <- function(emc, n_adapt, history = NULL, verbose = FALSE){
  if(legacy_sampler() || isFALSE(getOption("emc.adapt_converge")) || emc[[1]]$type == "single") return(TRUE)
  if(n_adapt < adapt_converge$min) return(structure(FALSE, rhat = history))
  h <- c(history, max(window_rhat(emc, kernel_window)))
  why <- adapt_stop(h, n_adapt)
  if(verbose) message(sprintf("  adapt, %d iterations: max Rhat of alpha over the last %d = %.3f%s", n_adapt, kernel_window, h[length(h)],
                              switch(why, limit = " (limit of adapt iterations reached)",
                                     stalled = sprintf(" (no improvement in %d checks)", adapt_converge$stall), "")))
  structure(why != "", rhat = h)
}

# The stop rule on the window Rhats h of this adapt run's checks (the last is
# the current one): "converged", "stalled", "limit" or "" (carry on).
adapt_stop <- function(h, n_adapt){
  k <- adapt_converge$stall; n <- length(h)
  if(h[n] < adapt_converge$rhat) return("converged")
  if(n > k && min(h[(n - k + 1):n]) >= min(h[1:(n - k)])) return("stalled")
  if(n_adapt >= adapt_converge$max) return("limit")
  ""
}

tune_sample_kernel <- function(emc, verbose, verboseProgress, fileName, particle_factor,
                               search_width, cores_per_chain, cores_for_chains, n_blocks,
                               on_singular, r_cores){
  if (verbose) message("Tuning the sample-stage kernel")
  n_cores <- cores_per_chain*cores_for_chains
  # a block of adapt-labelled iterations run with the sample-stage kernel
  run_block <- function(emc, iter){
    sub_emc <- subset(emc, filter = chain_n(emc)[1,"adapt"] - 1, stage = "adapt")
    sub_emc <- auto_mclapply(sub_emc, run_stages, stage = "adapt", kernel = "sample", iter = iter,
                             verbose = verbose, verboseProgress = verboseProgress,
                             particle_factor = particle_factor, search_width = search_width,
                             n_cores = cores_per_chain, mc.cores = cores_for_chains,
                             on_singular = on_singular, r_cores = r_cores)
    check_chain_failures(sub_emc, "adapt", fileName)
    class(sub_emc) <- "emc"
    if(cores_for_chains > 1) sub_emc <- fix_custom_kernel_pointers(sub_emc, emc)
    concat_emc(emc, sub_emc, iter, "adapt")
  }
  emc <- add_proposals(emc, "sample", n_cores, n_blocks, window = kernel_window)
  emc <- reset_kernel_tail(emc)
  emc <- run_block(emc, 100)
  # Freeze the local step size at the average over the tail, and the number
  # of particles where the average effective sample size meets its target
  emc <- map_pm_settings(emc, function(x){
    n_avg <- x$local_uses - 20             # local uses after the warm-up
    if(isTRUE(n_avg >= 10)) x$epsilon[1] <- exp(x$log_eps_sum / n_avg)
    if(isTRUE(x$log_ess_n >= 10)){
      n_particles <- round(x$ess_target / exp(x$log_ess_sum / x$log_ess_n))
      x$n_particles <- max(25, min(x$max_particles, n_particles))
    }
    x
  })
  map_scale_move(emc, scale_move_freeze)
}

# Restart the counters and averages that the tuning of the sample-stage kernel
# is read from (tune_sample_kernel)
reset_kernel_tail <- function(emc){
  emc <- map_pm_settings(emc, function(x){
    x <- reset_acc_counts(x); x$iter <- 25
    x$local_uses <- x$log_eps_sum <- 0
    x$log_ess_sum <- x$log_ess_n <- 0
    x
  })
  map_scale_move(emc, scale_move_reset_tail)
}

# f applied to every component of every subject's particle settings, or to
# the sweep's settings, of every chain
map_pm_settings <- function(emc, f){
  for(i in seq_along(emc)){
    pm <- attr(emc[[i]]$samples, "pm_settings")
    if(!is.null(pm)) attr(emc[[i]]$samples, "pm_settings") <- lapply(pm, function(x) lapply(x, f))
  }
  emc
}
map_scale_move <- function(emc, f){
  for(i in seq_along(emc)){
    sm <- attr(emc[[i]]$samples, "scale_move")
    if(!is.null(sm)) attr(emc[[i]]$samples, "scale_move") <- f(sm)
  }
  emc
}

reset_acc_counts <- function(x){
  x$proposal_counts <- rep(0, length(x$proposal_counts))
  x$acc_counts <- rep(0, length(x$proposal_counts))
  x
}

# Rhat (not split) of every subject x parameter alpha over the last n
# iterations; the largest entry finds the few subjects whose chains sit apart,
# which the mean (burn's criterion) hides.
window_rhat <- function(emc, n){
  idx <- emc[[1]]$samples$idx
  it <- max(1, idx - n + 1):idx
  n <- length(it)
  M <- lapply(emc, function(x) rowMeans(x$samples$alpha[, , it, drop = FALSE], dims = 2))
  V <- mapply(function(x, m) (rowMeans(x$samples$alpha[, , it, drop = FALSE]^2, dims = 2) - m^2) * n / (n - 1),
              emc, M, SIMPLIFY = FALSE)
  W <- Reduce(`+`, V) / length(V)
  Mbar <- Reduce(`+`, M) / length(M)
  B <- n * Reduce(`+`, lapply(M, function(m) (m - Mbar)^2)) / (length(M) - 1)
  r <- sqrt(((n - 1) / n * W + B / n) / W)
  r[!is.finite(r)] <- Inf
  r
}

# The interweaving sweep's gate (scale_move_gate), from the sweeps since its
# last call and the group-level draws of the last n_iter iterations
gate_scale_move <- function(emc, n_iter, decide = TRUE, verbose = FALSE){
  sm <- lapply(emc, function(x) attr(x$samples, "scale_move"))
  if(any(sapply(sm, is.null)) || emc[[1]]$type != "standard") return(emc)
  idx <- emc[[1]]$samples$idx; it <- max(2, idx - n_iter + 1):idx
  tv <- do.call(cbind, lapply(emc, function(x) apply(x$samples$theta_var[, , it, drop = FALSE], 3, diag)))
  p <- nrow(tv)
  mu <- if(is.null(emc[[1]]$group_designs) && nrow(emc[[1]]$samples$theta_mu) == p)
    do.call(cbind, lapply(emc, function(x) x$samples$theta_mu[, it, drop = FALSE])) else NULL
  sm <- scale_move_gate(sm, var_lsd = apply(.5 * log(tv), 1, stats::var),
                        var_mu = if(is.null(mu)) NULL else apply(mu, 1, stats::var),
                        mean_var = rowMeans(tv), n_subjects = emc[[1]]$n_subjects, decide = decide)
  for(i in seq_along(emc)) attr(emc[[i]]$samples, "scale_move") <- sm[[i]]
  if(verbose && decide){
    dropped <- p - sum(sm[[1]]$active_scale | sm[[1]]$active_loc)
    if(dropped > 0) message(sprintf("  sweep gate: %d of %d parameters left out of the sweep", dropped, p))
  }
  emc
}

run_stages <- function(sampler, stage = "preburn", iter=0, verbose = TRUE, verboseProgress = TRUE,
                       particle_factor=50, search_width= NULL, n_cores=1, on_singular = NULL, r_cores = 1,
                       kernel = NULL)
{
  particles <- round(particle_factor*sqrt(sampler$n_pars))
  if (!sampler$init) {
    sampler <- init(sampler, n_cores = n_cores, r_cores = r_cores)
  }
  if (iter == 0) return(sampler)
  tune <- list(search_width = search_width, kernel = kernel)
  sampler <- run_stage(sampler, stage = stage,iter = iter, particles = particles,
                       n_cores = n_cores, tune = tune, verbose = verbose,
                       verboseProgress = verboseProgress, on_singular = on_singular, r_cores = r_cores)
  return(sampler)
}

add_proposals <- function(emc, stage, n_cores, n_blocks, window = NULL){
  legacy <- legacy_sampler()
  # The sample stage's proposals are built once and then re-used at every step
  # and in every later call (see tune_sample_kernel); only the legacy sampler
  # re-estimates them from the sample draws.
  fixed <- stage == "sample" && !legacy
  if(fixed && length(emc[[1]]$sample_kernel$chains) == length(emc)){
    return(restore_sample_kernel(emc))
  }
  if(stage != "preburn"){
    # if(!is.null(emc[[1]]$g_map_fixed)){
    #   emc <- create_chain_proposals_lm(emc)
    # } else{    }
    emc <- create_chain_proposals(emc, do_block = stage != "sample")
    if(!is.null(n_blocks)){
      if(n_blocks > 1){
        components <- sub_blocking(emc, n_blocks)
        for(i in 1:length(emc)){
          attr(emc[[i]]$data, "components") <- components
        }
      }
    }
  }
  if(stage == "sample"){
    # if(!is.null(emc[[1]]$g_map_fixed)){
    #   emc <- create_eff_proposals_lm(emc, n_cores)
    # } else{    }
    emc <- create_eff_proposals(emc, n_cores, window = window)
  }
  if(stage %in% c("adapt", "sample") && !legacy && emc[[1]]$type != "single"){
    emc <- create_lik_prec(emc, n_cores)
  }
  if(fixed) emc <- store_sample_kernel(emc)
  return(emc)
}

# The sample-stage kernel's proposals are kept in the first chain's entry,
# which is the one strip_duplicates() preserves, so that they survive saving
# and a later call that adds samples.
sample_kernel_fields <- c("chains_var", "chains_mu", "eff_mu", "eff_var")

store_sample_kernel <- function(emc){
  emc[[1]]$sample_kernel <- list(
    chains = lapply(emc, function(x) c(x[sample_kernel_fields], list(prop_var = attr(x, "prop_var")))),
    lik_prec = emc[[1]]$lik_prec)
  return(emc)
}

restore_sample_kernel <- function(emc){
  kernel <- emc[[1]]$sample_kernel
  for(i in 1:length(emc)){
    for(nm in sample_kernel_fields) emc[[i]][[nm]] <- kernel$chains[[i]][[nm]]
    attr(emc[[i]], "prop_var") <- kernel$chains[[i]]$prop_var
    emc[[i]]$lik_prec <- kernel$lik_prec
  }
  return(emc)
}

# Likelihood precision of every subject from its recent draws, shared by the
# chains: by finite differences at the draws' mean (lik_precision(), the
# sweep's surrogate) and, for type "standard", from the draws' covariance
# (lik_precision_draws(), $post, used by the particle step where it exists),
# because one point's curvature is too narrow for a skewed posterior.
create_lik_prec <- function(emc, n_cores){
  idx <- emc[[1]]$samples$idx
  history_idx <- proposal_window(idx)
  alpha <- get_pars(emc, filter = history_idx, selection = "alpha",
                    stage = c('preburn', 'burn', 'adapt', 'sample'),
                    by_subject = T, merge_chains = T, return_mcmc = F,
                    remove_dup = F, remove_constants = F)
  group <- window_group_level(emc, history_idx, dim(alpha)[3])
  lik_prec <- auto_mclapply(1:emc[[1]]$n_subjects, function(sub){
    draws <- matrix(alpha[, sub, ], nrow = dim(alpha)[1], dimnames = list(dimnames(alpha)[[1]], NULL))
    lik <- tryCatch(lik_precision(rowMeans(draws), apply(draws, 1, stats::sd), emc[[1]]$data[[sub]], emc[[1]]$model),
                    error = function(e) NULL)
    if(!is.null(group) && !is.null(lik)){
      lik$post <- tryCatch(lik_precision_draws(draws, group$prec, group$mu[, sub], lik, n_chains = length(emc)),
                           error = function(e) NULL)
    }
    lik
  }, mc.cores = n_cores)
  for(i in 1:length(emc)) emc[[i]]$lik_prec <- lik_prec
  return(emc)
}

# The window's mean group precision and group means (p x n_subjects), or NULL
# where the draw-based precision does not apply: not "standard", nuisance
# parameters, < 100 draws, or a group precision that is not about constant
# over the window (largest CV >= lik_prec_draws_cv, e.g. a funnel), since
# lik_precision_draws() subtracts it.
lik_prec_draws_cv <- .5

window_group_level <- function(emc, history_idx, n_draws){
  x1 <- emc[[1]]
  if(x1$type != "standard" || any(x1$nuisance) || n_draws < 100) return(NULL)
  if(group_precision_cv(emc, history_idx) >= lik_prec_draws_cv) return(NULL)
  tryCatch({
    prec <- 0; mu <- 0
    for(x in emc){
      for(i in history_idx) prec <- prec + chol2inv(chol(x$samples$theta_var[, , i]))
      fs <- filtered_samples(x, history_idx, type = x1$type)
      mu <- mu + if(is.null(fs$subj_mu)) matrix(rowMeans(fs$theta_mu), nrow(fs$theta_mu), x1$n_subjects)
                 else apply(fs$subj_mu, 1:2, mean)
    }
    n <- length(emc) * length(history_idx)
    out <- list(prec = prec / n, mu = mu / length(emc))
    if(!all(is.finite(out$prec)) || !all(is.finite(out$mu)) || nrow(out$mu) != nrow(out$prec)) return(NULL)
    out
  }, error = function(e) NULL)
}

# The largest coefficient of variation, over the parameters, of the group
# precision 1 / sigma^2 over iterations history_idx of all chains (Inf if
# not computable)
group_precision_cv <- function(emc, history_idx){
  pr <- do.call(cbind, lapply(emc, function(x) 1 / apply(x$samples$theta_var[, , history_idx, drop = FALSE], 3, diag)))
  pr <- matrix(pr, ncol = length(emc) * length(history_idx))
  cv <- apply(pr, 1, stats::sd) / rowMeans(pr)
  if(!all(is.finite(cv))) return(Inf)
  max(cv)
}

# Likelihood precision of one subject from its draws: the inverse of their
# covariance V minus the group precision P of the same window. In coordinates
# where P is the identity the draws' precision is M = P^-1/2 V^-1 P^-1/2 and
# the likelihood's M - I. Directions where an eigenvalue of M - I exceeds tau
# (at least the largest eigenvalue's sampling noise for p parameters and the
# draws' ESS, (1 + sqrt(p / ESS))^2 - 1) take the draws' estimate; in the
# others the draws show the prior (a collapsed group level, or a likelihood
# that says little), so they take the finite-difference estimate `fd`, capped
# at tau since the draws rule out more than that.
# The linear term puts the conditional mean at the draws' mean when the group
# level is the window's (P, mu); in the directions left to `fd` it is the
# finite-difference quadratic's, taken with the other directions held at the
# draws' mean (without those cross terms the conditional proposal sits many
# posterior SDs off). draws: p x N, in n_chains blocks of equal length.
# Returns list(prec, lin, n_draws = directions taken from the draws,
# n_capped, tau, ess), or NULL (too few distinct draws, or not computable).
lik_precision_draws <- function(draws, P, mu, fd, tau = .5, n_chains = 1){
  p <- nrow(draws); N <- ncol(draws)
  if(length(unique(draws[1, ])) < max(50, 5 * p)) return(NULL)
  blk <- rep(seq_len(n_chains), each = N %/% n_chains, length.out = N)
  ess <- stats::median(sapply(seq_len(p), function(j) sum(tapply(draws[j, ], blk, function(x)
    tryCatch(as.numeric(coda::effectiveSize(x)), error = function(e) 1)))))
  tau <- max(tau, (1 + sqrt(p / max(ess, p)))^2 - 1)
  m <- rowMeans(draws)
  Vinv <- chol2inv(chol(stats::cov(t(draws))))
  e <- eigen(P, symmetric = TRUE)
  if(any(e$values <= 0)) return(NULL)
  Ph <- e$vectors %*% (sqrt(e$values) * t(e$vectors))        # P^1/2
  Pih <- e$vectors %*% (t(e$vectors) / sqrt(e$values))       # P^-1/2
  M <- Pih %*% Vinv %*% Pih
  eg <- eigen((M + t(M)) / 2, symmetric = TRUE)
  lam <- eg$values - 1
  trusted <- lam > tau
  U <- eg$vectors[, trusted, drop = FALSE]
  mt <- drop(Ph %*% m); mut <- drop(Ph %*% mu)
  Lw <- U %*% (lam[trusted] * t(U))
  lw <- drop(U %*% (eg$values[trusted] * drop(t(U) %*% mt) - drop(t(U) %*% mut)))
  n_capped <- 0
  if(!all(trusted)){
    C <- eg$vectors[, !trusted, drop = FALSE]                 # basis of the other directions
    Fw <- Pih %*% fd$prec %*% Pih
    Fc <- t(C) %*% Fw %*% C
    ef <- eigen((Fc + t(Fc)) / 2, symmetric = TRUE)
    Vc <- C %*% ef$vectors
    phi <- pmax(ef$values, 0)
    capped <- phi > tau
    n_capped <- sum(capped)
    phi[capped] <- tau
    Lw <- Lw + Vc %*% (phi * t(Vc))
    lfd <- drop(t(Vc) %*% (Pih %*% fd$lin - Fw %*% (U %*% drop(t(U) %*% mt))))
    lemp <- (phi + 1) * drop(t(Vc) %*% mt) - drop(t(Vc) %*% mut)
    lw <- lw + drop(Vc %*% ifelse(capped, lemp, lfd))
  }
  L <- Ph %*% Lw %*% Ph
  L <- (L + t(L)) / 2
  lin <- drop(Ph %*% lw)
  if(!all(is.finite(L)) || !all(is.finite(lin))) return(NULL)
  dimnames(L) <- dimnames(fd$prec); names(lin) <- names(fd$lin)
  list(prec = L, lin = lin, n_draws = sum(trusted), n_capped = n_capped)
}

check_progress <- function (emc, stage, iter, stop_criteria,
                            max_tries, step_size, n_cores, verbose, progress = NULL,
                            n_blocks)
{
  min_es <- stop_criteria$min_es
  if(is.null(min_es)) min_es <- 0
  selection <- stop_criteria$selection
  min_unique <- stop_criteria$min_unique
  total_iters_stage <- chain_n(emc)[, stage][1]
  if (is.null(progress)) {
    iters_total <- 0
    trys <- 0
  }
  else {
    iters_total <- progress$iters_total + step_size
    trys <- progress$trys + 1
    # use more informative message
    # if (verbose)
    #   message(trys, ": Iterations ", stage, " = ", total_iters_stage)
  }
  gd <- check_gd(emc, stage, stop_criteria[["max_gd"]], stop_criteria[["mean_gd"]], trys, verbose=FALSE,
                 iter = total_iters_stage, selection, omit_mpsrf = stop_criteria[["omit_mpsrf"]],
                 n_blocks, gd_quantile = stop_criteria[["gd_quantile"]])
  iter_done <- ifelse(is.null(iter) || length(iter) == 0, TRUE, total_iters_stage >= iter)
  if (min_es == 0) {
    es_done <- TRUE
  } else if (total_iters_stage != 0) {
    class(emc) <- "emc"
    curr_min_es <- Inf
    for(select in selection){
      curr_min_es <- min(c(ess_summary(emc, selection = select,
                                                stage = stage, stat_only = TRUE), curr_min_es))
    }
    # if (verbose)
    #   message("Smallest effective size = ", round(curr_min_es))
    es_done <- ifelse(!emc[[1]]$init, FALSE, curr_min_es >
                        min_es)
  }
  else {
    es_done <- FALSE
  }
  trys_done <- ifelse(is.null(max_tries), FALSE, trys >= max_tries)
  if (stage == "adapt") {
    samples_merged <- merge_chains(emc)
    test_samples <- extract_samples(samples_merged, stage = "adapt",
                                    samples_merged$samples$idx, n_chains = length(emc))
    # if(!is.null(emc[[1]]$g_map_fixed)){
    #   adapted <- test_adapted_lm(emc[[1]], test_samples, min_unique, n_cores, verbose)
    # } else{    }
    adapted <- test_adapted(emc[[1]], test_samples,
                            min_unique, n_cores, verbose)

  }
  else {
    adapted <- TRUE
  }
  # max_tries reached with only adapt's convergence rule unmet is not worth
  # the warning below: the kernel is then built from the draws there are
  enough_unique <- adapted
  adapt_rhat <- progress$adapt_rhat
  if(stage == "adapt" && adapted){
    adapted <- adapt_converged(emc, total_iters_stage, adapt_rhat, verbose)
    adapt_rhat <- attr(adapted, "rhat"); adapted <- as.vector(adapted)
  }
  done <- (es_done & iter_done & gd$gd_done & adapted) | (trys_done & iter_done)
  if(es_done & gd$gd_done & adapted & !iter_done){
    step_size <- min(step_size, abs(iter - total_iters_stage))[1]
  }
  if (trys_done & iter_done) {
    if(!(es_done & gd$gd_done & enough_unique)){
      warning("Max tries reached. If this happens in burn-in while trying to get
            gelman diagnostics small enough, you might have a particularly hard model.
            Make sure your model is well specified. If so, you can run adapt and
            sample, if run for long enough, sample usually converges eventually.")
    }
  }
  return(list(emc = gd$emc, done = done, step_size = step_size,
              trys = trys, n_blocks = gd$n_blocks, gd=gd,
              total_iters_stage=total_iters_stage, adapt_rhat = adapt_rhat,
              curr_min_es = if (min_es > 0 && total_iters_stage != 0) curr_min_es else NULL))
}

# The Rhats a stop rule reads (split_rhat()), subject-level ones from the
# stored array in one pass (gd_summary() costs one call per subject), in
# gd_summary()'s order.
stage_gds <- function(emc, selection, stage, omit_mpsrf = TRUE, filter = 0){
  gd_out <- c(); alpha <- NULL; alpha_all <- NULL; other <- c()
  for(select in selection){
    if(select == "alpha" && omit_mpsrf){
      its <- lapply(emc, function(x){
        it <- which(x$samples$stage[seq_len(x$samples$idx)] == stage)
        if(filter > 0) it <- it[-seq_len(min(filter, length(it)))]
        it
      })
      n <- min(lengths(its))
      p <- dim(emc[[1]]$samples$alpha)[1]; ns <- dim(emc[[1]]$samples$alpha)[2]
      if(n < 4){
        gd <- rep(NaN, p * ns); is_const <- rep(FALSE, p * ns)
      } else{
        X <- vapply(seq_along(emc), function(i){
          a <- emc[[i]]$samples$alpha[, , its[[i]][seq_len(n)], drop = FALSE]
          t(matrix(a, p * ns, n))                           # draws x (subjects, parameters within)
        }, matrix(0, n, p * ns))
        # an entry with one value in every draw of every chain is not a sampled
        # quantity: left out of the criterion, as get_pars() leaves it out
        is_const <- rowSums(colSums(abs(X - rep(X[1, , 1], each = n)), dims = 1)) == 0
        gd <- c(t(matrix(split_rhat(X), p, ns)))            # parameters, subjects within
        is_const <- c(t(matrix(is_const, p, ns)))
      }
      gd[is.na(gd)] <- Inf
      alpha_all <- gd
      gd <- gd[!is_const]
      alpha <- gd
    } else{
      # get_pars() keeps draws filter:n, subset() (the discard) drops the first
      # filter: one more here, so that both read the draws that would be kept
      gd <- unlist(gd_summary.emc(emc, selection = select, stage = stage, filter = if(filter > 0) filter + 1 else 0,
                                  omit_mpsrf = omit_mpsrf, stat = NULL, digits = 6))
      gd[is.na(gd)] <- Inf
      if(select == "alpha") alpha <- alpha_all <- gd else other <- c(other, c(gd))
    }
    gd_out <- c(gd_out, c(gd))
  }
  # alpha_all: every subject-level entry in its place, for set_tune_ess()
  list(gd = gd_out, alpha = alpha, other = other, alpha_all = alpha_all)
}

# max_gd's statistic: the largest Rhat, or with gd_quantile that quantile of
# the alpha Rhats and the largest of the rest (see ?fit).
gd_top <- function(g, gd_quantile = NULL){
  if(is.null(gd_quantile) || is.null(g$alpha) || !all(is.finite(g$alpha))) return(max(g$gd))
  max(as.numeric(stats::quantile(g$alpha, gd_quantile)), g$other)
}

check_gd <- function(emc, stage, max_gd, mean_gd, omit_mpsrf, trys, verbose,
                     selection, iter, n_blocks = 1, gd_quantile = NULL)
{
  if(is.null(max_gd) & is.null(mean_gd)) return(list(gd_done = TRUE, emc = emc))
  if(!emc[[1]]$init | !stage %in% emc[[1]]$samples$stage)
    return(list(gd_done = FALSE, emc = emc))
  if(is.null(omit_mpsrf)) omit_mpsrf <- TRUE
  gd_ok <- function(g){
    ok_max <- if(is.null(max_gd)) TRUE else all(is.finite(g$gd)) && gd_top(g, gd_quantile) < max_gd
    ok_mean <- if(is.null(mean_gd)) TRUE else all(is.finite(g$gd)) && mean(g$gd) < mean_gd
    ok_max & ok_mean
  }
  g <- stage_gds(emc, selection, stage, omit_mpsrf)
  ok_gd <- gd_ok(g)
  if(!ok_gd) {
    # Chains that are still moving into the posterior: if the diagnostic is
    # better without the first third of the stage's draws, those are dropped
    # for good. Decided on the largest Rhat whatever gd_quantile is: the
    # largest is the sensitive detector of what is left of a transient.
    n_remove <- round(chain_n(emc)[,stage][1]/3)
    g_short <- tryCatch(stage_gds(emc, selection, stage, omit_mpsrf, filter = n_remove), error = function(e) NULL)
    if(!is.null(g_short) &&
       ((is.null(max_gd) && mean(g_short$gd) < mean(g$gd)) || (!is.null(max_gd) && max(g_short$gd) < max(g$gd)))){
      emc_short <- try(subset.emc(emc, filter=n_remove,stage=stage, keep_stages = TRUE), silent = TRUE)
      if(!is(emc_short, "try-error")){
        g <- g_short
        emc <- emc_short
        ok_gd <- gd_ok(g)
      }
    }
  }
  gd <- g$gd
  if(verbose) {
    type <- "Rhat"
    if (!is.null(mean_gd)) message("Mean ",type," = ",round(mean(gd),3)) else
      if (!is.null(max_gd)) message("Max ",type," = ",round(max(gd),3))
  }
  emc <- set_tune_ess(emc, g$alpha_all, mean_gd, max_gd)
  class(emc) <- "emc"
  return(list(gd = gd, gd_done = ok_gd, emc = emc))
}

set_tune_ess <- function(emc, alpha_gd = NULL, mean_gd = NULL, max_gd = NULL){
  if(is.null(alpha_gd)){
    return(emc)
  }
  alpha_gd <- matrix(alpha_gd, emc[[1]]$n_subjects)
  if(!is.null(mean_gd)){
    mean_alpha_ok <- rowMeans(alpha_gd) < 1+(mean_gd-1)*.5
  } else{
    mean_alpha_ok <- rep(T, nrow(alpha_gd))
  }
  if(!is.null(max_gd)){
    max_alpha_ok <- apply(alpha_gd, 1, max) < 1+(max_gd-1)*.5
  } else{
    max_alpha_ok <- rep(T, nrow(alpha_gd))
  }
  alpha_ok <- max_alpha_ok & mean_alpha_ok
  emc <- lapply(emc, function(x){
    pm_settings <- attr(x$samples, "pm_settings")
    pm_settings <- mapply(pm_settings, alpha_ok, FUN = function(y, z){
      for(i in 1:length(y)){
        y[[i]]$gd_good <- z
      }
      return(list(y))
    })
    attr(x$samples, "pm_settings") <- pm_settings
    return(x)
  })
  return(emc)
}


# window: build them from the last `window` adapt iterations of each chain
# only (the draws adapt_converged() judged), or from as many more as the joint
# covariance of a subject's parameters and the group level needs (twice its
# dimension in draws).
create_eff_proposals <- function(emc, n_cores, window = NULL){
  samples_merged <- merge_chains(emc)
  full_filter <- NULL
  if(!is.null(window) && chain_n(emc)[1, "sample"] == 0){
    n_it <- emc[[1]]$samples$idx
    n_adapt <- chain_n(emc)[1, "adapt"]
    p <- sum(!emc[[1]]$nuisance)
    window <- max(window, ceiling(2 * (2 * p + p * (p + 1) / 2) / length(emc)))
    if(window < n_adapt){
      it <- (n_it - window + 1):n_it
      full_filter <- unlist(lapply(seq_along(emc), function(i) (i - 1) * n_it + it))
    }
  }
  test_samples <- extract_samples(samples_merged, stage = c("adapt", "sample"), max_n_sample = 750, n_chains = length(emc),
                                  full_filter = full_filter)

  type <- emc[[1]]$type
  components <- attr(emc[[1]]$data, "components")
  for(i in 1:length(emc)){
    iteration = round(test_samples$iteration * i/length(emc))
    n_pars <- emc[[1]]$n_pars
    nuisance <- emc[[1]]$nuisance
    n_subjects <- emc[[1]]$n_subjects
    eff_mu <- matrix(0, nrow = n_pars, ncol = n_subjects)
    eff_var <- array(0, dim = c(n_pars, n_pars, n_subjects))
    for(comp in unique(components)){
      idx <- comp == components
      nuis_idx <- nuisance[idx]
      if(any(nuis_idx)){
        type <- samples_merged$sampler_nuis$type
        conditionals <- auto_mclapply(X = 1:n_subjects,
                                      FUN = get_conditionals, samples = test_samples,
                                      n_pars = sum(idx[!nuisance]), iteration =  iteration, idx = idx[!nuisance],
                                      type = type ,
                                      mc.cores = n_cores)
        conditionals_nuis <- auto_mclapply(X = 1:n_subjects,
                                           FUN = get_conditionals, samples = test_samples$nuisance,
                                           n_pars = sum(idx[nuisance]), iteration =  iteration, idx = idx[nuisance],
                                           type = type,
                                           mc.cores = n_cores)
        conditionals <- array(unlist(conditionals), dim = c(sum(idx[!nuisance]), sum(idx[!nuisance]) + 1, n_subjects))
        conditionals_nuis <- array(unlist(conditionals_nuis), dim = c(sum(idx[nuisance]), sum(idx[nuisance]) + 1, n_subjects))
        eff_mu[idx & !nuisance,] <- conditionals[,1,]
        eff_var[idx & !nuisance, idx & !nuisance,] <- conditionals[,2:(sum(idx[!nuisance])+1),]
        eff_mu[idx & nuisance,] <- conditionals_nuis[,1,]
        eff_var[idx & nuisance,idx & nuisance,] <- conditionals_nuis[,2:(sum(idx[nuisance])+1),]
      } else{
        conditionals <- auto_mclapply(X = 1:n_subjects,
                                      FUN = get_conditionals, samples = test_samples,
                                      n_pars = sum(idx[!nuisance]), iteration =  iteration, idx = idx[!nuisance],
                                      type = type,
                                      mc.cores = n_cores)
        conditionals <- array(unlist(conditionals), dim = c(sum(idx[!nuisance]), sum(idx[!nuisance]) + 1, n_subjects))
        eff_mu[idx & !nuisance,] <- conditionals[,1,]
        eff_var[idx & !nuisance,idx & !nuisance,] <- conditionals[,2:(sum(idx[!nuisance])+1),]
      }

    }
    # eff_mu <- lapply(conditionals, FUN = function(x) x$eff_mu)
    # eff_var <- lapply(conditionals, FUN = function(x) x$eff_var)
    # eff_alpha <- lapply(conditionals, FUN = function(x) x$eff_alpha)
    # eff_tau <- lapply(conditionals, FUN = function(x) x$eff_tau)

    eff_mu <- split(eff_mu, col(eff_mu))
    eff_var <- apply(eff_var, 3, identity, simplify = F)
    emc[[i]]$eff_mu <- eff_mu
    emc[[i]]$eff_var <- eff_var
    # attr(emc[[i]], "eff_alpha") <- eff_alpha
    # attr(emc[[i]], "eff_tau") <- eff_tau
  }
  return(emc)
}


sub_blocking <- function(emc, n_blocks){
  covs <- lapply(emc, FUN = function(x){return(x$chains_var)})
  out <- array(0, dim = dim(covs[[1]][[1]]))
  for(i in 1:length(covs)){
    cov_tmp <- covs[[1]]
    for(j in 1:length(cov_tmp)){
      out <- out + cov2cor(cov_tmp[[1]])
    }
  }
  shared_ll_idx <- attr(emc[[1]]$data, "shared_ll_idx")
  min_comp <- 0
  components <- c()
  for(ll in unique(shared_ll_idx)){
    idx <- ll == shared_ll_idx
    distance <-as.dist(1- abs(out[idx, idx]/out[1,1]))
    clusts <- hclust(distance)
    sub_comps <- min_comp + cutree(clusts, k = n_blocks) # This could go wrong if one group has just one member
    min_comp <- max(sub_comps)
    components <- c(components, sub_comps)
  }
  return(components)
}

# The draws the proposals are built from: the last 250 iterations, or the last
# third once there are fewer than 375
proposal_window <- function(idx) unique(pmax(1, round(idx - min(250, idx / 1.5)):idx - 1))

create_chain_proposals <- function(emc, samples_idx = NULL, do_block = TRUE){
  n_subjects <- emc[[1]]$n_subjects
  n_chains <- length(emc)
  n_pars <- emc[[1]]$n_pars
  stage <- emc[[1]]$samples$stage[length(emc[[1]]$samples$stage)]
  history_idx <- if(is.null(samples_idx)) proposal_window(emc[[1]]$samples$idx) else unique(pmax(1, samples_idx - 1))
  LL <- get_pars(emc, filter = history_idx, selection = "LL",
                 stage = c('preburn', 'burn', 'adapt', 'sample'),
                 merge_chains = T, return_mcmc = F, remove_constants = F,
                 remove_dup = F)
  alpha <- get_pars(emc, filter = history_idx, selection = "alpha",
                    stage = c('preburn', 'burn', 'adapt', 'sample'),
                    by_subject = T, merge_chains = T, return_mcmc = F,
                    remove_dup = F, remove_constants = F)

  components <- attr(emc[[1]]$data, "components")
  block_idx <- block_variance_idx(components)
  for(j in 1:n_chains){
    # Take only a random half of all the samples for each chain
    rnd_index <- sample(1:ncol(LL), round(ncol(LL)/2))
    chains_var <- vector("list", n_subjects)
    chains_mu <- vector("list", n_subjects)
    for(sub in 1:n_subjects){
      moments <- weighted_moments(alpha[,sub,rnd_index], LL[sub, rnd_index])
      emp_covs <- moments$w_cov
      chains_mu[[sub]] <- moments$w_mu
      if(do_block) emp_covs[block_idx] <- 0
      if(!is.positive.definite(emp_covs)){
        # If not positive definite (e.g. too few distinct draws), do not use it
        next
      } else{
        chains_var[[sub]] <- emp_covs
      }
    }
    # Subjects without a usable covariance: use the mean of the other
    # subjects', else this chain's previous one, else a scaled-down prior
    # variance (its epsilon then adapts quickly). A fixed diag(.5) can be
    # orders of magnitude too wide for a concentrated posterior.
    null_idx <- sapply(chains_var, is.null)
    if(all(null_idx)){
      mean_chains_var <- if(legacy_sampler()) diag(n_pars) * .5 else NULL
    } else{
      mean_chains_var <- Reduce(`+`, chains_var[!null_idx]) / sum(!null_idx)
    }
    if(any(null_idx)){
      prev <- emc[[j]]$chains_var
      for(q in 1:n_subjects){
        if(null_idx[q]){
          chains_var[[q]] <- if(!is.null(mean_chains_var)) mean_chains_var
          else if(!is.null(prev) && !is.null(prev[[q]])) prev[[q]]
          else diag(diag(emc[[1]]$prior$theta_mu_var), n_pars) * .1
        }
      }
      if(is.null(mean_chains_var)) mean_chains_var <- Reduce(`+`, chains_var) / n_subjects
    }

    new_prop_var <- mean(diag(mean_chains_var))
    if(is.null(attr(emc[[j]], "prop_var"))){
      # This is in preburn case, group-level proposals are wider
      # So we scale the epsilon a bit to account for narrow individual proposals
      prop_var_ratio <- 2
    } else{
      # epsilon scales standard deviations, prop_var is a variance
      prop_var_ratio <- attr(emc[[j]], "prop_var")/new_prop_var
      if(!legacy_sampler()) prop_var_ratio <- sqrt(prop_var_ratio)
    }
    # Within adapt the hierarchical local kernel's epsilon is relative to
    # the likelihood + group precision (new_particle), not to chains_var
    lik_scaled <- stage == "adapt" && !legacy_sampler() && emc[[1]]$type != "single"
    if(stage != "sample" && !lik_scaled){
      emc[[j]] <- update_epsilon_scale(emc[[j]], prop_var_ratio)
    }
    attr(emc[[j]], "prop_var") <- new_prop_var
    emc[[j]]$chains_var <- chains_var
    emc[[j]]$chains_mu <- chains_mu
  }
  return(emc)
}

reset_pm_settings <- function(emc, stage){
  new_stage <- stage != get_last_stage(emc)
  legacy <- legacy_sampler()
  # Acceptance counts restart every step (step_size iterations), so the
  # epsilon and mixing-weight adaptation follows the recent window rather
  # than the whole stage's history.
  if(legacy && !(new_stage || stage == "burn")) return(emc)
  # The sample stage's kernel is fixed: it keeps the settings that the tail of
  # adapt tuned for it (tune_sample_kernel)
  if(!legacy && stage == "sample") return(emc)
  map_pm_settings(emc, function(x){
    x <- reset_acc_counts(x)
    if(new_stage){
      x$iter <- 25
      x$mix <- NULL
      # burn's epsilons belong to (prior-variance, chains_var) random
      # walks; adapt's first scaled component is the chains_var one, so
      # carry that epsilon (check_epsilon then pads the vector)
      if(stage == "adapt" && !legacy) x$epsilon <- x$epsilon[length(x$epsilon)]
    }
    x
  })
}

update_epsilon_scale <- function(pmwgs, prop_var_ratio){
  pm_settings <- attr(pmwgs$samples, "pm_settings")
  pm_settings <- lapply(pm_settings, function(x){
    for(i in 1:length(x)){
      x[[i]]$epsilon <- x[[i]]$epsilon * prop_var_ratio
    }
    return(x)
  })
  attr(pmwgs$samples, "pm_settings") <- pm_settings
  return(pmwgs)
}

test_adapted <- function(sampler, test_samples, min_unique, n_cores_conditional = 1,
                         verbose = FALSE)
{
  # Function used by run_adapt to check whether we can create the conditional.

  # Only need to check uniqueness for one parameter
  first_par <- as.matrix(test_samples$alpha[1, , ])
  if(ncol(first_par) == 1) first_par <- t(first_par)
  # Split the matrix into a list of vectors by subject
  # all subjects is greater than unq_vals
  n_unique_sub <- apply(first_par, 1, FUN = function(x) return(length(unique(x))))
  n_pars <- sampler$n_pars
  components <- attr(sampler$data, "components")
  nuisance <- sampler$nuisance

  if (length(n_unique_sub) != 0 & all(n_unique_sub > min_unique)) {
    if(verbose){
      message("Enough unique values detected: ", min_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
      for(comp in unique(components)){
        idx <- comp == components
        nuis_idx <- nuisance[idx]
        if(any(nuis_idx)){
          type <- sampler$sampler_nuis$type
          auto_mclapply(X = 1:sampler$n_subjects,
                        FUN = get_conditionals, samples = test_samples$nuisance,
                        n_pars = sum(idx[nuisance]), idx = idx[nuisance],
                        type = type,
                        mc.cores = n_cores_conditional)
        }
        auto_mclapply(X = 1:sampler$n_subjects,FUN = get_conditionals,samples = test_samples,
                      n_pars = sum(idx[!nuisance]), idx = idx[!nuisance], type = sampler$type,
                      mc.cores = n_cores_conditional)
      }
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
      if(verbose){
        message("Can't create efficient distribution yet")
        message("Increasing required unique values and continuing adaptation")
      }
      return(FALSE)
    }
    else {
      if(verbose) message("Successfully adapted - stopping adaptation")
      return(TRUE)
    }
  } else{
    return(FALSE) # Not enough unique particles found
  }
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# Apply a design's TC (truncation/censoring) list to data that have not been
# through make_missing(). No-op when TC requests nothing or the data already
# carry censoring columns.
apply_design_TC <- function(data, design, rt_resolution = NULL) {
  TC <- design$TC
  if (is.null(TC) || !is.data.frame(data)) return(data)
  if (any(c("LT", "UT", "LC", "UC", "missingness") %in% names(data))) return(data)
  is_default <- function(x, default) is.null(x) || (is.numeric(x) && all(x == default))
  if (is_default(TC$LT, 0) && is_default(TC$LC, 0) && is_default(TC$UT, Inf) &&
      is_default(TC$UC, Inf) && is_default(TC$pContaminant, 0)) return(data)
  message("Applying the design's truncation/censoring (TC) to the data")
  make_missing(data, LT = TC$LT, LC = TC$LC, UC = TC$UC, UT = TC$UT,
               LCresponse = TC$LCresponse, UCresponse = TC$UCresponse,
               LCdirection = TC$LCdirection, UCdirection = TC$UCdirection,
               pContaminant = TC$pContaminant,
               no_truncate = TC$no_truncate, no_censor = TC$no_censor,
               verbose = isTRUE(TC$verbose), rt_resolution = rt_resolution,
               digits = if (is.null(TC$digits)) 2 else TC$digits)
}


#' Make an emc Object
#'
#' Creates an emc object by combining the data, prior,
#' and model specification into a `emc` object that is needed in `fit()`.
#'
#' @param data A data frame, or a list of data frames. Needs to have the variable `subjects` as participant identifier.
#' @param design A list with a pre-specified design, the output of `design()`.
#' @param model A model list. If none is supplied, the model specified in `design()` is used.
#' @param type A string indicating whether to run a `standard` group-level, `blocked`, `diagonal`, `factor`, or `single` (i.e., non-hierarchical) model.
#' @param n_chains An integer. Specifies the number of mcmc chains to be run (has to be more than 1 to compute `rhat`).
#' @param compress A Boolean, if `TRUE` (i.e., the default), the data is compressed to speed up likelihood calculations.
#' @param rt_resolution A double. Used for compression, response times will be binned based on this resolution.
#' @param group_design A design for group-level mappings, made using `group_design()`.
#' @param par_groups A vector. Indicates which parameters are allowed to correlate. Could either be a list of character vectors of covariance blocks. Or
#' a numeric vector, e.g., `c(1,1,1,2,2)` means the covariances
#' of the first three and of the last two parameters are estimated as two separate blocks.
#' @param prior_list A named list containing the prior. Default prior created if `NULL`. For the default priors, see `?get_prior_{type}`.
#' @param memory_saver A Boolean. If `TRUE`, store a pooled design representation and drop per-parameter designs from data to reduce memory usage.
#' @param ... Additional, optional arguments.
#' @return An uninitialized emc object
#' @examples dat <- forstmann
#'
#' # function that takes the lR factor (named diff in the following function) and
#' # returns a logical defining the correct response for each stimulus. In this
#' # case the match is simply such that the S factor equals the latent response factor.
#' matchfun <- function(d)d$S==d$lR
#'
#' # design an "average and difference" contrast matrix
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))
#'
#' # specify design
#' design_LBABE <- design(data = dat,model=LBA,matchfun=matchfun,
#' formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
#' contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
#'
#' # specify priors
#' pmean <- c(v=1,v_lMdiff=1,sv_lMTRUE=log(.5), B=log(.5),B_Eneutral=log(1.5),
#'            B_Eaccuracy=log(2),B_lRright=0, A=log(0.25),t0=log(.2))
#' psd <- c(v=1,v_lMdiff=0.5,sv_lMTRUE=.5,
#'          B=0.3,B_Eneutral=0.3,B_Eaccuracy=0.3,B_lRright=0.3,A=0.4,t0=.5)
#' prior_LBABE <- prior(design_LBABE, type = 'standard',pmean=pmean,psd=psd)
#'
#' # create emc object
#' LBABE <- make_emc(dat,design_LBABE,type="standard",  prior=prior_LBABE,
#'                   compress = FALSE)
#'
#' @export

make_emc <- function(data,design,model=NULL,
                    type="standard",
                    n_chains=3,compress=TRUE,rt_resolution=1/60,
                    prior_list = NULL, group_design = NULL,
                    par_groups=NULL, memory_saver = FALSE, ...){
  # arguments for future compatibility
  n_factors <- NULL
  nuisance <- NULL
  nuisance_non_hyper <- NULL
  sem_settings <- NULL
  Lambda_mat <- NULL
  # overwrite those that were supplied
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if(!is.null(prior_list) & !is.null(prior_list$theta_mu_mean)){
    prior_list <- list(prior_list)
  }
  if(!is.null(prior_list)){
    type <- attr(prior_list[[1]], "type")
  }
  if(!is.null(group_design) && !type %in% c("standard", "diagonal", "blocked")){
    stop("group_design can only be used with standard, blocked or diagonal type")
  }
  if(type != "single" && length(unique(data$subjects)) == 1){
    stop("can only use type = `single` when there's only one subject in the data")
  }

  if (!(type %in% c("standard","diagonal","blocked","factor","single", "infnt_factor", "SEM", "diagonal-gamma")))
    stop("type must be one of: standard,diagonal,blocked,factor,infnt_factor, single")

  if(!is.null(nuisance) & !is.null(nuisance_non_hyper)){
    stop("You can only specify nuisance OR nuisance_non_hyper")
  }
  if (is(data, "data.frame")) data <- list(data)
  data <- lapply(data,function(d){
    d$subjects <- factor(d$subjects)
    d <- d[order(d$subjects),]
    LC <- attr(d,"LC")
    UC <- attr(d,"UC")
    LT <- attr(d,"LT")
    UT <- attr(d,"UT")
    d <- add_trials(d)
    attr(d,"LC") <- LC
    attr(d,"UC") <- UC
    attr(d,"LT") <- LT
    attr(d,"UT") <- UT
    d
  })
  if (!is.null(names(design)[1]) && names(design)[1]=="Flist"){
    design <- list(design)
  }
  checks <- sapply(design, function(x) is(x, "emc.design"))
  if (!all(checks)) stop("design must be a list of emc.design objects")
  if (length(design)!=length(data)){
    design <- rep(design,length(data))
  }
  if (is.null(model)) model <- lapply(design,function(x){x$model})
  if (any(unlist(lapply(model,is.null))))
    stop("Must supply model if model is not in all design components")
  if (!is.null(names(model)[1]) && names(model)[1]=="type")
    model <- list(model)
  if (length(model)!=length(data))
    model <- rep(model,length(data))

  ## SM: check for delta rules in trend, and override/turn off compression if the user supplied compress=TRUE.
  compress_passed <- compress
  compress <- rep(compress, length(model))
  has_delta_rule <- sapply(model, has_delta_rules)
  compress[has_delta_rule] <- FALSE
  if(compress_passed & any(has_delta_rule)) {
    if(length(model) == 1) message('Because the model contains a delta rule, data will not be compressed.')
    else message(paste0('Models ', which(has_delta_rule), ' contain a delta rule; the corresponding data will not be compressed.'))
    rt_resolution <- 0.001   # no need to downsample resolution when not compressing
  }
  ## SM END

  dadm_list <- vector(mode="list",length=length(data))
  rt_resolution <- rep(rt_resolution,length.out=length(data))
  for (i in 1:length(dadm_list)) {
    message("Processing data set ",i)
    if(is.null(attr(design[[i]], "custom_ll"))){
      # Censoring/truncation requested in design(TC = ...) but the data carry no
      # LT/UT/LC/UC/missingness columns (i.e. make_missing() was not applied,
      # as for real data): apply it now rather than silently fit without it.
      data[[i]] <- apply_design_TC(data[[i]], design[[i]], rt_resolution[i])
      dadm_list[[i]] <- design_model(data=data[[i]],design=design[[i]],
                                     compress=compress[[i]],model=model[[i]],rt_resolution=rt_resolution[i],
                                     memory_saver = memory_saver,
                                     check_identifiability = TRUE)
      sampled_p_names <- names(attr(design[[i]],"p_vector"))
    } else{
      if (memory_saver) {
        warning("memory_saver not supported for custom likelihoods; ignored")
      }
      dadm_list[[i]] <- design_model_custom_ll(data = data[[i]],
                                               design = design[[i]],model=model[[i]])
      sampled_p_names <- attr(design[[i]],"sampled_p_names")
    }
    if(length(prior_list) == length(data)){
      if(!is.null(prior_list[[i]])){
        prior_list[[i]] <- check_prior(prior_list[[i]], sampled_p_names, group_design)
      }
    }
  }
  # Warn before fitting about parameters (or linear combinations) that are
  # unidentified across the (joint) model, so a fit that would later diverge on
  # an ill-conditioned group covariance fails fast with a clear message instead.
  warn_identifiability(dadm_list, names(design))

  # Make sure class retains following changes
  class(design) <- "emc.design"
  prior_in <- merge_priors(prior_list)

  prior_in <- prior(design, type, update = prior_in, group_design = group_design, ...)
  attr(dadm_list[[1]], "prior") <- prior_in

  # if(!is.null(subject_covariates)) attr(dadm_list, "subject_covariates") <- subject_covariates
  if (type %in% c("single", "infnt_factor", "diagonal-gamma")) {
    out <- pmwgs(dadm_list, type, nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper,
                 n_factors = n_factors)
  } else if (type %in% c("standard", "blocked", "diagonal")) {
    if(type == "blocked"){
      if(is.null(par_groups)) stop("par_groups must be specified for blocked models")
    }
    if(type == "diagonal" && is.null(par_groups)){
      par_groups <- 1:length(sampled_pars(design))
    }
    if(type == "standard" && is.null(par_groups)){
      par_groups <- rep(1, length(sampled_pars(design)))
    }
    if(!is.null(par_groups)){
      if(is.character(par_groups)){
        par_groups <- list(par_groups)
      }
      if(is.list(par_groups)){
        par_names <- names(sampled_pars(design))
        new_par_groups <- rep(NA, length(par_names))
        for(i in 1:length(par_groups)){
          if(any(!par_groups[[i]] %in% par_names)) stop("Make sure you specified parameter names in par_groups correctly")
          new_par_groups[par_names %in% par_groups[[i]]] <- i
        }
        new_par_groups[is.na(new_par_groups)] <- (i+1):(i+sum(is.na(new_par_groups)))
        par_groups <- new_par_groups
      }
      if(length(par_groups) != length(sampled_pars(design))){
        stop("par_groups length does not match number of sampled parameters, make sure you specified par_groups correctly")
      }
    }
    if(type %in% c("diagonal", "blocked")) type <- "standard"
    out <- pmwgs(dadm_list, type, par_groups=par_groups,
                 group_design = group_design,
                 nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper)
  } else if (type == "factor") {
    out <- pmwgs(dadm_list, type, n_factors = n_factors,
                 nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper,
                 Lambda_mat = Lambda_mat)
  } else if (type == "SEM"){
    out <- pmwgs(dadm_list, type, sem_settings = sem_settings, nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper)
  }
  out$model <- lapply(design, function(x) x$model)
  # Only for joint models we need to keep a list of functions
  if(length(out$model) == 1) out$model <- out$model[[1]]
  out <- check_duplicate_designs(out)
  # replicate chains
  dadm_lists <- rep(list(out),n_chains)
  # For post predict
  class(dadm_lists) <- "emc"
  return(dadm_lists)
}

fix_fileName <- function(x){
  ext <- substr(x, nchar(x)-5, nchar(x))
  if(ext != ".RData" & ext != ".Rdata"){
    return(paste0(x, ".RData"))
  } else{
    return(x)
  }
}

check_duplicate_designs <- function(out){
  if(is.data.frame(out$data[[1]])) return(out)
  if(!is.null(attr(out$data[[1]][[1]], "design_pool"))) return(out)
  for(i in 1:length(out$data)){ # loop over subjects
    designs <- lapply(out$data[[i]], function(y) attr(y, "designs"))
    duplicacy <- duplicated(designs)
    unq_idx <-   sapply(seq_along(designs), function(i) {
      for (j in seq_along(designs)) {
        if (identical(designs[[i]], designs[[j]])) return(j)
      }
    })
    for(j in 1:length(out$data[[i]])){# Loop over data sets in this sub
      if(is.null(designs[[j]])) next
      if(duplicacy[j]){
        attr(out$data[[i]][[j]], "designs") <- unq_idx[j]
      }
    }
  }
  return(out)
}

extractDadms <- function(dadms, names = NULL){
  if(is.null(names)) names <- 1:length(dadms)
  N_models <- length(dadms)
  pars <- attr(dadms[[1]], "sampled_p_names")
  prior <- attr(dadms[[1]], "prior")
  subjects <- unique(factor(sapply(dadms, FUN = function(x) levels(x$subjects))))
  dadm_list <- dm_list(dadms[[1]])
  components <- rep(1, length(pars))
  if(N_models > 1){
    total_dadm_list <- vector("list", length = N_models)
    k <- 1
    pars <- paste(names[1], pars, sep = "|")
    dadm_list[as.character(which(!subjects %in% unique(dadms[[1]]$subjects)))] <- NA
    total_dadm_list[[1]] <- dadm_list
    for(dadm in dadms[-1]){
      k <- k + 1
      tmp_list <- vector("list", length = length(subjects))
      tmp_list[as.numeric(unique(dadm$subjects))] <- dm_list(dadm)
      total_dadm_list[[k]] <- tmp_list
      curr_pars <- attr(dadm, "sampled_p_names")
      components <- c(components, rep(k, length(curr_pars)))
    }
    dadm_list <- do.call(mapply, c(list, total_dadm_list, SIMPLIFY = F))
  }
  # subject_covariates_ok <- unlist(lapply(subject_covariates, FUN = function(x) length(x) == length(subjects)))
  # if(!is.null(subject_covariates_ok)) if(any(!subject_covariates_ok)) stop("subject_covariates must be as long as the number of subjects")
  attr(dadm_list, "components") <- components
  attr(dadm_list, "shared_ll_idx") <- components
  return(list(prior = prior,
              dadm_list = dadm_list, subjects = subjects))
}

auto_mclapply <- function(X, FUN, mc.cores, ...){
  if(Sys.info()[1] == "Windows") return(cluster_lapply(X, FUN, mc.cores, ...))
  parallel::mclapply(X, FUN, mc.cores = mc.cores, ...)
}

# Socket workers start without the session's options, and the sampler is
# switched by options(emc.*), so those are copied to the workers.
cluster_lapply <- function(X, FUN, cores, ...){
  cluster <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster, options, options()[grep("^emc\\.", names(options()))])
  parallel::parLapply(cl = cluster, X, FUN, ...)
}

# Turn a silent per-chain sampler failure into an informative error. mclapply
# puts a try-error object in a chain's slot when its worker errors (and a worker
# killed for memory returns something with no $samples); left unchecked this
# surfaces much later as "$ operator is invalid for atomic vectors" in chain_n().
check_chain_failures <- function(chains, stage, fileName = NULL){
  failed <- vapply(chains, function(x){
    inherits(x, "try-error") || is.null(x) || !is.list(x) || is.null(x$samples)
  }, logical(1))
  if(!any(failed)) return(invisible(NULL))
  reasons <- vapply(which(failed), function(i){
    x <- chains[[i]]
    cond <- attr(x, "condition")
    reason <- if(!is.null(cond)) conditionMessage(cond)
      else if(inherits(x, "try-error")) trimws(as.character(x))
      else "chain worker returned no samples (it likely crashed or ran out of memory)"
    paste0("  chain ", i, ": ", reason)
  }, character(1))
  saved <- if(!is.null(fileName)) paste0(" The fit up to the previous stage was saved to '",
                                         fix_fileName(fileName), "'.") else ""
  # Tailor the advice to what actually failed, rather than always blaming the
  # covariance. The per-chain reasons now carry context (iteration + on_singular
  # recovery so far + diverging parameters) added in run_stage.
  singular_like <- any(grepl("singular|reciprocal condition|not positive|leading minor|positive definite|ill-conditioned|Lapack",
                             reasons, ignore.case = TRUE))
  advice <- if (singular_like) {
    paste0("The group-level covariance became ill-conditioned, usually from an ",
           "unidentified or weakly-identified parameter (see the diverging parameters ",
           "listed above). Options: check for unidentified parameters, use a tighter ",
           "prior or fewer free parameters, or set the `on_singular` argument of fit() ",
           "to retry / carry_forward.")
  } else {
    paste0("This is not the known singular-covariance failure but some other error ",
           "(shown above) hit inside a chain. Read the underlying message and context; ",
           "if it looks like a bug rather than your model/data, please report it.")
  }
  # Wrap the advice to ~80 columns so it is readable; the saved-file note gets
  # its own line.
  advice <- paste(strwrap(advice, width = 80), collapse = "\n")
  if (nzchar(saved)) advice <- paste0(advice, "\n", trimws(saved))
  stop("Sampling failed in ", sum(failed), " of ", length(chains),
       " chain(s) during the '", stage, "' stage:\n",
       paste(reasons, collapse = "\n"),
       "\n\n", advice, call. = FALSE)
}

#' Strip all entries except samples from EMC list entries
#' @param emc A list of EMC objects
#' @return The same list with everything but samples removed from all but first entry
#' @noRd
strip_duplicates <- function(emc, incl_props = TRUE) {
  # Keep only samples in non-first entries
  for (i in 2:length(emc)) {
    samples <- emc[[i]]$samples
    prop_var <- attr(emc[[i]], "prop_var")
    emc[[i]] <- list(samples = samples)
    attr(emc[[i]], "prop_var") <- prop_var
  }
  if(incl_props){
    # Also remove eff_mu, eff_var, chains_cov
    for (i in 1:length(emc)) {
      emc[[i]]$eff_mu <- NULL
      emc[[i]]$eff_var <- NULL
      emc[[i]]$chains_cov <- NULL
      emc[[i]]$chains_mu <- NULL
      emc[[i]]$lik_prec <- NULL
    }
  }
  return(emc)
}

get_posterior_weights <- function(ll){
  max_ll <- max(ll)
  weights <- exp(ll - max_ll)
  weights <- weights / sum(weights)
  return(weights)
}

weighted_moments <- function(chain, ll = NULL) {
  # chain: a matrix with each col as a sample, and each row as a parameter
  # weights: an optional vector of weights. If not provided, equal weighting is assumed.
  if (is.null(dim(chain))) {
    chain <- matrix(chain, nrow = 1L)
  }
  n <- ncol(chain)
  d <- nrow(chain)

  # Use equal weights if none are provided.
  if (is.null(ll)) {
    weights <- rep(1 / n, n)
  } else {
    weights <- get_posterior_weights(ll)
  }

  # Compute the weighted mean of the chain.
  weighted_mean <- drop(chain %*% weights)
  # NIEK THIS CAUSES ERRORS
  # # Compute the weighted covariance matrix.
  # cov_matrix <- matrix(0, nrow = d, ncol = d)
  # for (i in 1:n) {
  #   diff <- chain[,i] - weighted_mean
  #   cov_matrix <- cov_matrix + weights[i] * (diff %*% t(diff))
  # }
  cov_matrix <- if (d == 1L) {
    matrix(stats::var(as.numeric(chain)), nrow = 1L, ncol = 1L)
  } else {
    cov(t(chain))
  }
  return(list(w_cov = cov_matrix, w_mu = weighted_mean))
}



#' Restore all entries to EMC list entries from first entry
#' @param emc A list of EMC objects with stripped duplicates
#' @return The same list with all fields restored from first entry except samples
#' @noRd
restore_duplicates <- function(emc) {
  # Restore everything except samples from first entry
  for (i in 2:length(emc)) {
    samples <- emc[[i]]$samples
    emc[[i]] <- emc[[1]]
    emc[[i]]$samples <- samples
  }
  return(emc)
}
