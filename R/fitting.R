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
  if(!is.null(stop_criteria[["min_es"]]) && is.null(stop_criteria$selection)){
    stop_criteria$selection <- if(type != "single") c('alpha', 'mu') else 'alpha'
  }
  # Chain culling (burn and sample) and the flat test (sample); see ?fit
  if(!is.null(stop_criteria$cull)){
    if(!stage %in% c("burn", "sample")) stop("`cull` applies to the burn and sample stages")
    given <- if(is.list(stop_criteria$cull)) names(stop_criteria$cull) else character(0)
    stop_criteria$cull <- cull_defaults(stop_criteria$cull)
    # In burn a false cull costs nothing but one chain's diversity, so the
    # comparison starts earlier there than in sample, where it restarts the stage.
    if(!is.null(stop_criteria$cull) && stage == "burn" && !"min_window" %in% given)
      stop_criteria$cull$min_window <- 250L
  }
  if(!is.null(stop_criteria$max_flat_loc)){
    if(stage != "sample") stop("`max_flat_loc` applies to the sample stage")
    if(is.null(stop_criteria$flat_selection)) stop_criteria$flat_selection <- c("mu", "alpha", "LL", "subj_ll")
    if(type == "single") stop_criteria$flat_selection <- setdiff(stop_criteria$flat_selection, "mu")
    if(is.null(stop_criteria$flat_p1)) stop_criteria$flat_p1 <- 1/3
    if(is.null(stop_criteria$flat_p2)) stop_criteria$flat_p2 <- 1/3
    if(is.null(stop_criteria$flat_quantile)) stop_criteria$flat_quantile <- c(mu = 1, alpha = .99, LL = 1, subj_ll = .99)
    if(is.null(stop_criteria$flat_mcse)) stop_criteria$flat_mcse <- 0
    if(stop_criteria$flat_p1 <= 0 || stop_criteria$flat_p1 > 1 || stop_criteria$flat_p2 <= 0 || stop_criteria$flat_p2 > 1)
      stop("flat_p1 and flat_p2 must be in (0, 1]")
    if(!all(stop_criteria$flat_selection %in% c("mu", "alpha", "LL", "subj_ll")))
      stop("flat_selection must be one or more of: mu, alpha, LL, subj_ll")
  }
  # Opt-in sample-stage criteria on the per-check statistics (sample_stop_stats()):
  # an effective-size floor on the group means and log group SDs, and vetoes on
  # the chains' log-likelihood gap and the subject spread. The statistics are
  # reported at every check whether or not these are set.
  for(k in c("min_ess", "max_ll_gap", "max_subj_spread")){
    v <- stop_criteria[[k]]
    if(is.null(v)) next
    if(stage != "sample") stop("`", k, "` applies to the sample stage")
    if(!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0) stop("`", k, "` must be a single non-negative number")
    if(k == "min_ess" && type == "single") stop("`min_ess` reads the group means: not for type 'single'")
  }
  if(!is.null(stop_criteria$max_sample_iter)){
    if(stage != "sample") stop("`max_sample_iter` applies to the sample stage")
    if(stop_criteria$max_sample_iter < 10) stop("max_sample_iter must be at least 10")
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
#' @param verboseProgress Logical. Whether to print progress within each sampling step.
#' @param fileName A string. If specified will autosave emc at this location on every iteration.
#' @param particle_factor An integer. `particle_factor` multiplied by the square root of the number of sampled parameters determines the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain.
#' Parallelizes across participant calculations using native threads.
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
#' covariance becomes computationally singular during sampling. `NULL`: up to 3 re-draws of the
#' group step, then an error.
#' A list may set `max_retries` (integer, re-draw the group step on failure).
#' @export
#' @param .context Internal compiled likelihood cache.
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
                    thin = FALSE, trim = TRUE, on_singular = NULL, .context = NULL){
  if(length(emc) == 1) stop("run_emc() currently requires n_chains > 1.")
  emc <- restore_duplicates(emc)
  if (is.null(.context)) .context <- sampler_contexts(emc, cores_per_chain)
  if (verbose) message(paste0("Running ", stage, " stage"))
  total_iters_stage <- chain_n(emc)[,stage][1]
  if(stage != "preburn"){
    iter <- stop_criteria[["iter"]] + total_iters_stage
    if(!is.null(stop_criteria$max_sample_iter) && iter > stop_criteria$max_sample_iter)
      stop("`max_sample_iter` is below the iterations requested, which could never be reached")
  } else{
    iter <- stop_criteria[["iter"]]
  }
  progress <- check_progress(emc, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains, verbose)
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
    sub_emc <- run_native_chains(sub_emc, stage = stage, iter= progress$step_size*max(1,cur_thin),
                             verbose=verbose,  verboseProgress = verboseProgress,
                             particle_factor=particle_factor,search_width=search_width,
                             n_cores=cores_per_chain, cores_for_chains = cores_for_chains,
                             on_singular = on_singular, context = .context)

    class(sub_emc) <- "emc"
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
    if(stage == "sample" && !is.null(stop_criteria$max_sample_iter)){
      emc <- trim_sample_stage(emc, stop_criteria$max_sample_iter)
    }
    progress <- check_progress(emc, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains,
                               verbose, progress)
    emc <- progress$emc
    progress <- progress[!names(progress) == 'emc']
    # Chain culling: a chain (or one subject in one chain) held in a worse
    # mode than the others is respawned from a better chain; in sample the
    # stage then restarts, with the tries counter, after a short re-burn.
    if(!progress$done && !is.null(stop_criteria$cull) && emc[[1]]$init){
      cull <- stop_criteria$cull
      n_culls <- if(is.null(progress$n_culls)) 0L else progress$n_culls
      if(n_culls < cull$max_culls){
        diag <- cull_diagnostics(emc, stage, cull)
        progress$cull_state <- update_cull_state(progress$cull_state, diag, length(emc))
        act <- apply_culls(emc, diag, stage, cull, progress$cull_state)
        if(length(act$done)){
          emc <- act$emc
          emc[[1]]$culls <- c(emc[[1]]$culls, act$done)
          progress$n_culls <- n_culls + 1L
          progress$cull_state <- NULL
          if(verbose) for(d in act$done) message(sprintf("Cull %d of %d: %s of chain %d respawned from chain %d (log posterior %.0f below, z = %.1f)",
            progress$n_culls, cull$max_culls, if(is.na(d$subject)) "the state" else paste("subject", d$subject), d$chain, d$donor, d$gap, d$z))
          if(stage == "sample"){
            emc <- cull_restart(emc, cull$reburn, verbose, verboseProgress, particle_factor, search_width,
                                cores_per_chain, cores_for_chains, n_blocks, on_singular, .context)
            iter <- stop_criteria[["iter"]] + chain_n(emc)[1, stage]
            progress$trys <- 0L; progress$iters_total <- 0
            if(verbose) message("Sample stage restarted after the cull")
          }
        }
      }
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
      if(!is.null(stop_criteria$max_gd)) gd_message <- rhat_message(progress$gd$rule, gd)
      # get min effective sample size (ESS)
      ess_message <- NULL
      if(stage == 'sample' & !is.null(progress$curr_min_es)) {
        ess_message <- paste0(" | min ESS=", round(progress$curr_min_es))
      }
      ess_message <- paste0(ess_message, stop_stats_message(progress$stop_stats))
      if(!is.null(progress$flat) && length(progress$flat)){
        ess_message <- paste0(ess_message, " | drift(se)=", paste(sprintf("%s %.2f(%.2f)", names(progress$flat),
          vapply(progress$flat, function(f) f$drift, 0), vapply(progress$flat, function(f) f$mcse, 0)), collapse = ", "))
        wess <- vapply(progress$flat, function(f) if(is.null(f$ess) || !is.finite(f$ess)) Inf else f$ess, 0)
        if(any(wess < 20 * length(emc))) ess_message <- paste0(ess_message, sprintf(" | window ESS %.0f: too short to tell a wave from a trend", min(wess)))
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

  # The fit's record of how the sample stage ended: the rules that were on, their
  # values at the last check, and whether the rules or max_tries ended it
  if(stage == "sample"){
    on <- c("max_gd", "mean_gd", "gd_quantile", "iter", "min_ess", "max_ll_gap", "max_subj_spread", "max_flat_loc", "flat_mcse", "min_es")
    emc[[1]]$stop_record <- list(criteria = stop_criteria[intersect(on, names(stop_criteria))], cull = stop_criteria$cull,
                                 met = isTRUE(progress$met), tries = progress$trys, rule = progress$gd$rule,
                                 stats = progress$stop_stats, culls = length(emc[[1]]$culls))
  }
  if(stage == "adapt" && is.null(emc[[1]]$sample_kernel) &&
     get_last_stage(emc) == "adapt"){
    emc <- tune_sample_kernel(emc, verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName, particle_factor = particle_factor, search_width = search_width,
                              cores_per_chain = cores_per_chain, cores_for_chains = cores_for_chains,
                              n_blocks = n_blocks, on_singular = on_singular, context = .context)
  }

  emc <- strip_duplicates(emc)
  class(emc) <- "emc"
  return(emc)
}

# Build the fixed sample kernel from adapt's tail, then tune it for 100 discarded draws.
kernel_window <- 250
adapt_converge <- list(min = 250, max = 1000, rhat = 1.2, stall = 3)

# Stop when the proposal window converges, stalls, or reaches the adapt limit.
adapt_converged <- function(emc, n_adapt, history = NULL){
  if(isFALSE(getOption("emc.adapt_converge")) || emc[[1]]$type == "single") return(TRUE)
  if(n_adapt < adapt_converge$min) return(structure(FALSE, rhat = history))
  h <- c(history, max(window_rhat(emc, kernel_window)))
  structure(adapt_stop(h, n_adapt), rhat = h)
}

adapt_stop <- function(h, n_adapt){
  k <- adapt_converge$stall; n <- length(h)
  h[n] < adapt_converge$rhat || n_adapt >= adapt_converge$max ||
    (n > k && min(tail(h, k)) >= min(head(h, -k)))
}

tune_sample_kernel <- function(emc, verbose, verboseProgress, fileName, particle_factor,
                               search_width, cores_per_chain, cores_for_chains, n_blocks,
                               on_singular, context){
  n_cores <- cores_per_chain*cores_for_chains
  emc <- add_proposals(emc, "sample", n_cores, n_blocks, window = kernel_window)
  emc <- map_pm_settings(emc, function(x){
    x <- reset_acc_counts(x); x$iter <- 25
    x$local_uses <- x$log_eps_sum <- x$log_ess_sum <- x$log_ess_n <- 0
    x
  })
  emc <- map_scale_move(emc, scale_move_reset_tail)
  sub_emc <- subset(emc, filter = chain_n(emc)[1,"adapt"] - 1, stage = "adapt")
  sub_emc <- run_native_chains(sub_emc, stage = "adapt", kernel = "sample", iter = 100,
                           verbose = verbose, verboseProgress = verboseProgress,
                           particle_factor = particle_factor, search_width = search_width,
                           n_cores = cores_per_chain, cores_for_chains = cores_for_chains,
                           on_singular = on_singular, context = context)
  class(sub_emc) <- "emc"
  emc <- concat_emc(emc, sub_emc, 100, "adapt")
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
  d <- dim(emc[[1]]$samples$alpha)
  r <- split_rhat(alpha_draws(emc, rep(list(it), length(emc))), split = FALSE)
  r[is.na(r)] <- Inf
  matrix(r, d[1], d[2])
}

# The alpha draws its[[i]] of chain i (the first min(lengths(its)) of each) as
# draws x (subjects, parameters within) x chains, split_rhat()'s layout
alpha_draws <- function(emc, its){
  n <- min(lengths(its)); d <- dim(emc[[1]]$samples$alpha)
  vapply(seq_along(emc), function(i){
    a <- emc[[i]]$samples$alpha[, , its[[i]][seq_len(n)], drop = FALSE]
    t(matrix(a, d[1] * d[2], n))
  }, matrix(0, n, d[1] * d[2]))
}

add_proposals <- function(emc, stage, n_cores, n_blocks, window = NULL){
  # Reuse the frozen proposals, including when a saved fit is resumed.
  fixed <- stage == "sample"
  if(fixed && length(emc[[1]]$sample_kernel$chains) == length(emc)){
    return(restore_sample_kernel(emc))
  }
  if(stage != "preburn"){
    emc <- create_chain_proposals(emc, do_block = stage != "sample")
    if(!is.null(n_blocks) && n_blocks > 1){
      components <- sub_blocking(emc, n_blocks)
      for(i in 1:length(emc)){
        attr(emc[[i]]$data, "components") <- components
      }
    }
  }
  if(stage == "sample"){
    emc <- create_eff_proposals(emc, n_cores, window = window)
  }
  if(stage %in% c("burn", "adapt", "sample") && emc[[1]]$type != "single"){
    # The joint location move's directions are rebuilt only when the stage
    # changes (and when the sample kernel is tuned, which starts a "sample"
    # kernel on an adapt object), so that their step sizes keep adapting
    # across the steps of a stage instead of restarting at every step.
    emc <- create_lik_prec(emc, n_cores, reset_directions = stage != get_last_stage(emc))
  }
  if(fixed) emc <- store_sample_kernel(emc)
  return(emc)
}

# The first chain survives strip_duplicates(), so store the fixed kernel there.
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

create_lik_prec <- function(emc, n_cores, reset_directions = TRUE){
  history_idx <- proposal_window(emc[[1]]$samples$idx)
  lik_prec <- native_lik_prec(emc, history_idx, native_likelihood_inputs(emc[[1]]), n_cores)
  for (i in seq_along(emc)) {
    emc[[i]]$lik_prec <- lik_prec
    sm <- attr(emc[[i]]$samples, "scale_move")
    if (reset_directions && !is.null(sm)) {
      sm$directions <- NULL
      attr(emc[[i]]$samples, "scale_move") <- sm
    }
  }
  emc
}

check_progress <- function (emc, stage, iter, stop_criteria,
                            max_tries, step_size, n_cores, verbose, progress = NULL)
{
  min_es <- stop_criteria[["min_es"]]
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
  }
  cull_state <- progress$cull_state; n_culls <- progress$n_culls
  gd <- check_gd(emc, stage, stop_criteria$max_gd, stop_criteria$mean_gd,
                 omit_mpsrf = stop_criteria$omit_mpsrf, selection = selection,
                 gd_quantile = stop_criteria$gd_quantile)
  # Everything below reads the draws the Rhat read: check_gd() may have dropped
  # the first third of the stage, and the draws that remain are what the run
  # keeps, so the iteration target, the ESS and the flat test count those
  emc <- gd$emc
  total_iters_stage <- chain_n(emc)[, stage][1]
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
    adapted <- adapt_converged(emc, total_iters_stage, adapt_rhat)
    adapt_rhat <- attr(adapted, "rhat"); adapted <- as.vector(adapted)
  }
  # The per-check statistics and the opt-in criteria on them (see sample_stop_stats())
  stop_stats <- if(stage == "sample") sample_stop_stats(emc, stage) else NULL
  crit_ok <- function(k, stat, below){
    v <- stop_criteria[[k]]
    if(is.null(v)) return(TRUE)
    s <- stop_stats[[stat]]
    if(is.null(s) || !is.finite(s)) return(FALSE)
    if(below) s <= v else s >= v
  }
  stats_ok <- crit_ok("min_ess", "ess_min", FALSE) & crit_ok("max_ll_gap", "ll_gap", TRUE) &
    crit_ok("max_subj_spread", "subj_spread", TRUE)
  # The flat test: the stage must not still be drifting between its first and last thirds
  flat <- NULL; flat_ok <- TRUE
  if(stage == "sample" && !is.null(stop_criteria$max_flat_loc) && total_iters_stage >= 30){
    flat <- flat_diagnostics(emc, stage, stop_criteria$flat_selection, stop_criteria$flat_p1,
                             stop_criteria$flat_p2, stop_criteria$flat_quantile)
    flat_ok <- flat_done(flat, stop_criteria$max_flat_loc, stop_criteria$flat_mcse)
  } else if(stage == "sample" && !is.null(stop_criteria$max_flat_loc)) flat_ok <- FALSE
  # In sample, max_tries binds even when the iteration target is unmet: a
  # discard can take the stage back below the target at every failing check,
  # so without this the run would continue past max_tries indefinitely
  done <- (es_done & iter_done & gd$gd_done & adapted & flat_ok & stats_ok) | (trys_done & (iter_done | stage == "sample"))
  if(es_done & gd$gd_done & adapted & flat_ok & stats_ok & !iter_done){
    step_size <- min(step_size, abs(iter - total_iters_stage))[1]
  }
  if (trys_done & (iter_done | stage == "sample")) {
    if(!(es_done & gd$gd_done & enough_unique & flat_ok & stats_ok & iter_done)){
      warning("Max tries reached. If this happens in burn-in while trying to get
            gelman diagnostics small enough, you might have a particularly hard model.
            Make sure your model is well specified. If so, you can run adapt and
            sample, if run for long enough, sample usually converges eventually.")
    }
  }
  return(list(emc = emc, done = done, step_size = step_size,
              trys = trys, gd=gd,
              total_iters_stage=total_iters_stage, adapt_rhat = adapt_rhat,
              curr_min_es = if (min_es > 0 && total_iters_stage != 0) curr_min_es else NULL,
              flat = flat, cull_state = cull_state, n_culls = n_culls, stop_stats = stop_stats,
              met = es_done & iter_done & gd$gd_done & adapted & flat_ok & stats_ok))
}

# The Rhats a stop rule reads (split_rhat()), subject-level ones from the
# stored array in one pass (gd_summary() costs one call per subject), in
# gd_summary()'s order.
stage_gds <- function(emc, selection, stage, omit_mpsrf = TRUE, filter = 0){
  gd_out <- c(); alpha <- NULL; other <- c()
  for(select in selection){
    if(select == "alpha"){
      its <- lapply(emc, function(x){
        it <- which(x$samples$stage[seq_len(x$samples$idx)] == stage)
        if(filter > 0) it <- it[-seq_len(min(filter, length(it)))]
        it
      })
      n <- min(lengths(its))
      p <- dim(emc[[1]]$samples$alpha)[1]; ns <- dim(emc[[1]]$samples$alpha)[2]
      if(n < 4){
        gd <- rep(NaN, p * ns); is_const <- rep(FALSE, p * ns)
        mpsrf <- rep(Inf, p)
      } else{
        X <- alpha_draws(emc, its)
        # an entry with one value in every draw of every chain is not a sampled
        # quantity: left out of the criterion, as get_pars() leaves it out
        is_const <- rowSums(colSums(abs(X - rep(X[1, , 1], each = n)), dims = 1)) == 0
        # One multivariate PSRF per parameter, across subjects.
        mpsrf <- if(omit_mpsrf) NULL else vapply(seq_len(p), function(j){
          cols <- (seq_len(ns) - 1) * p + j
          cols <- cols[!is_const[cols]]
          if(!length(cols)) return(NA_real_)
          gelman_diag_robust(lapply(seq_along(emc), function(i) coda::mcmc(matrix(X[, cols, i], n))),
                             omit_mpsrf = FALSE)[["mpsrf"]]
        }, 0)
        gd <- c(t(matrix(split_rhat(X), p, ns)))            # parameters, subjects within
        is_const <- c(t(matrix(is_const, p, ns)))
      }
      gd[is.na(gd)] <- Inf
      gd <- gd[!is_const]
      alpha <- gd
      if(!omit_mpsrf){
        mpsrf <- mpsrf[!is.na(mpsrf)]
        other <- c(other, mpsrf)
        gd <- c(gd, mpsrf)
      }
    } else{
      # get_pars() keeps draws filter:n, subset() (the discard) drops the first
      # filter: one more here, so that both read the draws that would be kept
      gd <- unlist(gd_summary.emc(emc, selection = select, stage = stage, filter = if(filter > 0) filter + 1 else 0,
                                  omit_mpsrf = omit_mpsrf, stat = NULL, digits = 6))
      gd[is.na(gd)] <- Inf
      other <- c(other, c(gd))
    }
    gd_out <- c(gd_out, c(gd))
  }
  list(gd = gd_out, alpha = alpha, other = other)
}

# The progress line's Rhat: the rule's value (what max_gd is compared with) and,
# when it differs, the plain maximum over every Rhat the rule read
rhat_message <- function(rule, gd){
  if(is.null(rule) || !length(gd)) return(NULL)
  m <- max(gd)
  if(isTRUE(all.equal(rule, m))) sprintf("Rhat=%.3f", rule) else sprintf("Rhat=%.3f (rule; max %.3f)", rule, m)
}

# Apply the quantile only to subject parameters; all other Rhats must pass.
gd_top <- function(g, gd_quantile = NULL){
  if(is.null(gd_quantile) || is.null(g$alpha) || !all(is.finite(g$alpha))) return(max(g$gd))
  max(stats::quantile(g$alpha, gd_quantile), g$other)
}

check_gd <- function(emc, stage, max_gd, mean_gd, omit_mpsrf,
                     selection, gd_quantile = NULL)
{
  if(is.null(max_gd) & is.null(mean_gd)) return(list(gd_done = TRUE, emc = emc))
  if(!emc[[1]]$init | !stage %in% emc[[1]]$samples$stage)
    return(list(gd_done = FALSE, emc = emc))
  if(is.null(omit_mpsrf)) omit_mpsrf <- TRUE
  gd_ok <- function(g){
    all(is.finite(g$gd)) &&
      (is.null(max_gd) || gd_top(g, gd_quantile) < max_gd) &&
      (is.null(mean_gd) || mean(g$gd) < mean_gd)
  }
  g <- stage_gds(emc, selection, stage, omit_mpsrf)
  ok_gd <- gd_ok(g)
  if(!ok_gd) {
    # Drop the first third if this improves the stopping diagnostic.
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
  class(emc) <- "emc"
  # the value the rule compares with max_gd (mean_gd): the maximum over the
  # non-subject Rhats and the gd_quantile quantile of the subject ones
  rule <- if(!is.null(max_gd)) gd_top(g, gd_quantile) else mean(gd)
  return(list(gd = gd, gd_done = ok_gd, emc = emc, rule = rule))
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
    p <- emc[[1]]$n_pars
    window <- max(window, ceiling(2 * (2 * p + p * (p + 1) / 2) / length(emc)))
    if(window < n_adapt){
      it <- (n_it - window + 1):n_it
      full_filter <- unlist(lapply(seq_along(emc), function(i) (i - 1) * n_it + it))
    }
  }
  test_samples <- extract_samples(samples_merged, stage = c("adapt", "sample"), max_n_sample = 750, n_chains = length(emc),
                                  full_filter = full_filter)
  native_eff_proposals(emc, test_samples)
}

sub_blocking <- function(emc, n_blocks){
  native_sub_blocking(emc, n_blocks)
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
  native_chain_proposals(emc, alpha, LL, do_block, stage)
}

reset_pm_settings <- function(emc, stage){
  new_stage <- stage != get_last_stage(emc)
  # Acceptance counts restart every step (step_size iterations), so the
  # epsilon and mixing-weight adaptation follows the recent window rather
  # than the whole stage's history.
  # The sample stage's kernel is fixed: it keeps the settings that the tail of
  # adapt tuned for it (tune_sample_kernel)
  if(stage == "sample") return(emc)
  map_pm_settings(emc, function(x){
    x <- reset_acc_counts(x)
    if(new_stage){
      x$iter <- 25
      x$mix <- NULL
      # burn's epsilons belong to (prior-variance, chains_var) random
      # walks; adapt's first scaled component is the chains_var one, so
      # carry that epsilon (check_epsilon then pads the vector)
      if(stage == "adapt") x$epsilon <- x$epsilon[length(x$epsilon)]
    }
    x
  })
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

  if (length(n_unique_sub) != 0 & all(n_unique_sub > min_unique)) {
    if(verbose){
      message("Enough unique values detected: ", min_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
      for(comp in unique(components)){
        idx <- comp == components
        lapply(seq_len(sampler$n_subjects), function(s)
          native_conditionals(test_samples, sampler$type, s, test_samples$iteration, idx))
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
#' @param type A string indicating whether to run a `standard`, `blocked`, `diagonal`, `diagonal-gamma`, `SEM`, `infnt_factor`, or `single` (non-hierarchical) model.
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
  sem_settings <- NULL
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

  if (!(type %in% c("standard","diagonal","blocked","single", "infnt_factor", "SEM", "diagonal-gamma")))
    stop("type must be one of: standard,diagonal,blocked,diagonal-gamma,SEM,infnt_factor,single")

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
    out <- pmwgs(dadm_list, type, n_factors = n_factors)
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
                 group_design = group_design)
  } else if (type == "SEM"){
    out <- pmwgs(dadm_list, type, sem_settings = sem_settings)
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
  if(mc.cores > 1 && length(X) > 1 && RNGkind()[1] == "L'Ecuyer-CMRG"){
    # Reserve separate streams for chains.
    # Advancing the parent also prevents repeated streams in later batches.
    seed <- get(".Random.seed", envir = .GlobalEnv)
    seeds <- vector("list", length(X))
    for(i in seq_along(X)) seeds[[i]] <- seed <- parallel::nextRNGStream(seed)
    assign(".Random.seed", parallel::nextRNGStream(seed), envir = .GlobalEnv)
    return(parallel::mclapply(setNames(seq_along(X), names(X)), function(i, ...){
      assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
      FUN(X[[i]], ...)
    }, ..., mc.cores = mc.cores, mc.set.seed = FALSE))
  }
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

#' Strip all entries except samples from EMC list entries
#' @param emc A list of EMC objects
#' @return The same list with everything but samples removed from all but first entry
#' @noRd
strip_duplicates <- function(emc, incl_props = TRUE) {
  # Keep only samples in non-first entries
  for (i in seq_along(emc)[-1L]) {
    samples <- emc[[i]]$samples
    prop_var <- attr(emc[[i]], "prop_var")
    emc[[i]] <- list(samples = samples)
    attr(emc[[i]], "prop_var") <- prop_var
  }
  if(incl_props){
    # Also remove proposal distributions.
    for (i in 1:length(emc)) {
      emc[[i]]$eff_mu <- NULL
      emc[[i]]$eff_var <- NULL
      emc[[i]]$chains_var <- NULL
      emc[[i]]$chains_mu <- NULL
      emc[[i]]$lik_prec <- NULL
    }
  }
  return(emc)
}



#' Restore all entries to EMC list entries from first entry
#' @param emc A list of EMC objects with stripped duplicates
#' @return The same list with all fields restored from first entry except samples
#' @noRd
restore_duplicates <- function(emc) {
  # Restore everything except samples from first entry
  for (i in seq_along(emc)[-1L]) {
    samples <- emc[[i]]$samples
    emc[[i]] <- emc[[1]]
    emc[[i]]$samples <- samples
  }
  return(emc)
}


# Keep at most max_sample_iter draws of the sample stage (the flat test's window)
trim_sample_stage <- function(emc, max_sample_iter){
  n <- chain_n(emc)[1, "sample"]
  if(n <= max_sample_iter) return(emc)
  subset(emc, stage = "sample", filter = n - max_sample_iter, keep_stages = TRUE)
}

# After a cull in the sample stage: keep only the respawned last draw, run
# `reburn` discarded iterations with the fixed sample kernel, and keep only
# the last draw again, so the stage starts afresh from settled states.
cull_restart <- function(emc, reburn, verbose, verboseProgress, particle_factor, search_width,
                         cores_per_chain, cores_for_chains, n_blocks, on_singular, context){
  last_only <- function(emc) subset(emc, stage = "sample", filter = chain_n(emc)[1, "sample"] - 1, keep_stages = TRUE)
  emc <- last_only(emc)
  if(reburn > 0){
    emc <- add_proposals(emc, "sample", cores_per_chain * cores_for_chains, n_blocks)
    sub_emc <- run_native_chains(emc, stage = "sample", iter = reburn, verbose = verbose,
                                 verboseProgress = verboseProgress, particle_factor = particle_factor,
                                 search_width = search_width, n_cores = cores_per_chain,
                                 cores_for_chains = cores_for_chains, on_singular = on_singular, context = context)
    class(sub_emc) <- "emc"
    emc <- concat_emc(emc, sub_emc, reburn, "sample")
    emc <- last_only(emc)
  }
  class(emc) <- "emc"
  emc
}
