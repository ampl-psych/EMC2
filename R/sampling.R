pmwgs <- function(dadm, type, pars = NULL, prior = NULL,
                  nuisance = NULL, nuisance_non_hyper = NULL, ...) {
  if(is.data.frame(dadm)) dadm <- list(dadm)
  if(is.null(pars)) pars <- names(sampled_pars(attr(dadm[[1]], "prior")))
  if(is.null(prior)) prior <- attr(dadm[[1]], "prior")
  dadm <- extractDadms(dadm)

  dadm_list <-dadm$dadm_list
  # Storage for the samples.
  subjects <- sort(as.numeric(unique(dadm$subjects)))
  if(!is.null(nuisance) & !is.numeric(nuisance)) nuisance <- which(pars %in% nuisance)
  if(!is.null(nuisance_non_hyper) & !is.numeric(nuisance_non_hyper)) nuisance_non_hyper <- which(pars %in% nuisance_non_hyper)

  if(!is.null(nuisance_non_hyper)){
    is_nuisance <- is.element(seq_len(length(pars)), nuisance_non_hyper)
    nuis_type <- "single"
  } else if(!is.null(nuisance)) {
    is_nuisance <- is.element(seq_len(length(pars)), nuisance)
    nuis_type <- "diagonal"
  } else{
    is_nuisance <- rep(F, length(pars))
  }


  sampler_nuis <- NULL
  if(any(is_nuisance)){
    sampler_nuis <- list(
      samples = sample_store(dadm, pars, nuis_type, integrate = F,
                                                    is_nuisance = !is_nuisance, ...),
      n_subjects = length(subjects),
      n_pars = sum(is_nuisance),
      nuisance = rep(F, sum(is_nuisance)),
      type = nuis_type
    )
    if(nuis_type == "single") sampler_nuis$samples <- NULL
    sampler_nuis <- add_info(sampler_nuis, prior$prior_nuis, nuis_type, ...)
  }
  samples <- sample_store(dadm, pars, type, is_nuisance = is_nuisance, ...)
  sampler <- list(
    data = dadm_list,
    par_names = pars,
    subjects = subjects,
    n_pars = length(pars),
    nuisance = is_nuisance,
    n_subjects = length(subjects),
    samples = samples,
    sampler_nuis = sampler_nuis,
    type = type,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler <- add_info(sampler, prior, type, ...)
  return(sampler)
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 verbose = FALSE, particles = 1000,
                 n_cores = 1, r_cores = 1) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  type <- pmwgs$type
  startpoints <-startpoints_comb <- get_startpoints(pmwgs, start_mu, start_var, type)
  if(any(pmwgs$nuisance)){
    type_nuis <- pmwgs$sampler_nuis$type
    startpoints_nuis <- get_startpoints(pmwgs$sampler_nuis, start_mu = NULL, start_var = NULL, type = type_nuis)
    startpoints_comb <- merge_group_level(startpoints$tmu, startpoints_nuis$tmu,
                                          startpoints$tvar, startpoints_nuis$tvar,
                                          pmwgs$nuisance, startpoints$subj_mu)
    pmwgs$sampler_nuis$samples <- fill_samples(samples = pmwgs$sampler_nuis$samples,
                                                                      group_level = startpoints_nuis,
                                                                      j = 1,
                                                                      proposals = NULL,
                                                                      n_pars = pmwgs$n_pars, type = type_nuis)
    pmwgs$sampler_nuis$samples$idx <- 1
  }
  proposals <- auto_mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,
                                  parameters = startpoints_comb, n_particles = particles,
                                  pmwgs = pmwgs, type = type,
                                  mc.cores = n_cores, r_cores = r_cores, rng_substream = TRUE)
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 1, pmwgs$n_subjects))

  # Sample the mixture variables' initial values.

  pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = startpoints, proposals = proposals,
                                             j = 1, n_pars = pmwgs$n_pars, type = type)
  pmwgs$init <- TRUE
  return(pmwgs)
}

#' Initialize Chains
#'
#' Adds a set of start points to each chain. These start points are sampled from a user-defined multivariate
#' normal across subjects.
#'
#' @param emc An emc object made by `make_emc()`
#' @param start_mu A vector. Mean of multivariate normal used in proposal distribution
#' @param start_var A matrix. Variance covariance matrix of multivariate normal used in proposal distribution.
#' Smaller values will lead to less deviation around the mean.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across participant calculations.
#' @param cores_for_chains An integer. How many cores to use to parallelize across chains. Default is the number of chains.
#' @param particles An integer. Number of starting values
#' @param ... optional additional arguments
#'
#' @return An emc object
#' @examples \donttest{
#' # Make a design and an emc object
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1),
#'                            constants=c(s=log(1)))
#'
#' DDMaE <- make_emc(forstmann, design_DDMaE, compress = FALSE)
#' # set up our mean starting points (same used across subjects).
#' mu <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'        t0=log(.2))
#' # Small variances to simulate start points from a tight range
#' var <- diag(0.05, length(mu))
#' # Initialize chains, 4 cores per chain, and parallelizing across our 3 chains as well
#' # so 4*3 cores used.
#' DDMaE <- init_chains(DDMaE, start_mu = mu, start_var = var,
#'                      cores_per_chain = 1, cores_for_chains = 1, particles = 3)
#' # Afterwards we can just use fit
#' # DDMaE <- fit(DDMaE, cores_per_chain = 4)
#' }
#' @export
init_chains <- function(emc, start_mu = NULL, start_var = NULL, particles = 1000,
                        cores_per_chain=1,cores_for_chains = length(emc),
                        ...)
{
  dots <- add_defaults(list(...),r_cores=1)
  emc <- auto_mclapply(emc,init,start_mu = start_mu, start_var = start_var,
           verbose = FALSE, particles = particles,r_cores=dots$r_cores,
           n_cores = cores_per_chain, mc.cores=cores_for_chains)
  class(emc) <- "emc"
  return(emc)
}

start_proposals <- function(s, parameters, n_particles, pmwgs, type, r_cores = 1){
  #Draw the first start point
  group_pars <- get_group_level(parameters, s, type)
  proposals <- particle_draws(n_particles, group_pars$mu, group_pars$var)
  colnames(proposals) <- rownames(pmwgs$samples$alpha) # preserve par names
  lw <- calc_ll_manager(proposals, dadm = pmwgs$data[[which(pmwgs$subjects == s)]],
                        model = pmwgs$model, r_cores = r_cores)
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}


check_tune_settings <- function(tune, n_pars, stage, particles){
  # Acceptance ratio tuning
  tune$alphaStar <- ifelse(stage == "sample", 2, 3)
  tune$p_accept <- set_p_accept(stage, tune$search_width)
  tune$local <- local_components(stage)
  # Potential blocking settings
  if(is.null(tune$components)) tune$components <- rep(1, n_pars)
  if(is.null(tune$shared_ll_idx)) tune$shared_ll_idx <- tune$components
  # Tuning of number of particles, might be a bit arbitrary
  if(is.null(tune$target_ESS)) tune$target_ESS <- 2.5*sqrt(n_pars)
  if(is.null(tune$ESS_scale)) tune$ESS_scale <- .05
  if(is.null(tune$max_particles)) tune$max_particles <- particles*1.2
  # Mix tuning settings
  if(is.null(tune$mix_adapt)) tune$mix_adapt <- .05
  # After n0 all the tuning kicks in
  tune$n0 <- 25
  return(tune)
}

check_sampling_settings <- function(pm_settings, stage, n_pars, particles){
  for(i in 1:length(pm_settings)){
    # Mix settings
    pm_settings[[i]]$mix <- check_mix(pm_settings[[i]]$mix, stage)
    # For p_accept
    pm_settings[[i]]$epsilon <- check_epsilon(pm_settings[[i]]$epsilon, n_pars, pm_settings[[i]]$mix)
    # Only local components have an adapted epsilon (see new_particle).
    if(!legacy_sampler()) pm_settings[[i]]$epsilon[!local_components(stage)[-1]] <- 1
    # For mix and p_accept tuning
    pm_settings[[i]]$proposal_counts <- check_prop_performance(pm_settings[[i]]$proposal_counts, stage)
    pm_settings[[i]]$acc_counts <- check_prop_performance(pm_settings[[i]]$acc_counts, stage)
    # Setting particles
    if(is.null(pm_settings[[i]]$n_particles)) pm_settings[[i]]$n_particles <- particles
    if(is.null(pm_settings[[i]]$iter)) pm_settings[[i]]$iter <- 1
    if(is.null(pm_settings[[i]]$gd_good)) pm_settings[[i]]$gd_good <- FALSE
  }
  return(pm_settings)
}

# Retry numerically singular group draws before applying on_exhausted.
resolve_on_singular <- function(on_singular) {
  defaults <- list(max_retries = 3, on_exhausted = "error",
                   max_carry_forward = 10)
  if (is.null(on_singular)) return(defaults)
  if (!is.list(on_singular)) stop("`on_singular` must be a list or NULL")
  bad <- setdiff(names(on_singular), names(defaults))
  if (length(bad)) stop("unknown `on_singular` field(s): ", paste(bad, collapse = ", "),
                        ". Valid fields: ", paste(names(defaults), collapse = ", "))
  out <- modifyList(defaults, on_singular)
  out$on_exhausted <- match.arg(out$on_exhausted, c("error", "carry_forward"))
  out
}

# Is this error the recoverable group-covariance failure? Covers solve()
# ("computationally singular") and chol() non-positive-definite failures. Note R
# >= 4.x reports the latter as "the leading minor of order N is not positive"
# (without "definite"), so match "not positive" and "leading minor" too.
is_singular_error <- function(e) {
  grepl("singular|reciprocal condition|not positive|leading minor|positive definite|Lapack|infinite or missing",
        conditionMessage(e), ignore.case = TRUE)
}

# Parameters with the largest group variance at iteration j, for give-up messages.
top_diverging_pars <- function(store, j, n = 6) {
  v <- tryCatch(diag(store$theta_var[, , j]), error = function(e) NULL)
  if (is.null(v) || all(!is.finite(v))) return("unavailable")
  nm <- rownames(store$theta_mu)
  paste0(nm[order(v, decreasing = TRUE)][seq_len(min(n, length(nm)))], collapse = ", ")
}

# One group (Gibbs) step with retry on a singular covariance. Returns the group
# parameters, or NULL if all retries were exhausted (caller then decides whether
# to carry forward or give up).
robust_gibbs_step <- function(sampler, alpha, type, on_singular) {
  attempt <- 0L
  repeat {
    res <- tryCatch(
      gibbs_step(sampler, alpha, type),
      error = function(e) if (is_singular_error(e))
        structure(list(), class = "gibbs_singular") else stop(e))
    if (!inherits(res, "gibbs_singular")) return(res)
    if (attempt >= on_singular$max_retries) return(NULL)
    attempt <- attempt + 1L
  }
}

run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 100,
                      n_cores = 1,
                      tune = NULL,
                      verbose = TRUE,
                      verboseProgress = TRUE,
                      on_singular = NULL,
                      r_cores = 1) {
  on_singular <- resolve_on_singular(on_singular)
  # Set defaults for NULL values
  # Set necessary local variables
  # Set stable (fixed) new_sample argument for this run
  n_pars <- pmwgs$n_pars
  tune$components <- attr(pmwgs$data, "components")
  tune$shared_ll_idx <- attr(pmwgs$data, "shared_ll_idx")

  pm_settings <- attr(pmwgs$samples, "pm_settings")
  # Intialize sampling tuning settings
  if(is.null(pm_settings)) pm_settings <- lapply(1:pmwgs$n_subjects, function(x) return(vector("list", length(unique(tune$components)))))
  # The tail of adapt tunes the sample kernel before it is frozen.
  kernel <- if(is.null(tune$kernel)) stage else tune$kernel
  tune <- check_tune_settings(tune, n_pars, kernel, particles)
  pm_settings <- lapply(pm_settings, FUN = check_sampling_settings,  stage = kernel, n_pars = n_pars, particles)
  tune$frozen <- stage == "sample" && !legacy_sampler()
  tune$lik_prec <- pmwgs$lik_prec
  tune$exact <- kernel %in% c("adapt", "sample") && !legacy_sampler()
  # Keep a local move on at least a quarter of hierarchical sample iterations.
  tune$min_local <- if(tune$exact && kernel == "sample" && pmwgs$type != "single") .25 else 0
  do_scale <- pmwgs$type == "standard" && kernel %in% scale_move_stages()
  scale_settings <- attr(pmwgs$samples, "scale_move")
  if(do_scale) scale_settings <- scale_move_init(scale_settings, pmwgs$par_names[!pmwgs$nuisance])

  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)

  # Add proposal distributions
  eff_mu <- pmwgs$eff_mu
  eff_var <- pmwgs$eff_var
  chains_var <- pmwgs$chains_var
  chains_mu <- pmwgs$chains_mu
  # Make sure that there's at least something to mapply over
  if(is.null(eff_mu)) eff_mu <- vector("list", pmwgs$n_subjects)
  if(is.null(chains_mu)) chains_mu <- vector("list", pmwgs$n_subjects)
  if(is.null(eff_var)) eff_var <- vector("list", pmwgs$n_subjects)
  if(is.null(chains_var)) chains_var <- vector("list", pmwgs$n_subjects)
  if (verboseProgress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx

  data <- pmwgs$data
  subjects <- pmwgs$subjects
  nuisance <- pmwgs$nuisance
  if(any(nuisance)){
    type <- pmwgs$sampler_nuis$type
    pmwgs$sampler_nuis$samples$idx <- pmwgs$samples$idx
  }
  # Group-covariance recovery bookkeeping (only active when on_singular is set)
  last_good_pars <- NULL; consec_cf <- 0L; n_cf <- 0L
  i <- 0L; j <- start_iter
  # Any error in the loop is enriched with where it happened and what the
  # on_singular recovery had done so far, so a chain failure is actually
  # diagnostic instead of a bare R message.
  tryCatch(
  # Main iteration loop
  for (i in 1:iter) {
    if (verboseProgress) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }
    j <- start_iter + i

    # Gibbs step (with optional recovery from a singular group covariance)
    pars <- robust_gibbs_step(pmwgs, pmwgs$samples$alpha[!nuisance,,j-1, drop = FALSE],
                              pmwgs$type, on_singular)
    if (is.null(pars)) {
      # Retries exhausted: either carry the previous group parameters forward, or
      # give up with a diagnosis of the diverging parameters.
      if (on_singular$on_exhausted == "carry_forward" && !is.null(last_good_pars)) {
        consec_cf <- consec_cf + 1L; n_cf <- n_cf + 1L
        # Terse message; the loop's error handler adds stage/iteration context and
        # the diverging parameters, and check_chain_failures adds the advice.
        if (consec_cf > on_singular$max_carry_forward) {
          stop("group covariance singular for ", consec_cf,
               " consecutive iterations; giving up (on_singular carry_forward)", call. = FALSE)
        }
        pars <- last_good_pars
        # gibbs_step returns alpha as a 2-D (p x n) matrix, so match that shape
        # here (the raw sample slice is 3-D and breaks the downstream indexing).
        a <- pmwgs$samples$alpha[!nuisance, , j - 1, drop = FALSE]
        pars$alpha <- matrix(a, nrow = dim(a)[1], ncol = dim(a)[2],
                             dimnames = list(pmwgs$par_names[!nuisance], NULL))
      } else {
        stop("group covariance became computationally singular", call. = FALSE)
      }
    } else {
      consec_cf <- 0L
      last_good_pars <- pars
    }
    alpha_full <- matrix(pmwgs$samples$alpha[, , j-1], nrow = n_pars, ncol = pmwgs$n_subjects,
                         dimnames = dimnames(pmwgs$samples$alpha)[1:2])
    prev_ll <- pmwgs$samples$subj_ll[, j-1]
    if(do_scale){
      sm <- scale_move_standard(pmwgs, pars, alpha_full, prev_ll, scale_settings, tune$lik_prec,
                                frozen = isTRUE(tune$frozen), n_cores = n_cores, r_cores = r_cores)
      pars <- sm$pars; alpha_full <- sm$alpha; prev_ll <- sm$ll; scale_settings <- sm$settings
    }
    pars_comb <- pars
    if(any(nuisance)){
      pars_nuis <- gibbs_step(pmwgs$sampler_nuis, pmwgs$samples$alpha[nuisance,,j-1, drop = FALSE], pmwgs$sampler_nuis$type)
      pars_comb <- merge_group_level(pars$tmu, pars_nuis$tmu, pars$tvar, pars_nuis$tvar, nuisance, pars$subj_mu)
      pars_comb$alpha <- alpha_full
      pmwgs$sampler_nuis$samples <- fill_samples(samples = pmwgs$sampler_nuis$samples,
                                                                        group_level = pars_nuis,
                                                                        j = j,
                                                                        proposals = NULL,
                                                                        n_pars = n_pars, type = pmwgs$sampler_nuis$type)
      pmwgs$sampler_nuis$samples$idx <- j
    }
    # Particle step
    proposals <- auto_mclapply(seq_len(pmwgs$n_subjects), function(s){
      new_particle(s, data[[s]], pm_settings[[s]], eff_mu[[s]], eff_var[[s]],
                   chains_mu[[s]], chains_var[[s]], prev_ll[s], pars_comb,
                   pmwgs$model, kernel, pmwgs$type, tune, r_cores)
    }, mc.cores = n_cores, rng_substream = TRUE)
    proposals <- do.call(cbind, proposals)
    pm_settings <- proposals[3,]
    proposals <- array(unlist(proposals[1:2,]), dim = c(pmwgs$n_pars + 1, pmwgs$n_subjects))

    #Fill samples
    pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = pars,
                                               proposals = proposals, j = j, n_pars = pmwgs$n_pars, type = pmwgs$type)
  }
  ,
  error = function(e) {
    recov <- if (n_cf > 0)
      sprintf("; on_singular recovery so far: %d carried forward (%d consecutive)",
              n_cf, consec_cf) else
      "; no on_singular recovery was active (see ?fit `on_singular`)"
    stop(conditionMessage(e),
         sprintf("\n  [context: '%s' stage, iteration %d of %d%s]", stage, i, iter, recov),
         "\n  parameters with the largest group variance:\n    ",
         top_diverging_pars(pmwgs$samples, max(start_iter, j - 1)),
         call. = FALSE)
  })
  attr(pmwgs$samples, "pm_settings") <- pm_settings
  if(do_scale) attr(pmwgs$samples, "scale_move") <- scale_settings
  if (verboseProgress) close(pb)
  if (verbose && n_cf > 0) {
    message("  [on_singular] group covariance recovered on ", n_cf, " carried forward",
            " of ", iter, " '", stage, "' iterations")
  }
  return(pmwgs)
}


# options(emc.sampler = "legacy") restores the previous particle step and tuning.
legacy_sampler <- function() identical(getOption("emc.sampler"), "legacy")

# Components centred on the current value, in new_particle's proposal order.
local_components <- function(stage){
  switch(stage,
         preburn = c(FALSE, TRUE),
         burn = c(FALSE, TRUE, TRUE),
         adapt = c(FALSE, TRUE, FALSE),
         c(FALSE, TRUE, FALSE, FALSE))
}

# Row-wise log(sum(exp(.))) of a matrix
log_sum_exp_rows <- function(x){
  m <- apply(x, 1, max)
  m + log(rowSums(exp(x - m)))
}

new_particle <- function (s, data, pm_settings, eff_mu = NULL,
                          eff_var = NULL, chains_mu = NULL,
                          chains_var = NULL, prev_ll,
                          parameters, model = NULL, stage,
                          type, tune, r_cores = 1)
{
  group_pars <- get_group_level(parameters, s, type)
  unq_components <- unique(tune$components)
  group_mu <- group_pars$mu
  group_var <- group_pars$var
  subj_mu <- parameters$alpha[,s]
  out_lls <- numeric(length(unique(tune$shared_ll_idx)))
  particle_multiplier <- 1
  if(stage == "preburn"){
    Mus <- list(group_mu, subj_mu)
    Sigmas <- list(group_var, group_var)
    # For preburn use a lot of proposals, to increase initial search a bit
    particle_multiplier <- 2
  } else if(stage == "burn"){ # Burn
    Mus <- list(group_mu, subj_mu, subj_mu)
    Sigmas <- list(group_var, group_var, chains_var)
  } else if(stage == "adapt"){
    Mus <- list(group_mu, subj_mu, chains_mu)
    Sigmas <- list(group_var, chains_var, chains_var)
  } else{ # Sample
    Mus <- list(group_mu, subj_mu, chains_mu, eff_mu)
    Sigmas <- list(group_var, chains_var, chains_var, eff_var)
  }
  n_proposals <- length(Mus)
  local <- local_components(stage)
  # Exact steps choose either local or independent proposals, never both.
  exact <- isTRUE(tune$exact)
  # Hierarchical proposals follow the current group precision.
  lik <- if(exact) tune$lik_prec[[s]] else NULL
  if(!is.null(lik$post)) lik <- lik$post
  lik_prec <- lik$prec
  prior_prec <- if(is.null(lik_prec)) NULL else group_precision(group_var)
  cond <- conditional_proposal(lik, prior_prec, group_mu)
  if(!is.null(cond)){
    Mus[[3]] <- cond$mu
    Sigmas[[3]] <- cond$var
  }

  for(i in unq_components){
    # Add 1 to epsilons such that prior/group-level proposals aren't scaled
    epsilons <- c(1, pm_settings[[i]]$epsilon)
    idx <- tune$components == i
    mix <- pm_settings[[i]]$mix
    if(exact){
      use_local <- runif(1) < sum(mix[local])
      active <- which(if(use_local) local else !local)
    } else{
      active <- seq_len(n_proposals)
    }
    mix_active <- mix[active] / sum(mix[active])
    particle_numbers <- numeric(n_proposals)
    particle_numbers[active] <- rmultinom(1, pm_settings[[i]]$n_particles*particle_multiplier, mix_active)
    if(!exact) particle_numbers[active] <- pmax(1, particle_numbers[active])
    # Auxiliary centres c_j ~ N(current, S_j) give the ensemble target
    # pi(theta) prod_j N(theta | c_j, S_j), whose marginal is pi.
    centres <- Mus
    covs <- vector("list", n_proposals)
    for(j in active){
      base_cov <- if(exact && local[j]) local_cov(lik_prec, prior_prec, idx) else NULL
      if(is.null(base_cov)) base_cov <- Sigmas[[j]][idx, idx, drop = FALSE]
      covs[[j]] <- base_cov * (epsilons[j]^2)
      if(exact && local[j]) centres[[j]][idx] <- particle_draws(1, Mus[[j]][idx], covs[[j]])
    }
    proposals <- vector("list", n_proposals +1)
    proposals[[1]] <- subj_mu[idx]
    for(j in active){
      proposals[[j + 1]] <- particle_draws(particle_numbers[j], centres[[j]][idx], covs[[j]])
    }
    proposals <- do.call(rbind, proposals)

    # Rejoin new proposals with current MCMC values for other components
    if(any(!idx)){
      proposals_other <- do.call(rbind, rep(list(subj_mu[!idx]), nrow(proposals)))
      colnames(proposals_other) <- names(subj_mu)[!idx]
      colnames(proposals) <- names(subj_mu)[idx]
      proposals <- cbind(proposals, proposals_other)
      proposals <- proposals[, names(subj_mu), drop = FALSE]
    } else{
      colnames(proposals) <- names(subj_mu)
    }

    # Multiple parameter blocks may share one model likelihood.
    shared_idx <- tune$shared_ll_idx[idx][1]
    is_shared <- shared_idx == tune$shared_ll_idx

    # Calculate likelihoods
    if(tune$components[length(tune$components)] > 1){
      lw <- calc_ll_manager(proposals[, is_shared, drop = FALSE], dadm = data, model,
                            component = shared_idx, r_cores = r_cores)
    } else{
      lw <- calc_ll_manager(proposals[, is_shared, drop = FALSE], dadm = data, model,
                            r_cores = r_cores)
    }
    lw_total <- lw + prev_ll - lw[1] # make sure lls from other components are included
    # Prior density
    lp <- mvtnorm::dmvnorm(
      x = proposals[, idx, drop = FALSE],
      mean = group_mu[idx],
      sigma = group_var[idx, idx, drop = FALSE],
      log = TRUE
    )
    if(length(unq_components) > 1){
      prior_density <- mvtnorm::dmvnorm(x = proposals, mean = group_mu, sigma = group_var, log = TRUE)
    } else{
      prior_density <- lp
    }
    # The exact local kernel has one component: its auxiliary target factor
    # cancels its proposal density. Independent/search mixtures need pi / q.
    l <- lw_total + prior_density
    if(!exact || !use_local){
      log_dens <- matrix(NA_real_, nrow(proposals), length(active))
      for(k in seq_along(active)){
        j <- active[k]
        dens <- if(j == 1) lp else mvtnorm::dmvnorm(
          x = proposals[, idx, drop = FALSE],
          mean = centres[[j]][idx],
          sigma = covs[[j]],
          log = TRUE
        )
        log_dens[, k] <- log(mix_active[k]) + dens
      }
      lm <- log_sum_exp_rows(log_dens)
      infnt_idx <- is.infinite(lm)
      lm[infnt_idx] <- min(lm[!infnt_idx])
      l <- l - lm
    }
    weights <- exp(l - max(l))
    # Do MH step and return everything
    idx_ll <- sample(x = sum(particle_numbers) + 1, size = 1, prob = weights)

    out_lls[shared_idx] <- lw[idx_ll]
    subj_mu[idx] <- proposals[idx_ll,idx]
    if(!isTRUE(tune$frozen)) pm_settings[[i]] <- update_pm_settings(pm_settings[[i]], weights, particle_numbers, tune, sum(idx))
  }
  return(list(proposal = unname(subj_mu), ll = sum(out_lls), pm_settings = pm_settings))
}

# Zero nuisance variances contribute no prior precision.
group_precision <- function(group_var){
  free <- diag(group_var) > 0
  P <- matrix(0, nrow(group_var), ncol(group_var), dimnames = dimnames(group_var))
  inv <- tryCatch(solve(group_var[free, free, drop = FALSE]), error = function(e) NULL)
  if(is.null(inv)) return(NULL)
  P[free, free] <- inv
  P
}

# NULL leaves the caller's chain covariance in place.
local_cov <- function(lik_prec, prior_prec, idx = seq_len(nrow(lik_prec))){
  if(is.null(lik_prec) || is.null(prior_prec)) return(NULL)
  tryCatch({
    S <- chol2inv(chol(lik_prec[idx, idx, drop = FALSE] + prior_prec[idx, idx, drop = FALSE]))
    if(all(is.finite(S))) S else NULL
  }, error = function(e) NULL)
}

# Combine the likelihood quadratic with the current Gaussian group prior.
conditional_proposal <- function(lik, prior_prec, group_mu){
  S <- local_cov(lik$prec, prior_prec)
  if(is.null(S)) return(NULL)
  mu <- drop(S %*% (lik$lin + prior_prec %*% group_mu))
  if(!all(is.finite(mu))) return(NULL)
  names(mu) <- names(group_mu)
  list(mu = mu, var = S)
}

# Central-difference quadratic: -1/2 x' prec x + lin' x. Search for steps
# giving a log-likelihood drop near target; bracketed log-scale bisection
# avoids bouncing between a flat stretch and a model bound. Unsettled
# directions get no precision or gradient, including their cross terms.
lik_precision <- function(centre, h, dadm, model, r_cores = 1, target = 1, max_rounds = 8){
  p <- length(centre)
  ll <- function(X){
    colnames(X) <- names(centre)
    as.vector(calc_ll_manager(X, dadm = dadm, model = model, r_cores = r_cores))
  }
  shift <- function(D) sweep(D, 2, centre, "+")
  f0 <- ll(matrix(centre, nrow = 1))
  if(!is.finite(f0)) return(NULL)
  h[!is.finite(h) | h <= 0] <- .1
  E <- diag(p)
  fp <- fm <- rep(NA_real_, p)
  todo <- rep(TRUE, p)
  lo <-rep(NA_real_, p); hi <- rep(NA_real_, p)
  for(r in seq_len(max_rounds)){
    k <- which(todo)
    D <- E[k, , drop = FALSE] * h[k]
    f <- ll(shift(rbind(D, -D)))
    fp[k] <- f[seq_along(k)]; fm[k] <- f[length(k) + seq_along(k)]
    drop_k <- f0 - (fp[k] + fm[k])/2
    ok <- is.finite(drop_k) & drop_k > target/3 & drop_k < target*3
    todo[k[ok]] <- FALSE
    if(!any(todo) || r == max_rounds) break
    kb <- k[!ok]; bad <- drop_k[!ok]
    small <- is.finite(bad) & bad <= target/3
    lo[kb[small]] <- pmax(lo[kb[small]], h[kb[small]], na.rm = TRUE)
    hi[kb[!small]] <- pmin(hi[kb[!small]], h[kb[!small]], na.rm = TRUE)
    h_new <- h[kb] * ifelse(!is.finite(bad), .25,
                            ifelse(bad <= 0, 4, pmin(10, pmax(.1, sqrt(target/bad)))))
    both <- is.finite(lo[kb]) & is.finite(hi[kb])
    h_new[both] <- sqrt(lo[kb[both]] * hi[kb[both]])
    h[kb] <- h_new
  }
  cliff <-todo | !is.finite(fp) | !is.finite(fm)
  H <- diag((2*f0 - fp - fm)/h^2, p)
  if(p > 1){
    pairs <- utils::combn(p, 2)
    D <- matrix(0, ncol(pairs), p)
    D[cbind(seq_len(ncol(pairs)), pairs[1,])] <- h[pairs[1,]]
    D[cbind(seq_len(ncol(pairs)), pairs[2,])] <- h[pairs[2,]]
    f <- ll(shift(rbind(D, -D)))
    fpp <- f[seq_len(ncol(pairs))]; fmm <- f[ncol(pairs) + seq_len(ncol(pairs))]
    off <- -(fpp + fmm - fp[pairs[1,]] - fm[pairs[1,]] - fp[pairs[2,]] - fm[pairs[2,]] + 2*f0) /
      (2*h[pairs[1,]]*h[pairs[2,]])
    H[t(pairs)] <- off
    H[t(pairs[2:1, , drop = FALSE])] <- off
  }
  H[!is.finite(H)] <- 0
  H[cliff, ] <- 0; H[, cliff] <- 0
  # Nearest positive semi-definite matrix: directions in which the likelihood
  # is flat or convex at the centre get no likelihood precision
  eig <- eigen(H, symmetric = TRUE)
  H <- eig$vectors %*% (pmax(eig$values, 0) * t(eig$vectors))
  H <- (H + t(H))/2
  dimnames(H) <- list(names(centre), names(centre))
  grad <- (fp - fm)/(2*h)
  grad[!is.finite(grad) | cliff] <- 0
  list(prec = H, lin = drop(H %*% centre) + grad)
}


update_pm_settings <- function(pm_settings, weights, particle_numbers,
                               tune, n_pars) {
  pm_settings$iter <- pm_settings$iter + 1
  if(pm_settings$iter <= tune$n0) return(pm_settings)

  pm_settings$proposal_counts <- pm_settings$proposal_counts + particle_numbers
  # Fraction of each proposal's particles that beat the current particle.
  offset <- 2
  rate_now <- rep(NA_real_, length(particle_numbers))
  for (j in which(particle_numbers > 0)) {
    n_j <- particle_numbers[j]
    better_j <- sum(weights[offset:(offset + n_j - 1)] > weights[1])
    pm_settings$acc_counts[j] <- pm_settings$acc_counts[j] + better_j
    rate_now[j] <- better_j / n_j
    offset <- offset + n_j
  }
  acc_rates <- ifelse(pm_settings$proposal_counts > 0, pm_settings$acc_counts / pm_settings$proposal_counts, 0)

  legacy <- legacy_sampler()
  clamp <- if(legacy) c(ifelse(length(pm_settings$mix) == 2, .1, ifelse(length(pm_settings$mix) == 3, .4, .6)), 5)
           else c(.01, 20)
  if(isTRUE(tune$exact)){
    # Tune only local proposals, from this iteration's rate. The adapt tail
    # averages later log steps before freezing the sample kernel.
    new_epsilon <- pm_settings$epsilon
    for(j in which(tune$local & !is.na(rate_now))){
      signal <- max(-1, min(3, (rate_now[j] - tune$p_accept[j - 1]) / tune$p_accept[j - 1]))
      new_epsilon[j - 1] <- min(clamp[2], max(clamp[1], exp(log(new_epsilon[j - 1]) + .1 * signal)))
      pm_settings$local_uses <- sum(pm_settings$local_uses, 1)
      if(pm_settings$local_uses > 20) pm_settings$log_eps_sum <- sum(pm_settings$log_eps_sum, log(new_epsilon[j - 1]))
    }
  } else {
    new_epsilon <- update_epsilon_continuous(
      epsilon   = pm_settings$epsilon,
      acceptance = acc_rates[-1],
      target     = tune$p_accept,
      iter       = pm_settings$iter,
      d          = n_pars,
      alphaStar  = tune$alphaStar,
      clamp      = clamp,
      relative   = !legacy
    )
  }
  pm_settings$epsilon <- new_epsilon

  if(length(pm_settings$mix) > 2){
    performance <- (acc_rates + 1e-12) / pm_settings$mix
    performance <- performance / c(mean(tune$p_accept), tune$p_accept)
    performance <- performance / sum(performance)
    new_mix <- (1 - tune$mix_adapt) * pm_settings$mix + tune$mix_adapt * performance
    new_mix <- pmax(new_mix, 0.02)
    new_mix <- new_mix / sum(new_mix)
    if(isTRUE(tune$min_local > 0) && sum(new_mix[tune$local]) < tune$min_local){
      new_mix[tune$local] <- new_mix[tune$local] * tune$min_local / sum(new_mix[tune$local])
      new_mix[!tune$local] <- new_mix[!tune$local] * (1 - tune$min_local) / sum(new_mix[!tune$local])
    }
    pm_settings$mix <- new_mix
  }
  # Exact sample-kernel tuning happens only in the adapt tail.
  if (length(pm_settings$mix) > 3 && (pm_settings$gd_good || isTRUE(tune$exact))) {
    ess <- sum(weights)^2 / sum(weights^2)
    scale_factor <- (tune$target_ESS / ess)^tune$ESS_scale
    new_num_particles <- round(pm_settings$n_particles * scale_factor)
    pm_settings$n_particles <- max(25, min(tune$max_particles, new_num_particles))
    if(isTRUE(tune$exact)){
      # Average ESS per particle for tune_sample_kernel().
      pm_settings$log_ess_sum <- sum(pm_settings$log_ess_sum, log(ess / (length(weights) - 1)))
      pm_settings$log_ess_n <- sum(pm_settings$log_ess_n, 1)
      pm_settings$ess_target <- tune$target_ESS
      pm_settings$max_particles <- tune$max_particles
    }
  }
  return(pm_settings)
}



# Utility functions for sampling below ------------------------------------
update_epsilon_continuous <- function(epsilon, acceptance, target, iter, d,
                                      alphaStar, damp = 100, clamp = c(0.6, 4),
                                      relative = TRUE) {
  c_term <- (1 - 1/d)*sqrt(2*pi)*exp(alphaStar^2/2)/(2*alphaStar) + 1/(d*target*(1-target))
  step_size <- c_term / max(damp, iter)
  # Relative error avoids slow shrinkage when the target acceptance is small.
  diff_accept <- if(relative) pmax(-1, pmin(1, (acceptance - target) / target)) else acceptance - target
  eps_new <- exp(log(epsilon) + step_size * diff_accept)
  pmin(clamp[2], pmax(eps_new, clamp[1]))
}

particle_draws <- function(n, mu, covar, alpha = NULL, tau= NULL) {
  if (n <= 0) {
    return(NULL)
  }
  if (is.null(dim(covar)) && length(covar) == 1L) {
    covar <- matrix(covar, nrow = 1L, ncol = 1L)
  }
  if(is.null(alpha)){
    return(mvtnorm::rmvnorm(n, mu, covar))
  }
}

extend_sampler <- function(sampler, n_samples, stage) {
  # This function takes the sampler and extends it along the intended number of
  # iterations, to ensure that we're not constantly increasing our sampled object
  # by 1. Big shout out to the rapply function
  sampler$samples$stage <- c(sampler$samples$stage, rep(stage, n_samples))
  last_theta_var_inv <- sampler$samples$last_theta_var_inv
  sampler$samples$last_theta_var_inv <- NULL
  if(any(sampler$nuisance)) {
    last_theta_var_inv_nuis <- sampler$sampler_nuis$samples$last_theta_var_inv
    sampler$sampler_nuis$samples$last_theta_var_inv <- NULL
    sampler$sampler_nuis$samples <- rapply(sampler$sampler_nuis$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
    sampler$sampler_nuis$samples$last_theta_var_inv <- last_theta_var_inv_nuis
  }
  sampler$samples <- rapply(sampler$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
  sampler$samples$last_theta_var_inv <- last_theta_var_inv
  return(sampler)
}

extend_obj <- function(obj, n_extend){
  old_dim <- dim(obj)
  n_dimensions <- length(old_dim)
  if(is.null(old_dim) | n_dimensions == 1) return(obj)
  if(n_dimensions == 2){
    if(nrow(obj) == ncol(obj)){
      if(nrow(obj) > 1){
        if(mean(abs(abs(rowSums(obj/max(obj))) - abs(colSums(obj/max(obj))))) < .01) return(obj)
      }
    }
  }
  new_dim <- c(rep(0, (n_dimensions -1)), n_extend)
  extended <- array(NA_real_, dim = old_dim +  new_dim, dimnames = dimnames(obj))
  extended[slice.index(extended,n_dimensions) <= old_dim[n_dimensions]] <- obj
  return(extended)
}

sample_store_base <- function(data, par_names, iters = 1, stage = "init", is_nuisance = rep(F, length(par_names)), ...) {
  subject_ids <- unique(data$subjects)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  samples <- list(
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL))
  )
}

block_variance_idx <- function(components){
  vars_out <- matrix(0, length(components), length(components))
  for(i in unique(components)){
    idx <- i == components
    vars_out[idx,idx] <- NA
  }
  return(vars_out == 0)
}

fill_samples_base <- function(samples, group_level, proposals, j = 1, n_pars){
  # Fill samples both group level and random effects
  samples$theta_mu[, j] <- group_level$tmu
  samples$theta_var[, , j] <- group_level$tvar
  if(!is.null(proposals)) samples <- fill_samples_RE(samples, proposals, j, n_pars)
  return(samples)
}



fill_samples_RE <- function(samples, proposals, j = 1, n_pars, ...){
  # Only for random effects, separated because group level sometimes differs.
  if(!is.null(proposals)){
    samples$alpha[, , j] <- proposals[1:n_pars,]
    samples$subj_ll[, j] <- proposals[n_pars + 1,]
    samples$idx <- j
  }
  return(samples)
}


set_p_accept <- function(stage, search_width){
  # Proposals distributions:
  # 1. Prior - unscaled: all stages
  # 2. Prev particle - scaled chain variance: all stages in preburn scaled by prior variance
  # 3. Chain mean - scaled chain variance: burn onwards
  # 4. Eff mean - scaled eff variance: sample onwards
  # adapt/sample: the local component targets ~3% of particles beating the
  # current value (where the ensemble move's ESS per iteration peaks);
  # independent components are not scaled.
  if(stage == "preburn") return(0.02 * (1/search_width))
  if(stage == "burn") return(c(0.02, 0.25)* (1/search_width))
  if(legacy_sampler()){
    if(stage == "adapt") return(c(0.2, 0.25)* (1/search_width))
    if(stage == "sample") return(c(0.3, 0.3, 0.3)* (1/search_width))
  }
  if(stage == "adapt") return(c(0.03, 0.25)* (1/search_width))
  if(stage == "sample") return(c(0.03, 0.3, 0.3)* (1/search_width))
}

get_default_mix <- function(stage){
  if (stage == "burn") {
    default_mix <- c(0.15, 0.35, 0.5)
  } else if(stage == "adapt"){
    default_mix <- c(0.1, 0.4, 0.4)
  } else if(stage == "sample"){
    default_mix <- c(0.05, 0.3, 0.3, 0.35)
  }  else{
    default_mix <- c(0.5, 0.5)
  }
  return(default_mix)
}


check_mix <- function(mix = NULL, stage) {
  default_mix <- get_default_mix(stage)
  if(is.null(mix)) mix <- default_mix
  if(length(mix) < length(default_mix)){
    if(length(default_mix) - length(mix) == 1){
      mix <- mix*(1-default_mix[length(default_mix)])
      mix <- c(mix, default_mix[length(default_mix)])
    } else{
      stop("mix settings are not compatible with stage")
    }
  }
  return(mix)
}

check_epsilon <- function(epsilon, n_pars, mix) {
  if (is.null(epsilon)) { # In preburn case, there's only one epsilon here
    if (n_pars > 15) {
      epsilon <- .5
    } else if (n_pars > 10) {
      epsilon <- .6
    } else {
      epsilon <- .7
    }
    # The first proposal is always unscaled
    # Every subsequent phase at most adds one epsilon
  } else if(length(epsilon) < (length(mix) -1)){
    epsilon <- c(epsilon, epsilon[length(epsilon)])
  }
  return(epsilon)
}

check_prop_performance <- function(prop_performance, stage){
  default_mix <- get_default_mix(stage)
  if(is.null(prop_performance) || stage == "adapt") prop_performance <- rep(0, length(default_mix))
  if(length(prop_performance) < length(default_mix)){
    if(length(default_mix) - length(prop_performance) == 1){
      n_total <- sum(prop_performance)
      prop_performance <- prop_performance*(1-default_mix[length(default_mix)])
      prop_performance <- c(prop_performance, default_mix[length(default_mix)]*n_total)
    } else{
      stop("prop_performance settings are not compatible with stage")
    }
  }
  return(round(prop_performance))
}

calc_ll_manager <- function(proposals, dadm, model, component = NULL, r_cores = 1, return_trialwise=FALSE){
  if(!is.data.frame(dadm)){
    lls <- log_likelihood_joint(proposals, dadm, model, component)
  } else{
    model <- model()
    if(is.null(model$c_name)){ # use the R implementation
      lls <- apply(proposals,1, calc_ll_R, model, dadm = dadm)
    } else {
      p_types <- names(model$p_types)
      designs <- get_designs_expanded(dadm, model)
      constants <- attr(dadm, "constants")
      if(is.null(constants)) constants <- NA

      backend   <- getOption("emc.ll_backend", default = "multiprocess")
      n_threads <- getOption("emc.n_threads", default = 1)

      if (backend == "multithreaded") {
        lls <- calc_ll_multithreaded(proposals, dadm, constants = constants, designs = designs,
                                     type = model$c_name, model$bound, model$transform,
                                     model$pre_transform, p_types = p_types,
                                     min_ll = log(1e-10), model$trend, n_threads = n_threads, return_trialwise=return_trialwise)
      } else {
        lls <- calc_ll(proposals, dadm, constants = constants, designs = designs,
                       type = model$c_name, model$bound, model$transform,
                       model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                       model$trend, return_trialwise=return_trialwise)
      }
    }
  }
  return(lls)
}

merge_group_level <- function(tmu, tmu_nuis, tvar, tvar_nuis, is_nuisance, subj_mu){
  n_pars <- length(is_nuisance)
  tmu_out <- numeric(n_pars)
  tmu_out[!is_nuisance] <- tmu
  tmu_out[is_nuisance] <- tmu_nuis
  tvar_out <- matrix(0, nrow = n_pars, ncol = n_pars)
  tvar_out[!is_nuisance, !is_nuisance] <- tvar

  subj_mu_out <- matrix(NA, ncol = ncol(subj_mu), nrow = length(tmu_out))
  subj_mu_out[is_nuisance,] <- do.call(cbind, rep(list(c(tmu_nuis)), ncol(subj_mu)))
  subj_mu_out[!is_nuisance,] <- subj_mu
  return(list(tmu = tmu_out, tvar = tvar_out, subj_mu = subj_mu_out))
}


#' Run a Group-level Model.
#'
#' Separate function for running only the group-level model. This can be useful in a
#' two-step analysis. Works similar in functionality to make_emc,
#' except also does the fitting and returns an emc object that works with
#' most posterior checking tests (but not the data generation/posterior predictives).
#'
#' @param prior an emc.prior object.
#' @param iter Number of MCMC samples to collect.
#' @inheritParams make_emc
#'
#' @returns an emc object with only group-level samples
#' @export
run_hyper <- function(type = "standard", data, prior = NULL, iter = 1000, n_chains =3, ...){
  args <- list(...)
  if(length(dim(data)) == 3){
    data_input <- data
    data <- as.data.frame(t(data_input[,,1]))
    data$subjects <- 1:nrow(data)
    iter <- dim(data_input)[3]
    is_mcmc <- T
    pars <- rownames(data_input)
  } else{
    data_input <- data[,colnames(data)!= "subjects"]
    is_mcmc <- F
    pars <- colnames(data_input)
  }
  emc <- list()
  for(j in 1:n_chains){
    samples <- sample_store(data = data ,par_names = pars, is_nuisance = rep(F, length(pars)), integrate = F, type = type, ...)
    subjects <- unique(data$subjects)
    sampler <- list(
      data = split(data, data$subjects),
      par_names = pars,
      subjects = subjects,
      n_pars = length(pars),
      nuisance = rep(F, length(pars)),
      n_subjects = length(subjects),
      samples = samples,
      init = TRUE
    )
    class(sampler) <- "pmwgs"
    sampler <- add_info(sampler, prior, type = type, ...)
    sampler$type <- type
    startpoints <- get_startpoints(sampler, start_mu = NULL, start_var = NULL, type = type)
    sampler$samples <- fill_samples(samples = sampler$samples, group_level = startpoints, proposals = NULL,
                                    j = 1, n_pars = sampler$n_pars, type = type)
    sampler$samples$idx <- 1
    sampler <- extend_sampler(sampler, iter-1, "sample")
    for(i in 2:iter){
      if(is_mcmc){
        group_pars <- gibbs_step(sampler, data_input[,,i], type = type)
      } else{
        group_pars <- gibbs_step(sampler, t(data_input), type = type)
      }
      sampler$samples$idx <- i
      sampler$samples <- fill_samples(samples = sampler$samples, group_level = group_pars, proposals = NULL,
                                                   j = i, n_pars = sampler$n_pars, type = type)
    }
    emc[[j]] <- sampler
  }
  emc[[1]]$type <- type
  class(emc) <- "emc"
  emc <- subset(emc, filter = 1)
  return(emc)
}

check_CR <- function(emc, p_vector, range = .2, N = 500){
  covs <- diag(length(p_vector)) * range
  props <- mvtnorm::rmvnorm(N, mean = p_vector, sigma = covs)
  model <- emc[[1]]$model
  if(is.null(model()$c_name)) stop("C not implemented yet for this model")
  dat <- emc[[1]]$data[[1]]
  modelRlist <- model()
  modelRlist$c_name <- NULL
  modelR <- function()return(modelRlist)
  t1 <- system.time(
    R <- calc_ll_manager(props, dat, modelR)
  )
  t2 <- system.time(
    C <- calc_ll_manager(props, dat, model)
  )
  print(paste0("C ", t1$elapsed/t2$elapsed, " times faster"))
  if(!identical(C, R)){
    warning("C and R results differ")
  }
  return(list(C = C, R = R))
}
