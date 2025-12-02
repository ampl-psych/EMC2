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
  proposals <- parallel::mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,
                                  parameters = startpoints_comb, n_particles = particles,
                                  pmwgs = pmwgs, type = type,
                                  mc.cores = n_cores, r_cores = r_cores)
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
  emc <- mclapply(emc,init,start_mu = start_mu, start_var = start_var,
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

run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 100,
                      n_cores = 1,
                      tune = NULL,
                      verbose = TRUE,
                      verboseProgress = TRUE,
                      r_cores = 1) {
  # Set defaults for NULL values
  # Set necessary local variables
  # Set stable (fixed) new_sample argument for this run
  n_pars <- pmwgs$n_pars
  tune$components <- attr(pmwgs$data, "components")
  tune$shared_ll_idx <- attr(pmwgs$data, "shared_ll_idx")

  pm_settings <- attr(pmwgs$samples, "pm_settings")
  # Intialize sampling tuning settings
  if(is.null(pm_settings)) pm_settings <- lapply(1:pmwgs$n_subjects, function(x) return(vector("list", length(unique(tune$components)))))
  tune <- check_tune_settings(tune, n_pars, stage, particles)
  pm_settings <- lapply(pm_settings, FUN = check_sampling_settings,  stage = stage, n_pars = n_pars, particles)

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
  block_idx <- block_variance_idx(tune$components)
  # Main iteration loop
  for (i in 1:iter) {
    if (verboseProgress) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }
    j <- start_iter + i

    # Gibbs step
    pars <- pars_comb <- gibbs_step(pmwgs, pmwgs$samples$alpha[!nuisance,,j-1], pmwgs$type)
    if(any(nuisance)){
      pars_nuis <- gibbs_step(pmwgs$sampler_nuis, pmwgs$samples$alpha[nuisance,,j-1], pmwgs$sampler_nuis$type)
      pars_comb <- merge_group_level(pars$tmu, pars_nuis$tmu, pars$tvar, pars_nuis$tvar, nuisance, pars$subj_mu)
      pars_comb$alpha <- pmwgs$samples$alpha[,,j-1]
      pmwgs$sampler_nuis$samples <- fill_samples(samples = pmwgs$sampler_nuis$samples,
                                                                        group_level = pars_nuis,
                                                                        j = j,
                                                                        proposals = NULL,
                                                                        n_pars = n_pars, type = pmwgs$sampler_nuis$type)
      pmwgs$sampler_nuis$samples$idx <- j
    }
    # Particle step
    proposals <- parallel::mcmapply(new_particle, 1:pmwgs$n_subjects, data, pm_settings, eff_mu, eff_var,
                                    chains_mu, chains_var, pmwgs$samples$subj_ll[,j-1],
                                    MoreArgs = list(pars_comb, pmwgs$model, stage,
                                                    pmwgs$type,
                                                    tune),
                                    mc.cores =n_cores, r_cores = r_cores)
    pm_settings <- proposals[3,]
    proposals <- array(unlist(proposals[1:2,]), dim = c(pmwgs$n_pars + 1, pmwgs$n_subjects))

    #Fill samples
    pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = pars,
                                               proposals = proposals, j = j, n_pars = pmwgs$n_pars, type = pmwgs$type)
  }
  attr(pmwgs$samples, "pm_settings") <- pm_settings
  if (verboseProgress) close(pb)
  return(pmwgs)
}


new_particle <- function (s, data, pm_settings, eff_mu = NULL,
                          eff_var = NULL, chains_mu = NULL,
                          chains_var = NULL, prev_ll,
                          parameters, model = NULL, stage,
                          type, tune, r_cores = 1)
{
  group_pars <- get_group_level(parameters, s, type)
  unq_components <- unique(tune$components)
  proposal_out <- numeric(length(group_pars$mu))
  group_mu <- group_pars$mu
  group_var <- group_pars$var
  subj_mu <- parameters$alpha[,s]
  out_lls <- numeric(length(unq_components))
  particle_multiplier <- 1
  # Set the proposals
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
  for(i in unq_components){
    # Add 1 to epsilons such that prior/group-level proposals aren't scaled
    epsilons <- c(1, pm_settings[[i]]$epsilon)
    idx <- tune$components == i
    # Draw new proposals for each component
    particle_numbers <- numbers_from_proportion(pm_settings[[i]]$mix, pm_settings[[i]]$n_particles*particle_multiplier)
    proposals <- vector("list", n_proposals +1)
    proposals[[1]] <- subj_mu[idx]
    for(j in 1:n_proposals){
      # Fill up the proposals
      proposals[[j + 1]] <- particle_draws(particle_numbers[j], Mus[[j]][idx], Sigmas[[j]][idx,idx] * (epsilons[j]^2))
    }
    proposals <- do.call(rbind, proposals)

    # Non -used proposals (for prior calculations)
    # Rejoin new proposals with current MCMC values for other components
    if(any(!idx)){
      proposals_other <- do.call(rbind, rep(list(subj_mu[!idx]), nrow(proposals)))
      colnames(proposals_other) <- names(subj_mu)[!idx]
      colnames(proposals) <- names(subj_mu)[idx]
      proposals <- cbind(proposals, proposals_other)
      proposals <- proposals[,names(subj_mu)]
    } else{
      colnames(proposals) <- names(subj_mu)
    }

    # Normally we assume that a component contains all the parameters to estimate the individual likelihood of a joint model
    # Sometimes we may also want to block within a model if it has very high dimensionality
    shared_idx <- tune$shared_ll_idx[idx][1]
    is_shared <- shared_idx == tune$shared_ll_idx

    # Calculate likelihoods
    if(tune$components[length(tune$components)] > 1){
      lw <- calc_ll_manager(proposals[,is_shared], dadm = data, model,
                            component = shared_idx, r_cores = r_cores)
    } else{
      lw <- calc_ll_manager(proposals[,is_shared], dadm = data, model,
                            r_cores = r_cores)
    }
    lw_total <- lw + prev_ll - lw[1] # make sure lls from other components are included
    # Prior density
    lp <- mvtnorm::dmvnorm(x = proposals[,idx], mean = group_mu[idx], sigma = group_var[idx,idx], log = TRUE)
    if(length(unq_components) > 1){
      prior_density <- mvtnorm::dmvnorm(x = proposals, mean = group_mu, sigma = group_var, log = TRUE)
    } else{
      prior_density <- lp
    }
    # We can start from 2, since first proposal is prior density
    lm <- pm_settings[[i]]$mix[1]*exp(lp)
    for(k in 2:length(Sigmas)){
      # Prior density is updated separately so start at 2
      lm <- lm + pm_settings[[i]]$mix[k]*mvtnorm::dmvnorm(x = proposals[,idx], mean = Mus[[k]][idx], sigma = Sigmas[[k]][idx,idx]*(epsilons[k]^2))
    }
    # Avoid infinite values
    lm <- log(lm)
    infnt_idx <- is.infinite(lm)
    lm[infnt_idx] <- min(lm[!infnt_idx])
    # Calculate weights and center
    l <- lw_total + prior_density - lm
    weights <- exp(l - max(l))
    # Do MH step and return everything
    idx_ll <- sample(x = sum(particle_numbers) + 1, size = 1, prob = weights)

    out_lls[i] <- lw[idx_ll]
    proposal_out[idx] <- proposals[idx_ll,idx]
    pm_settings[[i]] <- update_pm_settings(pm_settings[[i]], idx_ll, weights, particle_numbers, tune, sum(idx))
  }
  return(list(proposal = proposal_out, ll = sum(out_lls), pm_settings = pm_settings))
}


update_pm_settings <- function(pm_settings, chosen_idx, weights, particle_numbers,
                               tune, n_pars) {
  # 0) If we're past an initial burn-in, do the adaptation
  pm_settings$iter <- pm_settings$iter + 1
  if (pm_settings$iter > tune$n0) {
    # A) Update proposal_counts
    # -------------------------
    # Each proposal j was used "particle_numbers[j]" times
    pm_settings$proposal_counts <- pm_settings$proposal_counts + particle_numbers

    # B) Update acceptance counts: "percent better than old"
    # ------------------------------------------------------
    # weights[1] = old particle's weight
    # weights[2..(1+sum(particle_numbers))] = new draws' weights
    old_weight <- weights[1]

    # We'll parse out each proposal's chunk in weights[-1]
    # using the known counts in particle_numbers.
    # We're also tracking acceptance of group-level proposals, which is minorly wasteful
    offset <- 2  # start index in 'weights' for new proposals
    for (j in seq_along(particle_numbers)) {
      # The chunk of new weights for proposal j
      n_j <- particle_numbers[j]
      if (n_j > 0) {
        draws_j <- weights[offset:(offset + n_j - 1)]
        # Count how many draws_j exceed old_weight
        better_j <- sum(draws_j > old_weight)
        # Accumulate that in acceptance counts
        pm_settings$acc_counts[j] <- pm_settings$acc_counts[j] + better_j
        offset <- offset + n_j
      }
    }

    # C) Compute per-proposal acceptance rates
    # ----------------------------------------
    acc_rates <- ifelse(pm_settings$proposal_counts > 0, pm_settings$acc_counts / pm_settings$proposal_counts, 0)

    # .1 for preburn, .4 for burn and adapt and .6 for sample
    clamp_min <- ifelse(length(pm_settings$mix) == 2, .1, ifelse(length(pm_settings$mix) == 3, .4, .6))

    # D) Adapt epsilon via continuous approach
    # ----------------------------------------
    # pm_settings$epsilon is a vector, same length as pm_settings$mix
    # tune$p_accept is also a vector, e.g. c(0.2, 0.3, 0.6) for each proposal
    new_epsilon <- update_epsilon_continuous(
      epsilon   = pm_settings$epsilon,
      acceptance = acc_rates[-1],
      target     = tune$p_accept,
      iter       = pm_settings$iter,
      d          = n_pars,
      alphaStar  = tune$alphaStar,
      damp       = 100,          # Example
      clamp      = c(clamp_min, 5)    # Example range
    )
    pm_settings$epsilon <- new_epsilon

    # E) Adapt mixing weights based on acceptance vs. target
    # ------------------------------------------------------
    # If ratio_j > 1 (proposal j acceptance > target), its mix goes up;
    # if ratio_j < 1, mix goes down.

    if(length(pm_settings$mix) > 2){ # We're not in preburn
      eps_val <- 1e-12   # Avoid divide-by-zero

      # 1) Compute performance ~ (acceptance / old_mix), normalized
      performance <- (acc_rates + eps_val) / pm_settings$mix
      performance <- performance / sum(performance)

      # 2) Adjust the elements by expected acceptance (p_accept)
      # The first element is the group-level, so has no p_accept, just take the mean of the others
      # Bit hacky
      adj_factor <- c(mean(tune$p_accept), tune$p_accept)
      performance <- performance / adj_factor

      # 3) Normalize again
      performance <- performance / sum(performance)

      # 4) Blend with old mix: stable update
      new_mix <- (1 - tune$mix_adapt) * pm_settings$mix + tune$mix_adapt * performance

      # 5) Impose a floor, re-normalize
      new_mix <- pmax(new_mix, 0.02)
      new_mix <- new_mix / sum(new_mix)
      pm_settings$mix <- new_mix
    }


    # F) Adapt the number of particles (ESS logic)
    # -------------------------------------------------------
    # If length mix > 2, we're in sample stage
    # Only reduce number of particles when we're already converged
    if (length(pm_settings$mix) > 3 && pm_settings$gd_good) {
      ess <- sum(weights)^2 / sum(weights^2)
      desired_ess <- tune$target_ESS
      scale_factor <- (desired_ess / ess)^tune$ESS_scale
      new_num_particles <- round(pm_settings$n_particles * scale_factor)
      pm_settings$n_particles <- max(25, min(tune$max_particles, new_num_particles))
    }
  }

  return(pm_settings)
}



# Utility functions for sampling below ------------------------------------
update_epsilon_continuous <- function(
    epsilon,    # vector of current epsilons
    acceptance, # vector of acceptance rates, same length
    target,     # vector of target acceptance rates, same length
    iter,
    d,
    alphaStar,
    damp = 100,
    clamp = c(0.6, 4)
) {
  log_eps <- log(epsilon)
  # 2) define step size
  # We'll do one pass per element. If you want a single c_term, that's also fine.
  c_term <- (1 - 1/d)*sqrt(2*pi)*exp(alphaStar^2/2)/(2*alphaStar) + 1/(d*target*(1-target))
  step_size <- c_term / max(damp, iter)
  # 3) compute difference from target acceptance
  diff_accept <- acceptance - target
  # 4) update in log space (vectorized)
  log_eps_new <- log_eps + step_size * diff_accept
  # 5) exponentiate
  eps_new <- exp(log_eps_new)
  # 6) clamp
  eps_new <- pmin(clamp[2], pmax(eps_new, clamp[1]))

  return(eps_new)
}

update_epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c <- ((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta <- log(sqrt(epsilon2))
  Theta <- Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}


numbers_from_proportion <- function(mix_proportion, n_particles = 1000) {
  # Make sure each proposal has at least 1 particle
  return(pmax(1, rmultinom(1, n_particles, mix_proportion)))
}


particle_draws <- function(n, mu, covar, alpha = NULL, tau= NULL) {
  if (n <= 0) {
    return(NULL)
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
  if(any(sampler$nuisance)) sampler$sampler_nuis$samples <- rapply(sampler$sampler_nuis$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
  sampler$samples <- rapply(sampler$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
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
  if(stage == "preburn") return(0.02 * (1/search_width))
  if(stage == "burn") return(c(0.02, 0.25)* (1/search_width))
  if(stage == "adapt") return(c(0.2, 0.25)* (1/search_width))
  if(stage == "sample") return(c(0.3, 0.3, 0.3)* (1/search_width))
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

calc_ll_manager <- function(proposals, dadm, model, component = NULL, r_cores = 1){
  if(!is.data.frame(dadm)){
    lls <- log_likelihood_joint(proposals, dadm, model, component)
  } else{
    model <- model()
    if(is.null(model$c_name)){ # use the R implementation
      lls <- unlist(
        auto_mclapply(1:nrow(proposals),
          function(i) calc_ll_R(proposals[i,], model=model, dadm = dadm),
         mc.cores=r_cores))
    } else{
      p_types <- names(model$p_types)
      designs <- list()
      for(p in p_types){
        designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
      }
      constants <- attr(dadm, "constants")
      if(is.null(constants)) constants <- NA
      if (nrow(proposals) <= r_cores)
        lls <- calc_ll(proposals, dadm, constants = constants, designs = designs, type = model$c_name,
                     model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                     model$trend) else {
        idx <- rep(1:r_cores,each=1+(nrow(proposals) %/% r_cores))[1:nrow(proposals)]
        lls <- unlist(auto_mclapply(1:r_cores,function(i) {
          calc_ll(proposals[idx==i,,drop=FALSE], dadm, constants = constants,
            designs = designs, type = model$c_name, model$bound, model$transform,
            model$pre_transform, p_types = p_types, min_ll = log(1e-10),model$trend)
          },mc.cores=r_cores))
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

