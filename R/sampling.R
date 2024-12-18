pmwgs <- function(dadm, type, pars = NULL, prior = NULL,
                  nuisance = NULL, nuisance_non_hyper = NULL, ...) {
  if(is.data.frame(dadm)) dadm <- list(dadm)
  dadm <- extractDadms(dadm)
  if(is.null(pars)) pars <- dadm$pars
  if(is.null(prior)) prior <- dadm$prior
  dadm_list <-dadm$dadm_list
  # Storage for the samples.
  subjects <- sort(as.numeric(unique(dadm$subjects)))
  if(!is.null(nuisance) & !is.numeric(nuisance)) nuisance <- which(pars %in% nuisance)
  if(!is.null(nuisance_non_hyper) & !is.numeric(nuisance_non_hyper)) nuisance_non_hyper <- which(pars %in% nuisance_non_hyper)

  if(!is.null(nuisance_non_hyper)){
    is_nuisance <- is.element(seq_len(length(pars)), nuisance_non_hyper)
    type <- "single"
  } else if(!is.null(nuisance)) {
    is_nuisance <- is.element(seq_len(length(pars)), nuisance)
    type <- "diagonal"
  } else{
    is_nuisance <- rep(F, length(pars))
  }


  sampler_nuis <- NULL
  if(any(is_nuisance)){
    sampler_nuis <- list(
      samples = sample_store(dadm, pars, type, integrate = F,
                                                    is_nuisance = !is_nuisance, ...),
      n_subjects = length(subjects),
      n_pars = sum(is_nuisance),
      nuisance = rep(F, sum(is_nuisance)),
      type = type
    )
    if(type == "single") sampler_nuis$samples <- NULL
    sampler_nuis <- add_info(sampler_nuis, prior$prior_nuis, type, ...)
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
                 verbose = FALSE, particles = 1000, n_cores = 1) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  type <- pmwgs$type
  startpoints <-startpoints_comb <- get_startpoints(pmwgs, start_mu, start_var, type)
  if(any(pmwgs$nuisance)){
    type_nuis <- pmwgs$sampler_nuis$type
    startpoints_nuis <- get_startpoints(pmwgs$sampler_nuis, start_mu = NULL, start_var = NULL, type = type_nuis)
    startpoints_comb <- merge_group_level(startpoints$tmu, startpoints_nuis$tmu,
                                          startpoints$tvar, startpoints_nuis$tvar,
                                          pmwgs$nuisance)
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
                                  mc.cores = n_cores)
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))

  # Sample the mixture variables' initial values.

  pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = startpoints, proposals = proposals,
                                             j = 1, n_pars = pmwgs$n_pars, type = type)
  pmwgs$init <- TRUE
  return(pmwgs)
}

#' Initialize chains
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
#'
#' @return An emc object
#' @examples \donttest{
#' # Make a design and an emc object
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1),
#'                            constants=c(s=log(1)))
#'
#' DDMaE <- make_emc(forstmann, design_DDMaE)
#' # set up our mean starting points (same used across subjects).
#' mu <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'        t0=log(.2))
#' # Small variances to simulate start points from a tight range
#' var <- diag(0.05, length(mu))
#' # Initialize chains, 4 cores per chain, and parallelizing across our 3 chains as well
#' # so 4*3 cores used.
#' DDMaE <- init_chains(DDMaE, start_mu = mu, start_var = var,
#'                      cores_per_chain = 1, cores_for_chains = 1)
#' # Afterwards we can just use fit
#' # DDMaE <- fit(DDMaE, cores_per_chain = 4)
#' }
#' @export
init_chains <- function(emc, start_mu = NULL, start_var = NULL, particles = 1000,
                        cores_per_chain=1,cores_for_chains = length(emc))
{
  emc <- mclapply(emc,init,start_mu = start_mu, start_var = start_var,
           verbose = FALSE, particles = particles,
           n_cores = cores_per_chain, mc.cores=cores_for_chains)
  class(emc) <- "emc"
  return(emc)
}

start_proposals <- function(s, parameters, n_particles, pmwgs, type){
  #Draw the first start point
  group_pars <- get_group_level(parameters, s, type)
  proposals <- particle_draws(n_particles, group_pars$mu, group_pars$var)
  colnames(proposals) <- rownames(pmwgs$samples$alpha) # preserve par names
  lw <- calc_ll_manager(proposals, dadm = pmwgs$data[[which(pmwgs$subjects == s)]],
                        model = pmwgs$model)
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx]))
}


check_tune_settings <- function(tune, n_pars){
  if(is.null(tune$p_accept)) tune$p_accept <- .8
  if(is.null(tune$components)) tune$components <- rep(1, n_pars)
  if(is.null(tune$shared_ll_idx)) tune$shared_ll_idx <- tune$components
  if(is.null(tune$ESS_ratio)) tune$ESS_ratio <- .5
  if(is.null(tune$ESS_scale)) tune$ESS_scale <- .5
  if(is.null(tune$mix_adapt)) tune$mix_adapt <- .05
  tune$alphaStar <- -qnorm(tune$p_accept/2)
  tune$n0 <- 100
  return(tune)
}

check_sampling_settings <- function(pm_settings, stage, n_pars, particles){
  pm_settings$mix <- check_mix(pm_settings$mix, stage)
  pm_settings$epsilon <- check_epsilon(pm_settings$epsilon, n_pars)
  pm_settings$prop_performance <- check_prop_performance(pm_settings$prop_performance, stage)
  if(is.null(pm_settings$n_particles)) pm_settings$n_particles <- particles
  return(pm_settings)
}

run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 100,
                      n_cores = 1,
                      tune = NULL,
                      verbose = TRUE,
                      verboseProgress = TRUE) {
  # Set defaults for NULL values
  # Set necessary local variables
  # Set stable (fixed) new_sample argument for this run
  n_pars <- pmwgs$n_pars
  tune$components <- attr(pmwgs$data, "components")
  tune$shared_ll_idx <- attr(pmwgs$data, "shared_ll_idx")

  pm_settings <- pmwgs$samples$pm_settings
  # Intialize sampling tuning settings
  if(is.null(pm_settings)) pm_settings <- vector("list", pmwgs$n_subjects)
  tune <- check_tune_settings(tune, n_pars)
  pm_settings <- lapply(pm_settings, FUN = check_sampling_settings,  stage = stage, n_pars = n_pars, particles)

  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)

  # Add proposal distributions
  eff_mu <- pmwgs$eff_mu
  eff_var <- pmwgs$eff_var
  chains_var <- pmwgs$chains_var

  # Make sure that there's at least something to mapply over
  if(is.null(eff_mu)) eff_mu <- vector("list", pmwgs$n_subjects)
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
      pars_comb <- merge_group_level(pars$tmu, pars_nuis$tmu, pars$tvar, pars_nuis$tvar, nuisance)
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
                                    chains_var, pmwgs$samples$subj_ll[,j-1],
                                    MoreArgs = list(pars_comb, pmwgs$model, stage,
                                                    pmwgs$type,
                                                    tune, j),
                                    mc.cores =n_cores)
    pm_settings <- proposals[3,]
    proposals <- array(unlist(proposals[1:2,]), dim = c(pmwgs$n_pars + 1, pmwgs$n_subjects))

    #Fill samples
    pmwgs$samples <- fill_samples(samples = pmwgs$samples, group_level = pars,
                                               proposals = proposals, j = j, n_pars = pmwgs$n_pars, type = pmwgs$type)
    pmwgs$samples$pm_settings <- pm_settings
  }
  if (verboseProgress) close(pb)
  return(pmwgs)
}


new_particle <- function (s, data, pm_settings, eff_mu = NULL,
                          eff_var = NULL, chains_var, prev_ll,
                          parameters, model = NULL, stage,
                          type, tune, iter)
{
  group_pars <- get_group_level(parameters, s, type)
  unq_components <- unique(tune$components)
  proposal_out <- numeric(length(group_pars$mu))
  group_mu <- group_pars$mu
  group_var <- group_pars$var
  subj_mu <- parameters$alpha[,s]
  eff_var_old <- eff_var
  if(stage != "sample"){
    eff_mu <- subj_mu
    group_var_subj <- group_var
    if(length(unq_components) > 1){
      group_var_subj[block_idx] <- 0
    }
  }
  out_lls <- numeric(length(unq_components))
  for(i in unq_components){
    if(stage != "sample"){
      eff_var <- chains_var * pm_settings$epsilon^2
      var_subj <- group_var_subj *  pm_settings$epsilon^2
    } else{
      eff_var <- eff_var_old * pm_settings$epsilon^2
      var_subj <- chains_var *  pm_settings$epsilon^2
    }
    idx <- tune$components == i
    # Draw new proposals for this component
    particle_numbers <- numbers_from_proportion(pm_settings$mix, pm_settings$n_particles)
    pop_particles <- particle_draws(particle_numbers[1], group_mu[idx], group_var[idx, idx])
    ind_particles <- particle_draws(particle_numbers[2], subj_mu[idx], var_subj[idx, idx])
    if(length(pm_settings$mix) < 3){
      eff_particles <- NULL
    } else{
      eff_particles <- particle_draws(particle_numbers[3], eff_mu[idx], eff_var[idx, idx ])# eff_alpha, eff_tau)
    }
    # Rejoin new proposals with current MCMC values for other components
    proposals <- matrix(rep(subj_mu, pm_settings$n_particles + 1), nrow = pm_settings$n_particles + 1, byrow = T)
    colnames(proposals) <- names(subj_mu)
    proposals[2:(pm_settings$n_particles+1),idx] <- rbind(pop_particles, ind_particles, eff_particles)
    ll_proposals <- proposals
    # Normally we assume that a component contains all the parameters to estimate the individual likelihood of a joint model
    # Sometimes we may also want to block within a model if it has very high dimensionality
    shared_idx <- tune$shared_ll_idx[idx][1]
    is_shared <- shared_idx == tune$shared_ll_idx

    # Calculate likelihoods
    if(tune$components[length(tune$components)] > 1){
      lw <- calc_ll_manager(ll_proposals[,is_shared], dadm = data, model, component = shared_idx)
    } else{
      lw <- calc_ll_manager(ll_proposals[,is_shared], dadm = data, model)
    }
    lw_total <- lw + prev_ll - lw[1] # make sure lls from other components are included
    lp <- mvtnorm::dmvnorm(x = proposals[,idx], mean = group_mu[idx], sigma = group_var[idx,idx], log = TRUE)
    prop_density <- mvtnorm::dmvnorm(x = proposals[,idx], mean = subj_mu[idx], sigma = var_subj[idx,idx])
    if(length(unq_components) > 1){
      prior_density <- mvtnorm::dmvnorm(x = proposals, mean = group_mu, sigma = group_var, log = TRUE)
    } else{
      prior_density <- lp
    }

    if (length(pm_settings$mix) < 3) {
      eff_density <- 0
      lm <- log(pm_settings$mix[1] * exp(lp) + (pm_settings$mix[2] * prop_density))
    }
    else {
      eff_density <- mvtnorm::dmvnorm(x = proposals[,idx], mean = eff_mu[idx], sigma = eff_var[idx,idx])
      lm <- log(pm_settings$mix[1] * exp(lp) + (pm_settings$mix[2] * prop_density) + (pm_settings$mix[3] * eff_density))
    }

    infnt_idx <- is.infinite(lm)
    lm[infnt_idx] <- min(lm[!infnt_idx])
    # Calculate weights and center
    l <- lw_total + prior_density - lm
    weights <- exp(l - max(l))
    # Do MH step and return everything
    idx_ll <- sample(x = pm_settings$n_particles+1, size = 1, prob = weights)

    out_lls[i] <- lw[idx_ll]
    proposal_out[idx] <- proposals[idx_ll,idx]
    pm_settings <- update_pm_settings(pm_settings, idx_ll, weights, particle_numbers, tune, iter, length(subj_mu))
  } # Note that origin only contains last components origin, just used by devs anyway
  return(list(proposal = proposal_out, ll = sum(out_lls), pm_settings = pm_settings))
}


update_pm_settings <- function(pm_settings, chosen_idx, weights, particle_numbers, tune, iter, n_pars) {
  prop_performance <- pm_settings$prop_performance
  mix <- pm_settings$mix
  # Check if chosen_idx corresponds to a new proposal (not the previous particle)
  if (chosen_idx > 1) {
    # Identify which proposal block produced the chosen particle
    origin <- min(which(chosen_idx <=  (cumsum(particle_numbers) +1)))
    prop_performance[origin] <- prop_performance[origin] + 1
  }
  # First update prop_performance
  pm_settings$prop_performance <- prop_performance

  # After a small period of burn-in:
  if(iter > tune$n0){
    # 1. Update mixing weights continuously
    # Compute observed acceptance rates per proposal (scaled by mix)
    target_weights <- (prop_performance*(1-mix))/sum(prop_performance*(1-mix))
    # Smooth update
    mix_weights_new <- (1 - tune$mix_adapt)*pm_settings$mix + tune$mix_adapt*target_weights
    # Ensure minimum weight
    mix_weights_new <- pmax(mix_weights_new, 0.01)
    pm_settings$mix <- mix_weights_new / sum(mix_weights_new)

    # 2. Adapt number of particles based on ESS
    ess <- (sum(weights))^2 / sum(weights^2)
    desired_ess <- pm_settings$n_particles * tune$ESS_ratio
    scale_factor <- (desired_ess / ess)^tune$ESS_scale
    new_num_particles <- round(pm_settings$n_particles * scale_factor)
    pm_settings$n_particles <- max(10, min(200, new_num_particles))

    # 3. Update epsilon
    # If chosen_idx > 1 a new proposal was accepted
    pm_settings$epsilon <- update_epsilon(pm_settings$epsilon^2, chosen_idx > 1, tune$p_accept, iter, n_pars, tune$alphaStar)
  }
  return(pm_settings)
}


# Utility functions for sampling below ------------------------------------


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

check_mix <- function(mix, stage) {
  if (stage %in% c("burn", "adapt", "sample")) {
    default_mix <- c(0.1, 0.3, 0.6)
  } else{ #Preburn stage
    default_mix <- c(0.5, 0.5)
  }
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

check_epsilon <- function(epsilon, n_pars) {
  if (is.null(epsilon)) {
    if (n_pars > 15) {
      epsilon <- 0.5
    } else if (n_pars > 10) {
      epsilon <- 1
    } else {
      epsilon <- 1.5
    }
  }
  return(epsilon)
}

check_prop_performance <- function(prop_performance, stage){
  if (stage %in% c("burn", "adapt", "sample")) {
    default_mix <- c(0.1, 0.3, 0.6)
  } else{ #Preburn stage
    default_mix <- c(0.5, 0.5)
  }
  if(is.null(prop_performance)) prop_performance <- rep(0, length(default_mix))
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


condMVN <- function (mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE)
{
  if (missing(dependent.ind))
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given))
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(given.ind) == 0)
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(X.given) != length(given.ind))
    stop("lengths of `X.given' and `given.ind' must be same")
  if (check.sigma) {
    if (!isSymmetric(sigma))
      stop("sigma is not a symmetric matrix")
    eigenvalues <- eigen(sigma, only.values = TRUE)$values
    if (any(eigenvalues < 1e-08))
      stop("sigma is not positive-definite")
  }
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% chol2inv(chol(D))
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean = cMu, condVar = cVar)
}

calc_ll_manager <- function(proposals, dadm, model, component = NULL){
  if(!is.data.frame(dadm)){
    lls <- log_likelihood_joint(proposals, dadm, model, component)
  } else{
    model <- model()
    if(is.null(model$c_name)){ # use the R implementation
      lls <- apply(proposals,1, calc_ll_R, model, dadm = dadm)
    } else{
      p_types <- names(model$p_types)
      designs <- list()
      for(p in p_types){
        designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
      }
      constants <- attr(dadm, "constants")
      if(is.null(constants)) constants <- NA
      lls <- calc_ll(proposals, dadm, constants = constants, designs = designs, type = model$c_name,
                     model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                     model$trend)
    }
  }
  return(lls)
}

merge_group_level <- function(tmu, tmu_nuis, tvar, tvar_nuis, is_nuisance){
  n_pars <- length(is_nuisance)
  tmu_out <- numeric(n_pars)
  tmu_out[!is_nuisance] <- tmu
  tmu_out[is_nuisance] <- tmu_nuis
  tvar_out <- matrix(0, nrow = n_pars, ncol = n_pars)
  tvar_out[!is_nuisance, !is_nuisance] <- tvar
  tvar_out[is_nuisance, is_nuisance] <- tvar_nuis
  return(list(tmu = tmu_out, tvar = tvar_out))
}


run_hyper <- function(type, data, prior = NULL, iter = 1000, n_chains =3, ...){
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
  class(emc) <- "emc"
  return(emc)
}
