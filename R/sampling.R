
pmwgs <- function(dadm, variant_funs, pars = NULL, ll_func = NULL, prior = NULL,
                  nuisance = NULL, nuisance_non_hyper = NULL, grouped_pars = NULL, ...) {
  if(is.data.frame(dadm)) dadm <- list(dadm)
  dadm <- extractDadms(dadm)
  if(is.null(pars)) pars <- dadm$pars
  if(is.null(ll_func)) ll_func <- dadm$ll_func
  if(is.null(prior)) prior <- dadm$prior
  dadm_list <-dadm$dadm_list
  # Storage for the samples.
  subjects <- sort(as.numeric(unique(dadm$subjects)))
  if(!is.null(grouped_pars) & !is.numeric(grouped_pars)) grouped_pars <- which(pars %in% grouped_pars)
  if(!is.null(nuisance) & !is.numeric(nuisance)) nuisance <- which(pars %in% nuisance)
  if(!is.null(nuisance_non_hyper) & !is.numeric(nuisance_non_hyper)) nuisance_non_hyper <- which(pars %in% nuisance_non_hyper)

  if(!is.null(grouped_pars)){
    is_grouped <- is.element(seq_len(length(pars)), grouped_pars)
    group_pars <- array(NA_real_, dim = c(sum(is_grouped),1), dimnames = list(pars[is_grouped], NULL))
    epsilon_grouped <- .5
  } else{
    is_grouped <- rep(F, length(pars))
    group_pars <- NULL
    epsilon_grouped <- NULL
  }

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
      samples = get_variant_funs(type)$sample_store(dadm, pars, integrate = F,
                                                    is_nuisance = !is_nuisance, is_grouped = is_grouped, ...),
      n_subjects = length(subjects),
      n_pars = sum(is_nuisance),
      nuisance = rep(F, sum(is_nuisance)),
      grouped = rep(F, sum(is_nuisance)),
      type = type
    )
    sampler_nuis <- get_variant_funs(type)$add_info(sampler_nuis, prior$prior_nuis, ...)
  }
  samples <- variant_funs$sample_store(dadm, pars, is_nuisance = is_nuisance,
                                       is_grouped = is_grouped, ...)
  samples$grouped_pars <- group_pars
  samples$epsilon_grouped <- epsilon_grouped

  sampler <- list(
    data = dadm_list,
    par_names = pars,
    subjects = subjects,
    n_pars = length(pars),
    nuisance = is_nuisance,
    n_subjects = length(subjects),
    ll_func = ll_func,
    samples = samples,
    grouped = is_grouped,
    sampler_nuis = sampler_nuis,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler <- variant_funs$add_info(sampler, prior, ...)
  sampler$prior$prior_grouped <- get_prior_single(prior$prior_grouped,
                                                  n_pars = sum(is_grouped), sample = F)
  attr(sampler, "variant_funs") <- variant_funs
  return(sampler)
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 verbose = FALSE, particles = 1000, n_cores = 1, epsilon = NULL) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  variant_funs <- attr(pmwgs, "variant_funs")
  startpoints <-startpoints_comb <- variant_funs$get_startpoints(pmwgs, start_mu, start_var)
  if(is.null(epsilon)) epsilon <- rep(set_epsilon(pmwgs$n_pars), pmwgs$n_subjects)
  if(length(epsilon) == 1) epsilon <- rep(epsilon, pmwgs$n_subjects)
  if(any(pmwgs$nuisance)){
    type <- pmwgs$sampler_nuis$type
    startpoints_nuis <- get_variant_funs(type)$get_startpoints(pmwgs$sampler_nuis, start_mu = NULL, start_var = NULL)
    startpoints_comb <- merge_group_level(startpoints$tmu, startpoints_nuis$tmu,
                                          startpoints$tvar, startpoints_nuis$tvar,
                                          pmwgs$nuisance[!pmwgs$grouped])
    pmwgs$sampler_nuis$samples <- get_variant_funs(type)$fill_samples(samples = pmwgs$sampler_nuis$samples,
                                                                      group_level = startpoints_nuis,
                                                                      epsilon = epsilon, j = 1,
                                                                      proposals = NULL,
                                                                      n_pars = pmwgs$n_pars)
    pmwgs$sampler_nuis$samples$idx <- 1
  }
  if(any(pmwgs$grouped)){
    grouped_pars <- mvtnorm::rmvnorm(particles, pmwgs$prior$prior_grouped$theta_mu_mean,
                                     pmwgs$prior$prior_grouped$theta_mu_var)
    colnames(grouped_pars) <- pmwgs$par_names[pmwgs$grouped]
  } else{
    grouped_pars <- NULL
  }
  proposals <- parallel::mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,
                                  parameters = startpoints_comb, n_particles = particles,
                                  pmwgs = pmwgs, variant_funs = variant_funs, grouped_pars = grouped_pars[1,], is_grouped = pmwgs$grouped,
                                  mc.cores = n_cores)
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars - sum(pmwgs$grouped) + 2, pmwgs$n_subjects))

  # Sample the mixture variables' initial values.

  pmwgs$samples <- variant_funs$fill_samples(samples = pmwgs$samples, group_level = startpoints, proposals = proposals,
                                             epsilon = epsilon, j = 1, n_pars = sum(!pmwgs$grouped))
  if(any(pmwgs$grouped)){
    grouped_pars <- start_proposals_group(pmwgs$data, grouped_pars, pmwgs$samples$alpha[,,1], pmwgs$par_names,
                                          pmwgs$ll_func, pmwgs$grouped, variant_funs, pmwgs$subjects, n_cores)

    pmwgs$samples$grouped_pars[,1] <- grouped_pars
  }
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
#' @examples \dontrun{
#' # Make a design and an emc object
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#'
#' DDMaE <- make_emc(forstmann, design_DDMaE)
#' # set up our mean starting points (same used across subjects).
#' mu <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'        t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' # Small variances to simulate start points from a tight range
#' var <- diag(0.05, length(mu))
#' # Initialize chains, 4 cores per chain, and parallelizing across our 3 chains as well
#' # so 4*3 cores used.
#' DDMaE <- init_chains(DDMaE, start_mu = p_vector, start_var = var, cores_per_chain = 4)
#' # Afterwards we can just use fit
#' DDMaE <- fit(DDMaE, cores_per_chain = 4)
#' }
#' @export
init_chains <- function(emc, start_mu = NULL, start_var = NULL, particles = 1000,
                        cores_per_chain=1,cores_for_chains = length(emc))
{
  attributes <- get_attributes(emc)
  emc <- mclapply(emc,init,start_mu = start_mu, start_var = start_var,
           verbose = FALSE, particles = particles,
           n_cores = cores_per_chain, mc.cores=cores_for_chains)
  emc <- get_attributes(emc, attributes)
  class(emc) <- "emc"
  return(emc)
}

start_proposals_group <- function(data, group_pars, alpha, par_names,
                                  likelihood_func, is_grouped,
                                  variant_funs, subjects, n_cores){
  num_particles <- nrow(group_pars)
  n_subjects <- length(subjects)
  proposals_list <- vector("list", n_subjects)
  for(i in 1:n_subjects){
    proposals_list[[i]] <- bind_alpha(group_pars, alpha[,i], num_particles, is_grouped, par_names)
  }
  lws <- parallel::mcmapply(calc_ll_for_group, proposals_list, data, MoreArgs = list(ll = likelihood_func), mc.cores = n_cores)
  lw <- rowSums(lws)
  weight <- exp(lw - max(lw))
  idx <- sample(x = num_particles, size = 1, prob = weight)
  return(group_pars[idx,])
}

start_proposals <- function(s, parameters, n_particles, pmwgs, variant_funs, grouped_pars, is_grouped){
  #Draw the first start point
  group_pars <- variant_funs$get_group_level(parameters, s)
  proposals <- particle_draws(n_particles, group_pars$mu, group_pars$var)
  colnames(proposals) <- rownames(pmwgs$samples$alpha) # preserve par names
  if(any(is_grouped)){
    proposals <- update_proposals_grouped(proposals, grouped_pars, is_grouped,
                                          par_names = colnames(proposals))
  }
  lw <- calc_ll_manager(proposals, dadm = pmwgs$data[[which(pmwgs$subjects == s)]],
                        ll_func = pmwgs$ll_func)
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,!is_grouped], ll = lw[idx], origin = 2))
}

run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 100,
                      verbose = TRUE,
                      verboseProgress = TRUE,
                      force_prev_epsilon = TRUE,
                      n_cores = 1,
                      epsilon = NULL,
                      p_accept = NULL) {
  # Set defaults for NULL values
  # Set necessary local variables
  # Set stable (fixed) new_sample argument for this run
  n_pars <- pmwgs$n_pars
  components <- attr(pmwgs$data, "components")
  shared_ll_idx <- attr(pmwgs$data, "shared_ll_idx")

  alphaStar=-qnorm(p_accept/2) #Idk about this one
  n0=round(5/(p_accept*(1-p_accept)))[1] #Also not questioning this math for now

  epsilon <- fix_epsilon(pmwgs, epsilon, force_prev_epsilon, components)
  if(length(particles) == 1){
    particles <- rep(particles, max(pmwgs$n_subjects, 2)) # kluge to keep it as a vector
  }
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  eff_mu <- attr(pmwgs, "eff_mu")
  eff_var <- attr(pmwgs, "eff_var")
  # eff_alpha <- attr(pmwgs, "eff_alpha")
  # eff_tau <- attr(pmwgs, "eff_tau")
  if(is.null(eff_mu)) eff_mu <- vector("list", pmwgs$n_subjects)
  if(is.null(eff_var)) eff_var <- vector("list", pmwgs$n_subjects)
  # if(is.null(eff_alpha)) eff_alpha <- vector("list", pmwgs$n_subjects)
  # if(is.null(eff_tau)) eff_tau <- vector("list", pmwgs$n_subjects)
  chains_cov <- attr(pmwgs, "chains_cov")
  if(is.null(chains_cov)) chains_cov <- vector("list", pmwgs$n_subjects)

  chains_cov_grouped <- attr(pmwgs, "chains_cov_grouped")
  mix <- set_mix(stage)
  if (verboseProgress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx

  data <- pmwgs$data
  subjects <- pmwgs$subjects
  unq_components <- unique(components)
  variant_funs <- attr(pmwgs, "variant_funs")
  nuisance <- pmwgs$nuisance
  if(any(nuisance)){
    type <- pmwgs$sampler_nuis$type
    pmwgs$sampler_nuis$samples$idx <- pmwgs$samples$idx
  }
  grouped <- pmwgs$grouped
  if(any(grouped)){
    epsilon_grouped <- pmwgs$samples$epsilon_grouped
    mix_grouped <- c(0.1, 0.9, 0)
    particles_grouped <- round(sqrt(nrow(pmwgs$samples$grouped_pars))*40)
  }
  block_idx <- block_variance_idx(components[!grouped])
  # Main iteration loop
  for (i in 1:iter) {
    if (verboseProgress) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }
    j <- start_iter + i

    # Gibbs step
    pars <- pars_comb <- variant_funs$gibbs_step(pmwgs, pmwgs$samples$alpha[!nuisance[!grouped],,j-1])
    if(any(nuisance)){
      pars_nuis <- get_variant_funs(type)$gibbs_step(pmwgs$sampler_nuis, pmwgs$samples$alpha[nuisance[!grouped],,j-1])
      pars_comb <- merge_group_level(pars$tmu, pars_nuis$tmu, pars$tvar, pars_nuis$tvar, nuisance[!grouped])
      pars_comb$alpha <- pmwgs$samples$alpha[,,j-1]
      pmwgs$sampler_nuis$samples <- get_variant_funs(type)$fill_samples(samples = pmwgs$sampler_nuis$samples,
                                                                        group_level = pars_nuis,
                                                                        epsilon = epsilon, j = j,
                                                                        proposals = NULL,
                                                                        n_pars = n_pars)
      pmwgs$sampler_nuis$samples$idx <- j
    }
    if(any(grouped)){
      subj_prior <-sum(apply(pars_comb$alpha, 2, FUN = function(x) mvtnorm::dmvnorm(x, pars$tmu, sigma = pars$tvar, log = TRUE)))
      grouped_pars <- new_particle_group(pmwgs$data, particles_grouped, pmwgs$prior$prior_grouped,
                                         chains_cov_grouped, mix_grouped, epsilon_grouped,
                                         pmwgs$samples$grouped_pars[,j-1] , pars_comb$alpha, pmwgs$par_names,
                                         pmwgs$ll_func, pmwgs$grouped, stage, variant_funs, pmwgs$subjects, subj_prior, n_cores)
      pmwgs$samples$grouped_pars[,j] <- grouped_pars$proposal
      pmwgs$samples$epsilon_grouped <- epsilon_grouped
    } else{
      grouped_pars <- NULL
    }

    # Particle step
    proposals <- parallel::mcmapply(new_particle, 1:pmwgs$n_subjects, data, particles, eff_mu, eff_var,
                                    chains_cov,
                                    pmwgs$samples$subj_ll[,j-1],
                                    MoreArgs = list(pars_comb, mix, pmwgs$ll_func, epsilon, components, stage,
                                                    variant_funs$get_group_level, block_idx, shared_ll_idx, grouped_pars$proposal, grouped,
                                                    group_prior = grouped_pars$prior),
                                    mc.cores =n_cores)
    proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars - sum(grouped) + 2, pmwgs$n_subjects))

    #Fill samples
    pmwgs$samples <- variant_funs$fill_samples(samples = pmwgs$samples, group_level = pars,
                                               proposals = proposals, epsilon = rowMeans(epsilon), j = j, n_pars = sum(!pmwgs$grouped))

    # Update epsilon
    if(!is.null(p_accept)){
      if(j > n0){
        for(component in unq_components){
          idx <- components[!grouped] == component
          acc <-  pmwgs$samples$alpha[max(which(idx)),,j] != pmwgs$samples$alpha[max(which(idx)),,(j-1)]
          epsilon[,component] <-update_epsilon(epsilon[,component]^2, acc, p_accept, j, sum(idx), alphaStar)
        }
        if(any(grouped)){
          acc <-  pmwgs$samples$grouped_pars[1,j] !=  pmwgs$samples$grouped_pars[1,(j-1)]
          epsilon_grouped <-update_epsilon(epsilon_grouped^2, acc, mean(p_accept), j, sum(grouped), mean(alphaStar))
        }
      }
    }
  }
  if (verboseProgress) close(pb)
  attr(pmwgs, "epsilon") <- epsilon
  return(pmwgs)
}


new_particle <- function (s, data, num_particles, eff_mu = NULL,
                          eff_var = NULL, chains_cov, prev_ll,
                          parameters, mix_proportion = c(0.5, 0.5, 0),
                          likelihood_func = NULL, epsilon = NULL,
                          components, stage,  group_level_func,
                          block_idx, shared_ll_idx, grouped_pars, is_grouped,
                          group_prior)
{
  # if(stage == "sample"){
  #   if(rbinom(1, size = 1, prob = .5) == 1){
  #     components <- rep(1, length(components))
  #     shared_ll_idx <- rep(1, length(components))
  #   }
  # }
  group_pars <- group_level_func(parameters, s)
  unq_components <- unique(components)
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
      eff_var <- chains_cov * epsilon[s,i]^2
      var_subj <- group_var_subj *  epsilon[s,i]^2
    } else{
      eff_var <- eff_var_old * epsilon[s,i]^2
      var_subj <- chains_cov *  epsilon[s,i]^2
    }
    idx <- components[!is_grouped] == i
    # Draw new proposals for this component
    particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
    cumuNumbers <- cumsum(particle_numbers) + 1 # Include the particle from b4
    pop_particles <- particle_draws(particle_numbers[1], group_mu[idx], group_var[idx, idx])
    ind_particles <- particle_draws(particle_numbers[2], subj_mu[idx], var_subj[idx, idx])
    if(mix_proportion[3] == 0){
      eff_particles <- NULL
    } else{
      eff_particles <- particle_draws(particle_numbers[3], eff_mu[idx], eff_var[idx, idx ])# eff_alpha, eff_tau)
    }
    # Rejoin new proposals with current MCMC values for other components
    proposals <- matrix(rep(subj_mu, num_particles + 1), nrow = num_particles + 1, byrow = T)
    colnames(proposals) <- names(subj_mu)
    proposals[2:(num_particles+1),idx] <- rbind(pop_particles, ind_particles, eff_particles)
    ll_proposals <- proposals
    if(any(is_grouped)){
      ll_proposals <- update_proposals_grouped(proposals, grouped_pars, is_grouped,
                                               par_names = names(subj_mu))
    }
    # Normally we assume that a component contains all the parameters to estimate the individual likelihood of a joint model
    # Sometimes we may also want to block within a model if it has very high dimensionality
    shared_idx <- shared_ll_idx[idx][1]
    is_shared <- shared_idx == shared_ll_idx

    # Calculate likelihoods
    if(components[length(components)] > 1){
      lw <- calc_ll_manager(ll_proposals[,is_shared], dadm = data, likelihood_func, component = shared_idx)
    } else{
      lw <- calc_ll_manager(ll_proposals[,is_shared], dadm = data, likelihood_func)
    }
    lw_total <- lw + prev_ll - lw[1] # make sure lls from other components are included
    lp <- mvtnorm::dmvnorm(x = proposals[,idx], mean = group_mu[idx], sigma = group_var[idx,idx], log = TRUE)
    prop_density <- mvtnorm::dmvnorm(x = proposals[,idx], mean = subj_mu[idx], sigma = var_subj[idx,idx])
    if(length(unq_components) > 1){
      prior_density <- mvtnorm::dmvnorm(x = proposals, mean = group_mu, sigma = group_var, log = TRUE)
    } else{
      prior_density <- lp
    }
    if(any(is_grouped)){
      prior_density <- prior_density + group_prior
    }

    if (mix_proportion[3] == 0) {
      eff_density <- 0
    }
    else {
      #if(is.null(eff_alpha)){
      eff_density <- mvtnorm::dmvnorm(x = proposals[,idx], mean = eff_mu[idx], sigma = eff_var[idx,idx])
      # } else{
      #   eff_density <- sn::dmsn(x = proposals[,idx], xi = eff_mu[idx], Omega = eff_var[idx,idx],
      #                            alpha = eff_alpha, tau = eff_tau)
      # }
    }
    lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * prop_density) + (mix_proportion[3] * eff_density))
    infnt_idx <- is.infinite(lm)
    lm[infnt_idx] <- min(lm[!infnt_idx])
    # Calculate weights and center
    l <- lw_total + prior_density - lm
    weights <- exp(l - max(l))
    # Do MH step and return everything
    idx_ll <- sample(x = num_particles+1, size = 1, prob = weights)
    origin <- min(which(idx_ll <= cumuNumbers))
    out_lls[i] <- lw[idx_ll]
    proposal_out[idx] <- proposals[idx_ll,idx]
  } # Note that origin only contains last components origin, just used by devs anyway
  return(list(proposal = proposal_out, ll = sum(out_lls), origin = origin))
}

new_particle_group <- function(data, num_particles, prior,
                               chains_cov, mix_proportion = c(.1, .9), epsilon_grouped,
                               prev_mu, alpha, par_names, likelihood_func = NULL,
                               is_grouped, stage,
                               variant_funs, subjects, subj_prior, n_cores){
  prior_mu <- prior$theta_mu_mean
  prior_var <- prior$theta_mu_var
  if(stage == "preburn"){
    chains_cov <- prior_var
  }
  chains_cov <- chains_cov * epsilon_grouped^2
  particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
  prior_particles <- particle_draws(particle_numbers[1], prior_mu, prior_var)
  ind_particles <- particle_draws(particle_numbers[2], prev_mu, chains_cov)
  proposals <- rbind(prev_mu, prior_particles, ind_particles)
  n_subjects <- length(subjects)
  proposals_list <- vector("list", n_subjects)
  for(i in 1:n_subjects){
    proposals_list[[i]] <- bind_alpha(proposals, alpha[,i], num_particles + 1, is_grouped, par_names)
  }
  lws <- parallel::mcmapply(calc_ll_for_group, proposals_list, data, MoreArgs = list(ll = likelihood_func), mc.cores = n_cores)
  lw <- rowSums(lws)
  lp <- mvtnorm::dmvnorm(x = proposals, mean = prior_mu, sigma = prior_var, log = TRUE)
  prop_density <- mvtnorm::dmvnorm(x = proposals, mean = prev_mu, sigma = chains_cov)
  lm <- log(mix_proportion[1] * exp(lp) + mix_proportion[2] * prop_density)
  prior_density <- lp + subj_prior
  l <- lw + prior_density - lm
  weights <- exp(l - max(l))
  idx <- sample(x = num_particles + 1, size = 1, prob = weights)
  return(list(proposal = proposals[idx,], prior = lp[idx]/length(subjects)))
}




calc_ll_for_group <- function(proposals, data, ll){
  lw <- calc_ll_manager(proposals, dadm = data, ll)
}


bind_alpha <- function(proposals, alpha, num_particles, is_grouped, par_names){
  tmp <- matrix(0, ncol = length(is_grouped), nrow = num_particles)
  tmp[,!is_grouped] <- matrix(alpha, ncol = length(alpha), nrow = num_particles, byrow = T)
  tmp[,is_grouped] <- proposals
  colnames(tmp) <- par_names
  return(tmp)
}

# Utility functions for sampling below ------------------------------------


update_epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c <- ((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta <- log(sqrt(epsilon2))
  Theta <- Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}


numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
  numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_proportion)
  if (mix_proportion[3] == 0) {
    numbers[3] <- 0
    numbers[2] <- num_particles - numbers[1]
  } else {
    numbers[3] <- num_particles - sum(numbers)
  }
  return(numbers)
}


particle_draws <- function(n, mu, covar, alpha = NULL, tau= NULL) {
  if (n <= 0) {
    return(NULL)
  }
  if(is.null(alpha)){
    return(mvtnorm::rmvnorm(n, mu, covar))
  }
  # else{
  #   return(sn::rmsn(n, xi = mu, Omega = covar, alpha = alpha, tau = tau))
  # }
}

fix_epsilon <- function(pmwgs, epsilon, force_prev_epsilon, components){
  if(is.null(epsilon) | force_prev_epsilon){
    epsilon <- attr(pmwgs, "epsilon")
    if(is.null(epsilon)){
      epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
    }
  }
  if(is.matrix(epsilon)){
    if(ncol(epsilon) != max(components)){
      epsilon <- matrix(rowMeans(epsilon), nrow = pmwgs$n_subjects, ncol = max(components))
    }
  } else if (length(epsilon) == 1 | is.vector(epsilon)){
    epsilon <- matrix(epsilon, nrow = pmwgs$n_subjects, ncol = max(components))
  }
  return(epsilon)
}

set_epsilon <- function(n_pars) {
  if (n_pars > 15) {
    epsilon <- 0.5
  } else if (n_pars > 10) {
    epsilon <- 1
  } else {
    epsilon <- 1.5
  }
  return(epsilon)
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
    epsilon = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    origin = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
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

fill_samples_base <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  # Fill samples both group level and random effects
  samples$theta_mu[, j] <- group_level$tmu
  samples$theta_var[, , j] <- group_level$tvar
  if(!is.null(proposals)) samples <- fill_samples_RE(samples, proposals, epsilon,j, n_pars)
  return(samples)
}



fill_samples_RE <- function(samples, proposals, epsilon, j = 1, n_pars, ...){
  # Only for random effects, separated because group level sometimes differs.
  if(!is.null(proposals)){
    samples$alpha[, , j] <- proposals[1:n_pars,]
    samples$subj_ll[, j] <- proposals[n_pars + 1,]
    samples$origin[,j] <- proposals[n_pars + 2,]
    samples$idx <- j
    samples$epsilon[,j] <- epsilon
  }
  return(samples)
}

set_mix <- function(stage) {
  if (stage %in% c("burn", "adapt")) {
    mix <- c(0.15, 0.15, 0.7)
  } else if(stage == "sample"){
    mix <- c(0.1, 0.3, 0.6)
  } else{ #Preburn stage
    mix <- c(0.5, 0.5, 0)
  }
  return(mix)
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


# some sort of variants, used to be called via a list called variant_funs

get_variant_funs <- function(type = "standard") {
  if(type == "standard") {
    list_fun <- list(# store functions
      sample_store = sample_store_standard,
      get_prior = get_prior_standard,
      add_info = add_info_standard,
      get_startpoints = get_startpoints_standard,
      get_group_level = get_group_level_standard,
      fill_samples = fill_samples_standard,
      gibbs_step = gibbs_step_standard,
      group_IC = group__IC_standard,
      filtered_samples = filtered_samples_standard,
      get_conditionals = get_conditionals_standard,
      get_all_pars_IS2 = get_all_pars_standard,
      prior_dist_IS2 = prior_dist_standard,
      group_dist_IS2 = group_dist_standard,
      bridge_add_group = bridge_add_group_standard,
      bridge_add_info = bridge_add_info_standard,
      bridge_group_and_prior_and_jac = bridge_group_and_prior_and_jac_standard
    )
  } else if(type == "single"){
    list_fun <- list(# store functions
      sample_store = sample_store_base,
      get_prior = get_prior_single,
      add_info = add_info_single,
      get_startpoints = get_startpoints_single,
      get_group_level = get_group_level_single,
      fill_samples = fill_samples_RE,
      gibbs_step = gibbs_step_single,
      group_IC = group__IC_single,
      filtered_samples = filtered_samples_single,
      get_conditionals = get_conditionals_single,
      bridge_add_group = bridge_add_group_single,
      bridge_add_info = bridge_add_info_single,
      bridge_group_and_prior_and_jac = bridge_group_and_prior_and_jac_single
    )
  } else if(type == "blocked"){
    list_fun <- list(# store functions
      sample_store = sample_store_standard,
      add_info = add_info_blocked,
      get_prior = get_prior_blocked,
      get_startpoints = get_startpoints_blocked,
      get_group_level = get_group_level_standard,
      fill_samples = fill_samples_standard,
      gibbs_step = gibbs_step_blocked,
      group_IC = group__IC_standard,
      filtered_samples = filtered_samples_standard,
      get_conditionals = get_conditionals_blocked,
      get_all_pars_IS2 = get_all_pars_blocked,
      prior_dist_IS2 = prior_dist_blocked,
      group_dist_IS2 = group_dist_blocked,
      bridge_add_group = bridge_add_group_blocked,
      bridge_add_info = bridge_add_info_blocked,
      bridge_group_and_prior_and_jac = bridge_group_and_prior_and_jac_blocked
    )
  } else if(type == "diagonal"){
    list_fun <- list(# store functions
      sample_store = sample_store_standard,
      add_info = add_info_diag,
      get_startpoints = get_startpoints_diag,
      get_group_level = get_group_level_standard,
      fill_samples = fill_samples_standard,
      gibbs_step = gibbs_step_diag,
      group_IC = group__IC_standard,
      filtered_samples = filtered_samples_standard,
      get_conditionals = get_conditionals_diag,
      get_all_pars_IS2 = get_all_pars_standard,
      prior_dist_IS2 = prior_dist_diag,
      group_dist_IS2 = group_dist_diag,
      bridge_add_group = bridge_add_group_diag,
      bridge_add_info = bridge_add_info_diag,
      bridge_group_and_prior_and_jac = bridge_group_and_prior_and_jac_diag
    )
  } else if(type == "factor"){
    list_fun <- list(# store functions
      sample_store = sample_store_factor,
      add_info = add_info_factor,
      get_startpoints = get_startpoints_factor,
      get_prior = get_prior_factor,
      get_group_level = get_group_level_standard,
      fill_samples = fill_samples_factor,
      gibbs_step = gibbs_step_factor,
      group_IC = group__IC_standard,
      filtered_samples = filtered_samples_factor,
      get_conditionals = get_conditionals_factor,
      bridge_add_group = bridge_add_group_factor,
      bridge_add_info = bridge_add_info_factor,
      bridge_group_and_prior_and_jac = bridge_group_and_prior_and_jac_factor
    )
  } else if(type == "lm"){
    list_fun <- list(# store functions
      sample_store = sample_store_lm,
      add_info = add_info_lm,
      get_startpoints = get_startpoints_lm,
      get_group_level = get_group_level_lm,
      fill_samples = fill_samples_lm,
      gibbs_step = gibbs_step_lm,
      group_IC = group__IC_standard,
      filtered_samples = filtered_samples_lm,
      get_conditionals = get_conditionals_lm,
      get_all_pars_IS2 = get_all_pars_lm,
      prior_dist_IS2 = prior_dist_lm,
      group_dist_IS2 = group_dist_lm
    )
  }
  else if(type == "infnt_factor"){
    list_fun <- list(# store functions
      sample_store = sample_store_infnt_factor,
      add_info = add_info_infnt_factor,
      get_startpoints = get_startpoints_infnt_factor,
      get_prior = get_prior_factor,
      get_group_level = get_group_level_standard,
      fill_samples = fill_samples_infnt_factor,
      gibbs_step = gibbs_step_infnt_factor,
      group_IC = group__IC_standard,
      filtered_samples = filtered_samples_infnt_factor,
      get_conditionals = get_conditionals_infnt_factor,
      get_all_pars_IS2 = get_all_pars_infnt_factor,
      prior_dist_IS2 = prior_dist_infnt_factor,
      group_dist_IS2 = group_dist_infnt_factor
    )
  }
  else if(type == "SEM"){
    list_fun <- list(# store functions
      sample_store = sample_store_SEM,
      add_info = add_info_SEM,
      get_startpoints = get_startpoints_SEM,
      get_group_level = get_group_level_SEM,
      fill_samples = fill_samples_SEM,
      gibbs_step = gibbs_step_SEM,
      group_IC = group__IC_SEM,
      filtered_samples = filtered_samples_SEM,
      get_conditionals = get_conditionals_SEM,
      get_all_pars_IS2 = get_all_pars_infnt_factor,
      prior_dist_IS2 = prior_dist_infnt_factor,
      group_dist_IS2 = group_dist_infnt_factor,
      bridge_add_group = bridge_add_group_SEM,
      bridge_add_info = bridge_add_info_SEM,
      bridge_group_and_prior_and_jac = bridge_group_and_prior_and_jac_SEM
    )
  }
  list_fun$type <- type
  return(list_fun)
}

calc_ll_manager <- function(proposals, dadm, ll_func, component = NULL){
  if(!is.data.frame(dadm)){
    lls <- log_likelihood_joint(proposals, dadm, component)
  } else{
    c_name <- attr(dadm,"model")()$c_name
    if(is.null(c_name)){ # use the R implementation
      lls <- apply(proposals,1, ll_func,dadm = dadm)
    } else{
      p_types <- names(attr(dadm,"model")()$p_types)
      designs <- list()
      for(p in p_types){
        designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
      }
      constants <- attr(dadm, "constants")
      if(is.null(constants)) constants <- NA
      if(c_name == "DDM"){
        levels(dadm$R) <- c(0,1)
        pars <- get_pars_matrix(proposals[1,],dadm)
        pars <- cbind(pars, dadm$R)
        parameter_char <- apply(pars, 1, paste0, collapse = "\t")
        parameter_factor <- factor(parameter_char, levels = unique(parameter_char))
        parameter_indices <- split(seq_len(nrow(pars)), f = parameter_factor)
        names(parameter_indices) <- 1:length(parameter_indices)
      } else{
        parameter_indices <- list()
      }
      lls <- calc_ll(proposals, dadm, constants = constants, designs = designs, type = c_name,
                     p_types = p_types, min_ll = log(1e-10), group_idx = parameter_indices)
    }
  }
  return(lls)
}

update_proposals_grouped <- function(proposals, grouped_pars, is_grouped, par_names){
  proposals_full <- matrix(0, nrow = nrow(proposals), ncol = length(is_grouped))
  proposals_full[,!is_grouped] <- proposals
  proposals_full[,is_grouped] <- matrix(grouped_pars, ncol = length(grouped_pars), nrow = nrow(proposals), byrow = T)
  colnames(proposals_full)[!is_grouped] <- par_names
  colnames(proposals_full)[is_grouped] <- names(grouped_pars)
  return(proposals_full)
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


run_hyper <- function(type, data, prior = NULL, iter = 5000, ...){
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
  variant_funs <- get_variant_funs(type)
  samples <- variant_funs$sample_store(data = data ,par_names = pars, is_nuisance = rep(F, length(pars)), integrate = F, is_grouped =
                                         rep(F, length(pars)), ...)
  subjects <- unique(data$subjects)
  sampler <- list(
    data = data,
    par_names = pars,
    subjects = subjects,
    n_pars = length(pars),
    nuisance = rep(F, length(pars)),
    grouped = rep(F, length(pars)),
    n_subjects = length(subjects),
    samples = samples,
    init = TRUE
  )
  class(sampler) <- "pmwgs"
  sampler <- variant_funs$add_info(sampler, prior, ...)
  startpoints <- variant_funs$get_startpoints(sampler, start_mu = NULL, start_var = NULL)
  sampler$samples <- variant_funs$fill_samples(samples = sampler$samples, group_level = startpoints, proposals = NULL,
                                               epsilon = NA, j = 1, n_pars = sampler$n_pars)
  sampler$samples$idx <- 1
  sampler <- extend_sampler(sampler, iter-1, "sample")
  for(i in 2:iter){
    if(is_mcmc){
      pars <- variant_funs$gibbs_step(sampler, data_input[,,i])
    } else{
      pars <- variant_funs$gibbs_step(sampler, t(data_input))
    }
    sampler$samples$idx <- i
    sampler$samples <- variant_funs$fill_samples(samples = sampler$samples, group_level = pars, proposals = NULL,
                                                 epsilon = NA, j = i, n_pars = sampler$n_pars)
  }
  return(sampler)
}





