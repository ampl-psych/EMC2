#
sample_store_lm <- function(data, pars_fixed, pars_random, g_fixed, g_random,  subjects, iters = 1, stage = "init") {
  n_pars_fixed <- length(pars_fixed)
  n_pars_random <- length(pars_random)
  components <- c("between_subs", subjects)
  n_components <- length(components)
  samples <- list(
    fixed = array(NA_real_,dim = c(n_pars_fixed, iters), dimnames = list(pars_fixed, NULL)),
    random = array(NA_real_,dim = c(n_pars_random, iters),dimnames = list(pars_random, NULL)),
    g_fixed = array(NA_real_, c(length(unique(g_fixed)), iters), dimnames = list(unique(names(g_fixed)), NULL)),
    g_random = array(NA_real_, c(length(unique(g_random)), iters), dimnames = list(unique(names(g_random)), NULL)),
    a_half_fixed = array(NA_real_, c(length(unique(g_fixed)), iters), dimnames = list(unique(names(g_fixed)), NULL)),
    a_half_random = array(NA_real_, c(length(unique(g_random)), iters), dimnames = list(unique(names(g_random)), NULL)),
    epsilon = array(NA_real_,dim = c(n_components, iters),dimnames = list(components, NULL)),
    origin = array(NA_real_,dim = c(n_components, iters),dimnames = list(components, NULL)),
    stage = array(stage, iters),
    component_ll = array(NA_real_,dim = c(n_components, iters),dimnames = list(components, NULL))
  )
  return(samples)
}


make_dadm_lm <- function(data, DM_fixed, DM_random, full_data){
  s <- as.character(data$subjects[1])
  for(i in 1:length(DM_random)){
    curr_DM <- DM_random[[i]]
    if(!is.null(curr_DM)){
      curr_DM <- cbind(curr_DM[full_data$subject == s,!grepl("subjects", colnames(curr_DM)), drop = F],
                       curr_DM[full_data$subject == s,get_sub_idx(s, colnames(curr_DM)), drop = F])
      DM_random[[i]] <- curr_DM[,!apply(curr_DM, 2, sum) == 0, drop = F] # Remove empty columns that can come from nesting
    }
    DM_fixed[[i]] <- DM_fixed[[i]][full_data$subject == s,, drop = F]
  }
  DM <- mapply(cbind, DM_fixed, DM_random)
  for(i in 1:length(DM)){
    dm <- DM[[i]]
    cells <- apply(dm,1,paste,collapse="_")
    ass <- attr(dm,"assign")
    contr <- attr(dm,"contrasts")
    dups <- duplicated(cells)
    out <- as.matrix(dm[!dups,,drop=FALSE])
    attr(out,"expand") <- as.numeric(factor(cells,levels=unique(cells)))
    attr(out,"assign") <- ass
    attr(out,"contrasts") <- contr
    DM[[i]] <- out
  }
  attr(data, "designs") <- DM
  return(data)
}

dump_lm_attributes <- function(dadm){
  attr(dadm, "prior") <- NULL
  attr(dadm, "DM_fixed") <- NULL
  attr(dadm, "DM_random") <- NULL
  attr(dadm, "g_fixed") <- NULL
  attr(dadm, "g_random") <- NULL
  attr(dadm, "p_vector_fixed") <- NULL
  attr(dadm, "p_vector_random") <- NULL
  return(dadm)
}

pmwgs_lm <- function(dadm, pars = NULL, ll_func = NULL, prior = NULL, type, ...) {
  prior <- attr(dadm[[1]], "prior")
  DM_fixed <- attr(dadm[[1]], "DM_fixed")
  DM_random <- attr(dadm[[1]], "DM_random")
  g_fixed <- oneG(attr(dadm[[1]], "g_fixed"))
  g_random <- oneG(attr(dadm[[1]], "g_random"))
  pars_fixed <- attr(dadm[[1]], "p_vector_fixed")
  pars_random <- attr(dadm[[1]], "p_vector_random")
  pars_between <- c(pars_fixed, pars_random[!grepl("subjects", pars_random)])
  subjects <- unique(dadm[[1]]$subjects)
  ll_func <- attr(dadm[[1]], "model")()$log_likelihood
  dadm[[1]] <- dump_lm_attributes(dadm[[1]])
  split_dadm <- dm_list(dadm[[1]])
  split_dadm <- lapply(split_dadm, make_dadm_lm, DM_fixed, DM_random, dadm[[1]])
  is_intercept <- !grepl("_", pars_fixed)
  # Storage for the samples.
  samples <- sample_store_lm(dadm, pars_fixed, pars_random, g_fixed[!is_intercept], g_random, subjects)
  sampler <- list(
    data = split_dadm,
    dm_random = DM_random,
    dm_fixed = DM_fixed,
    pars_fixed = pars_fixed,
    pars_random = pars_random,
    g_map_fixed = g_fixed,
    g_map_random = g_random,
    is_intercept = is_intercept,
    subjects = subjects,
    n_subjects = length(subjects),
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler$prior <- get_prior_lm(prior, g_fixed, g_random, sample = F)
  return(sampler)
}

get_prior_lm <- function(prior, g_fixed, g_random, sample = F){
  intercepts <- names(g_fixed)[!grepl("_", names(g_fixed))]
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(prior$intercepts)){
    names(prior$intercepts_mu) <- intercepts # @Andrew how did you reorder the prior?
  } else{
    prior$intercepts_mu <- rep(0, length(intercepts))
    names(prior$intercepts_mu) <- intercepts
  }
  if(is.null(prior$intercepts_var)){
    prior$intercepts_var <-rep(1, length(intercepts))
  }
  if(is.null(prior$g_fixed_A)){
    prior$G_fixed_A <- rep(1, length(unique(g_fixed)))
  }
  if(is.null(prior$g_random_A)){
    prior$G_random_A <- rep(1, length(unique(g_random)))
  }
  if(is.null(prior$g_fixed_v)){
    prior$G_fixed_v <- 2
  }
  if(is.null(prior$g_random_v)){
    prior$G_random_v <- 2
  }
  return(prior)
}

#
get_startpoints_lm <- function(pmwgs){
  mu_fixed <- numeric(length(pmwgs$pars_fixed))
  var_fixed <- numeric(length(pmwgs$pars_fixed))
  mu_fixed[pmwgs$is_intercept] <- pmwgs$prior$intercepts_mu
  var_fixed[pmwgs$is_intercept] <- pmwgs$prior$intercepts_var
  # For now assume 0,1 normal for startpoints of effects
  mu_fixed[!pmwgs$is_intercept] <- 0
  var_fixed[!pmwgs$is_intercept] <- 1
  mu_random <- rep(0, length(pmwgs$pars_random))
  var_random <- rep(1, length(pmwgs$pars_random))
  names(mu_fixed) <- names(var_fixed) <- pmwgs$pars_fixed
  names(mu_random) <- names(var_random) <- pmwgs$pars_random
  return(list(mu_fixed = mu_fixed, var_fixed = var_fixed,
              mu_random = mu_random, var_random = var_random))
}

start_proposals_between <- function(n_particles, startpoints,
                                    pars_fixed, pars_random, data, subjects,
                                    ll_func, n_cores){
  #Draw the first start point
  within_idx <- grepl("subjects", pars_random)
  group_mean <- c(startpoints$mu_fixed, startpoints$mu_random[!within_idx])
  group_var <- c(startpoints$var_fixed, startpoints$var_random[!within_idx])
  proposals <- particle_draws(n_particles, group_mean, diag(group_var))
  proposals_list <- vector("list", length(subjects))
  colnames(proposals) <- c(pars_fixed, pars_random[!within_idx])
  random_pars <- numeric(length(pars_random))
  names(random_pars) <- pars_random
  sub_proposals<- lapply(subjects, add_within, proposals, random_pars)
  lw <- rowSums(parallel::mcmapply(calc_ll_manager, sub_proposals, data, MoreArgs = list(ll = ll_func), mc.cores = n_cores))
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx], origin = 2))
}

add_within <- function(s, proposals, random_pars){
  random_pars <-random_pars[get_sub_idx(s, names(random_pars))]
  random <- matrix(random_pars, ncol = length(random_pars), nrow = nrow(proposals), T)
  colnames(random) <- names(random_pars)
  return(cbind(proposals, random))
}

start_proposals_within <- function(pars, data, n_particles, startpoints, pars_between, ll_func){
  #Draw the first start point
  group_mean <- startpoints$mu_random[pars]
  group_var <-  startpoints$var_random[pars]
  proposals <- particle_draws(n_particles, group_mean, diag(group_var))
  colnames(proposals) <- pars
  between_pars <- matrix(pars_between, nrow = n_particles, ncol = length(pars_between), T)
  colnames(between_pars) <- names(pars_between)
  combined_proposals <- cbind(between_pars, proposals)
  lw <- calc_ll_manager(combined_proposals, dadm = data,ll_func)
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx], origin = 2))
}

unpack_within_proposals <- function(proposals_within){
  pars <- do.call(cbind, proposals_within[1,])
  pars <- c(t(pars)) # Drops the names, but I don't need them
  return(list(pars = pars,
              ll = do.call(c, proposals_within[2,]),
              origin = do.call(c, proposals_within[3,])))
}

fill_samples_lm <- function(samples, group_level, random_par_names,
                            proposals_between, proposals_within,
                            epsilon, j, stage){
  # Some bookkeeping to map the between and within pars to fixed and random again
  between_pars <- proposals_between$proposal
  rand_idx <- names(between_pars) %in% random_par_names
  fixed_pars <- between_pars[!rand_idx]
  random_pars <- numeric(length(random_par_names))
  between_idx <- random_par_names %in% names(between_pars[rand_idx])
  random_pars[between_idx] <- between_pars[rand_idx]
  random_pars[!between_idx] <- proposals_within$pars
  component_ll <- c(proposals_between$ll, proposals_within$ll)
  component_origin <- c(proposals_between$origin, proposals_within$origin)
  # And now we fill them in
  samples$fixed[,j] <- fixed_pars
  samples$random[,j] <- random_pars
  samples$g_fixed[,j] <- group_level$g_fixed
  samples$g_random[,j] <- group_level$g_random
  samples$a_half_fixed[,j] <- group_level$a_half_fixed
  samples$a_half_random[,j] <- group_level$a_half_random
  samples$epsilon[,j] <- epsilon
  samples$origin[,j] <- component_origin
  samples$component_ll[,j] <- component_ll
  samples$stage[j] <- stage
  samples$idx <- j
  return(samples)
}

#
init_lm <- function(pmwgs, verbose = FALSE, particles = 1000,
                 particles_fixed = round(particles*2/pmwgs$n_subjects), n_cores = 1, epsilon = NULL) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  startpoints <- get_startpoints_lm(pmwgs)
  if(is.null(epsilon)) epsilon <- rep(1, pmwgs$n_subjects + 1)
  if(length(epsilon) == 1) epsilon <- rep(epsilon, pmwgs$n_subjects + 1)
  proposals_between <- start_proposals_between(particles_fixed, startpoints, pmwgs$pars_fixed,
                                             pmwgs$pars_random, pmwgs$data, pmwgs$subjects,
                                             pmwgs$ll_func, n_cores)
  pars_per_subject <- lapply(pmwgs$subjects, FUN = function(x, pars) pars[get_sub_idx(x, pars)], pmwgs$pars_random)
  proposals_within <- parallel::mcmapply(start_proposals_within,pars_per_subject, pmwgs$data, MoreArgs = list(startpoints = startpoints,
                      n_particles = particles, pars_between = proposals_between$proposal, ll_func = pmwgs$ll_func), mc.cores = n_cores)
  proposals_within <- unpack_within_proposals(proposals_within)

  # Sample the mixture variables' initial values.

  pmwgs$samples <- fill_samples_lm(samples = pmwgs$samples, group_level = list(g_fixed = 1, g_random = 1, a_half_fixed = 1, a_half_random = 1),
                                   random_par_names = pmwgs$pars_random, proposals_between = proposals_between, proposals_within = proposals_within,
                                   epsilon = epsilon, j = 1, stage = "init")
  pmwgs$init <- TRUE
  return(pmwgs)
}

last_sample_lm <- function(sampler){
  list(
    g_fixed = sampler$g_fixed[, sampler$idx],
    g_random = sampler$g_random[, sampler$idx],
    a_half_fixed = sampler$a_half_fixed[,sampler$idx],
    a_half_random = sampler$a_half_random[, sampler$idx]
  )
}

get_par_sums_sq <- function(map, pars){
  out <- numeric(length(unique(map)))
  for(i in 1:length(unique(map))){
    idx <- map == map[i]
    out[i] <- sum(pars[idx]^2)
  }
  return(out)
}

get_n_obs <- function(map, pars){
  out <- numeric(length(unique(map)))
  for(i in 1:length(unique(map))){
    out[i] <- sum(map == map[i])
  }
  return(out)
}

gibbs_lm <- function(pmwgs, samples_fixed, samples_random){
  last <- last_sample_lm(pmwgs$samples)
  samples_fixed <- samples_fixed[!pmwgs$is_intercept]
  # Gibbs step for diagonal only
  prior <- pmwgs$prior
  sums_sq <- get_par_sums_sq(pmwgs$g_map_fixed[!pmwgs$is_intercept], samples_fixed)
  n_obs <- get_n_obs(pmwgs$g_map_fixed[!pmwgs$is_intercept], samples_fixed)
  tvinv_fixed = rgamma(n=length(sums_sq), shape=prior$G_fixed_v/2 + n_obs/2,
                 rate=prior$G_fixed_v/last$a_half_fixed + sums_sq / 2)
  tvar_fixed = 1/tvinv_fixed
  a_half_fixed <- 1 / rgamma(n = length(sums_sq), shape = (prior$G_fixed_v + 1) / 2,
                       rate = prior$G_fixed_v * tvinv_fixed + 1/(prior$G_fixed_A^2))
  # random
  sums_sq <- get_par_sums_sq(pmwgs$g_map_random, samples_random)
  n_obs <- get_n_obs(pmwgs$g_map_random, samples_random)
  tvinv_random = rgamma(n=length(sums_sq), shape=prior$G_random_v/2 + n_obs/2,
                       rate=prior$G_random_v/last$a_half_random + sums_sq / 2)
  tvar_random = 1/tvinv_random
  a_half_random <- 1 / rgamma(n = length(sums_sq), shape = (prior$G_random_v + 1) / 2,
                       rate = prior$G_random_v * tvinv_random + 1/(prior$G_random_A^2))
  return(list(g_fixed = tvar_fixed, a_half_fixed = a_half_fixed,
              g_random = tvar_random, a_half_random = a_half_random))
}

fix_epsilon_lm <- function(pmwgs, epsilon, force_prev_epsilon, components =1){
  if(is.null(epsilon) | force_prev_epsilon){
    epsilon <- attr(pmwgs, "epsilon")
    if(is.null(epsilon)){
      epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
    }
  }
  if(is.matrix(epsilon)){
    if(ncol(epsilon) != max(components)){
      epsilon <- matrix(rowMeans(epsilon), nrow = pmwgs$n_subjects +1, ncol = max(components))
    }
  } else if (length(epsilon) == 1 | is.vector(epsilon)){
    epsilon <- matrix(epsilon, nrow = pmwgs$n_subjects +1, ncol = max(components))
  }
  return(epsilon)
}

new_particle_between <- function(n_particles, hyper,
                                 pars_fixed, pars_random, data, subjects,
                                 ll_func, epsilon, stage, chains_cov = NULL,
                                 eff_props = NULL,
                                 mix_proportion = c(.1, .9, 0),
                                 n_cores){

  # Defining distributions
  within_idx <- grepl("subjects", names(pars_random))
  group_mean <- c(hyper$mu_fixed, hyper$mu_random[!within_idx])
  group_var <- c(hyper$var_fixed, hyper$var_random[!within_idx])
  prev_mu <- c(pars_fixed, pars_random[!within_idx])

  eff_mean <- eff_props$eff_mu
  eff_var <- eff_props$eff_var

  if(stage == "preburn"){
    eff_mean <- prev_mu
  }
  if(stage == "preburn"){
    eff_var_curr <- chains_cov * epsilon^2
    if(stage == "preburn"){
      ind_var <- diag(group_var) * epsilon^2
    } else{
      ind_var <- diag(group_var) *  mean(diag(eff_var_curr))/mean(diag(group_var))
    }
  } else{
    eff_var_curr <- eff_var * epsilon^2
    ind_var <- chains_cov *  epsilon^2
  }
  particle_numbers <- numbers_from_proportion(mix_proportion, n_particles)
  cumuNumbers <- cumsum(particle_numbers) + 1 # Include the particle from b4
  prior_particles <- particle_draws(particle_numbers[1], group_mean, diag(group_var))
  ind_particles <- particle_draws(particle_numbers[2], prev_mu, ind_var)
  if(mix_proportion[3] == 0){
    eff_particles <- NULL
  } else{
    eff_particles <- particle_draws(particle_numbers[3], eff_mean, eff_var_curr)
  }
  proposals <- rbind(prev_mu, prior_particles, ind_particles, eff_particles)
  proposals_list <- vector("list", length(subjects))
  colnames(proposals) <- c(names(pars_fixed), names(pars_random[!within_idx]))
  sub_proposals<- lapply(subjects, add_within, proposals, pars_random)
  lw <- rowSums(parallel::mcmapply(calc_ll_manager, sub_proposals, data, MoreArgs = list(ll = ll_func), mc.cores = n_cores))
  lp <- mvtnorm::dmvnorm(x = proposals, mean = group_mean, sigma = diag(group_var), log = TRUE)

  prop_density <- mvtnorm::dmvnorm(x = proposals, mean = prev_mu, sigma = ind_var)
  if (mix_proportion[3] == 0) {
    eff_density <- 0
  }
  else {
    eff_density <- mvtnorm::dmvnorm(x = proposals, mean = eff_mean, sigma = eff_var_curr)
  }
  lm <- log(mix_proportion[1] * exp(lp) + mix_proportion[2] * prop_density + (mix_proportion[3] * eff_density))
  l <- lw + lp - lm
  weight <- exp(l - max(l))
  idx <- sample(x = n_particles+1, size = 1, prob = weight)
  origin <- min(which(idx <= cumuNumbers))
  return(list(proposal = proposals[idx,], ll = lw[idx], origin = origin))
}

new_particle_within <- function(subj_mean, data, epsilon, chains_cov, eff_props, n_particles, hyper, pars_between, ll_func,
                                components = rep(1, length(subj_mean)),stage,
                                prev_ll, mix_proportion = c(.5, .5, 0)){
  #Draw the first start point
  unq_components <- unique(components)
  group_mean <- hyper$mu_random[names(subj_mean)]
  group_var <-  hyper$var_random[names(subj_mean)]

  eff_mean_sub <- eff_props$eff_mu
  eff_var <- eff_props$eff_var

  if(stage != "sample"){
    eff_mean_sub <- subj_mean
    group_var_subj <- group_var
  }
  out_lls <- numeric(length(unq_components))
  proposal_out <- numeric(length(subj_mean))

  for(i in unq_components){
    # In preburn use to proposals:
    # 1) person distribution with previous particle mean and scaled group variance
    # 2) group mean and group variance
    # In burn and sample add covariance proposals with
    # 3) previous particle mean and scaled chain covariance
    # In sample rather use:
    # 1) Previous particle mean and now scaled covariance
    # 2) Group mean and group variance (same as before)
    # 3) Conditional mean and covariance
    # Group distribution, with group mean and variance
    if(stage != "sample"){
      eff_var_curr <- chains_cov * epsilon[i]^2
      if(stage == "preburn"){
        var_subj <- diag(group_var_subj) * epsilon[i]^2
      } else{
        var_subj <- diag(group_var_subj) *  mean(diag(eff_var_curr))/mean(diag(group_var_subj))
      }
    } else{
      eff_var_curr <- eff_var
      var_subj <- chains_cov *  epsilon[i]^2
    }
    idx <- components == i
    # Draw new proposals for this component
    particle_numbers <- numbers_from_proportion(mix_proportion, n_particles)
    cumuNumbers <- cumsum(particle_numbers) + 1 # Include the particle from b4
    pop_particles <- particle_draws(particle_numbers[1], group_mean[idx], diag(group_var[idx]))
    ind_particles <- particle_draws(particle_numbers[2], subj_mean[idx], var_subj[idx, idx])
    if(mix_proportion[3] == 0){
      eff_particles <- NULL
    } else{
      eff_particles <- particle_draws(particle_numbers[3], eff_mean_sub[idx], eff_var_curr[idx, idx ])
    }
    # Rejoin new proposals with current MCMC values for other components
    proposals <- rbind(subj_mean, pop_particles, ind_particles, eff_particles)
    colnames(proposals) <- names(subj_mean)
    # We can safely add the between pars also the ones that aren't present in this subject's design matrix,
    # Since they are ignored when they're not present.
    between_pars <- matrix(pars_between, nrow = n_particles + 1, ncol = length(pars_between), T)
    colnames(between_pars) <- names(pars_between)
    combined_proposals <- cbind(between_pars, proposals)
    # Normally we assume that a component contains all the parameters to estimate the individual likelihood of a joint model
    # Sometimes we may also want to block within a model if it has very high dimensionality
    # Calculate likelihoods
    if(components[length(components)] > 1){
      lw <- calc_ll_manager(combined_proposals, dadm = data, ll_func, component = i)
    } else{
      lw <- calc_ll_manager(combined_proposals, dadm = data, ll_func)
    }
    lw_total <- lw + prev_ll - lw[1] # make sure lls from other components are included
    lp <- mvtnorm::dmvnorm(x = proposals[,idx], mean = group_mean[idx], sigma = diag(group_var[idx]), log = TRUE)
    prop_density <- mvtnorm::dmvnorm(x = proposals[,idx], mean = subj_mean[idx], sigma = var_subj[idx,idx])
    if (mix_proportion[3] == 0) {
      eff_density <- 0
    }
    else {
      eff_density <- mvtnorm::dmvnorm(x = proposals[,idx], mean = eff_mean_sub[idx], sigma = eff_var_curr[idx,idx])
    }
    lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * prop_density) + (mix_proportion[3] * eff_density))
    infnt_idx <- is.infinite(lm)
    lm[infnt_idx] <- min(lm[!infnt_idx])
    # Calculate weights and center
    l <- lw_total + lp - lm
    weights <- exp(l - max(l))
    # Do MH step and return everything
    idx_ll <- sample(x = n_particles+1, size = 1, prob = weights)
    origin <- min(which(idx_ll <= cumuNumbers))
    # ll <- ll + lw[idx_ll]
    out_lls[i] <- lw_total[idx_ll]
    proposal_out[idx] <- proposals[idx_ll,idx]
  } # Note that origin only contains last components origin, just used by devs anyway
  names(proposal_out) <- names(subj_mean)
  return(list(proposal = proposal_out, ll = mean(out_lls), origin = origin))
}


run_stage_lm <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 100,
                      verbose = TRUE,
                      verboseProgress = TRUE,
                      particles_fixed = round(particles*1.5),
                      n_cores = 1,
                      epsilon = NULL,
                      p_accept = NULL) {
  # Set defaults for NULL values
  # Set necessary local variables
  # Set stable (fixed) new_sample argument for this run
  components <- attr(pmwgs$data, "components")
  shared_ll_idx <- attr(pmwgs$data, "shared_ll_idx")
  # if(stage == "sample"){
  components <- rep(1, length(components))
  shared_ll_idx <- rep(1, length(shared_ll_idx))
  # }
  # Display stage to screen

  alphaStar=-qnorm(p_accept/2) #Idk about this one
  n0=round(5/(p_accept*(1-p_accept))) #Also not questioning this math for now

  epsilon <- fix_epsilon_lm(pmwgs, epsilon, force_prev_epsilon = T)
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  chains_cov_sub <- attr(pmwgs, "chains_cov_sub")
  eff_props <- attr(pmwgs, "eff_props")
  eff_props_between <-attr(pmwgs, "eff_props_between")
  if(is.null(chains_cov_sub)) chains_cov_sub <- vector("list", pmwgs$n_subjects)
  if(is.null(eff_props)) eff_props <- vector("list", pmwgs$n_subjects)
  chains_cov_between <- attr(pmwgs, "chains_cov_between")
  mix <- set_mix(stage)
  if (verboseProgress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx

  data <- pmwgs$data
  subjects <- pmwgs$subjects
  unq_components <- unique(components)
  block_idx <- block_variance_idx(components)
  g_intercept <-!grepl("_", unique(names(pmwgs$g_map_fixed)))
  g_fixed <- numeric(length(g_intercept))
  g_fixed[g_intercept] <- pmwgs$prior$intercepts_var

  # Main iteration loop
  for (i in 1:iter) {
    if (verboseProgress) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }
    j <- start_iter + i

    # Gibbs step
    pars <- gibbs_lm(pmwgs, pmwgs$samples$fixed[,j-1], pmwgs$samples$random[,j-1])

    # Create our normal distributions which are our priors on their parameters
    g_fixed[!g_intercept] <- pars$g_fixed
    hyper <- list(
      var_fixed = g_fixed[pmwgs$g_map_fixed],
      var_random = pars$g_random[pmwgs$g_map_random],
      mu_fixed = rep(0, length(pmwgs$g_map_fixed)),
      mu_random = rep(0, length(pmwgs$g_map_random))
    )

    names(hyper$mu_random) <- pmwgs$pars_random
    names(hyper$var_random) <- pmwgs$pars_random
    hyper$mu_fixed[pmwgs$is_intercept] = pmwgs$prior$intercepts_mu
    input_random <- pmwgs$samples$random[,j-1]
    if(stage == "preburn"){
      input_random[1:length(input_random)] <- 0
    } else{
      eff_props_between <- create_eff_proposals_between(pmwgs$samples)
    }
    # Particle step
    proposals_between <- new_particle_between(particles_fixed, hyper, pmwgs$samples$fixed[,j-1],
                                              input_random, pmwgs$data, pmwgs$subjects,
                                              pmwgs$ll_func, pmwgs$samples$epsilon[1, j -1], stage,  chains_cov = chains_cov_between,
                                              eff_props = eff_props_between, mix_proportion = mix,
                                              n_cores)
    pars_per_subject <- lapply(pmwgs$subjects, FUN = function(x, pars) pars[get_sub_idx(x, names(pars))], pmwgs$samples$random[,j-1])
    proposals_within <- parallel::mcmapply(new_particle_within,pars_per_subject, pmwgs$data, pmwgs$samples$epsilon[-1,j-1], chains_cov_sub, eff_props,
                                           MoreArgs = list(hyper = hyper, n_particles = particles, pars_between = proposals_between$proposal,
                                          ll_func = pmwgs$ll_func, stage = stage, prev_ll = pmwgs$samples$component_ll[1,j-1],
                                          mix_proportion = mix), mc.cores = n_cores)
    proposals_within <- unpack_within_proposals(proposals_within)

    # Sample the mixture variables' initial values.

    pmwgs$samples <- fill_samples_lm(samples = pmwgs$samples, group_level =  pars,
                                     random_par_names = pmwgs$pars_random, proposals_between = proposals_between, proposals_within = proposals_within,
                                     epsilon = epsilon, j = j, stage = stage)
    # Update epsilon
    if(!is.null(p_accept)){
        curr <- sapply(pmwgs$subjects, FUN = function(x, pars) pars[get_sub_idx(x, names(pars))][1], pmwgs$samples$random[,j])
        prev <- sapply(pmwgs$subjects, FUN = function(x, pars) pars[get_sub_idx(x, names(pars))][1], pmwgs$samples$random[,j-1])
        acc <- c(pmwgs$samples$fixed[1,j], curr) != c(pmwgs$samples$fixed[1,j-1], prev)
        epsilon <-update.epsilon(epsilon^2, acc, p_accept, j, length(pars_per_subject[[1]]), alphaStar)
    }
  }
  if (verboseProgress) close(pb)
  attr(pmwgs, "epsilon") <- epsilon
  return(pmwgs)
}
