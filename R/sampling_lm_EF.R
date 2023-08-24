#
pmwgs_ES <- function(dadm, pars = NULL, ll_func = NULL, prior = NULL, ...) {
  if(is.data.frame(dadm)) dadm <- list(dadm)
  dadm <- extractDadms(dadm)
  if(is.null(pars)) pars <- dadm$pars
  if(is.null(ll_func)) ll_func <- dadm$ll_func
  if(is.null(prior)) prior <- dadm$prior
  dadm_list <-dadm$dadm_list
  # Storage for the samples.
  subjects <- sort(as.numeric(unique(dadm$subjects)))
  samples <- sample_store_ES(dadm, random_effects, fixed_effects)
  sampler <- list(
    data = dadm_list,
    dm_random = dm_random,
    dm_fixed = dm_fixed,
    subjects = subjects,
    n_pars = length(pars),
    n_subjects = length(subjects),
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler <- add_info_ES(sampler, prior)
  return(sampler)
}
#
init_ES <- function(pmwgs, verbose = FALSE, particles_random = 1000,
                 particles_fixed = round(particles_random*2/pmwgs$n_subjects), n_cores = 1, epsilon = NULL) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  startpoints_fixed <- get_startpoints_fixed(pmwgs)
  startpoints_random <- get_startpoints_random(pmwgs)
  if(is.null(epsilon)) epsilon <- rep(set_epsilon(pmwgs$n_pars), pmwgs$n_subjects)
  if(length(epsilon) == 1) epsilon <- rep(epsilon, pmwgs$n_subjects)
  epsilon_fixed <-  epsilon[1]
  proposals <- parallel::mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals_random,
                                  parameters = startpoints_comb, n_particles = particles,
                                  pmwgs = pmwgs, variant_funs = variant_funs,
                                  mc.cores = n_cores)
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars - sum(pmwgs$grouped) + 2, pmwgs$n_subjects))

  # Sample the mixture variables' initial values.

  pmwgs$samples <- fill_samples_ES(samples = pmwgs$samples, group_level = startpoints, proposals = proposals,
                                             epsilon = epsilon, j = 1)
  pmwgs$init <- TRUE
  return(pmwgs)
}
#
start_proposals_fixed <- function(n_particles, pmwgs, fixed_pars, random_pars){
  #Draw the first start point
  proposals <- particle_draws(n_particles, rep(0, pmwgs$n_fixed), fixed_pars$var)
  colnames(proposals) <- rownames(pmwgs$samples$alpha) # preserve par names
  proposals <- add_random(proposals, fixed_pars, s)
  for(i in 1:n_subjects){
    proposals_list[[i]] <- add_random(proposals, random_pars, i)
  }
  lws <- parallel::mclapply(1:n_subjects, calc_ll_for_group, proposals_list, data, likelihood_func, subjects, mc.cores = n_cores)
  lw <- Reduce("+", lws)
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx], origin = 2))
}

start_proposals_random <- function(s, n_particles, pmwgs, fixed_pars, random_pars){
  #Draw the first start point
  proposals <- particle_draws(n_particles, rep(0, pmwgs$n_random), random_pars$var)
  colnames(proposals) <- rownames(pmwgs$samples$alpha) # preserve par names
  proposals <- add_fixed(proposals, fixed_pars, s)
  lw <- calc_ll_manager(proposals, dadm = pmwgs$data[[which(pmwgs$subjects == s)]],
                        ll_func = pmwgs$ll_func)
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx], origin = 2))
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
  components <- attr(pmwgs$data, "components")
  shared_ll_idx <- attr(pmwgs$data, "shared_ll_idx")
  # if(stage == "sample"){
  components <- rep(1, length(components))
  shared_ll_idx <- rep(1, length(shared_ll_idx))
  # }
  # Display stage to screen

  alphaStar=-qnorm(p_accept/2) #Idk about this one
  n0=round(5/(p_accept*(1-p_accept))) #Also not questioning this math for now

  epsilon <- fix_epsilon(pmwgs, epsilon, force_prev_epsilon, components)
  if(length(particles == 1)){
    particles <- rep(particles, max(pmwgs$n_subjects, 2)) # kluge to keep it as a vector
  }
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  eff_mu <- attr(pmwgs, "eff_mu")
  eff_var <- attr(pmwgs, "eff_var")
  chains_cov <- attr(pmwgs, "chains_cov")
  mix <- set_mix(stage)
  if (verboseProgress) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx

  data <- pmwgs$data
  subjects <- pmwgs$subjects
  unq_components <- unique(components)
  block_idx <- block_variance_idx(components[!grouped])
  # Main iteration loop
  for (i in 1:iter) {
    if (verboseProgress) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }
    j <- start_iter + i

    # Gibbs step
    pars <- variant_funs$gibbs_step(pmwgs, pmwgs$samples$samples_fixed[,j-1], pmwgs$samples$samples_random[,,j-1])

    # Particle step
    proposals_random <- mclapply(X=1:pmwgs$n_subjects,FUN = new_particle_random, data, particles_random, pars, eff_mu,
                       eff_var, mix, pmwgs$ll_func, epsilon, subjects, components,
                       prev_ll = pmwgs$samples$subj_ll[,j-1], stage, chains_cov,
                       block_idx, shared_ll_idx, pmwgs$random_names,
                       mc.cores =n_cores)
    proposals_random <- array(unlist(proposals), dim = c(pmwgs$n_random + 2, pmwgs$n_subjects))

    proposals_fixed <- new_particle_fixed(data, particles_fixed, pars, eff_mu,
                              eff_var, mix, pmwgs$ll_func, epsilon, subjects, components,
                              prev_ll = pmwgs$samples$subj_ll[,j-1], stage, chains_cov,
                              block_idx, shared_ll_idx, pmwgs$random_names,
                              mc.cores =n_cores)
    proposals_fixed <- array(unlist(proposals), dim = c(pmwgs$n_random + 2, pmwgs$n_subjects))

    #Fill samples
    pmwgs$samples <- variant_funs$fill_samples(samples = pmwgs$samples, group_level = pars,
                                               proposals = proposals, epsilon = rowMeans(epsilon), j = j, n_pars = sum(!pmwgs$grouped))

    # Update epsilon
    if(!is.null(p_accept)){
      if(j > n0){
        for(component in unq_components){
          idx <- components[!grouped] == component
          acc <-  pmwgs$samples$alpha[max(which(idx)),,j] != pmwgs$samples$alpha[max(which(idx)),,(j-1)]
          epsilon[,component] <-update.epsilon(epsilon[,component]^2, acc, p_accept, j, sum(idx), alphaStar)
        }
      }
    }
  }
  if (verboseProgress) close(pb)
  attr(pmwgs, "epsilon") <- epsilon
  return(pmwgs)
}
#
# new_particle_fixed <- function(data, num_particles, prior,
#                                chains_cov, mix_proportion = c(.1, .9), epsilon_grouped,
#                                prev_mu, alpha, par_names, likelihood_func = NULL,
#                                is_grouped, stage,
#                                variant_funs, subjects, n_cores){
#   prior_mu <- prior$theta_mu_mean
#   prior_var <- prior$theta_mu_var
#   if(stage == "preburn"){
#     chains_cov <- prior_var
#   }
#   chains_cov <- chains_cov * epsilon_grouped^2
#   particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
#   prior_particles <- particle_draws(particle_numbers[1], prior_mu, prior_var)
#   ind_particles <- particle_draws(particle_numbers[2], prev_mu, chains_cov)
#   proposals <- rbind(prev_mu, prior_particles, ind_particles)
#   n_subjects <- length(subjects)
#   proposals_list <- vector("list", n_subjects)
#   for(i in 1:n_subjects){
#     proposals_list[[i]] <- bind_alpha(proposals, alpha[,i], num_particles + 1, is_grouped, par_names)
#   }
#   lws <- parallel::mclapply(1:n_subjects, calc_ll_for_group, proposals_list, data, likelihood_func, subjects, mc.cores = n_cores)
#   lw <- Reduce("+", lws)
#   lp <- mvtnorm::dmvnorm(x = proposals, mean = prior_mu, sigma = prior_var, log = TRUE)
#   prop_density <- mvtnorm::dmvnorm(x = proposals, mean = prev_mu, sigma = chains_cov)
#   lm <- log(mix_proportion[1] * exp(lp) + mix_proportion[2] * prop_density)
#   l <- lw + lp - lm
#   weights <- exp(l - max(l))
#   idx <- sample(x = num_particles + 1, size = 1, prob = weights)
#   return(proposals[idx,])
# }
#
#
# new_particle <- function (s, data, num_particles, parameters, eff_mu = NULL,
#                           eff_var = NULL, mix_proportion = c(0.5, 0.5, 0),
#                           likelihood_func = NULL, epsilon = NULL, subjects,
#                           components, prev_ll, stage, chains_cov, variant_funs,
#                           block_idx, shared_ll_idx, grouped_pars, is_grouped,
#                           par_names)
# {
#   group_pars <- variant_funs$get_group_level(parameters, s)
#   unq_components <- unique(components)
#   proposal_out <- numeric(length(group_pars$mu))
#   ll <- 0
#   group_mu <- group_pars$mu
#   group_var <- group_pars$var
#   subj_mu <- parameters$alpha[,s]
#   eff_mu_sub <- eff_mu[,s]
#   num_particles <- num_particles[s]
#   if(stage != "sample"){
#     eff_mu_sub <- subj_mu
#     group_var_subj <- group_var
#   }
#   out_lls <- numeric(length(unq_components))
#   for(i in unq_components){
#     if(stage != "sample"){
#       # In preburn use to proposals:
#       # 1) person distribution with previous particle mean and scaled group variance
#       # 2) group mean and group variance
#       # In burn and sample add covariance proposals with
#       # 3) previous particle mean and scaled chain covariance
#       # In sample rather use:
#       # 1) Previous particle mean and now scaled covariance
#       # 2) Group mean and group variance (same as before)
#       # 3) Conditional mean and covariance
#       # Group distribution, with group mean and variance
#       eff_var_curr <- chains_cov[,,s] * epsilon[s,i]^2
#       if(stage == "preburn"){
#         var_subj <- group_var_subj * epsilon[s,i]^2
#       } else{
#         var_subj <- group_var_subj *  mean(diag(eff_var_curr))/mean(diag(group_var_subj))
#       }
#     } else{
#       eff_var_curr <- eff_var[,,s]
#       var_subj <- chains_cov[,,s] *  epsilon[s,i]^2
#     }
#     idx <- components[!is_grouped] == i
#     # Draw new proposals for this component
#     particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
#     cumuNumbers <- cumsum(particle_numbers) + 1 # Include the particle from b4
#     pop_particles <- particle_draws(particle_numbers[1], group_mu[idx], group_var[idx, idx])
#     ind_particles <- particle_draws(particle_numbers[2], subj_mu[idx], var_subj[idx, idx])
#     if(mix_proportion[3] == 0){
#       eff_particles <- NULL
#     } else{
#       eff_particles <- particle_draws(particle_numbers[3], eff_mu_sub[idx], eff_var_curr[idx, idx ])
#     }
#     # Rejoin new proposals with current MCMC values for other components
#     proposals <- matrix(rep(subj_mu, num_particles + 1), nrow = num_particles + 1, byrow = T)
#     colnames(proposals) <- names(subj_mu)
#     proposals[2:(num_particles+1),idx] <- rbind(pop_particles, ind_particles, eff_particles)
#     ll_proposals <- proposals
#     if(any(is_grouped)){
#       ll_proposals <- update_proposals_grouped(proposals, grouped_pars, is_grouped,
#                                                par_names = par_names)
#     }
#     # Normally we assume that a component contains all the parameters to estimate the individual likelihood of a joint model
#     # Sometimes we may also want to block within a model if it has very high dimensionality
#     shared_idx <- shared_ll_idx[idx][1]
#     is_shared <- shared_idx == shared_ll_idx
#
#     # Calculate likelihoods
#     if(components[length(components)] > 1){
#       lw <- calc_ll_manager(ll_proposals[,is_shared], dadm = data[[which(subjects == s)]], likelihood_func, component = shared_idx)
#     } else{
#       lw <- calc_ll_manager(ll_proposals[,is_shared], dadm = data[[which(subjects == s)]], likelihood_func)
#     }
#     lw_total <- lw + prev_ll[s] - lw[1] # make sure lls from other components are included
#     lp <- mvtnorm::dmvnorm(x = proposals[,idx], mean = group_mu[idx], sigma = group_var[idx,idx], log = TRUE)
#     prop_density <- mvtnorm::dmvnorm(x = proposals[,idx], mean = subj_mu[idx], sigma = var_subj[idx,idx])
#     if (mix_proportion[3] == 0) {
#       eff_density <- 0
#     }
#     else {
#       eff_density <- mvtnorm::dmvnorm(x = proposals[,idx], mean = eff_mu_sub[idx], sigma = eff_var_curr[idx,idx])
#     }
#     lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] * prop_density) + (mix_proportion[3] * eff_density))
#     infnt_idx <- is.infinite(lm)
#     lm[infnt_idx] <- min(lm[!infnt_idx])
#     # Calculate weights and center
#     l <- lw_total + lp - lm
#     weights <- exp(l - max(l))
#     # Do MH step and return everything
#     idx_ll <- sample(x = num_particles+1, size = 1, prob = weights)
#     origin <- min(which(idx_ll <= cumuNumbers))
#     # ll <- ll + lw[idx_ll]
#     out_lls[i] <- lw_total[idx_ll]
#     proposal_out[idx] <- proposals[idx_ll,idx]
#   } # Note that origin only contains last components origin, just used by devs anyway
#   return(list(proposal = proposal_out, ll = mean(out_lls), origin = origin))
# }
#
#
# calc_ll_for_group <- function(s, proposals, data, ll, subjects){
#   lw <- calc_ll_manager(proposals[[s]], dadm = data[[which(subjects == s)]], ll)
# }
#
# bind_alpha <- function(proposals, alpha, num_particles, is_grouped, par_names){
#   tmp <- matrix(0, ncol = length(is_grouped), nrow = num_particles)
#   tmp[,!is_grouped] <- matrix(alpha, ncol = length(alpha), nrow = num_particles, byrow = T)
#   tmp[,is_grouped] <- proposals
#   colnames(tmp) <- par_names
#   return(tmp)
# }
#
# # Utility functions for sampling below ------------------------------------
#
#
# update.epsilon<- function(epsilon2, acc, p, i, d, alpha) {
#   c <- ((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
#   Theta <- log(sqrt(epsilon2))
#   Theta <- Theta+c*(acc-p)/max(200, i/d)
#   return(exp(Theta))
# }
#
#
# numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
#   numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_proportion)
#   if (mix_proportion[3] == 0) {
#     numbers[3] <- 0
#     numbers[2] <- num_particles - numbers[1]
#   } else {
#     numbers[3] <- num_particles - sum(numbers)
#   }
#   return(numbers)
# }
#
#
# particle_draws <- function(n, mu, covar) {
#   if (n <= 0) {
#     return(NULL)
#   }
#   return(mvtnorm::rmvnorm(n, mu, covar))
# }
#
# fix_epsilon <- function(pmwgs, epsilon, force_prev_epsilon, components){
#   if(is.null(epsilon) | force_prev_epsilon){
#     epsilon <- attr(pmwgs, "epsilon")
#     if(is.null(epsilon)){
#       epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
#     }
#   }
#   if(is.matrix(epsilon)){
#     if(ncol(epsilon) != max(components)){
#       epsilon <- matrix(rowMeans(epsilon), nrow = pmwgs$n_subjects, ncol = max(components))
#     }
#   } else if (length(epsilon) == 1 | is.vector(epsilon)){
#     epsilon <- matrix(epsilon, nrow = pmwgs$n_subjects, ncol = max(components))
#   }
#   return(epsilon)
# }
#
# set_epsilon <- function(n_pars) {
#   if (n_pars > 15) {
#     epsilon <- 0.1
#   } else if (n_pars > 10) {
#     epsilon <- 0.3
#   } else {
#     epsilon <- 0.5
#   }
#   return(epsilon)
# }
#
#
# extend_sampler <- function(sampler, n_samples, stage) {
#   # This function takes the sampler and extends it along the intended number of
#   # iterations, to ensure that we're not constantly increasing our sampled object
#   # by 1. Big shout out to the rapply function
#   sampler$samples$stage <- c(sampler$samples$stage, rep(stage, n_samples))
#   if(any(sampler$nuisance)) sampler$sampler_nuis$samples <- rapply(sampler$sampler_nuis$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
#   sampler$samples <- rapply(sampler$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
#   return(sampler)
# }
#
# extend_obj <- function(obj, n_extend){
#   old_dim <- dim(obj)
#   n_dimensions <- length(old_dim)
#   if(is.null(old_dim) | n_dimensions == 1) return(obj)
#   if(n_dimensions == 2){
#     if(nrow(obj) == ncol(obj)){
#       if(nrow(obj) > 1){
#         if(mean(abs(abs(rowSums(obj/max(obj))) - abs(colSums(obj/max(obj))))) < .01) return(obj)
#       }
#     }
#   }
#   new_dim <- c(rep(0, (n_dimensions -1)), n_extend)
#   extended <- array(NA_real_, dim = old_dim +  new_dim, dimnames = dimnames(obj))
#   extended[slice.index(extended,n_dimensions) <= old_dim[n_dimensions]] <- obj
#   return(extended)
# }
#
# sample_store_base <- function(data, par_names, iters = 1, stage = "init", is_nuisance = rep(F, length(par_names)), ...) {
#   subject_ids <- unique(data$subjects)
#   n_pars <- length(par_names)
#   n_subjects <- length(subject_ids)
#   samples <- list(
#     epsilon = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
#     origin = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
#     alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
#     stage = array(stage, iters),
#     subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL))
#   )
# }
#
# block_variance_idx <- function(components){
#   vars_out <- matrix(0, length(components), length(components))
#   for(i in unique(components)){
#     idx <- i == components
#     vars_out[idx,idx] <- NA
#   }
#   return(vars_out == 0)
# }
#
# fill_samples_base <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
#   # Fill samples both group level and random effects
#   samples$theta_mu[, j] <- group_level$tmu
#   samples$theta_var[, , j] <- group_level$tvar
#   if(!is.null(proposals)) samples <- fill_samples_RE(samples, proposals, epsilon,j, n_pars)
#   return(samples)
# }
#
#
#
# fill_samples_RE <- function(samples, proposals, epsilon, j = 1, n_pars, ...){
#   # Only for random effects, separated because group level sometimes differs.
#   if(!is.null(proposals)){
#     samples$alpha[, , j] <- proposals[1:n_pars,]
#     samples$subj_ll[, j] <- proposals[n_pars + 1,]
#     samples$origin[,j] <- proposals[n_pars + 2,]
#     samples$idx <- j
#     samples$epsilon[,j] <- epsilon
#   }
#   return(samples)
# }
#
# set_mix <- function(stage) {
#   if (stage %in% c("burn", "adapt")) {
#     mix <- c(0.15, 0.15, 0.7)
#   } else if(stage == "sample"){
#     mix <- c(0.1, 0.3, 0.6)
#   } else{ #Preburn stage
#     mix <- c(0.5, 0.5, 0)
#   }
#   return(mix)
# }
#
# condMVN <- function (mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE)
# {
#   if (missing(dependent.ind))
#     return("You must specify the indices of dependent random variables in `dependent.ind'")
#   if (missing(given.ind) & missing(X.given))
#     return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
#                                                                           dependent.ind])))
#   if (length(given.ind) == 0)
#     return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
#                                                                           dependent.ind])))
#   if (length(X.given) != length(given.ind))
#     stop("lengths of `X.given' and `given.ind' must be same")
#   if (check.sigma) {
#     if (!isSymmetric(sigma))
#       stop("sigma is not a symmetric matrix")
#     eigenvalues <- eigen(sigma, only.values = TRUE)$values
#     if (any(eigenvalues < 1e-08))
#       stop("sigma is not positive-definite")
#   }
#   B <- sigma[dependent.ind, dependent.ind]
#   C <- sigma[dependent.ind, given.ind, drop = FALSE]
#   D <- sigma[given.ind, given.ind]
#   CDinv <- C %*% chol2inv(chol(D))
#   cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
#   cVar <- B - CDinv %*% t(C)
#   list(condMean = cMu, condVar = cVar)
# }
#
#
# # some sort of variants, used to be called via a list called variant_funs
#
# get_variant_funs <- function(type = "standard") {
#   if(type == "standard") {
#     list_fun <- list(# store functions
#       sample_store = sample_store_standard,
#       get_prior = get_prior_standard,
#       add_info = add_info_standard,
#       get_startpoints = get_startpoints_standard,
#       get_group_level = get_group_level_standard,
#       fill_samples = fill_samples_standard,
#       gibbs_step = gibbs_step_standard,
#       filtered_samples = filtered_samples_standard,
#       get_conditionals = get_conditionals_standard,
#       get_all_pars_IS2 = get_all_pars_standard,
#       prior_dist_IS2 = prior_dist_standard,
#       group_dist_IS2 = group_dist_standard
#     )
#   } else if(type == "single"){
#     list_fun <- list(# store functions
#       sample_store = sample_store_base,
#       get_prior = get_prior_single,
#       add_info = add_info_single,
#       get_startpoints = get_startpoints_single,
#       get_group_level = get_group_level_single,
#       fill_samples = fill_samples_RE,
#       gibbs_step = gibbs_step_single,
#       filtered_samples = filtered_samples_single,
#       get_conditionals = get_conditionals_single,
#       get_all_pars_IS2 = get_all_pars_single,
#       prior_dist_IS2 = prior_dist_single,
#       group_dist_IS2 = group_dist_single
#     )
#   } else if(type == "blocked"){
#     list_fun <- list(# store functions
#       sample_store = sample_store_standard,
#       add_info = add_info_blocked,
#       get_prior = get_prior_blocked,
#       get_startpoints = get_startpoints_blocked,
#       get_group_level = get_group_level_standard,
#       fill_samples = fill_samples_standard,
#       gibbs_step = gibbs_step_blocked,
#       filtered_samples = filtered_samples_standard,
#       get_conditionals = get_conditionals_blocked,
#       get_all_pars_IS2 = get_all_pars_blocked,
#       prior_dist_IS2 = prior_dist_blocked,
#       group_dist_IS2 = group_dist_blocked
#     )
#   } else if(type == "diagonal"){
#     list_fun <- list(# store functions
#       sample_store = sample_store_standard,
#       add_info = add_info_diag,
#       get_startpoints = get_startpoints_diag,
#       get_group_level = get_group_level_standard,
#       fill_samples = fill_samples_standard,
#       gibbs_step = gibbs_step_diag,
#       filtered_samples = filtered_samples_standard,
#       get_conditionals = get_conditionals_diag,
#       get_all_pars_IS2 = get_all_pars_diag,
#       prior_dist_IS2 = prior_dist_diag,
#       group_dist_IS2 = group_dist_diag
#     )
#   } else if(type == "factor"){
#     list_fun <- list(# store functions
#       sample_store = sample_store_factor,
#       add_info = add_info_factor,
#       get_startpoints = get_startpoints_factor,
#       get_group_level = get_group_level_standard,
#       fill_samples = fill_samples_factor,
#       gibbs_step = gibbs_step_factor,
#       filtered_samples = filtered_samples_factor,
#       get_conditionals = get_conditionals_factor,
#       get_all_pars_IS2 = get_all_pars_factor,
#       prior_dist_IS2 = prior_dist_factor,
#       group_dist_IS2 = group_dist_factor
#     )
#   } else if(type == "lm"){
#     list_fun <- list(# store functions
#       sample_store = sample_store_lm,
#       add_info = add_info_lm,
#       get_startpoints = get_startpoints_lm,
#       get_group_level = get_group_level_lm,
#       fill_samples = fill_samples_lm,
#       gibbs_step = gibbs_step_lm,
#       filtered_samples = filtered_samples_lm,
#       get_conditionals = get_conditionals_lm,
#       get_all_pars_IS2 = get_all_pars_lm,
#       prior_dist_IS2 = prior_dist_lm,
#       group_dist_IS2 = group_dist_lm
#     )
#   }
#   else if(type == "infnt_factor"){
#     list_fun <- list(# store functions
#       sample_store = sample_store_infnt_factor,
#       add_info = add_info_infnt_factor,
#       get_startpoints = get_startpoints_infnt_factor,
#       get_group_level = get_group_level_standard,
#       fill_samples = fill_samples_infnt_factor,
#       gibbs_step = gibbs_step_infnt_factor,
#       filtered_samples = filtered_samples_infnt_factor,
#       get_conditionals = get_conditionals_infnt_factor,
#       get_all_pars_IS2 = get_all_pars_infnt_factor,
#       prior_dist_IS2 = prior_dist_infnt_factor,
#       group_dist_IS2 = group_dist_infnt_factor
#     )
#   }
#   return(list_fun)
# }
#
# calc_ll_manager <- function(proposals, dadm, ll_func, component = NULL){
#   if(!is.data.frame(dadm)){
#     lls <- log_likelihood_joint(proposals, dadm, component)
#   } else{
#     c_name <- attr(dadm,"model")()$c_name
#     if(is.null(c_name)){ # use the R implementation
#       lls <- apply(proposals,1, ll_func,dadm = dadm)
#     } else{
#       p_types <- attr(dadm,"model")()$p_types
#       designs <- list()
#       for(p in p_types){
#         designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
#       }
#       constants <- attr(dadm, "constants")
#       if(is.null(constants)) constants <- NA
#       n_trials = nrow(dadm)
#       if(c_name == "DDM"){
#         levels(dadm$R) <- c(0,1)
#         pars <- get_pars(proposals[1,],dadm)
#         pars <- cbind(pars, dadm$R)
#         parameter_char <- apply(pars, 1, paste0, collapse = "\t")
#         parameter_factor <- factor(parameter_char, levels = unique(parameter_char))
#         parameter_indices <- split(seq_len(nrow(pars)), f = parameter_factor)
#         names(parameter_indices) <- 1:length(parameter_indices)
#       } else{
#         parameter_indices <- list()
#       }
#       lls <- calc_ll(proposals, dadm, constants = constants, n_trials = n_trials, designs = designs, type = c_name, p_types = p_types,
#                      min_ll = log(1e-10), winner = dadm$winner, expand = attr(dadm, "expand"), group_idx = parameter_indices)
#     }
#   }
#   return(lls)
# }
#
# update_proposals_grouped <- function(proposals, grouped_pars, is_grouped, par_names){
#   proposals_full <- matrix(0, nrow = nrow(proposals), ncol = length(is_grouped))
#   proposals_full[,!is_grouped] <- proposals
#   proposals_full[,is_grouped] <- matrix(grouped_pars, ncol = length(grouped_pars), nrow = nrow(proposals), byrow = T)
#   colnames(proposals_full) <- par_names
#   return(proposals_full)
# }
