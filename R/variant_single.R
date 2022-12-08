add_info_single <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, sampler$n_pars), theta_mu_var = diag(rep(1, sampler$n_pars)))
  }
  sampler$prior <- prior
  return(sampler)
}

get_startpoints_single <- function(pmwgs, start_mu, start_var){
  if (is.null(start_mu)) start_mu <- pmwgs$prior$theta_mu_mean
  # If no starting point for var just sample some
  if (is.null(start_var)) start_var <- pmwgs$prior$theta_mu_var
  return(list(tmu = start_mu, tvar = start_var))
}

get_group_level_single <- function(parameters, s){
  # "group level" is just prior here
  mu <- parameters$tmu
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

gibbs_step_single <- function(sampler, alpha){
  alpha_out <- matrix(alpha, nrow = sampler$n_pars, ncol = sampler$n_subjects)
  rownames(alpha_out) <- sampler$par_names
  return(list(tmu = sampler$prior$theta_mu_mean,tvar = sampler$prior$theta_mu_var,alpha = alpha_out))
}

get_conditionals_single <- function(s, samples, n_pars, iteration = NULL){
  stop("no adaptation/sample stage for single subject estimation")
}

filtered_samples_single <- function(sampler, filter){
  stop("no adaptation/sample stage for single subject estimation")
}

get_all_pars_single <- function(samples, idx, info){
  stop("no IS2 for single subject estimation yet")
}

group_dist_single <- function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  stop("no IS2 for single subject estimation yet")

}

prior_dist_single <- function(parameters, info){
  stop("no IS2 for single subject estimation yet")

}

