add_info_single <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  sampler$prior <- get_prior_single(prior, sampler$n_pars)
  return(sampler)
}

get_prior_single <- function(prior = NULL, n_pars = NULL, par_names = NULL, plot = F, plot_all_pars = F){
  if(is.null(prior)){
    prior <- list()
  }
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- diag(rep(1, n_pars))
  }
  if(plot){
    par(mfrow = c(3,3))
    N <- 1e6
    titles <- paste0("prior - ", par_names)
    for(i in 1:n_pars){
      data <- rnorm(N, prior$theta_mu_mean[i], prior$theta_mu_var[i,i])
      hist(data, prob = T, main = titles[i])
      lines(density(data))
    }
  }
  return(prior)
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
  n_pars <- sampler$n_pars-sum(sampler$nuisance) - sum(sampler$grouped)
  alpha_out <- matrix(alpha, nrow = n_pars, ncol = sampler$n_subjects)
  rownames(alpha_out) <- sampler$par_names[!(sampler$nuisance | sampler$grouped)]
  return(list(tmu = sampler$prior$theta_mu_mean,tvar = sampler$prior$theta_mu_var,alpha = alpha_out))
}

get_conditionals_single <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iters <- dim(samples$alpha)[3]
  iter_idx <- sample(1:iters, min(iters, 250))
  if(is.null(idx)) idx <- 1:n_pars
  all_samples <- samples$alpha[idx,s,iter_idx]
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- var(t(all_samples))
  return(list(eff_mu = mu_tilde, eff_var = var_tilde))
}

filtered_samples_single <- function(sampler, filter){
  out <- list(
    alpha = sampler$samples$alpha[, , filter],
    iteration = length(filter)
  )
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

