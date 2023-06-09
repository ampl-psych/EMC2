add_info_single <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  sampler$prior <- get_prior_single(prior, sampler$n_pars, sample = F)
  return(sampler)
}

#' Prior specification for single subject sampling.
#'
#' With this type of sampling, we assume that one, or multiple, subjects are estimated without any hierarchical link.
#' We need to specify a prior for each parameter. For now, the package assumes a multivariate normal prior on each parameter.
#' Thus you need to specify prior$theta_mu_mean a vector with an entry for each parameter. Furthermore you need a prior covariance matrix prior$theta_mu_var.
#' Default is: prior <- list(theta_mu_mean = rep(0, n_pars), theta_mu_var = diag(rep(0, n_pars)))
#'
#' Specific constraints, such as parameter x needs to be larger than 0, are enforced through transformations performed by the models.
#' Thus you need to be mindful that even though the prior might be normally distributed, the transformed prior might look very differently.
#' Set `sample = TRUE` and pass your design object to `design` to sample from the prior. You can then plot the samples using `plot_prior`
#'
#' @param prior A list of the prior containing the mean (theta_mu_mean) and variance (theta_mu_var).
#' @param n_pars Argument used by the sampler, best left NULL. In user case inferred from the design
#' @param sample Whether to sample from the prior. Default is TRUE
#' @param N How many samples to draw from the prior
#' @param type String specifying on what you want the prior. For single subject only "alpha" is a valid option.
#' @param design The design obtained from `make_design`
#'
#' @return A list of samples from the prior
#' @export
#'
#' @examples
get_prior_single <- function(prior = NULL, n_pars = NULL, sample = T, N = 1e5, type = "alpha", design = NULL){
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(attr(design, "p_vector"))
  }
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- diag(rep(1, n_pars))
  }
  if(sample){
    if(type != "alpha") stop("for variant single, only mu can be specified")
    samples <- mvtnorm::rmvnorm(N, prior$theta_mu_mean, prior$theta_mu_var)
    if(!is.null(design)){
      colnames(samples) <- names(attr(design, "p_vector"))
      samples[,colnames(samples) %in% design$model()$p_types] <- design$model()$Ntransform(samples[,colnames(samples) %in% design$model()$p_types])
    }
    return(list(alpha = samples))
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
    alpha = sampler$samples$alpha[, , filter, drop = F],
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

