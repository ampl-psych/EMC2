add_info_single <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  sampler$prior <- get_prior_single(prior, sampler$n_pars, sample = F)
  return(sampler)
}

#' Prior specification or prior sampling for single subject estimation
#'
#' With this type of estimation, one or multiple subjects are
#' estimated independently, without any hierarchical constraint.
#'
#' To specify a (multivariate normal) prior, `prior$theta_mu_mean` and `prior$theta_mu_var` an entry
#' is needed for each parameter.
#'
#'
#' @param prior A named list containing the prior mean (`theta_mu_mean`) and
#' variance (`theta_mu_var`). If `NULL`, the default prior is used.
#' @param n_pars Often inferred from the design, but if `design = NULL`, `n_pars`
#' will be used to determine the size of prior.
#' @param sample Boolean, defaults to `TRUE`, sample from the prior or simply
#' return the prior specifications?
#' @param map Boolean, defaults to `TRUE`. If `sample = TRUE`, the implied prior is sampled.
#' This includes back-transformations for naturally bounded parameters such as
#' the non-decision time and an inverse mapping from the design matrix back to the
#' cells of the design. If `FALSE`, the transformed, unmapped, parameters are used.
#' Note that `map` does not affect the prior used in the sampling process.
#' @param N How many samples to draw from the prior, the default is 1e5
#' @param design The design obtained from `design()`, required when `map = TRUE`
#' @param selection  Character. If `sample = TRUE`, what prior to sample from. Options: `"alpha"`.
#' @return A list with a single entry named `"alpha"` and samples from the prior (if `sample = TRUE`) or else a prior object
#' @examples
#' # First define a design for the model
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Now get the default prior
#' prior <- get_prior_single(design = design_DDMaE, sample = FALSE)
#' # We can change values in the default prior or use `prior`
#' # Then we can get samples from this prior e.g.
#' samples <- get_prior_single(prior = prior, design = design_DDMaE,
#'   sample = TRUE, selection = "alpha")
#'
#' @export

get_prior_single <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5,
                             selection = "alpha", design = NULL, map = FALSE){
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_p_vector(design, doMap = F))
  }
  if (!is.null(prior$theta_mu_mean)) {
    n_pars <- length(prior$theta_mu_mean)
  }
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- setNames(rep(0, n_pars),names(sampled_p_vector(design, doMap = F)))
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- diag(rep(1, n_pars))
  }
  attr(prior, "type") <- "single"
  out <- prior
  if(sample){
    out <- list()
    par_names <- names(sampled_p_vector(design, doMap = F))
    if(selection != "alpha") stop("for variant single, only alpha can be specified")
    out$alpha <- t(mvtnorm::rmvnorm(N, prior$theta_mu_mean, prior$theta_mu_var))
    out$alpha <- array(out$alpha, dim = c(length(par_names), 1, N))
    rownames(out$alpha) <- par_names
  }
  return(out)
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
#
# get_all_pars_single <- function(samples, idx, info){
#   n_subjects <- samples$n_subjects
#   n_iter = length(samples$samples$stage[idx])
#   # Exctract relevant objects
#   alpha <- samples$samples$alpha[,,idx, drop = F]
#   # Set up
#   n_params<- samples$n_pars
#   mu_tilde=array(dim = c(n_subjects,n_params))
#   var_tilde=array(dim = c(n_subjects,n_params,n_params))
#   for (j in 1:n_subjects){
#     # calculate the mean for re, mu and sigma
#     mu_tilde[j,] <- rowMeans(alpha[,j,])
#     # calculate the covariance matrix for random effects, mu and sigma
#     var_tilde[j,,] <- cov(t(alpha[,j,]))
#   }
#
#   for(i in 1:n_subjects){ #RJI_change: this bit makes sure that the sigma tilde is pos def
#     if(!corpcor::is.positive.definite(var_tilde[i,,], tol=1e-8)){
#       var_tilde[i,,]<-corpcor::make.positive.definite(var_tilde[i,,], tol=1e-6)
#     }
#   }
#   info$n_params <- n_params
#   return(list(mu_tilde = mu_tilde, var_tilde = var_tilde, info = info))
# }
#
# group_dist_single <- function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
#   # This is for the single case actually the prior distribution.
#   if (sample){
#     return(rmvnorm(n_samples, info$prior$theta_mu_mean, info$prior$theta_mu_var))
#   }else{
#     logw_second<-max(-5000*info$n_randeffect, dmvnorm(random_effect, info$prior$theta_mu_mean,info$prior$theta_mu_var,log=TRUE))
#     return(logw_second)
#   }
# }
#
# prior_dist_single <- function(parameters, info){
#   # This is quite confusing, but now the actual prior dist is the group dist.
#   # Just here for compatability with the other types.
#   return(0)
# }
#

# bridge sampling ---------------------------------------------------------

bridge_add_group_single <- function(all_samples, samples, idx){
  return(all_samples)
}

bridge_add_info_single <- function(info, samples){
  return(info)
}


bridge_group_and_prior_and_jac_single <- function(proposals_group, proposals_list, info){
  proposals <- do.call(rbind, proposals_list)
  prior_mu <- dmvnorm(proposals, mean = info$prior$theta_mu_mean, sigma = info$prior$theta_mu_var, log =T)
  n_iter <- nrow(proposals_list[[1]])
  prior_mu <- matrix(prior_mu, ncol = info$n_subjects)
  return(rowSums(prior_mu)) # Output is of length n_iter
}

group__IC_single <- function(emc, stage="sample",filter=0){
  return(list(mean_ll = NULL, Dmean = NULL,
              minD = NULL))
}
