add_info_single <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  sampler$prior <- get_prior_single(prior, sampler$n_pars, sample = F)
  return(sampler)
}

#' Prior specification or prior sampling for single subject estimation.
#'
#' With this type of estimation, we assume that one, or multiple, subjects are
#' estimated without any hierarchical constraint. We need to specify a prior
#' with a multivariate normal from, by providing specifying prior$theta_mu_mean
#' a vector with an entry for each parameter, and a prior covariance matrix
#' prior$theta_mu_var, with default list(theta_mu_mean = rep(0, n_pars),
#' theta_mu_var = diag(rep(0, n_pars))).
#'
#' @param prior A named list containing the prior mean (theta_mu_mean) and
#' variance (theta_mu_var). Default prior created if NULL
#' @param n_pars Argument used by the sampler, best left NULL. In user case inferred from the design
#' @param sample Whether to sample from the prior. Default is TRUE
#' @param map Boolean, default TRUE reverses malformation used by model to make
#' sampled parameters unbounded
#' @param N How many samples to draw from the prior, default 1e5
#' @param design The design obtained from `make_design`, required when map = TRUE
#' @param type  FIX ME
#'
#' @return A list with single entry named "alpha" of samples from the prior (if sample = TRUE) or else a prior object
#' @export

get_prior_single <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5,
                             type = "alpha", design = NULL, map = FALSE){
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(attr(design, "p_vector"))
  }
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- setNames(rep(0, n_pars),names(attr(design, "p_vector")))
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- diag(rep(1, n_pars))
  }
  if(sample){
    if(type != "alpha") stop("for variant single, only alpha can be specified")
    samples <- mvtnorm::rmvnorm(N, prior$theta_mu_mean, prior$theta_mu_var)
    if (map) {
      proot <- unlist(lapply(strsplit(colnames(samples),"_"),function(x)x[[1]]))
      isin <- proot %in% names(design$model()$p_types)
      fullnames <- colnames(samples)[isin]
      colnames(samples)[isin] <- proot
      samples[,isin] <- design$model()$Ntransform(samples[,isin])
      colnames(samples)[isin] <- fullnames
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

