
add_info_diag <- function(sampler, prior = NULL, ...){
  sampler$prior <- get_prior_diag(prior, sum(!(sampler$nuisance | sampler$grouped)), sample = F)
  return(sampler)
}

#' Prior specification or prior sampling for diagonal estimation
#'
#' To get the default prior for a created design: `get_prior_diag(design = design, sample = FALSE)`
#'
#'
#' For details see Huang, A., & Wand, M. P. (2013). Simple marginally noninformative
#' prior distributions for covariance matrices. *Bayesian Analysis*, 8, 439-452. https://doi.org/10.1214/13-BA815.
#'
#' Note that if `sample = FALSE`, prior$theta_mu_invar (the inverse of the prior covariance matrix on the group-level mean) is returned,
#' which is only used for computational efficiency.
#'
#' @inheritParams get_prior_standard
#' @return A list with a single entry of type of samples from the prior (if `sample = TRUE`) or else a prior object
#' @examples
#' # First define a design for the model
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Now get the default prior
#' prior <- get_prior_diag(design = design_DDMaE, sample = FALSE)
#' # We can change values in the default prior or use `prior`
#' # Then we can get samples from this prior e.g.
#' samples <- get_prior_diag(prior = prior, design = design_DDMaE,
#'   sample = TRUE, selection = "mu")
#'
#' @export
get_prior_diag <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL){
  # Checking and default priors
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
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- diag(rep(1, n_pars))
  }
  if(is.null(prior$v)){
    prior$v <- rep(2, n_pars)
  }
  if(is.null(prior$A)){
    prior$A <- rep(.3, n_pars)
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- ginv(prior$theta_mu_var) #Inverse of the matrix
  attr(prior, "type") <- "diagonal"
  out <- prior
  if(sample){
    par_names <- names(sampled_p_vector(design, doMap = F))
    samples <- list()
    if(selection %in% c("mu", "alpha")){
      mu <- t(mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
                               sigma = prior$theta_mu_var))
      rownames(mu) <- par_names
      if(selection %in% c("mu")){
        samples$theta_mu <- mu
      }
    }
    if(selection %in% c("sigma2", "covariance", "correlation", "Sigma", "alpha")) {
      vars <- array(NA_real_, dim = c(n_pars, n_pars, N))
      colnames(vars) <- rownames(vars) <- par_names
      for(i in 1:N){
        a_half <- 1 / rgamma(n = n_pars,shape = 1/2,
                             rate = 1/(prior$A^2))
        vars[,,i] <- diag(1/ rgamma(n = n_pars, shape = prior$v/2, rate = prior$v/a_half))
      }
      if(selection != "alpha") samples$theta_var <- vars
    }
    if(selection %in% "alpha"){
      samples$alpha <- get_alphas(mu, vars, "alpha")
    }
    out <- samples
  }
  return(out)
}
#   if(sample){
#     out <- list()
#     if(!selection %in% c("mu", "variance", "full_var")){
#       stop("for variant diagonal, you can only specify the prior on the mean or variance parameters")
#     }
#     if(selection == "mu"){
#       samples <- mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
#                                   sigma = prior$theta_mu_var)
#       if(!is.null(design)){
#         colnames(samples) <- par_names <- names(sampled_p_vector(design, doMap = F))
#         if(map){
#           samples <- map_mcmc(samples,design,design$model,include_constants=FALSE)
#         }
#       }
#       out$mu <- samples
#       return(out)
#     } else {
#
#       var <- array(NA_real_, dim = c(N, n_pars))
#       for(i in 1:N){
#         a_half <- 1 / rgamma(n = n_pars,shape = 1/2,
#                              rate = 1/(prior$A^2))
#         var[i,] <- 1 / rgamma(n = n_pars, shape = prior$v/2, rate = prior$v/a_half)
#       }
#       colnames(var) <- names(sampled_p_vector(design, doMap = F))
#       if (selection == "full_var"){
#         out$full_var <- var
#       } else{
#         out$variance <- var
#       }
#       return(out)
#     }
#   }
#   return(prior)
# }


get_startpoints_diag <- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!(pmwgs$nuisance | pmwgs$grouped))
  if (is.null(start_mu)) start_mu <- rnorm(n_pars, mean = pmwgs$prior$theta_mu_mean, sd = sqrt(diag(pmwgs$prior$theta_mu_var)))
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- diag(1/rgamma(n_pars, 10, 5), n_pars) #Bit stupid maybe as startpoint
  start_a_half <- 1 / rgamma(n = n_pars, shape = 2, rate = 1)
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var), a_half = start_a_half))
}

get_conditionals_diag <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  pts2_unwound <- log(apply(samples$theta_var[idx,idx,, drop = F],3,diag))
  all_samples <- rbind(samples$alpha[idx,s,],samples$theta_mu[idx,],pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- var(t(all_samples))
  if(n_pars == 1){
    X.given <- c(samples$theta_mu[idx,iteration], log(samples$theta_var[idx,idx,iteration]))
  } else{
    X.given <- c(samples$theta_mu[idx,iteration], log(diag(samples$theta_var[idx,idx,iteration])))
  }
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = X.given)
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}


gibbs_step_diag <- function(sampler, alpha){
  # Gibbs step for diagonal only
  # Get single iter versions, tmu = theta_mu, tvar = theta_var
  last <- last_sample_standard(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior
  last$tvinv <- diag(last$tvinv)
  n_pars <- sum(!(sampler$nuisance | sampler$grouped))
  alpha <- as.matrix(alpha)
  #Mu
  var_mu = 1.0 / (sampler$n_subjects * last$tvinv + diag(prior$theta_mu_invar))
  mean_mu = var_mu * ((apply(alpha, 1, sum) * last$tvinv + prior$theta_mu_mean * diag(prior$theta_mu_invar)))
  tmu <- rnorm(n_pars, mean_mu, sd = sqrt(var_mu))
  names(tmu) <- sampler$par_names[!(sampler$nuisance | sampler$grouped)]
  tvinv = rgamma(n=n_pars, shape=prior$v/2 + sampler$n_subjects/2, rate=prior$v/last$a_half +
                   rowSums( (alpha-tmu)^2 ) / 2)
  tvar = 1/tvinv
  #Contrary to standard pmwg I use shape, rate for IG()
  a_half <- 1 / rgamma(n = n_pars, shape = (prior$v + 1) / 2,
                       rate = prior$v * tvinv + 1/(prior$A^2))
  return(list(tmu = tmu, tvar = diag(tvar, n_pars), tvinv = diag(tvinv, n_pars), a_half = a_half, alpha = alpha))
}



unwind_diag_IS2 <- function(x,reverse=FALSE, diag = TRUE) {
  if (reverse) {
    if(diag){
      out <- diag(exp(x), nrow = length(x))
    } else{
      out <- exp(x)
    }
  } else {
    out <- log(diag(x))
  }
  return(out)
}

get_all_pars_diag <- function(samples, idx, info){
  n_subjects <- samples$n_subjects
  n_iter = length(samples$samples$stage[idx])
  # Exctract relevant objects
  alpha <- samples$samples$alpha[,,idx]
  theta_mu <- samples$samples$theta_mu[,idx]
  theta_var <- samples$samples$theta_var[,,idx]
  a_half <- log(samples$samples$a_half[,idx])
  theta_var.unwound = log(apply(samples$samples$theta_var[,,idx],3,diag))
  # Set up
  n_params<- samples$n_pars+samples$n_pars+samples$n_pars
  all_samples=array(dim=c(n_subjects,n_params,n_iter))
  mu_tilde=array(dim = c(n_subjects,n_params))
  var_tilde=array(dim = c(n_subjects,n_params,n_params))

  for (j in 1:n_subjects){
    all_samples[j,,] = rbind(alpha[,j,],theta_mu[,],theta_var.unwound[,])
    # calculate the mean for re, mu and sigma
    mu_tilde[j,] =apply(all_samples[j,,],1,mean)
    # calculate the covariance matrix for random effects, mu and sigma
    var_tilde[j,,] = cov(t(all_samples[j,,]))
  }

  for(i in 1:n_subjects){ #RJI_change: this bit makes sure that the sigma tilde is pos def
    if(!corpcor::is.positive.definite(var_tilde[i,,], tol=1e-8)){
      var_tilde[i,,]<-corpcor::make.positive.definite(var_tilde[i,,], tol=1e-6)
    }
  }
  X <- cbind(t(theta_mu),t(theta_var.unwound),t(a_half))
  info$n_params <- n_params
  info$given.ind <- (info$n_randeffect+1):n_params
  info$X.given_ind <- 1:(n_params-info$n_randeffect)
  return(list(X = X, mu_tilde = mu_tilde, var_tilde = var_tilde, info = info))
}

group_dist_diag = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  n_randeffect <- info$n_randeffect
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind_diag_IS2(param.theta.sig.unwound, reverse = TRUE)
  if (sample){
    return(mvtnorm::rmvnorm(n_samples, param.theta.mu,param.theta.sig2))
  }else{
    logw_second<-max(-5000*info$n_randeffect, mvtnorm::dmvnorm(random_effect, param.theta.mu,param.theta.sig2,log=TRUE))
    return(logw_second)
  }
}


prior_dist_diag = function(parameters, info){
  n_randeffect <- info$n_randeffect
  prior <- info$prior
  hyper <- info$hyper
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind_diag_IS2(param.theta.sig.unwound, reverse = TRUE, diag = FALSE)
  param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
  log_prior_mu=mvtnorm::dmvnorm(param.theta.mu, mean = prior$theta_mu_mean, sigma = prior$theta_mu_var, log =TRUE)
  log_prior_sigma = sum(logdinvGamma(param.theta.sig2, shape = prior$v/2, rate = prior$v/param.a))
  log_prior_a = sum(logdinvGamma(param.a,shape = 1/2,rate=1/(prior$A^2)))
  # These are Jacobian corrections for the transformations on these
  logw_den2 <- -sum(log(param.a))
  logw_den3 <- -sum(log(param.theta.sig2))
  return(log_prior_mu + log_prior_sigma + log_prior_a - logw_den3 - logw_den2)
}


# bridge_sampling ---------------------------------------------------------
bridge_group_and_prior_and_jac_diag <- function(proposals_group, proposals_list, info){
  prior <- info$prior
  proposals <- do.call(cbind, proposals_list)
  theta_mu <- proposals_group[,1:info$n_pars]
  theta_var <- proposals_group[,(info$n_pars + 1):(2*info$n_pars)]
  theta_a <- proposals_group[,(2*info$n_pars + 1):(3*info$n_pars)]
  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)
  for(i in 1:n_iter){ # these unfortunately can't be vectorized
    proposals_curr <- matrix(proposals[i,], ncol = info$n_pars, byrow = T)
    group_ll <- sum(dmvnorm(proposals_curr, theta_mu[i,], diag(exp(theta_var[i,]), nrow = info$n_pars), log = T))
    prior_var <- sum(logdinvGamma(exp(theta_var[i,]), shape = prior$v/2, rate = prior$v/exp(theta_a[i,])))
    prior_a <- sum(logdinvGamma(exp(theta_a[i,]), shape = 1/2,rate=1/(prior$A^2)))
    sum_out[i] <- group_ll + prior_var + prior_a
  }
  prior_mu <- dmvnorm(theta_mu, mean = prior$theta_mu_mean, sigma = prior$theta_mu_var, log =T)
  jac_var <- rowSums(theta_var)
  jac_a <- rowSums(theta_a)
  return(sum_out + prior_mu + jac_var + jac_a)
}


bridge_add_info_diag <- function(info, samples){
  info$group_idx <- (samples$n_pars*samples$n_subjects + 1):(samples$n_pars*samples$n_subjects + 3*samples$n_pars)
  return(info)
}

bridge_add_group_diag <- function(all_samples, samples, idx){
  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[,idx]))
  all_samples <- cbind(all_samples, t(log(apply(samples$samples$theta_var[,,idx], 3, diag))))
  all_samples <- cbind(all_samples, t(log(samples$samples$a_half[,idx])))
  return(all_samples)
}




