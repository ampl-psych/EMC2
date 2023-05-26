
add_info_diag <- function(sampler, prior = NULL, ...){
  n_pars <- sum(!(sampler$nuisance | sampler$grouped))
  # Checking and default priors
  if(is.null(prior)){
    prior <- list()
  }
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- rep(1, n_pars)
  }
  if(is.null(prior$v)){
    prior$v <- 2
  }
  if(is.null(prior$A)){
    prior$A <- rep(1, n_pars)
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- 1/prior$theta_mu_var #Inverse of the prior
  sampler$prior <- prior
  return(sampler)
}

get_startpoints_diag <- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!(pmwgs$nuisance | pmwgs$grouped))
  if (is.null(start_mu)) start_mu <- rnorm(n_pars, mean = pmwgs$prior$theta_mu_mean, sd = sqrt(pmwgs$prior$theta_mu_var))
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- diag(1/rgamma(n_pars, 10, 5)) #Bit stupid maybe as startpoint
  start_a_half <- 1 / rgamma(n = n_pars, shape = 2, rate = 1)
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var), a_half = start_a_half))
}

get_conditionals_diag <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  pts2_unwound <- log(apply(samples$theta_var[idx,idx,],3,diag))
  all_samples <- rbind(samples$alpha[idx,s,],samples$theta_mu[idx,],pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- var(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration], log(diag(samples$theta_var[idx,idx,iteration]))))
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
  #Mu
  var_mu = 1.0 / (sampler$n_subjects * last$tvinv + prior$theta_mu_invar)
  mean_mu = var_mu * ((apply(alpha, 1, sum) * last$tvinv + prior$theta_mu_mean * prior$theta_mu_invar))
  tmu <- rnorm(n_pars, mean_mu, sd = sqrt(var_mu))
  names(tmu) <- sampler$par_names[!(sampler$nuisance | sampler$grouped)]

  # if(!is.null(hyper$std_shape)){
  #   # InvGamma alternative (probably inferior) prior
  #   shape = hyper$std_shape + sampler$n_subjects / 2
  #   rate = hyper$std_rate + rowSums( (alpha-tmu)^2 ) / 2
  #   tvinv = rgamma(n=sampler$n_pars, shape=shape, rate=rate)
  #   tvar = 1/tvinv
  #   a_half <- NULL
  # } else {
  tvinv = rgamma(n=n_pars, shape=prior$v/2 + sampler$n_subjects/2, rate=prior$v/last$a_half +
                   rowSums( (alpha-tmu)^2 ) / 2)
  tvar = 1/tvinv
  #Contrary to standard pmwg I use shape, rate for IG()
  a_half <- 1 / rgamma(n = n_pars, shape = (prior$v + n_pars) / 2,
                       rate = prior$v * tvinv + 1/(prior$A^2))
  return(list(tmu = tmu, tvar = diag(tvar), tvinv = diag(tvinv), a_half = a_half, alpha = alpha))
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
