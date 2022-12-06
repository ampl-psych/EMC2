
sample_store_standard <- function(data, par_names, iters = 1, stage = "init", integrate = T, ...) {
  subject_ids <- unique(data$subjects)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    a_half = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

add_info_standard <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, sampler$n_pars), theta_mu_var = diag(rep(1, sampler$n_pars)))
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- MASS::ginv(prior$theta_mu_var) #Inverse of the matrix

  #Hyper parameters
  attr(sampler, "v_half") <- 2
  attr(sampler, "A_half") <- 1
  sampler$prior <- prior
  return(sampler)
}

get_startpoints_standard <- function(pmwgs, start_mu, start_var){
  if (is.null(start_mu)) start_mu <- mvtnorm::rmvnorm(1, mean = pmwgs$prior$theta_mu_mean, sigma = pmwgs$prior$theta_mu_var)
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- MCMCpack::riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
  start_a_half <- 1 / rgamma(n = pmwgs$n_pars, shape = 2, rate = 1)
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var), a_half = start_a_half))
}

get_group_level_standard <- function(parameters, s){
  # This function is modified for other versions
  mu <- parameters$tmu
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

fill_samples_standard <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$a_half[, j] <- group_level$a_half
  samples$last_theta_var_inv <- group_level$tvinv
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}



gibbs_step_standard <- function(sampler, alpha){
  # Gibbs step for group means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  last <- last_sample_standard(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  # Here mu is group mean, so we are getting mean and variance
  var_mu <- MASS::ginv(sampler$n_subjects * last$tvinv + prior$theta_mu_invar)
  mean_mu <- as.vector(var_mu %*% (last$tvinv %*% apply(alpha, 1, sum) +
                                     prior$theta_mu_invar %*% prior$theta_mu_mean))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  tmu <- mvtnorm::rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names

  # New values for group var
  theta_temp <- alpha - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  if(!is.null(hyper$std_df)){
    B_half <- hyper$std_scale * diag(1, nrow = sampler$n_pars) + cov_temp # nolint
    tvar <- MCMCpack::riwish(hyper$std_df + sampler$n_subjects, B_half) # New sample for group variance
    tvinv <- MASS::ginv(tvar)
    # Sample new mixing weights.
    a_half <- NULL
  } else{
    B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
    tvar <- MCMCpack::riwish(hyper$v_half + sampler$n_pars - 1 + sampler$n_subjects, B_half) # New sample for group variance
    tvinv <- MASS::ginv(tvar)

    # Sample new mixing weights.
    a_half <- 1 / rgamma(n = sampler$n_pars,shape = (hyper$v_half + sampler$n_pars) / 2,
                         rate = hyper$v_half * diag(tvinv) + 1/(hyper$A_half^2))
  }
  return(list(tmu = tmu,tvar = tvar,tvinv = tvinv,a_half = a_half,alpha = alpha))
}

get_conditionals_standard <- function(s, samples, n_pars, iteration = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  pts2_unwound <- apply(samples$theta_var,3,unwind)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- stats::cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[,iteration], unwind(samples$theta_var[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  y[lower.tri(y, diag = TRUE)]
}

last_sample_standard <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],
    tvar = store$theta_var[, , store$idx],
    tvinv = store$last_theta_var_inv,
    a_half = store$a_half[, store$idx]
  )
}

filtered_samples_standard <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter],
    theta_var = sampler$samples$theta_var[, , filter],
    alpha = sampler$samples$alpha[, , filter],
    iteration = length(filter)
  )
}

# Used by IS2
get_all_pars_standard <- function(samples, idx, info){
  n_subjects <- samples$n_subjects
  n_iter = length(samples$samples$stage[idx])
  # Exctract relevant objects
  alpha <- samples$samples$alpha[,,idx]
  theta_mu <- samples$samples$theta_mu[,idx]
  theta_var <- samples$samples$theta_var[,,idx]
  a_half <- log(samples$samples$a_half[,idx])
  theta_var.unwound = apply(theta_var,3,unwind_IS2)
  # Set up
  n_params<- samples$n_pars+nrow(theta_var.unwound)+samples$n_pars
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

# Used by IS2
robust_diwish <- function (W, v, S) { #RJI_change: this function is to protect against weird proposals in the diwish function, where sometimes matrices weren't pos def
  if (!is.matrix(S)) S <- matrix(S)
  if (!is.matrix(W)) W <- matrix(W)
  p <- nrow(S)
  gammapart <- sum(lgamma((v + 1 - 1:p)/2))
  ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
  if (corpcor::is.positive.definite(W, tol=1e-8)){
    cholW<-base::chol(W)
  }else{
    return(1e-10)
  }
  if (corpcor::is.positive.definite(S, tol=1e-8)){
    cholS <- base::chol(S)
  }else{
    return(1e-10)
  }
  halflogdetS <- sum(log(diag(cholS)))
  halflogdetW <- sum(log(diag(cholW)))
  invW <- chol2inv(cholW)
  exptrace <- sum(S * invW)
  lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * exptrace
  lpdf <- lnum - ldenom
  out <- exp(lpdf)
  if(!is.finite(out)) return(1e-100)
  if(out < 1e-10) return(1e-100)
  return(exp(lpdf))
}

# Used by IS2
unwind_IS2 <- function(x,reverse=FALSE) {

  if (reverse) {
    n=sqrt(2*length(x)+0.25)-0.5 ## Dim of matrix.
    out=array(0,dim=c(n,n))
    out[lower.tri(out,diag=TRUE)]=x
    diag(out)=exp(diag(out))
    out=out%*%t(out)
  } else {
    y=t(base::chol(x))
    diag(y)=log(diag(y))
    out=y[lower.tri(y,diag=TRUE)]
  }
  return(out)
}

# Used by IS2
group_dist_standard <- function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  n_randeffect <- info$n_randeffect
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind_IS2(param.theta.sig.unwound, reverse = TRUE)
  if (sample){
    return(mvtnorm::rmvnorm(n_samples, param.theta.mu,param.theta.sig2))
  }else{
    logw_second<-max(-5000*info$n_randeffect, mvtnorm::dmvnorm(random_effect, param.theta.mu,param.theta.sig2,log=TRUE))
    return(logw_second)
  }
}

# Used by IS2
prior_dist_standard <- function(parameters, info){
  n_randeffect <- info$n_randeffect
  prior <- info$prior
  hyper <- info$hyper
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind_IS2(param.theta.sig.unwound, reverse = TRUE)
  param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
  log_prior_mu=mvtnorm::dmvnorm(param.theta.mu, mean = prior$theta_mu_mean, sigma = prior$theta_mu_var, log =TRUE)
  log_prior_sigma = log(robust_diwish(param.theta.sig2, v=hyper$v_half+ n_randeffect-1, S = 2*hyper$v_half*diag(1/param.a)))
  log_prior_a = sum(logdinvGamma(param.a,shape = 1/2,rate=1/(hyper$A_half^2)))
  # These are Jacobian corrections for the transformations on these
  logw_den2 <- -sum(log(param.a))
  logw_den3 <- -(log(2^n_randeffect)+sum((n_randeffect:1+1)*log(diag(param.theta.sig2))))
  return(log_prior_mu + log_prior_sigma + log_prior_a - logw_den3 - logw_den2)
}

# Used by IS2
logdinvGamma <- function(x, shape, rate){
  alpha <- shape
  beta <- 1/rate
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha +
                                                        1) * log(x) - (beta/x)
  return(pmax(log.density, -50)) #Roughly equal to 1e-22 on real scale
}

