sample_store_diag_gamma <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, is_grouped, ...) {
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names[!is_grouped], iters, stage)
  par_names <- par_names[!is_nuisance & !is_grouped]
  n_pars <- length(par_names)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

fill_samples_diag_gamma <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$last_theta_var_inv <- group_level$tvinv
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}


add_info_diag_gamma <- function(sampler, prior = NULL, ...){
  sampler$prior <- get_prior_diag_gamma(prior, sum(!(sampler$nuisance | sampler$grouped)), sample = F)
  return(sampler)
}

get_prior_diag_gamma <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL){
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
  if(is.null(prior$shape)){
    prior$shape <- rep(2, n_pars)
  }
  if(is.null(prior$rate)){
    prior$rate <- rep(.3, n_pars)
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- ginv(prior$theta_mu_var) #Inverse of the matrix
  attr(prior, "type") <- "diagonal-gamma"
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
        vars[,,i] <- diag(1/ rgamma(n = n_pars, shape = prior$shape, rate = prior$rate))
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

get_startpoints_diag_gamma <- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!(pmwgs$nuisance | pmwgs$grouped))
  if (is.null(start_mu)) start_mu <- rnorm(n_pars, mean = pmwgs$prior$theta_mu_mean, sd = sqrt(diag(pmwgs$prior$theta_mu_var)))
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- diag(1/rgamma(n_pars, 10, 5), n_pars) #Bit stupid maybe as startpoint
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var)))
}

last_sample_diag_gamma <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],
    tvar = store$theta_var[, , store$idx],
    tvinv = store$last_theta_var_inv
  )
}


gibbs_step_diag_gamma <- function(sampler, alpha){
  # Gibbs step for diagonal only
  # Get single iter versions, tmu = theta_mu, tvar = theta_var
  last <- last_sample_diag_gamma(sampler$samples)
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


  # InvGamma alternative (probably inferior) prior
  shape = prior$shape + sampler$n_subjects / 2
  rate = prior$rate + rowSums( (alpha-tmu)^2 ) / 2
  tvinv = rgamma(n=sampler$n_pars, shape=shape, rate=rate)
  tvar = 1/tvinv
  return(list(tmu = tmu, tvar = diag(tvar, n_pars), tvinv = diag(tvinv, n_pars), alpha = alpha))
}


