add_info_infnt_factor <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  args <- list(...)
  max_factors <- args$n_factors
  if(is.null(max_factors)) max_factors <- 10
  attr(sampler, "max_factors") <- max_factors

  sampler$prior <- get_prior_infnt_factor(prior, sampler$n_pars, sample = F, n_factors = max_factors)
  return(sampler)
}

get_prior_infnt_factor <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL,
                                   map = FALSE, n_factors = 10){
  # Checking and default priors
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_pars(design, doMap = F))
  }
  if (!is.null(prior$theta_mu_mean)) {
    n_pars <- length(prior$theta_mu_mean)
  }
  if(is.null(prior$theta_mu_mean)){
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- rep(1, n_pars)
  }
  if(is.null(prior$as)){
    prior$as <- rep(5, n_pars) # shape prior on the error variances
  }
  if(is.null(prior$bs)){
    prior$bs <- rep(.25, n_pars) # rate prior on the error variances
  }
  if(is.null(prior$df)){
    prior$df <- 30 # Shape and rate prior on the global shrinkage
  }
  if(is.null(prior$ad1)){
    prior$ad1 <- 3 # Shape prior on first column
  }
  if(is.null(prior$bd1)){
    prior$bd1 <- 1 # Rate prior on first column
  }
  if(is.null(prior$ad2)){
    prior$ad2 <- 5 # Multiplicative prior on shape subsequent columns
  }
  if(is.null(prior$bd2)){
    prior$bd2 <- 2 # Multiplicative prior on rate of subsequent columns
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- 1/prior$theta_mu_var #Inverse of the matrix
  attr(prior, "type") <- "infnt_factor"
  out <- prior
  if(sample){
    samples <- list()
    par_names <- names(sampled_pars(design, doMap = F))
    if(!selection %in% c("mu", "sigma2", "covariance", "alpha", "correlation", "Sigma", "loadings", "std_loadings", "residuals")){
      stop("for variant infnt_factor, you can only specify the prior on the mean, variance, covariance, loadings, residuals, or the correlation of the parameters")
    }
    if(selection %in% c("mu", "alpha")){
      mu <- t(mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
                               sigma = diag(prior$theta_mu_var, n_pars)))
      rownames(mu) <- par_names
      if(selection %in% c("mu")){
        samples$theta_mu <- mu
      }
    }
    if(selection %in% c("loadings", "std_loadings", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      lambda_tmp <- matrix(0, nrow = n_pars, ncol = n_factors)
      lambda <- array(0, dim = c(n_pars, n_factors, N))
      for(i in 1:N){
        psi <- matrix(rgamma(n_factors*n_pars, prior$df/2, prior$df/2), nrow = n_pars, ncol = n_factors)
        delta <- numeric(n_factors)
        delta[1] <- rgamma(1, prior$ad1, prior$bd1)
        if(n_factors > 1){
          delta[2:n_factors] <- rgamma(n_factors -1, prior$ad2, prior$bd2)
        }
        tau <- cumprod(delta)
        for(j in 1:n_factors){
          lambda_tmp[,j] <- rnorm(n_pars, mean = 0, sd = sqrt(1/(psi[,j]*tau[j])))
        }
        lambda[,,i] <- lambda_tmp
      }
      rownames(lambda) <- par_names
      colnames(lambda) <- paste0("F", 1:n_factors)
      if(selection %in% c("loadings", "std_loadings")){
        samples$lambda <- lambda
      }
    }
    if(selection %in% c("residuals", "std_loadings", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      residuals <- t(matrix(rgamma(n_pars*N, shape = prior$as, rate = prior$bs),
                            ncol = n_pars, byrow = T))
      rownames(residuals) <- par_names
      if(selection %in% c("residuals", "std_loadings")){
        samples$epsilon_inv <- residuals
      }
    }
    if(selection %in% c("sigma2", "covariance", "correlation", "Sigma", "alpha")) {
      vars <- array(NA_real_, dim = c(n_pars, n_pars, N))
      colnames(vars) <- rownames(vars) <- par_names
      for(i in 1:N){
        sigma <- 1/residuals[,i]
        loadings <- matrix(lambda[,,i], n_pars, n_factors)
        vars[,,i] <- loadings %*% t(loadings) + diag(sigma, n_pars)
      }
      if(selection != "alpha") samples$theta_var <- vars
    }
    if(selection %in% "alpha"){
      samples$alpha <- get_alphas(mu, vars, design[[1]]$Ffactors$subjects)
    }
    out <- samples
  }
  return(out)
}



sample_store_infnt_factor <- function(data, par_names, iters = 1, stage = "init", integrate = T, ...) {
  n_factors <- list(...)$n_factors
  if(is.null(n_factors)) n_factors <- 10
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  n_pars <- length(par_names)
  FA_names <- paste0("F", 1:n_factors)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    lambda = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, FA_names, NULL)),
    # Psi here is the precision of the loadings, slightly different to standard factor model
    psi = array(NA_real_, dim = c(n_pars, n_factors, iters), dimnames = list(par_names, NULL, NULL)),
    delta = array(NA_real_, dim = c(n_factors, iters), dimnames = list(NULL, NULL)),
    epsilon_inv = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL)),
    eta = array(NA_real_, dim = c(n_subjects, n_factors, iters), dimnames = list(subject_ids, FA_names, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

filtered_samples_infnt_factor <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter, drop = FALSE],
    lambda = sampler$samples$lambda[, , filter, drop = F],
    epsilon_inv = sampler$samples$epsilon_inv[, filter, drop = FALSE],
    eta = sampler$samples$eta[, , filter, drop = F],
    alpha = sampler$samples$alpha[, , filter, drop = FALSE],
    n_factors = attr(sampler, "max_factors"),
    iteration = length(filter)
  )
}
