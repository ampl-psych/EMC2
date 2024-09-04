add_info_infnt_factor <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  args <- list(...)
  max_factors <- args$n_factors
  if(is.null(max_factors)) max_factors <- 10
  attr(sampler, "max_factors") <- max_factors

  sampler$prior <- get_prior_infnt_factor(prior, sum(!(sampler$nuisance | sampler$grouped)), sample = F, n_factors = max_factors)
  return(sampler)
}

get_prior_infnt_factor <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL,
                                   map = FALSE, n_factors = 10){
  # Checking and default priors
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_p_vector(design, doMap = F))
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
    prior$bs <- rep(.3, n_pars) # rate prior on the error variances
  }
  if(is.null(prior$df)){
    prior$df <- 30 # Shape and rate prior on the global shrinkage
  }
  if(is.null(prior$ad1)){
    prior$ad1 <- 3.5 # Shape prior on first column
  }
  if(is.null(prior$bd1)){
    prior$bd1 <- 1.7 # Rate prior on first column
  }
  if(is.null(prior$ad2)){
    prior$ad2 <- 6 # Multiplicative prior on shape subsequent columns
  }
  if(is.null(prior$bd2)){
    prior$bd2 <- 2.2 # Multiplicative prior on rate of subsequent columns
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- 1/prior$theta_mu_var #Inverse of the matrix
  attr(prior, "type") <- "infnt_factor"
  out <- prior
  if(sample){
    samples <- list()
    par_names <- names(sampled_p_vector(design, doMap = F))
    if(!selection %in% c("mu", "sigma2", "covariance", "alpha", "correlation", "Sigma", "loadings", "residuals")){
      stop("for variant factor, you can only specify the prior on the mean, variance, covariance, loadings, residuals, or the correlation of the parameters")
    }
    if(selection %in% c("mu", "alpha")){
      mu <- t(mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
                               sigma = diag(prior$theta_mu_var)))
      rownames(mu) <- par_names
      if(selection %in% c("mu")){
        samples$theta_mu <- mu
      }
    }
    if(selection %in% c("loadings", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      lambda_tmp <- matrix(0, nrow = n_pars, ncol = n_factors)
      lambda <- array(0, dim = c(n_pars, n_factors, N))
      for(i in 1:N){
        psi <- matrix(1/rgamma(n_factors*n_pars, prior$df/2, prior$df/2), nrow = n_pars, ncol = n_factors)
        delta <- numeric(n_factors)
        delta[1] <- 1/rgamma(1, prior$ad1, prior$bd1)
        if(n_factors > 1){
          delta[2:n_factors] <- 1/rgamma(n_factors -1, prior$ad2, prior$bd2)
        }
        tau <- cumprod(delta)
        for(j in 1:n_factors){
          lambda_tmp[,j] <- rnorm(n_pars, mean = 0, sd = psi[,j]*tau[j])
        }
        lambda[,,i] <- lambda_tmp
      }
      rownames(lambda) <- par_names
      colnames(lambda) <- paste0("F", 1:n_factors)
      if(selection %in% "loadings"){
        samples$theta_lambda <- lambda
      }
    }
    if(selection %in% c("residuals", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      residuals <- t(matrix(rgamma(n_pars*N, shape = prior$as, rate = prior$bs),
                          ncol = n_pars, byrow = T))
      rownames(residuals) <- par_names
      if(selection %in% "residuals"){
        samples$sig_err_inv <- residuals
      }
    }
    if(selection %in% c("sigma2", "covariance", "correlation", "Sigma", "alpha")) {
      vars <- array(NA_real_, dim = c(n_pars, n_pars, N))
      colnames(vars) <- rownames(vars) <- par_names
      for(i in 1:N){
        sigma <- 1/residuals[,i]
        loadings <- lambda[,,i]
        vars[,,i] <- loadings %*% t(loadings) + diag(sigma)
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



sample_store_infnt_factor <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance,
                                      is_grouped, ...) {
  n_factors <- list(...)$n_factors
  if(is.null(n_factors)) n_factors <- 10
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names[!is_grouped], iters, stage)
  par_names <- par_names[!is_nuisance & !is_grouped]
  n_pars <- length(par_names)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    theta_lambda = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, NULL, NULL)),
    theta_psi = array(NA_real_, dim = c(n_pars, n_factors, iters), dimnames = list(par_names, NULL, NULL)),
    theta_delta = array(NA_real_, dim = c(n_factors, iters), dimnames = list(NULL, NULL)),
    theta_sig_err_inv = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL)),
    theta_eta = array(NA_real_, dim = c(n_subjects, n_factors, iters), dimnames = list(subject_ids, NULL, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

get_startpoints_infnt_factor<- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!(pmwgs$nuisance | pmwgs$grouped))
  prior <- pmwgs$prior
  if (is.null(start_mu)) start_mu <- rnorm(n_pars, prior$theta_mu_mean, sd = sqrt(prior$theta_mu_var))
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(n_pars * 3,diag(n_pars))

  hyper <- attributes(pmwgs)

  start_sig_err_inv <- rgamma(n_pars, shape = prior$as, rate = prior$bs)
  start_lambda <- matrix(0, nrow = n_pars, ncol = hyper$max_factors)
  start_eta <- matrix(0, nrow = pmwgs$n_subjects, ncol = hyper$max_factors)
  start_psi <- matrix(rgamma(n_pars * hyper$max_factors,
                             shape = prior$df / 2, rate = prior$df/2),
                      n_pars, hyper$max_factors) # Local shrinkage
  start_delta <- rgamma(hyper$max_factors, shape=c(prior$ad1,
                                                   rep(prior$ad2, hyper$max_factors-1)),
                        scale=c(prior$bd1, rep(prior$bd2, hyper$max_factors-1))) # Global shrinkage
  return(list(tmu = start_mu, tvar = start_var, lambda = start_lambda,
              sig_err_inv = start_sig_err_inv, psi = start_psi,
              eta = start_eta, delta = start_delta))
}

fill_samples_infnt_factor <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$theta_lambda[,,j] <- group_level$lambda
  samples$theta_sig_err_inv[,j] <- group_level$sig_err_inv
  samples$theta_psi[,,j] <- group_level$psi
  samples$theta_eta[,,j] <- group_level$eta
  samples$theta_delta[,j] <- group_level$delta
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}

gibbs_step_infnt_factor <- function(sampler, alpha){
  # Gibbs step for group means with parameter expanded factor analysis from Ghosh & Dunson 2009
  # mu = theta_mu, var = theta_var
  last <- last_sample_infnt_factor(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  # extract previous values (for ease of reading)
  alpha_t <- t(alpha)
  n_subjects <- sampler$n_subjects
  n_pars <- sampler$n_pars - sum(sampler$nuisance) - sum(sampler$grouped)
  max_factors <- hyper$max_factors

  eta <- matrix(last$eta, n_subjects, max_factors)
  psi <- matrix(last$psi, n_pars, max_factors)
  sig_err_inv <- last$sig_err_inv
  lambda <- matrix(last$lambda, n_pars, max_factors)
  mu <- last$mu
  delta <- last$delta

  # Pre-compute
  tauh <- cumprod(delta) # global shrinkage coefficients
  Plam <- psi %*% diag(tauh, nrow = max_factors) # precision of loadings rows
  # Update mu
  # for(i in 1:n_pars){
  #   new_mu <- mu
  #   new_mu[i] <- new_mu[i] + rnorm(1, 0, sd = .05)
  #   ll_new <- sum(mvtnorm::dmvnorm(alpha_t, new_mu, lambda %*% t(lambda) + diag(1/sig_err_inv), log = T))
  #   ll_old <- sum(mvtnorm::dmvnorm(alpha_t, mu, lambda %*% t(lambda) + diag(1/sig_err_inv), log = T))
  #   ratio <- exp(ll_new - ll_old)
  #   if(is.na(ratio)) ratio <- 0 #Could be dividing 0 by 0
  #   #Then sometimes accept worse values
  #   if (runif(1) < ratio){
  #     #Accept!
  #     mu[i] <- new_mu[i]
  #   }
  # }
  mu_sig <- 1/(n_subjects * sig_err_inv + prior$theta_mu_invar)
  mu_mu <- mu_sig * (sig_err_inv * colSums(alpha_t - eta %*% t(lambda)) + prior$theta_mu_invar * prior$theta_mu_mean)
  mu <- rmvnorm(1, mu_mu, diag(mu_sig))
  colnames(mu) <- colnames(alpha_t)
  # calculate mean-centered observations
  alphatilde <- sweep(alpha_t, 2, mu)

  # Update eta
  Lmsg <- diag(sig_err_inv) %*% lambda
  Veta1 = diag(1, nrow = max_factors) + t(Lmsg)%*% lambda
  T_mat <- chol(Veta1)
  qrT <- qr(T_mat)
  R <- qr.R(qrT)
  S <- ginv(R)
  Veta <- tcrossprod(S)
  eta <- alphatilde %*% Lmsg %*% Veta + matrix(rnorm(n_subjects * max_factors),
                                               nrow = n_subjects, ncol = max_factors) %*% t(S)

  # Update lambda
  eta2 <- crossprod(eta)
  for(j in 1:n_pars){
    Qlam <- diag(Plam[j,], nrow = max_factors) + sig_err_inv[j] * eta2
    blam <- sig_err_inv[j] * (t(eta) %*% alphatilde[,j])
    Llam <- t(chol(Qlam))
    zlam <- rnorm(max_factors)
    vlam <- forwardsolve(Llam, blam)
    mlam <- backsolve(t(Llam), vlam)
    ylam <- backsolve(t(Llam), zlam)
    lambda[j,] <- t(ylam + mlam)
  }

  # Update psi_jh
  for(h in 1:max_factors){
    psi[,h] <- rgamma(n_pars, shape = (prior$df + 1)/2, rate = (prior$df + tauh[h] * lambda[,h]^2) / 2)
  }

  # Update delta & tau
  mat <- psi * lambda^2
  ad <- prior$ad1 + 0.5 * n_pars * max_factors
  bd <- prior$bd1 + 0.5 * (1 / delta[1])  * sum(tauh * colSums(mat))
  delta[1] <- rgamma(1, shape = ad, rate = bd)
  tauh <- cumprod(delta)
  if(!is.null(prior$ad2)){
    if(max_factors > 1){
      for(h in 2:max_factors) {
        ad <- prior$ad2 + 0.5 * n_pars * (max_factors - h + 1)
        bd <- prior$bd2 + 0.5 * (1 / delta[h]) * sum(tauh[h:max_factors] * colSums(as.matrix(mat[, h:max_factors])))
        delta[h] <- rgamma(1, shape = ad, rate = bd)
      }
    }
  } else if(max_factors > 1){
    delta[2:max_factors] <- 1
  }



  # Update sig_err
  sig_err_inv <- rgamma(n_pars ,shape = prior$as + .5*n_subjects, rate = prior$bs +
                          .5*colSums((alphatilde - eta %*% t(lambda))^2))

  # Give back
  var <- lambda %*% t(lambda) + diag(1/sig_err_inv)
  return(list(tmu = mu, tvar = var, lambda = lambda, eta = eta,
              sig_err_inv = sig_err_inv, psi = psi, alpha = alpha, delta = delta))
}

last_sample_infnt_factor <- function(store) {
  list(
    mu = store$theta_mu[, store$idx],
    eta = store$theta_eta[,,store$idx],
    delta = store$theta_delta[,store$idx],
    lambda = store$theta_lambda[,,store$idx],
    psi = store$theta_psi[,,store$idx],
    sig_err_inv = store$theta_sig_err_inv[,store$idx]
  )
}

get_conditionals_infnt_factor <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  sig_err <- log(samples$theta_sig_err_inv[idx,])
  eta <- matrix(samples$theta_eta[s,,], nrow = samples$n_factors)
  lambda <- apply(samples$theta_lambda[idx,,, drop = F], 3, as.numeric, samples$n_factors)
  theta_mu <- samples$theta_mu[idx,]
  all_samples <- rbind(samples$alpha[idx, s,],theta_mu, eta, sig_err, lambda)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration],
                                 samples$theta_eta[s,,iteration],
                                 log(samples$theta_sig_err_inv[idx, iteration]),
                                 as.numeric(samples$theta_lambda[idx,, iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

filtered_samples_infnt_factor <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter],
    theta_lambda = sampler$samples$theta_lambda[, , filter, drop = F],
    theta_sig_err_inv = sampler$samples$theta_sig_err_inv[, filter],
    theta_eta = sampler$samples$theta_eta[, , filter, drop = F],
    alpha = sampler$samples$alpha[, , filter],
    n_factors = attr(sampler, "max_factors"),
    iteration = length(filter)
  )
}


get_all_pars_infnt_factor <- function(samples, idx, info){
  stop("no IS2 for infinite factor estimation yet")
}

group_dist_infnt_factor <- function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  stop("no IS2 for infinite factor estimation yet")

}

prior_dist_infnt_factor <- function(parameters, info){
  stop("no IS2 for infinite factor estimation yet")

}
