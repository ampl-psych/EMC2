add_info_infnt_factor <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  args <- list(...)
  max_factors <- args$max_factors
  if(is.null(max_factors)) max_factors <- 10
  n_pars <- sampler$n_pars
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, n_pars),
                  theta_mu_var = rep(1, n_pars))
  }
  # Things I save rather than re-compute inside the loops.
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- 1/prior$theta_mu_var

  #Hyper parameters
  # for residual precision
  attr(sampler, "as") <- 1
  attr(sampler, "bs") <- .3
  # for t_{ij}
  attr(sampler, "df") <- 3
  # for delta_1, ad1 must be > 2
  attr(sampler, "ad1") <- 2.1
  attr(sampler, "bd1") <- 1
  # for delta_h, h >=2, ad2 must be > 3
  attr(sampler, "ad2") <- 3.1
  attr(sampler, "bd2") <- 2

  attr(sampler, "max_factors") <- max_factors

  sampler$prior <- prior
  return(sampler)
}

sample_store_infnt_factor <- function(data, par_names, iters = 1, stage = "init", integrate = T, ...) {
  n_factors <- list(...)$max_factors
  if(is.null(n_factors)) n_factors <- 10
  subject_ids <- unique(data$subjects)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
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
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$n_pars, sd = 1)
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))

  hyper <- attributes(pmwgs)

  start_sig_err_inv <- rgamma(pmwgs$n_pars, shape = hyper$as, rate = hyper$bs)
  start_lambda <- matrix(0, nrow = pmwgs$n_pars, ncol = hyper$max_factors)
  start_eta <- matrix(0, nrow = pmwgs$n_subjects, ncol = hyper$max_factors)
  start_psi <- matrix(rgamma(pmwgs$n_pars * hyper$max_factors,
                             shape = hyper$df / 2, rate = hyper$df/2),
                      pmwgs$n_pars, hyper$max_factors) # Local shrinkage
  start_delta <- rgamma(hyper$max_factors, shape=c(hyper$ad1,
                                                   rep(hyper$ad2, hyper$max_factors-1)),
                        scale=c(hyper$bd1, rep(hyper$bd2, hyper$max_factors-1))) # Global shrinkage
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
  n_pars <- sampler$n_pars
  max_factors <- hyper$max_factors

  eta <- matrix(last$eta, n_subjects, max_factors)
  psi <- matrix(last$psi, n_pars, max_factors)
  sig_err_inv <- last$sig_err_inv
  lambda <- matrix(last$lambda, n_pars, max_factors)
  mu <- last$mu
  delta <- last$delta

  # Pre-compute
  tauh <- cumprod(delta) # global shrinkage coefficients
  Plam <- matvec(psi, tauh) # precision of loadings rows


  #Update mu
  mu_sig <- 1/(n_subjects * sig_err_inv + prior$theta_mu_invar)
  mu_mu <- mu_sig * (sig_err_inv * colSums(alpha_t - eta %*% t(lambda)) + prior$theta_mu_invar * prior$theta_mu_mean)
  mu <- rmvnorm(1, mu_mu, diag(mu_sig))
  colnames(mu) <- colnames(alpha_t)
  # calculate mean-centered observations
  alphatilde <- sweep(alpha_t, 2, mu)

  # Update eta
  Lmsg <- vecmat(sig_err_inv, lambda)
  Veta1 = diag(max_factors) + t(Lmsg)%*% lambda
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
    Qlam <- diag(Plam[j,]) + sig_err_inv[j] * eta2
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
    psi[,h] <- rgamma(n_pars, shape = (hyper$df + 1)/2, rate = (hyper$df + tauh[h] * lambda[,h]^2) / 2)
  }

  # Update delta & tau
  mat <- psi * lambda^2
  ad <- hyper$ad1 + 0.5 * n_pars * max_factors
  bd <- hyper$bd1 + 0.5 * (1 / delta[1])  * sum(tauh * colSums(mat))
  delta[1] <- rgamma(1, shape = ad, rate = bd)
  tauh <- cumprod(delta)

  for(h in 2:max_factors) {
    ad <- hyper$ad2 + 0.5 * n_pars * (max_factors - h + 1)
    bd <- hyper$bd2 + 0.5 * (1 / delta[h]) * sum(tauh[h:max_factors] * colSums(as.matrix(mat[, h:max_factors])))
    delta[h] <- rgamma(1, shape = ad, rate = bd)
  }

  # Update sig_err
  sig_err_inv <- rgamma(n_pars ,shape = hyper$as + .5*n_subjects, rate = hyper$bs +
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
  sig_err <- log(samples$theta_sig_err_inv[idx])
  eta <- matrix(samples$theta_eta[s,,], nrow = samples$n_factors)
  lambda <- apply(samples$theta_lambda[idx,,], 3, as.numeric, samples$n_factors)
  theta_mu <- samples$theta_mu[idx,]
  all_samples <- rbind(samples$alpha[idx, s,],theta_mu, eta, sig_err, lambda)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(theta_mu[idx,iteration],
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

