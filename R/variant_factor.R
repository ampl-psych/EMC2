
sample_store_factor <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, ...) {
  n_factors <- list(...)$n_factors
  Lambda_mat <- list(...)$Lambda_mat
  if(is.null(n_factors)){
    n_factors <- ncol(Lambda_mat)
  }
  f_names <- colnames(Lambda_mat)
  if(is.null(f_names)){
    f_names <- paste0("F", 1:n_factors)
  }
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  par_names <- par_names[!is_nuisance]
  n_pars <- length(par_names)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    lambda = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, f_names, NULL)),
    lambda_untransf = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, f_names, NULL)),
    epsilon_inv = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL)),
    psi_inv = array(NA_real_, dim = c(n_factors, iters), dimnames = list(NULL, NULL)),
    eta = array(NA_real_, dim = c(n_subjects, n_factors, iters), dimnames = list(subject_ids, NULL, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}


add_info_factor <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  args <- list(...)
  n_factors <- args$n_factors
  Lambda_mat <- args$Lambda_mat
  if(is.null(n_factors)) n_factors <- ncol(Lambda_mat)
  n_pars <- sum(!sampler$nuisance)
  if(is.null(Lambda_mat)){
    Lambda_mat <- matrix(Inf, nrow = n_pars, ncol = n_factors)
    diag(Lambda_mat) <- 1
    Lambda_mat[upper.tri(Lambda_mat, diag = F)] <- 0
  }
  if(!is.null(args$signFix)){
    signFix <- args$signFix
  } else{
    signFix <- F
  }

  attr(sampler, "signFix") <- signFix
  attr(sampler, "Lambda_mat") <- Lambda_mat

  sampler$prior <- get_prior_factor(prior, sum(!sampler$nuisance), sample = F, n_factors = n_factors)
  sampler$n_factors <- n_factors
  return(sampler)
}

get_prior_factor <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL,
                             Lambda_mat = NULL, n_factors = NULL){


  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_pars(design, doMap = F))
  }
  if(is.null(Lambda_mat)) Lambda_mat <- attr(prior, "Lambda_mat")
  if(is.null(Lambda_mat)){
    Lambda_mat <- matrix(Inf, nrow = n_pars, ncol = n_factors)
    diag(Lambda_mat) <- 1
    Lambda_mat[upper.tri(Lambda_mat, diag = F)] <- 0
  }
  if(is.null(n_factors)) n_factors <- ncol(Lambda_mat)

  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- rep(1, n_pars)
  }
  if(is.null(prior$theta_lambda_var)){
    prior$theta_lambda_var <- rep(.2, n_pars)
  }
  if(is.null(prior$ap)){
    prior$ap <- 2
  }
  if(is.null(prior$bp)){
    prior$bp <- .5
  }
  if(is.null(prior$as)){
    prior$as <- rep(2.5, n_pars)
  }
  if(is.null(prior$bs)){
    prior$bs <- rep(.1, n_pars)
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- diag(1/prior$theta_mu_var)
  prior$theta_lambda_invar <-1/prior$theta_lambda_var
  # Things I save rather than re-compute inside the loops.
  attr(prior, "type") <- "factor"
  attr(prior, "Lambda_mat") <- Lambda_mat
  out <- prior
  if(sample){
    samples <- list()
    par_names <- names(sampled_pars(design, doMap = F))
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
    if(selection %in% c("loadings", "std_loadings", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      lambda <- array(0, dim = c(n_pars, n_factors, N))
      for(i in 1:n_factors){
        lambda[,i,] <- t(mvtnorm::rmvnorm(N, sigma = diag(prior$theta_lambda_var)))
      }
      lambda <- constrain_lambda(lambda, Lambda_mat)
      rownames(lambda) <- par_names
      if(is.null(colnames(Lambda_mat))){
        colnames(lambda) <- paste0("F", 1:n_factors)
      } else{
        colnames(lambda) <- colnames(Lambda_mat)
      }
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
        psi <- 1/rgamma(n_factors, prior$ap, prior$bp)
        loadings <- lambda[,,i]
        vars[,,i] <- loadings %*% diag(psi, n_factors) %*% t(loadings) + diag(sigma)
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

get_startpoints_factor<- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!pmwgs$nuisance)
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$prior$theta_mu_mean, sd = sqrt(pmwgs$prior$theta_mu_var))
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(n_pars * 3,diag(n_pars))
  start_psi_inv <- rep(1, pmwgs$n_factors)
  start_sig_err_inv <- rep(1, n_pars)
  start_lambda <- matrix(0, nrow = n_pars, ncol = pmwgs$n_factors)
  Lambda_mat <- attr(pmwgs, "Lambda_mat")
  start_lambda[Lambda_mat != Inf] <- Lambda_mat[Lambda_mat != Inf]
  start_eta <- matrix(0, nrow = pmwgs$n_subjects, ncol = pmwgs$n_factors)
  return(list(tmu = start_mu, tvar = start_var, lambda = start_lambda, lambda_untransf = start_lambda,
              sig_err_inv = start_sig_err_inv, psi_inv = start_psi_inv,
              eta = start_eta))
}

fill_samples_factor <- function(samples, group_level, proposals, j = 1, n_pars){
  samples$lambda[,,j] <- group_level$lambda
  samples$lambda_untransf[,,j] <- group_level$lambda_untransf
  samples$epsilon_inv[,j] <- group_level$sig_err_inv
  samples$psi_inv[,j] <- group_level$psi_inv
  samples$eta[,,j] <- group_level$eta
  samples <- fill_samples_base(samples, group_level, proposals, j = j, n_pars)
  return(samples)
}

gibbs_step_factor <- function(sampler, alpha){
  # Gibbs step for group means with parameter expanded factor analysis from Ghosh & Dunson 2009
  # mu = theta_mu, var = theta_var
  last <- last_sample_factor(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  # extract previous values (for ease of reading)
  alpha <- t(alpha)
  n_subjects <- sampler$n_subjects
  n_pars <- sum(!sampler$nuisance)
  n_factors <- sampler$n_factors

  # Save the original constraint matrix (with fixed values) for later rescaling.
  Lambda_constraints <- hyper$Lambda_mat

  # Create a binary indicator: TRUE for free parameters (Inf) and FALSE for fixed ones.
  Lambda_mat <- (Lambda_constraints == Inf)

  eta <- matrix(last$eta, n_subjects, n_factors)
  psi_inv <- diag(last$psi_inv, n_factors)
  sig_err_inv <- diag(last$sig_err_inv)
  lambda <- matrix(last$lambda, n_pars, n_factors)
  mu <- last$mu

  # Update mu
  mu_sig <- solve(n_subjects * sig_err_inv + prior$theta_mu_invar)
  mu_mu <- mu_sig %*% (sig_err_inv %*% colSums(alpha - eta %*% t(lambda)) +
                         prior$theta_mu_invar %*% prior$theta_mu_mean)
  mu <- rmvnorm(1, mu_mu, mu_sig)
  colnames(mu) <- colnames(alpha)

  # calculate mean-centered observations
  alphatilde <- sweep(alpha, 2, mu)

  # Update eta (latent factors)
  eta_sig <- solve(psi_inv + t(lambda) %*% sig_err_inv %*% lambda)
  eta_mu <- eta_sig %*% t(lambda) %*% sig_err_inv %*% t(alphatilde)
  eta[,] <- t(apply(eta_mu, 2, FUN = function(x){ rmvnorm(1, x, eta_sig) }))

  # Update sig_err (error precisions)
  sig_err_inv <- diag(rgamma(n_pars, shape = prior$as + n_subjects/2,
                             rate = prior$bs + colSums((alphatilde - eta %*% t(lambda))^2)/2))

  # Update lambda (factor loadings) for free entries only
  for (j in 1:n_pars) {
    constraint <- Lambda_mat[j,]  # TRUE for free parameters
    if(any(constraint)){  # Only update if there are free entries in row j
      etaS <- eta[, constraint]
      lambda_sig <- solve(sig_err_inv[j,j] * t(etaS) %*% etaS +
                            diag(prior$theta_lambda_invar[j], sum(constraint)))
      lambda_mu <- (lambda_sig * sig_err_inv[j,j]) %*% (t(etaS) %*% alphatilde[,j])
      lambda[j, constraint] <- rmvnorm(1, lambda_mu, lambda_sig)
    }
  }

  # Update psi_inv (latent factor precisions)
  psi_inv[,] <- diag(rgamma(n_factors, shape = prior$ap + n_subjects/2,
                            rate = prior$bp + colSums(eta^2)/2), n_factors)
  # Optionally update via inverse Wishart if desired:
  # psi_inv <- diag(n_factors) # or use riwish update

  # **** New Rescaling Step to Enforce the Marker Constraints ****
  # For each factor, locate the marker (fixed value) and rescale the factor loadings and scores.
  # for (j in 1:n_factors) {
  #   marker_idx <- which(!is.infinite(Lambda_constraints[, j]) & Lambda_constraints[,j] != 0)
  #   if(length(marker_idx) == 0){
  #     stop(sprintf("No constraint found for factor %d", j))
  #   }
  #   # Choose the first marker in the column as the anchor.
  #   marker_row <- marker_idx[1]
  #   fixed_val <- Lambda_constraints[marker_row, j]  # e.g., should be 1 or another constant
  #   # Compute scale factor: how far is the current loading from the fixed value?
  #   scale_factor <- fixed_val / lambda[marker_row, j]
  #   # Rescale the entire j-th column of lambda and adjust eta accordingly.
  #   lambda[, j] <- lambda[, j] * scale_factor
  #   eta[, j] <- eta[, j] / scale_factor
  # }
  # **** End Rescaling Step ****

  # The rest of the code (e.g., signFix) can follow if desired.
  # for(l in 1:n_factors){
  #   mult <- ifelse(lambda[l, l] < 0, -1, 1)
  #   lambda[,l] <- mult * lambda[, l]
  # }

  var <- lambda %*% solve(psi_inv) %*% t(lambda) + diag(1/diag(sig_err_inv))
  lambda <- lambda %*% matrix(diag(sqrt(1/diag(psi_inv)), n_factors), nrow = n_factors)

  return(list(tmu = mu, tvar = var, lambda_untransf = lambda,
              lambda = lambda, eta = eta,
              sig_err_inv = diag(sig_err_inv), psi_inv = diag(psi_inv), alpha = t(alpha)))
}

last_sample_factor <- function(store) {
  list(
    mu = store$theta_mu[, store$idx],
    eta = store$eta[,,store$idx],
    lambda = store$lambda_untransf[,,store$idx],
    psi_inv = store$psi_inv[,store$idx],
    sig_err_inv = store$epsilon_inv[,store$idx]
  )
}

get_conditionals_factor <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  sig_err <- log(samples$epsilon_inv[idx,])
  psi <- log(samples$psi_inv)
  eta <- matrix(samples$eta[s,,], nrow = samples$n_factors)
  lambda <- apply(samples$lambda_untransf[idx,,,drop = F], 3, unwind_lambda, samples$Lambda_mat[idx,])
  theta_mu <- samples$theta_mu[idx,]
  all_samples <- rbind(samples$alpha[idx, s,],theta_mu, eta, sig_err, psi, lambda)#, sig_err, psi, lambda)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration],
                                 samples$eta[s,,iteration],
                                 log(samples$epsilon_inv[idx, iteration]),
                                 log(samples$psi_inv[,iteration, drop = F]),
                                 unwind_lambda(samples$lambda_untransf[idx,, iteration], samples$Lambda_mat[idx,])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

filtered_samples_factor <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter],
    lambda_untransf = sampler$samples$lambda_untransf[, , filter, drop = F],
    psi_inv = sampler$samples$psi_inv[, filter, drop = F],
    epsilon_inv = sampler$samples$epsilon_inv[, filter],
    eta = sampler$samples$eta[, , filter, drop = F],
    theta_var = sampler$samples$theta_var[,,filter],
    alpha = sampler$samples$alpha[, , filter],
    Lambda_mat = attributes(sampler)$Lambda_mat,
    n_factors = sampler$n_factors,
    iteration = length(filter)
  )
}
unwind_lambda <- function(lambda, constraintMat, reverse = F){
  if(reverse){
    out <- constraintMat
    out[constraintMat == Inf] <- lambda
  } else{
    out <- as.numeric(lambda[constraintMat == Inf])
  }
  return(out)
}

constrain_lambda <- function(lambda, constraintMat){
  for(i in 1:dim(lambda)[3]){
    tmp <- lambda[,,i]
    tmp[constraintMat != Inf] <- 0
    lambda[,,i] <- tmp
  }
  return(lambda)
}

# bridge_sampling ---------------------------------------------------------

bridge_add_info_factor <- function(info, samples){
  info$n_factors <- samples$n_factors
  info$Lambda_mat <- attr(samples, "Lambda_mat")
  # Free loadings + group-level mean + factor variances + residual variances
  # note not factor scores, since we use the marginal model as adviced by Merkle et al. 2023
  info$group_idx <- (samples$n_pars*samples$n_subjects + 1):(samples$n_pars*samples$n_subjects + sum(info$Lambda_mat == Inf) + samples$n_pars + samples$n_factors + samples$n_pars)
  return(info)
}

bridge_add_group_factor <- function(all_samples, samples, idx){
  Lambda_mat <- attr(samples, "Lambda_mat")
  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[,idx]))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$lambda_untransf[,,idx,drop = F], 3, unwind_lambda, Lambda_mat), ncol = nrow(all_samples))))
  all_samples <- cbind(all_samples, t(log(samples$samples$epsilon_inv[,idx])))
  all_samples <- cbind(all_samples, t(log(samples$samples$psi_inv[,idx, drop = F])))
  return(all_samples)
}

bridge_group_and_prior_and_jac_factor <- function(proposals_group, proposals_list, info){
  prior <- info$prior
  proposals <- do.call(cbind, proposals_list)
  theta_mu <- proposals_group[,1:info$n_pars]
  lambda <- proposals_group[,(info$n_pars +1):(info$n_pars + sum(info$Lambda_mat == Inf))]
  theta_epsilon_inv <- proposals_group[,(1 + info$n_pars + sum(info$Lambda_mat == Inf)): (info$n_pars + sum(info$Lambda_mat == Inf) + info$n_pars)]
  psi_inv <- proposals_group[,(1 + info$n_pars + sum(info$Lambda_mat == Inf) + info$n_pars):
                                     (info$n_pars + sum(info$Lambda_mat == Inf) + info$n_pars + info$n_factors), drop = F]

  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)
  for(i in 1:n_iter){ # these unfortunately can't be vectorized
    lambda_curr <- unwind_lambda(lambda[i,], info$Lambda_mat, reverse = T)
    epsilon_curr <- diag(1/exp(theta_epsilon_inv[i,]))
    psi_curr <- diag(1/exp(psi_inv[i,]), info$n_factors)
    theta_var_curr <- lambda_curr %*% psi_curr %*% t(lambda_curr) + epsilon_curr
    proposals_curr <- matrix(proposals[i,], ncol = info$n_pars, byrow = T)
    group_ll <- sum(dmvnorm(proposals_curr, theta_mu[i,], theta_var_curr, log = T))
    prior_epsilon <- sum(logdinvGamma(1/exp(theta_epsilon_inv[i,]), shape = prior$as, rate = prior$bs))
    prior_psi <- sum(logdinvGamma(1/exp(psi_inv[i,]), prior$ap, rate = prior$bp))
    sum_out[i] <- group_ll + prior_epsilon + prior_psi
  }

  prior_lambda <- dmvnorm(lambda, mean = rep(0, ncol(lambda)),
                          sigma = diag(prior$theta_lambda_var, ncol(lambda)), log = T)
  prior_mu <- dmvnorm(theta_mu, mean = prior$theta_mu_mean, sigma = diag(prior$theta_mu_var), log =T)

  jac_psi <- rowSums(theta_epsilon_inv)
  jac_epsilon <- rowSums(psi_inv)
  return(sum_out + prior_mu + prior_lambda + jac_psi + jac_epsilon) # Output is of length nrow(proposals)
}

