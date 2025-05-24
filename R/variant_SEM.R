sample_store_SEM <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, ...) {
  args <- list(...)
  n_factors <- ncol(args$Lambda_mat)
  covariates <- as.matrix(args$covariates)
  n_cov <- ncol(covariates)
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  par_names <- par_names[!is_nuisance]
  n_pars <- length(par_names)
  x_names <- colnames(covariates)
  factor_names <- paste0("F", 1:n_factors)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    lambda = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, factor_names, NULL)),
    B = array(NA_real_,dim = c(n_factors, n_factors, iters),dimnames = list(factor_names, factor_names, NULL)),
    epsilon_inv = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    delta_inv = array(NA_real_, dim = c(n_factors, n_factors, iters), dimnames = list(factor_names, factor_names, NULL)),
    eta = array(NA_real_, dim = c(n_subjects, n_factors, iters), dimnames = list(subject_ids, factor_names, NULL)),
    K = array(NA_real_, dim = c(n_pars, n_cov, iters), dimnames = list(par_names, x_names, NULL)),
    G = array(NA_real_, dim = c(n_factors, n_cov, iters), dimnames = list(factor_names, x_names, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}


add_info_SEM <- function(sampler, prior = NULL, ...){
  args <- list(...)
  sem_settings <- if (!is.null(args$sem_settings)) args$sem_settings else list()

  # Extract matrices and covariates from sem_settings
  Lambda_mat <- sem_settings$Lambda_mat
  n_pars <- sum(!sampler$nuisance)
  n_factors <- if (!is.null(Lambda_mat)) ncol(Lambda_mat) else {
    if (!is.null(sem_settings$B_mat)) ncol(sem_settings$B_mat)
    else if (!is.null(sem_settings$G_mat)) nrow(sem_settings$G_mat)
    else if (!is.null(sem_settings$factor_names)) length(sem_settings$factor_names) # Fallback if only factor_names provided
    else stop("n_factors cannot be determined if Lambda_mat, B_mat, G_mat, and factor_names are all NULL in sem_settings.")
  }
  if(is.null(Lambda_mat) && n_factors > 0){
    Lambda_mat <- matrix(0, nrow = n_pars, ncol = n_factors)
  }
  sem_settings$Lambda_mat <- Lambda_mat

  B_mat <- sem_settings$B_mat
  if(is.null(B_mat) && n_factors > 0){
    B_mat <- matrix(0, nrow = n_factors, ncol = n_factors)
  }
  sem_settings$B_mat <- B_mat

  # Get covariates from sem_settings
  covariates_from_sem <- sem_settings$covariates
  if (is.null(covariates_from_sem)) {
    # If covariates are not in sem_settings (e.g. old usage or direct provision),
    # try to get them from args for backward compatibility or other scenarios.
    # However, the primary expectation is they come from define_sem_structure via sem_settings.
    if (!is.null(args$covariates)) {
        covariates_matrix <- as.matrix(args$covariates)
        warning("Covariates were taken from '...' args instead of 'sem_settings$covariates'. Ensure 'sem_settings' is correctly populated by 'define_sem_structure'.")
    } else {
        covariates_matrix <- matrix(0, nrow = sampler$n_subjects, ncol = 0) # Default to 0 covariates if none found
    }
  } else {
    covariates_matrix <- as.matrix(covariates_from_sem)
  }
  n_cov <- ncol(covariates_matrix)
  # Store the definitive covariates and n_cov in sem_settings for consistency
  sem_settings$covariates <- covariates_matrix
  sem_settings$n_cov <- n_cov

  K_mat <- sem_settings$K_mat
  if(is.null(K_mat)){
    K_mat <- matrix(0, nrow = n_pars, ncol = n_cov)
  }
  sem_settings$K_mat <- K_mat

  G_mat <- sem_settings$G_mat
  if(is.null(G_mat) && n_factors > 0){
    G_mat <- matrix(0, nrow = n_factors, ncol = n_cov)
  }
  sem_settings$G_mat <- G_mat

  factor_groups <- sem_settings$factor_groups
  if (is.null(factor_groups)) {
    if (n_factors > 0) factor_groups <- 1:n_factors else factor_groups <- integer(0)
  } else {
    if (length(factor_groups) != n_factors) {
      stop("Length of factor_groups in sem_settings must be equal to n_factors.")
    }
  }
  sem_settings$factor_groups <- factor_groups
  sem_settings$n_factors <- n_factors # Ensure n_factors is also in sem_settings

  sampler$sem_settings <- sem_settings
  # Remove legacy direct storage if they existed
  sampler$covariates <- NULL
  sampler$n_cov <- NULL

  sampler$prior <- get_prior_SEM(prior, n_pars = n_pars, sample = F,
                                 sem_settings = sampler$sem_settings, # Pass the whole sem_settings
                                 covariates = sampler$sem_settings$covariates) # Explicitly pass covariates from settings
  sampler$n_factors <- n_factors # Retain for direct access if needed elsewhere, though primary is via sem_settings
  return(sampler)
}

get_prior_SEM <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL,
                          sem_settings = NULL,
                          covariates = NULL){

  Lambda_mat <- if (!is.null(sem_settings)) sem_settings$Lambda_mat else NULL
  B_mat <- if (!is.null(sem_settings)) sem_settings$B_mat else NULL
  K_mat <- if (!is.null(sem_settings)) sem_settings$K_mat else NULL
  factor_groups_prior <- if (!is.null(sem_settings)) sem_settings$factor_groups else NULL

  n_factors <- if(!is.null(Lambda_mat)) ncol(Lambda_mat) else if(!is.null(B_mat)) ncol(B_mat) else 0
  n_cov <- if(!is.null(covariates)) ncol(covariates) else if(!is.null(K_mat)) ncol(K_mat) else 0


  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_pars(design, doMap = F))
  }
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- rep(1, n_pars)
  }
  if(is.null(prior$lambda_var)){
    prior$lambda_var <- rep(.7, n_pars)
  }
  if(is.null(prior$K_var)){
    prior$K_var <- rep(1, n_pars)
  }
  if(is.null(prior$B_var)){
    prior$B_var <- rep(.5, n_factors)
  }
  if(is.null(prior$G_var)){
    prior$G_var <- rep(.5, n_factors)
  }
  if(is.null(prior$a_d)){
    prior$a_d <- 5
  }
  if(is.null(prior$b_d)){
    prior$b_d <- .3
  }
  if(is.null(prior$a_e)){
    prior$a_e <- rep(5, n_pars)
  }
  if(is.null(prior$b_e)){
    prior$b_e <- rep(.3, n_pars)
  }
  attr(prior, "type") <- "SEM"
  out <- prior
  if(sample){
    x_mu <- colMeans(covariates)
    x_var <- cov(covariates)
    isFree_B <- B_mat == Inf
    if (is.null(factor_groups_prior)) {
        if (!is.null(design) && !is.null(attr(design, "factor_groups"))) {
          factor_groups_prior <- attr(design, "factor_groups")
        } else if (n_factors > 0) {
          factor_groups_prior <- 1:n_factors
        } else {
          factor_groups_prior <- integer(0)
        }
    }
    unique_factor_groups_prior <- unique(factor_groups_prior)

    factor_names <- if (n_factors > 0) paste0("F", 1:n_factors) else character(0)
    samples <- list()
    par_names <- names(sampled_pars(design, doMap = F))
    if(selection %in% c("mu", "alpha", "mu_implied")){
      mu <- t(mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
                               sigma = diag(prior$theta_mu_var)))
      rownames(mu) <- par_names
      if(selection %in% c("mu")){
        samples$theta_mu <- mu
      }
    }
    if(selection %in% c("regressors", "alpha", "mu_implied", "Sigma", "correlation", "covariance", "sigma2")){
      K <- array(0, dim = c(n_pars, n_cov, N))
      for(i in 1:n_cov){
        K[,i,] <- t(mvtnorm::rmvnorm(N, sigma = diag(prior$K_var)))
      }
      K <- constrain_lambda(K, K_mat)
      rownames(K) <- par_names
      colnames(K) <- colnames(covariates)
      if(selection %in% 'regressors'){
        samples$K <- K
      }
    }
    if(selection %in% c("factor_regressors", "alpha", "mu_implied", "Sigma", "correlation", "covariance", "sigma2")){
      G <- array(0, dim = c(n_factors, n_cov, N))
      for(i in 1:n_cov){
        G[,i,] <- t(mvtnorm::rmvnorm(N, sigma = diag(prior$G_var)))
      }
      G <- constrain_lambda(G, G_mat)

      rownames(G) <- factor_names
      colnames(G) <- colnames(covariates)
      if(selection %in% 'factor_regressors'){
        samples$G <- G
      }
    }
    if(selection %in% c("structural_regressors", "alpha", "mu_implied", "Sigma", "correlation", "covariance", "sigma2")){
      B <- array(0, dim = c(n_factors, n_factors, N))
      for(i in 1:n_factors){
        B[,i,] <- t(mvtnorm::rmvnorm(N, sigma = diag(prior$B_var)))
      }
      B <- constrain_lambda(B, B_mat)

      rownames(B) <- colnames(B) <- factor_names
      if(selection %in% 'structural_regressors'){
        samples$B <- B
      }
    }
    if(selection %in% c("loadings", "alpha", "mu_implied", "Sigma", "correlation", "covariance", "sigma2")){
      lambda <- array(0, dim = c(n_pars, n_factors, N))
      for(i in 1:n_factors){
        lambda[,i,] <- t(mvtnorm::rmvnorm(N, sigma = diag(prior$lambda_var)))
      }
      lambda <- constrain_lambda(lambda, Lambda_mat)
      rownames(lambda) <- par_names
      colnames(lambda) <- factor_names
      if(selection %in% "loadings"){
        samples$lambda <- lambda
      }
    }
    if(selection %in% c("residuals", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      epsilon_inv <- t(matrix(rgamma(n_pars*N, shape = prior$a_e, rate = prior$b_e),
                            ncol = n_pars, byrow = T))
      rownames(epsilon_inv) <- par_names
      if(selection %in% "residuals"){
        samples$epsilon_inv <- epsilon_inv
      }
    }
    if(selection %in% c("factor_residuals", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      delta_inv <- array(0, dim = c(n_factors, n_factors, N))
      for(i in 1:N){
        for (group_id in unique_factor_groups_prior) {
          group_idx <- which(factor_groups_prior == group_id)
          d_block <- length(group_idx)
          if (d_block > 1) {
            S_iw_prior <- diag(prior$b_d, d_block)
            df_iw_prior <- max(prior$a_d, d_block + 2)
            sampled_cov_block_prior <- riwish(df_iw_prior, S_iw_prior)
            delta_inv[group_idx, group_idx, i] <- solve(sampled_cov_block_prior)
          } else if (d_block == 1) {
            delta_inv[group_idx, group_idx, i] <- rgamma(1, shape = prior$a_d, rate = prior$b_d)
          }
        }
      }
      rownames(delta_inv) <- colnames(delta_inv) <- factor_names
      if(selection %in% "factor_residuals"){
        samples$delta_inv <- delta_inv
      }
    }
    if(selection %in% c("mu_implied", "alpha")) {
      mu_implied <- mu
      for(i in 1:N){
        B_0_inv <- solve(diag(n_factors) - B[,,i])
        mu_implied[,i] <- c(mu[,i] + lambda[,,i] %*% B_0_inv %*% G[,,i] %*% x_mu + K[,,i] %*% x_mu)
      }
      if(selection != "alpha") samples$mu_implied <- mu_implied
    }
    if(selection %in% c("sigma2", "covariance", "correlation", "Sigma", "alpha")) {
      vars <- array(NA_real_, dim = c(n_pars, n_pars, N))
      colnames(vars) <- rownames(vars) <- par_names
      for(i in 1:N){
        B_0_inv <- solve(diag(n_factors) - as.matrix(B[,,i]))
        vars[,,i] <- as.matrix(lambda[,,i]) %*% B_0_inv %*% (as.matrix(G[,,i]) %*% x_var %*% t(as.matrix(G[,,i])) +
                     solve(as.matrix(delta_inv[,,i]))) %*% t(B_0_inv) %*% t(as.matrix(lambda[,,i])) +
                     as.matrix(K[,,i]) %*% x_var %*% t(as.matrix(K[,,i])) + diag(1/epsilon_inv[,i])
      }
      if(selection != "alpha") samples$theta_var <- vars
    }
    if(selection %in% "alpha"){
      samples$alpha <- get_alphas(mu_implied, vars, "alpha")
    }
    out <- samples
  }
  return(out)
}

get_startpoints_SEM<- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!pmwgs$nuisance)
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$prior$theta_mu_mean, sd = sqrt(pmwgs$prior$theta_mu_var))
  if (is.null(start_var)) start_var <- riwish(n_pars * 3,diag(n_pars))
  start_delta_inv <- diag(1, pmwgs$n_factors)
  start_epsilon_inv <- diag(1, n_pars)
  start_eta <- matrix(0, nrow = pmwgs$n_subjects, ncol = pmwgs$n_factors)

  start_lambda <- matrix(0, nrow = n_pars, ncol = pmwgs$n_factors)
  start_B <- matrix(0, nrow = pmwgs$n_factors, ncol = pmwgs$n_factors)
  start_K <- matrix(0, nrow = n_pars, ncol = pmwgs$n_cov)
  start_G <- matrix(0, nrow = pmwgs$n_factors, ncol = pmwgs$n_cov)

  sem_settings <- pmwgs$sem_settings
  Lambda_mat <- sem_settings$Lambda_mat
  B_mat <- sem_settings$B_mat
  K_mat <- sem_settings$K_mat
  G_mat <- sem_settings$G_mat

  start_lambda[Lambda_mat != Inf] <- Lambda_mat[Lambda_mat != Inf]
  start_B[B_mat != Inf] <- B_mat[B_mat != Inf]
  start_K[K_mat != Inf] <- K_mat[K_mat != Inf]
  start_G[G_mat != Inf] <- G_mat[G_mat != Inf]
  return(list(tmu = start_mu, tvar = start_var, lambda = start_lambda, B = start_B,
              K = start_K, G = start_G,
              epsilon_inv = start_epsilon_inv, delta_inv = start_delta_inv,
              eta = start_eta, sub_mu = start_mu))
}

fill_samples_SEM <- function(samples, group_level, proposals, j = 1, n_pars){
  samples$lambda[,,j] <- group_level$lambda
  samples$B[,,j] <- group_level$B
  samples$K[,,j] <- group_level$K
  samples$G[,,j] <- group_level$G
  samples$epsilon_inv[,,j] <- group_level$epsilon_inv
  samples$delta_inv[,,j] <- group_level$delta_inv
  samples$eta[,,j] <- group_level$eta
  samples <- fill_samples_base(samples, group_level, proposals, j = j, n_pars)
  return(samples)
}

gibbs_step_SEM <- function(sampler, alpha){
  last          <- last_sample_SEM(sampler$samples)
  sem_settings  <- sampler$sem_settings
  prior         <- sampler$prior

  y             <- t(alpha)                     # subjects × variables
  n_subjects    <- sampler$n_subjects
  n_pars        <- sum(!sampler$nuisance)
  n_factors     <- sampler$n_factors
  n_cov         <- sampler$n_cov
  covariates    <- sampler$covariates

  ## free‑parameter masks ----------------------------------------------------
  isFree_Lambda <- sem_settings$Lambda_mat == Inf
  isFree_B      <- sem_settings$B_mat     == Inf
  isFree_K      <- sem_settings$K_mat     == Inf
  isFree_G      <- sem_settings$G_mat     == Inf

  factor_groups <- sem_settings$factor_groups
  unique_factor_groups <- unique(factor_groups)

  ## current state -----------------------------------------------------------
  eta         <- matrix(last$eta,     n_subjects, n_factors)
  delta_inv   <- matrix(last$delta_inv, n_factors, n_factors)
  epsilon_inv <- last$epsilon_inv
  lambda      <- matrix(last$lambda,  n_pars,    n_factors)
  B           <- matrix(last$B,       n_factors, n_factors)
  K           <- matrix(last$K,       n_pars,    n_cov)
  G           <- matrix(last$G,       n_factors, n_cov)
  mu          <- last$mu

  ## ---- update mu ----------------------------------------------------------
  mu_sig <- solve(n_subjects * epsilon_inv +
                  diag(1 / prior$theta_mu_var, n_pars))
  mu_mu  <- mu_sig %*%
            (epsilon_inv %*%
               colSums(y - covariates %*% t(K) - eta %*% t(lambda)) +
             diag(1/prior$theta_mu_var, n_pars) %*% prior$theta_mu_mean)
  mu <- drop(rmvnorm(1, mu_mu, mu_sig))
  ytilde <- sweep(y, 2, mu)

  ## ---- update eta ---------------------------------------------------------
  B0_inv   <- solve(diag(n_factors) - B)
  Psi0_inv <- solve(B0_inv %*% solve(delta_inv) %*% t(B0_inv))
  eta_sig  <- solve(Psi0_inv + t(lambda) %*% epsilon_inv %*% lambda)
  for (i in seq_len(n_subjects)){
    eta_mean <- eta_sig %*%
      ( t(lambda) %*% epsilon_inv %*%
          (ytilde[i,] - K %*% covariates[i,]) +
        Psi0_inv %*% B0_inv %*% (G %*% covariates[i,]) )
    eta[i,] <- rmvnorm(1, eta_mean, eta_sig)
  }

  ## ---- update epsilon (item precision) ------------------------------------
  epsilon_inv <- diag(rgamma(
    n_pars,
    shape = prior$a_e + n_subjects/2,
    rate  = prior$b_e + 0.5 * colSums((ytilde -
                                        covariates %*% t(K) -
                                        eta %*% t(lambda))^2)))

  ## ---- update loadings: lambda and K --------------------------------------
  lambda_y       <- cbind(K, lambda)
  lambda_y_prior <- cbind(matrix(prior$K_var,     n_pars, n_cov),
                          matrix(prior$lambda_var,n_pars, n_factors))
  for (j in seq_len(n_pars)){
    isFree <- c(isFree_K[j,], isFree_Lambda[j,])
    if(any(isFree)){
      etaS    <- cbind(covariates, eta)[, isFree, drop = FALSE]
      lam_sig <- solve(epsilon_inv[j,j] * crossprod(etaS) +
                       diag(1/lambda_y_prior[j,isFree], sum(isFree)))
      lam_mu  <- lam_sig %*% (epsilon_inv[j,j] * crossprod(etaS, ytilde[,j]))
      lambda_y[j, isFree] <- rmvnorm(1, lam_mu, lam_sig)
    }
  }
  K      <- lambda_y[, seq_len(n_cov),        drop = FALSE]
  lambda <- lambda_y[, -(seq_len(n_cov)),     drop = FALSE]

  ## ---- correlated update for G and B --------------------------------------
  G_new <- G
  B_new <- B

  # Design matrix shared across rows
  Z_full <- cbind(covariates, eta)   # n_subjects × (n_cov + n_factors)

  # Start with residuals under current coefficients
  eta_residuals <- eta - eta %*% t(B_new) - covariates %*% t(G_new)

  for (p in seq_len(n_factors)){
    free_idx <- c(isFree_G[p,], isFree_B[p,])
    if(!any(free_idx)) next

    Z_p <- Z_full[, free_idx, drop = FALSE]
    delta_pp <- delta_inv[p,p]

    # cross‑residual term capturing correlation with other factors
    c_vec <- if (n_factors > 1){
      eta_residuals[, -p, drop = FALSE] %*% delta_inv[-p, p, drop = FALSE]
    } else {
      matrix(0, n_subjects, 1)
    }
    y_star <- delta_pp * eta[, p, drop = FALSE] + c_vec

    B_prior_vec <- c(rep(prior$G_var, n_cov),
                     rep(prior$B_var, n_factors))[free_idx]

    B_sig <- solve(delta_pp * crossprod(Z_p) +
                   diag(1 / B_prior_vec, length(B_prior_vec)))
    B_mu  <- B_sig %*% crossprod(Z_p, y_star)

    coef_sample <- rmvnorm(1, B_mu, B_sig)

    # Write back sampled coefficients
    if(any(isFree_G[p,]))
      G_new[p, isFree_G[p,]] <- coef_sample[seq_len(n_cov)][isFree_G[p,]]
    if(any(isFree_B[p,]))
      B_new[p, isFree_B[p,]] <- coef_sample[-seq_len(n_cov)][isFree_B[p,]]

    # refresh residuals with new row before next iteration
    eta_residuals[,p] <- eta[,p] - Z_full[, free_idx, drop = FALSE] %*% t(coef_sample)
  }
  G <- G_new
  B <- B_new

  ## ---- update delta_inv ---------------------------------------------------
  eta_residuals <- eta - eta %*% t(B) - covariates %*% t(G)

  for (group_id in unique_factor_groups){
    group_idx <- which(factor_groups == group_id)
    d_block   <- length(group_idx)
    resid_blk <- eta_residuals[, group_idx, drop = FALSE]
    cov_block <- crossprod(resid_blk)

    S_iw  <- diag(prior$b_d, d_block)
    df_iw <- prior$a_d + n_subjects
    if(df_iw <= d_block - 1)
      stop("Inverse‑Wishart degrees‑of‑freedom too small for factor group ", group_id)

    if(d_block > 1){
      sampled_cov <- riwish(df_iw, S_iw + cov_block)
      delta_inv[group_idx, group_idx] <- solve(sampled_cov)
    } else {
      delta_inv[group_idx, group_idx] <- rgamma(
        1,
        shape = prior$a_d + n_subjects/2,
        rate  = prior$b_d + 0.5 * cov_block)
    }
  }

  ## ---- derived population moments ----------------------------------------
  x_mu  <- colMeans(covariates)
  x_var <- cov(covariates)

  B_0_inv  <- solve(diag(n_factors) - B)
  pop_mean <- drop(mu + lambda %*% B_0_inv %*% G %*% x_mu + K %*% x_mu)
  pop_var  <- lambda %*% B_0_inv %*% (G %*% x_var %*% t(G) +
                                      solve(delta_inv)) %*%
              t(B_0_inv) %*% t(lambda) +
              K %*% x_var %*% t(K) +
              solve(epsilon_inv)

  ## ---- return -------------------------------------------------------------
  list(mu          = mu,
       lambda      = lambda,
       eta         = eta,
       B           = B,
       K           = K,
       G           = G,
       epsilon_inv = epsilon_inv,
       delta_inv   = delta_inv,
       alpha       = t(y),
       tmu         = pop_mean,
       tvar        = pop_var)
}



last_sample_SEM <- function(store) {
  list(
    mu = store$theta_mu[, store$idx],
    eta = store$eta[,,store$idx],
    lambda = store$lambda[,,store$idx],
    B = store$B[,,store$idx],
    K = store$K[,,store$idx],
    G = store$G[,,store$idx],
    delta_inv = store$delta_inv[,,store$idx],
    epsilon_inv = store$epsilon_inv[,,store$idx]
  )
}

get_group_level_SEM <- function(parameters, s){
  mu <- parameters$sub_mu
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}


get_conditionals_SEM <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  epsilon_inv <- log(apply(samples$epsilon_inv[idx,idx,],3 , diag))
  eta <- matrix(samples$eta[s,,], nrow = samples$n_factors)
  current_Lambda_mat <- if(!is.null(samples$sem_settings)) samples$sem_settings$Lambda_mat else samples$Lambda_mat
  lambda <- apply(samples$lambda[idx,,,drop = F], 3, unwind_lambda, current_Lambda_mat[idx,])
  theta_mu <- samples$theta_mu[idx,]
  all_samples <- rbind(samples$alpha[idx, s,],theta_mu, eta, epsilon_inv, lambda)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration],
                                 samples$eta[s,,iteration],
                                 log(diag(samples$epsilon_inv[idx,idx, iteration])),
                                 unwind_lambda(samples$lambda[idx,, iteration], current_Lambda_mat[idx,])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

filtered_samples_SEM <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter],
    lambda = sampler$samples$lambda[, , filter, drop = F],
    epsilon_inv = sampler$samples$epsilon_inv[,, filter],
    eta = sampler$samples$eta[, , filter, drop = F],
    alpha = sampler$samples$alpha[, , filter],
    n_factors = sampler$n_factors,
    iteration = length(filter),
    sem_settings = sampler$sem_settings
  )
}

group__IC_SEM <- function(emc, stage="sample",filter=NULL, ...){
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    return_mcmc = FALSE, merge_chains = TRUE)
  theta_mu <- get_pars(emc, selection = "mu_implied", stage = stage, filter = filter,
                       return_mcmc = FALSE, merge_chains = TRUE)
  theta_var <- get_pars(emc, selection = "Sigma", stage = stage, filter = filter,
                        return_mcmc = FALSE, merge_chains = TRUE, remove_constants = F)
  mean_alpha <- apply(alpha, 1:2, mean)
  mean_mu <- rowMeans(theta_mu)
  mean_var <- apply(theta_var, 1:2, mean)

  N <- ncol(theta_mu)
  lls <- numeric(N)
  if(list(...)$for_WAIC){
    lls <- matrix(NA, nrow = ncol(mean_alpha), ncol = N)
    for(i in 1:N){
      lls[,i] <- dmvnorm(t(alpha[,,i]), theta_mu[,i], theta_var[,,i], logd = T)
    }
    return(lls)
  }
  for(i in 1:N){
    lls[i] <- sum(dmvnorm(t(alpha[,,i]), theta_mu[,i], theta_var[,,i], logd = T))
  }
  minD <- -2*max(lls)
  mean_ll <- mean(lls)
  mean_pars_ll <-  sum(dmvnorm(t(mean_alpha), mean_mu, mean_var, logd = TRUE))
  Dmean <- -2*mean_pars_ll
  return(list(mean_ll = mean_ll, Dmean = Dmean,
              minD = minD))
}

get_mu_implied <- function(x, idx){
  mu <- x$samples$theta_mu[,idx]
  B <- x$samples$B[,,idx, drop = F]
  G <- x$samples$G[,,idx, drop = F]
  K <- x$samples$K[,,idx, drop = F]
  loadings <- x$samples$lambda[,,idx, drop = F]
  n_factors <- ncol(loadings)
  x_mu <- colMeans(x$covariates)
  for(i in 1:ncol(mu)){
    B_0_inv <- solve(diag(n_factors) - as.matrix(B[,,i]))
    mu[,i] <- as.matrix(mu[,i]) + as.matrix(loadings[,,i]) %*% B_0_inv %*% as.matrix(G[,,i]) %*% x_mu
    + as.matrix(K[,,i]) %*% x_mu
  }
  return(mu)
}


bridge_add_info_SEM <- function(info, samples){
  sem_settings <- samples$sem_settings
  if (is.null(sem_settings)) stop("sem_settings not found in samples object for bridge sampling.")

  info$Lambda_mat <- sem_settings$Lambda_mat
  info$B_mat <- sem_settings$B_mat
  info$K_mat <- sem_settings$K_mat
  info$G_mat <- sem_settings$G_mat
  info$factor_groups_bs <- sem_settings$factor_groups

  info$n_factors <- samples$n_factors
  info$n_cov <- samples$n_cov
  info$covariates <- samples$covariates
  free_regrs <- sum(info$Lambda_mat == Inf) + sum(info$B_mat == Inf) + sum(info$K_mat == Inf) + sum(info$G_mat == Inf)
  other <- samples$n_pars + samples$n_pars

  if (is.null(info$factor_groups_bs)) {
      info$factor_groups_bs <- 1:info$n_factors
  }
  unique_f_groups_bs <- unique(info$factor_groups_bs)
  n_factor_params_delta <- 0
  for (fg_id in unique_f_groups_bs) {
    fg_idx <- which(info$factor_groups_bs == fg_id)
    d_block_fg <- length(fg_idx)
    if (d_block_fg > 1) {
      n_factor_params_delta <- n_factor_params_delta + (d_block_fg * (d_block_fg + 1)) / 2
    } else {
      n_factor_params_delta <- n_factor_params_delta + 1
    }
  }
  other <- other + n_factor_params_delta

  info$group_idx <- (samples$n_pars*samples$n_subjects + 1):
                    (samples$n_pars*samples$n_subjects + free_regrs + other)
  return(info)
}


bridge_add_group_SEM <- function(all_samples, samples, idx){
  sem_settings <- samples$sem_settings
  if (is.null(sem_settings)) stop("sem_settings not found in samples object for bridge sampling.")

  Lambda_mat <- sem_settings$Lambda_mat
  B_mat <- sem_settings$B_mat
  K_mat <- sem_settings$K_mat
  G_mat <- sem_settings$G_mat
  factor_groups_bs_add <- sem_settings$factor_groups

  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[,idx]))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$lambda[,,idx,drop = F], 3, unwind_lambda, Lambda_mat), ncol = nrow(all_samples))))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$B[,,idx,drop = F], 3, unwind_lambda, B_mat), ncol = nrow(all_samples))))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$K[,,idx,drop = F], 3, unwind_lambda, K_mat), ncol = nrow(all_samples))))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$G[,,idx,drop = F], 3, unwind_lambda, G_mat), ncol = nrow(all_samples))))

  all_samples <- cbind(all_samples, t(log(matrix(apply(samples$samples$epsilon_inv[,,idx, drop = F], 3, diag), ncol = nrow(all_samples)))))

  if (is.null(factor_groups_bs_add)) {
      factor_groups_bs_add <- 1:samples$n_factors
  }
  unique_f_groups_bs_add <- unique(factor_groups_bs_add)

  for (fg_id in unique_f_groups_bs_add) {
    fg_idx <- which(factor_groups_bs_add == fg_id)
    d_block_fg <- length(fg_idx)
    current_delta_block_samples <- samples$samples$delta_inv[fg_idx, fg_idx, idx, drop = F]

    if (d_block_fg > 1) {
      all_samples <- cbind(all_samples, t(matrix(apply(current_delta_block_samples, 3, unwind_chol), ncol = nrow(all_samples))))
    } else {
      all_samples <- cbind(all_samples, t(log(matrix(apply(current_delta_block_samples, 3, diag), ncol = nrow(all_samples)))))
    }
  }

  return(all_samples)
}




bridge_group_and_prior_and_jac_SEM <- function(proposals_group,
                                               proposals_list,
                                               info)
{
  ## --------------------------------------------------------------------- ##
  ## 0.  CONSTANTS & SHORT-HANDS                                           ##
  ## --------------------------------------------------------------------- ##
  prior        <- info$prior
  n_pars       <- info$n_pars
  n_factors    <- info$n_factors
  n_cov        <- info$n_cov
  Lambda_mat   <- info$Lambda_mat
  B_mat        <- info$B_mat
  K_mat        <- info$K_mat
  G_mat        <- info$G_mat
  factor_groups <- info$factor_groups_bs   # guaranteed non-NULL by add_info helper
  unique_fg    <- unique(factor_groups)

  ## covariate moments used in population moments
  x_mu  <- colMeans(info$covariates)
  x_var <- cov(info$covariates)

  ## posterior samples of subject-level α stacked row-wise
  proposals <- do.call(cbind, proposals_list)

  ## --------------------------------------------------------------------- ##
  ## 1.  SLICE THE PROPOSAL VECTOR INTO NAMED BLOCKS                       ##
  ##     (each row = one draw)                                             ##
  ## --------------------------------------------------------------------- ##
  pos          <- 0
  grab <- function(len)
  {
    if (len == 0) return(NULL)
    out <- proposals_group[ , (pos + 1):(pos + len), drop = FALSE]
    pos <<- pos + len
    out
  }

  theta_mu <- grab(n_pars)

  v_lambda <- grab(sum(Lambda_mat == Inf))
  v_B      <- grab(sum(B_mat     == Inf))
  v_K      <- grab(sum(K_mat     == Inf))
  v_G      <- grab(sum(G_mat     == Inf))

  eps_log  <- grab(n_pars)           # log-precision for manifest residuals

  ## δ-blocks: list indexed by factor-group id
  delta_vec_list <- lapply(unique_fg, function(fg_id)
  {
    d  <- sum(factor_groups == fg_id)
    len <- if (d > 1) d * (d + 1) / 2 else 1
    grab(len)
  })
  names(delta_vec_list) <- as.character(unique_fg)

  stopifnot(pos == ncol(proposals_group))  # sanity

  n_iter <- nrow(proposals_group)
  log_post_plus_jac <- numeric(n_iter)

  ## --------------------------------------------------------------------- ##
  ## 2.  MAIN LOOP OVER ITERATIONS                                         ##
  ## --------------------------------------------------------------------- ##
  for (i in seq_len(n_iter)) {

    ## -------- 2a. Reconstruct free matrices (fallback = fixed) ----------
    lambda_i <- if (!is.null(v_lambda))
      unwind_lambda(v_lambda[i,], Lambda_mat, reverse = TRUE)
    else Lambda_mat

    B_i      <- if (!is.null(v_B))
      unwind_lambda(v_B[i,], B_mat, reverse = TRUE)
    else B_mat

    K_i      <- if (!is.null(v_K))
      unwind_lambda(v_K[i,], K_mat, reverse = TRUE)
    else K_mat

    G_i      <- if (!is.null(v_G))
      unwind_lambda(v_G[i,], G_mat, reverse = TRUE)
    else G_mat

    ## manifest residual precisions
    eps_prec_i <- exp( eps_log[i,] )          # size n_pars

    ## -------- 2b. Factor-level precision matrix δ_inv -------------------
    delta_inv_i <- matrix(0, n_factors, n_factors)
    lp_delta    <- 0          # prior part
    jac_delta   <- 0          # Jacobian part

    for (fg_id in unique_fg) {

      fg_idx <- which(factor_groups == fg_id)
      d_blk  <- length(fg_idx)
      v_blk  <- delta_vec_list[[as.character(fg_id)]][i, ]

      if (d_blk > 1) {

        ##   multivariate block – v_blk holds the *lower-triangular*
        ##   elements (diag on log scale) of the precision Cholesky.
        L_prec          <- unwind_chol(v_blk, reverse = TRUE)  # precision
        delta_inv_i[fg_idx, fg_idx] <- L_prec

        ##   prior on covariance (inverse-Wishart)  + Jacobian
        cov_blk   <- solve(L_prec)
        lp_delta  <- lp_delta +
          log( robust_diwish(W = cov_blk,
                             v = prior$a_d,
                             S = diag(prior$b_d, d_blk)) )
        jac_delta <- jac_delta +
          calc_log_jac_chol(v_blk)

      } else {  # d_blk == 1  (Gamma on precision)

        log_prec   <- v_blk[1]
        prec_val   <- exp(log_prec)
        delta_inv_i[fg_idx, fg_idx] <- prec_val

        lp_delta  <- lp_delta +
          dgamma(prec_val,
                 shape = prior$a_d,
                 rate  = prior$b_d,
                 log   = TRUE)
        jac_delta <- jac_delta + log_prec
      }
    }

    ## -------- 2c. Priors for other free blocks --------------------------
    lp_theta_mu <- dmvnorm(theta_mu[i, ],
                           mean   = prior$theta_mu_mean,
                           sigma  = diag(prior$theta_mu_var),
                           log    = TRUE)

    lp_lambda <- if (!is.null(v_lambda))
      dmvnorm(v_lambda[i, ],
              mean = rep(0, length(v_lambda[i, ])),
              sigma = diag(prior$lambda_var,
                           length(v_lambda[i, ])),
              log = TRUE) else 0

    lp_B <- if (!is.null(v_B))
      dmvnorm(v_B[i, ],
              mean  = rep(0, length(v_B[i, ])),
              sigma = diag(prior$B_var,
                           length(v_B[i, ])),
              log   = TRUE) else 0

    lp_K <- if (!is.null(v_K))
      dmvnorm(v_K[i, ],
              mean  = rep(0, length(v_K[i, ])),
              sigma = diag(prior$K_var,
                           length(v_K[i, ])),
              log   = TRUE) else 0

    lp_G <- if (!is.null(v_G))
      dmvnorm(v_G[i, ],
              mean  = rep(0, length(v_G[i, ])),
              sigma = diag(prior$G_var,
                           length(v_G[i, ])),
              log   = TRUE) else 0

    ## ε-prior + Jacobian on log scale
    lp_eps <- sum( dgamma(eps_prec_i,
                          shape = prior$a_e,
                          rate  = prior$b_e,
                          log   = TRUE) )
    jac_eps <- sum(eps_log[i, ])

    ## -------- 2d. Group-level log-likelihood ---------------------------
    B0_inv   <- solve(diag(n_factors) - B_i)
    pop_mean <- drop(theta_mu[i, ] +
                       lambda_i %*% B0_inv %*% G_i %*% x_mu +
                       K_i      %*%            x_mu)

    pop_var  <- lambda_i %*% B0_inv %*%
      (G_i %*% x_var %*% t(G_i) +
         solve(delta_inv_i)) %*%
      t(B0_inv) %*% t(lambda_i) +
      K_i %*% x_var %*% t(K_i) +
      diag(1 / eps_prec_i)

    ##  subject-level rows  ×  parameter columns
    alpha_mat <- matrix(proposals[i, ],
                        ncol = n_pars,
                        byrow = TRUE)

    ll_group <- sum( mvtnorm::dmvnorm(alpha_mat,
                                      mean   = pop_mean,
                                      sigma  = pop_var,
                                      log    = TRUE) )

    ## -------- 2e. Combine all pieces -----------------------------------
    log_post_plus_jac[i] <- ll_group     +
      lp_theta_mu  + lp_lambda + lp_B +
      lp_K         + lp_G      +
      lp_eps       +
      lp_delta     +
      jac_eps      + jac_delta
  }

  return(log_post_plus_jac)
}

#' Define Structural Equation Model (SEM) Matrices
#'
#' @description
#' This function helps create the specification matrices (Lambda, B, K, G) for an SEM.
#' It takes a design object, data, factor names, covariate column names, and list-based
#' specifications for the paths to be estimated.
#' The manifest variable names for Lambda_mat and K_mat rows are derived from `sampled_pars(design)`.
#' It validates that covariates are consistent per subject (subject column in `data` must be named "subjects")
#' and includes an aggregated subject-level covariate data frame named `covariates` in the output list.
#' For identifiability, the first parameter listed in `lambda_specs` for each factor is fixed to 1.
#'
#' @param data A data frame containing subject identifiers (in a column named "subjects")
#'   and any covariate columns specified in `covariate_cols`.
#' @param design An emc.design object, as created by the `design()` function.
#'   The parameter names for the SEM are derived from `names(sampled_pars(design))`.
#' @param factor_names Character vector. Names of the latent factors for the SEM.
#' @param covariate_cols Character vector or NULL. Column names in `data` to be used
#'   as covariates for K_mat and G_mat. If NULL, no covariates are processed.
#' @param lambda_specs A list defining `Lambda_mat` (factor loadings).
#'   The list names should be factor names (from `factor_names`), and each element should be a
#'   character vector of parameter names (from `names(sampled_pars(design))`) that load onto that factor.
#'   The first parameter listed for each factor will be fixed to 1 for identifiability.
#'   Example: `list(Factor1 = c("v_Sleft", "a_Eneutral"), Factor2 = "t0")`
#'   Here, `Lambda_mat["v_Sleft", "Factor1"]` would be 1.
#' @param b_specs A list defining `B_mat` (regressions among factors).
#'   List names are outcome factors, elements are character vectors of predictor factors.
#'   Example: `list(Factor2 = "Factor1", Factor3 = c("Factor1", "Factor2"))`
#' @param k_specs A list defining `K_mat` (covariate effects on manifest design parameters).
#'   List names are parameter names (from `names(sampled_pars(design))`), elements are character vectors of covariate names
#'   (must be present in `covariate_cols` and thus in the processed `covariates` data frame).
#'   Example: `list(v_Sleft = "cov1", a_Eneutral = c("cov1", "cov2"))`
#' @param g_specs A list defining `G_mat` (covariate effects on factors).
#'   List names are factor names, elements are character vectors of covariate names.
#'   Example: `list(Factor1 = "cov1", Factor2 = c("cov1", "cov2"))`
#' @param fixed_value Numeric. The value used for fixed paths in the matrices that
#'   are not set to 1 for identifiability or `Inf` for estimation. Default is 0.
#'
#' @return A list (intended to be used as `sem_settings`) containing:
#'   - `Lambda_mat`: The factor loading matrix.
#'   - `B_mat`: The matrix of regressions among factors.
#'   - `K_mat`: The matrix of covariate effects on manifest design parameters.
#'   - `G_mat`: The matrix of covariate effects on factors.
#'   - `par_names`: The manifest design parameter names derived from `sampled_pars(design)`.
#'   - `factor_names`: The provided SEM factor names.
#'   - `covariates`: A data frame with one row per unique subject (ordered by first
#'     appearance in `data[["subjects"]]`) and columns for each validated covariate,
#'     containing the unique subject-level values. Column names are the covariate names.
#'     If no covariates are processed, this will be a data frame with 0 columns and rows for each subject.
#'
#' @export
#'
#' @examples
#' # Create a design object (simplified from design.R example)
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))
#' matchfun_example <- function(d) d$S==d$lR # Example match function
#'
#' example_design_obj <- design(
#'   data = forstmann,
#'   model= LBA,
#'   matchfun=matchfun_example,
#'   formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
#'   contrasts=list(v=list(lM=ADmat)),
#'   constants=c(sv=log(1)),
#'   report_p_vector = FALSE
#' )
#'
#' # SEM Factor names
#' sem_factor_names <- c("Speed", "Caution")
#'
#' # Make a copy of forstmann for example modification
#' forstmann_mod <- forstmann
#' set.seed(123) # for reproducibility
#' subj_trait_values <- stats::setNames(rnorm(length(levels(forstmann_mod$subjects))),
#'                                     levels(forstmann_mod$subjects))
#' forstmann_mod$SubjTrait <- subj_trait_values[forstmann_mod$subjects]
#'
#' my_cov_cols <- c("SubjTrait")
#'
#' lambda_example_specs <- list(
#'   Speed = c("v", "v_lMdiff"), # "v" will be fixed to 1
#'   Caution = c("B", "B_Eneutral", "B_Eaccuracy", "B_lRright", "A") # "B" fixed to 1
#' )
#' b_example_specs <- list(Caution = "Speed")
#' k_example_specs <- list(t0 = "SubjTrait") # "SubjTrait" must be in my_cov_cols
#' g_example_specs <- list(Speed = "SubjTrait")
#'
#' sem_settings_definition <- define_sem_structure(
#'   data = forstmann_mod,
#'   design = example_design_obj,
#'   factor_names = sem_factor_names,
#'   covariate_cols = my_cov_cols,
#'   lambda_specs = lambda_example_specs,
#'   b_specs = b_example_specs,
#'   k_specs = k_example_specs,
#'   g_specs = g_example_specs
#' )
#'
#' print(sem_settings_definition$Lambda_mat)
#' print(sem_settings_definition$B_mat)
#' print(sem_settings_definition$K_mat)
#' print(sem_settings_definition$G_mat)
#' print(head(sem_settings_definition$covariates))
#'
define_sem_structure <- function(data,
                                 design,
                                 factor_names,
                                 covariate_cols = NULL,
                                 lambda_specs = NULL,
                                 b_specs = NULL,
                                 k_specs = NULL,
                                 g_specs = NULL,
                                 fixed_value = 0) {

  subjects_col_internal <- "subjects"
  free_value_internal <- Inf

  if (!inherits(design, "emc.design")) stop("'design' must be an emc.design object.")
  par_names <- names(sampled_pars(design))
  if (length(par_names) == 0) stop("Could not extract parameter names from the design object. Ensure sampled_pars(design) works.")

  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  if (!subjects_col_internal %in% colnames(data)) stop(paste0("Subjects column '", subjects_col_internal, "' not found in data."))
  if (!is.character(factor_names) || length(factor_names) == 0) stop("'factor_names' must be a non-empty character vector for the SEM factors.")
  if (!is.null(covariate_cols) && !is.character(covariate_cols)) stop("'covariate_cols' must be a character vector or NULL.")

  n_pars <- length(par_names)
  n_factors <- length(factor_names)

  processed_covariate_names <- character(0)
  processed_covariates_df <- NULL
  unique_subject_ids <- unique(data[[subjects_col_internal]])

  if (!is.null(covariate_cols) && length(covariate_cols) > 0) {
    subject_cov_list <- vector("list", length(covariate_cols))
    names(subject_cov_list) <- covariate_cols
    temp_processed_cov_names <- character(length(covariate_cols))

    for (j in seq_along(covariate_cols)) {
      cov_name <- covariate_cols[j]
      if (!cov_name %in% colnames(data)) {
        stop(paste0("Specified covariate column '", cov_name, "' not found in data."))
      }

      current_cov_values_for_subjects <- sapply(unique_subject_ids, function(id) {
        subject_data_for_cov <- data[data[[subjects_col_internal]] == id, cov_name]
        if (length(subject_data_for_cov) == 0) {
            stop(paste0("No data for covariate '", cov_name, "' for subject '", id, "'."))
        }
        unique_vals <- unique(subject_data_for_cov)
        if (length(unique_vals) > 1) {
          stop(paste0("Covariate '", cov_name, "' has multiple different values for subject '", id, "'. Covariates must be consistent per subject."))
        }
        return(unique_vals[1])
      })
      subject_cov_list[[cov_name]] <- current_cov_values_for_subjects
      temp_processed_cov_names[j] <- cov_name
    }
    processed_covariate_names <- temp_processed_cov_names
    processed_covariates_df <- as.data.frame(subject_cov_list, stringsAsFactors = FALSE)
  } else {
    processed_covariates_df <- data.frame(matrix(ncol = 0, nrow = length(unique_subject_ids)))
    if (length(unique_subject_ids) > 0 && !is.null(rownames(data))) {
        # This part is tricky if unique_subject_ids isn't directly usable as rownames
        # For simplicity, if rownames are needed, user should ensure `data` has them or handle post-hoc
    }
  }
  n_cov <- length(processed_covariate_names)


  Lambda_mat <- matrix(fixed_value, nrow = n_pars, ncol = n_factors,
                       dimnames = list(par_names, factor_names))
  B_mat <- matrix(fixed_value, nrow = n_factors, ncol = n_factors,
                  dimnames = list(factor_names, factor_names))

  if (n_cov > 0) {
    K_mat <- matrix(fixed_value, nrow = n_pars, ncol = n_cov,
                    dimnames = list(par_names, processed_covariate_names))
    G_mat <- matrix(fixed_value, nrow = n_factors, ncol = n_cov,
                    dimnames = list(factor_names, processed_covariate_names))
  } else {
    K_mat <- matrix(fixed_value, nrow = n_pars, ncol = 0, dimnames = list(par_names, NULL))
    G_mat <- matrix(fixed_value, nrow = n_factors, ncol = 0, dimnames = list(factor_names, NULL))
  }

  if (!is.null(lambda_specs)) {
    if (!is.list(lambda_specs)) stop("'lambda_specs' must be a list.")
    for (f_name in names(lambda_specs)) {
      if (!f_name %in% factor_names) stop(paste0("Factor '", f_name, "' in lambda_specs not in SEM factor_names."))
      p_names_for_f <- lambda_specs[[f_name]]
      if (!is.character(p_names_for_f) || length(p_names_for_f) == 0) {
        stop(paste0("Values in lambda_specs for factor '", f_name, "' must be a non-empty character vector of design parameter names."))
      }
      
      # Identifiability constraint: first parameter fixed to 1
      first_p_name_spec <- p_names_for_f[1]
      if (!first_p_name_spec %in% par_names) stop(paste0("Parameter '", first_p_name_spec, "' for factor '", f_name, "' in lambda_specs (for identifiability) not in names derived from sampled_pars(design)."))
      Lambda_mat[first_p_name_spec, f_name] <- 1
      
      # Other specified parameters set to free
      if (length(p_names_for_f) > 1) {
        for (p_name_spec in p_names_for_f[-1]) {
          if (!p_name_spec %in% par_names) stop(paste0("Parameter '", p_name_spec, "' for factor '", f_name, "' in lambda_specs not in names derived from sampled_pars(design)."))
          Lambda_mat[p_name_spec, f_name] <- free_value_internal
        }
      }
    }
  }

  if (!is.null(b_specs)) {
    if (!is.list(b_specs)) stop("'b_specs' must be a list.")
    for (to_f_name in names(b_specs)) {
      if (!to_f_name %in% factor_names) stop(paste0("Factor '", to_f_name, "' (as outcome) in b_specs not in SEM factor_names."))
      from_f_names_for_f <- b_specs[[to_f_name]]
      if (!is.character(from_f_names_for_f)) stop(paste0("Values in b_specs for factor '", to_f_name, "' must be character vectors of SEM factor_names."))
      for (from_f_name in from_f_names_for_f) {
        if (!from_f_name %in% factor_names) stop(paste0("Factor '", from_f_name, "' (as predictor) for '",to_f_name,"' in b_specs not in SEM factor_names."))
        if (from_f_name == to_f_name) stop(paste0("Factor '", from_f_name, "' cannot predict itself in b_specs (no self-loops)."))
        B_mat[to_f_name, from_f_name] <- free_value_internal
      }
    }
  }

  if (n_cov > 0 && !is.null(k_specs)) {
    if (!is.list(k_specs)) stop("'k_specs' must be a list.")
    for (p_name_spec in names(k_specs)) {
      if (!p_name_spec %in% par_names) stop(paste0("Parameter '", p_name_spec, "' in k_specs not in names derived from sampled_pars(design)."))
      cov_names_for_p <- k_specs[[p_name_spec]]
      if (!is.character(cov_names_for_p)) stop(paste0("Values in k_specs for parameter '", p_name_spec, "' must be character vectors of covariate names."))
      for (cov_name in cov_names_for_p) {
        if (!cov_name %in% processed_covariate_names) stop(paste0("Covariate '", cov_name, "' for parameter '", p_name_spec,"' in k_specs not in processed_covariate_names (i.e., not in covariate_cols or not found in data)."))
        K_mat[p_name_spec, cov_name] <- free_value_internal
      }
    }
  }

  if (n_cov > 0 && !is.null(g_specs)) {
    if (!is.list(g_specs)) stop("'g_specs' must be a list.")
    for (f_name in names(g_specs)) {
      if (!f_name %in% factor_names) stop(paste0("Factor '", f_name, "' in g_specs not in SEM factor_names."))
      cov_names_for_f <- g_specs[[f_name]]
      if (!is.character(cov_names_for_f)) stop(paste0("Values in g_specs for factor '", f_name, "' must be character vectors of covariate names."))
      for (cov_name in cov_names_for_f) {
        if (!cov_name %in% processed_covariate_names) stop(paste0("Covariate '", cov_name, "' for factor '", f_name,"' in g_specs not in processed_covariate_names (i.e., not in covariate_cols or not found in data)."))
        G_mat[f_name, cov_name] <- free_value_internal
      }
    }
  }

  return(list(
    Lambda_mat = Lambda_mat,
    B_mat = B_mat,
    K_mat = K_mat,
    G_mat = G_mat,
    par_names = par_names,
    factor_names = factor_names,
    covariates = processed_covariates_df
  ))
}
