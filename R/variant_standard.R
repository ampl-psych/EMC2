
sample_store_standard <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance,...) {
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  par_names <- par_names[!is_nuisance]
  n_pars <- length(par_names)
  betas <- list(...)$betas
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    a_half = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL))
  )
  if(!is.null(betas)){
    samples$betas <- array(NA_real_,dim = c(length(betas), iters),dimnames = list(betas, NULL))
  }
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

add_info_standard <- function(sampler, prior = NULL, ...){
  sampler$prior <- get_prior_standard(prior, sum(!sampler$nuisance), sample = F, betas = list(...)$betas)
  return(sampler)
}


get_prior_standard <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL, betas = NULL){
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

  if(!is.null(betas) > 0 & is.null(prior$beta_mean)){
    prior$beta_mean <- rep(0, length(betas))
  }
  if(!is.null(betas) & is.null(prior$beta_var)){
    prior$beta_var <- diag(rep(1, length(betas)))
  }

  if(is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- diag(rep(1, n_pars))
  }
  if(is.null(prior$v)){
    prior$v <- 2
  }
  if(is.null(prior$A)){
    prior$A <- rep(.3, n_pars)
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- ginv(prior$theta_mu_var) #Inverse of the matrix
  attr(prior, "type") <- "standard"
  out <- prior
  if(sample){
    par_names <- names(sampled_pars(design, doMap = F))
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
        attempt <- tryCatch({
          vars[,,i] <- riwish(prior$v + n_pars - 1, 2 * prior$v * diag(1 / a_half))
        },error=function(e) e, warning=function(w) w)
        if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
          sample_idx <- sample(1:(i-1),1)
          vars[,,i] <- vars[,,sample_idx]
        }
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

get_startpoints_standard <- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!pmwgs$nuisance)
  if (is.null(start_mu)) start_mu <- rmvnorm(1, mean = pmwgs$prior$theta_mu_mean, sigma = pmwgs$prior$theta_mu_var)
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(n_pars * 3,diag(n_pars))
  start_a_half <- 1 / rgamma(n = n_pars, shape = 2, rate = 1)
  if(!is.null(pmwgs$betas)){
    betas <- rmvnorm(1, mean = pmwgs$prior$theta_mu_mean, sigma = pmwgs$prior$theta_mu_var)
  } else{
    betas <- NULL
  }
  return(list(tmu = start_mu, tvar = start_var, tvinv = ginv(start_var), a_half = start_a_half,
              beta = betas))
}

get_group_level_standard <- function(parameters, s){
  # This function is modified for other versions
  mu <- parameters$tmu
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

fill_samples_standard <- function(samples, group_level, proposals,  j = 1, n_pars){
  samples$a_half[, j] <- group_level$a_half
  samples$beta <- group_level$beta
  samples$last_theta_var_inv <- group_level$tvinv
  samples <- fill_samples_base(samples, group_level, proposals, j = j, n_pars)
  return(samples)
}

gibbs_step_standard <- function(sampler, alpha) {
  #
  # alpha:         (p x n) subject-level parameter matrix
  #
  # sampler$design_mats:   list of length p, each (n x m_k)
  # sampler$beta_indices:  list of length p, each integer indices in Beta for row k
  # sampler$is_blocked:    length p logical; is_blocked[k] = TRUE => param k is in some covariance block
  # sampler$par_group:     length p integer group label for each param (only matters if is_blocked=TRUE)
  #
  # sampler$prior contains:
  #   - theta_mu_mean, theta_mu_invar  => prior mean & precision for Beta (length M and MxM)
  #   - v, A => IW/Gamma hyperparams
  #
  #
  # This function returns an updated (beta, tvar, tvinv, a_half),
  # doing partial block updates of tvar:
  #   => If is_blocked[k] is TRUE, param k shares off-diagonal cov with its group in par_group.
  #   => If is_blocked[k] is FALSE, param k is updated as diagonal-only with a Gamma draw.


  design_mats   <- sampler$design_mats
  beta_indices  <- sampler$beta_indices
  is_blocked    <- sampler$is_blocked
  par_group     <- sampler$par_group
  prior         <- sampler$prior

  last    <- last_sample_standard(sampler$samples)
  beta    <- last$beta      # (M)
  tvar    <- last$tvar      # (p x p)
  tvinv   <- last$tvinv     # (p x p)
  a_half  <- last$a_half    # (p)

  p <- nrow(alpha)          # number of parameters
  n <- ncol(alpha)          # number of subjects

  ##------------------------ 1) Update Beta (jointly) ------------------------##
  M <- length(beta)         # total dimension of Beta
  prec_data <- matrix(0, M, M)
  mean_data <- numeric(M)

  # Build data-based precision & mean
  for (i in seq_len(n)) {
    # M_i = (p x M)
    M_i <- matrix(0, nrow=p, ncol=M)
    for (k in seq_len(p)) {
      cols_k <- beta_indices[[k]]
      x_ik   <- design_mats[[k]][i, , drop=FALSE]
      M_i[k, cols_k] <- x_ik
    }
    alpha_i <- alpha[, i, drop=FALSE]

    prec_data <- prec_data + crossprod(M_i, tvinv %*% M_i)
    mean_data <- mean_data + crossprod(M_i, tvinv %*% alpha_i)
  }

  prec_post <- prior$theta_mu_invar + prec_data
  cov_post  <- solve(prec_post)
  mean_post <- cov_post %*% ( prior$theta_mu_invar %*% prior$theta_mu_mean + mean_data )

  # Sample new Beta
  L         <- t(chol(cov_post))
  z         <- rnorm(M)
  beta_new  <- as.vector(mean_post + L %*% z)


  ##---------------------- 2) Compute residuals for tvar update -------------##
  resid <- matrix(0, nrow=p, ncol=n)
  for (i in seq_len(n)) {
    alpha_i <- alpha[, i, drop=FALSE]

    M_i <- matrix(0, nrow=p, ncol=M)
    for (k in seq_len(p)) {
      cols_k <- beta_indices[[k]]
      x_ik   <- design_mats[[k]][i, , drop=FALSE]
      M_i[k, cols_k] <- x_ik
    }

    mu_i <- M_i %*% beta_new
    resid[, i] <- alpha_i[,1] - mu_i
  }


  ##---------------- 3) Build a brand-new tvar_new partially blocked ---------##
  tvar_new <- matrix(0, p, p)

  # Identify the subset of parameters that are blocked vs unblocked
  blocked_idx   <- which(is_blocked)
  unblocked_idx <- which(!is_blocked)

  #--- 3a) Multi-parameter block IW for each group among blocked -------------
  if (length(blocked_idx) > 0) {
    block_groups <- unique(par_group[blocked_idx])
    for (g in block_groups) {
      group_idx <- blocked_idx[ par_group[blocked_idx] == g ]
      d <- length(group_idx)  # dimension of this block

      resid_block <- resid[group_idx, , drop=FALSE]
      cov_block   <- resid_block %*% t(resid_block)

      B_half_block <- 2 * prior$v * diag(1 / a_half[group_idx], d) + cov_block
      df_block <- prior$v + d - 1 + n

      Sigma_block <- riwish(df_block, B_half_block)

      # place the block into tvar_new
      tvar_new[group_idx, group_idx] <- Sigma_block
    }
  }

  #--- 3b) Vectorized gamma update for unblocked (diagonal-only) -------------
  if (length(unblocked_idx) > 0) {
    # sum of squares for each unblocked row
    # rowSums( (resid[k,])^2 ) is the total squared residual for param k
    sum_sq_vec <- rowSums(resid[unblocked_idx, , drop=FALSE]^2)

    # shape = v/2 + n/2
    shape_vec  <- prior$v/2 + n/2
    # rate = prior$v / a_half[k] + 0.5 * sum_sq
    rate_vec   <- prior$v / a_half[unblocked_idx] + 0.5 * sum_sq_vec

    # Draw the precision for each unblocked param
    tvinv_diag <- rgamma(length(unblocked_idx), shape=shape_vec, rate=rate_vec)
    tvar_diag  <- 1 / tvinv_diag

    # Fill these diagonal entries into tvar_new
    tvar_new[cbind(unblocked_idx, unblocked_idx)] <- tvar_diag
  }

  ##------------------ 4) Invert tvar_new to get tvinv_new -------------------##
  tvinv_new <- solve(tvar_new)


  ##--------- 5) Update mixing weights a_half for all p parameters ----------##
  #
  # We'll allow the shape to depend on the dimension of each parameter's block:
  #   shape_k = (v + block_dim[k]) / 2
  # Where block_dim[k] = 1 if unblocked, or = size of that block if blocked.
  block_dim <- integer(p)
  block_dim[unblocked_idx] <- 1
  # For each block group g, set block_dim to the block size:
  if (length(blocked_idx) > 0) {
    block_groups <- unique(par_group[blocked_idx])
    for (g in block_groups) {
      group_idx <- blocked_idx[ par_group[blocked_idx] == g ]
      block_dim[group_idx] <- length(group_idx)
    }
  }

  shape_vec <- (prior$v + block_dim) / 2
  rate_vec  <- prior$v * diag(tvinv_new) + 1 / (prior$A^2)

  a_half_new <- 1 / rgamma(p, shape=shape_vec, rate=rate_vec)


  ##---------------------- 6) Return updated values -------------------------##
  out <- list(
    beta   = beta_new,
    tmu = tmu,
    tvar   = tvar_new,
    tvinv  = tvinv_new,
    a_half = a_half_new,
    alpha = alpha
  )
  return(out)
}


# For now don't consider betas in conditionals
get_conditionals_standard <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  pts2_unwound <- apply(samples$theta_var[idx,idx,],3,unwind)
  all_samples <- rbind(samples$alpha[idx, s,],samples$theta_mu[idx,],pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- stats::cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration], unwind(samples$theta_var[idx,idx,iteration])))
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

filtered_samples_standard <- function(sampler, filter, ...){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter, drop = F],
    theta_var = sampler$samples$theta_var[, , filter, drop = F],
    alpha = sampler$samples$alpha[, , filter, drop = F],
    iteration = length(filter)
  )
}


calc_log_jac_chol <- function(x) {
  # x: a vector of length n*(n+1)/2 representing the lower-triangular
  # matrix with its diagonal stored on the log-scale.
  # Determine n from the length of x.
  n <- floor(sqrt(2 * length(x) + 0.25) - 0.5)

  # Reconstruct the lower-triangular matrix in "parameter space"
  mat <- matrix(NA, n, n)
  mat[lower.tri(mat, diag = TRUE)] <- x

  # The diagonal of the actual Cholesky factor is exp(x_diag),
  # so log(diag(L)) = x_diag.
  # The log-Jacobian is then:
  #   sum_{i=1}^{n} (n - i + 2) * x_diag[i]
  log_jac <- sum( (n - seq_len(n) + 2) * diag(mat) )

  return(log_jac)
}

unwind_chol <- function(x,reverse=FALSE) {

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

# bridge_sampling ---------------------------------------------------------
bridge_group_and_prior_and_jac_blocked <- function(proposals_group, proposals_list, info){
  prior <- info$prior
  par_groups <- info$par_groups
  has_cov <- par_groups %in% which(table(par_groups) > 1)
  proposals <- do.call(cbind, proposals_list)
  theta_mu <- proposals_group[,1:info$n_pars]
  theta_a <- proposals_group[,(info$n_pars + 1):(2*info$n_pars)]
  theta_var1 <- proposals_group[,(2*info$n_pars + 1):(2*info$n_pars + sum(!has_cov))]
  theta_var2 <- proposals_group[,(2*info$n_pars + sum(!has_cov) + 1):(2*info$n_pars + sum(!has_cov) + (sum(has_cov) * (sum(has_cov) +1))/2)]

  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)
  var_curr <- matrix(0, nrow = info$n_pars, ncol = info$n_pars)

  for(i in 1:n_iter){ # these unfortunately can't be vectorized
    # prior_delta2 <- log(robust_diwish(solve(delta2_curr), v=prior$a_d, S = diag(prior$b_d, sum(!info$is_structured))))
    var2curr <- unwind_chol(theta_var2[i,], reverse = T)
    var_curr[!has_cov, !has_cov] <- diag(exp(theta_var1[i,]), sum(!has_cov))
    var_curr[has_cov, has_cov] <- var2curr
    proposals_curr <- matrix(proposals[i,], ncol = info$n_pars, byrow = T)
    group_ll <- sum(dmvnorm(proposals_curr, theta_mu[i,], var_curr, log = T))
    prior_var1 <- sum(logdinvGamma(exp(theta_var1[i,]), shape = prior$v/2, rate = prior$v/exp(theta_a[i,!has_cov])))
    prior_var2 <- log(robust_diwish(var2curr, v=prior$v+ sum(has_cov)-1, S = 2*prior$v*diag(1/theta_a[i,has_cov])))
    prior_a <- sum(logdinvGamma(exp(theta_a[i,]), shape = 1/2,rate=1/(prior$A^2)))
    jac_var2 <- calc_log_jac_chol(theta_var2[i, ])
    sum_out[i] <- group_ll + prior_var1 + prior_var2 + jac_var2 + prior_a
  }
  prior_mu <- dmvnorm(theta_mu, mean = prior$theta_mu_mean, sigma = prior$theta_mu_var, log =T)
  jac_var1 <- rowSums(theta_var1)
  jac_a <- rowSums(theta_a)
  return(sum_out + prior_mu + jac_var1 + jac_a)
}


bridge_add_info_blocked <- function(info, samples){
  par_groups <- samples$par_groups
  has_cov <- par_groups %in% which(table(par_groups) > 1)
  info$par_groups <- par_groups
  info$group_idx <- (samples$n_pars*samples$n_subjects + 1):(samples$n_pars*samples$n_subjects + 2*samples$n_pars +
                                                               sum(!has_cov) + (sum(has_cov) * (sum(has_cov) +1))/2)
  return(info)
}

bridge_add_group_blocked <- function(all_samples, samples, idx){
  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[,idx]))
  all_samples <- cbind(all_samples, t(log(samples$samples$a_half[,idx])))
  par_groups <- samples$par_groups
  has_cov <- par_groups %in% which(table(par_groups) > 1)
  all_samples <- cbind(all_samples, t(log(matrix(apply(samples$samples$theta_var[!has_cov,!has_cov,idx, drop = F], 3, diag), ncol = nrow(all_samples)))))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$theta_var[has_cov,has_cov,idx, drop = F], 3, unwind_chol), ncol = nrow(all_samples))))
  return(all_samples)
}


# for IC ------------------------------------------------------------------

group__IC_standard <- function(emc, stage="sample",filter=NULL){
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                       return_mcmc = FALSE, merge_chains = TRUE)
  theta_mu <- get_pars(emc, selection = "mu", stage = stage, filter = filter,
                          return_mcmc = FALSE, merge_chains = TRUE)
  theta_var <- get_pars(emc, selection = "Sigma", stage = stage, filter = filter,
                           return_mcmc = FALSE, merge_chains = TRUE, remove_constants = F)
  mean_alpha <- apply(alpha, 1:2, mean)
  mean_mu <- rowMeans(theta_mu)
  mean_var <- apply(theta_var, 1:2, mean)

  N <- ncol(theta_mu)
  lls <- numeric(N)
  for(i in 1:N){
    lls[i] <- sum(dmvnorm(t(alpha[,,i]), theta_mu[,i], theta_var[,,i], log = T))
  }
  minD <- -2*max(lls)
  mean_ll <- mean(lls)
  mean_pars_ll <-  sum(dmvnorm(t(mean_alpha), mean_mu, mean_var, log = TRUE))
  Dmean <- -2*mean_pars_ll
  return(list(mean_ll = mean_ll, Dmean = Dmean,
              minD = minD))
}

