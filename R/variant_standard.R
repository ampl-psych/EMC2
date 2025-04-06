
sample_store_standard <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance,...) {
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  par_names <- par_names[!is_nuisance]
  n_pars <- length(par_names)
  betas <- get_betas(list(...)$group_design)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    a_half = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL))
  )
  if(!is.null(betas)){
    samples$beta <- array(NA_real_,dim = c(length(betas), iters),dimnames = list(betas, NULL))
  }
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

add_info_standard <- function(sampler, prior = NULL, ...){
  n_pars <-sum(!sampler$nuisance)
  betas <- get_betas(list(...)$group_design)
  if(!is.null(betas)){
    sampler$betas <- betas
    sampler$group_designs <- list(...)$group_design
  }
  sampler$par_group <- list(...)$par_groups
  sampler$is_blocked <- sampler$par_group %in% which(table(sampler$par_group) > 1)
  sampler$prior <- get_prior_standard(prior, n_pars, sample = F, betas = betas)
  return(sampler)
}

get_betas <- function(group_design){
  if(is.null(group_design)) return(NULL)
  betas <- unlist(Map(function(name, values) paste0(name, "_", values), names(group_design),
                      lapply(group_design, function(x) colnames(x))))
  names(betas) <- NULL
  return(betas)
}

get_prior_standard <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL, betas = NULL,
                               par_groups = NULL){
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

  if(!is.null(betas) & is.null(prior$beta_mean)){
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
    if(selection %in% "beta"){
      beta <- t(mvtnorm::rmvnorm(N, mean = prior$beta_mean,
                               sigma = prior$beta_var))
      rownames(beta) <- betas
      samples$beta <- beta
    }
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
      if(is.null(par_groups)) par_groups <- rep(1, n_pars)
      constraintMat <- matrix(0, n_pars, n_pars)
      for(i in 1:length(unique(par_groups))){
        idx <- par_groups == i
        constraintMat[idx, idx] <- Inf
      }
      vars <- constrain_lambda(vars, constraintMat)
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
    betas <- rmvnorm(1, mean = pmwgs$prior$beta_mean, sigma = pmwgs$prior$beta_var)
  } else{
    betas <- NULL
  }
  subj_mu <- do.call(cbind, rep(list(c(start_mu)), pmwgs$n_subjects))
  return(list(tmu = start_mu, tvar = start_var, tvinv = ginv(start_var),
              a_half = start_a_half, subj_mu = subj_mu, beta = betas))
}

get_group_level_standard <- function(parameters, s){
  # This function is modified for other versions
  mu <- parameters$subj_mu[,s]
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

fill_samples_standard <- function(samples, group_level, proposals,  j = 1, n_pars){
  samples$a_half[, j] <- group_level$a_half
  if(!is.null(group_level$beta)){
    samples$beta[,j] <- group_level$beta
  }
  samples$last_theta_var_inv <- group_level$tvinv
  samples <- fill_samples_base(samples, group_level, proposals, j = j, n_pars)
  return(samples)
}

add_intercepts <- function(par_names, group_designs = NULL, n_subjects) {
  #
  # par_names:   character vector of parameter names (length p)
  # group_designs: named list of existing design matrices, e.g. $"param1" => (n x m1), etc.
  #                    some par_names might not appear here.
  # n_subjects:    integer, number of subjects (rows).
  #
  # Returns a new list, same length as par_names, each entry an (n_subjects x m_k) matrix
  # that includes an intercept column as the FIRST column. If param was not in group_designs,
  # we create a single column of 1’s of dimension (n_subjects x 1).

  out_list <- vector("list", length(par_names))
  names(out_list) <- par_names

  for (k in seq_along(par_names)) {
    pname <- par_names[k]
    if (!is.null(group_designs[[pname]])) {
      # existing design => prepend an intercept column
      mat_k <- group_designs[[pname]]
      # cbind(1, mat_k)
      out_list[[k]] <- cbind(1, mat_k)
    } else {
      # no existing design => just a column of 1's
      out_list[[k]] <- matrix(1, nrow = n_subjects, ncol = 1)
    }
  }

  return(out_list)
}

gibbs_step_standard <- function(sampler, alpha) {
  #
  # alpha:  (p x n) subject-level parameter matrix
  #
  # sampler$group_designs:   list of length p, each (n x m_k)
  # sampler$is_blocked:    length p logical; is_blocked[k] = TRUE => param k is in some covariance block
  # sampler$par_group:     length p integer group labels (only used if is_blocked=TRUE)
  #
  # sampler$prior contains at least:
  #   - theta_mu_mean (length p) and theta_mu_invar (p x p)
  #   - beta_mean (length of "all slopes") and beta_var (matching dimension)
  #   - v, A => hyperparams for the IW/Gamma updates
  #
  # last_sample_standard(...) must provide:
  #   - last$tmu (p-vector) of group means
  #   - last$beta (the additional slopes, if any)
  #   - last$tvar, last$tvinv (p x p)
  #   - last$a_half (p)
  #
  # The function returns an updated (tmu, beta, tvar, tvinv, a_half, alpha),
  # doing partial block updates for the covariance tvar, and a single joint
  # update for [mu + beta] in one big vector.


  group_designs   <- sampler$group_designs
  is_blocked    <- sampler$is_blocked
  par_group     <- sampler$par_group
  prior         <- sampler$prior



  last    <- last_sample_standard(sampler$samples)
  tmu     <- last$tmu        # group means, dimension p
  beta    <- last$beta       # additional slopes
  tvar    <- last$tvar       # (p x p)
  tvinv   <- last$tvinv
  a_half  <- last$a_half     # length p

  p <- nrow(alpha)           # must match length(tmu)
  n <- ncol(alpha)           # number of subjects

  # Some backwards compatibility
  if(is.null(is_blocked)){
    is_blocked <- rep(T, p)
  }
  if(is.null(par_group)){
    par_group <- rep(1, p)
  }



  group_designs <- add_intercepts(names(tmu), group_designs, n)
  ##--------------------------------------------------
  ## 1) Combine tmu & beta into one big vector
  ##--------------------------------------------------
  # "mu" is length p, "beta" is length(...). So total dimension:
  M <- length(tmu) + length(beta)

  # We will build a single prior mean vector 'prior_mean' (length M)
  # and a single prior precision 'prior_invar' (MxM).
  prior_mean    <- numeric(M)
  prior_invar   <- matrix(0, nrow=M, ncol=M)  # we'll fill the diagonal blocks

  # Decide which indices in {1..M} correspond to mu, which to beta
  # That set of columns is treated as the group means (the intercept dimension).
  mean_index <- unlist(lapply(group_designs, \(x) c(TRUE, rep(FALSE, ncol(x) - 1))))

  # Fill prior_mean for the group means vs. slopes
  # e.g. prior$theta_mu_mean is length p, prior$beta_mean is length(M - p)
  prior_mean[mean_index]  <- prior$theta_mu_mean
  prior_mean[!mean_index] <- prior$beta_mean

  # Fill prior_invar for the group means (block) vs. slopes (block)
  # e.g. prior$theta_mu_invar is (p x p), prior$beta_var is (some dim x some dim)
  # Fill the diagonal blocks with the correct inverses
  # group means block:
  prior_invar[mean_index, mean_index] <- prior$theta_mu_invar
  # slopes block:
  if(any(!mean_index)){
    beta_invar <- solve(prior$beta_var)
    prior_invar[!mean_index, !mean_index] <- beta_invar
  }


  ##--------------------------------------------------
  ## 2) Build data-based precision & mean, for the
  ##    big vector of dimension M
  ##--------------------------------------------------
  prec_data <- matrix(0, M, M)
  mean_data <- numeric(M)

  # For each subject i, we build an (p x M) matrix that maps [mu+beta] -> predicted alpha_i
  for (i in seq_len(n)) {
    par_idx <- 0
    M_i <- matrix(0, nrow=p, ncol=M)
    for (k in seq_len(p)) {
      x_ik   <- group_designs[[k]][i, , drop=FALSE]
      # place them in row k:
      M_i[k, par_idx + 1:ncol(group_designs[[k]])] <- x_ik
      par_idx <- par_idx + ncol(group_designs[[k]])
    }
    # Accumulate
    prec_data <- prec_data + crossprod(M_i, tvinv %*% M_i)
    mean_data <- mean_data + crossprod(M_i, tvinv %*% alpha[, i, drop=FALSE])
  }

  prec_post <- prior_invar + prec_data
  cov_post  <- solve(prec_post)
  mean_post <- cov_post %*% (prior_invar %*% prior_mean + mean_data)

  # Draw new big vector
  L          <- t(chol(cov_post))
  z          <- rnorm(M)
  mean_new   <- as.vector(mean_post + L %*% z)

  ##--------------------------------------------------
  ## 3) Compute residuals = alpha - (X * [mu+beta])
  ##--------------------------------------------------
  resid <- matrix(0, nrow=p, ncol=n)
  subj_mu <- matrix(0, nrow=p, ncol=n)
  for (i in seq_len(n)) {
    par_idx <- 0
    M_i <- matrix(0, nrow=p, ncol=M)
    for (k in seq_len(p)) {
      x_ik   <- group_designs[[k]][i, , drop=FALSE]
      M_i[k, par_idx + 1:ncol(group_designs[[k]])] <- x_ik
      par_idx <- par_idx + ncol(group_designs[[k]])
    }

    mu_i <- M_i %*% mean_new
    subj_mu[,i] <- mu_i
    resid[, i] <- alpha[,i] - mu_i
  }

  ##--------------------------------------------------
  ## 4) Partial-block update of tvar_new
  ##--------------------------------------------------
  tvar_new <- matrix(0, p, p)

  blocked_idx   <- which(is_blocked)
  unblocked_idx <- which(!is_blocked)

  # 4a) For each block among the "blocked" subset, do an inverse-Wishart
  if (length(blocked_idx) > 0) {
    block_groups <- unique(par_group[blocked_idx])
    for (g in block_groups) {
      group_idx <- blocked_idx[ par_group[blocked_idx] == g ]
      d <- length(group_idx)
      resid_block <- resid[group_idx, , drop=FALSE]
      cov_block   <- resid_block %*% t(resid_block)

      B_half_block <- 2 * prior$v * diag(1 / a_half[group_idx], d) + cov_block
      df_block     <- prior$v + d - 1 + n

      Sigma_block  <- riwish(df_block, B_half_block)
      tvar_new[group_idx, group_idx] <- Sigma_block
    }
  }

  # 4b) For each "unblocked" param, do a vectorized Gamma for diagonal-only
  if (length(unblocked_idx) > 0) {
    sum_sq_vec <- rowSums(resid[unblocked_idx, , drop=FALSE]^2)
    # shape = v/2 + n/2
    shape_vec <- prior$v/2 + n/2
    # rate = v/a_half + 0.5 * sum_sq
    rate_vec  <- prior$v / a_half[unblocked_idx] + 0.5 * sum_sq_vec

    tvinv_diag <- rgamma(length(unblocked_idx), shape=shape_vec, rate=rate_vec)
    tvar_diag  <- 1 / tvinv_diag

    tvar_new[cbind(unblocked_idx, unblocked_idx)] <- tvar_diag
  }

  # Invert
  tvinv_new <- solve(tvar_new)

  ##--------------------------------------------------
  ## 5) Update a_half
  ##--------------------------------------------------
  # shape depends on block dimension: (v + block_dim[k]) / 2
  block_dim <- integer(p)
  block_dim[unblocked_idx] <- 1
  if (length(blocked_idx) > 0) {
    block_groups <- unique(par_group[blocked_idx])
    for (g in block_groups) {
      group_idx <- blocked_idx[ par_group[blocked_idx] == g ]
      block_dim[group_idx] <- length(group_idx)
    }
  }

  shape_vec <- (prior$v + block_dim) / 2
  rate_vec  <- prior$v * diag(tvinv_new) + 1/(prior$A^2)
  a_half_new <- 1 / rgamma(p, shape=shape_vec, rate=rate_vec)

  ##--------------------------------------------------
  ## 6) Separate out the updated tmu & beta
  ##--------------------------------------------------
  # mean_index = those columns for 'group means' => new tmu
  # !mean_index = those columns for slopes => new beta
  tmu_new  <- mean_new[ mean_index ]
  if(any(!mean_index)){
    beta_new <- mean_new[ !mean_index ]
  } else{
    beta_new <- NULL
  }
  ##--------------------------------------------------
  ## 7) Return updated values
  ##--------------------------------------------------
  out <- list(
    tmu    = tmu_new,          # updated group means
    beta   = beta_new,         # updated slopes
    tvar   = tvar_new,
    tvinv  = tvinv_new,
    a_half = a_half_new,
    subj_mu = subj_mu,
    alpha  = alpha
  )
  return(out)
}


last_sample_standard <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],
    beta = store$beta[,store$idx],
    tvar = store$theta_var[, , store$idx],
    tvinv = store$last_theta_var_inv,
    a_half = store$a_half[, store$idx]
  )
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
bridge_group_and_prior_and_jac_standard <- function(
    proposals_group,   # (n_iter x [theta_mu, theta_a, var1, var2, (beta)?])
    proposals_list,    # list of length n_iter, each item is the subject-level alpha draws
    info
) {
  # 'info' should contain:
  #   info$n_pars        = p        (# of rows in alpha)
  #   info$n_subjects    = n_subj   (# of subjects)
  #   info$is_blocked    = logical p  (which rows are blocked)
  #   info$par_group     = integer p  (group IDs for blocked rows)
  #   info$group_designs   = list of length p, each (n_subj x m_k)
  #   info$prior:
  #       - $theta_mu_mean (p), $theta_mu_var (pxp)
  #       - $beta_mean (maybe NULL or length=?), $beta_var (matching dimension)
  #       - $v, $A (IW & gamma hyperparams)
  #
  # The 'proposals_group' columns are arranged:
  #   1) theta_mu (p columns)
  #   2) theta_a  (p columns) -> log(a_half)
  #   3) theta_var1 (sum(!has_cov) columns) -> log of unblocked diagonal
  #   4) theta_var2 ((d*(d+1))/2 columns)    -> cholesky for the blocked subset
  #   5) (optional) betas  (length(prior$beta_mean) columns) -> could be log or direct
  #
  # 'proposals_list' is a list of length n_iter, each item is an (n_subj * p)-vector or
  #   something that we can reshape into (n_subj x p), i.e. the subject-level alphas.

  par_group <- info$par_group
  has_cov   <- info$is_blocked
  block_groups <- unique(par_group[has_cov])
  prior     <- info$prior

  p         <- info$n_pars        # dimension of each alpha_i
  n_subj    <- info$n_subjects
  B <- length(info$betas)
  M <- p + B
  group_designs <- add_intercepts(info$par_names, info$group_designs, n_subj)
  mean_index <- unlist(lapply(group_designs, \(x) c(TRUE, rep(FALSE, ncol(x) - 1))))

  # Extract columns:
  theta_mu   <- proposals_group[, seq_len(p),   drop=FALSE]  # (n_iter x p)
  theta_a    <- proposals_group[, p + seq_len(p),    drop=FALSE]  # (n_iter x p)
  if(B > 0){
    theta_beta <- proposals_group[,(2*p + 1):(2*p+B), drop = FALSE]
  }
  if(any(!has_cov)){
    theta_var1 <- proposals_group[, (2*p + 1 + B) : (2*p + B + sum(!has_cov)),       drop=FALSE]  # (n_iter x sum(!has_cov))
  }
  if(any(has_cov)){
    min_idx <- (2*p + B + sum(!has_cov))
    theta_var2_list <- list()
    for(block in block_groups){
      cur_idx <- has_cov[par_group == block]
      max_idx <- min_idx + (sum(cur_idx)*(sum(cur_idx)+1)) / 2
      theta_var2_list[[block]] <- proposals_group[, (min_idx + 1) :
                                      max_idx,drop=FALSE]  # (n_iter x (d*(d+1))/2)
      min_idx <- max_idx
    }

  }
  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)

  # Precompute the prior on mu, beta outside the loop (vectorized)
  prior_mu_log <- dmvnorm(theta_mu,
                          mean  = prior$theta_mu_mean,
                          sigma = prior$theta_mu_var,
                          log   = TRUE)
  prior_beta_log <- 0
  if (B > 0) {
    prior_beta_log <- dmvnorm(theta_beta,
                              mean  = prior$beta_mean,
                              sigma = prior$beta_var,
                              log   = TRUE)
  }

  # For the Jacobian of var1 => rowSums(theta_var1)
  # For the Jacobian of a    => rowSums(theta_a)
  if(any(!has_cov)){
    jac_var1 <- rowSums(theta_var1)
  } else{
    jac_var1 <- 0
  }
  jac_a_vec    <- rowSums(theta_a)

  #--- Main loop over MCMC draws ---
  # We reconstruct var_curr & compute the group-likelihood, plus prior on var1/var2/a
  # then we store in sum_out[i].  After the loop, we add prior_mu_log + prior_beta_log + jac_var1 + jac_a.
  for (i in seq_len(n_iter)) {
    # 1) Reconstruct partial-block var_curr
    var_curr <- matrix(0, p, p)

    # Unblocked diagonal => exp(theta_var1[i,])
    if (any(!has_cov)) {
      var_curr[!has_cov, !has_cov] <- diag(exp(theta_var1[i, ]), sum(!has_cov))
    }

    # Jacobian for var2
    jac_var2 <- 0
    prior_var2 <- 0
    # Blocked => unwind_chol
    if (any(has_cov)) {
      # Fill in theta_var2
      for(block in block_groups){
        cur_idx <- par_group == block
        var_block <- unwind_chol(theta_var2_list[[block]][i, ], reverse=TRUE)
        var_curr[cur_idx, cur_idx] <- var_block
        # Prior influence
        d_block   <- sum(cur_idx)
        prior_var2 <- prior_var2 + log( robust_diwish(
          W = var_block,
          v = prior$v + d_block - 1,
          S = 2 * prior$v * diag( 1 / exp(theta_a[i,cur_idx]) )
        ))
        # Jacobian influence
        jac_var2 <- jac_var2 + calc_log_jac_chol(theta_var2_list[[block]][i, ])
      }
    }
    # 3) Compute group likelihood = sum_{s=1..n_subj} dmvnorm(alpha_s, mu_s, var_curr)
    regressors_i <- numeric(M)
    regressors_i[mean_index]  <- theta_mu[i, ]   # fill the intercept positions
    if(B > 0){
      regressors_i[!mean_index] <- theta_beta[i, ] # fill the slope positions
    }
    group_ll <- 0
    for (s in seq_len(n_subj)) {
      alpha_s <- proposals_list[[s]][i,]  # length p
      mu_s <- numeric(p)            # will hold subject s's mean for each row k=1..p
      par_idx <- 0
      for (k in seq_len(p)) {
        # The row-k design vector for subject s is group_designs[[k]][s, ] => e.g. (1 x m_k).
        x_sk <- group_designs[[k]][s, , drop = FALSE]
        # Then the mean for row k is the dot product of x_sk with the relevant slice of 'regressors_i'.
        mu_s[k] <- x_sk %*% regressors_i[par_idx + 1:ncol(group_designs[[k]])]
        par_idx <- par_idx + ncol(group_designs[[k]])
      }
      # Now alpha_s ~ N(mu_s, var_curr)
      group_ll <- group_ll + dmvnorm(alpha_s, mu_s, var_curr, log = TRUE)
    }

    # 4) Prior on var1, var2, and a => same partial-block logic
    # "var1" => Inverse-Gamma with shape=v/2, rate=v/exp(theta_a[i, !has_cov])
    # "var2" => robust_diwish( var_block, v= v + sum(has_cov)-1, S= 2 v diag(1/ a_half ) ), but a_half=exp(theta_a)
    # "a" => a_half = exp(theta_a), shape=1/2, rate=1/(A^2)
    prior_var1 <- 0
    if (any(!has_cov)) {
      # sum of logdinvGamma for each unblocked row
      # x = exp(theta_var1), shape= v/2, rate= v / exp(theta_a)
      # note that 'theta_a[i, !has_cov]' picks the relevant subset
      prior_var1 <- sum( logdinvGamma(
        x     = exp(theta_var1[i, ]),
        shape = prior$v/2,
        rate  = prior$v / exp( theta_a[i, !has_cov, drop=FALSE] )
      ))
    }

    prior_a <- sum(logdinvGamma(
      x     = exp(theta_a[i, ]),
      shape = 1/2,
      rate  = 1/(prior$A^2)
    ))
    # Combine
    sum_out[i] <- group_ll + prior_var1 + prior_var2 + jac_var2 + prior_a
  } # end loop i

  #--- 5) Add the vectorized prior on theta_mu, beta, plus Jacobians var1, a
  # prior_mu_log, prior_beta_log each length=n_iter
  sum_out <- sum_out + prior_mu_log + prior_beta_log + jac_var1 + jac_a_vec

  return(sum_out)
}


bridge_add_info_standard <- function(info, samples){
  if(is.null(samples$is_blocked)){
    has_cov <- rep(T, samples$n_pars)
  } else{
    has_cov <- samples$is_blocked
  }
  info$par_group <- samples$par_group
  if(is.null(info$par_group)){
    info$par_group <- rep(1, samples$n_pars)
  }
  info$is_blocked <- has_cov
  info$betas <- samples$betas
  info$group_designs   <- samples$group_designs
  info$group_idx <- (samples$n_pars*samples$n_subjects + 1):(samples$n_pars*samples$n_subjects + 2*samples$n_pars + length(samples$betas) +
                                                               sum(!has_cov) + (sum(has_cov) * (sum(has_cov) +1))/2)
  return(info)
}

bridge_add_group_standard <- function(all_samples, samples, idx){
  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[,idx]))
  all_samples <- cbind(all_samples, t(log(samples$samples$a_half[,idx])))
  if(!is.null(samples$betas)){
    all_samples <- cbind(all_samples, t(samples$samples$beta[,idx]))
  }
  par_group <- samples$par_group
  has_cov <- samples$is_blocked
  if(is.null(has_cov)){
    has_cov <- rep(T, samples$n_pars)
  }
  par_group <- samples$par_group
  if(is.null(par_group)){
    par_group <- rep(1, samples$n_pars)
  }

  if(any(!has_cov)){
    all_samples <- cbind(all_samples, t(log(matrix(apply(samples$samples$theta_var[!has_cov,!has_cov,idx, drop = F], 3, diag), ncol = nrow(all_samples)))))
  }
  if(any(has_cov)){
    for(block in unique(par_group[has_cov])){
      cur_idx <- par_group == block
      all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$theta_var[cur_idx,cur_idx,idx, drop = F], 3, unwind_chol), ncol = nrow(all_samples))))
    }
  }
  return(all_samples)
}


# for IC ------------------------------------------------------------------

group__IC_standard <- function(emc, stage="sample", filter=NULL) {
  # 1) Retrieve draws
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    return_mcmc = FALSE, merge_chains = TRUE)
  theta_mu <- get_pars(emc, selection = "mu", stage = stage, filter = filter,
                       return_mcmc = FALSE, merge_chains = TRUE)
  theta_var <- get_pars(emc, selection = "Sigma", stage = stage, filter = filter,
                        return_mcmc = FALSE, merge_chains = TRUE, remove_constants = F)
  if(!is.null(emc[[1]]$betas)){
    theta_beta <- get_pars(emc, selection = "beta", stage = stage, filter = filter,
                         return_mcmc = FALSE, merge_chains = TRUE)               # (B, N) or NULL
  }


  p <- dim(alpha)[1]
  n_subj <- dim(alpha)[2]
  N <- dim(alpha)[3]
  has_betas <- !is.null(emc[[1]]$betas)

  # 2) Averages
  mean_alpha <- apply(alpha, c(1,2), mean)
  mean_mu <- rowMeans(theta_mu)
  mean_var <- apply(theta_var, c(1,2), mean)
  mean_beta <- if(has_betas) rowMeans(theta_beta) else NULL

  # 3) Build bridging design if betas
  if(has_betas) {
    group_designs <- add_intercepts(emc[[1]]$par_names, emc[[1]]$group_designs, n_subj)
    # figure out which columns go to mu vs beta
    M <- p + length(mean_beta)
    mean_index <- unlist(lapply(group_designs, function(x)
      c(TRUE, rep(FALSE, ncol(x) - 1))))
  }

  # 4) Per-draw log-likelihood
  lls <- numeric(N)
  for(i in seq_len(N)) {
    var_i <- theta_var[,, i]
    if(!has_betas) {
      # no betas => old approach
      lls[i] <- sum(dmvnorm(t(alpha[,, i]), theta_mu[, i], var_i, log=TRUE))
    } else {
      # bridging approach
      regressors_i <- numeric(p + length(mean_beta))
      regressors_i[mean_index] <- theta_mu[, i]
      regressors_i[!mean_index] <- theta_beta[, i]
      ll_i <- 0
      for(s in seq_len(n_subj)) {
        # build M_i row by row
        M_i <- numeric(p) # we'll store the mean for each row k
        par_start <- 1
        for(k in seq_len(p)) {
          x_sk <- group_designs[[k]][s,,drop=FALSE]
          ncols_k <- ncol(group_designs[[k]])
          M_i[k] <- x_sk %*% regressors_i[par_start:(par_start + ncols_k - 1)]
          par_start <- par_start + ncols_k
        }
        alpha_s <- alpha[, s, i]
        ll_i <- ll_i + dmvnorm(alpha_s, M_i, var_i, log=TRUE)
      }
      lls[i] <- ll_i
    }
  }

  minD <- -2*max(lls)
  mean_ll <- mean(lls)

  # 5) Likelihood at posterior mean
  mean_pars_ll <- 0
  if(!has_betas) {
    for(s in seq_len(n_subj)) {
      mean_pars_ll <- mean_pars_ll + dmvnorm(mean_alpha[, s], mean_mu, mean_var, log=TRUE)
    }
  } else {
    # bridging with mean mu, mean beta
    regressors_mean <- numeric(p + length(mean_beta))
    regressors_mean[mean_index] <- mean_mu
    regressors_mean[!mean_index] <- mean_beta
    for(s in seq_len(n_subj)) {
      M_s <- numeric(p)
      par_start <- 1
      for(k in seq_len(p)) {
        x_sk <- group_designs[[k]][s,,drop=FALSE]
        nc_k <- ncol(group_designs[[k]])
        M_s[k] <- x_sk %*% regressors_mean[par_start:(par_start+nc_k-1)]
        par_start <- par_start + nc_k
      }
      mean_pars_ll <- mean_pars_ll +
        dmvnorm(mean_alpha[, s], M_s, mean_var, log=TRUE)
    }
  }
  Dmean <- -2*mean_pars_ll

  list(mean_ll = mean_ll, Dmean = Dmean, minD = minD)
}
