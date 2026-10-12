sample_store_SEM <- function(data, par_names, iters = 1, stage = "init", integrate = T, ...) {
  args <- list(...)
  sem_settings <- args$sem_settings
  if (is.null(sem_settings)) stop("sem_settings is required")

  n_factors <- if (!is.null(sem_settings$Lambda_mat)) ncol(sem_settings$Lambda_mat) else 0
  covariates <- sem_settings$covariates
  if (is.null(covariates)) {
    n_cov <- 0
  } else {
    covariates <- as.matrix(covariates)
    n_cov <- ncol(covariates)
  }

  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  n_pars <- length(par_names)
  x_names <- if (n_cov > 0) colnames(covariates) else character(0)
  factor_names <- if (n_factors > 0) paste0("F", seq_len(n_factors)) else character(0)
  if(!is.null(colnames(sem_settings$Lambda_mat)))  factor_names <- colnames(sem_settings$Lambda_mat)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    lambda = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, factor_names, NULL)),
    B = array(NA_real_,dim = c(n_factors, n_factors, iters),dimnames = list(factor_names, factor_names, NULL)),
    epsilon_inv = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL)),
    delta_inv = array(NA_real_, dim = c(n_factors, n_factors, iters), dimnames = list(factor_names, factor_names, NULL)),
    eta = array(NA_real_, dim = c(n_subjects, n_factors, iters), dimnames = list(subject_ids, factor_names, NULL)),
    K = array(NA_real_, dim = c(n_pars, n_cov, iters), dimnames = list(par_names, x_names, NULL)),
    G = array(NA_real_, dim = c(n_factors, n_cov, iters), dimnames = list(factor_names, x_names, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}


# Lambda_mat and K_mat are applied to the parameter vector by position, so their
# rows must be in the same order as the sampled parameters. A joint design list
# supplied to make_sem_structure in a different order than the one being fitted
# yields a row permutation that would otherwise be applied silently.
align_sem_rows <- function(mat, par_names, mat_name){
  if (is.null(mat)) return(NULL)
  if (nrow(mat) != length(par_names)) {
    stop(mat_name, " has ", nrow(mat), " rows but the model has ",
         length(par_names), " sampled parameters.")
  }
  row_names <- rownames(mat)
  if (is.null(row_names)) {
    rownames(mat) <- par_names
    return(mat)
  }
  if (identical(row_names, par_names)) return(mat)
  if (!setequal(row_names, par_names)) {
    stop(mat_name, " rownames do not match the sampled parameters.\n",
         "  Missing from ", mat_name, ": ",
         paste(setdiff(par_names, row_names), collapse = ", "), "\n",
         "  Not a parameter: ",
         paste(setdiff(row_names, par_names), collapse = ", "), "\n",
         "Was make_sem_structure given the same design as the one being fitted?")
  }
  mat[par_names, , drop = FALSE]
}

add_info_SEM <- function(sampler, prior = NULL, ...){
  args <- list(...)
  sem_settings <- args$sem_settings
  if (is.null(sem_settings)) stop("sem_settings is required")

  n_pars <- sampler$n_pars
  par_names <- sampler$par_names
  sem_settings$Lambda_mat <- align_sem_rows(sem_settings$Lambda_mat, par_names, "Lambda_mat")
  sem_settings$K_mat <- align_sem_rows(sem_settings$K_mat, par_names, "K_mat")

  # Create backup matrices if NULL
  if (is.null(sem_settings$Lambda_mat)) {
    n_factors <- if (!is.null(sem_settings$B_mat)) ncol(sem_settings$B_mat) else 0
    sem_settings$Lambda_mat <- matrix(0, nrow = n_pars, ncol = n_factors)
  } else {
    n_factors <- ncol(sem_settings$Lambda_mat)
  }

  if (is.null(sem_settings$B_mat) && n_factors > 0) {
    sem_settings$B_mat <- matrix(0, nrow = n_factors, ncol = n_factors)
  }

  n_cov <- if (!is.null(sem_settings$covariates)) ncol(sem_settings$covariates) else 0

  if (is.null(sem_settings$K_mat)) {
    sem_settings$K_mat <- matrix(0, nrow = n_pars, ncol = n_cov)
  }

  if (is.null(sem_settings$G_mat) && n_factors > 0) {
    sem_settings$G_mat <- matrix(0, nrow = n_factors, ncol = n_cov)
  }

  if (is.null(sem_settings$factor_groups)) {
    sem_settings$factor_groups <- if (n_factors > 0) 1:n_factors else integer(0)
  }

  sampler$sem_settings <- sem_settings
  sampler$prior <- get_prior_SEM(prior, n_pars = n_pars, sample = F, sem_settings = sem_settings)
  sampler$n_factors <- n_factors
  normalise_SEM_prior(sampler)
}

normalise_factor_var_prior <- function(x, n_factors, factor_names) {
  if (length(x) == 1L) return(rep(unname(x), n_factors))
  if (!is.null(names(x))) x <- x[factor_names]
  unname(x)
}

normalise_lambda_matrix_prior <- function(x, n_pars, n_factors,
                                          par_names = NULL, factor_names = NULL) {
  if (is.matrix(x)) {
    if (!is.null(rownames(x)) && !is.null(par_names)) x <- x[par_names, , drop = FALSE]
    if (!is.null(colnames(x)) && !is.null(factor_names)) x <- x[, factor_names, drop = FALSE]
    return(unname(x))
  }
  if (length(x) > 1L && !is.null(names(x))) x <- x[factor_names]
  matrix(x, n_pars, n_factors, byrow = TRUE)
}

normalise_SEM_prior <- function(sampler) {
  settings <- sampler$sem_settings
  settings$covariates <- if (is.null(settings$covariates)) {
    matrix(numeric(0), sampler$n_subjects, 0)
  } else as.matrix(settings$covariates)
  factor_names <- colnames(settings$Lambda_mat)
  if (is.null(factor_names)) factor_names <- if (sampler$n_factors > 0) paste0("F", seq_len(sampler$n_factors)) else character(0)
  for (name in c("lambda_var", "lambda_mean")) {
    sampler$prior[[name]] <- normalise_lambda_matrix_prior(
      sampler$prior[[name]], sampler$n_pars, sampler$n_factors,
      sampler$par_names, factor_names)
  }
  sampler$sem_settings <- settings
  sampler
}

sem_moments <- function(mu, lambda, B, K, G, epsilon_inv, delta_inv,
                        covariates, population = TRUE) {
  mu <- as.vector(mu)
  p <- length(mu)
  q <- sqrt(length(B))
  x <- as.matrix(covariates)
  effect <- matrix(K, p, ncol(x))
  variance <- diag(1 / as.vector(epsilon_inv), p)
  if (q > 0) {
    loading <- matrix(lambda, p, q) %*% solve(diag(q) - matrix(B, q, q))
    effect <- effect + loading %*% matrix(G, q, ncol(x))
    variance <- variance + loading %*% solve(matrix(delta_inv, q, q)) %*% t(loading)
  }
  if (population) {
    list(mu = drop(mu + effect %*% colMeans(x)),
         var = variance + effect %*% cov(x) %*% t(effect))
  } else {
    list(mu = sweep(effect %*% t(x), 1, mu, "+"), var = variance)
  }
}

get_prior_SEM <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL, sem_settings = NULL){

  if (is.null(sem_settings)) stop("sem_settings is required")

  Lambda_mat <- sem_settings$Lambda_mat
  B_mat <- sem_settings$B_mat
  K_mat <- sem_settings$K_mat
  G_mat <- sem_settings$G_mat
  factor_groups_prior <- sem_settings$factor_groups
  covariates <- sem_settings$covariates

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
    prior$lambda_var <- rep(.7, n_factors)
  }
  if(is.null(prior$lambda_mean)){
    prior$lambda_mean <- 0
  }
  if(is.null(prior$K_var)){
    prior$K_var <- rep(1, n_cov)
  }
  if(is.null(prior$B_var)){
    prior$B_var <- rep(.5, n_factors)
  }
  if(is.null(prior$G_var)){
    prior$G_var <- rep(.5, n_cov)
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
  fn_norm <- if (!is.null(colnames(Lambda_mat))) colnames(Lambda_mat) else if (n_factors > 0) paste0("F", seq_len(n_factors)) else character(0)
  prior$a_d <- normalise_factor_var_prior(prior$a_d, n_factors, fn_norm)
  prior$b_d <- normalise_factor_var_prior(prior$b_d, n_factors, fn_norm)
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

    if(is.null(colnames(Lambda_mat))){
      factor_names <- if (n_factors > 0) paste0("F", seq_len(n_factors)) else character(0)
    } else{
      factor_names <- colnames(Lambda_mat)
    }
    samples <- list()
    par_names <- names(sampled_pars(design, doMap = F))
    if(selection %in% c("mu", "alpha", "mu_implied", "std_loadings")){
      mu <- t(mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
                               sigma = diag(prior$theta_mu_var)))
      rownames(mu) <- par_names
      if(selection %in% c("mu")){
        samples$theta_mu <- mu
      }
    }
    if(selection %in% c("regressors","std_loadings", "alpha", "mu_implied", "Sigma", "correlation", "covariance", "sigma2")){
      K <- array(0, dim = c(n_pars, n_cov, N))
      for(i in seq_len(n_cov)){
        K[,i,] <- t(mvtnorm::rmvnorm(N, sigma = diag(prior$K_var[i], n_pars)))
      }
      K <- constrain_lambda(K, K_mat)
      rownames(K) <- par_names
      colnames(K) <- colnames(covariates)
      if(selection %in% 'regressors'){
        samples$K <- K
      }
    }
    if(selection %in% c("factor_regressors","std_loadings", "alpha", "mu_implied", "Sigma", "correlation", "covariance", "sigma2")){
      G <- array(0, dim = c(n_factors, n_cov, N))
      if (n_factors > 0) for(i in seq_len(n_cov)){
        G[,i,] <- t(mvtnorm::rmvnorm(N, sigma = diag(prior$G_var[i], n_factors)))
      }
      G <- constrain_lambda(G, G_mat)

      rownames(G) <- factor_names
      colnames(G) <- colnames(covariates)
      if(selection %in% 'factor_regressors'){
        samples$G <- G
      }
    }
    if(selection %in% c("structural_regressors", "std_loadings", "alpha", "mu_implied", "Sigma", "correlation", "covariance", "sigma2")){
      B <- array(0, dim = c(n_factors, n_factors, N))
      for(i in seq_len(n_factors)){
        B[,i,] <- t(mvtnorm::rmvnorm(N, sigma = diag(prior$B_var[i], n_factors)))
      }
      B <- constrain_lambda(B, B_mat)

      rownames(B) <- colnames(B) <- factor_names
      if(selection %in% 'structural_regressors'){
        samples$B <- B
      }
    }
    if(selection %in% c("loadings", "std_loadings", "alpha", "mu_implied", "Sigma", "correlation", "covariance", "sigma2")){
      Lambda_var_mat  <- normalise_lambda_matrix_prior(prior$lambda_var, n_pars, n_factors, par_names, factor_names)
      Lambda_mean_mat <- normalise_lambda_matrix_prior(prior$lambda_mean, n_pars, n_factors, par_names, factor_names)
      lambda <- array(0, dim = c(n_pars, n_factors, N))
      for(i in seq_len(n_factors)){
        lambda[,i,] <- t(mvtnorm::rmvnorm(N, mean = Lambda_mean_mat[, i], sigma = diag(Lambda_var_mat[, i], nrow = n_pars)))
      }
      lambda <- constrain_lambda(lambda, Lambda_mat)
      rownames(lambda) <- par_names
      colnames(lambda) <- factor_names
      if(selection %in% c("loadings", "std_loadings")){
        samples$lambda <- lambda
      }
    }
    if(selection %in% c("residuals","std_loadings", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      epsilon_inv <- t(matrix(rgamma(n_pars*N, shape = prior$a_e, rate = prior$b_e),
                            ncol = n_pars, byrow = T))
      rownames(epsilon_inv) <- par_names
      if(selection %in% c("residuals")){
        samples$epsilon_inv <- epsilon_inv
      }
    }
    if(selection %in% c("factor_residuals", "std_loadings", "alpha", "correlation", "Sigma", "covariance", "sigma2")) {
      delta_inv <- array(0, dim = c(n_factors, n_factors, N))
      for(i in 1:N){
        for (group_id in unique_factor_groups_prior) {
          group_idx <- which(factor_groups_prior == group_id)
          d_block <- length(group_idx)
          if (d_block > 1) {
            S_iw_prior <- diag(prior$b_d[group_idx], d_block)
            df_iw_prior <- prior$a_d[group_idx][1]
            sampled_cov_block_prior <- riwish(df_iw_prior, S_iw_prior)
            delta_inv[group_idx, group_idx, i] <- solve(sampled_cov_block_prior)
          } else if (d_block == 1) {
            delta_inv[group_idx, group_idx, i] <- rgamma(1, shape = prior$a_d[group_idx], rate = prior$b_d[group_idx])
          }
        }
      }
      rownames(delta_inv) <- colnames(delta_inv) <- factor_names
      if(selection %in% c("factor_residuals", "std_loadings")){
        samples$delta_inv <- delta_inv
      }
    }
    if (selection == "mu_implied") {
      samples$mu_implied <- mu
      for (i in seq_len(N)) {
        effect <- matrix(K[,,i], n_pars, n_cov)
        if (n_factors > 0) {
          effect <- effect + matrix(lambda[,,i], n_pars, n_factors) %*%
            solve(diag(n_factors) - matrix(B[,,i], n_factors, n_factors)) %*%
            matrix(G[,,i], n_factors, n_cov)
        }
        samples$mu_implied[,i] <- mu[,i] + effect %*% x_mu
      }
    }
    if(selection %in% c("alpha", "sigma2", "covariance", "correlation", "Sigma", "std_loadings")) {
      vars <- array(NA_real_, c(n_pars, n_pars, N), list(par_names, par_names, NULL))
      if (selection == "alpha") {
        subjects <- if (is.null(design$Ffactors)) design[[1]]$Ffactors$subjects else design$Ffactors$subjects
        samples$alpha <- array(NA_real_, c(n_pars, length(subjects), N), list(par_names, subjects, NULL))
      }
      for(i in seq_len(N)) {
        moments <- sem_moments(if (selection %in% c("mu_implied", "alpha", "std_loadings")) mu[,i] else rep(0, n_pars),
          lambda[,,i], B[,,i], K[,,i], G[,,i], epsilon_inv[,i], delta_inv[,,i],
          covariates, population = selection != "alpha")
        vars[,,i] <- moments$var
        if (selection == "alpha") {
          samples$alpha[,,i] <- moments$mu + t(mvtnorm::rmvnorm(length(subjects), sigma = moments$var))
        }
      }
      if(selection %in% c("sigma2", "covariance", "correlation", "Sigma", "std_loadings")) samples$theta_var <- vars
    }
    out <- samples
  }
  return(out)
}

filtered_samples_SEM <- function(sampler, filter) {
  samples <- sampler$samples
  list(theta_mu = samples$theta_mu[, filter, drop = FALSE],
       lambda = samples$lambda[, , filter, drop = FALSE],
       epsilon_inv = samples$epsilon_inv[, filter, drop = FALSE],
       eta = samples$eta[, , filter, drop = FALSE],
       alpha = samples$alpha[, , filter, drop = FALSE],
       n_factors = sampler$n_factors, iteration = length(filter),
       sem_settings = sampler$sem_settings)
}

group__IC_SEM <- function(emc, stage = "sample", filter = NULL, ...) {
  p <- emc[[1]]$n_pars
  n <- emc[[1]]$n_subjects
  mean_alpha <- mean_mu <- matrix(0, p, n)
  mean_var <- matrix(0, p, p)
  lls <- list()
  for (chain in emc) {
    samples <- chain$samples
    idx <- drop(filter_emc(matrix(get_idx(list(chain), stage), 1), filter = filter))
    for (i in idx) {
      moments <- sem_moments(samples$theta_mu[,i], samples$lambda[,,i],
        samples$B[,,i], samples$K[,,i], samples$G[,,i], samples$epsilon_inv[,i],
        samples$delta_inv[,,i], chain$sem_settings$covariates, population = FALSE)
      alpha <- matrix(samples$alpha[,,i], p, n)
      lls[[length(lls) + 1L]] <- mvtnorm::dmvnorm(t(alpha - moments$mu), sigma = moments$var, log = TRUE)
      mean_alpha <- mean_alpha + alpha
      mean_mu <- mean_mu + moments$mu
      mean_var <- mean_var + moments$var
    }
  }
  N <- length(lls)
  lls <- do.call(cbind, lls)
  if (isTRUE(list(...)$for_WAIC)) return(lls)
  list(mean_ll = mean(colSums(lls)), minD = -2 * max(colSums(lls)),
       Dmean = -2 * sum(mvtnorm::dmvnorm(t((mean_alpha - mean_mu) / N),
                                      sigma = mean_var / N, log = TRUE)))
}

get_mu_implied <- function(x, idx) {
  samples <- x$samples
  mu <- samples$theta_mu[,idx, drop = FALSE]
  for (i in seq_along(idx)) {
    j <- idx[i]
    mu[,i] <- sem_moments(samples$theta_mu[,j], samples$lambda[,,j],
      samples$B[,,j], samples$K[,,j], samples$G[,,j], samples$epsilon_inv[,j],
      samples$delta_inv[,,j], x$sem_settings$covariates)$mu
  }
  mu
}


bridge_add_info_SEM <- function(info, samples){
  sem_settings <- samples$sem_settings
  if (is.null(sem_settings)) stop("sem_settings not found in samples object for bridge sampling.")

  info$sem_settings <- sem_settings
  info$n_factors <- samples$n_factors

  Lambda_mat <- sem_settings$Lambda_mat
  B_mat <- sem_settings$B_mat
  K_mat <- sem_settings$K_mat
  G_mat <- sem_settings$G_mat
  factor_groups_bs <- sem_settings$factor_groups

  free_regrs <- sum(Lambda_mat == Inf) + sum(B_mat == Inf) + sum(K_mat == Inf) + sum(G_mat == Inf)
  other <- samples$n_pars + samples$n_pars

  if (is.null(factor_groups_bs)) {
      factor_groups_bs <- 1:info$n_factors
  }
  unique_f_groups_bs <- unique(factor_groups_bs)
  n_factor_params_delta <- 0
  for (fg_id in unique_f_groups_bs) {
    fg_idx <- which(factor_groups_bs == fg_id)
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

  all_samples <- cbind(all_samples, t(log(samples$samples$epsilon_inv[,idx])))

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




bridge_group_and_prior_and_jac_SEM <- function(proposals_group, proposals_list, info) {
  prior <- info$prior
  p <- info$n_pars
  q <- info$n_factors
  settings <- info$sem_settings
  factor_names <- colnames(settings$Lambda_mat)
  if (is.null(factor_names)) factor_names <- if (q > 0) paste0("F", seq_len(q)) else character(0)
  a_d <- normalise_factor_var_prior(prior$a_d, q, factor_names)
  b_d <- normalise_factor_var_prior(prior$b_d, q, factor_names)
  groups <- split(seq_len(q), factor(settings$factor_groups, levels = unique(settings$factor_groups)))
  matrices <- settings[c("Lambda_mat", "B_mat", "K_mat", "G_mat")]
  lambda_mean <- normalise_lambda_matrix_prior(prior$lambda_mean, p, q,
    rownames(settings$Lambda_mat), factor_names)
  lambda_var <- normalise_lambda_matrix_prior(prior$lambda_var, p, q,
    rownames(settings$Lambda_mat), factor_names)
  variances <- list(lambda_var,
    matrix(prior$B_var, q, q, byrow = TRUE),
    matrix(prior$K_var, p, ncol(settings$K_mat), byrow = TRUE),
    matrix(prior$G_var, q, ncol(settings$G_mat), byrow = TRUE))
  means <- list(lambda_mean, matrix(0, q, q),
    matrix(0, p, ncol(settings$K_mat)), matrix(0, q, ncol(settings$G_mat)))
  alpha <- do.call(cbind, proposals_list)
  out <- numeric(nrow(proposals_group))
  for (i in seq_along(out)) {
    pos <- p
    mu <- proposals_group[i,seq_len(p)]
    lp <- sum(dnorm(mu, prior$theta_mu_mean, sqrt(prior$theta_mu_var), log = TRUE))
    values <- matrices
    for (j in seq_along(matrices)) {
      free <- matrices[[j]] == Inf
      k <- sum(free)
      if (k > 0) {
        value <- proposals_group[i,pos + seq_len(k)]
        values[[j]][free] <- value
        lp <- lp + sum(dnorm(value, means[[j]][free], sqrt(variances[[j]][free]), log = TRUE))
      }
      pos <- pos + k
    }
    log_eps <- proposals_group[i,pos + seq_len(p)]
    eps <- exp(log_eps)
    pos <- pos + p
    lp <- lp + sum(dgamma(eps, prior$a_e, rate = prior$b_e, log = TRUE) + log_eps)
    delta_inv <- matrix(0, q, q)
    for (idx in groups) {
      d <- length(idx)
      k <- d * (d + 1) / 2
      value <- proposals_group[i,pos + seq_len(k)]
      pos <- pos + k
      if (d == 1) {
        delta_inv[idx,idx] <- exp(value)
        lp <- lp + dgamma(exp(value), a_d[idx], rate = b_d[idx], log = TRUE) + value
      } else {
        precision <- unwind_chol(value, reverse = TRUE)
        delta_inv[idx,idx] <- precision
        lp <- lp + log(robust_diwish(solve(precision), a_d[idx[1]], diag(b_d[idx], d))) +
          calc_log_jac_chol(value) - (d + 1) * as.numeric(determinant(precision, logarithm = TRUE)$modulus)
      }
    }
    moments <- sem_moments(mu, values[[1]], values[[2]], values[[3]], values[[4]],
      eps, delta_inv, settings$covariates, population = FALSE)
    residual <- matrix(alpha[i,], ncol = p, byrow = TRUE) - t(moments$mu)
    out[i] <- lp + sum(mvtnorm::dmvnorm(residual, sigma = moments$var, log = TRUE))
  }
  out
}

#' Define Structural Equation Model (SEM) Matrices
#'
#' @description
#' This function helps create the specification matrices (Lambda, B, K, G) for an SEM.
#' It takes a design object, data, factor names, covariate column names, and list-based
#' specifications for the paths to be estimated.
#' The subject-level parameter names for Lambda_mat and K_mat rows are derived from `sampled_pars(design)`.
#' It validates that covariates are consistent per subject (subject column in `data` must be named "subjects")
#' and includes an aggregated subject-level covariate data frame named `covariates` in the output list.
#' For identifiability, when a factor's loadings are given as a character vector, the first parameter
#' listed in `lambda_specs` for that factor is fixed to 1.
#'
#' @param data A data frame containing a column named "subjects"
#'   and any covariate columns specified in `covariate_cols`.
#' @param design An emc.design object, as created by the `design()` function.
#'   The parameter names for the SEM are derived from `names(sampled_pars(design))`.
#' @param covariate_cols Character vector or NULL. Column names in `data` to be used
#'   as covariates for K_mat and G_mat. If NULL, no covariates are processed.
#' @param lambda_specs A list defining factor loadings. The list names should be factor names.
#'   Each element may be given in one of two forms:
#'   \enumerate{
#'     \item A character vector of parameter names (from `names(sampled_pars(design))`) that load
#'       onto that factor. The first parameter listed is fixed to 1 for identifiability and the
#'       remaining ones are estimated.
#'       Example: `list(Factor1 = c("v_Sleft", "a_Eneutral"), Factor2 = "t0")`.
#'       Here, `Lambda_mat["v_Sleft", "Factor1"]` would be 1.
#'     \item A named numeric vector, where the names are parameter names and the values are written
#'       directly into `Lambda_mat`. A value of `Inf` marks an estimated loading; any finite value
#'       fixes the loading to that value. No automatic identifiability constraint is applied in this
#'       mode, so at least one loading should typically be fixed.
#'       Example: `list(Factor1 = c(v_Sleft = 1, a_Eneutral = Inf))`.
#'   }
#' @param b_specs A list defining regressions among factors.
#'   List names are outcome factors, elements are character vectors of predictor factors.
#'   Example: `list(Factor2 = "Factor1", Factor3 = c("Factor1", "Factor2"))`
#' @param k_specs A list defining covariate effects on subject-level parameters.
#'   List names are parameter names (from `names(sampled_pars(design))`), elements are character vectors of covariate names
#'   (must be present in `covariate_cols` and thus in the processed `covariates` data frame).
#'   Example: `list(v_Sleft = "cov1", a_Eneutral = c("cov1", "cov2"))`
#' @param g_specs A list defining covariate effects on factors.
#'   List names are factor names, elements are character vectors of covariate names.
#'   Example: `list(Factor1 = "cov1", Factor2 = c("cov1", "cov2"))`
#' @param factor_groups Optional factor-correlation block specification.
#'   Can be either (1) a character vector of factor names defining one correlated block, or
#'   (2) a list of character vectors where each element is one correlated block.
#'   Factors not listed are treated as independent singleton blocks.
#'   Example: `c("Speed", "Caution", "Ability")` or
#'   `list(c("Speed", "Caution"), c("Ability", "Strategy"))`.
#' @param fixed_value Numeric. The value used for fixed paths in the matrices that
#'   are not set to 1 for identifiability or `Inf` for estimation. Default is 0.
#'
#' @return A list containing:
#'   - `Lambda_mat`: The factor loading matrix.
#'   - `B_mat`: The matrix of regressions among factors.
#'   - `K_mat`: The matrix of covariate effects on subject-level parameters.
#'   - `G_mat`: The matrix of covariate effects on factors.
#'   - `par_names`: The subject-level parameter names derived from `sampled_pars(design)`.
#'   - `factor_names`: The provided SEM factor names.
#'   - `factor_groups`: Integer vector indicating factor-correlation blocks.
#'   - `covariates`: A data frame with one row per unique subject and columns for each covariate,
#'     containing the unique subject-level values. Column names are the covariate names.
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
#' )
#'
#' # SEM Factor names
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
#' sem_settings_definition <- make_sem_structure(
#'   data = forstmann_mod,
#'   design = example_design_obj,
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
make_sem_structure <- function(data = NULL,
                                 design,
                                 covariate_cols = NULL,
                                 lambda_specs = NULL,
                                 b_specs = NULL,
                                 k_specs = NULL,
                                 g_specs = NULL,
                                 factor_groups = NULL,
                                 fixed_value = 0) {

  subjects_col_internal <- "subjects"
  free_value_internal <- Inf
  par_names <- names(sampled_pars(design))

  if (!is.null(covariate_cols) && !is.character(covariate_cols)) stop("'covariate_cols' must be a character vector or NULL.")
  if (!is.null(lambda_specs) && !is.list(lambda_specs)) stop("'lambda_specs' must be a list.")
  factor_names <- if (is.null(lambda_specs)) character(0) else names(lambda_specs)
  if (!is.null(lambda_specs) && (is.null(factor_names) || any(factor_names == ""))) {
    stop("'lambda_specs' must be a named list of factor definitions.")
  }
  if (anyDuplicated(factor_names)) stop("Factor names in 'lambda_specs' must be unique.")
  n_pars <- length(par_names)
  n_factors <- length(factor_names)

  factor_groups_out <- if (n_factors > 0) seq_len(n_factors) else integer(0)
  names(factor_groups_out) <- factor_names
  if (!is.null(factor_groups)) {
    if (n_factors == 0) {
      stop("'factor_groups' supplied, but no factors were defined in 'lambda_specs'.")
    }
    groups_spec <- if (is.character(factor_groups)) {
      list(factor_groups)
    } else if (is.list(factor_groups)) {
      factor_groups
    } else {
      stop("'factor_groups' must be a character vector, a list of character vectors, or NULL.")
    }
    if (length(groups_spec) == 0) stop("'factor_groups' must not be empty.")

    parsed_groups <- rep(NA_integer_, n_factors)
    names(parsed_groups) <- factor_names
    seen_factors <- character(0)
    next_group_id <- 1L
    for (group_vec in groups_spec) {
      if (!is.character(group_vec) || length(group_vec) == 0) {
        stop("Each element in 'factor_groups' must be a non-empty character vector.")
      }
      if (any(is.na(group_vec)) || any(group_vec == "")) {
        stop("'factor_groups' contains missing or empty factor names.")
      }
      if (anyDuplicated(group_vec)) {
        stop("A factor was listed more than once within a 'factor_groups' block.")
      }
      unknown_factors <- setdiff(group_vec, factor_names)
      if (length(unknown_factors) > 0) {
        stop("Unknown factor name(s) in 'factor_groups': ",
             paste(unknown_factors, collapse = ", "), ".")
      }
      duplicate_across_blocks <- intersect(group_vec, seen_factors)
      if (length(duplicate_across_blocks) > 0) {
        stop("Factor name(s) repeated across 'factor_groups' blocks: ",
             paste(duplicate_across_blocks, collapse = ", "), ".")
      }
      parsed_groups[group_vec] <- next_group_id
      seen_factors <- c(seen_factors, group_vec)
      next_group_id <- next_group_id + 1L
    }
    remaining <- which(is.na(parsed_groups))
    if (length(remaining) > 0) {
      for (idx in remaining) {
        parsed_groups[idx] <- next_group_id
        next_group_id <- next_group_id + 1L
      }
    }
    factor_groups_out <- parsed_groups
  }

  processed_covariate_names <- character(0)
  processed_covariates_df <- NULL

  if (!is.null(covariate_cols) && length(covariate_cols) > 0) {
    if (!is.data.frame(data)) stop("'data' must be a data frame.")
    unique_subject_ids <- unique(data[[subjects_col_internal]])

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
    if(is.null(design$Ffactors$subjects)){
      n_sub <- length(design[[1]]$Ffactors$subjects)
    } else{
      n_sub <- length(design$Ffactors$subjects)
    }
    processed_covariates_df <- data.frame(matrix(ncol = 0, nrow = n_sub))
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
    for (f_name in names(lambda_specs)) {
      if (!f_name %in% factor_names) stop(paste0("Factor '", f_name, "' in lambda_specs not in SEM factor_names."))
      spec_for_f <- lambda_specs[[f_name]]
      if (length(spec_for_f) == 0) {
        stop(paste0("Values in lambda_specs for factor '", f_name, "' must be a non-empty character vector or a named numeric vector."))
      }

      if (is.numeric(spec_for_f)) {
        # Named numeric vector: values are written into Lambda_mat directly.
        # Inf entries are estimated, finite entries are fixed to that value.
        # No automatic identifiability constraint is applied in this mode.
        p_names_for_f <- names(spec_for_f)
        if (is.null(p_names_for_f) || any(is.na(p_names_for_f)) || any(p_names_for_f == "")) {
          stop(paste0("lambda_specs for factor '", f_name, "' is numeric but has missing or empty names; names must be parameter names from sampled_pars(design)."))
        }
        if (anyDuplicated(p_names_for_f)) {
          stop(paste0("lambda_specs for factor '", f_name, "' has duplicate parameter names."))
        }
        if (any(is.na(spec_for_f))) {
          stop(paste0("lambda_specs for factor '", f_name, "' contains NA values; use Inf for estimated loadings or a finite value for fixed loadings."))
        }
        for (p_name_spec in p_names_for_f) {
          if (!p_name_spec %in% par_names) stop(paste0("Parameter '", p_name_spec, "' for factor '", f_name, "' in lambda_specs not in names derived from sampled_pars(design)."))
        }
        Lambda_mat[p_names_for_f, f_name] <- unname(spec_for_f)
      } else if (is.character(spec_for_f)) {
        p_names_for_f <- spec_for_f

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
      } else {
        stop(paste0("Values in lambda_specs for factor '", f_name, "' must be a non-empty character vector or a named numeric vector."))
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

  if(n_cov == 0 && !is.null(k_specs)){
    stop("k_specs supplied, but no covariates_cols specified")
  }
  if(n_cov == 0 && !is.null(g_specs)){
    stop("g_specs supplied, but no covariates_cols specified")
  }

  if (n_cov > 0 && !is.null(k_specs)) {
    if (!is.list(k_specs)) stop("'k_specs' must be a list.")
    for (p_name_spec in names(k_specs)) {
      if (!p_name_spec %in% par_names) stop(paste0("Parameter '", p_name_spec, "' in k_specs not in names derived from sampled_pars(design)."))
      cov_names_for_p <- k_specs[[p_name_spec]]
      if (!is.character(cov_names_for_p)) stop(paste0("Values in k_specs for parameter '", p_name_spec, "' must be character vectors of covariate names."))
      for (cov_name in cov_names_for_p) {
        if (!cov_name %in% processed_covariate_names) stop(paste0("Covariate '", cov_name, "' for parameter '", p_name_spec,"' in k_specs not in covariate_cols"))
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
        if (!cov_name %in% processed_covariate_names) stop(paste0("Covariate '", cov_name, "' for factor '", f_name,"' in g_specs not in covariate_cols"))
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
    factor_groups = factor_groups_out,
    covariates = processed_covariates_df
  ))
}
