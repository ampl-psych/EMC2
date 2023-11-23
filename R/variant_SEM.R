
sample_store_SEM <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, is_grouped, ...) {
  args <- list(...)
  n_factors <- ncol(args$Lambda_mat)
  xy <- as.matrix(args$xy)
  xeta <- as.matrix(args$xeta)
  n_cov_y <- ncol(xy)
  n_cov_eta <- ncol(xeta)
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names[!is_grouped], iters, stage)
  par_names <- par_names[!is_nuisance & !is_grouped]
  n_pars <- length(par_names)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    lambda = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, NULL, NULL)),
    B = array(NA_real_,dim = c(n_factors, n_factors, iters),dimnames = list(NULL, NULL, NULL)),
    epsilon_inv = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    delta_inv = array(NA_real_, dim = c(n_factors, n_factors, iters), dimnames = list(NULL, NULL, NULL)),
    eta = array(NA_real_, dim = c(n_subjects, n_factors, iters), dimnames = list(subject_ids, NULL, NULL)),
    K = array(NA_real_, dim = c(n_pars, n_cov_y, iters), dimnames = list(NULL, NULL, NULL)),
    G = array(NA_real_, dim = c(n_factors, n_cov_eta, iters), dimnames = list(NULL, NULL, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}


add_info_SEM <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  args <- list(...)
  n_factors <- ncol(args$Lambda_mat)
  Lambda_mat <- args$Lambda_mat
  B_mat <- args$B_mat
  K_mat <- args$K_mat
  G_mat <- args$G_mat

  xy <- as.matrix(args$xy)
  xeta <- as.matrix(args$xeta)
  n_cov_y <- ncol(xy)
  n_cov_eta <- ncol(xeta)

  n_pars <- sum(!(sampler$nuisance | sampler$grouped))
  if(is.null(prior)){
    prior <- list()
  }
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- rep(1, n_pars)
  }
  if(is.null(prior$lambda_var)){
    prior$lambda_var <- 1
  }
  if(is.null(prior$B_var)){
    prior$B_var <- 1
  }
  if(is.null(prior$a_p)){
    prior$a_d <- 5
  }
  if(is.null(prior$b_p)){
    prior$b_d <- 5
  }
  if(is.null(prior$a_e)){
    prior$a_e <- 5
  }
  if(is.null(prior$b_e)){
    prior$b_e <- 5
  }

  if(is.null(Lambda_mat)){
    Lambda_mat <- matrix(0, nrow = n_pars, ncol = n_factors)
  }

  if(is.null(B_mat)){
    B_mat <- matrix(0, nrow = n_factors, ncol = n_factors)
  }
  if(is.null(K_mat)){
    K_mat <- matrix(0, nrow = n_pars, ncol = n_cov_y)
  }
  if(is.null(G_mat)){
    G_mat <- matrix(0, nrow = n_factors, ncol = n_cov_eta)
  }


  attr(sampler, "Lambda_mat") <- Lambda_mat
  attr(sampler, "B_mat") <- B_mat
  attr(sampler, "K_mat") <- K_mat
  attr(sampler, "G_mat") <- G_mat

  sampler$xy <- xy
  sampler$xeta <- xeta
  sampler$n_cov_y <- n_cov_y
  sampler$n_cov_eta <- n_cov_eta
  sampler$prior <- prior
  sampler$n_factors <- n_factors
  return(sampler)
}

get_startpoints_SEM<- function(pmwgs, start_mu, start_var){
  n_pars <- sum(!(pmwgs$nuisance | pmwgs$grouped))
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$prior$theta_mu_mean, sd = sqrt(pmwgs$prior$theta_mu_var))
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(n_pars * 3,diag(n_pars))
  start_delta_inv <- diag(1, pmwgs$n_factors)
  start_epsilon_inv <- diag(1, n_pars)
  start_eta <- matrix(0, nrow = pmwgs$n_subjects, ncol = pmwgs$n_factors)

  start_lambda <- matrix(0, nrow = n_pars, ncol = pmwgs$n_factors)
  start_B <- matrix(0, nrow = pmwgs$n_factors, ncol = pmwgs$n_factors)
  start_K <- matrix(0, nrow = pmwgs$n_pars, ncol = pmwgs$n_cov_y)
  start_G <- matrix(0, nrow = pmwgs$n_factors, ncol = pmwgs$n_cov_eta)

  Lambda_mat <- attr(pmwgs, "Lambda_mat")
  B_mat <- attr(pmwgs, "B_mat")
  K_mat <- attr(pmwgs, "K_mat")
  G_mat <- attr(pmwgs, "G_mat")
  start_lambda[Lambda_mat != Inf] <- Lambda_mat[Lambda_mat != Inf]
  start_B[B_mat != Inf] <- B_mat[B_mat != Inf]
  start_K[K_mat != Inf] <- K_mat[K_mat != Inf]
  start_G[G_mat != Inf] <- G_mat[G_mat != Inf]
  sub_mu <- matrix(mu, nrow = pmwgs$n_subjects, ncol = n_pars, byrow = T)
  return(list(tmu = start_mu, tvar = start_var, lambda = start_lambda, B = start_B,
              K = start_K, G = start_G,
              epsilon_inv = start_epsilon_inv, delta_inv = start_delta_inv,
              eta = start_eta, sub_mu = sub_mu))
}

fill_samples_SEM <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$lambda[,,j] <- group_level$lambda
  samples$B[,,j] <- group_level$B
  samples$K[,,j] <- group_level$K
  samples$G[,,j] <- group_level$G
  samples$epsilon_inv[,,j] <- group_level$epsilon_inv
  samples$delta_inv[,,j] <- group_level$delta_inv
  samples$eta[,,j] <- group_level$eta
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}

gibbs_step_SEM <- function(sampler, alpha){
  # First some bookkeeping
  last <- last_sample_SEM(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  # Just some ease of reading
  y <- t(alpha)
  n_subjects <- sampler$n_subjects
  n_pars <- sampler$n_pars-sum(sampler$nuisance) - sum(sampler$grouped)
  n_factors <- sampler$n_factors
  n_cov_y <- sampler$n_cov_y
  n_cov_eta <- sampler$n_cov_eta
  xy <- sampler$xy
  xeta <- sampler$xeta

  # Update regression matrices
  isFree_Lambda <- hyper$Lambda_mat == Inf #For indexing
  isFree_B <- hyper$B_mat == Inf #For indexing
  isFree_K <- hyper$K_mat == Inf
  isFree_G <- hyper$G_mat == Inf
  # Keeps track of which factors have a latent structure and which don't (thus can be estimated with a covariance)
  is_structured <- rowSums(isFree_B) != 0

  # Get previous values
  eta <- matrix(last$eta, n_subjects, n_factors)
  delta_inv <- matrix(last$delta_inv, n_factors)
  epsilon_inv <- last$epsilon_inv
  lambda <- matrix(last$lambda, n_pars, n_factors)
  B <- matrix(last$B, n_factors, n_factors)
  K <- matrix(last$K, n_pars, n_cov_y)
  G <- matrix(last$G, n_factors, n_cov_eta)
  mu <- last$mu
  # Start of the Gibbs steps
  #Update mu
  mu_sig <- solve(n_subjects * epsilon_inv + diag(1/prior$theta_mu_var, nrow = n_pars))
  mu_mu <- mu_sig %*% (epsilon_inv %*% colSums(y - xy %*% t(K) - eta %*% t(lambda)) + diag(1/prior$theta_mu_var, nrow = n_pars)%*% prior$theta_mu_mean)
  mu <- rmvnorm(1, mu_mu, mu_sig)
  colnames(mu) <- colnames(y)
  # calculate mean-centered observations
  ytilde <- sweep(y, 2, mu)

  B0_inv <- solve(diag(n_factors) - B)
  Psi0_inv <- solve(B0_inv %*% solve(delta_inv) %*% t(B0_inv))
  eta_sig <- solve(Psi0_inv + t(lambda) %*% epsilon_inv %*% lambda)
  for(i in 1:n_subjects){
    eta[i,] <- rmvnorm(1, eta_sig%*% (t(lambda) %*% epsilon_inv %*% (ytilde[i,] - K %*% xy[i,]) + Psi0_inv %*% B0_inv %*% (G %*% xeta[i,])), eta_sig)
  }
  epsilon_inv <- diag(rgamma(n_pars,shape=prior$a_e+n_subjects/2, rate=prior$b_e + 0.5*colSums((ytilde - xy %*% t(K) - eta %*% t(lambda))^2)))

  # #Update lambda and K, these can be updated together since they are both regressions on y (technically mu as well, but we have a mean prior on that one)
  lambda_y <- cbind(K, lambda)
  for (j in 1:n_pars) {
    isFree <- c(isFree_K[j,], isFree_Lambda[j,])
    if(any(isFree)){ #Don't do this if there are no free entries in lambda
      etaS <- cbind(xy, eta)[,isFree]
      lambda_sig <- solve(epsilon_inv[j,j] * t(etaS) %*% etaS + solve(prior$lambda_var * diag(1,sum(isFree))))
      lambda_mu <- (lambda_sig * epsilon_inv[j,j]) %*% (t(etaS) %*% ytilde[,j]) # assume 0 prior on mean
      lambda_y[j,isFree] <- rmvnorm(1,lambda_mu,lambda_sig)
    }
  }
  K <- lambda_y[,1:n_cov_y]
  lambda <- lambda_y[,((n_cov_y + 1):ncol(lambda_y))]
  #Update B and G, these can also be updated together since they are both regressions on eta
  B_eta <- cbind(G, B)
  for(p in 1:n_factors){
    isFree <- c(isFree_G[p,], isFree_B[p,])
    if(any(isFree)){
      etaS <- cbind(xeta, eta)[,isFree]
      B_sig <- solve(delta_inv[p,p] * t(etaS) %*% etaS + solve(prior$B_var * diag(1,sum(isFree))))
      B_mu <- (B_sig * delta_inv[p,p]) %*% (t(etaS) %*% eta[,p])
      B_eta[p,isFree] <- rmvnorm(1,B_mu,B_sig)
    }
  }
  G <- B_eta[,1:n_cov_eta]
  B <- B_eta[,((n_cov_eta + 1):ncol(B_eta))]

  #Update delta_inv, using diagonal entries for structural factors and covariances between non-structural entries
  eta_sq <- t(eta - eta %*% t(B) - xeta %*% t(G)) %*% (eta - eta %*% t(B) - xeta %*% t(G))
  if(any(is_structured)){
    delta_inv[is_structured,is_structured] <- diag(rgamma(sum(is_structured) ,shape=prior$a_d+n_subjects/2,rate=prior$b_d+ 0.5*diag(eta_sq)[is_structured]), sum(is_structured))
  }
  if(any(!is_structured)){
    delta_inv[!is_structured, !is_structured] <- solve(riwish(n_subjects + prior$a_d, diag(prior$b_d, nrow = sum(!is_structured)) + eta_sq[!is_structured, !is_structured]))
  }
  sub_mu <- matrix(mu, n_subjects, n_pars, byrow = T) + xy %*% t(K) + eta %*% t(lambda)
  var <- lambda %*% t(lambda) + solve(epsilon_inv)
  return(list(tmu = mu, tvar = var, lambda = lambda,  eta = eta, B = B, K=K, G = G,
              epsilon_inv = epsilon_inv, delta_inv = delta_inv, alpha = t(y),
              sub_mu = sub_mu))
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
  mu <- parameters$sub_mu[s,]
  var <- solve(parameters$epsilon_inv)
  return(list(mu = mu, var = var))
}


get_conditionals_SEM <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  epsilon_inv <- log(apply(samples$epsilon_inv[idx,idx,],3 , diag))
  eta <- matrix(samples$eta[s,,], nrow = samples$n_factors)
  lambda <- apply(samples$lambda[idx,,,drop = F], 3, unwind_lambda, samples$Lambda_mat[idx,])
  theta_mu <- samples$theta_mu[idx,]
  all_samples <- rbind(samples$alpha[idx, s,],theta_mu, eta, epsilon_inv, lambda)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration],
                                 samples$eta[s,,iteration],
                                 log(diag(samples$epsilon_inv[idx,idx, iteration])),
                                 unwind_lambda(samples$lambda[idx,, iteration], samples$Lambda_mat[idx,])))
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
    Lambda_mat = attributes(sampler)$Lambda_mat
  )
}



# bridge_sampling ---------------------------------------------------------
bridge_add_info_SEM <- function(info, samples, idx){
  info$Lambda_mat <- attr(samples, "Lambda_mat")
  info$B_mat <- attr(samples, "B_mat")
  info$K_mat <- attr(samples, "K_mat")
  info$G_mat <- attr(samples, "G_mat")
  info$n_factors <- samples$n_factors
  info$n_cov_y <- samples$n_cov_y
  info$n_cov_eta <- samples$n_cov_eta
  info$xy <- samples$xy
  # How many free regressors do we have
  free_regrs <- sum(info$Lambda_mat == Inf) + sum(info$B_mat == Inf) + sum(info$K_mat == Inf) + sum(info$G_mat == Inf)
  # Also group_level mean and residual parameter variances
  other <- samples$n_pars + samples$n_pars
  # Now we split residual factor variances in structured and unstructured
  is_structured <- rowSums(info$B_mat == Inf) != 0
  # We only get one parameter for the structured and a cholesky decomp number of parameters for the unstructured
  other <- other + sum(is_structured) + (sum(!is_structured) * (sum(!is_structured) +1))/2
  info$is_structured <- is_structured
  info$group_idx <- (samples$n_pars*samples$n_subjects + 1):(samples$n_pars*samples$n_subjects + free_regrs + other)
  # add factor scores here, they're not a parameter with a prior on them and I treat them as such
  info$eta <- samples$samples$eta[,,idx, drop = F]
  return(info)
}


bridge_add_group_SEM <- function(all_samples, samples, idx){
  Lambda_mat <- attr(samples, "Lambda_mat")
  B_mat <- attr(samples, "B_mat")
  K_mat <- attr(samples, "K_mat")
  G_mat <- attr(samples, "G_mat")

  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[,idx]))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$lambda[,,idx,drop = F], 3, unwind_lambda, Lambda_mat), ncol = nrow(all_samples))))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$B[,,idx,drop = F], 3, unwind_lambda, B_mat), ncol = nrow(all_samples))))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$K[,,idx,drop = F], 3, unwind_lambda, K_mat), ncol = nrow(all_samples))))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$G[,,idx,drop = F], 3, unwind_lambda, G_mat), ncol = nrow(all_samples))))

  all_samples <- cbind(all_samples, t(log(apply(samples$samples$epsilon_inv[,,idx], 3, diag))))
  # For delta we split it in structured and unstructured parts of the covariance matrix
  # The unstructured parts get covariances and thus have to be decomposed with cholesky
  is_structured <- rowSums(B_mat == Inf) != 0

  all_samples <- cbind(all_samples, t(log(matrix(apply(samples$samples$delta_inv[is_structured,is_structured,idx, drop = F], 3, diag), ncol = nrow(all_samples)))))
  all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$delta_inv[!is_structured,!is_structured,idx, drop = F], 3, unwind_chol), ncol = nrow(all_samples))))

  return(all_samples)
}




bridge_group_and_prior_and_jac_SEM <- function(proposals_group, proposals_list, info){
  prior <- info$prior
  proposals <- do.call(cbind, proposals_list)
  # Keep an index of where we are to clean things up, start with mu, so much bookkeeping!
  prev_end <- new_end <- info$n_pars
  theta_mu <- proposals_group[,1:new_end, drop = F]
  new_end <- prev_end + sum(info$Lambda_mat == Inf)

  # Regressors
  lambda <- proposals_group[,(prev_end+1):new_end, drop = F]
  prev_end <- new_end
  new_end <- prev_end + sum(info$B_mat == Inf)
  B <- proposals_group[,(prev_end+1):new_end, drop = F]
  prev_end <- new_end
  new_end <- prev_end + sum(info$K_mat == Inf)
  K <- proposals_group[,(prev_end+1):new_end, drop = F]
  prev_end <- new_end
  new_end <- prev_end +  sum(info$G_mat == Inf)
  G <- proposals_group[,(prev_end+1):new_end, drop = F]

  # Others
  prev_end <- new_end
  new_end <- prev_end + info$n_pars
  epsilon_inv <- proposals_group[,(prev_end+1):new_end, drop = F]

  # prev_end <- new_end
  # new_end <- prev_end + info$n_factors*info$n_subjects
  # eta <- proposals_group[,(prev_end+1):new_end, drop = F]
  #
  # Get the deltas separately
  prev_end <- new_end
  new_end <- prev_end + sum(info$is_structured)
  delta_inv1 <- proposals_group[,(prev_end+1):new_end, drop = F]

  prev_end <- new_end
  new_end <- prev_end + (sum(!info$is_structured) * (sum(!info$is_structured) +1))/2
  delta_inv2 <- proposals_group[,(prev_end+1):new_end, drop = F]

  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)
  for(i in 1:n_iter){ # these unfortunately can't be vectorized
    # Put all our stuff back together
    lambda_curr <- unwind_lambda(lambda[i,], info$Lambda_mat, reverse = T)
    # B_curr <- unwind_lambda(B[i,], info$B_mat, reverse = T)
    K_curr <- unwind_lambda(K[i,], info$K_mat,  reverse = T)
    # G_curr <- unwind_lambda(G[i,], info$G_mat, reverse = T)
    eta_curr <- info$eta[,,i]

    proposals_curr <- matrix(proposals[i,], ncol = info$n_pars, byrow = T)

    group_means <- matrix(theta_mu[i,], nrow = info$n_subjects, ncol = info$n_pars, byrow = T) + info$xy %*% t(K_curr) + eta_curr %*% t(lambda_curr)
    group_ll <- sum(dnorm(t(proposals_curr), t(group_means), exp(1/epsilon_inv[i,]), log = T))

    delta2_curr <- unwind_chol(delta_inv2[i,], reverse = T)
    prior_delta1 <- sum(logdinvGamma(exp(delta_inv1[i,]), shape = prior$a_d, rate = prior$b_d))
    prior_delta2 <- log(robust_diwish(delta2_curr, v=prior$a_d, S = diag(prior$b_d, sum(!info$is_structured))))
    prior_epsilon_inv <- sum(logdinvGamma(exp(epsilon_inv[i,]), shape = prior$a_e, rate = prior$b_e))
    jac_delta2 <- log(2^sum(!info$is_structured))+sum((sum(!info$is_structured) + 1)*log(diag(delta2_curr))) # Log of derivative of cholesky transformation
    sum_out[i] <- group_ll + prior_epsilon_inv + prior_delta1 + prior_delta2 + jac_delta2
  }
  prior_mu <- dmvnorm(theta_mu, mean = prior$theta_mu_mean, sigma = diag(prior$theta_mu_var), log =T)
  prior_lambda <- dmvnorm(lambda, mean = rep(0, ncol(lambda)), sigma = diag(prior$lambda_var, ncol(lambda)), log = T)
  prior_B <- dmvnorm(B, mean = rep(0, ncol(B)), sigma = diag(prior$B_var, ncol(B)), log = T)
  prior_K <- dmvnorm(K, mean = rep(0, ncol(K)), sigma = diag(prior$lambda_var, ncol(K)), log = T)
  prior_G <- dmvnorm(G, mean = rep(0, ncol(G)), sigma = diag(prior$B_var, ncol(G)), log = T)

  jac_delta1 <- rowSums(delta_inv1)
  jac_epsilon_inv <- rowSums(epsilon_inv)
  return(sum_out + prior_mu + prior_lambda + prior_B + prior_K + prior_G + jac_delta1 + jac_epsilon_inv) # Output is of length nrow(proposals)
}
