
sample_store_factor <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, is_grouped, ...) {
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
  base_samples <- sample_store_base(data, par_names[!is_grouped], iters, stage)
  par_names <- par_names[!is_nuisance & !is_grouped]
  n_pars <- length(par_names)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    theta_lambda = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, f_names, NULL)),
    lambda_untransf = array(NA_real_,dim = c(n_pars, n_factors, iters),dimnames = list(par_names, f_names, NULL)),
    theta_sig_err_inv = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL)),
    theta_psi_inv = array(NA_real_, dim = c(n_factors, iters), dimnames = list(NULL, NULL)),
    theta_eta = array(NA_real_, dim = c(n_subjects, n_factors, iters), dimnames = list(subject_ids, NULL, NULL))
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
  n_pars <- sum(!(sampler$nuisance | sampler$grouped))
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

  sampler$prior <- get_prior_factor(prior, sum(!(sampler$nuisance | sampler$grouped)), sample = F, n_factors = n_factors)
  sampler$n_factors <- n_factors
  return(sampler)
}


#' Prior specification and prior sampling for factor estimation
#'
#' To get the default priors for a given design: `get_prior_factor(design = design, sample = FALSE)`
#'
#' For details see Ghosh, J., & Dunson, D. B. (2009).
#' Default prior distributions and efficient posterior computation in Bayesian factor analysis.
#' *Journal of Computational and Graphical Statistics*, 18, 306-320. or
#' Stevenson, N., Innes, R. J., Gronau, Q. F., Miletic, S., Heathcote, A., PhD,
#' Forstmann, B., & Brown, S. (2024). Using group level factor models to resolve
#' high dimensionality in model-based sampling. https://doi.org/10.31234/osf.io/pn3wv.
#'
#' Note that if `sample = FALSE`, prior$theta_mu_invar (the inverse of the prior covariance matrix on the group-level mean) is returned,
#' which is only used for computational efficiency
#'
#' @param prior A named list that can contain the prior mean (`theta_mu_mean`) and
#' variance (`theta_mu_var`) on the group-level mean; the variance of the loadings (`theta_lambda_var`);
#' shape and rate of the factor variances (`ap` and `bp`) and shape and rate of the residual variances
#' (`as` and `bs`). For `NULL` entries, the default prior is is used.
#' @param n_pars Often inferred from the design, but if `design = NULL`, `n_pars`
#' will be used to determine the size of prior.
#' @param sample Whether to sample from the prior or to simply return the prior. Default is TRUE,
#' @param N How many samples to draw from the prior, the default is 1e5
#' @param design The design obtained from `design()`, required when `map = TRUE`
#' @param selection  Character. If `sample = TRUE`, what priors to sample from.
#' @param n_factors Integer. The number of factors.
#' @param Lambda_mat The loadings constraint matrix.
#'
#' @return A list with a single entry of type of samples from the prior (if `sample = TRUE`) or else a prior object
#' @examples
#' # First define a design for the model
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Now get the default prior
#' prior <- get_prior_factor(design = design_DDMaE, sample = FALSE, n_factors = 3)
#' # We can change values in the default prior or use `prior`
#' # Then we can get samples from this prior e.g.
#' samples <- get_prior_factor(prior = prior, design = design_DDMaE,
#'   sample = TRUE, selection = "mu", n_factors = 3)
#'
#' @export
get_prior_factor <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL,
                             Lambda_mat = NULL, n_factors = NULL){


  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_p_vector(design, doMap = F))
  }
  if(is.null(n_factors)) n_factors <- ncol(Lambda_mat)
  if(is.null(Lambda_mat)){
    Lambda_mat <- matrix(Inf, nrow = n_pars, ncol = n_factors)
    diag(Lambda_mat) <- 1
    Lambda_mat[upper.tri(Lambda_mat, diag = F)] <- 0
  }
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars)
  }
  if(is.null(prior$theta_mu_var)){
    prior$theta_mu_var <- rep(1, n_pars)
  }
  if(is.null(prior$theta_lambda_var)){
    prior$theta_lambda_var <- rep(.7, n_pars)
  }
  if(is.null(prior$ap)){
    prior$ap <- 2
  }
  if(is.null(prior$bp)){
    prior$bp <- .5
  }
  if(is.null(prior$as)){
    prior$as <- rep(2, n_pars)
  }
  if(is.null(prior$bs)){
    prior$bs <- rep(.1, n_pars)
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- diag(1/prior$theta_mu_var)
  prior$theta_lambda_invar <-1/prior$theta_lambda_var
  # Things I save rather than re-compute inside the loops.
  attr(prior, "type") <- "factor"
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
  n_pars <- sum(!(pmwgs$nuisance | pmwgs$grouped))
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

fill_samples_factor <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$theta_lambda[,,j] <- group_level$lambda
  samples$lambda_untransf[,,j] <- group_level$lambda_untransf
  samples$theta_sig_err_inv[,j] <- group_level$sig_err_inv
  samples$theta_psi_inv[,j] <- group_level$psi_inv
  samples$theta_eta[,,j] <- group_level$eta
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}

gibbs_step_factor <- function(sampler, alpha){
  # Gibbs step for group means with parameter expanded factor analysis from Ghosh & Dunson 2009
  # mu = theta_mu, var = theta_var
  last <- last_sample_factor(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  #extract previous values (for ease of reading)

  alpha <- t(alpha)
  n_subjects <- sampler$n_subjects
  n_pars <- sampler$n_pars-sum(sampler$nuisance) - sum(sampler$grouped)
  n_factors <- sampler$n_factors
  Lambda_mat <- hyper$Lambda_mat
  Lambda_mat <- Lambda_mat == Inf #For indexing

  eta <- matrix(last$eta, n_subjects, n_factors)
  psi_inv <- diag(last$psi_inv, n_factors)
  sig_err_inv <- diag(last$sig_err_inv)
  lambda <- matrix(last$lambda, n_pars, n_factors)
  mu <- last$mu

  #Update mu
  mu_sig <- solve(n_subjects * sig_err_inv + prior$theta_mu_invar)
  mu_mu <- mu_sig %*% (sig_err_inv %*% colSums(alpha - eta %*% t(lambda)) + prior$theta_mu_invar %*% prior$theta_mu_mean)
  mu <- rmvnorm(1, mu_mu, mu_sig)
  colnames(mu) <- colnames(alpha)
  # calculate mean-centered observations
  alphatilde <- sweep(alpha, 2, mu)

  #Update eta, I do this one first since I don't want to save eta
  eta_sig <- solve(psi_inv + t(lambda) %*% sig_err_inv %*% lambda)
  eta_mu <- eta_sig %*% t(lambda) %*% sig_err_inv %*% t(alphatilde)
  eta[,] <- t(apply(eta_mu, 2, FUN = function(x){rmvnorm(1, x, eta_sig)}))

  # for(p in 1:n_pars){
  #   constraint <- Lambda_mat[p,] == Inf
  #   for(j in 1:n_factors){
  #     alphatilde[,p] <- alphatilde[,p] - lambda[p,j] * eta[,j] * (1-constraint[j])
  #   }
  # }

  #Update sig_err
  sig_err_inv <- diag(rgamma(n_pars,shape=prior$as+n_subjects/2, rate= prior$bs + colSums((alphatilde - eta %*% t(lambda))^2)/2))

  #Update lambda
  for (j in 1:n_pars) {
    constraint <- Lambda_mat[j,] #T if item is not constraint (bit confusing tbh)
    if(any(constraint)){ #Don't do this if there are no free entries in lambda
      etaS <- eta[,constraint]
      lambda_sig <- solve(sig_err_inv[j,j] * t(etaS) %*% etaS + diag(prior$theta_lambda_invar[j], sum(constraint)))
      lambda_mu <- (lambda_sig * sig_err_inv[j,j]) %*% (t(etaS) %*% alphatilde[,j])
      lambda[j,constraint] <- rmvnorm(1,lambda_mu,lambda_sig)
    }
  }

  #Update psi_inv
  psi_inv[,] <- diag(rgamma(n_factors ,shape=prior$ap+n_subjects/2,rate=prior$bp+colSums(eta^2)/2), n_factors)
  # psi_inv <- diag(n_factors)#solve(riwish(n_subjects + prior$rho_0, t(eta) %*% eta + solve(prior$R_0)))

  lambda_orig <- lambda
  #If the diagonals of lambda aren't constrained to be 1, we should fix the signs
  if(hyper$signFix){
    for(l in 1:n_factors){
      mult <- ifelse(lambda[l, l] < 0, -1, 1) #definitely a more clever function for this
      lambda_orig[,l] <- mult * lambda[, l]
    }
  }

  var <- lambda_orig %*% solve(psi_inv) %*% t(lambda_orig) + diag(1/diag((sig_err_inv)))
  lambda_orig <- lambda_orig %*% matrix(diag(sqrt(1/diag(psi_inv)), n_factors), nrow = n_factors)
  return(list(tmu = mu, tvar = var, lambda_untransf = lambda, lambda = lambda_orig, eta = eta,
              sig_err_inv = diag(sig_err_inv), psi_inv = diag(psi_inv), alpha = t(alpha)))
}

last_sample_factor <- function(store) {
  list(
    mu = store$theta_mu[, store$idx],
    eta = store$theta_eta[,,store$idx],
    lambda = store$lambda_untransf[,,store$idx],
    psi_inv = store$theta_psi_inv[,store$idx],
    sig_err_inv = store$theta_sig_err_inv[,store$idx]
  )
}

get_conditionals_factor <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  sig_err <- log(samples$theta_sig_err_inv[idx,])
  psi <- log(samples$theta_psi_inv)
  eta <- matrix(samples$theta_eta[s,,], nrow = samples$n_factors)
  lambda <- apply(samples$lambda_untransf[idx,,,drop = F], 3, unwind_lambda, samples$Lambda_mat[idx,])
  theta_mu <- samples$theta_mu[idx,]
  all_samples <- rbind(samples$alpha[idx, s,],theta_mu, eta, sig_err, psi, lambda)#, sig_err, psi, lambda)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration],
                                 samples$theta_eta[s,,iteration],
                                 log(samples$theta_sig_err_inv[idx, iteration]),
                                 log(samples$theta_psi_inv[,iteration, drop = F]),
                                 unwind_lambda(samples$lambda_untransf[idx,, iteration], samples$Lambda_mat[idx,])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

filtered_samples_factor <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter],
    lambda_untransf = sampler$samples$lambda_untransf[, , filter, drop = F],
    theta_psi_inv = sampler$samples$theta_psi_inv[, filter, drop = F],
    theta_sig_err_inv = sampler$samples$theta_sig_err_inv[, filter],
    theta_eta = sampler$samples$theta_eta[, , filter, drop = F],
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
  all_samples <- cbind(all_samples, t(log(samples$samples$theta_sig_err_inv[,idx])))
  all_samples <- cbind(all_samples, t(log(samples$samples$theta_psi_inv[,idx, drop = F])))
  return(all_samples)
}

bridge_group_and_prior_and_jac_factor <- function(proposals_group, proposals_list, info){
  prior <- info$prior
  proposals <- do.call(cbind, proposals_list)
  theta_mu <- proposals_group[,1:info$n_pars]
  theta_lambda <- proposals_group[,(info$n_pars +1):(info$n_pars + sum(info$Lambda_mat == Inf))]
  theta_epsilon_inv <- proposals_group[,(1 + info$n_pars + sum(info$Lambda_mat == Inf)): (info$n_pars + sum(info$Lambda_mat == Inf) + info$n_pars)]
  theta_psi_inv <- proposals_group[,(1 + info$n_pars + sum(info$Lambda_mat == Inf) + info$n_pars):
                                     (info$n_pars + sum(info$Lambda_mat == Inf) + info$n_pars + info$n_factors), drop = F]

  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)
  for(i in 1:n_iter){ # these unfortunately can't be vectorized
    lambda_curr <- unwind_lambda(theta_lambda[i,], info$Lambda_mat, reverse = T)
    epsilon_curr <- diag(1/exp(theta_epsilon_inv[i,]))
    psi_curr <- diag(1/exp(theta_psi_inv[i,]), info$n_factors)
    theta_var_curr <- lambda_curr %*% psi_curr %*% t(lambda_curr) + epsilon_curr
    proposals_curr <- matrix(proposals[i,], ncol = info$n_pars, byrow = T)
    group_ll <- sum(dmvnorm(proposals_curr, theta_mu[i,], theta_var_curr, log = T))
    prior_epsilon <- sum(logdinvGamma(1/exp(theta_epsilon_inv[i,]), shape = prior$as, rate = prior$bs))
    prior_psi <- sum(logdinvGamma(1/exp(theta_psi_inv[i,]), prior$ap, rate = prior$bp))
    sum_out[i] <- group_ll + prior_epsilon + prior_psi
  }
  prior_lambda <- dmvnorm(theta_lambda, mean = rep(0, ncol(theta_lambda)),
                          sigma = diag(prior$theta_lambda_var, ncol(theta_lambda)), log = T)
  prior_mu <- dmvnorm(theta_mu, mean = prior$theta_mu_mean, sigma = diag(prior$theta_mu_var), log =T)
  jac_psi <- rowSums(theta_epsilon_inv)
  jac_epsilon <- rowSums(theta_psi_inv)
  return(sum_out + prior_mu + prior_lambda + jac_psi + jac_epsilon) # Output is of length nrow(proposals)
}



# get_all_pars_factor <- function(samples, filter, info){
#   n_subjects <- samples$n_subjects
#   n_iter = length(samples$samples$stage[samples$samples$stage== filter])
#   # Extract relevant objects
#   alpha <- samples$samples$alpha[,,samples$samples$stage== filter]
#   theta_mu <- samples$samples$theta_mu[,samples$samples$stage== filter]
#   lambda <- samples$samples$lambda_untransf[,,samples$samples$stage==filter, drop = F]
#   psi_inv <- samples$samples$theta_psi_inv[,,samples$samples$stage==filter, drop = F]
#   sig_err_inv <- samples$samples$theta_sig_err_inv[,,samples$samples$stage==filter]
#
#   Lambda_mat <- info$hyper$Lambda_mat
#   n_factors <- samples$n_factors
#
#   lambda.unwound <- apply(lambda,3,unwind_lambda, Lambda_mat)
#   sig_err_inv.diag <- log(apply(sig_err_inv, 3, diag))
#   psi_inv.diag <- matrix(log(apply(psi_inv, 3, diag)), nrow = n_factors)
#
#   # Set up
#   n_params<- nrow(alpha) + nrow(theta_mu) + nrow(sig_err_inv.diag) + nrow(psi_inv.diag) + nrow(lambda.unwound)
#   all_samples=array(dim=c(n_subjects,n_params,n_iter))
#   mu_tilde=array(dim = c(n_subjects,n_params))
#   var_tilde=array(dim = c(n_subjects,n_params,n_params))
#
#   for (j in 1:n_subjects){
#     all_samples[j,,] = rbind(alpha[,j,],theta_mu[,],sig_err_inv.diag[,],psi_inv.diag[,],lambda.unwound[,])
#     # calculate the mean for re, mu and sigma
#     mu_tilde[j,] =apply(all_samples[j,,],1,mean)
#     # calculate the covariance matrix for random effects, mu and sigma
#     var_tilde[j,,] = cov(t(all_samples[j,,]))
#   }
#
#   for(i in 1:n_subjects){ #RJI_change: this bit makes sure that the sigma tilde is pos def
#     if(!corpcor::is.positive.definite(var_tilde[i,,], tol=1e-8)){
#       var_tilde[i,,]<-corpcor::make.positive.definite(var_tilde[i,,], tol=1e-6)
#     }
#   }
#
#   X <- cbind(t(theta_mu),t(sig_err_inv.diag),t(psi_inv.diag), t(lambda.unwound))
#   info$n_params <- n_params
#   info$n_factors <- n_factors
#   info$given.ind <- (info$n_randeffect+1):n_params
#   info$X.given_ind <- 1:length(info$given.ind)
#   return(list(X = X, mu_tilde = mu_tilde, var_tilde = var_tilde, info = info))
# }
#

# group_dist_factor = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
#   n_randeffect <- info$n_randeffect
#   n_factors <- info$n_factors
#   param.theta_mu <- parameters[1:n_randeffect]
#   param.sig_err_inv <- exp(parameters[(n_randeffect+1):(n_randeffect + n_randeffect)])
#   param.psi_inv <- exp(parameters[(n_randeffect+n_randeffect+1):(n_randeffect + n_randeffect+ n_factors)])
#   param.lambda.unwound <- parameters[(n_randeffect+n_randeffect+n_factors+1):length(parameters)]
#   param.lambda <- unwind_lambda(param.lambda.unwound, info$hyper$Lambda_mat, reverse = T)
#   param.var <- param.lambda %*% diag(1/param.psi_inv, length(param.psi_inv)) %*% t(param.lambda) + diag(1/param.sig_err_inv)
#   if (sample){
#     return(mvtnorm::rmvnorm(n_samples, param.theta_mu, param.var))
#   }else{
#     logw_second<-max(-5000*info$n_randeffect, mvtnorm::dmvnorm(random_effect, param.theta_mu,param.var,log=TRUE))
#     return(logw_second)
#   }
# }
#
# prior_dist_factor = function(parameters, info){
#   n_randeffect <- info$n_randeffect
#   n_factors <- info$n_factors
#   prior <- info$prior
#   hyper <- info$hyper
#   #Extract and when necessary transform back
#   param.theta_mu <- parameters[1:n_randeffect]
#   param.sig_err_inv <- exp(parameters[(n_randeffect+1):(n_randeffect + n_randeffect)])
#   param.psi_inv <- exp(parameters[(n_randeffect+n_randeffect+1):(n_randeffect + n_randeffect+ n_factors)])
#   param.lambda.unwound <- parameters[(n_randeffect+n_randeffect+n_factors+1):length(parameters)]
#
#   log_prior_mu=sum(dnorm(param.theta_mu, mean = prior$theta_mu_mean, sd = sqrt(prior$theta_mu_var), log =TRUE))
#   log_prior_sig_err_inv = sum(pmax(-1000, dgamma(param.sig_err_inv, shape = prior$nu/2, rate = (prior$s2*prior$nu)/2, log=TRUE)))
#   log_prior_psi_inv = sum(pmax(-1000, dgamma(param.psi_inv, shape = prior$ap, rate = prior$bp, log=TRUE)))
#   log_prior_lambda=sum(dnorm(param.lambda.unwound, mean = 0, sd = sqrt(prior$theta_lambda_var), log =TRUE))
#
#   jac_sig_err_inv <- -sum(log(param.sig_err_inv)) # Jacobian determinant of transformation of log of the sig_err_inv
#   jac_psi_inv <- -sum(log(param.psi_inv)) # Jacobian determinant of transformation of log of the psi_inv
#   # Jacobians are actually part of the denominator (dnorm(prop_theta)) since transformations of the data (rather than parameters),
#   # warrant a jacobian added. But we add the jacobians here for ease of calculations.
#   return(log_prior_mu + log_prior_sig_err_inv + log_prior_psi_inv + log_prior_lambda - jac_psi_inv - jac_sig_err_inv)
# }
