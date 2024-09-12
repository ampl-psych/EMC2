
sample_store_standard <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, is_grouped,...) {
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names[!is_grouped], iters, stage)
  par_names <- par_names[!is_nuisance & !is_grouped]
  n_pars <- length(par_names)

  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    a_half = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

add_info_standard <- function(sampler, prior = NULL, ...){
  sampler$prior <- get_prior_standard(prior, sum(!(sampler$nuisance | sampler$grouped)), sample = F)
  return(sampler)
}

#' Prior specification or prior sampling for standard estimation.
#'
#' To get the default prior for a created design: get_prior_standard(design = design, sample = FALSE)
#'
#' For details see Huang, A., & Wand, M. P. (2013). Simple marginally noninformative
#' prior distributions for covariance matrices. *Bayesian Analysis*, 8, 439-452. https://doi.org/10.1214/13-BA815.
#'
#' Note that if `sample = FALSE`, prior$theta_mu_invar (the inverse of the prior covariance matrix on the group-level mean)
#' is also returned, which is only used for computational efficiency
#'
#' @param prior A named list that can contain the prior mean (`theta_mu_mean`) and
#' variance (`theta_mu_var`) on the group-level mean, or the scale (`A`), or degrees of freedom (`v`)
#' for the group-level variance-covariance matrix. For `NULL` entries, a default prior gets created.
#' @param n_pars Often inferred from the design, but if `design = NULL`, `n_pars`
#' will be used to determine the size of prior.
#' @param sample Boolean, defaults to `TRUE`, sample from the prior or simply return the prior specifications?
#' @param N How many samples to draw from the prior, the default is 1e5
#' @param design The design obtained from `design()`, required when `map = TRUE`
#' @param selection Character. If `sample = TRUE`, what prior to sample from. Options:
#' `"mu"`, `"sigma2"`, `"covariance"` `"Sigma"`, `"alpha", "correlation"`.
#'
#' @return A list with a single entry of type of samples from the prior (if sample = TRUE) or else a prior object
#' @examples
#' # First define a design for the model
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Now get the default prior
#' prior <- get_prior_standard(design = design_DDMaE, sample = FALSE)
#' # We can change values in the default prior or use `prior`
#' # Then we can get samples from this prior e.g.
#' samples <- get_prior_standard(prior = prior, design = design_DDMaE,
#'   sample = TRUE, selection = "mu")
#' @export
get_prior_standard <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL){
  # Checking and default priors
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(sampled_p_vector(design, doMap = F))
  }
  if (!is.null(prior$theta_mu_mean)) {
    n_pars <- length(prior$theta_mu_mean)
  }
  if (is.null(prior$theta_mu_mean)) {
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
    par_names <- names(sampled_p_vector(design, doMap = F))
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
  n_pars <- sum(!(pmwgs$nuisance | pmwgs$grouped))
  if (is.null(start_mu)) start_mu <- rmvnorm(1, mean = pmwgs$prior$theta_mu_mean, sigma = pmwgs$prior$theta_mu_var)
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(n_pars * 3,diag(n_pars))
  start_a_half <- 1 / rgamma(n = n_pars, shape = 2, rate = 1)
  return(list(tmu = start_mu, tvar = start_var, tvinv = ginv(start_var), a_half = start_a_half))
}

get_group_level_standard <- function(parameters, s){
  # This function is modified for other versions
  mu <- parameters$tmu
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

fill_samples_standard <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$a_half[, j] <- group_level$a_half
  samples$last_theta_var_inv <- group_level$tvinv
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}



gibbs_step_standard <- function(sampler, alpha){
  # Gibbs step for group means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  last <- last_sample_standard(sampler$samples)
  prior <- sampler$prior

  n_pars <- sampler$n_pars-sum(sampler$nuisance) - sum(sampler$grouped)
  # Here mu is group mean, so we are getting mean and variance
  var_mu <- ginv(sampler$n_subjects * last$tvinv + prior$theta_mu_invar)
  mean_mu <- as.vector(var_mu %*% (last$tvinv %*% apply(alpha, 1, sum) +
                                     prior$theta_mu_invar %*% prior$theta_mu_mean))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  tmu <- rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names[!(sampler$nuisance | sampler$grouped)]

  # New values for group var
  theta_temp <- alpha - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  B_half <- 2 * prior$v * diag(1 / last$a_half) + cov_temp # nolint
  tvar <- riwish(prior$v + n_pars - 1 + sampler$n_subjects, B_half) # New sample for group variance
  tvinv <- ginv(tvar)

  # Sample new mixing weights.
  a_half <- 1 / rgamma(n = n_pars,shape = (prior$v + n_pars) / 2,
                       rate = prior$v * diag(tvinv) + 1/(prior$A^2))
  return(list(tmu = tmu,tvar = tvar,tvinv = tvinv,a_half = a_half,alpha = alpha))
}

# conditionalSECdistr <- function (object, fixed.comp, fixed.values, name, drop = TRUE)
# {
#   family <- slot(object, "family")
#   if (!(family %in% c("SN", "ESN")))
#     stop("family must be either SN or ESN")
#   dp <- slot(object, "dp")
#   xi <- dp$xi
#   Omega <- dp$Omega
#   alpha <- dp$alpha
#   tau <- if (family == "SN")
#     0
#   else dp$tau
#   d <- length(alpha)
#   fix <- fixed.comp
#   h <- length(fix)
#   if (any(fix != round(fix)) | !all(fix %in% 1:d) | h == d)
#     stop("fixed.comp makes no sense")
#   if (length(fixed.values) != h)
#     stop("length(fixed.comp) != lenght(fixed.values)")
#   compNames <- slot(object, "compNames")
#   if (missing(name)) {
#     basename <- if (object@name != "")
#       object@name
#     else deparse(substitute(object))
#     name <- paste(basename, "|(", paste(compNames[fix], collapse = ","),
#                   ")=(", paste(format(fixed.values), collapse = ","),
#                   ")", sep = "")
#   }
#   else name <- as.character(name)[1]
#   omega <- sqrt(diag(Omega))
#   omega1 <- omega[fix]
#   omega2 <- omega[-fix]
#   R <- cov2cor(Omega)
#   R11 <- R[fix, fix, drop = FALSE]
#   R12 <- R[fix, -fix, drop = FALSE]
#   R21 <- R[-fix, fix, drop = FALSE]
#   R22 <- R[-fix, -fix, drop = FALSE]
#   alpha1 <- matrix(alpha[fix], ncol = 1)
#   alpha2 <- matrix(alpha[-fix], ncol = 1)
#   iR11 <- mnormt::pd.solve(R11)
#   R22.1 <- R22 - R21 %*% iR11 %*% R12
#   a.sum <- as.vector(t(alpha2) %*% R22.1 %*% alpha2)
#   alpha1_2 <- as.vector(alpha1 + iR11 %*% R12 %*% alpha2)/sqrt(1 +
#                                                                  a.sum)
#   tau2.1 <- (tau * sqrt(1 + sum(alpha1_2 * as.vector(iR11 %*%
#                                                        alpha1_2))) + sum(alpha1_2 * (fixed.values - xi[fix])/omega1))
#   O11 <- Omega[fix, fix, drop = FALSE]
#   O12 <- Omega[fix, -fix, drop = FALSE]
#   O21 <- Omega[-fix, fix, drop = FALSE]
#   O22 <- Omega[-fix, -fix, drop = FALSE]
#   iO11 <- (1/omega1) * iR11 * rep(1/omega1, each = h)
#   reg <- O21 %*% iO11
#   xi2.1 <- as.vector(xi[-fix] + reg %*% (fixed.values - xi[fix]))
#   O22.1 <- O22 - reg %*% O12
#   omega22.1 <- sqrt(diag(O22.1))
#   alpha2.1 <- as.vector((omega22.1/omega2) * alpha2)
#   dp2.1 <- list(xi = xi2.1, Omega = O22.1, alpha = alpha2.1,
#                 tau = tau2.1)
#   return(dp2.1)
# }

get_conditionals_standard <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  pts2_unwound <- apply(samples$theta_var[idx,idx,],3,unwind)
  all_samples <- rbind(samples$alpha[idx, s,],samples$theta_mu[idx,],pts2_unwound)
  # moments <- msn.mle(y = t(all_samples))
  # sndist <- makeSECdistr(dp=list(xi = moments$dp$beta, Omega = moments$dp$Omega, alpha = moments$dp$alpha), family="SN")
  # condsn <- conditionalSECdistr(sndist, fixed.comp = (n_pars + 1):nrow(all_samples),
  #                               fixed.values = c(samples$theta_mu[idx,iteration], unwind(samples$theta_var[idx,idx,iteration])))

  mu_tilde <- rowMeans(all_samples)
  var_tilde <- stats::cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration], unwind(samples$theta_var[idx,idx,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
  #
  # return(list(eff_mu = condsn$xi, eff_var = condsn$Omega
  #             eff_alpha = condsn$alpha, eff_tau = condsn$tau))
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

get_all_pars_standard <- function(samples, idx, info){
  n_subjects <- samples$n_subjects
  n_iter = length(samples$samples$stage[idx])
  # Exctract relevant objects
  alpha <- samples$samples$alpha[,,idx]
  theta_mu <- samples$samples$theta_mu[,idx]
  theta_var <- samples$samples$theta_var[,,idx]
  a_half <- log(samples$samples$a_half[,idx])
  theta_var.unwound = apply(theta_var,3,unwind_chol)
  # Set up
  n_params<- samples$n_pars+nrow(theta_var.unwound)+samples$n_pars
  all_samples=array(dim=c(n_subjects,n_params,n_iter))
  mu_tilde=array(dim = c(n_subjects,n_params))
  var_tilde=array(dim = c(n_subjects,n_params,n_params))

  for (j in 1:n_subjects){
    all_samples[j,,] = rbind(alpha[,j,],theta_mu[,],theta_var.unwound[,])
    # calculate the mean for re, mu and sigma
    mu_tilde[j,] =apply(all_samples[j,,],1,mean)
    # calculate the covariance matrix for random effects, mu and sigma
    var_tilde[j,,] = cov(t(all_samples[j,,]))
  }

  for(i in 1:n_subjects){ #RJI_change: this bit makes sure that the sigma tilde is pos def
    if(!corpcor::is.positive.definite(var_tilde[i,,], tol=1e-8)){
      var_tilde[i,,]<-corpcor::make.positive.definite(var_tilde[i,,], tol=1e-6)
    }
  }
  X <- cbind(t(theta_mu),t(theta_var.unwound),t(a_half))
  info$n_params <- n_params
  info$given.ind <- (info$n_randeffect+1):n_params
  info$X.given_ind <- 1:(n_params-info$n_randeffect)
  return(list(X = X, mu_tilde = mu_tilde, var_tilde = var_tilde, info = info))
}
#
# robust_diwish <- function (W, v, S) { #RJI_change: this function is to protect against weird proposals in the diwish function, where sometimes matrices weren't pos def
#   if (!is.matrix(S)) S <- matrix(S)
#   if (!is.matrix(W)) W <- matrix(W)
#   p <- nrow(S)
#   gammapart <- sum(lgamma((v + 1 - 1:p)/2))
#   ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
#   if (corpcor::is.positive.definite(W, tol=1e-8)){
#     cholW<-base::chol(W)
#   }else{
#     return(1e-10)
#   }
#   if (corpcor::is.positive.definite(S, tol=1e-8)){
#     cholS <- base::chol(S)
#   }else{
#     return(1e-10)
#   }
#   halflogdetS <- sum(log(diag(cholS)))
#   halflogdetW <- sum(log(diag(cholW)))
#   invW <- chol2inv(cholW)
#   exptrace <- sum(S * invW)
#   lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * exptrace
#   lpdf <- lnum - ldenom
#   out <- exp(lpdf)
#   if(!is.finite(out)) return(1e-100)
#   if(out < 1e-10) return(1e-100)
#   return(exp(lpdf))
# }

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


group_dist_standard <- function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  n_randeffect <- info$n_randeffect
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind_chol(param.theta.sig.unwound, reverse = TRUE)
  if (sample){
    return(rmvnorm(n_samples, param.theta.mu,param.theta.sig2))
  }else{
    logw_second<-max(-5000*info$n_randeffect, dmvnorm(random_effect, param.theta.mu,param.theta.sig2,log=TRUE))
    return(logw_second)
  }
}

prior_dist_standard <- function(parameters, info){
  n_randeffect <- info$n_randeffect
  prior <- info$prior
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind_chol(param.theta.sig.unwound, reverse = TRUE)
  param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
  log_prior_mu=dmvnorm(param.theta.mu, mean = prior$theta_mu_mean, sigma = prior$theta_mu_var, log =TRUE)
  log_prior_sigma = log(robust_diwish(param.theta.sig2, v=prior$v+ n_randeffect-1, S = 2*prior$v*diag(1/param.a)))
  log_prior_a = sum(logdinvGamma(param.a,shape = 1/2,rate=1/(prior$A^2)))
  # These are Jacobian corrections for the transformations on these
  logw_den2 <- -sum(log(param.a))
  logw_den3 <- -(log(2^n_randeffect)+sum((n_randeffect:1+1)*log(diag(param.theta.sig2))))
  return(log_prior_mu + log_prior_sigma + log_prior_a - logw_den3 - logw_den2)
}



# bridge_sampling ---------------------------------------------------------
bridge_add_group_standard <- function(all_samples, samples, idx){
  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[,idx]))
  all_samples <- cbind(all_samples, t(log(samples$samples$a_half[,idx])))
  all_samples <- cbind(all_samples, t(apply(samples$samples$theta_var[,,idx], 3, unwind_chol)))
  return(all_samples)
}

bridge_add_info_standard <- function(info, samples){
  info$group_idx <- (samples$n_pars*samples$n_subjects + 1):(samples$n_pars*samples$n_subjects + 2*samples$n_pars + (samples$n_pars * (samples$n_pars +1))/2)
  return(info)
}


bridge_group_and_prior_and_jac_standard <- function(proposals_group, proposals_list, info){
  prior <- info$prior
  proposals <- do.call(cbind, proposals_list)
  theta_mu <- proposals_group[,1:info$n_pars]
  theta_a <- proposals_group[,(info$n_pars + 1):(2*info$n_pars)]
  theta_var <- proposals_group[,(2*info$n_pars + 1):(2*info$n_pars + info$n_pars*(info$n_pars + 1)/2)]
  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)
  for(i in 1:n_iter){ # these unfortunately can't be vectorized
    theta_var_curr <- unwind_chol(theta_var[i,], reverse = T)
    proposals_curr <- matrix(proposals[i,], ncol = info$n_pars, byrow = T)
    group_ll <- sum(dmvnorm(proposals_curr, theta_mu[i,], theta_var_curr, log = T))
    prior_var <- log(robust_diwish(theta_var_curr, v=prior$v+ info$n_pars-1, S = 2*prior$v*diag(1/theta_a[i,])))
    prior_a <- sum(logdinvGamma(exp(theta_a[i,]), shape = 1/2,rate=1/(prior$A^2)))
    jac_var <- log(2^info$n_pars)+sum((info$n_pars + 1)*log(diag(theta_var_curr))) # Log of derivative of cholesky transformation
    sum_out[i] <- group_ll + prior_var + prior_a + jac_var
  }
  prior_mu <- dmvnorm(theta_mu, mean = prior$theta_mu_mean, sigma = prior$theta_mu_var, log =T)
  jac_a <- rowSums(theta_a)
  return(sum_out + prior_mu + jac_a) # Output is of length nrow(proposals)
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

