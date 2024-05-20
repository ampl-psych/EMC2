add_info_blocked <- function(sampler, prior = NULL, ...){
  # blocked specific attributes
  n_pars <- sum(!(sampler$nuisance | sampler$grouped))
  par_groups <- list(...)$par_groups
  sampler$par_groups <- par_groups
  sampler$n_par_groups <- length(unique(par_groups))
  sampler$prior <- get_prior_blocked(prior, sum(!(sampler$nuisance | sampler$grouped)), sample = F)
  return(sampler)
}


#' Prior specification or prior sampling for blocked estimation
#'
#' Works analogous to `get_prior_standard`. Blocks of the covariance matrix to estimate
#' are only considered in sampling. To get the default prior for a created design:
#' `get_prior_diag(design = design, sample = FALSE)`
#'
#' For details see Huang, A., & Wand, M. P. (2013). Simple marginally noninformative
#' prior distributions for covariance matrices. *Bayesian Analysis*, 8, 439-452. https://doi.org/10.1214/13-BA815.
#'
#' Note that if `sample = FALSE`, prior$theta_mu_invar (the inverse of the prior covariance matrix on the group-level mean) is returned,
#' which is only used for computational efficiency
#'
#' @param prior A named list that can contain the prior mean (`theta_mu_mean`) and
#' variance (`theta_mu_var`) on the group-level mean, or the scale (`A`), or degrees of freedom (`v`)
#' for the group-level variance-covariance matrix. For `NULL` entries, the default prior is created
#' @param n_pars Often inferred from the design, but if `design = NULL`, `n_pars`
#' will be used to determine the size of prior.
#' @param sample Boolean, defaults to `TRUE`, sample from the prior or simply
#' return the prior specifications?
#' @param map Boolean, defaults to `TRUE`. If `sample = TRUE`, the implied prior is sampled.
#' This includes back-transformations for naturally bounded parameters such as
#' the non-decision time and an inverse mapping from the design matrix back to the
#' cells of the design. If `FALSE`, the transformed, unmapped, parameters are used.
#' Note that `map` does not affect the prior used in the sampling process.
#' @param N How many samples to draw from the prior, default 1e5
#' @param design The design obtained from `make_design`, required when map = TRUE
#' @param type  character. If sample = TRUE, what prior to sample from. Options: "mu", "variance", "covariance" "full_var", "alpha".
#'
#' @return A list with a single entry of type of samples from the prior (if sample = TRUE) or else a prior object
#' @examples \dontrun{
#' # First define a design for the model
#' design_DDMaE <- make_design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Now get the default prior
#' prior <- get_prior_blocked(design = design_DDMaE, sample = FALSE)
#' # We can change values in the default prior or use make_prior
#' # Then we can get samples from this prior e.g.
#' samples <- get_prior_blocked(prior = prior, design = design_DDMaE,
#'   sample = TRUE, type = "mu")
#' }
#' @export

get_prior_blocked <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL,
                               map = FALSE){
  # Checking and default priors
  if(is.null(prior)){
    prior <- list()
  }
  if(!is.null(design)){
    n_pars <- length(attr(design, "p_vector"))
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
  if(sample){
    out <- list()
    if(!selection %in% c("mu", "variance", "covariance", "correlation", "full_var")){
      stop("for variant standard, you can only specify the prior on the mean, variance, covariance or the correlation of the parameters")
    }
    if(selection == "mu"){
      samples <- mvtnorm::rmvnorm(N, mean = prior$theta_mu_mean,
                                  sigma = prior$theta_mu_var)
      if(!is.null(design)){
        colnames(samples) <- par_names <- names(attr(design, "p_vector"))
        if(map){
          samples <- map_mcmc(samples,design,design$model,include_constants=FALSE)
        }
      }
      out$mu <- samples
      return(out)
    } else {
      var <- array(NA_real_, dim = c(n_pars, n_pars, N))
      for(i in 1:N){
        a_half <- 1 / rgamma(n = n_pars,shape = 1/2,
                             rate = 1/(prior$A^2))
        attempt <- tryCatch({
          var[,,i] <- riwish(prior$v + n_pars - 1, 2 * prior$v * diag(1 / a_half))
        },error=function(e) e, warning=function(w) w)
        if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
          sample_idx <- sample(1:(i-1),1)
          var[,,i] <- var[,,sample_idx]
        }
      }
      if (selection == "variance") {
        vars_only <- t(apply(var,3,diag))
        if(!is.null(design)){
          colnames(vars_only) <- names(attr(design, "p_vector"))
        }
        out$variance <- vars_only
      }
      lt <- lower.tri(var[,,1])
      if (selection == "correlation"){
        corrs <- array(apply(var,3,cov2cor),dim=dim(var),dimnames=dimnames(var))
        out$correlation <- t(apply(corrs,3,function(x){x[lt]}))
      }
      if(selection == "covariance"){
        out$covariance <- t(apply(var,3,function(x){x[lt]}))
      }
      if (selection == "full_var"){
        out$full_var <- t(apply(var, 3, c))
      }
      return(out)
    }
  }
  return(prior)
}


get_startpoints_blocked <- function(pmwgs, start_mu, start_var){
  if (is.null(start_mu)) start_mu <- rmvnorm(1, mean = pmwgs$prior$theta_mu_mean, sigma = pmwgs$prior$theta_mu_var)
  # If no starting point for group var just sample some
  if (is.null(start_var)) {
    start_var <- matrix(nrow = 0, ncol = 0)
    for(i in 1:pmwgs$n_par_groups){
      # Check how many parameters are in the current group
      n_pars_group <- sum(pmwgs$par_groups == i)
      # Make a subblock of start variances for those parameters
      start_var <- adiag(start_var, riwish(n_pars_group * 3,diag(n_pars_group)))
    }
  }
  start_a_half <- 1 / rgamma(n = pmwgs$n_pars, shape = 2, rate = 1)
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var), a_half = start_a_half))
}

gibbs_step_blocked <- function(sampler, alpha){
  # Gibbs step for group means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  last <- last_sample_standard(sampler$samples)
  prior <- sampler$prior
  n_pars_total <- sampler$n_pars-sum(sampler$nuisance) - sum(sampler$grouped)
  tmu_out <- numeric()
  a_half_out <- numeric()
  tvar_out <- matrix(0, nrow = n_pars_total, ncol = n_pars_total)
  tvinv_out <- matrix(0, nrow = n_pars_total, ncol = n_pars_total)
  idx <- 0
  for(group in 1:sampler$n_par_groups){
    group_idx <- sampler$par_groups == group
    n_pars <- sum(group_idx)
    if(n_pars > 1){
      var_mu <- ginv(sampler$n_subjects * last$tvinv[group_idx, group_idx] + prior$theta_mu_invar[group_idx, group_idx])
      mean_mu <- as.vector(var_mu %*% (last$tvinv[group_idx, group_idx] %*% apply(alpha[group_idx,], 1, sum) +
                                         prior$theta_mu_invar[group_idx, group_idx] %*% prior$theta_mu_mean[group_idx]))
      chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
      # New sample for mu.
      tmu <- rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
      # New values for group var
      theta_temp <- alpha[group_idx,] - tmu
      cov_temp <- (theta_temp) %*% (t(theta_temp))
      B_half <- 2 * prior$v * diag(1 / last$a_half[group_idx]) + cov_temp # nolint
      tvar <- riwish(prior$v + n_pars - 1 + sampler$n_subjects, B_half) # New sample for group variance
      tvinv <- ginv(tvar)

      # Sample new mixing weights.
      a_half <- 1 / rgamma(n = n_pars,shape = (prior$v + n_pars) / 2,
                           rate = prior$v * diag(tvinv) + 1/(prior$v^2))
    } else{
      var_mu = 1.0 / (sampler$n_subjects * last$tvinv[group_idx, group_idx] + prior$theta_mu_invar[group_idx, group_idx])
      mean_mu = var_mu * (sum(alpha[group_idx,]) * last$tvinv[group_idx, group_idx] + prior$theta_mu_invar[group_idx, group_idx] * prior$theta_mu_mean[group_idx])
      tmu <- rnorm(n_pars, mean_mu, sd = sqrt(var_mu))
      tvinv = rgamma(n=n_pars, shape=prior$v/2 + sampler$n_subjects/2, rate=prior$v/last$a_half +
                       rowSums( (alpha-tmu)^2 ) / 2)
      tvar = 1/tvinv
      #Contrary to standard pmwg I use shape, rate for IG()
      a_half <- 1 / rgamma(n = n_pars, shape = (prior$v + 1) / 2,
                           rate = prior$v * tvinv + 1/(prior$v^2))
    }

    tmu_out <- c(tmu_out, tmu)
    a_half_out <- c(a_half_out, a_half)

    tvar_out[(idx+1):(idx+n_pars), (idx+1):(idx+n_pars)] <- tvar
    tvinv_out[(idx+1):(idx+n_pars), (idx+1):(idx+n_pars)] <- tvinv
    idx <- idx + n_pars
  }

  names(tmu_out) <- sampler$par_names

  return(list(tmu = tmu_out,tvar = tvar_out,tvinv = tvinv_out,a_half = a_half_out,alpha = alpha))
}

get_conditionals_blocked <- function(s, samples, n_pars, iteration = NULL, idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  pts2_unwound <- apply(samples$theta_var[idx,idx,],3,unwind)
  index <- rowMeans(pts2_unwound == 0) == 0 #remove all 0 elements
  pts2_unwound <- pts2_unwound[index,]
  all_samples <- rbind(samples$alpha[idx,s,],samples$theta_mu[idx,],pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- var(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[idx,iteration], unwind(samples$theta_var[idx,idx,iteration])[index]))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

get_all_pars_blocked <- function(samples, idx, info){
  stop("no IS2 for blocked estimation yet")
}

group_dist_blocked <- function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  stop("no IS2 for blocked estimation yet")

}

prior_dist_blocked <- function(parameters, info){
  stop("no IS2 for blocked estimation yet")

}
