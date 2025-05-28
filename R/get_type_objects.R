get_objects <- function(type, selection = NULL, sample_prior = F, design = NULL, sampler = NULL,
                        prior = NULL, stage = 'sample', N = 1e5, ...){
  return_prior <- ifelse(is.null(sampler), TRUE, FALSE)
  if(type %in% c("standard", "blocked", "diagonal")){
    out <- get_objects_standard(selection, sample_prior, return_prior, design, prior, stage, N = N,
                                sampler, ...)
  }
  else if(type == "single"){
    out <- get_objects_single(selection, sample_prior, return_prior, design, prior, stage, N = N,
                                sampler, ...)
  }
  else if(type == "diagonal-gamma"){
    out <- get_objects_diag_gamma(selection, sample_prior, return_prior, design, prior, stage, N = N,
                            sampler, ...)
  }
  else if(type == "infnt_factor"){
    out <- get_objects_infnt_factor(selection, sample_prior, return_prior, design, prior, stage, N = N,
                               sampler, ...)
  }
  else if(type == "factor"){
    out <- get_objects_factor(selection, sample_prior, return_prior, design, prior, stage, N = N,
                                    sampler,...)
  }
  else if(type == "SEM"){
    out <- get_objects_SEM(selection, sample_prior, return_prior, design, prior, stage, N = N,
                              sampler,...)
  }
  else{
    stop("make sure type is supported!")
  }
  return(out)
}

add_prior_names <- function(prior, design, ...){
  dots <- list(...)
  pnames <- names(sampled_pars(design))
  # If there is a group_design add its names
  gnames <- add_group_par_names(pnames, list(...)$group_design)
  nuisance <- is.element(seq_len(length(pnames)), unique(c(dots$nuisance, dots$nuisance_non_hyper)))
  # This might be a weird place to set the nuisance parameters in the prior,
  # but it's the most efficient for now
  for(pri in names(prior)){
    if(is.null(dim(prior[[pri]]))){
      if(length(pnames) == length(prior[[pri]])){
        prior[[pri]] <- prior[[pri]][!nuisance]
        names(prior[[pri]]) <- pnames[!nuisance]
      }
      if(length(gnames) == length(prior[[pri]])){
        prior[[pri]] <- prior[[pri]][!nuisance]
        names(prior[[pri]]) <- gnames[!nuisance]
      }
    } else{
      if(length(pnames) == nrow(prior[[pri]])){
        prior[[pri]] <- prior[[pri]][!nuisance,,drop=FALSE]
        rownames(prior[[pri]]) <- pnames[!nuisance]
      }
      if(length(pnames) == ncol(prior[[pri]])){
        prior[[pri]] <- prior[[pri]][,!nuisance,drop=FALSE]
        colnames(prior[[pri]]) <- pnames[!nuisance]
      }
      if(length(gnames) == nrow(prior[[pri]])){
        prior[[pri]] <- prior[[pri]][!nuisance,,drop=FALSE]
        rownames(prior[[pri]]) <- gnames[!nuisance]
      }
      if(length(gnames) == ncol(prior[[pri]])){
        prior[[pri]] <- prior[[pri]][,!nuisance,drop=FALSE]
        colnames(prior[[pri]]) <- gnames[!nuisance]
      }
    }
  }
  return(prior)
}

get_objects_diag_gamma <- function(selection, sample_prior, return_prior, design = NULL,
                             prior = NULL, stage = 'sample', N = 1e5, sampler = NULL,...){
  acc_selection <- c("mu", "sigma2", "alpha", "LL", "Sigma")
  if(return_prior & !sample_prior){
    if(is.null(list(...)$return_info)) prior$prior <- get_prior_diag_gamma(design = design, sample = F, prior = prior)
    prior$descriptions <- list(
      theta_mu_mean = "mean of the group-level mean prior",
      theta_mu_var = "variance of the group-level mean prior",
      shape = "shape of the group-level variance prior",
      rate = "rate of the group-level variance prior"
    )
    prior$types <- list(
      mu = c("theta_mu_mean", "theta_mu_var"),
      Sigma = c("shape", "rate")
    )
    prior$type_descriptions <- list(
      mu = "Group-level mean",
      Sigma = 'Group-level covariance matrix'
    )
    if(!is.null(list(...)$return_info)) return(prior[c("types", "type_descriptions", "descriptions")])
    prior$prior <- add_prior_names(prior$prior, design, ...)
    return(prior)
  } else{
    if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
    if(sample_prior){
      if(selection == "alpha" & !is.null(sampler)){
        mu <- get_pars(sampler, selection = "mu", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        var <- get_pars(sampler, selection = "Sigma", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        sub_names <- names(sampler[[1]]$data)
        sampler <- list(list(samples =  list(alpha = get_alphas(mu, var, sub_names))))
      } else{
        sampler <- list(list(samples = get_prior_diag_gamma(prior = prior, design = design, selection = selection,N = N)))
      }
      sampler[[1]]$prior <- prior
      class(sampler) <- "emc"
      return(sampler)
    }
    idx <- get_idx(sampler, stage)
    return(get_base(sampler, idx, selection))
  }
}


get_objects_standard <- function(selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, stage = 'sample', N = 1e5, sampler = NULL, ...){
  acc_selection <- c("mu", "sigma2", "beta", "covariance", "correlation", "alpha", "Sigma", "LL")
  if(return_prior & !sample_prior){
    if(is.null(list(...)$return_info)) prior$prior <- get_prior_standard(design = design, sample = F, prior = prior, group_des = list(...)$group_design)
    prior$descriptions <- list(
      theta_mu_mean = "mean of the group-level mean prior",
      theta_mu_var = "variance of the group-level mean prior",
      beta_mean = "mean of the group-level regressors",
      beta_var = "variance of the group-level regressors",
      v = "degrees of freedom on the group-level (co-)variance prior, 2 leads to uniform correlations. Single value",
      A = "scale on the group-level variance prior, larger values lead to larger variances"
    )
    prior$types <- list(
      mu = c("theta_mu_mean", "theta_mu_var"),
      Sigma = c("v", "A"),
      beta = c("beta_mean", "beta_var")
    )
    prior$type_descriptions <- list(
      mu = "Group-level mean",
      Sigma = 'Group-level covariance matrix',
      beta = "Group-level regressors"
    )
    if(!is.null(list(...)$return_info)) return(prior[c("types", "type_descriptions", "descriptions")])
    prior$prior <- add_prior_names(prior$prior, design, ...)
    return(prior)
  } else{
    if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
    if(sample_prior){
      if(selection == "alpha" & !is.null(sampler)){
        mu <- get_pars(sampler, selection = "mu", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        var <- get_pars(sampler, selection = "Sigma", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        sub_names <- names(sampler[[1]]$data)
        sampler <- list(list(samples =  list(alpha = get_alphas(mu, var, sub_names))))
      } else{
        sampler <- list(list(samples = get_prior_standard(prior = prior, design = design, selection = selection,N = N,
                                                          par_groups = sampler[[1]]$par_groups, group_des = sampler[[1]]$group_design)))
      }
      sampler[[1]]$prior <- prior
      class(sampler) <- "emc"
      return(sampler)
    }
    idx <- get_idx(sampler, stage)
    if(selection == "mu"){
      return(lapply(sampler, implied_mean, idx))
    }
    return(get_base(sampler, idx, selection))
  }
}

implied_mean <- function(sampler, idx){
  beta <- sampler$samples$theta_mu[,idx]
  N <- ncol(beta)
  group_des <- sampler$group_designs
  par_names <- sampler$par_names
  if(is.null(par_names)){
    group_des <- sampler$samples$group_designs
    par_names <- sampler$samples$par_names
  }
  n_pars <- length(par_names)
  # Calculate mu (implied means) if we have group designs
  if(!is.null(group_des)) {
    # Prepare group designs in proper format
    group_designs <- add_group_design(par_names, group_des, N)

    # Calculate mean design matrices
    mean_designs <- calculate_mean_design(group_designs, n_pars)

    # Calculate implied means for each sample
    mu <- matrix(0, nrow = n_pars, ncol = N)
    for (i in 1:N) {
      mu[, i] <- calculate_implied_means(mean_designs, beta[, i], n_pars)
    }
    rownames(mu) <- par_names
  } else {
    # If no group designs, mu and beta are the same
    mu <- beta
  }
  return(mu)
}



get_idx <- function(sampler, stage){
  if(is.null(sampler[[1]]$samples$stage)){
    dims <- dim(sampler[[1]][[1]][[1]])
    idx <- 1:(dims[length(dims)])
  } else{
    if(sum(chain_n(sampler)) == 0) stop("Run fit to collect MCMC samples before using this function")
    idx <- which(sampler[[1]]$samples$stage %in% stage)
  }
  if(length(idx) == 0) stop("Make sure there are already samples of the selected stage")
  return(idx)
}

get_objects_single <- function(selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, stage = 'sample', N = 1e5, sampler = NULL,...){
  acc_selection <- c("alpha", "LL")
  if(return_prior & !sample_prior){
    if(is.null(list(...)$return_info)) prior$prior <- get_prior_single(design = design, sample = F, prior = prior)
    prior$descriptions <- list(
      theta_mu_mean = "mean of the prior",
      theta_mu_var = "variance of the prior"
    )
    prior$types <- list(
      alpha = c("theta_mu_mean", "theta_mu_var")
    )
    prior$type_descriptions <- list(
      alpha = "Subject-level prior"
    )
    if(!is.null(list(...)$return_info)) return(prior[c("types", "type_descriptions", "descriptions")])
    prior$prior <- add_prior_names(prior$prior, design, ...)
    return(prior)
  } else{
    if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
    if(sample_prior){
      sampler <- list(list(samples = get_prior_single(prior = prior, design = design, selection = selection,N = N)))
      sampler[[1]]$prior <- prior
      class(sampler) <- "emc"
      return(sampler)
    }
    idx <- get_idx(sampler, stage)
    return(get_base(sampler, idx, selection))
  }
}

get_objects_factor <- function(selection, sample_prior, return_prior, design = NULL,
                                     prior = NULL, stage = 'sample', N = 1e5, sampler = NULL, ...){
  acc_selection <- c("mu", "sigma2", "covariance", "correlation", "alpha", "Sigma", "loadings", "residuals", "LL")
  if(return_prior & !sample_prior){
    if(is.null(list(...)$return_info)) prior$prior <- do.call(get_prior_factor, c(list(design = design, sample = F, prior = prior), fix_dots(list(...), get_prior_factor)))
    prior$descriptions <- list(
      theta_mu_mean = "mean of the group-level mean prior",
      theta_mu_var = "variance of the group-level mean prior",
      theta_lambda_var = "variance of the factor loadings",
      as = "shape of inverse-gamma prior on the residual variances",
      bs = "rate of inverse-gamma prior on the residual variances",
      ap = "shape prior of inverse gamma on factor variances",
      bp = "rate prior of inverse gamma on factor variances"
    )
    prior$types <- list(
      mu = c("theta_mu_mean", "theta_mu_var"),
      loadings = c("theta_lambda_var", "ap", "bp"),
      residuals = c("as", "bs")
    )
    prior$type_descriptions <- list(
      mu = "Group-level mean",
      loadings = "Factor loadings",
      residuals = "Residual errors on the variances"
    )
    if(!is.null(list(...)$return_info)) return(prior[c("types", "type_descriptions", "descriptions")])
    prior$prior <- add_prior_names(prior$prior, design, ...)
    return(prior)
  } else{
    if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
    if(sample_prior){
      if(selection == "alpha" & !is.null(sampler)){
        mu <- get_pars(sampler, selection = "mu", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        var <- get_pars(sampler, selection = "Sigma", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        sub_names <- names(sampler[[1]]$data)
        sampler <- list(list(samples =  list(alpha = get_alphas(mu, var, sub_names))))
      } else{
        dots <- list(...)
        if(!is.null(sampler)){
          dots <- add_defaults(dots, Lambda_mat = attr(sampler[[1]], "Lambda_mat"))
        }
        sampler <- list(list(samples = do.call(get_prior_factor,
                                               c(list(prior = prior, design = design,
                                                      selection = selection,N = N), fix_dots(dots, get_prior_factor)))))
      }
      sampler[[1]]$prior <- prior
      class(sampler) <- "emc"
      return(sampler)
    }
    idx <- get_idx(sampler, stage)
    if(selection == "loadings"){
      return(lapply(sampler, FUN = function(x) return(x$samples$theta_lambda[,,idx])))
    }
    if(selection == "residuals"){
      return(lapply(sampler, FUN = function(x) return(1/x$samples$theta_sig_err_inv[,idx])))
    }
    return(get_base(sampler, idx, selection))
  }
}



get_objects_infnt_factor <- function(selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, stage = 'sample', N = 1e5, sampler = NULL, ...){
  acc_selection <- c("mu", "sigma2", "covariance", "correlation", "alpha", "Sigma", "loadings", "residuals", "LL")
  if(return_prior & !sample_prior){
    if(is.null(list(...)$return_info)) prior$prior <- do.call(get_prior_infnt_factor, c(list(design = design, sample = F, prior = prior), fix_dots(list(...), get_prior_infnt_factor)))
    prior$descriptions <- list(
      theta_mu_mean = "mean of the group-level mean prior",
      theta_mu_var = "variance of the group-level mean prior",
      as = "shape of inverse-gamma prior on the residual variances",
      bs = "rate of inverse-gamma prior on the residual variances",
      df = "shape and rate prior on loadings (local) shrinkage parameter",
      ad1 = "shape prior on factor loading variances of first column",
      bd1 = "rate prior on factor loading variances of first column",
      ad2 = "multiplicative shape prior on factor loading variances of subsequent columns",
      bd2 = "multiplicative rate prior on factor loading variances of subsequent columns"
    )
    prior$types <- list(
      mu = c("theta_mu_mean", "theta_mu_var"),
      loadings = c("df", "ad1", "bd1",  "ad2", "bd2"),
      residuals = c("as", "bs")
    )
    prior$type_descriptions <- list(
      mu = "Group-level mean",
      loadings = "Factor loadings",
      residuals = "Residual errors on the variances"
    )
    if(!is.null(list(...)$return_info)) return(prior[c("types", "type_descriptions", "descriptions")])
    prior$prior <- add_prior_names(prior$prior, design, ...)
    return(prior)
  } else{
    if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
    if(sample_prior){
      if(selection == "alpha" & !is.null(sampler)){
        mu <- get_pars(sampler, selection = "mu", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        var <- get_pars(sampler, selection = "Sigma", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        sub_names <- names(sampler[[1]]$data)
        sampler <- list(list(samples =  list(alpha = get_alphas(mu, var, sub_names))))
      } else{
        sampler <- list(list(samples = do.call(get_prior_infnt_factor,
                                               c(list(prior = prior, design = design,
                                                      selection = selection,N = N), fix_dots(list(...), get_prior_infnt_factor)))))
      }
      sampler[[1]]$prior <- prior
      class(sampler) <- "emc"
      return(sampler)
    }
    idx <- get_idx(sampler, stage)
    if(selection == "loadings"){
      return(lapply(sampler, FUN = function(x) return(x$samples$theta_lambda[,,idx])))
    }
    if(selection == "residuals"){
      return(lapply(sampler, FUN = function(x) return(1/x$samples$theta_sig_err_inv[,idx])))
    }
    return(get_base(sampler, idx, selection))
  }
}


get_objects_SEM <- function(selection, sample_prior, return_prior, design = NULL,
                               prior = NULL, stage = 'sample', N = 1e5, sampler = NULL, ...){
  acc_selection <- c("mu", "sigma2", "covariance", "alpha", "correlation", "Sigma", "loadings", "residuals",
                     "factor_residuals", "regressors", "factor_regressors", "structural_regressors",
                     "mu_implied", "LL")
  if(return_prior & !sample_prior){
    if(is.null(list(...)$return_info)) prior$prior <- do.call(get_prior_SEM, c(list(design = design, sample = F, prior = prior), fix_dots(list(...), get_prior_SEM)))
    prior$descriptions <- list(
      theta_mu_mean = "mean of the group-level mean prior",
      theta_mu_var = "variance of the group-level mean prior",
      lambda_var = "variance of the factor loadings",
      K_var = "variance of the parameter regressors",
      G_var = "variance of the factor regressors",
      B_var = "variance of structural regressors",
      a_d = "shape prior of inverse gamma/inverse wishart on factor variances",
      b_d = "rate prior of inverse gamma/inverse wishart on factor variances",
      a_e = "shape prior of inverse gamma on residuals",
      b_e = "rate prior of inverse gamma on residuals"
    )
    prior$types <- list(
      mu = c("theta_mu_mean", "theta_mu_var"),
      loadings = c("theta_lambda_var"),
      residuals = c("a_e", "b_e"),
      factor_residuals = c("a_d", "b_d"),
      regressors = c("K_var"),
      factor_regressors = c("G_var"),
      structural_regressors = c("B_var")
    )
    prior$type_descriptions <- list(
      mu = "group-level mean",
      loadings = "factor variances",
      residuals = "residuals on parameter variances",
      factor_residuals = "residuals on factor variances",
      regressors = "regressors on parameters",
      factor_regressors = "regressors on factors",
      structural_regressors = "structural regressors between factors"
    )
    if(!is.null(list(...)$return_info)) return(prior[c("types", "type_descriptions", "descriptions")])
    prior$prior <- add_prior_names(prior$prior, design, ...)
    return(prior)
  } else{
    if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
    if(sample_prior){
      if(selection == "alpha" & !is.null(sampler)){
        mu <- get_pars(sampler, selection = "mu_implied", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        var <- get_pars(sampler, selection = "Sigma", stage = stage, map = FALSE, return_mcmc = FALSE, merge_chains = TRUE, ...)
        sub_names <- names(sampler[[1]]$data)
        sampler <- list(list(samples =  list(alpha = get_alphas(mu, var, sub_names))))
      } else{
        dots <- list(...)
        if(!is.null(sampler)){
          dots <- add_defaults(dots, K_mat = attr(sampler[[1]], "K_mat"), B_mat = attr(sampler[[1]], "B_mat"),
                               covariates = sampler[[1]]$covariates,
                               Lambda_mat = attr(sampler[[1]], "Lambda_mat"), G_mat = attr(sampler[[1]], "G_mat"))
        }
        sampler <- list(list(samples = do.call(get_prior_SEM,
                                               c(list(prior = prior, design = design,
                                                      selection = selection,N = N), fix_dots(dots, get_prior_SEM)))))
      }
      sampler[[1]]$prior <- prior
      class(sampler) <- "emc"
      return(sampler)
    }
    idx <- get_idx(sampler, stage)
    if(selection == "residuals"){
      return(lapply(sampler, FUN = function(x){
        resids <- x$samples$epsilon_inv[,,idx]
        for(i in 1:dim(resids)[3]){
          resids[,,i] <- solve(resids[,,i])
        }
        return(resids)
      }))
    }
    if(selection == "factor_residuals"){
      return(lapply(sampler, FUN = function(x){
        resids <- x$samples$delta_inv[,,idx]
        for(i in 1:dim(resids)[3]){
          resids[,,i] <- solve(resids[,,i])
        }
        return(resids)
      }))
    }
    if(selection == "loadings"){
      return(lapply(sampler, FUN = function(x) return(x$samples$lambda[,,idx, drop = F])))
    }
    if(selection == "regressors"){
      return(lapply(sampler, FUN = function(x) return(x$samples$K[,,idx, drop = F])))
    }
    if(selection == "factor_regressors"){
      return(lapply(sampler, FUN = function(x) return(x$samples$G[,,idx, drop = F])))
    }
    if(selection == "structural_regressors"){
      return(lapply(sampler, FUN = function(x) return(x$samples$B[,,idx, drop = F])))
    }
    if(selection == "mu_implied"){
      return(lapply(sampler, FUN = get_mu_implied, idx))
    }
    return(get_base(sampler, idx, selection))
  }
}

get_base <- function(sampler, idx, selection){
  if(selection == "alpha"){
    return(lapply(sampler, FUN = function(x) return(x$samples$alpha[,,idx, drop = F])))
  } else if(selection == "LL"){
    return(lapply(sampler, FUN = function(x) return(x$samples$subj_ll[,idx, drop = F])))
  } else if(selection == "mu"){
    return(lapply(sampler, FUN = function(x) return(x$samples$theta_mu[,idx, drop = F])))
  } else if(selection == "beta"){
    return(lapply(sampler, FUN = function(x) {
      if (!is.null(x$samples$theta_beta)) {
        return(x$samples$theta_beta[,idx, drop = F])
      } else {
        return(x$samples$theta_mu[,idx, drop = F])
      }
    }))
  } else if(selection == "covariance"){
    return(lapply(sampler, FUN = function(x){
      out <- x$samples$theta_var[,,idx, drop = F]
      for(i in 1:dim(out)[3]){
        diag(out[,,i]) <- 0
      }
      return(out)
    }))
  } else if(selection == "sigma2"){
    return(lapply(sampler, FUN = function(x){
      out <- x$samples$theta_var[,,idx, drop = F]
      out <- apply(out,3,diag)
      return(out)
    }))
  }
  else if(selection == "Sigma"){
    return(lapply(sampler, FUN = function(x) return(x$samples$theta_var[,,idx, drop = F])))
  }
  else if(selection == "correlation"){
    return(lapply(sampler, FUN = function(x) return(
      array(apply(x$samples$theta_var[,,idx],3,cov2cor),dim=dim(x$samples$theta_var[,,idx, drop = F]),
            dimnames=dimnames(x$samples$theta_var)))))
  }
}


get_alphas <- function(mu, var, sub_names, N = ncol(mu)){
  n_pars <- nrow(mu)
  alpha <- array(NA_real_, dim = c(n_pars, length(sub_names), N))
  for(i in 1:N){
    alpha[,,i] <- t(rmvnorm(1, mu[,i], var[,,i]))
  }
  rownames(alpha) <- rownames(mu)
  colnames(alpha) <- sub_names
  return(alpha)
}




