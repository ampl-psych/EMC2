get_objects <- function(type, selection = NULL, sample_prior = F, design = NULL, sampler = NULL,
                        prior = NULL, mapped = F, N = 1e5, ...){
  return_prior <- ifelse(is.null(sampler), TRUE, FALSE)
  if(type == "standard"){
    out <- get_objects_standard(type, selection, sample_prior, return_prior, design, prior, mapped, N = N,
                                sampler, ...)
  }
  else if(type == "single"){
    out <- get_objects_single(type, selection, sample_prior, return_prior, design, prior, mapped, N = N,
                                sampler, ...)
  }
  else if(type == "diagonal"){
    out <- get_objects_diag(type, selection, sample_prior, return_prior, design, prior, mapped, N = N,
                              sampler, ...)
  }
  else if(type == "blocked"){
    out <- get_objects_blocked(type, selection, sample_prior, return_prior, design, prior, mapped, N = N,
                            sampler, ...)
  }
  else if(type == "infnt_factor"){
    out <- get_objects_infnt_factor(type, selection, sample_prior, return_prior, design, prior, mapped, N = N,
                               sampler, ...)
  }
  else if(type == "factor"){
    out <- get_objects_factor(type, selection, sample_prior, return_prior, design, prior, mapped, N = N,
                                    sampler,...)
  }
  else{
    stop("make sure type is supported!")
  }
  return(out)
}

add_prior_names <- function(prior, design){
  pnames <- names(sampled_p_vector(design))
  for(pri in names(prior)){
    if(is.null(dim(prior[[pri]]))){
      if(length(pnames) == length(prior[[pri]])){
        names(prior[[pri]]) <- pnames
      }
    } else{
      if(length(pnames) == nrow(prior[[pri]])){
        rownames(prior[[pri]]) <- pnames
      }
      if(length(pnames) == ncol(prior[[pri]])){
        colnames(prior[[pri]]) <- pnames
      }
    }
  }
  return(prior)
}

get_objects_diag <- function(type, selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, mapped = F, N = 1e5, sampler = NULL,...){
  acc_selection <- c("mu", "variance", "alpha", "full_var")
  if(return_prior){
    if(sample_prior){
      if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
      prior <- get_prior_diag(prior = prior, design = design, selection = selection, map = mapped,
                                  N = N)[[1]]
    } else{
      prior$prior <- get_prior_diag(design = design, sample = F, prior = prior)
      prior$descriptions <- list(
        theta_mu_mean = "mean of the group-level mean prior",
        theta_mu_var = "variance of the group-level mean prior",
        v = "degrees of freedom on the group-level variance prior",
        A = "scale on the group-level variance prior, larger values lead to larger variances"
      )
      prior$prior <- add_prior_names(prior$prior, design)
    }
    return(prior)
  }
}

get_objects_standard <- function(type, selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, mapped = F, N = 1e5, sampler = NULL, ...){
  acc_selection <- c("mu", "variance", "covariance", "correlation", "alpha", "full_var")
  if(return_prior){
    if(sample_prior){
      if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
      prior <- get_prior_standard(prior = prior, design = design, selection = selection, map = mapped,
                                  N = N)[[1]]
    } else{
      prior$prior <- get_prior_standard(design = design, sample = F, prior = prior)
      prior$descriptions <- list(
        theta_mu_mean = "mean of the group-level mean prior",
        theta_mu_var = "variance of the group-level mean prior",
        v = "degrees of freedom on the group-level (co-)variance prior, 2 leads to uniform correlations. Single value",
        A = "scale on the group-level variance prior, larger values lead to larger variances"
      )
      prior$prior <- add_prior_names(prior$prior, design)
    }
    return(prior)
  } else{
    if(selection == "mu"){
      return(lapply(sampler, FUN = function(x) return(x$samples$theta_mu)))
    } else if(selection == "alpha"){
      return(lapply(sampler, FUN = function(x) return(x$samples$alpha)))
    } else if(selection == "covariance"){
      return(lapply(sampler, FUN = function(x){
        out <- x$samples$theta_var
        for(i in 1:dim(out)[3]){
          diag(out[,,i]) <- 0
        }
        return(out)
      }))
    } else if(selection == "variance"){
      return(lapply(sampler, FUN = function(x){
        out <- x$samples$theta_var
        out <- apply(out,3,diag)
        return(out)
      }))
    }

  }
}


get_objects_blocked <- function(type, selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, mapped = F, N = 1e5, sampler = NULL,...){
  acc_selection <- c("mu", "variance", "covariance", "correlation", "alpha", "full_var")
  if(return_prior){
    if(sample_prior){
      if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
      prior <- get_prior_blocked(prior = prior, design = design, selection = selection, map = mapped,
                                  N = N, ...)[[1]]
    } else{
      prior$prior <- get_prior_blocked(design = design, sample = F, prior = prior, ...)
      prior$descriptions <- list(
        theta_mu_mean = "mean of the group-level mean prior",
        theta_mu_var = "variance of the group-level mean prior",
        v = "degrees of freedom on the group-level (co-)variance prior, 2 leads to uniform correlations. Single value",
        A = "scale on the group-level variance prior, larger values lead to larger variances"
      )
      prior$prior <- add_prior_names(prior$prior, design)
    }
    return(prior)
  }
}



get_objects_single <- function(type, selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, mapped = F, N = 1e5, sampler = NULL,...){
  acc_selection <- "alpha"
  if(return_prior){
    if(sample_prior){
      if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
      prior <- get_prior_single(prior = prior, design = design, selection = selection, map = mapped,
                                  N = N)[[1]]
    } else{
      prior$prior <- get_prior_single(design = design, sample = F, prior = prior)
      prior$descriptions <- list(
        theta_mu_mean = "mean of the prior",
        theta_mu_var = "variance of the prior"
      )
      prior$prior <- add_prior_names(prior$prior, design)
    }
    return(prior)
  }
}

get_objects_factor <- function(type, selection, sample_prior, return_prior, design = NULL,
                                     prior = NULL, mapped = F, N = 1e5, sampler = NULL, ...){
  acc_selection <- c("mu", "variance", "covariance", "correlation", "alpha", "full_var", "loadings", "residuals")
  if(return_prior){
    if(sample_prior){
      if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
      prior <- get_prior_factor(prior = prior, design = design, selection = selection, map = mapped,
                                      N = N, ...)[[1]]
    } else{
      prior$prior <- get_prior_factor(design = design, sample = F, prior = prior, ...)
      prior$descriptions <- list(
        theta_mu_mean = "mean of the group-level mean prior",
        theta_mu_var = "variance of the group-level mean prior",
        theta_lambda_var = "variance of the factor loadings",
        as = "shape of inverse-gamma prior on the residual variances",
        bs = "rate of inverse-gamma prior on the residual variances",
        ap = "shape prior of inverse gamma on factor variances",
        bp = "rate prior of inverse gamma on factor variances"
      )
      prior$prior <- add_prior_names(prior$prior, design)
    }
    return(prior)
  } else{
    if(selection == "mu"){
      return(lapply(sampler, FUN = function(x) return(x$samples$theta_mu)))
    }

  }
}



get_objects_infnt_factor <- function(type, selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, mapped = F, N = 1e5, sampler = NULL, ...){
  acc_selection <- c("mu", "variance", "covariance", "correlation", "alpha", "full_var", "loadings", "residuals")
  if(return_prior){
    if(sample_prior){
      if(!selection %in% acc_selection) stop(paste0("selection must be in : ", paste(acc_selection, collapse = ", ")))
      prior <- get_prior_infnt_factor(prior = prior, design = design, selection = selection, map = mapped,
                                  N = N, ...)[[1]]
    } else{
      prior$prior <- get_prior_infnt_factor(design = design, sample = F, prior = prior, ...)
      prior$descriptions <- list(
        theta_mu_mean = "mean of the group-level mean prior",
        theta_mu_var = "variance of the group-level mean prior",
        as = "shape of inverse-gamma prior on the residual variances",
        bs = "rate of inverse-gamma prior on the residual variances",
        df = "shape and rate prior on cross-loadings (local) shrinkage parameter",
        ad1 = "shape prior on factor loading variances of first column",
        bd1 = "rate prior on factor loading variances of first column",
        ad2 = "multiplicative shape prior on factor loading variances of subsequent columns",
        bd2 = "multiplicative rate prior on factor loading variances of subsequent columns"
      )
      prior$prior <- add_prior_names(prior$prior, design)
    }
    return(prior)
  } else{
    if(selection == "mu"){
      return(lapply(sampler, FUN = function(x) return(x$samples$theta_mu)))
    }

  }
}




