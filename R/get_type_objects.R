get_objects <- function(type, selection = NULL, sample_prior = F, return_prior = T, design = NULL, sampler = NULL,
                        prior = NULL, mapped = F, ...){
  if(type == "standard"){
    return(get_objects_standard(type, selection, sample_prior, return_prior, design, prior, mapped, ...))
  }
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

get_objects_standard <- function(type, selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, mapped = F, ...){
  if(return_prior){
    if(sample_prior){
      prior <- get_prior_standard(prior = prior, design = design, selection = selection, map = mapped)[[1]]
    } else{
      prior$prior <- get_prior_standard(design = design, sample = F)
      prior$descriptions <- list(
        theta_mu_mean = "Mean of the group-level mean prior",
        theta_mu_var = "Variance-covariance matrix of the group-level mean prior",
        v = "degrees of freedom on the group-level variance prior, 2 leads to uniform correlations. Single value",
        A = "scale on the group-level variance prior, larger values lead to larger variances"
      )
      prior$prior <- add_prior_names(prior$prior, design)
    }
    return(prior)
  }
}
