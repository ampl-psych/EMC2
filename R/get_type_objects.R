get_objects <- function(type, selection, sample_prior = F, return_prior = T, design = NULL,
                        prior = NULL, mapped = F){
  if(type == "standard"){
    return(get_objects_standard(type, selection, sample_prior, return_prior, design, prior, mapped))
  }
}

get_objects_standard <- function(type, selection, sample_prior, return_prior, design = NULL,
                                 prior = NULL, mapped = F){
  if(return_prior){
    if(sample_prior){
      prior <- get_prior_standard(prior = prior, design = design, selection = selection, mapped = mapped)[[1]]
    } else{
      prior <- get_prior_standard(design = design, sample = F)
      prior$descriptions <- list(
        theta_mu_mean = "Mean of the group-level mean prior",
        theta_mu_var = "Variance of the group-level mean prior",
        v = "degrees of freedom on the group-level variance prior, 2 leads to uniform correlations",
        A = "scale on the group-level variance prior, larger values lead to larger variances"
      )
    }
    return(prior)
  }
}
