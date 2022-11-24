
single_out_joint <- function(joint_samples_list, i){
  single_samples_list <- lapply(joint_samples_list, return_single_sampler, i)
  attr(single_samples_list, "data_list") <-   attr(joint_samples_list,"data_list")[i]
  attr(single_samples_list, "design_list") <-   attr(joint_samples_list,"design_list")[i]
  attr(single_samples_list, "model_list") <-   attr(joint_samples_list,"model_list")[i]
  return(single_samples_list)
}

return_single_sampler <- function(joint_samples, i){
  par_names <- joint_samples$par_names
  prefix <- unique(gsub("[|].*", "", par_names))[i]
  idx <- grep(paste0(prefix, "|"), par_names, fixed = T)
  current_pars <- par_names[idx]
  replacement <-  gsub(".*[|]", "", current_pars)
  single_samples <- joint_samples
  single_samples$samples <- base::rapply(joint_samples$samples, f = function(x) fix_single_object(x, prefix, current_pars, replacement), how = "replace")
  single_samples$par_names <- replacement
  single_samples$data <- lapply(joint_samples$data, FUN = function(x) return(x[[i]]))
  single_samples$ll_func <- attr(single_samples$data[[1]], "model")$log_likelihood
  single_samples$prior$theta_mu_mean <- single_samples$prior$theta_mu_mean[idx]
  single_samples$prior$theta_mu_var <- single_samples$prior$theta_mu_var[idx, idx]
  return(single_samples)
}


fix_single_object <- function(object, prefix, current_pars, replacement){
  # Recursively go through every object in the list
  # Check whether the names match the prefix
  # Unleash all my R powers
  dims <- dim(object)
  n_dims <- length(dims)
  if(n_dims > 1){
    dim_names <- dimnames(object)
    for(i in 1:n_dims){
      tmp_idx <- dim_names[[i]] %in% current_pars
      if(any(tmp_idx)){
        object <- object[slice.index(object, i) %in% which(tmp_idx)]
        dims[i] <- sum(tmp_idx)
        dim(object) <- dims
        dim_names[[i]] <- replacement
        dimnames(object) <- dim_names # Give back to the community
      }
    }
  }
  return(object)
}
