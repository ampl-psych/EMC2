
single_out_joint <- function(joint_samples_list, i){
  single_samples_list <- lapply(joint_samples_list, return_single_sampler, i)
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
  # single_samples$model <- list(joint_samples$model[[i]])
  single_samples$model <- lapply(joint_samples$model, FUN = function(x) return(x[[i]]))
  single_samples$prior <- fix_single_prior(single_samples$prior, idx)
  return(single_samples)
}

fix_single_prior <- function(prior, idx){
  out <- prior
  prior_names <- names(prior)
  for(name in prior_names){
    obj <- prior[[name]]
    if(length(obj) > 1){
      if(length(dim(obj)) == 2){
        out[[name]] <- obj[idx, idx]
      } else{
        out[[name]] <- obj[idx]
      }
    }
  }
  return(out)
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
