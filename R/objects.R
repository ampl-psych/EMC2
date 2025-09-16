filter_obj <- function(obj, idx){
  dims <- dim(obj)
  dim_names <- dimnames(obj)
  if(is.null(dims)) return(obj)
  if(length(dims) == 2){
    if(nrow(obj) == ncol(obj)){
      if(nrow(obj) > 1){
        if(mean(abs(abs(rowSums(obj/max(obj))) - abs(colSums(obj/max(obj))))) < .01) return(obj)
      }
    }
  }
  obj <- obj[slice.index(obj, length(dims)) %in% idx]
  dims[length(dims)] <- length(idx)
  dim(obj) <- dims
  dimnames(obj) <- dim_names # Give back to the community
  return(obj)
}

remove_samples <- function(samples, stage = "sample", filter = NULL, thin = 1,
                           length.out = NULL, keep_stages = FALSE){
  if(!any(samples$samples$stage %in% stage)) stop("No samples from the chosen stage are present")
  filter_idx <- which(samples$samples$stage %in% stage)
  max_dim <- max(filter_idx); min_dim <- min(filter_idx)
  if(!is.null(filter)){
    filter <- min_dim + filter
    if(length(filter) == 1){
      if(filter > max_dim) stop("filter larger than available iterations")
      filter_idx <- filter:max_dim
    } else{
      if(length(filter) > max_dim) stop("filter larger than available iterations")
      filter_idx <- filter
    }
  }
  if(!is.null(length.out)){
    if(length.out > length(filter_idx)) stop("Length.out longer than available samples")
    filter_idx <- round(seq(min_dim, max_dim, length.out = length.out))
  } else if(thin > 1){
    filter_idx <- round(seq(min_dim, max_dim, by = thin))
  }
  if(keep_stages) filter_idx <- c(which(!(samples$samples$stage %in% stage)), filter_idx)

  if(any(samples$nuisance)) {
    samples$sampler_nuis$samples <- base::rapply(samples$sampler_nuis$samples, f = function(x) filter_obj(x, filter_idx), how = "replace")
    samples$sampler_nuis$samples$idx <- length(filter_idx)
  }
  stage <- samples$samples$stage
  samples$samples <- base::rapply(samples$samples, f = function(x) filter_obj(x, filter_idx), how = "replace")
  samples$samples$stage <- stage[filter_idx]
  samples$samples$idx <- length(filter_idx)
  return(samples)
}

#' Merge Samples
#'
#' Merges samples from all chains as one unlisted object.
#'
#' Note that all sampling stages are included in the merged output,
#' including iterations from the `preburn`, `burn`, and `adapt` stages.
#' `merge_chains(emc)$samples$stage` shows the corresponding sampling stages.
#'
#' @param emc An emc object, commonly the output of `fit()`
#'
#' @return An unlisted emc object with all chains merged
#' @export
#'
merge_chains <- function(emc){
  out_samples <- emc[[1]]
  # Only thing that differs between the chains is the samples$samples
  sampled_objects <- lapply(emc, FUN = function(x) return(x$samples))
  keys <- unique(unlist(lapply(sampled_objects, names)))
  sampled_objects <- setNames(do.call(mapply, c(abind, lapply(sampled_objects, '[', keys))), keys)
  sampled_objects$idx <- sum(sampled_objects$idx)
  out_samples$samples <- sampled_objects
  if(any(out_samples$nuisance)){
    sampled_objects <- lapply(emc, FUN = function(x) return(x$sampler_nuis$samples))
    keys <- unique(unlist(lapply(sampled_objects, names)))
    sampled_objects <- setNames(do.call(mapply, c(abind, lapply(sampled_objects, '[', keys))), keys)
    out_samples$sampler_nuis$samples <- sampled_objects
  }
  return(out_samples)
}

#' MCMC Chain Iterations
#'
#' Returns a matrix with the number of samples per chain for each stage that is present
#' in the emc object (i.e., `preburn`, `burn`, `adapt`,
#' `sample`). The number of rows of the matrix reflects the number of chains
#' and the number of columns the number of sampling stages.
#'
#' @param emc A list, the output of `fit()`.
#'
#' @return A matrix
#' @examples
#' chain_n(samples_LNR)
#'
#' @export

chain_n <- function(emc)
{
  do.call(rbind,lapply(emc, function(x){
    table(factor(x$samples$stage,levels=c("preburn", "burn","adapt","sample")))
  }))
}

extract_samples <- function(sampler, stage = c("adapt", "sample"), max_n_sample = NULL, n_chains) {
  type <- sampler$type
  samples <- sampler$samples
  nuis_type <- sampler$sampler_nuis$type
  if("sample" %in% stage & !is.null(max_n_sample)){
    sample_filter <- which(samples$stage %in% "sample" & seq_along(samples$stage) <= samples$idx)
    adapt_filter <- which(samples$stage %in% "adapt" & seq_along(samples$stage) <= samples$idx)
    if(length(sample_filter) > max_n_sample){
      chain_idx <- seq(1, length(sample_filter), by = length(sample_filter)/n_chains)
      prob_filter <- sample_filter - rep(sample_filter[chain_idx], each = length(sample_filter)/n_chains)
      sample_filter <- sample(sample_filter, max_n_sample, prob = prob_filter/max(prob_filter))
    }
    if(length(sample_filter) > length(adapt_filter)){
      full_filter <- sample_filter
    } else{
      full_filter <- c(adapt_filter, sample_filter)
    }
  } else{
    full_filter <- which(samples$stage %in% stage & seq_along(samples$stage) <= samples$idx)
  }
  if(any(sampler$nuisance)){
    nuisance <- sampler$nuisance
    suppressWarnings(sampler$sampler_nuis$samples$alpha <- sampler$samples$alpha[nuisance,,,drop = F])
    sampler$samples$alpha <- sampler$samples$alpha[!nuisance,,]
    out <- filtered_samples(sampler, full_filter, type = type)
    out$nuisance <- filtered_samples(sampler$sampler_nuis, full_filter, type = nuis_type)
  } else{
    out <-  filtered_samples(sampler, full_filter, type = type)
  }
  return(out)
}

filter_emc <- function(samples, thin = 1, length.out = NULL, filter = NULL){
  max_dim <- dim(samples)[length(dim(samples))]
  if(!is.null(filter)){
    if(length(filter) == 1){
      if(filter > max_dim) stop("filter larger than available iterations")
      keep <- filter:max_dim
    } else{
      if(length(filter) > max_dim) stop("filter larger than available iterations")
      keep <- filter
    }
    if(length(dim(samples)) > 2){
      samples <- samples[,,keep, drop = F]
    } else{
      samples <- samples[,keep, drop = F]
    }
  }

  max_dim <- dim(samples)[length(dim(samples))]
  if(!is.null(length.out)){
    if(length.out > max_dim) stop("Length.out longer than available samples")
    keep <- round(seq(1, max_dim, length.out = length.out))
  } else if(thin > 1){
    keep <- seq(1, max_dim, by = thin)
  }
  if(thin > 1 || !is.null(length.out)){
    if(length(dim(samples)) > 2){
      samples <- samples[,,keep, drop = F]
    } else{
      samples <- samples[,keep, drop = F]
    }
  }

  return(samples)
}


# filter for constants or duplicates
filter_const_and_dup <- function(samples, remove_dup = TRUE, remove_constants = TRUE){
  if(!remove_dup & !remove_constants) return(samples)
  # We only need the first list entry to calculate the idx
  if(length(dim(samples[[1]])) == 2){
    x <- do.call(cbind, samples)
    is_constant <- apply(x[,1:ncol(x), drop = F], 1, sd) == 0
    is_duplicate <- duplicated(round(rowSums(x[,1:ncol(x), drop = F]), 8))
    if(remove_dup){
      filter <- is_duplicate | is_constant
    } else{
      filter <- is_constant
    }
    samples <- lapply(samples, function(x) x[!filter,,drop = F])
  } else{
    x <- do.call(abind, samples)
    is_constant <- apply(x[,,1:dim(x)[3], drop = F], 1:2, sd) == 0
    # Add all the samples together, if any of them are the same (up till 8 digits)
    # We should probably assume that they are a duplicated entry
    # This is useful for correlations and such (on which samples are usually mirrored)
    all_sums <- c(apply(x[,,1:dim(x)[3], drop = F], 1:2, sum))
    is_duplicate <- duplicated(round(all_sums/mean(all_sums, na.rm = TRUE), 8))
    if(remove_dup){
      filter <- is_duplicate | is_constant
    } else{
      filter <- is_constant
    }
    samples <- lapply(samples, FUN = function(x){
      out <- x
      for(i in 1:dim(out)[3]){
        tmp <- out[,,i, drop = F]
        tmp[filter] <- NA
        out[,,i] <- tmp
      }
      return(out)
    })
  }

  return(samples)
}

# Two horrible functions that generalize very well though!
make_3dim_mcmc <- function(x, split){
  samples <- lapply(x, FUN = function(y){
    lapply(split_along_dim(y, split), FUN = function(z){
      idx <- !is.na(z[,1])
      coda::mcmc(t(z[idx,, drop = F]))
    })
  })
  out <- setNames(lapply(1:length(samples[[1]]), function(i) coda::as.mcmc.list(
    lapply(samples, function(m) m[[i]]))),  names(samples[[1]]))
  out <- lapply(out, function(x){
    if(ncol(x[[1]]) == 0) return(NULL) else return(x)
  })
  out[sapply(out, is.null)] <- NULL
  return(out)
}

# Function to split array a into list along any dimension n
split_along_dim <- function(a, n){
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
}

filter_sub_and_par <- function(obj, sub, sub_names, par){
  par_names <- c(colnames(obj), rownames(obj))
  if(is.character(sub)){
    if(any(!(sub %in% sub_names))){
      stop("Make sure you specify the subject names correctly")
    }
  }
  if(!is.null(par)){
    if(any(!(par %in% par_names))){
      stop("Make sure you specify the parameter names correctly")
    }
  }

  if(length(dim(obj)) > 2){
    is_sub_idx <- sapply(dimnames(obj), FUN = function(x) identical(x, sub_names))
    if(any(is_sub_idx)){
      if(is_sub_idx[2]){ # Second dimension has subject names
        if(is.character(sub)){
          obj <- obj[,sub_names %in% sub,,drop = F]
        }
        if(is.numeric(sub)){
          obj <- obj[,sub,,drop = F]
        }
      } else{ # Must be the first dimension
        if(is.character(sub)){
          obj <- obj[sub_names %in% sub,,,drop = F]
        }
        if(is.numeric(sub)){
          obj <- obj[sub,,,drop = F]
        }
      }
    }
    if(!is.null(par)){
      if(any(colnames(obj) %in% par)){
        obj <- obj[,colnames(obj) %in% par,,drop = F]
      } else if(any(rownames(obj) %in% par)){
        obj <- obj[rownames(obj) %in% par,,,drop = F]
      }
    }
  } else{
    if(all(rownames(obj) %in% sub_names)){
      if(is.character(sub)){
        obj <- obj[sub_names %in% sub,, drop = F]
      }
      if(is.numeric(sub)){
        obj <- obj[sub,,drop = F]
      }
    }
    if(!is.null(par)){
      obj <- obj[rownames(obj) %in% par,,drop = F]
    }
  }
  return(obj)
}

# Returns the last ran stage of an emc object
get_last_stage <- function(emc){
  nstage <- colSums(chain_n(emc))
  if(all(nstage == 0)){
    stage <- "preburn"
  } else{
    has_ran <- nstage[nstage != 0]
    stage <- names(has_ran)[length(has_ran)]
  }
  return(stage)
}


add_defaults <- function(dots, ...){
  extra_args <- list(...)
  arg_names <- names(extra_args)
  for(i in 1:length(extra_args)){
    cur_name <- arg_names[i]
    if(is.null(dots[[cur_name]])) dots[[cur_name]] <- extra_args[[i]]
  }
  return(dots)
}

fix_dots_plot <- function(dots, exclude = ""){
  dots <- dots[!(names(dots) %in% exclude)]
  dots <- dots[names(dots) %in% c(names(par()), names(formals(arrows)), names(formals(plot.default)))]
  return(dots)
}

fix_dots <- function(dots, fun, exclude = "", consider_dots = TRUE){
  dots <- dots[!(names(dots) %in% exclude)]
  fun_args <- names(formals(fun))
  if("..." %in% fun_args & consider_dots){
    return(dots)
  } else{
    return(dots[names(dots) %in% fun_args])
  }
}

#' Filter/Manipulate Parameters from emc Object
#'
#' Underlying function used in most plotting and object handling functions in
#' EMC2. Can for example be used to filter/thin a parameter type
#' (i.e, group-level means ``mu``) and convert to an mcmc.list.
#'
#' @param emc an emc object.
#' @param selection A Character string. Indicates which parameter type to select (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param stage A character string. Indicates from which sampling stage(s) to take the samples from (i.e. `preburn`, `burn`, `adapt`, `sample`)
#' @param thin An integer. By how much to thin the chains
#' @param filter Integer or numeric vector. If an integer is supplied, iterations
#' up until that integer are removed. If a vector is supplied, the iterations
#' within the range are kept.
#' @param length.out Integer. Alternatively to thinning, you can also select a
#' desired length of the MCMC chains, which will be thinned appropriately.
#' @param map Boolean. If `TRUE` parameters will be mapped back to the cells of
#' the experimental design using the design matrices.
#' Otherwise the sampled parameters are returned.
#' Only works for `selection = mu` or `selection = alpha`.
#' @param add_recalculated Boolean. If `TRUE` will also add recalculated parameters,
#' such as b in the LBA (b = B + A; see `?LBA`), or z in the DDM z = Z*A (see `?DDM`)
#' only works when `map = TRUE`
#' @param by_subject Boolean. If `TRUE` for selections that include subject parameters (e.g. `alpha`),
#' plot/stats are organized by subject, otherwise by parameter.
#' @param return_mcmc Boolean. If `TRUE` returns an mcmc.list object, otherwise a matrix/array with the parameter type.
#' @param merge_chains Boolean. If `TRUE` returns parameter type merged across chains.
#' @param subject Integer (vector) or character (vector). If an integer will select the 'x'th subject(s),
#' if a character it should match subject names in the data which will be selected.
#' @param flatten Boolean. If `FALSE` for 3-dimensional samples (e.g., correlations: n-pars x n-pars x iterations).
#' organizes by the dimension containing parameter names, otherwise collapses names across the first and second dimension.
#' Does not apply for `selection = "alpha"`
#' @param remove_dup Boolean. If `TRUE` removes duplicate values from the samples. Automatically set to `TRUE` if `flatten = TRUE`
#' @param remove_constants Boolean. If `TRUE` removes constant values from the samples (e.g. 0s in the covariance matrix).
#' @param use_par Character (vector). If specified, only these parameters are returned. Should match the parameter names
#' (i.e. these are collapsed when `flatten = TRUE` and use_par should also be collapsed names).
#' @param type Character indicating the group-level model selected. Only necessary if sampler isn't specified.
#' @param true_pars Set of `true_parameters` can be specified to apply flatten or use_par on a set of true parameters
#' @param covariates Only needed with `plot` for priors and covariates in the design
#' @param chain Integer. Which of the chain(s) to return
#'
#' @return An mcmc.list object of the selected parameter types with the specified manipulations
#' @export
#'
#' @examples
#' # E.g. get the group-level mean parameters mapped back to the design
#' get_pars(samples_LNR, stage = "sample", map = TRUE, selection = "mu")
#'
#' # Or return the flattened correlation, with 10 iterations per chain
#' get_pars(samples_LNR, stage = "sample", selection = "correlation", flatten = TRUE, length.out = 10)
get_pars <- function(emc,selection= "mu", stage=get_last_stage(emc),thin=1,filter=0,
                    map = FALSE, add_recalculated = FALSE, length.out = NULL,
                    by_subject = FALSE, return_mcmc = TRUE, merge_chains = FALSE,
                    subject = NULL, flatten = FALSE, remove_dup = FALSE,
                    remove_constants = TRUE, use_par = NULL, type = NULL,
                    true_pars = NULL, chain = NULL, covariates = NULL)
{
  if(!is.null(emc[[1]]$n_subjects)) emc <- restore_duplicates(emc)
  if(add_recalculated) map <- TRUE
  if(!(selection %in% c("mu", "alpha"))) map <- FALSE
  if(is.null(type)) type <- emc[[1]]$type
  # Centralized guard: mapped population-level parameters not yet correct for hierarchical models
  if(type != "single" && isTRUE(map) && selection %in% c("mu")){
    warning("map=TRUE is not yet correctly implemented for population-level parameters (i.e., `mu`, `sigma2`, `covariances` and `correlation`). See https://github.com/ampl-psych/EMC2/issues/119 Please use map = FALSE until this issue is fixed.")
  }
  if(type == "single" & !(selection %in% c("LL", "alpha"))) selection <- "alpha"
  samples <- get_objects(type = type, sampler = emc, stage = stage,
                         selection = selection)

  if(map){
    samples <- lapply(samples, do_map, get_design(emc), include_constants = FALSE,
                      add_recalculated = add_recalculated, covariates = covariates,
                      emc = emc)
  }
  if(flatten) remove_dup <- TRUE
  if(!is.null(true_pars)){ # Kluge to make sure the right object dimensions/filtering is performed on simulated parameters
    samples <- samples[1]
    if(length(dim(samples[[1]])) > 2){
      if(is.vector(true_pars)){
        samples[[1]][,,1] <- true_pars
      } else if(ncol(true_pars) == nrow(samples[[1]])){
        samples[[1]][,,1] <- t(true_pars)
      } else{
        samples[[1]][,,1] <- true_pars
      }
    } else{
      samples[[1]][,1] <- true_pars
    }
  }
  if(selection == "alpha"){
    flatten <- FALSE
    subnames <- colnames(samples[[1]])
  } else{
    subnames <- names(emc[[1]]$data)
  }
  if(length(dim(samples[[1]])) > 2 & flatten){
    if(is.null(rownames(samples[[1]]))){
      pnams <- colnames(samples[[1]])
    } else{
      pnams <- c(outer(rownames(samples[[1]]), colnames(samples[[1]]), paste, sep = "."))
    }
    samples <- lapply(samples, function(x) out <- apply(x,3,function(y){c(y)}))
    samples <- lapply(samples, function(x) {rownames(x) <- pnams; return(x)})
  }
  samples <- filter_const_and_dup(samples, remove_dup, remove_constants)


  samples <- lapply(samples, filter_sub_and_par, subject, subnames, use_par)
  samples <- lapply(samples, filter_emc, thin, length.out, filter)
  if(!is.null(chain)){
    if(any(!(chain %in% 1:length(samples)))) stop("chain selection exceeds number of chains")
    samples <- samples[chain]
  }
  if(merge_chains){
    if(length(dim(samples[[1]])) == 2){
      samples <- do.call(cbind, samples)
    } else{
      samples <- do.call(abind, samples)
    }
    if(return_mcmc) samples <- list(samples)
  }
  if(!return_mcmc) return(samples)
  if(length(dim(samples[[1]])) > 2){
    is_sub_idx <- sapply(dimnames(samples[[1]]), FUN = function(x) any(x %in% subnames))
    if(any(is_sub_idx) & by_subject){
      samples <- make_3dim_mcmc(samples, which(is_sub_idx))
    } else if(any(is_sub_idx)){
      samples <- make_3dim_mcmc(samples, which(!is_sub_idx)[1])
    } else{
      samples <- make_3dim_mcmc(samples, 2)
    }
  } else{
    samples <- list(coda::as.mcmc.list(lapply(samples, FUN = function(x) coda::mcmc(t(x)))))
    names(samples) <- selection
  }
  if(!is.null(true_pars)){ # Kluge to make sure the right object dimensions/filtering is performed on simulated parameters
    samples <- lapply(samples, FUN = function(x) return(list(x[[1]][1,, drop = F])))
  }
  return(samples)
}

concat_emc <- function(emc1, emc2, step_size, stage){
  out_samples <- emc2
  remove_idx <- ifelse(chain_n(emc1)[1,stage] != 0, 1, 0)
  emc2 <- subset(emc2, stage = c("preburn", "burn", "adapt", "sample"), filter = 1:(chain_n(emc2)[1,stage] - remove_idx))
  for(i in 1:length(emc1)){
    sampled_objects <- list(emc1[[i]]$samples, emc2[[i]]$samples)
    keys <- unique(unlist(lapply(sampled_objects, names)))
    sampled_objects <- do.call(mapply, c(abind, lapply(sampled_objects, '[', keys)))
    sampled_objects$idx <- sum(sampled_objects$idx)
    attr(sampled_objects, "pm_settings") <- attr(emc2[[i]]$samples, "pm_settings")
    out_samples[[i]]$samples <- sampled_objects
    out_samples[[i]]$samples$last_theta_var_inv <- emc2[[i]]$samples$last_theta_var_inv
    if(any(out_samples[[1]]$nuisance)){
      sampled_objects <- list(emc1[[i]]$sampler_nuis$samples, emc2[[i]]$sampler_nuis$samples)
      keys <- unique(unlist(lapply(sampled_objects, names)))
      sampled_objects <- do.call(mapply, c(abind, lapply(sampled_objects, '[', keys)))
      if(!is.list(sampled_objects)){ # deals with niche case where there's only item and this drops list making
        sampled_objects <- list(idx= sum(sampled_objects))
      } else{
        sampled_objects$idx <- sum(sampled_objects$idx)
      }
      out_samples[[i]]$sampler_nuis$samples <- sampled_objects
      out_samples[[i]]$sampler_nuis$samples$last_theta_var_inv <- emc2[[i]]$sampler_nuis$samples$last_theta_var_inv
    }
  }
  rm(emc1); rm(emc2)
  return(out_samples)
}


fit_remove_samples <- function(emc){
  if(chain_n(emc)[1,4] > chain_n(emc)[1,3]){
    emc <- subset(emc, stage = "sample")
  } else if (chain_n(emc)[1,4] > 0){
    emc <- subset(emc, stage = c("adapt", "sample"))
  } else if(chain_n(emc)[1,3] > 0){
    emc <- subset(emc, stage = c("burn", "adapt"))
  }
  return(emc)
}





