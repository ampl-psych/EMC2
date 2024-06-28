

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

remove_samples <- function(samples, filter = "sample", subfilter = NULL, thin = 1,
                           length.out = NULL){

  if(!any(samples$samples$stage %in% filter)) stop("No samples from the chosen filter are present")
  filter_idx <- which(samples$samples$stage %in% filter)
  max_dim <- max(filter_idx); min_dim <- min(filter_idx)
  if(!is.null(subfilter)){
    subfilter <- min_dim + subfilter
    if(length(subfilter) == 1){
      if(subfilter > max_dim) stop("Subfilter larger than available iterations")
      filter_idx <- subfilter:max_dim
    } else{
      if(length(subfilter) > max_dim) stop("Subfilter larger than available iterations")
      filter_idx <- subfilter
    }
  }
  if(!is.null(length.out)){
    if(length.out > length(filter_idx)) stop("Length.out longer than available samples")
    filter_idx <- round(seq(min_dim, max_dim, length.out = length.out))
  } else if(thin > 1){
    filter_idx <- round(seq(min_dim, max_dim, by = thin))
  }
  if(any(samples$nuisance)) {
    samples$sampler_nuis$samples <- base::rapply(pmwg$sampler_nuis$samples, f = function(x) filter_obj(x, filter_idx), how = "replace")
    samples$sampler_nuis$samples$idx <- length(filter_idx)
  }
  samples$samples <- base::rapply(samples$samples, f = function(x) filter_obj(x, filter_idx), how = "replace")
  samples$samples$stage <- samples$samples$stage[filter_idx]
  samples$samples$idx <- length(filter_idx)
  return(samples)
}

subset_emc <- function(samplers, filter = "sample", subfilter = NULL, thin = 1,
                              length.out = NULL){
  design_list <- attr(samplers, "design_list")
  data_list <- attr(samplers, "data_list")
  model_list <- attr(samplers, "model_list")
  samplers <- lapply(samplers, remove_samples, filter = filter, subfilter = subfilter,
                     thin = thin, length.out = length.out)
  attr(samplers, "design_list") <- design_list
  attr(samplers, "data_list") <- data_list
  attr(samplers, "model_list") <- model_list
  class(samplers) <- "emc"
  return(samplers)
}




#' Merge samples
#'
#' Merges samples from all chains as one unlisted object.
#'
#' Note that all sampling stages are included in the merged output,
#' including iterations from the `preburn`, `burn`, and `adapt` stages.
#' `merge_samples(samplers)$samples$stage` shows the corresponding sampling stages.
#'
#' @param samplers A samplers object, commonly the output of `run_emc()`
#'
#' @return An unlisted samplers object with all chains merged
#' @export
#'
merge_samples <- function(samplers){
  out_samples <- samplers[[1]]
  # Only thing that differs between the chains is the samples$samples
  sampled_objects <- lapply(samplers, FUN = function(x) return(x$samples))
  keys <- unique(unlist(lapply(sampled_objects, names)))
  sampled_objects <- setNames(do.call(mapply, c(abind, lapply(sampled_objects, '[', keys))), keys)
  sampled_objects$idx <- sum(sampled_objects$idx)
  out_samples$samples <- sampled_objects
  if(any(out_samples$nuisance)){
    sampled_objects <- lapply(samplers, FUN = function(x) return(x$sampler_nuis$samples))
    keys <- unique(unlist(lapply(sampled_objects, names)))
    sampled_objects <- setNames(do.call(mapply, c(abind, lapply(sampled_objects, '[', keys))), keys)
    out_samples$sampler_nuis$samples <- sampled_objects
  }
  return(out_samples)
}

#' chain_n()
#'
#' Returns a matrix with the number of samples per chain for each stage that is present
#' in the sampling object (i.e., `preburn`, `burn`, `adapt`,
#' `sample`). The number of rows of the matrix reflects the number of chains
#' and the number of columns the number of sampling stages.
#'
#' @param samplers A list, the output of `run_emc()`.
#'
#' @return A matrix
#' @examples \dontrun{
#' chain_n(samplers)
#' }
#' @export

chain_n <- function(samplers)
{
  do.call(rbind,lapply(samplers, function(x){
    table(factor(x$samples$stage,levels=c("preburn", "burn","adapt","sample")))
  }))
}

extract_samples <- function(sampler, stage = c("adapt", "sample"), max_n_sample = NULL, n_chains) {
  variant_funs <- attr(sampler, "variant_funs")
  samples <- sampler$samples
  type <- sampler$sampler_nuis$type
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
    nuisance <- sampler$nuisance[!sampler$grouped]
    sampler$sampler_nuis$samples$alpha <- sampler$samples$alpha[nuisance,,,drop = F]
    sampler$samples$alpha <- sampler$samples$alpha[!nuisance,,]
    out <- variant_funs$filtered_samples(sampler, full_filter)
    out$nuisance <- get_variant_funs(type)$filtered_samples(sampler$sampler_nuis, full_filter)
  } else{
    if(!is.null(sampler$g_map_fixed)){
      out <- 1# filtered_samples_lm(sampler, full_filter)
    } else{
      out <- variant_funs$filtered_samples(sampler, full_filter)
    }
  }
  return(out)
}

#' get data
#'
#' Extracts data from an emc object
#'
#' @param samplers an emc object
#'
#' @return a data frame
get_data <- function(samplers) {
  dat <- do.call(rbind,lapply(samplers[[1]]$data,\(x) x[attr(x,"expand"),]))
  row.names(dat) <- NULL
  dat <- dat[dat$lR == levels(dat$lR)[1],]
  dat <- dat[,!(colnames(dat) %in% c("trials","lR","lM","winner"))]
  dat
}

#' Converts a pmwgs object (or list of such) to a data frame.
#'
#' @param samples A list of samplers or samplers converted to mcmc objects.
#' @param selection String designating parameter type (e.g. mu, sigma2, correlation, alpha)
#' @param filter A string. Specifies which stage you want to plot.
#' @param thin An integer. Keep only iterations that are a multiple of thin.
#' @param subfilter An integer or vector. If integer it will exclude up until
#' that integer. If vector it will include everything in that range.
#' @param mapped Boolean (default FALSE) if TRUE plot parameters mapped to design
#' otherwise sampled parameters
#' @param include_constants Include parameters that are not sampled (i.e., constants)
#' @param stat A function that is applied to each column of the data frame
#' @param subject_n Number of alpha iterations to retain for each subject (default NULL = all)
#'
#' @return A data frame with one row for each sample (with a subjects column if selection = "alpha")
#' @export

parameters_data_frame <- function(samples,selection = "mu", N = 1000,
                                  stat=NULL,...)
  # extracts and stacks chains into a matrix
{
  dots <- list(...)
  dots$merge_chains <- TRUE ; dots$return_mcmc <- FALSE
  dots$flatten <- TRUE; dots$length.out <- N/3
  out <- do.call(as_mcmc_new, c(list(samples,selection=selection), fix_dots(dots, as_mcmc_new)))
  if(selection == "alpha"){
    out <- aperm(out, c(3,1,2))
    out <- apply(out, 3, identity, simplify=FALSE)
    for(i in 1:length(out)){
      out[[i]] <- as.data.frame(out[[i]])
      out[[i]] <- cbind(names(out)[i], out[[i]])
      colnames(out[[i]])[1] <- "subjects"
      out[[i]] <- out[[i]][1:N,]
    }
    out <- do.call(rbind, out)
    rownames(out) <- NULL
  } else{
    out <- as.data.frame(t(out))
    out <- out[1:N,]
  }
  out
}

filter_emc <- function(samples, thin = 1, length.out = NULL, subfilter = NULL){
  max_dim <- dim(samples)[length(dim(samples))]
  if(!is.null(subfilter)){
    if(length(subfilter) == 1){
      if(subfilter > max_dim) stop("Subfilter larger than available iterations")
      keep <- subfilter:max_dim
    } else{
      if(length(subfilter) > max_dim) stop("Subfilter larger than available iterations")
      keep <- subfilter
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
filter_const_and_dup <- function(samples, remove_dup = TRUE){
  # We only need the first list entry to calculate the idx
  x <- samples[[1]]
  if(length(dim(x)) == 2){
    is_constant <- apply(x[,1:(min(50, ncol(x)))], 1, sd) == 0
    is_duplicate <- duplicated(round(rowSums(x[,1:(min(50, ncol(x)))]), 4))
    if(remove_dup){
      filter <- is_duplicate | is_constant
    } else{
      filter <- is_constant
    }
    samples <- lapply(samples, function(x) x[!filter,,drop = F])
  } else{
    is_constant <- apply(x[,,1:(min(50, ncol(x))), drop = F], 1:2, sd) == 0
    all_sums <- c(apply(x[,,1:(min(50, ncol(x))), drop = F], 1:2, sum))
    is_duplicate <- duplicated(round(all_sums/mean(all_sums), 4))
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


#' Filter/manipulate parameters from emc object
#'
#' Underlying function used in most plotting and object handling functions in
#' EMC2. Can for example be used to filter/thin a parameter type
#' (i.e, group-level means ``mu``) and convert to an mcmc.list.
#'
#' @param sampler an emc object.
#' @param selection A Character string. Indicates which parameter type to select (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param stage A character string. Indicates from which sampling stage(s) to take the samples from (i.e. `preburn`, `burn`, `adapt`, `sample`)
#' @param thin An integer. By how much to thin the chains
#' @param subfilter Integer or numeric vector. If an integer is supplied, iterations
#' up until that integer are removed. If a vector is supplied, the iterations
#' within the range are kept.
#' @param length.out Integer. Alternatively to thinning, you can also select a
#' desired length of the MCMC chains, which will be thinned appropriately.
#' @param map Boolean. If `TRUE` parameters will be mapped back to the cells of
#' the experimental design using the design matrices.
#' Otherwise the sampled parameters are returned.
#' Only works for `selection = mu` or `selection = alpha`.
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
#' @param use_par Character (vector). If specified, only these parameters are returned. Should match the parameter names
#' (i.e. these are collapsed when `flatten = TRUE` and use_par should also be collapsed names).
#' @param type Character indicating the group-level model selected. Only necessary if sampler isn't specified.
#' @param true_pars Set of `true_parameters` can be specified to apply flatten or use_par on a set of true parameters
#' @param chain Integer. Which of the chain(s) to return
#'
#' @return
#' @export
#'
#' @examples
#' # E.g. get the group-level mean parameters mapped back to the design
#' as_mcmc_new(samplers_LNR, filter = "sample", map = TRUE, selection = "mu")
#'
#' # Or return the flattened correlation, with 10 iterations per chain
#' as_mcmc_new(samplers_LNR, filter = "sample", selection = "correlation", flatten = TRUE, length.out = 10)
as_mcmc_new <- function(sampler,selection= "mu", filter="sample",thin=1,subfilter=0,
                        map = FALSE, length.out = NULL, by_subject = FALSE,
                        return_mcmc = TRUE, merge_chains = FALSE,
                        subject = NULL, flatten = FALSE, remove_dup = FALSE,
                        use_par = NULL, type = NULL,
                        true_pars = NULL, chain = NULL)
{
  # filter <- which(sampler$samples$stage %in% filter)
  if(!(selection %in% c("mu", "alpha"))) map <- FALSE
  if(is.null(type)) type <- attr(sampler[[1]], "variant_funs")$type
  samples <- get_objects(type = type, sampler = sampler, filter = filter,
                         selection = selection)

  if(map){
    samples <- lapply(samples, map_mcmc, attr(sampler,"design_list")[[1]], include_constants = FALSE)
  }
  if(flatten) remove_dup <- TRUE
  if(!is.null(true_pars)){ # Kluge to make sure the right object dimensions/filtering is performed on simulated parameters
    samples <- samples[1]
    if(length(dim(samples[[1]])) > 2){
      samples[[1]][,,1] <- true_pars
    } else{
      samples[[1]][,1] <- true_pars
    }
  }

  subnames <- names(sampler[[1]]$data)
  if(selection == "alpha") flatten = FALSE

  if(length(dim(samples[[1]])) > 2 & flatten){
    if(is.null(rownames(samples[[1]]))){
      pnams <- colnames(samples[[1]])
    } else{
      pnams <- c(outer(rownames(samples[[1]]), colnames(samples[[1]]), paste, sep = "."))
    }
    samples <- lapply(samples, function(x) out <- apply(x,3,function(y){c(y)}))
    samples <- lapply(samples, function(x) {rownames(x) <- pnams; return(x)})
  }
  samples <- filter_const_and_dup(samples, remove_dup)
  subnames <- names(sampler[[1]]$data)

  samples <- lapply(samples, filter_sub_and_par, subject, subnames, use_par)
  samples <- lapply(samples, filter_emc, thin, length.out, subfilter)
  if(!is.null(chain)){
    if(any(!(chain %in% 1:length(samples)))) stop("chain selection exceeds number of chains")
    samples <- samples[chain]
  }
  if(is.null(attr(sampler[[1]], "variant_funs")$type)) subnames <- "alpha"
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




