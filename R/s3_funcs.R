#' @export
print.emc <- function(x, ...){
  n_chain <- chain_n(x)
  cat("Iterations: \n")
  print(chain_n(x))
  sub_names <- names(x[[1]]$data)
  cat("\n")
  cat("Subjects: \n")
  print(sub_names)
  par_names <- x[[1]]$par_names
  cat("\n")
  cat("Parameters: \n")
  print(par_names)
  return(invisible(list(iterations = n_chain, subjects = sub_names, parameters = par_names)))
}

#' Summary statistics for EMC2 samplers objects
#'
#' Computes quantiles, `Rhat` and `ESS` for selected model parameters.
#'
#' Note that if `selection = alpha` is used, summary statistics are computed
#' at the individual level. Only the first subject's summary output is printed
#' to the console but summary statistics for all subjects are returned by the
#' function.
#'
#' If `subfilter` is set to an integer other than 0 (i.e., the default), all iterations
#' up to that integer are ignored. If a vector is supplied, only those iterations
#' are considered.
#'
#' @param object An object of class `emc`
#' @param selection A character string indicating the parameter type
#' Defaults to `mu`, `sigma2`, and `alpha`.
#' (i.e., individual-level parameters). See below for more information.
#' @param probs The quantiles to be computed. Defaults to the the 2.5%, 50% and 97.5% quantiles.
#' @param filter A character string indicating the sampling stage to be summarized.
#' Can be `preburn`, `burn`, `adapt`, or `sample`.
#' @param subfilter An integer or vector, defaults to 0. See below
#' for more details.
#' @param digits An integer specifying rounding of output.
#' @param ... Further optional arguments.
#'
#' @return A list of summary output.
#' @export
summary.emc <- function(object, selection = c("mu", "sigma2", "alpha"), probs = c(0.025, .5, 0.975),
                        digits = 3, ...){
  args <- list(...)
  for (name in names(args) ) {
    assign(name, args[[name]])
  }
  if(attr(object[[1]], "variant_funs")$type == "single"){
    selection <- "alpha"
  }
  out_list <- list()
  for(select in selection){
    stats <- get_summary_stat(object, fun = c(get_posterior_quantiles, gelman_diag_robust, effectiveSize), probs = probs,
                     stat_name = c(paste0(probs*100, "%"), "Rhat", "ESS"), selection = select, ...)
    for(i in 1:length(stats)){
      if(length(stats) > 1){
        cat("\n", paste0(select, " ", names(stats)[i]), "\n")
      } else{
        cat("\n", names(stats)[i], "\n")
      }
      stat <- round(stats[[i]], digits)
      stat[,ncol(stat)] <- round(stat[,ncol(stat)])
      print(stat)
    }
    out_list[[select]] <- stat
  }
  return(invisible(out_list))
}


#' Plot function for EMC2 samplers objects
#'
#' Plots posterior densities and traceplots for EMC2 model parameters.
#'
#' Note that if `subfilter` is set to an integer other than 0 (i.e., the default), all iterations
#' up to that integer are ignored. If a vector is supplied, only those iterations
#' are considered.
#'
#' @param x An object of class `emc`
#' @param filter A character string indicating the sampling stage to be summarized.
#' Can be `preburn`, `burn`, `adapt`, or `sample`.
#' @param selection A character string indicating the parameter group.
#' Defaults to `mu`, `sigma2`, and `alpha`.
#' @param subfilter An integer or vector, defaults to 0. See below for more details.
#' @param subject An integer or character vector. Plot samples selected subjects.
#' @param thin An integer, only consider samples that are a multiple of `thin`.
#' @param ... Further optional arguments.
#' @export
plot.emc <- function(x, filter = "sample", selection = c("mu", "sigma2", "alpha"), ...){
  if(attr(x[[1]], "variant_funs")$type == "single"){
    selection <- "alpha"
  }
  for(select in selection){
    plot_chains_new(x, filter = filter, selection = select, subfilter = subfilter,
                    subject = subject, thin = thin, ...)
  }
}
