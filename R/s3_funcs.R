#' @export
print.emc <- function(x, ...){
  n_chain <- chain_n(x)
  print(chain_n(x))
  return(invisible(n_chain))
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
#' @param selection A character string indicating the parameter group.
#' Defaults to `mu`, `variance`, and `correlation`. Other options are `covariance`, and `alpha`
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
summary.emc <- function(object, selection = c("mu", "variance", "correlation"), probs = c(0.025, .5, 0.975),
                        filter = "sample", subfilter = 0, digits = 3, ...){
  args <- list(...)
  for (name in names(args) ) {
    assign(name, args[[name]])
  }
  if(attr(object[[1]], "variant_funs")$type == "single"){
    selection <- "alpha"
  }

  out_list <- list()
  for(select in selection){
    par_df <- parameters_data_frame(object, filter = filter, subfilter = subfilter, selection = select)
    Rhat <- gd_pmwg(object, filter = filter, subfilter = subfilter, selection = select, print_summary = F, omit_mpsrf = T)
    ESS <- es_pmwg(object, filter = filter, subfilter = subfilter, selection = select, print_summary = F, summary_alpha = NULL)
    if(select == "alpha"){
      unq_subs <- as.character(unique(par_df$subjects))
      for(i in 1:length(unq_subs)){
        tmp_df <- par_df[par_df$subjects == unq_subs[i],-1]
        quants <- t(apply(tmp_df, 2, quantile, probs))
        combined <- cbind(quants, Rhat[i,], ESS[i,])
        colnames(combined)[c(4,5)] <- c("Rhat", "ESS")
        out_list[[unq_subs[i]]] <- combined
        if(i == 1){
          cat("\n", unq_subs[i], "\n")
          print(round(combined,digits))
        }
      }
    } else{
      quants <- t(apply(par_df, 2, quantile, probs))
      combined <- cbind(quants, Rhat, ESS)
      cat("\n", select, "\n")
      print(round(combined,digits))
      out_list[[select]] <- combined
    }
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
#' Defaults to `mu`, `variance`, and `correlation`. Other options are `covariance`, and `alpha`
#' (i.e., individual-level parameters).
#' @param subfilter An integer or vector, defaults to 0. See below for more details.
#' @param subject An integer or character vector. Plot samples selected subjects.
#' @param thin An integer, only consider samples that are a multiple of `thin`.
#' @param ... Further optional arguments.
#' @export
plot.emc <- function(x, filter = "sample", selection = c("mu", "variance", "correlation"),
                     subfilter = 0, subject = NULL, thin = 1,  ...){

  layout <- NULL
  ylim <- NULL
  args <- list(...)
  for (name in names(args) ) {
    assign(name, args[[name]])
  }

  subject <- ifelse(is.null(subject), 1, subject)
  input <- x
  if(attr(x[[1]], "variant_funs")$type == "single") selection <- "alpha"
  for(select in selection){
    pmwg_mcmc <- as_mcmc.list(input,selection=select,filter=filter,
                              thin=thin,subfilter=subfilter)
    auto.layout <- is.null(layout)
    if (!auto.layout) par(mfrow=layout)
    if (attr(pmwg_mcmc,"selection")=="alpha" || attr(pmwg_mcmc,"selection")=="random") {
      snams <- names(pmwg_mcmc)
      if(is.null(subject)) subject <- snams[1]
      if (is.numeric(subject)) subject <- snams[subject]
      if (!all(subject %in% names(pmwg_mcmc)))
        stop("Subject not present\n")
      for (i in subject) {
        plot(pmwg_mcmc[[i]],auto.layout=auto.layout,density=TRUE,
             xlab=paste0("Iterations ",i),ask=FALSE,trace=TRUE,
             ylab=attr(pmwg_mcmc,"selection"),smooth=FALSE,ylim=ylim)
      }
    } else if (attr(pmwg_mcmc,"selection") %in% c("LL","epsilon")) {
      if (any(is.na(subject))) subject <- names(pmwg_mcmc)
      if (!all(subject %in% names(pmwg_mcmc)))
        stop("Subject not present\n")
      for (i in subject) {
            plot(pmwg_mcmc[[i]],auto.layout=auto.layout,density=TRUE,
                 xlab=paste0("Iterations ",i),ask=FALSE,trace=TRUE,
                 ylab=attr(pmwg_mcmc,"selection"),smooth=FALSE,ylim=ylim)
      }
    } else {
          plot(pmwg_mcmc,auto.layout=auto.layout,density=TRUE,
               ask=FALSE,trace=TRUE,ylim=ylim,
               ylab=attr(pmwg_mcmc,"selection"),smooth=FALSE)
    }
  }
}
