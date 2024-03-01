#' @export
print.emc <- function(x, ...){
  n_chain <- chain_n(x)
  print(chain_n(x))
  return(invisible(n_chain))
}

#' @export
summary.emc <- function(x, ...){
  args <- list(...)
  filter <- ifelse(is.null(args$filter), "sample", args$filter)
  subfilter <- ifelse(is.null(args$subfilter), 0, args$subfilter)
  if(is.null(args$selection)){
    selection <- "mu"
  } else{
    selection <- args$selection
  }
  if(is.null(args$probs)){
    probs <- c(0.025, .5, 0.975)
  } else{
    probs <- args$probs
  }
  out_list <- list()
  for(select in selection){
    par_df <- parameters_data_frame(x, filter = filter, subfilter = subfilter, selection = select)
    Rhat <- gd_pmwg(x, filter = filter, subfilter = subfilter, selection = select, print_summary = F, omit_mpsrf = T)
    ESS <- es_pmwg(x, filter = filter, subfilter = subfilter, selection = select, print_summary = F, summary_alpha = NULL)
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
          print(combined)
        }
      }
    } else{
      quants <- t(apply(par_df, 2, quantile, probs))
      combined <- cbind(quants, Rhat, ESS)
      cat("\n", select, "\n")
      print(combined)
      out_list[[select]] <- combined
    }
  }
  return(invisible(out_list))
}


#' @export
plot.emc <- function(x, ...){
  args <- list(...)
  filter <- ifelse(is.null(args$filter), "sample", args$filter)
  subfilter <- ifelse(is.null(args$subfilter), 0, args$subfilter)
  layout <- args$layout
  plot_acf <- ifelse(is.null(args$plot_acf), FALSE, args$plot_acf)
  subject <- ifelse(is.null(args$subject), 1, args$subject)
  thin <- ifelse(is.null(args$thin), 1, args$thin)
  ylim <- args$ylim
  if(is.null(args$selection)){
    selection <- "mu"
  } else{
    selection <- args$selection
  }
  input <- x
  for(select in selection){
    pmwg_mcmc <- as_mcmc.list(input,selection=select,filter=filter,
                              thin=thin,subfilter=subfilter)
    auto.layout <- any(is.na(layout))
    no_layout <- is.null(layout)
    if (!auto.layout & !no_layout) par(mfrow=layout)
    if (attr(pmwg_mcmc,"selection")=="alpha" || attr(pmwg_mcmc,"selection")=="random") {
      snams <- names(pmwg_mcmc)
      if(is.null(subject)) subject <- snams[1]
      if (is.numeric(subject)) subject <- snams[subject]
      if (!all(subject %in% names(pmwg_mcmc)))
        stop("Subject not present\n")
      for (i in subject) {
        if (!auto.layout & !no_layout) par(mfrow=layout)
        if (plot_acf) for (j in dimnames(pmwg_mcmc[[i]])[[2]]) {
          acf(pmwg_mcmc[[i]][,j],main=paste0("Chain ",acf_chain,": ","s",i,": ",j))
        } else {
          if (!auto.layout & !no_layout) par(mfrow=layout)
          plot(pmwg_mcmc[[i]],auto.layout=auto.layout,density=TRUE,
               xlab=paste0("Iterations ",i),ask=FALSE,trace=TRUE,
               ylab=attr(pmwg_mcmc,"selection"),smooth=FALSE,ylim=ylim)
        }
      }
    } else if (attr(pmwg_mcmc,"selection") %in% c("LL","epsilon")) {
      if (any(is.na(subject))) subject <- names(pmwg_mcmc)
      if (!all(subject %in% names(pmwg_mcmc)))
        stop("Subject not present\n")
      for (i in subject) {
        if (plot_acf)
          acf(pmwg_mcmc[[i]],main=paste0("Chain ",acf_chain,": ","LL ",i)) else
            plot(pmwg_mcmc[[i]],auto.layout=auto.layout,density=TRUE,
                 xlab=paste0("Iterations ",i),ask=FALSE,trace=TRUE,
                 ylab=attr(pmwg_mcmc,"selection"),smooth=FALSE,ylim=ylim)
      }
    } else {
      if (plot_acf) for (j in dimnames(pmwg_mcmc)[[2]])
        acf(pmwg_mcmc[,j],main=paste0("Chain ",acf_chain,": ",j)) else
          plot(pmwg_mcmc,auto.layout=auto.layout,density=TRUE,
               ask=FALSE,trace=TRUE,ylim=ylim,
               ylab=attr(pmwg_mcmc,"selection"),smooth=FALSE)
    }
  }
}
