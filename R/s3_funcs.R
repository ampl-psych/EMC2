print.emc <- function(x, ...){
  print(chain_n(x))
}

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
  return(out_list)
}
