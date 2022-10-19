std_error_IS2 <- function(IS_samples, n_bootstrap = 50000){
  log_marglik_boot= array(dim = n_bootstrap)
  for (i in 1:n_bootstrap){
    log_weight_boot = sample(IS_samples, length(IS_samples), replace = TRUE) #resample with replacement from the lw
    log_marglik_boot[i] <- median(log_weight_boot)
  }
  return(sd(log_marglik_boot))
}

es_pmwg <- function(pmwg_mcmc,selection="alpha",summary_alpha=mean,
                    print_summary=TRUE,sort_print=TRUE,
                    filter="sample",thin=1,subfilter=NULL)
  # Effective size
{
  if (!(class(pmwg_mcmc[[1]]) %in% c("mcmc","mcmc.list"))) {
    if (class(pmwg_mcmc)=="pmwgs")
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter) else
                             pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                                       thin=thin,subfilter=subfilter)
  }
  if (attr(pmwg_mcmc,"selection")=="LL")
    stop("Effective size not sensible for LL\n")
  out <- do.call(rbind,lapply(pmwg_mcmc,effectiveSize))
  if (attr(pmwg_mcmc,"selection")=="alpha") {
    if (!is.null(summary_alpha)) out <- apply(out,2,summary_alpha)
    if (print_summary) if (sort_print) print(round(sort(out))) else
      print(round(out))
    invisible(out)
  } else {
    out <- apply(out,2,sum)
    if (print_summary) if (sort_print) print(round(sort(out))) else
      print(round(out))
    invisible(out)
  }
}

