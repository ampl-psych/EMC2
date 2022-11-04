plot_chains <- function(pmwg_mcmc,layout=NA,subject=NA,ylim=NULL,
                        selection="alpha",filter="sample",thin=1,subfilter=0,
                        plot_acf=FALSE,acf_chain=1, verbose=TRUE) # ,use_par=NA
  # Plots chains  (if alpha, LL or epsilon can do individual subject, all by default)
{
  if (!(class(pmwg_mcmc) %in% c("mcmc","mcmc.list"))) {
    if (is(pmwg_mcmc, "pmwgs"))
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter) else
                             pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                                       thin=thin,subfilter=subfilter)
  }
  auto.layout <- any(is.na(layout))
  no_layout <- is.null(layout)
  if (!auto.layout & !no_layout) par(mfrow=layout)
  if (attr(pmwg_mcmc,"selection")=="alpha") {
    if (any(is.na(subject))) subject <- names(pmwg_mcmc)
    if (!all(subject %in% names(pmwg_mcmc)))
      stop("Subject not present\n")
    for (i in subject) {
      if (!auto.layout & !no_layout) par(mfrow=layout)
      if (plot_acf) for (j in dimnames(pmwg_mcmc[[i]])[[2]]) {
        acf(pmwg_mcmc[[i]][,j],main=paste0("Chain ",acf_chain,": ","s",i,": ",j))
      } else {
        if (!auto.layout & !no_layout) par(mfrow=layout)
        plot(pmwg_mcmc[[i]],auto.layout=auto.layout,density=FALSE,
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
          plot(pmwg_mcmc[[i]],auto.layout=auto.layout,density=FALSE,
               xlab=paste0("Iterations ",i),ask=FALSE,trace=TRUE,
               ylab=attr(pmwg_mcmc,"selection"),smooth=FALSE,ylim=ylim)
    }
  } else {
    if (plot_acf) for (j in dimnames(pmwg_mcmc)[[2]])
      acf(pmwg_mcmc[,j],main=paste0("Chain ",acf_chain,": ",j)) else
        plot(pmwg_mcmc,auto.layout=auto.layout,density=FALSE,
             ask=FALSE,trace=TRUE,ylim=ylim,
             ylab=attr(pmwg_mcmc,"selection"),smooth=FALSE)
  }
}

plot_acfs <- function(samples,layout=NULL,subject=1,
                      selection="alpha",filter="sample",subfilter=0)
  # Plots acf for all chains
{
  if (is(samples, "pmwgs")) samples <- list(samples)
  if (selection=="alpha") {
    snams <- names(samples[[1]]$data)
    if (is.numeric(subject)) subject <- snams[subject]
    if (!(subject %in% snams)) stop("Subject not present\n")
    message("Plotting chains for subject ",subject)
  }
  for (i in 1:length(samples)) {
    if (selection=="alpha") {
      plot_chains(samples[[i]],selection=selection,filter=filter,subfilter=subfilter,
                  layout=layout,plot_acf=TRUE,acf_chain=i,verbose=FALSE,subject=subject)
    } else plot_chains(samples[[i]],selection=selection,filter=filter,subfilter=subfilter,
                       layout=layout,plot_acf=TRUE,acf_chain=i,verbose=FALSE)
  }
}

check_run <- function(samples,pdf_name="check_run.pdf",interactive=TRUE,
                      filter="sample",subfilter=0,
                      layout=c(3,4),width=NULL,height=NULL) {
  print(chain_n(samples))
  pdf(pdf_name,width=width,height=height)
  plot_chains(samples,selection="LL",layout=layout,filter=filter,subfilter=subfilter)
  if (interactive) readline("Enter for mu check")
  cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("\nR-hat\n")
  print(round(gd_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter,print_summary = FALSE),2))
  cat("\nIntegrated autocorrelation time\n")
  iat_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter)
  cat("\nEffectvie Size\n")
  es_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter,summary==es_stat)
  plot_chains(samples,selection="mu",layout=layout,filter=filter,subfilter=subfilter)
  if (interactive) readline("Enter for variance check")
  cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!! VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("\nR-hat\n")
  print(round(gd_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter,print_summary = FALSE),2))
  cat("\nIntegrated autocorrelation time\n")
  iat_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter)
  cat("\nEffectvie Size\n")
  round(es_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter,summary==es_stat))
  plot_chains(samples,selection="variance",layout=layout,filter=filter,subfilter=subfilter)
  if (interactive) readline("Enter for correlation check")
  cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!! CORRELATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("\nR-hat\n")
  print(round(gd_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter,print_summary = FALSE),2))
  if (interactive) readline("Enter for next correlation check")
  cat("\nIntegrated autocorrelation time\n")
  iat_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter)
  if (interactive) readline("Enter for next correlation check")
  cat("\nEffectvie Size\n")
  round(es_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter,summary==es_stat))
  plot_chains(samples,selection="correlation",layout=layout,ylim=c(-1,1),filter=filter,subfilter=subfilter)
  if (interactive) readline("Enter for alpha check")
  cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALPHA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("\nR-hat\n")
  print(round(gd_pmwg(samples,selection="alpha",filter=filter,subfilter=subfilter,print_summary = FALSE),2))
  cat("\nIntegrated autocorrelation time\n")
  iat_pmwg(samples,selection="alpha",filter=filter,subfilter=subfilter)
  cat("\nEffectvie Size (minimum)\n")
  round(es_pmwg(samples,selection="alpha",summary=min,filter=filter,subfilter=subfilter))
  cat("\nEffectvie Size (mean)\n")
  round(es_pmwg(samples,selection="alpha",summary=mean,filter=filter,subfilter=subfilter))
  plot_chains(samples,selection="alpha",layout=layout,filter=filter,subfilter=subfilter)
  dev.off()
  message("\n\nGraphical checks available in ",pdf_name)
}

