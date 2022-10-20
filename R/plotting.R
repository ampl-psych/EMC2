plot_chains <- function(pmwg_mcmc,layout=NA,subject=NA,ylim=NULL,
                        selection="alpha",filter="sample",thin=1,subfilter=0,
                        plot_acf=FALSE,acf_chain=1, verbose=TRUE) # ,use_par=NA
  # Plots chains  (if alpha, LL or epsilon can do individual subject, all by default)
{
  if (!(class(pmwg_mcmc) %in% c("mcmc","mcmc.list"))) {
    if (class(pmwg_mcmc)=="pmwgs")
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
  if (class(samples)=="pmwgs") samples <- list(samples)
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
