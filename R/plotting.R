
#' Plot MCMC chains
#'
#' Plots the trace plots of the MCMC chains on top of each other to visualize convergence
#' and chain stability.
#'
#' @param samplers A list of samplers.
#' @param layout A vector specifying the layout as in par(mfrow = layout).
#' If NA (default) will use CODA defaults (unless ``plot_acf = TRUE``), if NULL use current.
#' @param subject Integer or character vector. Only applicable if selection = "alpha". Will plot only these subject(s).
#' NA (default) will plot all.
#' @param ylim A vector. The y limits of the chain plot.
#' @param selection String. Which parameter type to plot ("alpha", "mu", "variance", "covariance", "correlation").
#' "LL" will plot the log-likelihood chains.
#' @param filter A string. Specifies from which stage you want to plot the MCMC chains ("preburn", "burn", "adapt", "sample")                        .
#' @param thin An integer. Will keep only iterations of the MCMC chains that are a multiple of ``thin``.
#' @param subfilter An integer or vector. If integer it will exclude up until
#' that integer. If vector it will include everything in that range.
#' @param plot_acf Boolean. If TRUE will plot autocorrelation
#' for the first MCMC chain (of the default three chains).
#' @return None
#' @examples \dontrun{
#' # For a set of samplers run using run_emc:
#' plot_chains(samplers)
#' # Or for the second subject:
#' plot_chains(samplers, subject = 2)
#'
#' # Can also plot the autocorrelation of for example the group-level mean:
#' plot_chains(samplers, selection = "mu", plot_acf = TRUE)
#' }
#' @export

plot_chains <- function(samplers,layout=NA,subject=NA,ylim=NULL,
                        selection="mu",filter="sample",thin=1,subfilter=0,
                        plot_acf=FALSE,acf_chain=1) # ,use_par=NA
  # Plots chains  (if alpha, LL or epsilon can do individual subject, all by default)
{
  MCMC_samples <- samplers
  if (!(inherits(MCMC_samples,  c("mcmc","mcmc.list")))) {
    if (inherits(MCMC_samples, "pmwgs"))
      MCMC_samples <- as_Mcmc(MCMC_samples,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter) else
                             MCMC_samples <- as_mcmc.list(MCMC_samples,selection=selection,filter=filter,
                                                       thin=thin,subfilter=subfilter)
  }
  auto.layout <- any(is.na(layout))
  no_layout <- is.null(layout)
  if (!auto.layout & !no_layout) par(mfrow=layout)
  if (attr(MCMC_samples,"selection")=="alpha" || attr(MCMC_samples,"selection")=="random") {
    snams <- names(MCMC_samples)
    if (any(is.na(subject))) subject <- snams
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% names(MCMC_samples)))
      stop("Subject not present\n")
    for (i in subject) {
      if (!auto.layout & !no_layout) par(mfrow=layout)
      if (plot_acf){
        for (j in colnames(MCMC_samples[[i]][[1]])){
          acf(MCMC_samples[[i]][[1]][,j],main=paste0("Chain ",acf_chain,": ","s",i,": ",j))
        }
      } else {
        plot(MCMC_samples[[i]],auto.layout=auto.layout,density=FALSE,
             xlab=paste0("Iterations ",i),ask=FALSE,trace=TRUE,
             ylab=attr(MCMC_samples,"selection"),smooth=FALSE,ylim=ylim)
      }
    }
  } else if (attr(MCMC_samples,"selection") %in% c("LL","epsilon")) {
    if (any(is.na(subject))) subject <- names(MCMC_samples)
    if (!all(subject %in% names(MCMC_samples)))
      stop("Subject not present\n")
    for (i in subject) {
      if (plot_acf){
        acf(MCMC_samples[[i]][[1]],main=paste0("Chain ",acf_chain,": ","LL ",i))
      } else{
        plot(MCMC_samples[[i]],auto.layout=auto.layout,density=FALSE,
             xlab=paste0("Iterations ",i),ask=FALSE,trace=TRUE,
             ylab=attr(MCMC_samples,"selection"),smooth=FALSE,ylim=ylim)
      }
    }
  } else {
    if (plot_acf) {
      for (j in colnames(MCMC_samples[[1]])){
        acf(MCMC_samples[[1]][,j],main=paste0("Chain ",acf_chain,": ",j))
      }
    } else {
      plot(MCMC_samples,auto.layout=auto.layout,density=FALSE,
           ask=FALSE,trace=TRUE,ylim=ylim,
           ylab=attr(MCMC_samples,"selection"),smooth=FALSE)
    }
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
    if (!all(subject %in% snams)) stop("Subject not present\n")
    message("Plotting chains for subject(s) ",paste(subject,collapse=" "))
  }
  for (i in 1:length(samples)) {
    if (selection=="alpha") {
      plot_chains(samples[[i]],selection=selection,filter=filter,subfilter=subfilter,
                  layout=layout,plot_acf=TRUE,acf_chain=i,subject=subject)
    } else plot_chains(samples[[i]],selection=selection,filter=filter,subfilter=subfilter,
                       layout=layout,plot_acf=TRUE,acf_chain=i)
  }
}

#' plot_alpha_recovery
#' Uses output from plot_pars to plot true vs. estimated (median with CI) alpha
#' parameters for each subject (plot density must be called with true values passed
#' through the pars argument).
#'
#' @param tabs Tables of actual and estimated alpha parameters (with CIs) from plot_pars
#' @param layout A 2-vector specifying the layout as in par(mfrow = layout)
#' @param do_ci Boolean (Default TRUE). Add CIs to plot?
#' @param ci_col Color of CI.
#' @param cap Width of CI cap (passed to arrows)
#' @param do_rmse Boolean (default FALSE) Add root-mean squared error to plot
#' instead of default Pearson correlation
#' @param r_pos Position of Pearson/RMSE (passed to legend)
#' @param rmse_digits Digits for RMSE
#' @param pearson_digits Digits for Pearson correlation
#' @param do_coverage Boolean (default TRUE) add coverage percentage estimate (otherwise
#' add Spearman correlation)
#' @param coverage_pos Position of Coverage/Pearson (passed to legend)
#' @param coverage_digits Digits for coverage
#' @param spearman_digits Digits for Spearman correlation
#'
#' @return Invisible list with RMSE, coverage and Pearson and Spearman correlations.
#' @export

plot_alpha_recovery <- function(tabs,layout=c(2,3),
                                do_ci = TRUE,ci_col="grey",cap=.05,
                                do_rmse=FALSE,r_pos="topleft",
                                rmse_digits=3,pearson_digits=2,
                                do_coverage=TRUE,coverage_pos="bottomright",
                                coverage_digits=1,spearman_digits=2)
  # Takes tables output by plot_pars with par and plots recovery
{
  par(mfrow=layout)
  pnams <- dimnames(tabs[[1]])[[2]]
  if (!do_rmse) {
    rmse <- NULL
    pearson <- setNames(numeric(length(pnams)),pnams)
  } else {
    pearson <- NULL
    rmse <- setNames(numeric(length(pnams)),pnams)
  }
  if (do_coverage) {
    coverage <- NULL
    spearman <- setNames(numeric(length(pnams)),pnams)
  } else {
    spearman=NULL
    coverage <- setNames(numeric(length(pnams)),pnams)
  }
  for (p in pnams) {
    xy <- do.call(rbind,lapply(tabs,function(x){x[,p]}))
    ylim <- c(min(xy),max(xy))
    plot(xy[,"true"],xy[,"50%"],ylim=ylim,xlim=ylim,main=p,
         xlab="Generating",ylab="Estimated")
    abline(a=0,b=1,lty=3)
    if (do_ci) for (i in 1:dim(xy)[1]) {
      arrows(xy[i,"true"],xy[i,"50%"],xy[i,"true"],xy[i,"97.5%"],
             col=ci_col,angle=90,length=cap)
      arrows(xy[i,"true"],xy[i,"50%"],xy[i,"true"],xy[i,"2.5%"],
             col=ci_col,angle=90,length=cap)
    }
    if (do_rmse) {
      rmse[p] <- sqrt(mean((xy[,"true"] - xy[,"50%"])^2))
      legend(r_pos,paste("RMSE = ",round(rmse[p],rmse_digits)),bty="n")
    } else {
      pearson[p] <- cor(xy[,"true"],xy[,"50%"],method="pearson")
      legend(r_pos,paste("r(pearson) = ",round(pearson[p],pearson_digits)),bty="n")
    }
    if (do_coverage) {
      coverage[p] = 100*mean((xy[,"97.5%"] > xy[,"true"])  & (xy[,"2.5%"] < xy[,"true"]))
      legend(coverage_pos,paste("95% Coverage = ",round(coverage[p],coverage_digits)),bty="n")
    } else {
      spearman[p] <- cor(xy[,"true"],xy[,"50%"],method="spearman")
      legend(coverage_pos,paste("r(spearman) = ",round(spearman[p],spearman_digits)),bty="n")
    }
  }
  invisible(list(RMSE = rmse,COVERAGE = coverage,PEARSON=pearson,SPEARMAN=spearman))
}

#' Plot defective densities for each subject and cell.
#'
#' Plots multiple panels that contain a set of densities for each response option in the data.
#' These densities are defective; their areas are relative to their response's proportion.
#' Across all responses the area sums to 1.
#'
#' @param data A dataframe. The experimental data in EMC format with at least the subjects factor,
#'  R (response factor) and rt (response time) columns,
#' and optionally other factor columns of the design.
#' @param subject Integer or string selecting a subject from the data. If specified will only plot that subject
#' (default NULL = all).
#' @param factors character vector of factor names in design to aggregated across (default NULL = all).
#' @param layout vector specifying plot window layout: par(mfrow = layout), default NULL uses current plot window layout.
#' @param xlim x-axis limit for all cells (default NULL = scale per cell).
#' @param bw number or string bandwidth for density (default "nrd0"). See ``?density``.
#' @param adjust Numeric. Density function bandwidth adjust parameter. See ``?density``.
#' @param correct_fun If specified will calculate the accuracy for each subject using the supplied function and
#' invisibly returns an accuracy vector for each subject.
#' @param rt legend function position string for mean RT (default "top)
#' @param accuracy legend function position string for accuracy (default "topright")
#'
#' @return If correct_fun is specified, will invisibly return a subject accuracy vector
#' @examples
#' # First for each subject and the factor combination in the design:
#' plot_defective_density(forstmann)
#' # Now collapsing across subjects:
#' plot_defective_density(forstmann, factors = c("S", "E"))
#' # If your data is response coded it always makes sense to include the "S" factor
#' # because EMC will plot the "R" factor automatically. This way you can see how often
#' # "S" matches "R".
#' # We can also return each subject's accuracy using a custom function:
#' print(plot_defective_density(forstmann, correct_fun = function(d) d$R == d$S))
#'
#' @export

plot_defective_density <- function(data,subject=NULL,factors=NULL,
                                   layout=NULL,
                                   xlim=NULL,bw = "nrd0",adjust=1,
                                   correct_fun=NULL,rt="top",accuracy="topright")
{
  if (!is.null(subject)) {
    snams <- levels(data$subjects)
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% snams)) stop("Subject not present\n")
    dat <- data[data$subjects==subject,]
    fnams <- names(dat)[!(names(dat) %in% c("subjects","trials","R","rt"))]
  } else {
    dat <- data
    fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))]
  }
  if (!is.null(factors)) {
    if (!all(factors %in% fnams))
      stop("factors must name factors in data")
    fnams <- factors
  }
  alldat <- dat
  cellsall <- alldat[,fnams,drop=FALSE]
  for (i in fnams) cellsall[,i] <- paste(i,cellsall[,i],sep="=")
  cellsall <- apply(cellsall,1,paste,collapse=" ")
  # Remove missing
  dat <- dat[is.finite(dat$rt),]
  cells <- dat[,fnams,drop=FALSE]
  for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
  cells <- apply(cells,1,paste,collapse=" ")
  if (!is.null(layout))
    par(mfrow=layout)
  R <- levels(dat$R)
  for (i in sort(unique(cells))) {
    pR <- table(alldat$R[cellsall==i])/dim(alldat[cellsall==i,])[1]
    dati <- dat[cells==i,]
    mrt <- tapply(dati$rt,dati$R,median)
    dens <- setNames(vector(mode="list",length=length(R)),R)
    for (j in R) if (length(dati$rt[dati$R==j])>1) {
      dens[[j]] <- density(dati$rt[dati$R==j],bw=bw,adjust=adjust)
      dens[[j]]$y <- dens[[j]]$y*pR[j]
    } else dens[[j]] <- NULL
    rx <- do.call(rbind,lapply(dens,function(x){range(x$x)}))
    if (is.null(xlim)) xlimi <- c(min(rx[,1]),max(rx[,2])) else xlimi <- xlim
    ylim <- c(0,max(unlist(lapply(dens,function(x){max(x$y)}))))
    ltys <- c(1:length(mrt))[!is.na(mrt)]
    plot(dens[[1]],xlim=xlimi,ylim=ylim,lty=ltys[1],main=i,xlab="RT")
    if (length(dens)>1) for (j in 2:length(dens))
      lines(dens[[j]],lty=ltys[j])
    if (!is.null(accuracy) && length(R) > 1) {
      ltys <- 1:length(pR)
      ltys[is.na(mrt)] <- NA
      legend(accuracy,paste(names(pR),round(pR,2),sep="="),lty=ltys,bty="n",title="p(R)")
    }
    if (!is.null(rt))
      legend(rt,paste(names(pR),round(mrt,3),sep="="),bty="n",title="Med(RT)")
  }
  if (!is.null(correct_fun))
    invisible(tapply(correct_fun(dat),dat$subjects,mean))
}



#' Plots density for parameters.
#'
#' Plots the posterior and prior density for selected parameters of a model.
#'
#' @param samplers A list of samplers.
#' @param layout A vector specifying the layout as in par(mfrow = layout).
#' If NA or NULL use current plot layout.
#' @param selection A string. Which parameter type to plot ("alpha", "mu", "variance", "covariance", "correlation").
#' @param use_par Character vector of names of parameters to plot (default NULL = plot all)
#' @param filter A string. Specifies from which stage you want to plot the densities ("preburn", "burn", "adapt", "sample")
#' @param thin An integer. Will keep only iterations of the MCMC chains that are a multiple of ``thin``.
#' @param subfilter An integer or vector. If integer it will exclude up until
#' that integer. If vector it will include everything in that range.
#' @param mapped Boolean. If TRUE plots the parameters mapped back to experimental design
#' otherwise plots the sampled parameters.
#' @param plot_prior Boolean. Overlay prior density in the plot (in red)
#' @param xlim x-axis plot limit. If a vector is supplied will use the same axes for all.
#  Alternatively a matrix can be supplied with one row for each parameter.
#' @param ylim y-axis plot limit. If a vector is supplied will use the same axes for all.
#  Alternatively a matrix can be supplied with one row for each parameter.
#' @param prior_xlim A vector giving upper and lower quantiles of prior when choosing
#' xlim if ``plot_prior = TRUE``. If set to NULL xlim is used instead.
#' @param show_chains Boolean (default FALSE) plot separate density for each chain.
#' @param do_plot Boolean. Set to ``FALSE`` to only return the parameter credible intervals and omit the plots.
#' @param subject Integer or character vector. Only applicable if ``selection = "alpha"``. Will plot only these subject(s).
#' NA (default) will plot all.
#' @param add_means Boolean. Whether to add parameter means as an attribute
#' to return invisibly
#' @param pars Named vector or matrix of known/simulated parameters, or the output of
#' plot_pars, in which case the posterior medians are extracted. If xlim is not supplied,
#  the plot will be adjusted to include the parameter values, which are plotted as vertical lines.
#' @param probs Vector. The quantiles of the selected parameters to return invisibly.
#' @param bw number or string bandwidth for density (default "nrd0"). See ``?density``.
#' @param adjust Numeric. Adjustment for density plot. See ``?density``
#' @param lpos Character. Position of the contraction in the plot
#' @param digits Integer. Rounding of contraction
#'
#' @return invisibly return quantiles for the selected parameters
#' @examples \dontrun{
#' # For a set of samplers plot prior and posterior densities:
#' plot_pars(samplers)
#' # Or just get the quantiles and omit the plot
#' quants <- plot_pars(samplers, selection = "variance", do_plot = FALSE)
#' print(quants)
#' }
#' @export
plot_pars <- function(samplers,layout=c(2,3),use_par=NULL,
                      selection="mu",filter="sample",thin=1,subfilter=0,mapped=FALSE,
                      plot_prior=TRUE,xlim=NULL,ylim=NULL,prior_xlim=c(.1,.9),
                      show_chains=FALSE,do_plot=TRUE,subject=NA,add_means=FALSE,
                      pars=NULL,probs=c(.025,.5,.975),bw = "nrd0", adjust = 1,
                      lpos="topright",digits=3)

{

  robust_density <- function(ps,r,bw,adjust,use_robust=FALSE)
    # density estimate over range r (to avoid outlier influence, used
    # here for hyper co/variance which can have long tails)
  {
    if (use_robust) {
      # isin <- ps>r[1]*.9 & ps < r[2]*1.1
      isin <- ps> 0 & ps < r[2]*1.1
      p <- mean(isin)
      psc <- ps[isin]
      if (length(psc)<10) {
        dens <- density(ps)
        dens$y[1:length(dens$y)] <- 0
      } else {
        dens <- density(psc,bw=bw,adjust=adjust)
        dens$y <- p*dens$y
      }
      dens
    } else density(ps,bw=bw,adjust=adjust)
  }
  n_prior <- 1e3
  pmwg_mcmc <- samplers
  if (mapped & !(selection %in% c("mu","alpha")))
    stop("mapped only available for mu and alpha")

  if (show_chains & plot_prior)
    warning("Prior plots not implemented for show_chains=TRUE")
  do_contraction <- TRUE
  if (do_contraction & !plot_prior) do_contraction <- FALSE

  if (is.list(pars)) pars <- do.call(rbind,lapply(pars,function(x)x[2,]))

  if (!(inherits(pmwg_mcmc, c("mcmc","mcmc.list")))) {
    if (plot_prior) {
      psamples <- get_prior_samples(pmwg_mcmc,selection,filter,thin,subfilter,n_prior)
      if (mapped) psamples <- map_mcmc(psamples,design=attr(pmwg_mcmc,"design_list")[[1]],
                                       model=attr(pmwg_mcmc,"model_list")[[1]])
      if (!is.null(attr(psamples,"isConstant")))
        psamples <- psamples[,!attr(psamples,"isConstant"),drop=FALSE]
    }
    if (inherits(pmwg_mcmc, "pmwgs")){
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc, selection=selection,filter=filter,thin=thin,subfilter=subfilter)
    } else {
      if (mapped & !is.null(pars)) {
        pars <- map_mcmc(pars,design=attr(pmwg_mcmc,"design_list")[[1]],
                         model=attr(pmwg_mcmc,"model_list")[[1]])
        if (!is.null(attr(pars,"isConstant")))
          pars <- pars[,!attr(pars,"isConstant"),drop=FALSE]
      }
      pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                thin=thin,subfilter=subfilter,mapped=mapped)
    }
  } else plot_prior <- FALSE
  if (attr(pmwg_mcmc,"selection")=="LL")
    stop("No density plots for LL\n")
  no_layout <- any(is.na(layout)) | is.null(layout)
  chains <- 0


  if (attr(pmwg_mcmc,"selection") == "alpha") {
    snams <- names(pmwg_mcmc)
    if (any(is.na(subject))) subject <- snams
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% snams))
      stop("Subject not present\n")
    if (!is.null(pars)) {
      if (!is.null(dim(pars))) pars <- t(pars) else {
        if (length(subject)!=1) reps <- length(subject) else reps <- 1
        pars <- matrix(rep(pars,each=reps),ncol=1,dimnames=list(names(pars),subject) )
      }
    }
    if (!is.null(pars) && !is.matrix(pars))
      stop("pars must be a matrix for alpha")
    if (inherits(pmwg_mcmc[[1]], "mcmc.list")) {
      selection <- attr(pmwg_mcmc,"selection")
      pmwg_mcmc_combined <- lapply(pmwg_mcmc,function(x){do.call(rbind,x)})
      attr(pmwg_mcmc_combined,"selection") <- selection
      if (show_chains) chains <- length(pmwg_mcmc[[1]]) else
        pmwg_mcmc <- pmwg_mcmc_combined
    } else pmwg_mcmc_combined <- pmwg_mcmc
    tabs <- lapply(pmwg_mcmc_combined,function(x){apply(x,2,quantile,probs=probs)})
    # if (plot_prior) dimnames(psamples) <- list(NULL,colnames(pmwg_mcmc_combined[[1]]))
    if (do_contraction)
      contraction <- setNames(vector(mode="list",length=length(subject)),subject)
    if (!is.null(use_par)) {
      ok <- colnames(pmwg_mcmc_combined[[1]]) %in% use_par
      if (!any(ok)) stop("use_par did not specify parameters that are present")
      tabs <- lapply(tabs,function(x) x[,ok,drop=FALSE])
    } else ok <- rep(TRUE,length(colnames(pmwg_mcmc_combined[[1]])))
    for (i in subject) {
      if (do_plot) {
        if (!no_layout) par(mfrow=layout)
        if (do_contraction)
          contractioni <- setNames(numeric(length(colnames(pmwg_mcmc_combined[[i]]))),
                                   colnames(pmwg_mcmc_combined[[i]]))
        for (j in  colnames(pmwg_mcmc_combined[[i]])[ok] ) {
          if (chains>0) {
            dens <- lapply(pmwg_mcmc[[i]],function(x){density(x[,j],bw=bw,adjust=adjust)})
            if (!is.null(xlim)) {
              if (!is.matrix(xlim)) xlimi <- xlim else xlimi <- xlim[j,]
            } else {
              xlimi <- c(min(unlist(lapply(dens,function(x){min(x$x)}))),
                         max(unlist(lapply(dens,function(x){max(x$x)}))))
              if (!is.null(pars)) xlimi <- c(min(c(pars[j,i]-abs(pars[j,i])/10,xlimi[1])),
                                             max(c(pars[j,i]+abs(pars[j,i])/10,xlimi[2])))
            }
            if (!is.null(ylim)) {
              if (!is.matrix(ylim)) ylimi <- ylim else ylimi <- ylim[j,]
            } else ylimi <- c(0,max(unlist(lapply(dens,function(x){max(x$y)}))))
            main <- paste0(attr(pmwg_mcmc,"selection")," s",i)
            plot(dens[[1]],xlab=j,main=main,col=1,xlim=xlimi,ylim=ylimi)
            if (chains>1) for (k in 2:chains) lines(dens[[k]],col=k)
          } else {
            dens <- density(pmwg_mcmc[[i]][,j],bw=bw,adjust=adjust)
            if (plot_prior) pdens <- density(psamples[,j],bw=bw,adjust=adjust)
            if (!is.null(xlim)) {
              if (!is.matrix(xlim)) xlimi <- xlim else xlimi <- xlim[j,]
            } else if (plot_prior & !is.null(prior_xlim)) {
              xlimi <- c(quantile(pdens$x,probs=prior_xlim[1]),
                         quantile(pdens$x,probs=prior_xlim[2]))
            } else {
              xlimi <- c(min(dens$x),max(dens$x))
              if (!is.null(pars)) xlimi <- c(min(c(pars[j,i]-abs(pars[j,i])/10,xlimi[1])),
                                             max(c(pars[j,i]+abs(pars[j,i])/10,xlimi[2])))
            }
            if (!is.null(ylim)) {
              if (!is.matrix(ylim)) ylimi <- ylim else ylimi <- ylim[j,]
            } else {
              ylimi <- c(0,max(dens$y))
              if (plot_prior) ylimi[2] <- max(c(ylimi[2],pdens$y))
            }
            main <- paste0(attr(pmwg_mcmc,"selection")," s",i)
            plot(dens,xlab=j,xlim=xlimi,ylim=ylimi,main=main)
            if (plot_prior) lines(pdens,col="red")
            if (do_contraction) {
              contractioni[j] <- 1-(var(pmwg_mcmc[[i]][,j])/var(psamples[,j]))
              legend(lpos,legend=round(contractioni[j],digits),bty="n",title="Contraction")
            }
          }
          if (!is.null(pars)) abline(v=pars[j,i])
        }
      }
      if (do_plot & do_contraction) contraction[[i]] <- contractioni
      if (!is.null(pars)) {
        tabs[[i]] <- rbind(true=pars[dimnames(tabs[[i]])[[2]],i],tabs[[i]])
        tabs[[i]] <- rbind(tabs[[i]],Miss=tabs[[i]][3,]-tabs[[i]][1,])
      }
    }
    tabs <- tabs[as.character(subject)]
    if (add_means)
      attr(tabs,"mean") <- lapply(pmwg_mcmc_combined,function(x){apply(x,2,mean)})
    if (do_contraction & exists("contraction"))
      attr(tabs,"contraction") <- do.call(rbind,contraction)
    invisible(tabs)
  } else {
    if (!is.na(subject)) warning("Selecting a subject has no effect on population parameters")
    if (inherits(pmwg_mcmc, "mcmc.list")) {
      selection <- attr(pmwg_mcmc,"selection")
      pmwg_mcmc_combined <- do.call(rbind,pmwg_mcmc)
      attr(pmwg_mcmc_combined,"selection") <- selection
      if (show_chains) chains <- length(pmwg_mcmc) else
        pmwg_mcmc <- pmwg_mcmc_combined
    } else pmwg_mcmc_combined <- pmwg_mcmc
    if(!is.null(use_par)){
      ok <- colnames(pmwg_mcmc_combined) %in% use_par
      if (!any(ok)) stop("use_par did not specify parameters that are present")
    }

    if (!is.null(pars)) {
      if ( !is.vector(pars) ) {
        if (attr(pmwg_mcmc,"selection") == "mu")
          stop("pars must be a vector for mu")
        if (!is.matrix(pars))
          stop("pars must be a vector or matrix")
        if (attr(pmwg_mcmc,"selection") == "variance")
          pars <- diag(pars) else
            pars <- pars[lower.tri(pars)]
      }
      if(is.null(use_par)){
        if (length(pars) != dim(pmwg_mcmc_combined)[2])
          stop("pars is wrong length")
        names(pars) <- colnames(pmwg_mcmc_combined)
      } else{
        if (length(pars) != length(use_par))
          stop("pars is wrong length")
        names(pars) <- use_par
      }
    }
    tabs <- rbind(true=pars,apply(pmwg_mcmc_combined,2,quantile,probs=probs))
    if (!is.null(pars)) tabs <- rbind(tabs,Miss=tabs[3,]-tabs[1,])
    if (add_means) attr(tabs,"mean") <- apply(pmwg_mcmc_combined,2,mean)
    if (!no_layout) par(mfrow=layout)
    if (plot_prior) dimnames(psamples) <- list(NULL,colnames(pmwg_mcmc_combined))
    if (do_contraction)
      contraction <- setNames(numeric(length(colnames(pmwg_mcmc_combined))),colnames(pmwg_mcmc_combined))
    if (!is.null(use_par)) {
      tabs <- tabs[,ok,drop=FALSE]
    } else ok <- rep(TRUE,length(colnames(pmwg_mcmc_combined)))
    if (do_plot) for (j in colnames(pmwg_mcmc_combined)[ok] ) {
      if (chains > 0) {
        dens <- lapply(pmwg_mcmc,function(x){density(x[,j],bw=bw,adjust=adjust)})
        if (!is.null(xlim)) xlimi <- xlim else {
          xlimi <- c(min(unlist(lapply(dens,function(x){min(x$x)}))),
                     max(unlist(lapply(dens,function(x){max(x$x)}))))
          if (!is.null(pars)) xlimi <- c(min(c(pars[j,i]-abs(pars[j,i])/10,xlimi[1])),
                                         max(c(pars[j,i]+abs(pars[j,i])/10,xlimi[2])))
        }
        if (!is.null(ylim)) ylimi <- ylim else
          ylimi <- c(0,max(unlist(lapply(dens,function(x){max(x$y)}))))
        main <- j
        plot(dens[[1]],xlab=attr(pmwg_mcmc,"selection"),main=main,
             col=1,xlim=xlimi,ylim=ylimi)
        if (chains>1) for (k in 2:chains) lines(dens[[k]],col=k)
      } else {
        dens <- density(pmwg_mcmc[,j],bw=bw,adjust=adjust)
        if (plot_prior)  pdens <- robust_density(psamples[,j],range(pmwg_mcmc[,j]),
                                                 bw=bw,adjust=adjust,use_robust=!(attr(pmwg_mcmc,"selection") %in% c("mu","correlation")))
        if (!is.null(xlim)) xlimi <- xlim else if (plot_prior & !is.null(prior_xlim))
          xlimi <- c(quantile(pdens$x,probs=prior_xlim[1]),
                     quantile(pdens$x,probs=prior_xlim[2])) else {
                       xlimi <- c(min(dens$x),max(dens$x))
                       if (!is.null(pars)) xlimi <- c(min(c(pars[j,i]-abs(pars[j,i])/10,xlimi[1])),
                                                      max(c(pars[j,i]+abs(pars[j,i])/10,xlimi[2])))
                     }
        if (!is.null(ylim)) ylimi <- ylim else {
          ylimi <- c(0,max(dens$y))
          if (plot_prior) ylimi[2] <- max(c(ylimi[2],pdens$y))
        }
        main <- j
        plot(dens,xlim=xlimi,ylim=ylimi,xlab=attr(pmwg_mcmc,"selection"),main=main)
        if (plot_prior) {
          lines(pdens,col="red")
          if (do_contraction) {
            contraction[j] <- 1-(var(pmwg_mcmc[,j])/var(psamples[,j]))
            legend(lpos,legend=round(contraction[j],digits),bty="n",title="Contraction")
          }
        }
      }
      if (!is.null(pars)) abline(v=pars[j])
    }
    if (do_contraction & exists("contraction")) attr(tabs,"contraction") <- contraction
    invisible(tabs)
  }
}


plot_roc <- function(data,signalFactor="S",zROC=FALSE,qfun=NULL,main="",lim=NULL)

{
  tab <- table(data$R,data[[signalFactor]])
  ctab <- 1-apply(t(tab)/apply(tab,2,sum),1,cumsum)[-dim(tab)[1],] # p(Signal)
  if (!zROC) {
    if (!is.null(lim)) {xlim <- lim; ylim <- lim} else
    {xlim <- c(0,1); ylim <- c(0,1)}
    plot(ctab[,1],ctab[,2],xlab="p(FA)",ylab="p(H)",xlim=xlim,ylim=ylim,main=main)
    lines(ctab[,1],ctab[,2])
    abline(a=0,b=1,lty=3)
  } else {
    ctab <- qfun(ctab)
    ctab <- ctab[apply(ctab,1,function(x){all(is.finite(x))}),]
    if (is.null(lim)) lim <- c(min(ctab),max(ctab))
    plot(ctab[,1],ctab[,2],main=main,
         xlab="z(FA)",ylab="z(H)",xlim=lim,ylim=lim)
    lines(ctab[,1],ctab[,2])
    abline(a=0,b=1,lty=3)
  }
  invisible(ctab)
}

#' Plots overlaying observed and fitted data.
#'
#' If rt  available plots defective cumulative distributions
#' functions (CDFs), if not plots region under curve (ROC).
#' CDFs plot the probability of a response, p(R) as a function of response time
#' (RT) for data (black lines and points at qpoints quantiles) and posterior-
#' predictive simulations (grey lines and points). Large grey points show the
#' average of replicate quantiles and small grey points percentiles for
#' individual replicates, providing a representation of uncertainty in the model
#' predictions.
#'
#' If the stat argument (which calculates
#' a statistic based on the data) is supplied fit is plotted as a density with
#' a vertical line at the position of the data statistic. If more than one
#' subject is included data and fits are aggregated over subjects. If data or
#' fit contains NA in responses or rt, or is.infinite(rt) these are treated as
#' missing and defective cdfs sum to the probability of non-missing.
#'

#'
#' @param data A dataframe. The experimental data in EMC format with at least the subjects factor,
#' R (response factor) and rt (response time) columns,
#' and optionally other factor columns of the design.
#' @param pp Posterior predictives created by ``post_predict``
#' @param subject Integer or string selecting a subject from the data. If specified will only plot that subject
#' (default NULL = all).
#' @param factors Character vector of factors in data to display separately. If
#' NULL (default) use names of all columns in data except "trials","R", and "rt".
#' Omitted factors are aggregated over. If NA treats entire data set as a single cell.
#' Must be NA or NULL when using stat argument.
#' @param functions A named list of functions that create new factors which can then be
#' used by the factors and stat arguments.
#' @param stat A function that takes a the data and returns a single value.
#' @param stat_name A string naming what the stat argument calculates.
#' @param quants A vector. Quantiles to return when
#' stat argument is supplied.
#' @param do_plot Boolean. Set to ``FALSE`` to only return the quantiles and omit the plots.
#' @param xlim x-axis plot limit. If a vector is supplied will use the same axes for all.
#  Alternatively a matrix can be supplied with one row for each parameter.
#' @param ylim y-axis plot limit. If a vector is supplied will use the same axes for all.
#  Alternatively a matrix can be supplied with one row for each parameter.
#' @param layout A vector specifying the layout as in par(mfrow = layout).
#' If NA or NULL use current.
#' @param probs Vector of probabilities at which to calculate cdf
#' @param data_lwd Integer. Line width for data in cdf
#' @param fit_lwd Integer. Line width for posterior predictives in cdf
#' @param qp_cex Numeric. Cex for data quantile points in cdf
#' @param q_points Vector. Quantile points to plot in cdf
#' @param pqp_cex Numeric. Cex for predicted quantile points in cdf
#' @param lpos Character. Legend position, see ``?legend``.
#' @param signalFactor The name of the "signal" factor in an ROC.
#' @param zROC Boolean. Wheter to plot ROC on transformed scale.
#' @param qfun A function. Scale transform function for zROC (default qnorm)
#' @param lim Vector. Limits for x and y limit of ROC
#' @param rocfit_cex Numeric. Cex for predicted ROC points
#' @param adjust Numeric. Density function bandwidth adjust parameter. See ``?density`
#' @param main Text title, pasted before cell name.
#'
#' @return If stat argument is provided a table of observed values and predicted quantiles
#' @export
plot_fit <- function(data,pp,subject=NULL,factors=NULL,functions=NULL,
                     stat=NULL,stat_name="",adjust=1,
                     quants=c(.025,.5,.975),do_plot=TRUE,
                     xlim=NULL,ylim=NULL,main="",
                     layout=NULL,
                     probs=c(1:99)/100,
                     data_lwd=2,fit_lwd=1,qp_cex=1,
                     q_points=c(.1,.3,.5,.7,.9),pqp_cex=.5,lpos="topleft",
                     signalFactor="S",zROC=FALSE,qfun=qnorm,lim=NULL,rocfit_cex=.5)
{
  if (!is.null(stat) & (!all(is.na(factors))|is.null(factors))) {
    if (is.null(factors)) factors <- NA else {
      warning("factors must be NA or NULL when using stat, set to NA")
      factors <- NA
    }
  }
  if (!is.null(subject)) {
    snams <- levels(data$subjects)
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% snams)) stop("Subject(s) not present\n")
    dat <- as.data.frame(droplevels(data[data$subjects %in% subject,]))
    pp <- droplevels(pp[pp$subjects %in% subject,])
    if (length(subject>1))
      fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))] else
        fnams <- names(dat)[!(names(dat) %in% c("subjects","trials","R","rt"))]
  } else {
    dat <- as.data.frame(data)
    fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))]
  }

  okd <- !is.na(dat$R) & is.finite(dat$rt)
  okpp <- !is.na(pp$R) & is.finite(pp$rt)

  if (!is.null(functions)) for (i in 1:length(functions)) {
    dat <- cbind.data.frame(functions[[i]](dat),dat)
    names(dat)[1] <- names(functions)[i]
    pp <- cbind.data.frame(functions[[i]](pp),pp)
    names(pp)[1] <- names(functions)[i]
    fnams <- c(names(functions)[i],fnams)
  }

  if (!is.null(factors)) {
    if (any(is.na(factors))) fnams <- NA else {
      if (!all(factors %in% fnams))
        stop("factors must name factors in data")
      fnams <- factors
    }
  }
  if (!any(is.na(layout))) if (!is.null(layout))
      par(mfrow=layout)
  if (all(is.na(data$rt))) {  # type=SDT
    if (length(levels(data$R))==2 & is.null(stat))
      stop("No plots for binary responses, use an accuracy function in stat arguement.")
    if (!is.null(stat)) { # statistic
      if (!any(is.na(fnams))) {
        cells <- dat[,fnams,drop=FALSE]
        for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
        cells <- apply(cells,1,paste,collapse=" ")
        pp_cells <- pp[,fnams,drop=FALSE]
        for (i in fnams) pp_cells[,i] <- paste(i,pp_cells[,i],sep="=")
        pp_cells <- apply(pp_cells,1,paste,collapse=" ")
        postn <- unique(pp$postn)
        ucells <- sort(unique(cells))
      } else ucells <- ""
      tab <- matrix(nrow=length(ucells),ncol=4,
                    dimnames=list(ucells,c("Observed",names(quantile(1:5,quants)))))
      for (i in ucells) {
        if (i=="") {
          obs <- stat(dat)
          ppi <- pp
          pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
          tab[1,] <- c(obs,quantile(pred,quants))
        } else {
          obs <- stat(dat[cells==i,])
          ppi <- pp[pp_cells==i,]
          pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
          tab[i,] <- c(obs,quantile(pred,quants))
        }
        if (do_plot) {
          dens <- density(pred,adjust=adjust)
          if (!is.null(xlim)) xlimi <- xlim else
            xlimi <- c(pmin(obs,min(dens$x)),pmax(obs,max(dens$x)))
          plot(dens,main=paste(main,i),xlab=stat_name,xlim=xlimi)
          abline(v=obs)
        }

      }
      invisible(tab)
    } else {
      if (!any(fnams==signalFactor))
        stop("Data does not have a column specified in the signalFactor argument: ",signalFactor)
      if (length(levels(data[[signalFactor]]))!=2)
        stop("signalFactor must have exactly two levels for an ROC plot")
      if (zROC & is.null(qfun))
        stop("Must supply qfun for zROC")
      fnams <- fnams[fnams != signalFactor]
      if (!any(is.na(fnams))) {
        cells <- dat[,fnams,drop=FALSE]
        for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
        cells <- apply(cells,1,paste,collapse=" ")
        pp_cells <- pp[,fnams,drop=FALSE]
        for (i in fnams) pp_cells[,i] <- paste(i,pp_cells[,i],sep="=")
        pp_cells <- apply(pp_cells,1,paste,collapse=" ")
        ucells <- sort(unique(cells))
      } else ucells <- ""
      postn <- unique(pp$postn)
      for (i in ucells) {
        if (i=="") {
          dati <- dat
          ppi <- pp
        } else {
          dati <- dat[cells==i,]
          ppi <- pp[pp_cells==i,]
        }
        dpts <- plot_roc(dati,zROC=zROC,qfun=qfun,lim=lim,main=paste(main,i),signalFactor=signalFactor)
        tab <- table(ppi$postn,ppi$R,ppi[[signalFactor]])
        ctab <- apply(tab,1,function(x){list(1-apply(t(x)/apply(x,2,sum),1,cumsum)[-dim(x)[1],])})
        if (!zROC) lapply(ctab,function(x){
          points(x[[1]][,1],x[[1]][,2],col="grey",pch=16,cex=rocfit_cex)
        }) else ctab <- lapply(ctab,function(x){
          x[[1]] <- qnorm(x[[1]])
          points(x[[1]][row.names(dpts),1],x[[1]][row.names(dpts),2],col="grey",pch=16,cex=rocfit_cex)
        })
        points(dpts[,1],dpts[,2])
        lines(dpts[,1],dpts[,2])
      }
    }
  } else {
    if (!any(is.na(fnams))) {
      cells <- dat[,fnams,drop=FALSE]
      for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
      cells <- apply(cells,1,paste,collapse=" ")
      pp_cells <- pp[,fnams,drop=FALSE]
      for (i in fnams) pp_cells[,i] <- paste(i,pp_cells[,i],sep="=")
      pp_cells <- apply(pp_cells,1,paste,collapse=" ")
    }
    if (!is.null(stat)) { # statistic
      postn <- unique(pp$postn)
      if (any(is.na(fnams))) ucells <- "" else ucells <- sort(unique(cells))
      tab <- matrix(nrow=length(ucells),ncol=4,
                    dimnames=list(ucells,c("Observed",names(quantile(1:5,quants)))))
      for (i in ucells) {
        if (i=="") {
          dati <- dat
          ppi <- pp
          obs <- stat(dati)
          pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
          tab[1,] <- c(obs,quantile(pred,quants))
        } else {
          dati <- dat[cells==i,]
          ppi <- pp[pp_cells==i,]
          obs <- stat(dati)
          pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
          tab[i,] <- c(obs,quantile(pred,quants))
        }
        if (do_plot) {
          dens <- density(pred,adjust=adjust)
          if (!is.null(xlim)) xlimi <- xlim else
            xlimi <- c(pmin(obs,min(dens$x)),pmax(obs,max(dens$x)))
          plot(dens,main=paste(main,i),xlab=stat_name,xlim=xlimi)
          abline(v=obs)
        }
      }
      invisible(tab)
    } else { # cdf
      if (any(is.na(fnams))) cells <- ""
      pok <- probs %in% q_points
      R <- levels(dat$R)
      if (is.null(ylim)) ylim <- c(0,1)
      # set common xlim
      if (is.null(xlim)) {
        xlim <- c(Inf,-Inf)
        for (i in sort(unique(cells))) {
          if (i=="") {
            dati <- dat
            ppi <- pp
          } else {
            dati <- dat[cells==i & okd,]
            ppi <- pp[pp_cells==i & okpp,]
          }
          pqs <- pq <- qs <- setNames(vector(mode="list",length=length(R)),R)
          for (j in R) if (length(dati$rt[dati$R==j])>=length(q_points)) {
            qs[[j]] <- quantile(dati$rt[dati$R==j],probs=probs)
            pq[[j]] <- quantile(ppi$rt[ppi$R==j],probs=probs)
            pqs[[j]] <- tapply(ppi$rt[ppi$R==j],ppi$postn[ppi$R==j],
                               quantile,probs=probs[pok])
          } else qs[[j]] <- pq[[j]] <- pqs[[j]] <- NA
          rx <- cbind(do.call(rbind,lapply(qs,function(x){x[c(1,length(probs))]})),
                      do.call(rbind,lapply(pq,function(x){x[c(1,length(probs))]})))
          xlimi <- c(min(rx,na.rm=TRUE),max(rx,na.rm=TRUE))
          if (!any(is.na(xlimi))) {
            xlim[1] <- pmin(xlim[1],xlimi[1])
            xlim[2] <- pmax(xlim[2],xlimi[2])
          }
        }
      }
      for (i in sort(unique(cells))) {
        if (i=="") {
          dati <- dat
          ppi <- pp
          okdi <- okd
          okppi <- okpp
        } else {
          dati <- dat[cells==i,]
          ppi <- pp[pp_cells==i,]
          okdi <- okd[cells==i]
          okppi <- okpp[pp_cells==i]
        }
        pR <- tapply(okdi,dati$R,sum)/dim(dati)[1]
        ppR <- tapply(okppi,ppi$R,sum)/dim(ppi)[1]
        dati <- dati[okdi,]
        ppi <- ppi[okppi,]
        pqs <- pq <- qs <- setNames(vector(mode="list",length=length(R)),R)
        for (j in R) if (length(dati$rt[dati$R==j])>=length(q_points)) {
          isj <- ppi$R==j
          qs[[j]] <- quantile(dati$rt[dati$R==j],probs=probs)
          pq[[j]] <- quantile(ppi$rt[isj],probs=probs)
          pqs[[j]] <- tapply(ppi$rt[isj],ppi$postn[isj],quantile,probs=probs[pok])
        } else qs[[j]] <- pq[[j]] <- pqs[[j]] <- NA
        if ( !any(is.na(pq[[1]])) ) {
          plot(pq[[1]],probs*ppR[1],xlim=xlim,ylim=ylim,main=paste(main,i),xlab="RT",type="l",
               lwd=fit_lwd,ylab="p(R)",lty=1)
          tmp=lapply(pqs[[1]],function(x){
            points(x,probs[pok]*ppR[1],col="grey",pch=16,cex=pqp_cex)})
          points(pq[[1]][pok],probs[pok]*ppR[1],cex=pqp_cex*3,pch=16,col="grey")
          lines(qs[[1]],probs*pR[1],lwd=data_lwd,lty=1)
          points(qs[[1]][pok],probs[pok]*pR[1],cex=qp_cex,pch=16)
          do_plot=FALSE
        } else do_plot=TRUE
        if (length(qs)>1) {
          for (j in 2:length(qs)) if (!any(is.na(pq[[j]]))) {
            if (do_plot) {
              plot(pq[[j]],probs*ppR[j],xlim=xlim,ylim=ylim,main=paste(main,i),xlab="RT",type="l",
                   lwd=fit_lwd,ylab="p(R)",lty=j)
              do_plot <- FALSE
            } else lines(pq[[j]],probs*ppR[j],lwd=fit_lwd,lty=j)
            tmp=lapply(pqs[[j]],function(x){
              points(x,probs[pok]*ppR[j],col="grey",pch=16,cex=pqp_cex)})
            points(pq[[j]][pok],probs[pok]*ppR[j],cex=pqp_cex*3,pch=16,col="grey")
            lines(qs[[j]],probs*pR[j],lwd=data_lwd,lty=j)
            points(qs[[j]][pok],probs[pok]*pR[j],cex=qp_cex,pch=16)
          }
        }
        legend(lpos,as.character(R),lty=1:length(R),bty="n",title="Response")
      }
    }
  }
}


plot_trials <- function(data,pp=NULL,subject=NULL,factors=NULL,Fcovariates=NULL,
                        layout=NULL,mfcol=TRUE)
  # Plots rt as a function of trials, returns correlation with average predicted
{
  if (!any(names(data)=="trials")) stop("data must have trials column")
  if (!is.null(pp)) {
    if (!any(names(pp)=="trials")) stop("posterior predictives must have trials column")
  } else OvsP <- FALSE
  if (!is.null(subject)) {
    snams <- levels(data$subjects)
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% snams)) stop("Subject(s) not present\n")
    dat <- data[data$subjects %in% subject,]
    dat$subjects <- factor(as.character(dat$subjects))
    pp <- pp[pp$subjects %in% subject,]
    if (length(subject>1))
      fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))] else
        fnams <- names(dat)[!(names(dat) %in% c("subjects","trials","R","rt"))]
  } else {
    dat <- data
    fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt",Fcovariates))]
  }
  if (!is.null(factors)) {
    if (!all(factors %in% fnams)) stop("factors must name factors in data")
    if (!any(factors=="subjects")) stop("factors must include subjects")
    fnams <- factors
  }
  if (!is.null(layout)) if (mfcol) par(mfcol=layout) else par(mfrow=layout)
  postn <- unique(pp$postn)
  rlist <- vector(mode="list",length=length(levels(dat$subjects)))
  names(rlist) <- levels(dat$subjects)
  for (s in levels(dat$subjects)) {
    sdat <- dat[dat$subjects==s,]
    cells <- sdat[,fnams,drop=FALSE]
    for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
    cells <- apply(cells,1,paste,collapse=" ")
    if (!is.null(pp)) {
      spp <- pp[pp$subjects==s,]
      pp_cells <- spp[,fnams,drop=FALSE]
      for (i in fnams) pp_cells[,i] <- paste(i,pp_cells[,i],sep="=")
      pp_cells <- apply(pp_cells,1,paste,collapse=" ")
    }
    ucells <- sort(unique(cells))
    rlist[[s]] <- setNames(vector(mode="numeric",length=length(ucells)),ucells)
    for (i in ucells) {
      obs <- sdat[cells==i,]
      obs <- obs[order(obs$trials),]
      plot(obs$trials,obs$rt,type="l",xlab="Trials",ylab="RT (s)",main=i)
      if (!is.null(pp)) {
        ppi <- spp[pp_cells==i,]
        ppi <- ppi[order(ppi$trials),]
        avpred <- tapply(ppi$rt,ppi$trials,mean)
        lines(obs$trials,avpred,col="red")
        r <- cor(obs$rt,avpred)
        legend("top",paste("r = ",round(r,2)),bty="n")
        rlist[[s]][i] <- r
      }
    }
  }
  invisible(rlist)
}




#' Convergence checks for an EMC2 samplers object
#'
#' Runs a series of convergence checks, prints statistics to the console, and
#' saves traceplots to a PDF file.
#'
#' Note that the `Rhat` is calculated by doubling the number of chains by
#' first splitting chains into first and second half, so it also a test of
#' stationarity.
#'
#' Efficiency of sampling is (optionally) indicated by integrated
#' autocorrelation time (from the `LaplacesDemon` R package) and by the effective
#' number of samples (from the `coda` R package).
#'
#' @param samples A list of samplers or samplers converted to mcmc objects.
#' @param pdf_name A character string, indicating the path to and name of the PDF file.
#' @param interactive A Boolean, defaults to `TRUE`. Indicates whether console output between
#' hierarchical parameter types should be paused or not.
#' @param filter A character string, indicating the sampling stage to examine.
#' Can be `preburn`, `burn`, `adapt`, or `sample`.
#' @param subfilter An integer or vector. If an integer is supplied, iterations up until
#' that integer are examined. If a vector is supplied, everything in that range is examined.
#' @param thin An integer. Keeps only iterations that are a multiple of thin.
#' @param layout A vector specifying the layout as in `par(mfrow = layout)`.
#' If `NA` (i.e., the default), it uses `coda` defaults, if it is `NULL`, it uses
#' the current layout.
#' @param width PDF page width in inches.
#' @param height PDF page height in inches.
#' @param print_IAT Include a print out of the Integrated Autocorrelation Time.
#' @param subject An integer or character string denoting a particular subject.
#' Creates convergence diagnostics only for that subject within a hierarchical model.
#' @return A PDF written to `pdf_name`.
#' @examples \dontrun{
#' check_run(samplers)
#' }
#' @export

check_run <- function(samples,pdf_name="check_run.pdf",interactive=TRUE,
                      filter="sample",subfilter=0,thin=1,print_IAT=FALSE,subject=NULL,
                      layout=c(3,4),width=NULL,height=NULL) {
  print(chain_n(samples))
  pdf(pdf_name,width=width,height=height)
  plot_chains(samples,selection="LL",layout=layout,filter=filter,
              subfilter=subfilter,thin=thin)

  if (is.null(subject) & any(names(samples[[1]]$samples)=="theta_mu")) {
    if (interactive) readline("Enter for mu check")
    cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("\nR-hat\n")
    print(round(gd_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter,
                        print_summary = FALSE,thin=thin),2))
    if (print_IAT) {
      cat("\nIntegrated autocorrelation time\n")
      iat_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter,thin=thin)
    }
    cat("\nEffective Size\n")
    es_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter,thin=thin)
    plot_chains(samples,selection="mu",layout=layout,filter=filter,
                subfilter=subfilter,thin=thin)
    if (interactive) readline("Enter for variance check")
    cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!! VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("\nR-hat\n")
    print(round(gd_pmwg(samples,selection="variance",filter=filter,thin=thin,
                        subfilter=subfilter,print_summary = FALSE),2))
    if (print_IAT) {
      cat("\nIntegrated autocorrelation time\n")
      iat_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter,thin=thin)
    }
    cat("\nEffective Size\n")
    round(es_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter,thin=thin))
    plot_chains(samples,selection="variance",layout=layout,filter=filter,
                subfilter=subfilter,thin=thin)
    if (interactive) readline("Enter for correlation check")
    cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!! CORRELATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("\nR-hat\n")
    print(round(gd_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter,
                        print_summary = FALSE,thin=thin),2))
    if (print_IAT) {
      if (interactive) readline("Enter for next correlation check")
      cat("\nIntegrated autocorrelation time\n")
      iat_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter,thin=thin)
    }
    if (interactive) readline("Enter for next correlation check")
    cat("\nEffective Size\n")
    round(es_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter,thin=thin))
    plot_chains(samples,selection="correlation",layout=layout,ylim=c(-1,1),filter=filter,
                subfilter=subfilter,thin=thin)
    if (interactive) readline("Enter for alpha check")
  }
  cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALPHA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("\nR-hat\n")
  gd <- gd_pmwg(samples,selection="alpha",filter=filter,subfilter=subfilter,
                      print_summary = FALSE,thin=thin)
  if (!is.null(subject)) {
    if (is.numeric(subject)) {
      if (subject>dim(gd)[1] | subject<1) stop("subject value not in range")
      subject <- dimnames(gd)[[1]][round(subject)]
    }
    if (!(subject %in% dimnames(gd)[[1]])) stop("not a subject name")
    ok <- rownames(gd) == subject
  } else ok <- rep(TRUE,dim(gd)[1])
  print(round(gd[ok,],2))
  if (print_IAT) {
    if (is.null(subject)) {
      cat("\nMean Integrated autocorrelation time\n")
      iat_pmwg(samples,selection="alpha",filter=filter,subfilter=subfilter,thin=thin)
    } else {
      cat("\nIntegrated autocorrelation time\n")
      print(round(iat_pmwg(samples,selection="alpha",filter=filter,subfilter=subfilter,
               thin=thin,summary_alpha=NULL,print_summary=FALSE)[ok,],2))
    }
  }
  if (any(names(samples[[1]]$samples)=="theta_mu")) {
    if (is.null(subject)) {
      cat("\nEffective Size (minimum)\n")
      es_pmwg(samples,selection="alpha",summary_alpha=min,filter=filter,
                    subfilter=subfilter,thin=thin)
      cat("\nEffective Size (mean)\n")
      es_pmwg(samples,selection="alpha",summary_alpha=mean,filter=filter,
                  subfilter=subfilter,thin=thin)
    } else {
      cat("\nEffective Size\n")
      print(round(es_pmwg(samples,selection="alpha",summary_alpha=NULL,filter=filter,
                    print_summary=FALSE,subfilter=subfilter,thin=thin)[ok,]))
    }
  } else {
    cat("\nEffective Size\n")
    es_pmwg(samples,selection="alpha",summary_alpha=mean,filter=filter,
                    subfilter=subfilter,thin=thin)
  }
  plot_chains(samples,selection="alpha",layout=layout,filter=filter,
              subfilter=subfilter,thin=thin)
  dev.off()
  message("\n\nGraphical checks available in ",pdf_name)
}


#' Plots matrix of parameter correlations (upper triangle) and corresponding
#' scatterplots (lower triangle)
#'
#' @param samples List of pmwgs objects
#' @param filter A string. Specifies which stage you want to plot.
#' @param thin An integer. Keep only iterations that are a multiple of thin.
#' @param subfilter An integer or vector. If integer it will exclude up until
#' that integer. If vector it will include everything in that range.
#' @param selection String designating parameter type (mu, variance, correlation, alpha = default)
#' @param mapped Boolean (default FALSE) if TRUE plot parameters mapped to design
#' otherwise sampled parameters
#' @param scale.subjects Boolean (default TRUE) when plotting alphas standardize each participant
#' @param use Integer or character vector specifying subset of parameters to use
#' @param do_plot Boolean, do plot
#' @param maxp Integer for maximum number of iterations used (default 500).
#'
#' @return None
#' @export

pairs_posterior <- function(samples,filter="sample",thin=1,subfilter=0,mapped=FALSE,
                            selection=c("alpha","mu","variance","covariance","correlation")[1],
                            scale.subjects=TRUE,use=NA,do_plot=TRUE,maxp=500)
{

  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
  }

  ## put correlations on the upper panels,
  panel.cor <- function(x, y, prefix = "", ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y),3)
    txt <- format(c(r, 0.123456789), digits = 3)[1]
    txt <- paste0(prefix, txt)
    text(0.5, 0.5, txt, cex = 2)
  }

  do_scale <- function(df) {
    for (i in levels(df[,1])) {
      isin <- df[,1]==i
      df[isin,2] <- (df[isin,2]-mean(df[isin,2]))/sd(df[isin,2])
    }
    df[,2]
  }

  pmat <- parameters_data_frame(samples,filter=filter,thin=thin,subfilter=subfilter,
                                mapped=mapped,selection=selection)
  if (!any(is.na(use))) {
    if (is.numeric(use)) {
      if (any(use<1) || any(use>dim(pmat))) stop("stat outside parameter range")
      use <- names(pmat)[use]
    }
    if (!all(use %in% names(pmat))) stop("stat has a name not in parameters")
    if (length(use)==1) stop("must select more than one parameter")
    pmat <- pmat[,use]
  }
  if (selection=="alpha") {
    if (length(levels(pmat$subjects))>1 && scale.subjects)
      for (i in names(pmat)[-1]) pmat[,i] <- do_scale(pmat[,c("subjects",i)])
    pmat <- pmat[,-1]
  }
  if (dim(pmat)[1]>maxp) pmat <- pmat[sample(dim(pmat)[1],maxp),]
  if (do_plot) suppressWarnings(
    pairs(pmat,diag.panel = panel.hist,upper.panel = panel.cor))
  rs <- cor(pmat)
  r.names <- outer(dimnames(rs)[[1]],dimnames(rs)[[2]],paste,sep="~")[upper.tri(rs)]
  rs <- rs[upper.tri(rs)]
  names(rs) <- r.names
  invisible(rs)
}

#' Creates a likelihood profile plot from a dadm object (see design_model) by
#' varying one model parameter while holding all others constant.
#'
#' @param pname Name of parameter to profile
#' @param p Named vector of parameter values (typically created with sampled_p_vector)
#' @param p_min Minimum value of profile range
#' @param p_max Maximum value of profile range
#' @param dadm Data augmented design and model (created by design_model)
#' @param n_point Number of evenly spaced points at which to calculate likelihood
#' @param main Plot title
#' @param cores Number of likelihood points to calculate in parallel
#'
#' @return vector with value of p(pname), highest likelihood point and p(pname)
#' minus the parameter values at that point
#' @export

profile_plot <- function(pname,p,p_min,p_max,dadm,n_point=100,main="",cores=1)

{

  lfun <- function(i,x,p,pname,dadm) {
    p[pname] <- x[i]
    attr(dadm,"model")()$log_likelihood(p,dadm)
  }

  x <- seq(p_min,p_max,length.out=n_point)
  ll <- unlist(mclapply(1:n_point,lfun,dadm=dadm,x=x,p=p,pname=pname,mc.cores = cores))
  plot(x,ll,type="l",xlab=pname,ylab="LL",main=main)
  abline(v=p[pname])
  c(true=p[pname],max=x[which.max(ll)],miss=p[pname]-x[which.max(ll)])
}
