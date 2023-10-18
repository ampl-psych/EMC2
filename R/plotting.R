
#' Plot MCMC chains
#'
#' @param pmwg_mcmc A list of samplers or samplers converted to mcmc objects.
#' @param layout A vector specifying the layout as in par(mfrow = layout).
#' If NA (defualt) use CODA defaults, if NULL use current.
#' @param subject Integer or character vector, if selection = "alpha" picks out subject(s)
#' (defualt NA plot all)
#' @param ylim The y limits of the chain plot.
#' @param selection String designating parameter type (mu, variance, correlation, alpha = default)
#' @param filter A string. Specifies which stage you want to plot.
#' @param thin An integer. Keep only iterations that are a multiple of thin.
#' @param subfilter An integer or vector. If integer it will exclude up until
#' that integer. If vector it will include everything in that range.
#' @param plot_acf Boolean (default FALSE). If TRUE will also plot autocorrelation
#' for the chain specified in acf_chain.
#' @param acf_chain Integer. For which chain to plot the acf, if plot_acf = TRUE.
#' @return None
#' @export

plot_chains <- function(pmwg_mcmc,layout=NA,subject=NA,ylim=NULL,
                        selection="alpha",filter="sample",thin=1,subfilter=0,
                        plot_acf=FALSE,acf_chain=1) # ,use_par=NA
  # Plots chains  (if alpha, LL or epsilon can do individual subject, all by default)
{
  if (!(inherits(pmwg_mcmc,  c("mcmc","mcmc.list")))) {
    if (inherits(pmwg_mcmc, "pmwgs"))
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter) else
                             pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                                       thin=thin,subfilter=subfilter)
  }
  auto.layout <- any(is.na(layout))
  no_layout <- is.null(layout)
  if (!auto.layout & !no_layout) par(mfrow=layout)
  if (attr(pmwg_mcmc,"selection")=="alpha" || attr(pmwg_mcmc,"selection")=="random") {
    snams <- names(pmwg_mcmc)
    if (any(is.na(subject))) subject <- snams
    if (is.numeric(subject)) subject <- snams[subject]
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

#' Calls plot_chains to plot auto-correlation functions for each parameter
#'
#' @param samples Single pmwgs object (in which case ACF for one chain is plotted)
#' or list of pmwgs objects (in which case ACFs for each chain plotted)
#' @param layout A 2-vector specifying the layout handled as in plot_chains
#' @param subject integer or character vector, if selection = "alpha" picks out subjects(s)
#' @param selection String designating parameter type (mu, variance, correlation, alpha = default)
#' @param filter which stage to plot (default "sample")
#' @param subfilter an integer or vector. If integer it will exclude up until that integer.
#' If vector it will include everything in that range.
#'
#' @return None
#' @export
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
#' Each panel contains a set of densities (i.e., densities for each possible
#' possible response) that are defective (i.e., have areas potentially less
#' than 1, where for all responses the area sums to 1).
#'
#' @param data data frame with at least subjects (subjects factor) R (response factor)
#' and rt (response time) columns,and optionally other factor columns with any name except
#' subjects, R, rt or trials.
#' @param subject Integer or string selecting a subject (default NULL = all).
#' @param factors character vector of factor names in design (default NULL = all).
#' @param layout 2-vector specifying par(mfrow) or par(mfcol) (default NULL use current).
#' @param mfcol Boolean, default TRUE use mfcol else mfrow.
#' @param xlim x-axis limit for all cells (default NULL = scale per cell).
#' @param bw number or string bandwidth for density (default "nrd0").
#' @param adjust density function bandwidth adjust parameter.
#' @param correct_fun function scoring accuracy using columns in data.
#' @param rt legend function position string for mean RT (default "top)
#' @param accuracy legend function position string for accuracy (default "topright")
#'
#' @return Invisibly if correct_fun specified a subject accuracy vector
#' @export

plot_defective_density <- function(data,subject=NULL,factors=NULL,
                                   layout=NULL,mfcol=FALSE,
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
    if (mfcol) par(mfcol=layout) else par(mfrow=layout)
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


# layout=c(2,3);selection="alpha";filter="sample";thin=1;subfilter=0;mapped=FALSE
# plot_prior=TRUE;n_prior=1e3;xlim=NULL;ylim=NULL;
# show_chains=FALSE;do_plot=TRUE;subject=NA;add_means=FALSE;
# pars=NULL;probs=c(.025,.5,.975);bw = "nrd0";adjust = 1
# subject=1; filter="burn"; subfilter=300

#' Plots density for parameter estimates.
#'
#' @param pmwg_mcmc A list of samplers or samplers converted to mcmc objects.
#' @param layout A 2-vector specifying the layout as in par(mfrow = layout).
#' If NA or NULL use current.
#' @param selection String designating parameter type (mu, variance, correlation, alpha = default)
#' @param use_par Character vector of names of parameters to plot (default NULL = plot all)
#' @param filter A string. Specifies which stage you want to plot.
#' @param thin An integer. Keep only iterations that are a multiple of thin.
#' @param subfilter An integer or vector. If integer it will exclude up until
#' that integer. If vector it will include everything in that range.
#' @param mapped Boolean (default FALSE) if TRUE plot parameters mapped to design
#' otherwise sampled parameters
#' @param plot_prior Boolean. Add prior distribution to plot (in red)
#' @param n_prior Number of samples to approximate prior (default = 1e3)
#' @param xlim x-axis plot limit, 2-vector (same for all) or matrix (one row for each parameter)
#' @param ylim y-axis plot limit, 2-vector (same for all) or matrix (one row for each parameter)
#' @param prior_xlim A vector giving upper and lower quantiles of prior when choosing xlim.
#' @param show_chains Boolean (default FALSE) plot separate density for each chain.
#' @param do_plot Boolean (default TRUE) do plot
#' @param subject Integer or character vector, if selection = "alpha" picks out
#' subject(s) (default NA plots all).
#' @param add_means Boolean (default FALSE) add parameter means as an attribute
#' to return
#' @param pars Named vector or matrix of true parameters, or the output of
#' plot_pars, in which case the posterior medians are extracted. If supplied with
#' no xlim the plot will be adjusted to include the parameter values, which are
#' plotted as vertical lines.
#' @param probs Vector (default c(.025,.5,.975)) for CI and central tendency of return
#' @param bw Bandwidth for density plot (see density)
#' @param adjust Adjustment for density plot (see density)
#' @param do_contraction Print the prior to posterior contraction (1 - var(prior)/var(posterior))
#' @param lpos Position of contraction in graph
#' @param digits Rounding of contraction
#'
#' @return Invisibly returns tables of true and 95% CIs (for all chains combined
#'no matter what show_chains is), if do_contraction with a "contraction" attribute.
#'
#' @export
plot_pars <- function(pmwg_mcmc,layout=c(2,3),use_par=NULL,
  selection="alpha",filter="sample",thin=1,subfilter=0,mapped=FALSE,
  plot_prior=TRUE,n_prior=1e3,xlim=NULL,ylim=NULL,prior_xlim=NULL,
  show_chains=FALSE,do_plot=TRUE,subject=NA,add_means=FALSE,
  pars=NULL,probs=c(.025,.5,.975),bw = "nrd0", adjust = 1,
  do_contraction=TRUE,lpos="topright",digits=3)
  #  (if alpha can do individual subject, all by default)
  # If show_chains superimposes densities for each chain on same plot
  #
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

  if (mapped & !(selection %in% c("mu","alpha")))
    stop("mapped only available for mu and alpha")

  if (show_chains & plot_prior)
    warning("Prior plots not implemented for show_chains=TRUE")

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
            plot(dens[[1]],xlab=j,main=paste0(attr(pmwg_mcmc,"selection")," s",i),
                 col=1,xlim=xlimi,ylim=ylimi)
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
            plot(dens,xlab=j,xlim=xlimi,ylim=ylimi,
                 main=paste0(attr(pmwg_mcmc,"selection")," s",i))
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
    if (inherits(pmwg_mcmc, "mcmc.list")) {
      selection <- attr(pmwg_mcmc,"selection")
      pmwg_mcmc_combined <- do.call(rbind,pmwg_mcmc)
      attr(pmwg_mcmc_combined,"selection") <- selection
      if (show_chains) chains <- length(pmwg_mcmc) else
        pmwg_mcmc <- pmwg_mcmc_combined
    } else pmwg_mcmc_combined <- pmwg_mcmc
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
      if (length(pars) != dim(pmwg_mcmc_combined)[2])
        stop("pars is wrong length")
      names(pars) <- colnames(pmwg_mcmc_combined)
    }
    tabs <- rbind(true=pars,apply(pmwg_mcmc_combined,2,quantile,probs=probs))
    if (!is.null(pars)) tabs <- rbind(tabs,Miss=tabs[3,]-tabs[1,])
    if (add_means) attr(tabs,"mean") <- apply(pmwg_mcmc_combined,2,mean)
    if (!no_layout) par(mfrow=layout)
    if (plot_prior) dimnames(psamples) <- list(NULL,colnames(pmwg_mcmc_combined))
    if (do_contraction)
      contraction <- setNames(numeric(length(colnames(pmwg_mcmc_combined))),colnames(pmwg_mcmc_combined))
    if (!is.null(use_par)) {
      ok <- colnames(pmwg_mcmc_combined) %in% use_par
      if (!any(ok)) stop("use_par did not specify parameters that are present")
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
        plot(dens[[1]],xlab=attr(pmwg_mcmc,"selection"),main=j,
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
        plot(dens,xlim=xlimi,ylim=ylimi,xlab=attr(pmwg_mcmc,"selection"),main=j)
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


#' Plots data and fit. If rt  available plots cdf, if not plots ROC. If stat
#' argument (which calculates a statistics based on the data) is supplied
#' instead fit is plotted as a density with a vertical line at the position of
#' the data statistic. If more than one subject is included data and fits are
#' aggregated over subjects. If data or fit contains NA in responses or rt, or
#' is.infinite(rt) these are treating as missing and defective cdfs sum to the
#' probability of non-missing.
#'
#' @param data Data frame with subjects and R factors, and possibly other factors
#' and an rt column
#' @param pp Posterior predictives created by post_predict
#' @param subject Integer or string picking out subject(s).
#' @param factors Character vector of factors in data to display separately. If
#' NULL (default) use names of all columns in data except "trials","R", and "rt".
#' Omitted factors are aggregated over.
#' @param stat A function that takes a the data and returns a single value.
#' @param stat_name A string naming what the stat argument calculates.
#' @param ci Credible interval and central tendency quantiles for return when
#' stat argument is supplied (default c(.025,.5,.975))
#' @param do_plot Boolean (default TRUE) for making a plot
#' @param xlim x-axis plot limit, 2-vector (same for all) or matrix (one row for each paramter)
#' @param ylim y-axis plot limit, 2-vector (same for all) or matrix (one row for each paramter)
#' @param layout 2-vector specifying par(mfrow) or par(mfcol) (default NULL use current).
#' @param mfcol Boolean, default TRUE use mfcol else mfrow.
#' @param probs Vector of probabilities at which to calculate cdf (default percentiles)
#' @param data_lwd Integer line width for data in cdf (default = 2)
#' @param fit_lwd Integer line widht for fit in cdf (default = 1)
#' @param qp_cex cex for data quantile points in cdf
#' @param q_points Quantile points to plot in cdf (default c(.1,.3,.5,.7,.9))
#' @param pqp_cex cex for predicted quantile points in cdf
#' @param lpos Legend position (see legend)
#' @param signalFactor The name of the "signal" factor in an ROC.
#' @param zROC Boolean, (default FALSE) to plot ROC on transformed scale.
#' @param qfun Scale transform function for zROC (default qnorm)
#' @param lim 2-vector for x and y limit of ROC
#' @param rocfit_cex cex for predicted ROC points (default = 0.5.
#' @param adjust Control of smoothing in density plots
#'
#' @return If stat argument is provided a table of observed values and predicted quantiles
#' @export
plot_fit <- function(data,pp,subject=NULL,factors=NULL,
                     stat=NULL,stat_name="",adjust=1,
                     ci=c(.025,.5,.975),do_plot=TRUE,
                     xlim=NULL,ylim=NULL,
                     layout=NULL,mfcol=TRUE,
                     probs=c(1:99)/100,
                     data_lwd=2,fit_lwd=1,qp_cex=1,
                     q_points=c(.1,.3,.5,.7,.9),pqp_cex=.5,lpos="topleft",
                     signalFactor="S",zROC=FALSE,qfun=qnorm,lim=NULL,rocfit_cex=.5)
{
  if (!is.null(subject)) {
    snams <- levels(data$subjects)
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% snams)) stop("Subject(s) not present\n")
    dat <- droplevels(data[data$subjects %in% subject,])
    pp <- droplevels(pp[pp$subjects %in% subject,])
    if (length(subject>1))
      fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))] else
      fnams <- names(dat)[!(names(dat) %in% c("subjects","trials","R","rt"))]
  } else {
    dat <- data
    fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))]
  }

  okd <- !is.na(dat$R) & is.finite(dat$rt)
  okpp <- !is.na(pp$R) & is.finite(pp$rt)

  if (!is.null(factors)) {
    if (!all(factors %in% fnams))
      stop("factors must name factors in data")
    fnams <- factors
  }
  if (!is.null(layout)) if (mfcol) par(mfcol=layout) else par(mfrow=layout)
  if (all(is.na(data$rt))) {  # type=SDT
    if (length(levels(data$R))==2 & is.null(stat))
      stop("No plots for binary responses, use an accuracy function in stat arguement.")
    if (!is.null(stat)) { # statistic
      cells <- dat[,fnams,drop=FALSE]
      for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
      cells <- apply(cells,1,paste,collapse=" ")
      pp_cells <- pp[,fnams,drop=FALSE]
      for (i in fnams) pp_cells[,i] <- paste(i,pp_cells[,i],sep="=")
      pp_cells <- apply(pp_cells,1,paste,collapse=" ")
      postn <- unique(pp$postn)
      ucells <- sort(unique(cells))
      tab <- matrix(nrow=length(ucells),ncol=4,
                    dimnames=list(ucells,c("Observed",names(quantile(1:5,ci)))))
      for (i in ucells) {
        obs <- stat(dat[cells==i,])
        ppi <- pp[pp_cells==i,]
        pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
        if (do_plot) {
          dens <- density(pred,adjust=adjust)
          if (!is.null(xlim)) xlimi <- xlim else
            xlimi <- c(pmin(obs,min(dens$x)),pmax(obs,max(dens$x)))
          plot(dens,main=i,xlab=stat_name,xlim=xlimi)
          abline(v=obs)
        }
        tab[i,] <- c(obs,quantile(pred,ci))
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
      cells <- dat[,fnams,drop=FALSE]
      for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
      cells <- apply(cells,1,paste,collapse=" ")
      pp_cells <- pp[,fnams,drop=FALSE]
      for (i in fnams) pp_cells[,i] <- paste(i,pp_cells[,i],sep="=")
      pp_cells <- apply(pp_cells,1,paste,collapse=" ")
      postn <- unique(pp$postn)
      ucells <- sort(unique(cells))
      for (i in ucells) {
        dpts <- plot_roc(dat[cells==i,],zROC=zROC,qfun=qfun,lim=lim,main=i,signalFactor=signalFactor)
        tab <- table(pp[pp_cells==i,]$postn,pp[pp_cells==i,]$R,pp[pp_cells==i,][[signalFactor]])
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
    cells <- dat[,fnams,drop=FALSE]
    for (i in fnams) cells[,i] <- paste(i,cells[,i],sep="=")
    cells <- apply(cells,1,paste,collapse=" ")
    pp_cells <- pp[,fnams,drop=FALSE]
    for (i in fnams) pp_cells[,i] <- paste(i,pp_cells[,i],sep="=")
    pp_cells <- apply(pp_cells,1,paste,collapse=" ")
    if (!is.null(stat)) { # statistic
      postn <- unique(pp$postn)
      ucells <- sort(unique(cells))
      tab <- matrix(nrow=length(ucells),ncol=4,
                    dimnames=list(ucells,c("Observed",names(quantile(1:5,ci)))))
      for (i in ucells) {
        obs <- stat(dat[cells==i,])
        ppi <- pp[pp_cells==i,]
        pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
        if (do_plot) {
          dens <- density(pred,adjust=adjust)
          if (!is.null(xlim)) xlimi <- xlim else
            xlimi <- c(pmin(obs,min(dens$x)),pmax(obs,max(dens$x)))
          plot(dens,main=i,xlab=stat_name,xlim=xlimi)
          abline(v=obs)
        }
        tab[i,] <- c(obs,quantile(pred,ci))
      }
      invisible(tab)
    } else { # cdf
      pok <- probs %in% q_points
      R <- levels(dat$R)
      if (is.null(ylim)) ylim <- c(0,1)
      # set common xlim
      if (is.null(xlim)) {
        xlim <- c(Inf,-Inf)
        for (i in sort(unique(cells))) {
        dati <- dat[cells==i & okd,]
          ppi <- pp[pp_cells==i & okpp,]
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
        dati <- dat[cells==i,]
        ppi <- pp[pp_cells==i,]
        okdi <- okd[cells==i]
        okppi <- okpp[pp_cells==i]
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
          plot(pq[[1]],probs*ppR[1],xlim=xlim,ylim=ylim,main=i,xlab="RT",type="l",
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
              plot(pq[[j]],probs*ppR[j],xlim=xlim,ylim=ylim,main=i,xlab="RT",type="l",
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


# subject=NULL;factors=NULL;
# pp=pprdm_B_MT; data=attr(rdm_B_MT,"data_list")[[1]]; Fcovariates="MT", layout=c(2,3)
# pp=pp_MT; data=attr(ddm_a_MT,"data_list")[[1]]; Fcovariates="MT"; layout=c(2,3)
# pp=pp_st0; data=attr(ddm_a_st0,"data_list")[[1]]; Fcovariates="MT"; layout=c(2,3)
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




# pdf_name="check_run.pdf";interactive=TRUE;
#                       filter="sample";subfilter=0
#                       layout=c(3,4);width=NULL;height=NULL
# subfilter=2000

#' Runs a series of convergence checks, printing statistics to the console and
#' saving plots to a pdf. Note that the R_hat (psrf and mpsrf, i.e., gelman_diag
#' from the coda package) is calculated by doubling the number of chains by
#' first splitting chains into first and second half so it also a test of
#' stationarity. Efficiency of sampling is indicated by integrated autocorrelation
#' time (from the LaplacesDemon package) and the effective numebr of samples
#' (from coda).
#'
#' @param samples A list of samplers or samplers converted to mcmc objects.
#' @param pdf_name The name of the plot save file
#' @param interactive Pause console output between hierarchical parameter types
#' @param filter Sampling stage to examine.
#' @param subfilter An integer or vector. If integer it will exclude up until
#' that integer. If vector it will include everything in that range.
#' @param thin An integer. Keep only iterations that are a multiple of thin.
#' @param layout A vector specifying the layout as in par(mfrow = layout).
#' If NA (defualt) use CODA defaults, if NULL use current.
#' @param width PDF page width
#' @param height PDF page height
#'
#' @return None
#' @export

check_run <- function(samples,pdf_name="check_run.pdf",interactive=TRUE,
                      filter="sample",subfilter=0,thin=1,
                      layout=c(3,4),width=NULL,height=NULL) {
  print(chain_n(samples))
  pdf(pdf_name,width=width,height=height)
  plot_chains(samples,selection="LL",layout=layout,filter=filter,
              subfilter=subfilter,thin=thin)
  if (any(names(samples[[1]]$samples)=="theta_mu")) {
  if (interactive) readline("Enter for mu check")
    cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("\nR-hat\n")
    print(round(gd_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter,
                        print_summary = FALSE,thin=thin),2))
    cat("\nIntegrated autocorrelation time\n")
    iat_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter,thin=thin)
    cat("\nEffective Size\n")
    es_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter,thin=thin)
    plot_chains(samples,selection="mu",layout=layout,filter=filter,
                subfilter=subfilter,thin=thin)
    if (interactive) readline("Enter for variance check")
    cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!! VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("\nR-hat\n")
    print(round(gd_pmwg(samples,selection="variance",filter=filter,thin=thin,
                        subfilter=subfilter,print_summary = FALSE),2))
    cat("\nIntegrated autocorrelation time\n")
    iat_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter,thin=thin)
    cat("\nEffective Size\n")
    round(es_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter,thin=thin))
    plot_chains(samples,selection="variance",layout=layout,filter=filter,
                subfilter=subfilter,thin=thin)
    if (interactive) readline("Enter for correlation check")
    cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!! CORRELATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("\nR-hat\n")
    print(round(gd_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter,
                        print_summary = FALSE,thin=thin),2))
    if (interactive) readline("Enter for next correlation check")
    cat("\nIntegrated autocorrelation time\n")
    iat_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter,thin=thin)
    if (interactive) readline("Enter for next correlation check")
    cat("\nEffective Size\n")
    round(es_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter,thin=thin))
    plot_chains(samples,selection="correlation",layout=layout,ylim=c(-1,1),filter=filter,
                subfilter=subfilter,thin=thin)
    if (interactive) readline("Enter for alpha check")
  }
  cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALPHA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("\nR-hat\n")
  print(round(gd_pmwg(samples,selection="alpha",filter=filter,subfilter=subfilter,
                      print_summary = FALSE,thin=thin),2))
  cat("\nIntegrated autocorrelation time\n")
  iat_pmwg(samples,selection="alpha",filter=filter,subfilter=subfilter,thin=thin)
  if (any(names(samples[[1]]$samples)=="theta_mu")) {
    cat("\nEffective Size (minimum)\n")
    round(es_pmwg(samples,selection="alpha",summary_alpha=min,filter=filter,
                  subfilter=subfilter,thin=thin))
    cat("\nEffectvie Size (mean)\n")
    round(es_pmwg(samples,selection="alpha",summary_alpha=mean,filter=filter,
                  subfilter=subfilter,thin=thin))
  } else {
    cat("\nEffective Size\n")
    round(es_pmwg(samples,selection="alpha",summary_alpha=mean,filter=filter,
                  subfilter=subfilter,thin=thin))
  }
  plot_chains(samples,selection="alpha",layout=layout,filter=filter,
              subfilter=subfilter,thin=thin)
  dev.off()
  message("\n\nGraphical checks available in ",pdf_name)
}


# filter="sample";thin=1;subfilter=0;mapped=FALSE
# selection=c("alpha","mu","variance","covariance","correlation")[1]; stat=NULL
# collapse.subjects=TRUE;scale.subjects=TRUE;use=NA;do_plot=TRUE

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
    if (is.numeric(stat)) {
      if (any(stat<1) || any(stat>dim(pmat))) stop("stat outside parameter range")
      stat <- names(pmat)[stat]
    }
    if (!all(stat %in% names(pmat))) stop("stat has a name not in parameters")
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

