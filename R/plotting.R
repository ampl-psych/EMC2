
#' Plot MCMC chains
#'
#' @param pmwg_mcmc A list of samplers or samplers converted to mcmc objects.
#' @param layout a vector or matrix specifying the layout as in par(mfrow = layout).
#' @param subject a string or integer which you can specify to only plot chains of that subject.
#' @param ylim an intenger. The y limits of the chain plot.
#' @param selection a string. Specifies which chains you want to plot.
#' @param filter a string. Specifies which stage you want to plot.
#' @param thin an integer. Specify if you want to thin the chains by a factor n.
#' @param subfilter an integer or vector. If integer it will exclude up until that integer. If vector it will include everything in that range.
#' @param plot_acf bool. If TRUE will also plot autocorrelation for the chain specified in acf_chain.
#' @param acf_chain an integer. For which chain to plot the acf, if plot_acf = TRUE.
#' @return
#' @export
#'
#' @examples
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
                  layout=layout,plot_acf=TRUE,acf_chain=i,subject=subject)
    } else plot_chains(samples[[i]],selection=selection,filter=filter,subfilter=subfilter,
                       layout=layout,plot_acf=TRUE,acf_chain=i)
  }
}

plot_alpha_recovery <- function(tabs,layout=c(2,3),
                                do_ci = TRUE,ci_col="grey",cap=.05,
                                do_rmse=FALSE,rmse_pos="topleft",
                                rmse_digits=3,pearson_digits=2,
                                do_coverage=FALSE,coverage_pos="bottomright",
                                coverage_digits=1,spearman_digits=2)
  # Takes tables output by plot_density with par and plots recovery
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
      legend(rmse_pos,paste("RMSE = ",round(rmse[p],rmse_digits)),bty="n")
    } else {
      pearson[p] <- cor(xy[,"true"],xy[,"50%"],method="pearson")
      legend(rmse_pos,paste("r(pearson) = ",round(pearson[p],pearson_digits)),bty="n")
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
#' @param subject string selecting a subject (default NULL = all).
#' @param factors character vector of factor names in design (default NULL = all).
#' @param layout 2-vector specifying par(mfrow) or par(mfcol) (default NULL use current).
#' @param mfcol boolean, default TRUE use mfcol else mfrow.
#' @param xlim x-axis limit for all cells (default NULL = scale per cell).
#' @param bw number or string bandwidth for density (default "nrd0").
#' @param adjust density function bandwidth adjust parameter.
#' @param correct_fun function scoring accuracy using columns in data.
#' @param rt legend function position string for mean RT (default "top)
#' @param accuracy legend function position string for accuracy (default "topright")
#'
#' @return Invisibly if correct_fun specified a subject accuracy vector
#' @export
#'
#' @examples
plot_defective_density <- function(data,subject=NULL,factors=NULL,
                                   layout=NULL,mfcol=TRUE,
                                   xlim=NULL,bw = "nrd0",adjust=1,
                                   correct_fun=NULL,rt="top",accuracy="topright")
{
  if (!is.null(subject)) {
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
plot_density <- function(pmwg_mcmc,layout=c(2,3),
                         selection="alpha",filter="sample",thin=1,subfilter=0,mapped=FALSE,
                         plot_prior=TRUE,n_prior=1e3,xlim=NULL,ylim=NULL,
                         show_chains=FALSE,do_plot=TRUE,subject=NA,add_means=FALSE,
                         pars=NULL,probs=c(.025,.5,.975),bw = "nrd0", adjust = 1)
  # Plots density (if alpha can do individual subject, all by default)
  # If show_chains superimposes destinies for each chain on same plot
  # invisibly returns tables of true and 95% CIs (for all chains combined
  # no matter what show_chains is)
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

  if (show_chains & plot_prior)
    warning("Prior plots not implemented for show_chains=TRUE")
  if (!is.na(subject) & plot_prior) {
    warning("Can't plot prior for single subject.")
    plot_prior <- FALSE
  }
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
    if (any(is.na(subject))) subject <- names(pmwg_mcmc)
    if (!all(subject %in% names(pmwg_mcmc)))
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
    for (i in subject) {
      if (do_plot) {
        if (!no_layout) par(mfrow=layout)
        for (j in colnames(pmwg_mcmc_combined[[i]])) {
          if (chains>0) {
            dens <- lapply(pmwg_mcmc[[i]],function(x){density(x[,j],bw=bw,adjust=adjust)})
            if (!is.null(xlim)) {
              if (!is.matrix(xlim)) xlimi <- xlim else xlimi <- xlim[j,]
            } else xlimi <- c(min(unlist(lapply(dens,function(x){min(x$x)}))),
                              max(unlist(lapply(dens,function(x){max(x$x)}))))
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
            } else xlimi <- c(min(dens$x),max(dens$x))
            if (!is.null(ylim)) {
              if (!is.matrix(ylim)) ylimi <- ylim else ylimi <- ylim[j,]
            } else {
              ylimi <- c(0,max(dens$y))
              if (plot_prior) ylimi[2] <- max(c(ylimi[2],pdens$y))
            }
            plot(dens,xlab=j,xlim=xlimi,ylim=ylimi,
                 main=paste0(attr(pmwg_mcmc,"selection")," s",i))
            if (plot_prior) lines(pdens,col="red")
          }
          if (!is.null(pars)) abline(v=pars[j,i])
        }
      }
      if (!is.null(pars)) tabs[[i]] <- rbind(true=pars[dimnames(tabs[[i]])[[2]],i],tabs[[i]])
    }
    tabs <- tabs[as.character(subject)]
    if (add_means)
      attr(tabs,"mean") <- lapply(pmwg_mcmc_combined,function(x){apply(x,2,mean)})
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
    if (add_means) attr(tabs,"mean") <- apply(pmwg_mcmc_combined,2,mean)
    if (!no_layout) par(mfrow=layout)
    if (plot_prior) dimnames(psamples) <- list(NULL,colnames(pmwg_mcmc_combined))
    if (do_plot) for (j in colnames(pmwg_mcmc_combined)) {
      if (chains > 0) {
        dens <- lapply(pmwg_mcmc,function(x){density(x[,j],bw=bw,adjust=adjust)})
        if (!is.null(xlim)) xlimi <- xlim else
          xlimi <- c(min(unlist(lapply(dens,function(x){min(x$x)}))),
                     max(unlist(lapply(dens,function(x){max(x$x)}))))
        if (!is.null(ylim)) ylimi <- ylim else
          ylimi <- c(0,max(unlist(lapply(dens,function(x){max(x$y)}))))
        plot(dens[[1]],xlab=attr(pmwg_mcmc,"selection"),main=j,
             col=1,xlim=xlimi,ylim=ylimi)
        if (chains>1) for (k in 2:chains) lines(dens[[k]],col=k)
      } else {
        dens <- density(pmwg_mcmc[,j],bw=bw,adjust=adjust)
        if (plot_prior)  pdens <- robust_density(psamples[,j],range(pmwg_mcmc[,j]),
                                                 bw=bw,adjust=adjust,
                                                 use_robust=!(attr(pmwg_mcmc,"selection") %in% c("mu","correlation")))
        if (!is.null(xlim)) xlimi <- xlim else
          xlimi <- c(min(dens$x),max(dens$x))
        if (!is.null(ylim)) ylimi <- ylim else {
          ylimi <- c(0,max(dens$y))
          if (plot_prior) ylimi[2] <- max(c(ylimi[2],pdens$y))
        }
        plot(dens,xlim=xlimi,ylim=ylimi,xlab=attr(pmwg_mcmc,"selection"),main=j)
        if (plot_prior) lines(pdens,col="red")
      }
      if (!is.null(pars)) abline(v=pars[j])
    }
    invisible(tabs)
  }
}

plot_fit <- function(data,pp,subject=NULL,factors=NULL,
                     stat=NULL,stat_name="",
                     ci=c(.025,.5,.975),do_plot=TRUE,
                     xlim=NULL,ylim=NULL,
                     layout=NULL,mfcol=TRUE,
                     probs=c(1:99)/100,
                     data_lwd=2,fit_lwd=1,qp_cex=1,
                     q_points=c(.1,.3,.5,.7,.9),pqp_cex=.5,lpos="topleft",
                     signalFactor="S",zROC=FALSE,qfun=NULL,lim=NULL,rocfit_cex=.5)
  # Plots data and fit or if stat supplied returns a table of statistics.
  # If rt available plots cdf, if not plots ROC (but if only 2 levels of
  # response factor returns an accuracy table).
{
  if (!is.null(subject)) {
    dat <- data[data$subjects==subject,]
    pp <- pp[pp$subjects==subject,]
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
          dens <- density(pred)
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
        dpts <- plot_roc(dat[cells==i,],zROC=zROC,qfun=qfun,lim=lim,main=i)
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
          dens <- density(pred)
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
          dati <- dat[cells==i,]
          ppi <- pp[pp_cells==i,]
          pR <- table(dati$R)/dim(dati)[1]
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
        # R <- sort(unique(dati$R))
        pR <- table(dati$R)/dim(dati)[1]
        pqs <- pq <- qs <- setNames(vector(mode="list",length=length(R)),R)
        for (j in R) if (length(dati$rt[dati$R==j])>=length(q_points)) {
          qs[[j]] <- quantile(dati$rt[dati$R==j],probs=probs)
          pq[[j]] <- quantile(ppi$rt[ppi$R==j],probs=probs)
          pqs[[j]] <- tapply(ppi$rt[ppi$R==j],ppi$postn[ppi$R==j],
                             quantile,probs=probs[pok])
        } else qs[[j]] <- pq[[j]] <- pqs[[j]] <- NA
        if ( !any(is.na(pq[[1]])) ) {
          plot(pq[[1]],probs*pR[1],xlim=xlim,ylim=ylim,main=i,xlab="RT",type="l",
               lwd=fit_lwd,ylab="p(R)",lty=1)
          tmp=lapply(pqs[[1]],function(x){
            points(x,probs[pok]*pR[1],col="grey",pch=16,cex=pqp_cex)})
          points(pq[[1]][pok],probs[pok]*pR[1],cex=pqp_cex*3,pch=16,col="grey")
          lines(qs[[1]],probs*pR[1],lwd=data_lwd,lty=1)
          points(qs[[1]][pok],probs[pok]*pR[1],cex=qp_cex,pch=16)
          do_plot=FALSE
        } else do_plot=TRUE
        if (length(qs)>1) {
          for (j in 2:length(qs)) if (!any(is.na(pq[[j]]))) {
            if (do_plot) {
              plot(pq[[j]],probs*pR[j],xlim=xlim,ylim=ylim,main=i,xlab="RT",type="l",
                   lwd=fit_lwd,ylab="p(R)",lty=j)
              do_plot <- FALSE
            } else lines(pq[[j]],probs*pR[j],lwd=fit_lwd,lty=j)
            tmp=lapply(pqs[[j]],function(x){
              points(x,probs[pok]*pR[j],col="grey",pch=16,cex=pqp_cex)})
            points(pq[[j]][pok],probs[pok]*pR[j],cex=pqp_cex*3,pch=16,col="grey")
            lines(qs[[j]],probs*pR[j],lwd=data_lwd,lty=j)
            points(qs[[j]][pok],probs[pok]*pR[j],cex=qp_cex,pch=16)
          }
        }
        legend(lpos,as.character(R),lty=1:length(R),bty="n",title="Response")
      }
    }
  }
}


plot_fits <- function(data,pp,factors=NULL,
                      stat=NULL,stat_name="",ci=c(.025,.5,.975),do_plot=TRUE,
                      xlim=NULL,ylim=NULL,
                      layout=NULL,mfcol=TRUE,
                      probs=c(1:99)/100,
                      data_lwd=2,fit_lwd=1,qp_cex=1,
                      q_points=c(.1,.3,.5,.7,.9),pqp_cex=.5,lpos="topleft",
                      signalFactor="S",zROC=FALSE,qfun=NULL,lim=NULL,rocfit_cex=.5){
  # as for plot_fits but does it per subject.
  for (i in levels(data$subjects))
    plot_fit(data,pp,subject=i,factors,stat,stat_name,ci,do_plot,xlim,ylim,layout,mfcol,
             probs,data_lwd,fit_lwd,qp_cex,q_points,pqp_cex,lpos,
             signalFactor,zROC=qfun,lim,rocfit_cex)
}

plot_roc <- function(data,zROC=FALSE,qfun=NULL,main="",lim=NULL)

{
  tab <- table(data$R,data$S)
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
    dat <- data[data$subjects==subject,]
    dat$subjects <- factor(as.character(dat$subjects))
    pp <- pp[pp$subjects==subject,]
    fnams <- names(dat)[!(names(dat) %in% c("subjects","trials","R","rt",Fcovariates))]
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
  es_pmwg(samples,selection="mu",filter=filter,subfilter=subfilter)
  plot_chains(samples,selection="mu",layout=layout,filter=filter,subfilter=subfilter)
  if (interactive) readline("Enter for variance check")
  cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!! VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("\nR-hat\n")
  print(round(gd_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter,print_summary = FALSE),2))
  cat("\nIntegrated autocorrelation time\n")
  iat_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter)
  cat("\nEffectvie Size\n")
  round(es_pmwg(samples,selection="variance",filter=filter,subfilter=subfilter))
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
  round(es_pmwg(samples,selection="correlation",filter=filter,subfilter=subfilter))
  plot_chains(samples,selection="correlation",layout=layout,ylim=c(-1,1),filter=filter,subfilter=subfilter)
  if (interactive) readline("Enter for alpha check")
  cat("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALPHA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("\nR-hat\n")
  print(round(gd_pmwg(samples,selection="alpha",filter=filter,subfilter=subfilter,print_summary = FALSE),2))
  cat("\nIntegrated autocorrelation time\n")
  iat_pmwg(samples,selection="alpha",filter=filter,subfilter=subfilter)
  cat("\nEffectvie Size (minimum)\n")
  round(es_pmwg(samples,selection="alpha",summary_alpha=min,filter=filter,subfilter=subfilter))
  cat("\nEffectvie Size (mean)\n")
  round(es_pmwg(samples,selection="alpha",summary_alpha=mean,filter=filter,subfilter=subfilter))
  plot_chains(samples,selection="alpha",layout=layout,filter=filter,subfilter=subfilter)
  dev.off()
  message("\n\nGraphical checks available in ",pdf_name)
}

profile_pmwg <- function(pname,p,p_min,p_max,dadm,n_point=100,main="",cores=1)

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

