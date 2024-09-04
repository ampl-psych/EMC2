#' Plot defective densities for each subject and cell
#'
#' Plots panels that contain a set of densities for each response option in the data.
#' These densities are defective; their areas are relative to the respective
#' response proportion. Across all responses, the area sums to 1.
#'
#' @param data A data frame. The experimental data in EMC2 format with at least `subject` (i.e., the
#' subject factor), `R` (i.e., the response factor) and `rt` (i.e., response time) variable.
#' Additional factor variables of the design are optional.
#' @param subject An integer or character string selecting a subject from the data.
#' If specified, only that subject is plotted. Per default (i.e., `NULL`), all subjects
#' are plotted.
#' @param factors A character vector of the factor names in the design to aggregate across
#' Defaults to all (i.e., `NULL`).
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param correct_fun If specified, the accuracy for each subject is calculated, using the supplied function and
#' an accuracy vector for each subject is returned invisibly.
#' @param rt_pos legend function position character string for the mean response time (defaults to `top`)
#' @param accuracy legend function position character string for accuracy (defaults to `topright`)
#' @param ... Optional arguments that can be passed to `get_pars`, `density`, or `plot.default` (see `par()`)
#' @return If `correct_fun` is specified, a subject accuracy vector is returned invisibly
#' @examples
#' # First for each subject and the factor combination in the design:
#' plot_defective_density(forstmann)
#' # Now collapsing across subjects:
#' plot_defective_density(forstmann, factors = c("S", "E"))
#' # If the data is response coded, it generally makes sense to include the "S" factor
#' # because EMC2 will plot the "R" factor automatically. This way, choice accuracy can
#' # be examined
#' # Each subject's accuracy can be returned using a custom function:
#' print(plot_defective_density(forstmann, correct_fun = function(d) d$R == d$S))
#'
#' @export

plot_defective_density <- function(data,subject=NULL,factors=NULL,
                                   layout=NA, correct_fun=NULL,
                                   rt_pos="top",accuracy="topright",
                                   ...)
{
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  dots <- list(...)
  if (!is.null(subject)) {
    snams <- levels(data$subjects)
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% snams)) stop("Subject not present\n")
    dat <- data[data$subjects==subject,]
    dat$subjects <- droplevels(dat$subjects)
    fnams <- names(dat)[!(names(dat) %in% c("subjects","trials","R","rt"))]
  } else {
    dat <- data
    fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt", "subjects"))]
  }
  if (!is.null(factors)) {
    if (!all(factors %in% c(fnams, "subjects")))
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
  R <- levels(dat$R)
  if(any(is.na(layout))){
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(unique(cells)),
                              nplots = 1))
  } else{par(mfrow=layout)}
  for (i in sort(unique(cells))) {
    pR <- table(alldat$R[cellsall==i])/dim(alldat[cellsall==i,])[1]
    dati <- dat[cells==i,]
    mrt <- tapply(dati$rt,dati$R,median)
    dens <- setNames(vector(mode="list",length=length(R)),R)
    for (j in R) if (length(dati$rt[dati$R==j])>1) {
      dens[[j]] <- do.call(density, c(list(dati$rt[dati$R==j]), fix_dots(dots, density.default, consider_dots = FALSE)))
      dens[[j]]$y <- dens[[j]]$y*pR[j]
    } else dens[[j]] <- NULL
    rx <- do.call(rbind,lapply(dens,function(x){range(x$x)}))
    cur_dots <- add_defaults(dots, xlim = c(min(rx[,1]),max(rx[,2])),
                             ylim = c(0,max(unlist(lapply(dens,function(x){max(x$y)})))),
                             main = i, xlab = "RT")
    ltys <- c(1:length(mrt))[!is.na(mrt)]
    do.call(plot, c(list(dens[[1]]), fix_dots_plot(cur_dots)))
    if (length(dens)>1) for (j in 2:length(dens))
      lines(dens[[j]],lty=ltys[j])
    if (!is.null(accuracy) && length(R) > 1) {
      ltys <- 1:length(pR)
      ltys[is.na(mrt)] <- NA
      legend(accuracy,paste(names(pR),round(pR,2),sep="="),lty=ltys,bty="n",title="p(R)")
    }
    if (!is.null(rt_pos))
      legend(rt_pos,paste(names(pR),round(mrt,3),sep="="),bty="n",title="Med(RT)")
  }
  if (!is.null(correct_fun))
    invisible(tapply(correct_fun(dat),dat$subjects,mean))
}

robust_hist <- function(ps, breaks = 50, cutoff = 0.0015, do_plot = TRUE, ...){
  ps <- ps[abs(ps) < 1000]
  for(i in 1:log10(length(ps))){
    cuts <- cut(ps, breaks = breaks)
    tab_cuts <- table(cuts)
    good_cuts <- tab_cuts/length(ps) > cutoff
    ps <- ps[cuts %in% names(good_cuts)[good_cuts]]
  }
  if(do_plot){
    hist(ps, breaks = breaks, ...)
  } else{
    return(ps)
  }

}

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

#' Plots choice data
#'
#' Plots choice data with no response times.
#'
#' @param data A data frame. The experimental data in EMC2 format with at least `subject` (i.e., the
#' subject factor), `R` (i.e., the response factor) and `rt` (i.e., response time) variable.
#' Additional factor variables of the design are optional.
#' @param pp Posterior predictives created by `predict()`
#' @param subject Integer or string picking out subject(s).
#' @param factors Character vector of factors in data to display separately. If
#' `NULL` (i.e., the default), use names of all columns in data except `trials`,`R`, and `rt`.
#' Omitted factors are aggregated over. If `NA`, treats entire data set as a single cell.
#' If `stat` is used, the default is changed to `NA`.
#' @param functions A named list of functions that create new factors which can then be
#' used by the `factors` and `stat` arguments.
#' @param stat A function that takes a the data and returns a single value.
#' @param stat_name A string naming what the `stat` argument calculates.
#' @param adjust Control of smoothing in density plots
#' @param ci Credible interval and central tendency quantiles for return when
#' stat argument is supplied (defaults to the 2.5\\%, the 50\\% and the 97.5\\%
#' quantiles)
#' @param do_plot Boolean (defaults to `TRUE`) whether a plot should be created or not
#' @param xlim x-axis plot limit, 2-vector (same for all) or matrix (one row for each paramter)
#' @param ylim y-axis plot limit, 2-vector (same for all) or matrix (one row for each paramter)
#' @param main Text title, pasted before cell name.
#' @param layout 2-vector specifying `par(mfrow)` or `par(mfcol)`. The default `NULL` uses current,
#' `NA` keeps `par` currently active.
#' @param mfcol Boolean for `layout` settings, the default `TRUE` uses `mfcol`, else `mfrow`.
#' @param signalFactor Character name of factor for the signal
#' @param zROC Boolean, plot Z transformed ROC (defaults to `FALSE`)
#' @param qfun Type of Z transform (defaults to probit)
#' @param lim `x` = `y` limit for ROC plots
#' @param rocfit_cex Size of points in ROC plot (default 0.5)
#'
#' @return If stat argument is provided a matrix of observed values and predicted quantiles
#' is returned
plot_fit_choice <- function(data,pp,subject=NULL,factors=NULL,functions=NULL,
                            stat=NULL,stat_name="",adjust=1,
                            ci=c(.025,.5,.975),do_plot=TRUE,
                            xlim=NULL,ylim=NULL,main="",
                            layout=NULL,mfcol=TRUE,
                            signalFactor="S",zROC=FALSE,qfun=qnorm,lim=NULL,rocfit_cex=.5)
{
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  if (!is.null(stat) & is.null(factors)) factors <- NA
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
    if("subjects" %in% factors){
      fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))]
    } else{
      fnams <- names(dat)[!(names(dat) %in% c("subjects", "trials","R","rt"))]
    }
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
    if (mfcol) par(mfcol=layout) else par(mfrow=layout)

  if (!all(is.na(data$rt))) stop("Use plot_fit for rt data")

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
    } else {
      ucells <- ""
      postn <- unique(pp$postn)
    }
    tab <- matrix(nrow=length(ucells),ncol=4,
                  dimnames=list(ucells,c("Observed",names(quantile(1:5,ci)))))
    for (i in ucells) {
      if (i=="") {
        obs <- stat(dat)
        ppi <- pp
        pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
        tab[1,] <- c(obs,quantile(pred,ci))
      } else {
        obs <- stat(dat[cells==i,])
        ppi <- pp[pp_cells==i,]
        pred <- sapply(postn,function(x){stat(ppi[ppi$postn==x,])})
        tab[i,] <- c(obs,quantile(pred,ci))
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
}


#' Posterior predictive checks
#'
#' Plot (defective) cumulative density functions of the observed data and data from the
#' posterior predictive distribution: the probability of a response,
#' p(R) as a function of response time for the experimental data and posterior
#' predictive simulations.
#'
#' The data is plotted in black. Large grey points show the average quantiles across the posterior predictives.
#' The small grey points represent the predicted quantile of an individual replicate,
#' providing a representation of uncertainty in the model predictions.
#'
#' If the stat argument is supplied (which calculates a statistic based on the data),
#' the posterior predictives are plotted as a density over the different replicates.
#' A vertical line is plotted at the value of that statistic for the experimental data.
#'
#' If more than one subject is included, the data and fits are aggregated across subjects by default.
#'
#' Also see `?plot_defective_density()` for more details.

#' @param data A data frame. The experimental data in EMC2 format with at least `subject` (i.e., the
#' subject factor), `R` (i.e., the response factor) and `rt` (i.e., response time) variable.
#' Additional factor variables of the design are optional.
#' @param pp A data frame. Posterior predictives created by ``predict()``
#' @param subject Integer or string selecting a subject from the data. If specified only that subject
#' is plotted. `NULL` (i.e., the default), will plot all subjects.
#' @param factors Character vector of factors in data to display separately. If
#' NULL (default) use names of all columns in data except "trials","R", and "rt".
#' Omitted factors are aggregated over. If NA treats entire data set as a single cell.
#' Must be NA or NULL when using stat argument.
#' @param functions A named list of functions that create new factors which can then be
#' used by the factors and stat arguments.
#' @param stat A function that takes the data/the posterior predictives and returns a single value.
#' For the posterior predictives it will use a single value per replicate, which are then plotted as a density.
#' @param stat_name A string naming what the stat argument calculates, used in labeling the x-axis of the plot.
#' @param adjust Numeric. Density function bandwidth adjust parameter. See ``?density`
#' @param quants A vector. Quantiles of the posterior predictives to return when stat argument is supplied.
#' @param do_plot Boolean. Set to ``FALSE`` to only return the quantiles and omit the plots.
#' @param xlim A numeric vector. x-axis plot limit.
#' @param ylim A numeric vector. y-axis plot limit.
#' @param layout A vector specifying the layout as in `par(mfrow = layout)`.
#' If `NA` or `NULL` uses current plot window layout.
#' @param mfcol Boolean.  If ``TRUE`` uses par(mfrow = layout), otherwise uses par(mfcol = layout)
#' @param probs Vector of probabilities at which to calculate cumulative density function
#' @param data_lwd Integer. Line width for data
#' @param fit_lwd Integer. Line width for posterior predictives
#' @param q_points Vector. Quantile points to plot
#' @param qp_cex Numeric. Cex for data quantile points
#' @param pqp_cex Numeric. Cex for predicted quantile points
#' @param lpos Character. Legend position, see ``?legend()``.
#' @param main Character. Pasted before the plot title, especially useful when specifying a stat argument.
#' @return If stat argument is provided, a vector of observed values and predicted quantiles
#' is returned
#' @examples \donttest{
#' # First generate posterior predictives based on an emc object run with run_emc
#' pp <- predict(samples_LNR, n_cores = 1, n_post = 10)
#' # Then visualize the model fit
#' plot_fit(forstmann, pp, factors = c("S", "E"), layout = c(2,3))
#'
#' # Specific statistics on the posterior predictives can also be specified
#' # This function calculates the difference in rt between two S levels.
#' # It takes the data (or the posterior predictives) as an argument
#' drt <- function(data) diff(tapply(data$rt,data[,c("S")],mean))
#' plot_fit(forstmann, pp, stat=drt,stat_name="Rt difference",
#'          main=("Left vs Right"))
#' }
#' @export
plot_fit <- function(data,pp,subject=NULL,factors=NULL,functions=NULL,
                     stat=NULL,stat_name="",adjust=1,
                     quants=c(.025,.5,.975),do_plot=TRUE,
                     xlim=NULL,ylim=NULL,
                     layout=NULL,mfcol=FALSE,
                     probs=c(1:99)/100,
                     data_lwd=2,fit_lwd=1,
                     q_points=c(.1,.3,.5,.7,.9),
                     qp_cex=1,pqp_cex=.5,lpos="right", main = "")
{
  if (!is.null(stat) & is.null(factors)) factors <- NA
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
    if("subjects" %in% factors){
      fnams <- names(dat)[!(names(dat) %in% c("trials","R","rt"))]
    } else{
      fnams <- names(dat)[!(names(dat) %in% c("subjects", "trials","R","rt"))]
    }
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
    if (mfcol) par(mfcol=layout) else par(mfrow=layout)

  if (all(is.na(data$rt))) stop("Use plot_fit_choice if no rt data")

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





#' Plot within-chain correlations
#'
#' Plots within-chain parameter correlations (upper triangle) and corresponding
#' scatterplots (lower triangle) to visualize parameter sloppiness.
#'
#' If ``selection = alpha`` the parameter chains are concatenated across participants,
#' (after standardizing if ``scale_subjects = TRUE``) and then correlated.
#'
#' @param emc An emc object
#' @param selection A Character string. Indicates which parameter type to
#' plot (`alpha`, `mu`, `variance`, `covariance`, `correlation`).
#' @param scale_subjects Boolean. To standardize each participant with ``selection = "alpha"``,
#'  by subtracting the mean and divding by the standard deviation. This ensures the plot has every participant on the same scale.
#' @param do_plot Boolean. Whether to plot the pairs plot, if ``FALSE``, only the correlations
#' are returned.
#' @param N Integer for maximum number of iterations used (defaults to 500).
#' If number of samples in stage or selection exceeds N, a random subset will be taken of size N
#' @param ... Optional arguments that can be passed to `get_pars`
#' @return Invisibly returns a matrix with the correlations between the parameters.
#' @examples \donttest{
#' # Plot the sloppiness for the individual-level subjects
#' pairs_posterior(samples_LNR, selection = "alpha")
#'
#' # We can also choose group-level parameters and subsets of the parameter space
#' pairs_posterior(samples_LNR, use_par = c("m", "t0"), selection = "sigma2")
#' }
#' @export

pairs_posterior <- function(emc, selection="alpha", scale_subjects=TRUE,
                            do_plot=TRUE,N=500, ...)
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
  dots <- add_defaults(list(...), stage = "sample")
  use_par <- dots$use_par
  N <- min(chain_n(emc)[,dots$stage], N)
  pmat <- do.call(parameters, c(list(emc,selection=selection, N=N), fix_dots(dots, get_pars)))
  if (selection=="alpha") {
    if (length(levels(pmat$subjects))>1 && scale_subjects)
      for (i in names(pmat)[-1]) pmat[,i] <- do_scale(pmat[,c("subjects",i)])
    pmat <- pmat[,-1]
  }
  if (do_plot) suppressWarnings(
    pairs(pmat,diag.panel = panel.hist,upper.panel = panel.cor))
  rs <- cor(pmat)
  r.names <- outer(dimnames(rs)[[1]],dimnames(rs)[[2]],paste,sep="~")[upper.tri(rs)]
  rs <- rs[upper.tri(rs)]
  names(rs) <- r.names
  invisible(rs)
}

#' Likelihood profile plots
#'
#' Creates likelihood profile plots from a design and the experimental data by
#' varying one model parameter while holding all others constant.
#'
#' @param data A dataframe. Experimental data used, needed for the design mapping
#' @param design A design list. Created using ``design``.
#' @param p_vector Named vector of parameter values (typically created with ``sampled_p_vector(design)``)
#' @param range Numeric. The max and min will be p_vector + range/2 and p_vector - range/2, unless specified in p_min or p_max.
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param p_min Named vector. If specified will instead use these values for minimum range of the selected parameters.
#' @param p_max Named vector. If specified will instead use these values for maximum range of the selected parameters.
#' @param use_par Character vector. If specified will only plot the profiles for the specified parameters.
#' @param n_point Integer. Number of evenly spaced points at which to calculate likelihood
#' @param n_cores Number of likelihood points evenly spaced between the minimum and maximum likelihood range.
#' @param true_plot_args A list. Optional additional arguments that can be passed to plot.default for the plotting of the true vertical line.
#' @param round Integer. To how many digits will the output be rounded.
#' @param ... Optional additional arguments that can be passed to plot.default.
#' @return Vector with highest likelihood point, input and mismatch between true and highest point
#' @examples \donttest{
#' # First create a design
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                       formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                       constants=c(s=log(1)))
#' # Then create a p_vector:
#' p_vector=c(v_Sleft=-2,v_Sright=2,a=log(.95),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'           t0=log(.25),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' # Make a profile plot for some parameters. Specifying a custom range for t0.
#' profile_plot(p_vector = p_vector, p_min = c(t0 = -1.35),
#'              p_max = c(t0 = -1.45), use_par = c("a", "t0", "SZ"),
#'              data = forstmann, design = design_DDMaE, n_point = 10)
#' }
#' @export

profile_plot <- function(data, design, p_vector, range = .5, layout = NA,
                         p_min = NULL,p_max = NULL, use_par = NULL,
                         n_point=100,n_cores=1, round = 3,
                         true_plot_args = list(),
                         ...)

{
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  dots <- list(...)
  lfun <- function(i,x,p_vector,pname,dadm) {
    p_vector[pname] <- x[i]
    attr(dadm,"model")()$log_likelihood(p_vector,dadm)
  }
  if(!identical(names(p_min), names(p_max))) stop("p_min and p_max should be specified for the same parameters")
  if(!is.null(names(p_min)) & length(p_min) == length(use_par)) names(p_min) <- use_par
  if(!is.null(names(p_max)) & length(p_max) == length(use_par)) names(p_max) <- use_par
  if(is.null(use_par)) use_par <- names(p_vector)
  if(any(is.na(layout))){
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(use_par),
                              nplots = 1))
  } else{par(mfrow=layout)}
  if(is.null(dots$dadm)){
    dadm <- design_model(data, design, verbose = FALSE)
  } else{
    dadm <- dots$dadm
  }
  out <- data.frame(true = rep(NA, length(use_par)), max = rep(NA, length(use_par)), miss = rep(NA, length(use_par)))
  rownames(out) <- use_par
  for(p in 1:length(p_vector)){
    cur_name <- names(p_vector)[p]
    if(cur_name %in% use_par){
      cur_par <- p_vector[p]
      pmax_cur <- cur_par + range/2
      pmin_cur <- cur_par - range/2
      if(!is.null(p_min)){
        if(!is.na(p_min[cur_name])){
          pmin_cur <- p_min[cur_name]
        }
      }
      if(!is.null(p_max)){
        if(!is.na(p_max[cur_name])){
          pmax_cur <- p_max[cur_name]
        }
      }
      x <- seq(pmin_cur,pmax_cur,length.out=n_point)
      x <- c(x, cur_par)
      x <- unique(sort(x))
      ll <- unlist(mclapply(1:length(x),lfun,dadm=dadm,x=x,p_vector=p_vector,pname=cur_name,mc.cores = n_cores))
      do.call(plot, c(list(x,ll), fix_dots_plot(add_defaults(dots, type="l",xlab=cur_name,ylab="LL"))))
      do.call(abline, c(list(v=cur_par), fix_dots_plot(add_defaults(true_plot_args, lty = 2))))
      out[cur_name,] <- c(p_vector[cur_name], x[which.max(ll)], p_vector[cur_name] - x[which.max(ll)])
    }
  }
  return(round(out, 3))
}

#' Plots density for parameters
#'
#' Plots the posterior and prior density for selected parameters of a model.
#' Full range of samples manipulations described in `get_pars`.
#'
#' @param emc An emc object
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param selection A Character string. Indicates which parameter type to use (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param show_chains Boolean (defaults to `FALSE`) plots a separate density for each chain.
#' @param plot_prior Boolean. If ``TRUE`` will overlay prior density in the plot (default in red)
#' @param N Integer. How many prior samples to draw
#' @param use_prior_lim Boolean. If `TRUE` will use xlimits based on prior density, otherwise based on posterior density.
#' @param lpos Character. Where to plot the contraction statistic.
#' @param true_pars A vector or emc object. Can be used to visualize recovery.
#' If a vector will plot a vertical line for each parameter at the appropriate place.
#' If an emc object will plot the densities of the object as well, assumed to be the data-generating posteriors.
#' @param all_subjects Boolean. Will plot the densities of all (selected) subjects overlaid with the group-level distribution
#' @param prior_plot_args A list. Optional additional arguments to be passed to plot.default for the plotting of the prior density (see `par()`)
#' @param true_plot_args A list. Optional additional arguments to be passed to plot.default for the plotting of the true parameters (see `par()`)
#' @param ... Optional arguments that can be passed to `get_pars`, `density`, or `plot.default` (see `par()`)
#'
#' @return An invisible return of the contraction statistics for the selected parameter type
#' @export
#'
#' @examples
#' # Full range of possibilities described in get_pars
#' plot_pars(samples_LNR)
#' # Or plot all subjects
#' plot_pars(samples_LNR, all_subjects = TRUE, col = 'purple')
#' # Or plot recovery
#' true_emc <- samples_LNR # This would normally be the data-generating samples
#' plot_pars(samples_LNR, true_pars = true_emc, true_plot_args = list(col = 'blue'), adjust = 2)
plot_pars <- function(emc,layout=NA, selection="mu", show_chains = FALSE, plot_prior = TRUE, N = 1e4,
                      use_prior_lim = !all_subjects, lpos = "topright", true_pars = NULL, all_subjects = FALSE,
                      prior_plot_args = list(), true_plot_args = list(), ...)
{
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  dots <- list(...)
  type <- attr(emc[[1]], "variant_funs")$type
  if(length(dots$subject) == 1 || emc[[1]]$n_subjects == 1) dots$by_subject <- TRUE
  if(type == "single" & !(selection %in% c("LL", "alpha"))) selection <- "alpha"
  if(all_subjects){
    dots$by_subject <- TRUE; show_chains <- TRUE; selection <- "alpha"
    true_pars <- NULL
    if(is.null(dots$lwd)) dots$lwd <- .3
  }

  MCMC_samples <- do.call(get_pars, c(list(emc, selection = selection), fix_dots(dots, get_pars)))
  if(all_subjects) {
    MCMC_samples <- list(lapply(MCMC_samples, function(x) do.call(rbind, x)))
    names(MCMC_samples) <- "subjects"
  }
  psamples <-  get_objects(sampler = emc, design = attr(emc,"design_list"),
                           type = type, sample_prior = T,
                           selection = selection, N = N,
                           prior = emc[[1]]$prior)
  pMCMC_samples <- do.call(get_pars, c(list(psamples, selection = selection, type = type),
                                       fix_dots(dots, get_pars, exclude = c("thin", "filter", "chain", "subject"))))
  if(length(pMCMC_samples) != length(MCMC_samples)) pMCMC_samples <- rep(pMCMC_samples, length(MCMC_samples))

  true_MCMC_samples <- NULL
  if(!is.null(true_pars)){
    if(!is(true_pars, "emc")){
      true_pars <- do.call(get_pars, c(list(emc, selection = selection, type = type, true_pars = true_pars),
                                       fix_dots(dots, get_pars, exclude = c("thin", "filter", "chain"))))
    } else{
      true_MCMC_samples <- do.call(get_pars, c(list(true_pars, selection = selection), fix_dots(dots, get_pars)))
      true_pars <- NULL
    }
  }


  n_objects <- length(MCMC_samples)
  contraction_list <- list()
  for(i in 1:n_objects){
    cur_mcmc <- MCMC_samples[[i]]
    merged <- do.call(rbind, cur_mcmc) # Need for contraction calculation
    if(!show_chains) {
      cur_mcmc <- list(merged)
      if(!is.null(true_MCMC_samples)) true_MCMC_samples[[i]] <- list(do.call(rbind, true_MCMC_samples[[i]]))
    }
    n_chains <- length(cur_mcmc)
    n_pars <- ncol(merged)
    cur_contraction <- setNames(numeric(n_pars), colnames(merged))
    if(any(is.na(layout))){
      par(mfrow = coda_setmfrow(Nchains = n_chains, Nparms = n_pars, nplots = 1))
    } else{
      par(mfrow=layout)
    }
    for(l in 1:n_pars){
      denses <- lapply(cur_mcmc, function(x) do.call(density, c(list(x[,l]), fix_dots(dots, density.default, consider_dots = FALSE))))
      xlim <- range(c(sapply(cur_mcmc, function(x) return(range(x[,l]))), true_pars[[i]][[1]][,l]))
      ylim <- range(sapply(denses, function(x) return(range(x$y))))
      if(ncol(pMCMC_samples[[i]][[1]]) != n_pars){
        p_idx <- 1
      } else{p_idx <- l}
      if(plot_prior){
        if(selection != "alpha"){
          p_curr <- robust_hist(pMCMC_samples[[i]][[1]][,p_idx], do_plot = F)
        } else{ p_curr <-  pMCMC_samples[[i]][[1]][,p_idx]}
        pdenses <- do.call(density, c(list(p_curr), fix_dots(dots, density.default, consider_dots = FALSE)))
        if(use_prior_lim){
          xlim <- range(c(xlim, p_curr))
          ylim <- range(c(ylim, pdenses$y))
        }
      }
      if(!is.null(true_MCMC_samples)){
        true_denses <- lapply(true_MCMC_samples[[i]], function(x) do.call(density, c(list(x[,l]), fix_dots(dots, density.default, consider_dots = FALSE))))
        xlim <- range(c(sapply(true_MCMC_samples[[i]], function(x) return(range(x[,l]))), xlim))
        ylim <- range(c(sapply(true_denses, function(x) return(range(x$y))), ylim))
      }
      for(k in 1:n_chains){
        x_name <- ifelse(n_objects == 1, names(MCMC_samples)[i], paste0(selection, ": ", names(MCMC_samples)[i]))
        if(k == 1){
          cur_dots <- add_defaults(dots, xlim = xlim, ylim = ylim, ylab = "Density",
                                   xlab = x_name, main = colnames(cur_mcmc[[k]])[l])
          do.call(plot, c(list(denses[[k]]), fix_dots_plot(cur_dots)))
        } else{
          do.call(lines, c(list(denses[[k]]), fix_dots_plot(cur_dots)))
        }
        if(!is.null(true_MCMC_samples)){
          cur_true_plot_args <- add_defaults(true_plot_args, col = "darkgreen")
          do.call(lines, c(list(true_denses[[k]]), fix_dots_plot(cur_true_plot_args)))
        }
      }
      if(!is.null(true_pars)){
        cur_true_plot_args <- add_defaults(true_plot_args, col = "darkgreen", lwd = 1.5, lty = 2)
        do.call(abline, c(list(v = true_pars[[i]][[1]][,l]), fix_dots_plot(cur_true_plot_args)))
      }
      cur_contraction[l] <- 1-(var(merged[,l])/var(pMCMC_samples[[i]][[1]][,p_idx]))
      if(plot_prior){
        cur_prior_plot_args <- add_defaults(prior_plot_args, col = "red")
        do.call(lines, c(list(pdenses), fix_dots_plot(cur_prior_plot_args)))
        if(!all_subjects) legend(lpos,legend=round(cur_contraction[l],3),bty="n",title="Contraction")
      }
    }
    contraction_list[[x_name]] <- cur_contraction
  }
  return(invisible(contraction_list))
}

get_recovery_stats <- function(MCMC, true_MCMC, true_pars, CI){
  quantiles = c(.5 - CI/2, .5, .5 + CI/2)
  if(!is.null(true_MCMC)){
    quants_true <- get_posterior_quantiles(true_MCMC, probs = quantiles)
    true_pars <- quants_true[,"50%"]
  }
  if(is.list(true_pars)) true_pars <- unlist(true_pars)
  quants <- get_posterior_quantiles(MCMC, probs = quantiles)
  pearson <- cor(true_pars,quants[,"50%"],method="pearson")
  spearman <- cor(true_pars,quants[,"50%"],method="spearman")
  rmse <- sqrt(mean((true_pars - quants[,"50%"])^2))
  coverage = 100*mean((quants[,3] > true_pars)  & (quants[,1] < true_pars))
  if(is.null(true_MCMC)){
    quants_true <- matrix(true_pars)
    colnames(quants_true) <- "50%"
  }
  return(list(true = quants_true, recovered = quants,
              pearson = pearson, spearman = spearman, rmse = rmse, coverage = coverage))
}

make_recov_summary <- function(stats){
  out <- list()
  rec <- stats$recovered
  true <- stats$true
  if(ncol(true) == 1){
    miss <- rec[,"50%"] - true
    out$quantiles <- cbind(rec, true, miss)
    colnames(out$quantiles)[4:5] <- c("true", "miss")
  } else{
    out$recovered <- rec
    out$true <- true
    out$miss <- rec-true
  }
  out$stats <- setNames(c(stats$pearson, stats$spearman, stats$rmse, stats$coverage),
                        c("pearson", "spearman", "rmse", "coverage"))
  return(out)
}


#' Plot MCMC
#'
#' Uses the coda plot functions that are applied per chain
#'
#' @param emc An emc object
#' @param selection A Character string. Indicates which parameter type to plot (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param fun A plot function that takes a vector/mcmc object as input, e.g. cumuplot, acf
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param plot_type type argument passed on to coda fun.
#' @param chain Integer, which chain to include, if more than 1 will make separate plots per chain.
#' @param ... Optional arguments that can be passed to `get_pars`,
#' the chosen coda plot function, or `plot.default` (see `par()`)
#'
#' @return A coda plot
plot_mcmc <- function(emc, selection = "mu", fun = 'cumuplot', layout=NA, chain = 1,
                      plot_type = NULL, ...)
{
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  dots <- list(...)
  if(length(dots$subject) == 1 || emc[[1]]$n_subjects == 1) dots$by_subject <- TRUE
  MCMC_samples <- do.call(get_pars, c(list(emc, selection = selection, chain = chain),
                                      fix_dots(dots, get_pars)))

  for(i in 1:length(MCMC_samples)){
    for(k in 1:length(MCMC_samples[[i]])){
      if(any(is.na(layout))){
        par(mfrow = coda_setmfrow(Nchains = length(MCMC_samples[[1]]), Nparms = max(ncol(MCMC_samples[[1]][[1]]), length(MCMC_samples)),
                                  nplots = 1))
      } else{par(mfrow=layout)}
      for(j in 1:ncol(MCMC_samples[[i]][[1]])){
        cur_dots <- dots
        x_name <- ifelse(length(MCMC_samples) == 1, names(MCMC_samples)[i], paste0(selection, ": ", names(MCMC_samples)[i]))
        cur_dots <- add_defaults(cur_dots, ylab = paste0("Chain : ", chain[k]),
                                 xlab = x_name, main = colnames(MCMC_samples[[i]][[1]])[j],
                                 type = plot_type, auto.layout = F, ask = F)
        # do.call(fun, c(list(MCMC_samples[[i]][[k]][,j], ask = F, auto.layout =F), fix_dots(cur_dots, get(fun))))
        do.call(fun, c(list(MCMC_samples[[i]][[k]][,j]), fix_dots_plot(cur_dots),
                       fix_dots(cur_dots, get(fun), consider_dots = F, exclude =
                                  c(names(par()), names(formals(arrows)), names(formals(plot.default))))))
      }


    }
  }
}

#' Plot MCMC.list
#'
#' Uses the coda plot functions that are applied across chain
#'
#' @param emc An emc object
#' @param selection A Character string. Indicates which parameter type to plot (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param fun A coda plot function choice from
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param ... Optional arguments that can be passed to `get_pars`,
#' the chosen coda plot function, or `plot.default` (see `par()`)
#'
#' @return A coda plot
plot_mcmc_list <- function(emc, selection = "mu", fun = 'traceplot', layout=NA, ...)
{
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  dots <- list(...)
  if(length(dots$subject) == 1 || emc[[1]]$n_subjects == 1) dots$by_subject <- TRUE
  MCMC_samples <- do.call(get_pars, c(list(emc, selection = selection), fix_dots(dots, get_pars)))

  for(i in 1:length(MCMC_samples)){
    if(any(is.na(layout))){
      par(mfrow = coda_setmfrow(Nchains = length(MCMC_samples[[1]]), Nparms = max(ncol(MCMC_samples[[1]][[1]]), length(MCMC_samples)),
                                nplots = 1))
    } else{par(mfrow=layout)}
    cur_dots <- dots
    x_name <- ifelse(length(MCMC_samples) == 1, names(MCMC_samples)[i], paste0(selection, ": ", names(MCMC_samples)[i]))
    cur_dots <- add_defaults(cur_dots, xlab = x_name)
    do.call(fun, c(list(MCMC_samples[[i]]), fix_dots_plot(cur_dots),
                        fix_dots(cur_dots, get(fun), consider_dots = F, exclude =
                                   c(names(par()), names(formals(arrows)), names(formals(plot.default))))))
  }
}






coda_setmfrow <- function (Nchains = 1, Nparms = 1, nplots = 1, sepplot = FALSE)
{
  mfrow <- if (sepplot && Nchains > 1 && nplots == 1) {
    if (Nchains == 2) {
      switch(min(Nparms, 5), c(1, 2), c(2, 2), c(3, 2),
             c(4, 2), c(3, 2))
    }
    else if (Nchains == 3) {
      switch(min(Nparms, 5), c(2, 2), c(2, 3), c(3, 3),
             c(2, 3), c(3, 3))
    }
    else if (Nchains == 4) {
      if (Nparms == 1)
        c(2, 2)
      else c(4, 2)
    }
    else if (any(Nchains == c(5, 6, 10, 11, 12)))
      c(3, 2)
    else if (any(Nchains == c(7, 8, 9)) || Nchains >= 13)
      c(3, 3)
  }
  else {
    if (nplots == 1) {
      mfrow <- switch(min(Nparms, 13), c(1, 1), c(1, 2),
                      c(2, 2), c(2, 2), c(3, 2), c(3, 2), c(3, 3),
                      c(3, 3), c(3, 3), c(3, 2), c(3, 2), c(3, 2),
                      c(3, 3))
    }
    else {
      mfrow <- switch(min(Nparms, 13), c(1, 2), c(2, 2),
                      c(3, 2), c(4, 2), c(3, 2), c(3, 2), c(4, 2),
                      c(4, 2), c(4, 2), c(3, 2), c(3, 2), c(3, 2),
                      c(4, 2))
    }
  }
  return(mfrow)
}
