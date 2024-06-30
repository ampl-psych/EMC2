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
#' @param layout A vector specifying the layout as in `par(mfrow = layout)`.
#' The default `NULL` uses the current layout.
#' @param xlim x-axis limit for all cells (default NULL will scale each plot such that the xlimits encompass the densities).
#' @param bw number or string bandwidth for density (defaults to `nrd0`). See ``?density()``.
#' @param adjust Numeric. Density function bandwidth adjustment parameter. See ``?density()``.
#' @param correct_fun If specified, the accuracy for each subject is calculated, using the supplied function and
#' an accuracy vector for each subject is returned invisibly.
#' @param rt legend function position character string for the mean response time (defaults to `top`)
#' @param accuracy legend function position character string for accuracy (defaults to `topright`)
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
                                   layout=NULL,
                                   xlim=NULL,bw = "nrd0",adjust=1,
                                   correct_fun=NULL,rt="top",accuracy="topright")
{
  if (!is.null(subject)) {
    snams <- levels(data$subjects)
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% snams)) stop("Subject not present\n")
    dat <- data[data$subjects==subject,]
    dat$subjects <- droplevels(dat$subjects)
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

robust_hist <- function(ps, breaks = 50, cutoff = 0.0015, main = "",
                        cex.lab = 1, cex.main = 1.5, prob = TRUE,
                        xlab = "", do_plot = TRUE, ...){
  ps <- ps[abs(ps) < 1000]
  for(i in 1:log10(length(ps))){
    cuts <- cut(ps, breaks = breaks)
    tab_cuts <- table(cuts)
    good_cuts <- tab_cuts/length(ps) > cutoff
    ps <- ps[cuts %in% names(good_cuts)[good_cuts]]
  }
  if(do_plot){
    hist(ps, breaks = breaks, main = main, cex.lab = cex.lab, cex.main = cex.main,
         prob = prob, xlab = xlab, ...)
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
#' @param pp Posterior predictives created by `post_predict()`
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
#' @param pp A data frame. Posterior predictives created by ``post_predict()``
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
#' @examples \dontrun{
#' # First generate posterior predictives based on an emc object run with run_emc
#' pp <- post_predict(samplers_LNR, n_cores = 4)
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
                     qp_cex=1,pqp_cex=.5,lpos="topleft", main = "")
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
#' (after standardizing if ``scale.subjects = TRUE``) and then correlated.
#'
#' @param emc An emc object
#' @param filter A character. Specifies which sampling stage should be plotted.
#' Defaults to the sampling stage `sample`.
#' @param thin An integer. Only iterations that are a multiple of `thin` are kept.
#' @param subfilter Integer or numeric vector. If an integer is supplied, iterations
#' up until that integer within the sampling stage `filter` are kept. If a vector is supplied, the iterations
#' within the range are kept.
#' @param selection A Character string. Indicates which parameter type to plot (`alpha`, `mu`, `variance`, `covariance`, `correlation`).
#' `LL` will plot the log-likelihood chains.
#' @param mapped Boolean. If ``TRUE``, plots the parameters as mapped back to the factor levels of the design,
#' if `FALSE`, the sampled parameters are plotted
#' @param scale.subjects Boolean. To standardize each participant with ``selection = "alpha"``,
#'  by subtracting the mean and divding by the standard deviation. This ensures the plot has every participant on the same scale.
#' @param use_par Character vector of names of parameters to plot (default ``NULL`` plots everything)
#' @param do_plot Boolean. Whether to plot the pairs plot, if ``FALSE``, only the correlations
#' are returned.
#' @param maxp Integer for maximum number of iterations used (defaults to 500).
#' If number of samples in stage or selection exceeds maxp, a random subset will be taken of size maxp.
#' @return Invisibly returns a matrix with the correlations between the parameters.
#' @examples \dontrun{
#' # Plot the sloppiness for the individual-level subjects
#' pairs_posterior(samplers_LNR, selection = "alpha")
#'
#' # We can also choose group-level parameters and subsets of the parameter space
#' pairs_posterior(samplers_LNR, use = c("v", "B", "t0"), selection = "sigma2")
#' }
#' @export

pairs_posterior <- function(emc,filter="sample",thin=1,subfilter=0,mapped=FALSE,
                            selection="mu",
                            scale.subjects=TRUE,use_par=NULL,do_plot=TRUE,maxp=500, ...)
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

  pmat <- parameters_data_frame(emc,filter=filter,thin=thin,subfilter=subfilter,
                                mapped=mapped,selection=selection, ...)
  if (!is.null(use_par)) {
    if (is.numeric(use_par)) {
      if (any(use_par<1) || any(use_par>dim(pmat))) stop("use_par outside parameter range")
      use_par <- names(pmat)[use_par]
    }
    if (!all(use_par %in% names(pmat))) stop("use_par has a name not in parameters")
    if (length(use_par)==1) stop("must select more than one parameter")
    pmat <- pmat[,use_par]
  }
  if (selection=="alpha" & is.null(use_par)) {
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

#' Likelihood profile plots
#'
#' Creates likelihood profile plots from a design and the experimental data by
#' varying one model parameter while holding all others constant.
#'
#' @param pname String. Name of parameter to profile
#' @param p_vector Named vector of parameter values (typically created with ``sampled_p_vector(design)``)
#' @param p_min Numeric. Minimum value of profile range
#' @param p_max Numeric. Maximum of profile range
#' @param n_point Integer. Number of evenly spaced points at which to calculate likelihood
#' @param main Plot title
#' @param cores Number of likelihood points to calculate in parallel
#' @param data A dataframe. Experimental data used, needed for the design mapping
#' @param design A design list. Created using ``make_design``.
#'
#' @return Vector with highest likelihood point, input and mismatch between true and highest point
#' @examples
#' # First create a design
#' design_DDMaE <- make_design(data = forstmann,model=DDM,
#'                 formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                 constants=c(s=log(1)))
#' # Then create a p_vector:
#' p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'           t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' # Make a profile plot for the boundary separation parameter
#' profile_plot(pname = "a", p_vector = p_vector, p_min = -1, p_max = 1,
#'             data = forstmann, design = design_DDMaE, n_point = 10)
#'
#' @export

profile_plot <- function(data, design, p_vector, layout = NA,
                         p_min = NULL,p_max = NULL,
                         n_point=100,main="",n_cores=1,
                         true_plot_args = list(),
                         ...)

{
  dots <- list(...)
  lfun <- function(i,x,p_vector,pname,dadm) {
    p_vector[pname] <- x[i]
    attr(dadm,"model")()$log_likelihood(p_vector,dadm)
  }
  if(!identical(names(p_min), names(p_max))) stop("p_min and p_max should be specified for the same parameters")
  pars_to_use <- p_vector
  if(!is.null(names(p_min))){
    pars_to_use <- pars_to_use[names(pars_to_use) %in% names(p_min)]
  } else{
    if(!is.null(p_min)){
      if(length(p_vector) != length(p_min) & length(p_min) != 1) stop("Without names you can either specify p_min/p_max for one parameter or for all parameters")
    }
  }
  if(any(is.na(layout))){
    par(mfrow = coda:::set.mfrow(Nchains = 1, Nparms = length(pars_to_use),
                                 nplots = 1))
  } else{par(mfrow=layout)}
  if(is.null(dots$dadm)) dadm <- design_model(data, design, verbose = FALSE)
  for(p in 1:length(pars_to_use)){
    cur_name <- names(pars_to_use)[p]
    cur_par <- pars_to_use[p]
    if(!is.null(p_min) & !is.null(p_max)){
      pmin_cur <- p_min[p]
      pmax_cur <- p_max[p]
    } else{
      pmin_cur <- cur_par - .5
      pmax_cur <- cur_par + .5
    }
    x <- seq(pmin_cur,pmax_cur,length.out=n_point)
    ll <- unlist(mclapply(1:n_point,lfun,dadm=dadm,x=x,p_vector=p_vector,pname=cur_name,mc.cores = n_cores))
    do.call(plot, c(list(x,ll), fix_dots_plot(add_defaults(dots, type="l",xlab=cur_name,ylab="LL",main=main))))
    do.call(abline, c(list(v=cur_par), fix_dots_plot(add_defaults(true_plot_args, lty = 2))))
  }
  # c(true=p_vector[pname],max=x[which.max(ll)],miss=p_vector[pname]-x[which.max(ll)])
}

#' Plot MCMC chains
#'
#' Plots the trace plots of the MCMC chains on top of each other. Visualizes convergence
#' and chain stability.
#'
#' @param emc An emc object
#' @param selection A Character string. Indicates which parameter type to plot (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param plot_acf Boolean. If `FALSE` will make trace plots. If `TRUE` will plot the autocorrelation for a chain.
#' By default plots acf for the first chain, but can be changed, by setting i.e. chain = 2; see `?as_mcmc_new`.
#' @param ... A list of optional arguments that can be passed to `as_mcmc_new` or `plot.default` (see `par()`)
#'
#' @return A trace/acf plot of the selected MCMC chains
#' @export
#'
#' @examples
#' # Full range of possibilities described in as_mcmc_new.
#' plot_chains(samplers_LNR)
#' # Or trace plots for the second subject:
#' plot_chains(samplers_LNR, subject = 2, selection = "alpha", plot_acf = TRUE)
#'
#' # Can also plot the trace of for example the group-level correlation:
#' plot_chains(emc, selection = "correlation", plot_acf = TRUE, colors = c("green", "purple", "orange"), lwd = 2)

plot_chains <- function(emc, selection = "mu", layout=NA, plot_acf=FALSE, ...)
{
  dots <- list(...)
  if(plot_acf) dots <- add_defaults(dots, chain = 1)
  MCMC_samples <- do.call(as_mcmc_new, c(list(emc, selection = selection), fix_dots(dots, as_mcmc_new)))

  for(i in 1:length(MCMC_samples)){
    if(any(is.na(layout))){
      par(mfrow = coda:::set.mfrow(Nchains = length(MCMC_samples[[1]]), Nparms = max(ncol(MCMC_samples[[1]][[1]]), length(MCMC_samples)),
                                   nplots = 1))
    } else{par(mfrow=layout)}
    if (plot_acf){
      for(k in 1:length(MCMC_samples[[i]])){
        for (j in colnames(MCMC_samples[[i]][[1]])){
          cur_dots <- dots
          x_name <- ifelse(length(MCMC_samples) == 1, names(MCMC_samples)[i], paste0(selection, ": ", names(MCMC_samples)[i]))
          cur_dots <- add_defaults(cur_dots, ylab = "acf", main = paste0("Chain ",dots$chain[k],": ",j), xlab = x_name)
          do.call(acf, c(list(MCMC_samples[[i]][[k]][,j]), fix_dots_plot(cur_dots)))
        }
      }

    } else {
      cur_dots <- dots
      x_name <- ifelse(length(MCMC_samples) == 1, names(MCMC_samples)[i], paste0(selection, ": ", names(MCMC_samples)[i]))
      cur_dots <- add_defaults(cur_dots, xlab = x_name)
      do.call(plot, c(list(MCMC_samples[[i]], auto.layout = FALSE, density = FALSE, ask = FALSE,smooth = FALSE),
                      fix_dots_plot(cur_dots)))
    }
  }
}

#' Plots density for parameters
#'
#' Plots the posterior and prior density for selected parameters of a model.
#'
#' @param emc An emc object
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param selection A Character string. Indicates which parameter type to plot (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
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
#' @param ... A list of optional arguments that can be passed to `as_mcmc_new`, `density`, or `plot.default` (see `par()`)
#'
#' @return An invisible return of the contraction statistics for the selected parameter type
#' @export
#'
#' @examples
#' # Full range of possibilities described in as_mcmc_new
#' plot_pars(samplers_LNR)
#' # Or plot all subjects
#' plot_pars(samplers_LNR, all_subjects = TRUE, col = 'purple')
#' # Or plot recovery
#' true_emc <- samplers_LNR # This would normally be the data-generating samples
#' plot_pars(samplers_LNR, true_pars = true_emc, true_plot_args = list(col = 'blue'), adjust = 2)
plot_pars <- function(emc,layout=NA, selection="mu", show_chains = FALSE, plot_prior = TRUE, N = 1e4,
                          use_prior_lim = !all_subjects, lpos = "topright", true_pars = NULL, all_subjects = FALSE,
                          prior_plot_args = list(), true_plot_args = list(), ...)
{
  dots <- list(...)
  type <- attr(emc[[1]], "variant_funs")$type
  if(length(dots$subject) == 1) dots$by_subject <- TRUE
  if(all_subjects){
    dots$by_subject <- TRUE; show_chains <- TRUE; selection <- "alpha"
    true_pars <- NULL
    if(is.null(dots$lwd)) dots$lwd <- .3
  }

  MCMC_samples <- do.call(as_mcmc_new, c(list(emc, selection = selection), fix_dots(dots, as_mcmc_new)))
  if(all_subjects) {
    MCMC_samples <- list(lapply(MCMC_samples, function(x) do.call(rbind, x)))
    names(MCMC_samples) <- "test"
  }
  if(plot_prior){
    psamples <-  get_objects(sampler = emc, design = attr(emc,"design_list")[[1]],
                             type = type, sample_prior = T,
                             selection = selection, N = N)
    pMCMC_samples <- do.call(as_mcmc_new, c(list(psamples, selection = selection, type = type),
                                           fix_dots(dots, as_mcmc_new, exclude = c("thin", "subfilter", "chain", "subject"))))
    if(length(pMCMC_samples) != length(MCMC_samples)) pMCMC_samples <- rep(pMCMC_samples, length(MCMC_samples))
  }

  true_MCMC_samples <- NULL
  if(!is.null(true_pars)){
    if(!is(true_pars, "emc")){
      if(selection == "sigma2" & !is.matrix(true_pars)) true_pars <- diag(true_pars)
      true_pars <- do.call(as_mcmc_new, c(list(emc, selection = selection, type = type, true_pars = true_pars),
                                          fix_dots(dots, as_mcmc_new, exclude = c("thin", "subfilter", "chain"))))
    } else{
      true_MCMC_samples <- do.call(as_mcmc_new, c(list(true_pars, selection = selection), fix_dots(dots, as_mcmc_new)))
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
      par(mfrow = coda:::set.mfrow(Nchains = n_chains, Nparms = n_pars, nplots = 1))
    } else{
      par(mfrow=layout)
    }
    for(l in 1:n_pars){
      denses <- lapply(cur_mcmc, function(x) do.call(density, c(list(x[,l]), fix_dots(dots, density.default, consider_dots = FALSE))))
      xlim <- range(c(sapply(cur_mcmc, function(x) return(range(x[,l]))), true_pars[[i]][[1]][,l]))
      ylim <- range(sapply(denses, function(x) return(range(x$y))))
      if(plot_prior){
        if(ncol(pMCMC_samples[[i]][[1]]) != n_pars){
          p_idx <- 1
        } else{p_idx <- l}
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


#' Convergence checks for an emc object
#'
#' Runs a series of convergence checks, prints statistics to the console, and
#' makes traceplots of the worst converged parameter per selection.
#'
#' Note that the `Rhat` is calculated by doubling the number of chains by
#' first splitting chains into first and second half, so it also a test of
#' stationarity.
#'
#' Efficiency of sampling is indicated by the effective
#' sample size (ESS) (from the `coda` R package).
#'
#' @param emc An emc object
#' @param selection A Character vector. Indicates which parameter types to check (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param digits Integer. How many digits to round the ESS and Rhat to in the plots
#' @param plot_worst Boolean. If `TRUE` also plots the chain plots for the worst parameter
#' @param ... A list of optional arguments that can be passed to `as_mcmc_new` or `plot.default` (see `par()`)
#'
#' @return a list with the statistics for the worst converged parameter per selection
#' @export
#' @examples
#' check_run(samplers_LNR)
#'
check_run <- function(emc, selection = c('mu', 'sigma2', 'alpha'), digits = 3,
                          plot_worst = TRUE, ...){
  dots <- list(...)
  out_list <- list()
  cat("Iterations:\n")
  print(chain_n(emc))
  if(attr(emc[[1]], "variant_funs")$type == "single") selection <- "alpha"
  if(plot_worst){
    mfrow <- coda:::set.mfrow(Nchains = length(emc), Nparms = length(selection),nplots = 1)
    oldpar <- par(mfrow = mfrow)
  }
  for(select in selection){
    dots$flatten <- ifelse(select == "alpha", FALSE, TRUE)
    dots$by_subject <- TRUE
    ESS <- do.call(es_summary_new, c(list(emc, selection = select, stat= NULL), fix_dots(dots, es_summary_new)))
    gds <- do.call(gd_summary_new, c(list(emc, selection = select, stat= NULL), fix_dots(dots, gd_summary_new)))
    out <- list()
    max_gd <- -Inf
    for(name in names(ESS)){
      combined <- rbind(round(gds[[name]], digits), round(ESS[[name]]))
      rownames(combined) <- c("Rhat", "ESS")
      out[[name]] <- combined
      if(max(gds[[name]]) > max_gd){
        max_gd <- max(gds[[name]])
        cur_max <- name
        max_par <- names(gds[[name]])[which.max(gds[[name]])]
      }
    }
    if(length(ESS) > 1){
      cat("\n", paste0(select, " highest Rhat : ", cur_max), "\n")
    } else{
      cat("\n", cur_max, "\n")
    }
    print(out[[cur_max]])
    if(plot_worst){
      cur_dots <- dots
      if(select == "alpha"){
        cur_dots$subject <- cur_max
        cur_dots$use_par <- max_par
        cur_dots$by_subject <- TRUE
        MCMCs <- do.call(as_mcmc_new, c(list(emc, selection = select), fix_dots(cur_dots, as_mcmc_new)))
        names(MCMCs) <- paste0("alpha : ", names(MCMCs))
      } else{
        cur_dots$use_par <- max_par
        MCMCs <- do.call(as_mcmc_new, c(list(emc, selection = select), fix_dots(cur_dots, as_mcmc_new)))
      }
      cur_dots <- add_defaults(cur_dots, xlab = names(MCMCs)[1], ylab = "Highest Rhat parameter")
      do.call(plot, c(list(MCMCs[[1]], auto.layout = FALSE, density = FALSE, ask = FALSE,smooth = FALSE),
                      fix_dots_plot(cur_dots)))
      legend("topleft",legend=paste0("Rhat : ",round(max_gd,digits)), bty = "n")
      legend("topright",legend=paste0("ESS : ", round(ESS[[cur_max]][max_par])), bty = "n")
    }
    out_list[[select]] <- out
  }
  return(invisible(out_list))
}

get_recovery_stats <- function(MCMC, true_MCMC, true_pars, CI){
  quantiles = c(.5 - CI/2, .5, .5 + CI/2)
  if(!is.null(true_MCMC)){
    quants_true <- get_posterior_quantiles(true_MCMC, probs = quantiles)
    true_pars <- quants_true[,"50%"]
  }
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

#' Plots recovery of data generating parameters/samples
#'
#' @param emc An emc object
#' @param true_pars A vector of data-generating parameters or and emc object with data-generating samples
#' @param selection A Character vector. Indicates which parameter types to plot (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param do_CI Boolean. If `TRUE` will also include bars representing the credible intervals
#' @param correlation Character. Which correlation to include in the plot. Options are either `pearson` or `spearman`
#' @param stat Character. Which statistic to include in the plot. Options are either `rmse` or `coverage`
#' @param digits Integer. How many digits to round the statistic and correlation in the plot to
#' @param CI Numeric. The size of the credible intervals. Default is .95 (95%).
#' @param ci_plot_args A list. Optional additional arguments to be passed to plot.default for the plotting of the credible intervals (see `par()`)
#' @param ... A list of optional arguments that can be passed to `as_mcmc_new` or `plot.default` (see `par()`)
#'
#' @return Invisible list with RMSE, coverage, and Pearson and Spearman correlations.
#' @export
#' @examples
#' # Make up some values that resemble posterior samples
#' # Normally this would be true values that were used to simulate the data
#' pmat <- matrix(rnorm(12, mean = c(-1, -.6, -.4, -1.5), sd = .01), ncol = 4, byrow = TRUE)
#' # Conventionally this would be created before one makes data with true values
#' plot_recovery(samplers_LNR, pmat, correlation = "pearson", stat = "rmse")
#'
#' # Similarly we can plot other group-level parameters with a set of true samples
#' true_samples <- samplers_LNR # Normally this would be data-generating samples
#' recovery_plot(samplers_LNR, true_samples, correlation = "spearman", stat = "coverage", col = "red")
plot_recovery <- function(emc, true_pars,
                          selection = "mu",
                          layout=NA,
                          do_CI = TRUE,
                          correlation = "pearson",
                          stat = "rmse",
                          digits = 3,
                          CI = .95, ci_plot_args = list(), ...)
{
  dots <- list(...)
  type <- attr(emc[[1]], "variant_funs")$type
  if(length(dots$subject) == 1) dots$by_subject <- TRUE
  dots$merge_chains <- TRUE
  MCMC_samples <- do.call(as_mcmc_new, c(list(emc, selection = selection), fix_dots(dots, as_mcmc_new)))
  true_MCMC_samples <- NULL
  if(!is(true_pars, "emc")){
    if(selection == "sigma2" & !is.matrix(true_pars)) true_pars <- diag(true_pars)
    true_pars <- do.call(as_mcmc_new, c(list(emc, selection = selection, type = type, true_pars = true_pars),
                                        fix_dots(dots, as_mcmc_new, exclude = c("thin", "subfilter"))))
  } else{
    true_MCMC_samples <- do.call(as_mcmc_new, c(list(true_pars, selection = selection), fix_dots(dots, as_mcmc_new)))
    true_pars <- NULL
  }
  # pearson <- spearman <- rmse <- coverage <- setNames(numeric(length(MCMC_samples)),names(MCMC_samples))
  stats_list <- list()
  if(any(is.na(layout))){
    par(mfrow = coda:::set.mfrow(Nchains = 1, Nparms = length(MCMC_samples),
                                 nplots = 1))
  } else{par(mfrow=layout)}
  for(i in 1:length(MCMC_samples)){

    cur_name <- names(MCMC_samples)[i]
    stats <- get_recovery_stats(MCMC_samples[[i]], true_MCMC_samples[[i]],
                                                          true_pars[[i]], CI)
    ylim <- range(c(stats$true, stats$recovered))
    main_name <- ifelse(length(MCMC_samples) == 1, cur_name, paste0(selection, ": ", cur_name))
    cur_dots <- add_defaults(dots, main = main_name, ylim = ylim, xlim = ylim, xlab = "Generated", ylab = "Estimated")
    do.call(plot, c(list(stats$true[,"50%"],stats$recovered[,"50%"]), fix_dots_plot(cur_dots)))
    abline(a=0,b=1,lty=3)
    if(do_CI){
      cur_ci_plot_args <- add_defaults(ci_plot_args, col = "grey", angle = 90, length = .05)
      do.call(arrows, c(list(stats$true[,"50%"],stats$recovered[,2],stats$true[,"50%"],stats$recovered[,3]),
                        fix_dots_plot(cur_ci_plot_args)))
      do.call(arrows, c(list(stats$true[,"50%"],stats$recovered[,2],stats$true[,"50%"],stats$recovered[,1]),
                        fix_dots_plot(cur_ci_plot_args)))
      if(is.null(true_pars)){
        do.call(arrows, c(list(stats$true[,"50%"],stats$recovered[,2],stats$true[,3],stats$recovered[,2]),
                          fix_dots_plot(cur_ci_plot_args)))
        do.call(arrows, c(list(stats$true[,"50%"],stats$recovered[,2],stats$true[,1],stats$recovered[,2]),
                          fix_dots_plot(cur_ci_plot_args)))
      }
    }
    if(correlation == "pearson"){
      legend("topleft",paste("r(pearson) = ",round(stats$pearson,digits)),bty="n")
    } else if(correlation == "spearman"){
      legend("topleft",paste("r(spearman) = ",round(stats$spearman,digits)),bty="n")
    }
    if (stat == "rmse") {
      legend("bottomright",paste("RMSE = ",round(stats$rmse,digits)),bty="n")
    } else if(stat == "coverage") {
      legend("bottomright",paste("95% Coverage = ",round(stats$coverage,digits)),bty="n")
    }
    stats <- make_recov_summary(stats)
    stats_list[[cur_name]] <- stats
  }
  return(invisible(stats_list))
}



