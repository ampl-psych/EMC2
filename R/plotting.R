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

# #' Plots choice data
# #'
# #' Plots choice data with no response times.
# #'
# #' @param data A data frame. The experimental data in EMC2 format with at least `subject` (i.e., the
# #' subject factor), `R` (i.e., the response factor) and `rt` (i.e., response time) variable.
# #' Additional factor variables of the design are optional.
# #' @param pp Posterior predictives created by `predict()`
# #' @param subject Integer or string picking out subject(s).
# #' @param factors Character vector of factors in data to display separately. If
# #' `NULL` (i.e., the default), use names of all columns in data except `trials`,`R`, and `rt`.
# #' Omitted factors are aggregated over. If `NA`, treats entire data set as a single cell.
# #' If `stat` is used, the default is changed to `NA`.
# #' @param functions A named list of functions that create new factors which can then be
# #' used by the `factors` and `stat` arguments.
# #' @param stat A function that takes a the data and returns a single value.
# #' @param stat_name A string naming what the `stat` argument calculates.
# #' @param adjust Control of smoothing in density plots
# #' @param ci Credible interval and central tendency quantiles for return when
# #' stat argument is supplied (defaults to the 2.5\\%, the 50\\% and the 97.5\\%
# #' quantiles)
# #' @param do_plot Boolean (defaults to `TRUE`) whether a plot should be created or not
# #' @param xlim x-axis plot limit, 2-vector (same for all) or matrix (one row for each paramter)
# #' @param ylim y-axis plot limit, 2-vector (same for all) or matrix (one row for each paramter)
# #' @param main Text title, pasted before cell name.
# #' @param layout 2-vector specifying `par(mfrow)` or `par(mfcol)`. The default `NULL` uses current,
# #' `NA` keeps `par` currently active.
# #' @param mfcol Boolean for `layout` settings, the default `TRUE` uses `mfcol`, else `mfrow`.
# #' @param signalFactor Character name of factor for the signal
# #' @param zROC Boolean, plot Z transformed ROC (defaults to `FALSE`)
# #' @param qfun Type of Z transform (defaults to probit)
# #' @param lim `x` = `y` limit for ROC plots
# #' @param rocfit_cex Size of points in ROC plot (default 0.5)
# #'
# #' @return If stat argument is provided a matrix of observed values and predicted quantiles
# #' is returned
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

#' Plot Within-Chain Correlations
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
  if(!is(emc, "emc")) stop("input must be an emc object")
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

#' Likelihood Profile Plots
#'
#' Creates likelihood profile plots from a design and the experimental data by
#' varying one model parameter while holding all others constant.
#'
#' @param data A dataframe. Experimental data used, needed for the design mapping
#' @param design A design list. Created using ``design``.
#' @param p_vector Named vector of parameter values (typically created with ``sampled_pars(design)``)
#' @param range Numeric. The max and min will be p_vector + range/2 and p_vector - range/2, unless specified in p_min or p_max.
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param p_min Named vector. If specified will instead use these values for minimum range of the selected parameters.
#' @param p_max Named vector. If specified will instead use these values for maximum range of the selected parameters.
#' @param use_par Character vector. If specified will only plot the profiles for the specified parameters.
#' @param n_point Integer. Number of evenly spaced points at which to calculate likelihood
#' @param n_cores Number of likelihood points evenly spaced between the minimum and maximum likelihood range.
#' @param true_args A list. Optional additional arguments that can be passed to plot.default for the plotting of the true vertical line.
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
                         true_args = list(),
                         ...)

{
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  dots <- list(...)
  lfun <- function(i,x,p_vector,pname,dadm) {
    p_vector[pname] <- x[i]
    calc_ll_R(p_vector, attr(dadm, "model")(), dadm)
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
      do.call(abline, c(list(v=cur_par), fix_dots_plot(add_defaults(true_args, lty = 2))))
      out[cur_name,] <- c(p_vector[cur_name], x[which.max(ll)], p_vector[cur_name] - x[which.max(ll)])
    }
  }
  return(round(out, 3))
}

#' Plots Density for Parameters
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
#' @param prior_args A list. Optional additional arguments to be passed to plot.default for the plotting of the prior density (see `par()`)
#' @param true_args A list. Optional additional arguments to be passed to plot.default for the plotting of the true parameters (see `par()`)
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
#' plot_pars(samples_LNR, true_pars = true_emc, true_args = list(col = 'blue'), adjust = 2)
plot_pars <- function(emc,layout=NA, selection="mu", show_chains = FALSE, plot_prior = TRUE, N = 1e4,
                      use_prior_lim = !all_subjects, lpos = "topright", true_pars = NULL, all_subjects = FALSE,
                      prior_args = list(), true_args = list(), ...)
{
  if(!is(emc, "emc")) stop("input must be an emc object")
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  dots <- list(...)
  type <- emc[[1]]$type
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
  psamples <-  do.call(get_objects, c(list(sampler = emc, design = get_design(emc),
                           type = type, sample_prior = T,
                           selection = selection, N = N,
                           prior = get_prior(emc)), fix_dots(dots, get_objects,consider_dots = F)))
  pMCMC_samples <- do.call(get_pars, c(list(psamples, selection = selection, type = type),
                                       fix_dots(dots, get_pars, exclude = c("thin", "filter", "chain", "subject"))))
  if(length(pMCMC_samples) < length(MCMC_samples)) pMCMC_samples <- rep(pMCMC_samples, length(MCMC_samples))

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
    all_cols <- unique(unlist(lapply(cur_mcmc, colnames)))
    cur_mcmc <- lapply(cur_mcmc, function(m) {
      m <- cbind(m, matrix(NA_real_, nrow(m), length(setdiff(all_cols, colnames(m))),
                           dimnames = list(NULL, setdiff(all_cols, colnames(m)))))
      m[, all_cols, drop = FALSE]
    })
    merged <- do.call(rbind, cur_mcmc)
    # merged <- do.call(rbind, cur_mcmc_filled)
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
      denses <- lapply(cur_mcmc, function(x){
        if(is.na(x[1,l])) return(NULL)
        do.call(density, c(list(x[,l]), fix_dots(dots, density.default, consider_dots = FALSE)))
      })
      denses <- denses[!sapply(denses, is.null)]
      xlim <- range(c(sapply(cur_mcmc, function(x) {
        if(is.na(x[1,l])) return(NA)
        range(x[,l])
        }), true_pars[[i]][[1]][,l]), na.rm = TRUE)
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
        true_denses <- lapply(true_MCMC_samples[[i]], function(x){
          if(is.na(x[1,l])) return(NULL)
          do.call(density, c(list(x[,l]), fix_dots(dots, density.default, consider_dots = FALSE)))
        })
        xlim <- range(c(sapply(true_MCMC_samples[[i]], function(x) {
          if(is.na(x[1,l])) return(NULL)
          range(x[,l])
        }), true_pars[[i]][[1]][,l]))
        ylim <- range(c(sapply(true_denses, function(x) return(range(x$y))), ylim))
      }
      n_chains <- length(denses)
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
          cur_true_args <- add_defaults(true_args, col = "darkgreen")
          do.call(lines, c(list(true_denses[[k]]), fix_dots_plot(cur_true_args)))
        }
      }
      if(!is.null(true_pars)){
        cur_true_args <- add_defaults(true_args, col = "darkgreen", lwd = 1.5, lty = 2)
        do.call(abline, c(list(v = true_pars[[i]][[1]][,l]), fix_dots_plot(cur_true_args)))
      }
      cur_contraction[l] <- 1-(var(merged[,l], na.rm = TRUE)/var(pMCMC_samples[[i]][[1]][,p_idx]))
      if(plot_prior){
        cur_prior_args <- add_defaults(prior_args, col = "red")
        do.call(lines, c(list(pdenses), fix_dots_plot(cur_prior_args)))
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


# #' Plot MCMC
# #'
# #' Uses the coda plot functions that are applied per chain
# #'
# #' @param emc An emc object
# #' @param selection A Character string. Indicates which parameter type to plot (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
# #' @param fun A plot function that takes a vector/mcmc object as input, e.g. cumuplot, acf
# #' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
# #' @param plot_type type argument passed on to coda fun.
# #' @param chain Integer, which chain to include, if more than 1 will make separate plots per chain.
# #' @param ... Optional arguments that can be passed to `get_pars`,
# #' the chosen coda plot function, or `plot.default` (see `par()`)
# #'
# #' @return A coda plot
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

# #' Plot MCMC.list
# #'
# #' Uses the coda plot functions that are applied across chain
# #'
# #' @param emc An emc object
# #' @param selection A Character string. Indicates which parameter type to plot (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
# #' @param fun A coda plot function choice from
# #' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
# #' @param ... Optional arguments that can be passed to `get_pars`,
# #' the chosen coda plot function, or `plot.default` (see `par()`)
# #'
# #' @return A coda plot
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


#' Plots trends over time
#'
#' Plots the trend for selected parameters of a model. Can be used either with a p_vector, or trial-wise parameters or covariates obtained
#' from predict()
#'
#' @param input_data a p_vector or posterior predictives compatible with the provided emc object
#' @param emc An emc object
#' @param par_name Parameter name (or covariate name) to plot
#' @param subject Subject number to plot
#' @param filter Optional function that takes a data frame and returns a logical vector indicating which rows to include in the plot
#' @param on_x_axis Column name in the `dadm` to plot on the x-axis. By default 'trials'.
#' @param pp_shaded Boolean. If `TRUE` will plot 95% credible interval as a shaded area. Otherwise plots separate lines for each iteration of the posterior predictives. Only applicable if `input_data` are posterior predictives.
#' @param ... Optional arguments that can be passed to `plots`.
#'
#' @return A trend plot
#' @export
#'
#' @examples
#' dat <- EMC2:::add_trials(forstmann)
#' dat$trials2 <- dat$trials/1000
#'
#' lin_trend <- make_trend(cov_names='trials2',
#'                         kernels = 'exp_incr',
#'                         par_names='B',
#'                         bases='lin',
#'                         phase = "premap")
#'
#' design_RDM_lin_B <- design(model=RDM,
#'                            data=dat,
#'                            covariates='trials2',   # specify relevant covariate columns
#'                            matchfun=function(d) d$S==d$lR,
#'                            transform=list(func=c('B'='identity')),
#'                            formula=list(B ~ 1, v ~ lM, t0 ~ 1),
#'                            trend=lin_trend)       # add trend
#'
#' emc <- make_emc(dat, design=design_RDM_lin_B)
#' p_vector <- c('B'=1, 'v'=1, 'v_lMTRUE'=1, 't0'=0.1, 'B.w'=1, 'B.d_ei'=1)
#'
#' # Visualize trend
#' plot_trend(p_vector, emc=emc,
#'            par_name='B', subject='as1t',
#'            filter=function(d) d$lR=='right', main='Threshold for right')
plot_trend <- function(input_data, emc, par_name, subject=1,
                       filter=NULL, on_x_axis='trials',
                       pp_shaded=TRUE,
                       ...) {

  if(!is.list(input_data)) {
    # user supplied p_vector
    dadm <- emc[[1]]$data[[subject]]
  } else {
    dadm <- do.call(rbind, emc[[1]]$data)
  }
  row_filter <- dadm$subjects==subject
  if(!is.null(filter)) {
    if(is.function(filter)) {
      row_filter <- row_filter & filter(dadm)
    } else {
      stop("filter must be a function that takes a data frame and returns a logical vector")
    }
  }
  if(!is.list(input_data)) {
    updated <- get_pars_matrix(p_vector=input_data,
                               dadm=dadm,
                               model=emc[[1]]$model())
    trend <- updated[row_filter, par_name]
    credible_interval <- NULL
    ylim <- range(trend)
    x <- dadm[row_filter, on_x_axis]
    trend <- trend[order(x)]
    x <- x[order(x)]
  } else {
    # user supplied posterior predictives
    updated <- lapply(seq_along(input_data), function(i) data.frame(input_data[[i]][row_filter,par_name,drop=FALSE], x=dadm[row_filter,on_x_axis]))
    credible_interval <- aggregate(.~x, do.call(rbind, updated), quantile, c(0.025, .5, .9))

    ## order
    credible_interval <- credible_interval[order(credible_interval[[1]]),]
    trend <- credible_interval[[2]][,2]  # median
    x <- credible_interval[[1]]
    if(pp_shaded) {
      ylim <- range(credible_interval[[2]])
    } else {
      ylim <- range(sapply(updated, function(x) range(x[,par_name])))
    }
  }

  # parse dots for plotting arguments
  dots <- list(...)
  if(!'xlab' %in% names(dots)) dots$xlab=on_x_axis
  if(!'ylim' %in% names(dots)) dots$ylim=ylim

  full_args <- c(list(x=x,
                      y=trend,
                      type='l', ylab=par_name), dots)
  do.call(plot, full_args)

  if(!is.null(credible_interval)) {
    if(pp_shaded) {
      polygon(x=c(x, rev(x)),
              y=c(credible_interval[[2]][,1], rev(credible_interval[[2]][,3])), col=adjustcolor(3, alpha.f = .3), border=adjustcolor(1, alpha.f=.4))
    } else {
      for(i in updated) lines(i$x, i[,par_name], lwd=.2)
    }
  }
}


# Spectrum plotting -------------------------------------------------------
order_pp <- function(df) {
  df[order(df$subjects, df$postn, df$trials), ]
}

## Returns a function that interpolates power at common frequencies
make_interpolator <- function(freqs, power) {
  approxfun(x = freqs, y = power, rule = 2)
}

## Given a list of (freq, power) per subject, align in common frequency space
interpolate_to_common_grid <- function(freqs_list, power_list) {
  # use the subject with fewest frequency bins
  idx_min <- which.min(lengths(freqs_list))
  common_f <- freqs_list[[idx_min]]

  interpolators <- mapply(make_interpolator, freqs_list, power_list, SIMPLIFY = FALSE)
  interpolated  <- t(sapply(interpolators, function(f) f(common_f)))

  list(freq = common_f, power = interpolated)
}

#' Compute Power Spectra With Optional Subject-Level Aggregation
#'
#' Computes power spectral density estimates using \code{\link{spectrum}},
#' optionally aggregated across subjects or posterior predictive samples.
#' All arguments intended for the underlying spectral estimator should be
#' supplied through \code{spectrum.args}.
#'
#' @param data A data frame with reaction time data. Must contain
#'   \code{subjects} and \code{rt}, and for posterior predictive data
#'   optionally \code{postn} and \code{trials}.
#'
#' @param by.postn Logical. If \code{TRUE}, compute a separate spectrum
#'   for each posterior predictive draw and each posterior sample index.
#'
#' @param spectrum.args A named list of arguments passed directly to
#'   \code{\link{spectrum}}. These override the defaults internally used
#'   in this function. Useful for customizing smoothing spans, detrending,
#'   tapering, and so on.
#'   Defaults: `list(spans=c(3, 5), detrend=TRUE, demean=FALSE, log=FALSE)`.
#'   By default, we run `spectrum` without `log`, and log-transform while plotting
#'
#' @details
#' The function organizes the data by subject (and optionally posterior
#' sample index), computes spectra individually, interpolates spectra to
#' a common frequency grid if needed, and averages them appropriately.
#'
#' @importFrom utils modifyList
#' @importFrom stats spec.pgram
#' @return
#' Either a data frame with columns \code{freq} and \code{power}, or (if
#' \code{by.postn = TRUE}) a list with frequency vector and a matrix of
#' spectra across posterior samples.
#'
#' @export
get_power_spectra <- function(data,
                              by.postn = FALSE,
                              spectrum.args = list()) {

  # Default spectrum parameters (user may override any of them)
  # Default from Wagenmakers et al., (2004), except demean,
  # which we set to FALSE to capture the power at f=0.
  # We don't get the log, we transform the log later.
  default_spec_args <- list(
    spans   = c(3, 5),
    detrend = TRUE,
    demean  = FALSE,
    log     = FALSE,
    taper   = 0.1  # default from spec.pgram
    fast    = TRUE
  )

  spec_args <- modifyList(default_spec_args, spectrum.args)

  # Helper to compute a spectrum for vector x
  compute_spec <- function(x, spec_args) {
    spec.pgram(
      x = x,
      deamean = spec_args$demean,
      spans   = spec_args$spans,
      taper   = spec_args$taper,
      detrend = spec_args$detrend,
      log     = spec_args$log,
      plot    = FALSE,
      fast    = spec_args$fast   # FFT acceleration
    )
  }

  # --------------------------
  # Case 1: Mean posterior predictive spectra
  # --------------------------
  # if (mean.pp) {
  #   pp <- data[order(data$subjects, data$postn, data$trials), ]
  #
  #   power_df <- aggregate(rt ~ postn * subjects, pp, function(x) compute_spec(x, spec_args=spec_args)$spec)
  #   freq_df  <- aggregate(rt ~ postn * subjects, pp, function(x) compute_spec(x, spec_args=spec_args)$freq)
  #
  #   # Unlist spectra & frequencies into matrices
  #   if (is.list(power_df$rt)) power_df$rt <- do.call(rbind, power_df$rt)
  #   if (is.list(freq_df$rt))  freq_df$rt  <- do.call(rbind, freq_df$rt)
  #
  #   power_by_subj <- lapply(unique(power_df$subjects),
  #                           function(s) colMeans(power_df[power_df$subjects == s, "rt"]))
  #
  #   freq_by_subj  <- lapply(unique(freq_df$subjects),
  #                           function(s) colMeans(freq_df[freq_df$subjects == s, "rt"]))
  #
  #   mean_power <- colMeans(do.call(rbind, power_by_subj))
  #   mean_freq  <- colMeans(do.call(rbind, freq_by_subj))
  #
  #   return(data.frame(freq = mean_freq, power = mean_power))
  # }

  # --------------------------
  # Case 2: Posterior predictive spectra by postn
  # --------------------------
  if (by.postn) {
    pp <- data[order(data$subjects, data$postn, data$trials), ]
    subjects <- unique(pp$subjects)

    spectra_by_subj <- lapply(subjects, function(s) {
      lapply(unique(pp[pp$subjects == s, "postn"]),
             function(pn) compute_spec(pp[pp$subjects == s & pp$postn == pn, "rt"], spec_args=spec_args))
    })

    # Extract frequencies (one per subject)
    freq_by_subj <- lapply(spectra_by_subj, function(lst) lst[[1]]$freq)
    power_by_subj <- lapply(spectra_by_subj,
                            function(lst) lapply(lst, function(s) s$spec))

    n_postn <- length(power_by_subj[[1]])
    n_subj  <- length(power_by_subj)

    # Interpolate if subjects differ in trial count
    min_len <- min(sapply(freq_by_subj, length))
    freq_grid <- freq_by_subj[[which.min(sapply(freq_by_subj, length))]]

    interp <- function(freqs, powers) {
      sapply(powers, function(pwr) approxfun(freqs, pwr)(freq_grid))
    }

    power_interp <- lapply(1:n_subj, function(i) {
      interp(freq_by_subj[[i]], power_by_subj[[i]])
    })

    mean_power <- do.call(rbind,
                          lapply(1:n_postn, function(pn) {
                            colMeans(do.call(rbind, lapply(power_interp,
                                                           function(m) m[, pn])))
                          }))

    return(list(freq = freq_grid, power = mean_power))
  }

  # --------------------------
  # Case 3: Simple subject-averaged spectrum
  # --------------------------
  spectra <- lapply(unique(data$subjects),
                    function(s) compute_spec(data[data$subjects == s, "rt"], spec_args=spec_args))

  freqs <- lapply(spectra, function(s) s$freq)
  specs <- lapply(spectra, function(s) s$spec)

  # Interpolate to common frequency grid (shortest)
  res <- interpolate_to_common_grid(freqs, specs)

  freq_grid  <- res$freq
  interp_mat <- res$power         # matrix: subjects × freq bins
  mean_power <- colMeans(interp_mat)

  data.frame(freq = freq_grid, power = mean_power)
}



#' Plot Empirical and Posterior Predictive Power Spectra
#'
#' Computes and plots the empirical power spectrum, with optional overlay
#' of spectra from posterior predictive simulations. All customization of
#' the spectral estimator is done through \code{spectrum.args}, while plot
#' appearance is controlled via \code{...}.
#'
#' @param dat A data frame containing empirical reaction time data, with
#'   at least columns \code{subjects} and \code{rt}.
#'
#' @param pp Optional posterior predictive data in the same format as
#'   \code{dat}, including \code{subjects}, \code{postn}, and \code{trials}.
#'
#' @param plot.log Logical. Whether to log-transform frequencies and
#'   power before plotting. This does not affect the call to
#'   \code{\link{spectrum}}: use \code{spectrum.args$list(log = TRUE)}
#'   to request log spectral density from the estimator itself.
#'
#' @param spectrum.args A named list of arguments forwarded directly to
#'   \code{\link{spectrum}} inside \code{\link{get_power_spectra}}.
#'
#' @param trial_duration Optional duration of a trial in seconds. If
#'   supplied, the x-axis is labeled in human-readable time units. Otherwise
#'   the x-axis is in (log) frequencies of 1/trial.
#'
#' @param ... Additional graphical parameters passed to \code{\link{plot}}.
#'
#' @importFrom graphics axis
#' @return Invisibly returns a list containing the empirical spectrum
#'   and, if posterior predictive data is supplied, the posterior
#'   predictive spectra and their mean.
#'
#' @export
plot_spectrum <- function(dat, pp = NULL,
                          plot.log = TRUE,
                          spectrum.args = list(),
                          trial_duration = NULL,
                          ...) {

  # Compute spectrum
  sp_dat <- get_power_spectra(dat, spectrum.args = spectrum.args)

  f <- if (plot.log) log else identity
  freq_t <- f(sp_dat$freq)
  pow_t  <- f(sp_dat$power)

  if(plot.log) {
    dots <- add_defaults(list(...),
                         'xlab'='Log frequency (1/trial)',
                         'ylab'='Log power',
                         'type'='n')
  } else {
    dots <- add_defaults(list(...),
                         'xlab'='Frequency (1/trial)',
                         'ylab'='Power',
                         'type'='n')
  }

  # Base plot
  if (is.null(trial_duration)) {

    plot_args <- dots
    plot_args$x <- freq_t
    plot_args$y <- pow_t
    do.call(plot, plot_args)

  } else {
    ticks_sec <- c(7200, 3600, 1800, 900, 300, 120, 60, 30, 5, 1)
    ticks_rel <- (1 / ticks_sec) * trial_duration
    tick_x    <- f(ticks_rel)

    labels <- c("2 hr", "1 hr", "30 m", "15 m", "5 m",
                "2 m", "1 m", "30 s", "5 s", "1 s")

    if(!'xaxt' %in% names(dots)) dots$xaxt <- 'n'
    plot_args <- dots
    plot_args$x <- freq_t
    plot_args$y <- pow_t
    do.call(plot, dots)

    axis(1,
         at     = tick_x,
         labels = labels,
         las    = 2)
  }
  result <- list(dat = sp_dat)

  # Overlay posterior predictive spectra
  if (!is.null(pp)) {
    sp_pp <- get_power_spectra(pp, by.postn = TRUE,
                               spectrum.args = spectrum.args)

    freqs_pp <- f(sp_pp$freq)
    power_pp <- sp_pp$power

    # Individual PP spectra
    for (i in 1:nrow(power_pp)) {
      lines(freqs_pp, f(power_pp[i, ]),
            col = adjustcolor("darkgreen", alpha.f = 0.1))
    }

    # Posterior predictive mean
    lines(freqs_pp, f(colMeans(power_pp)),
          col = "darkgreen")

    result$pp      <- sp_pp
    result$pp_mean <- data.frame(freq = sp_pp$freq,
                                 power = colMeans(sp_pp$power))
  }

  # Empirical spectrum trace
  lines(freq_t, pow_t, col = par()$col)

  invisible(result)
}


