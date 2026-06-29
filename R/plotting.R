robust_hist <- function(ps, breaks = 50, cutoff = 0.0015, do_plot = TRUE, prob = FALSE, ...){
  ps <- ps[abs(ps) < 1000]
  for(i in 1:log10(length(ps))){
    cuts <- cut(ps, breaks = breaks)
    tab_cuts <- table(cuts)
    good_cuts <- tab_cuts/length(ps) > cutoff
    ps <- ps[cuts %in% names(good_cuts)[good_cuts]]
  }
  if(do_plot){
    hist(ps, breaks = breaks, freq = !prob, ...)
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
  ctab <- .choice_roc_coords(data, signalFactor = signalFactor, zROC = zROC,
                             qfun = if (is.null(qfun)) qnorm else qfun)
  if (!zROC) {
    if (!is.null(lim)) {xlim <- lim; ylim <- lim} else
    {xlim <- c(0,1); ylim <- c(0,1)}
    plot(ctab[, "x"],ctab[, "y"],xlab="p(FA)",ylab="p(H)",xlim=xlim,ylim=ylim,main=main)
    lines(ctab[, "x"],ctab[, "y"])
    abline(a=0,b=1,lty=3)
  } else {
    if (is.null(lim)) lim <- c(min(ctab),max(ctab))
    plot(ctab[, "x"],ctab[, "y"],main=main,
         xlab="z(FA)",ylab="z(H)",xlim=lim,ylim=lim)
    lines(ctab[, "x"],ctab[, "y"])
    abline(a=0,b=1,lty=3)
  }
  invisible(ctab)
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
#' # profile_plot(p_vector = p_vector, p_min = c(t0 = -1.35),
#' #              p_max = c(t0 = -1.45), use_par = c("a", "t0", "SZ"),
#' #              data = forstmann, design = design_DDMaE, n_point = 10)
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
  # allow for setting par outside of this function (needed because par(mfrow=...) always overwrites the cex argument as well)
  set.par <- !isFALSE(dots$set.par)
  dots$set.par <- NULL
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
    if(set.par) {
      if(any(is.na(layout))){
        par(mfrow = coda_setmfrow(Nchains = n_chains, Nparms = n_pars, nplots = 1))
      } else{
        par(mfrow=layout)
      }
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
#' Plots the trend for selected parameters of a model. Can be used either with a p_vector,
#' or trial-wise parameters or covariates obtained from predict(). Optionally overlays
#' a ground-truth trajectory for recovery evaluation.
#'
#' @param input_data A p_vector or posterior predictives compatible with the provided emc object.
#'   If posterior predictives, should be the output of \code{predict(..., return_trialwise_parameters = TRUE)}.
#' @param emc An emc object.
#' @param par_name Parameter name (or covariate name) to plot.
#' @param true_trajectories Optional. A data frame with columns \code{subjects}, the column named
#'   by \code{on_x_axis}, and at least \code{par_name}. If supplied, the true trajectory is
#'   overlaid as a dashed black line, and \code{ylim} expands to include it.
#' @param subject Subject to plot. If NULL, plots all subjects separately.
#'   Ignored when \code{group_average = TRUE}.
#' @param filter Optional function that takes a data frame and returns a logical vector
#'   indicating which rows to include. For race models, should select one accumulator,
#'   e.g. \code{function(d) d$lR == 'left'}.
#' @param on_x_axis Column name to plot on the x-axis. Default \code{'trials'}.
#' @param pp_shaded Logical. If TRUE, plots 95% credible interval as a shaded area.
#'   Otherwise plots separate lines for each posterior predictive iteration.
#'   Only applicable when \code{input_data} is posterior predictives.
#' @param group_average Logical. If TRUE, plots a single panel averaged across subjects
#'   rather than one panel per subject.
#' @param layout Optional integer vector of length 2 to override the automatic mfrow layout.
#'   Ignored when \code{group_average = TRUE}.
#' @param ... Optional arguments passed to \code{plot()}.
#'
#' @return Invisibly returns the credible interval data frame when \code{input_data} is
#'   posterior predictives, otherwise NULL.
#' @export
plot_trend <- function(input_data, emc, par_name,
                       true_trajectories = NULL,
                       subject           = NULL,
                       filter            = NULL,
                       on_x_axis         = 'trials',
                       pp_shaded         = TRUE,
                       group_average     = FALSE,
                       layout            = NULL,
                       ...) {

  all_subjects <- names(emc[[1]]$data)
  dots <- list(...)

  # --- Subject validation ---
  if (!is.null(subject) && !group_average) {
    if (!subject %in% all_subjects)
      stop(sprintf("Subject '%s' not found. Available: %s",
                   subject, paste(all_subjects, collapse = ", ")))
    subjects_to_plot <- subject
    multi <- FALSE
  } else {
    subjects_to_plot <- all_subjects
    multi <- !group_average
  }

  # --- Prepare posterior predictives if input_data is a list ---
  # Also merge dadm columns (e.g. lR) needed for filtering
  if (is.list(input_data)) {
    # validate existence of trialwise pars
    if('trialwise_parameters' %in% names(attributes(input_data))) {
      # full pp is provided
      input_data_pp <- lapply(attr(input_data, 'trialwise_parameters'), data.frame)
    } else {
      # user extracted trialwise parameters
      input_data_pp <- lapply(input_data, data.frame)
    }

    dadm_all       <- do.call(rbind, emc[[1]]$data)

    # add postn number and dadm columns for aggregation
    for (post_n in seq_along(input_data_pp)) {
      input_data_pp[[post_n]] <- cbind(
        postn=post_n,
        input_data_pp[[post_n]],
        dadm_all[, setdiff(names(dadm_all), names(input_data_pp[[post_n]])), drop = FALSE]
      )
    }
    if (!is.null(true_trajectories)) {
      true_trajectories <- data.frame(true_trajectories)
      true_trajectories <- cbind(
        true_trajectories,
        dadm_all[, setdiff(names(dadm_all), names(true_trajectories)), drop = FALSE]
      )
    }
  } else {
    input_data_pp <- NULL
  }

  # --- Aggregate posterior predictives (per-subject path only) ---
  if (!is.null(input_data_pp) && !group_average) {
    combined <- do.call(rbind, input_data_pp)
    if (!is.null(filter)) combined <- combined[filter(combined), ]
    agg <- aggregate(
      stats::as.formula(paste0(par_name, " ~ ", on_x_axis, " * subjects")),
      data = combined, stats::quantile, c(0.025, 0.5, 0.975)
    )
    agg <- agg[order(agg$subjects, agg[[on_x_axis]]), ]
  } else {
    agg <- NULL
  }

  # --- Group average path ---
  if (group_average) {
    if (!is.null(input_data_pp)) {
      # posterior predictives are proviced. First filter
      combined <- do.call(rbind, input_data_pp)
      if (!is.null(filter)) combined <- combined[filter(combined), ]

      # aggregate across subjects first, then get 95% CI and median
      agg_grp <- aggregate(
        stats::as.formula(paste0(par_name, " ~ ", on_x_axis, " * postn")), combined, mean)
      agg_ci <- aggregate(
        stats::as.formula(paste0(par_name, " ~ ", on_x_axis)), agg_grp, stats::quantile, c(.025, .5, .975)
      )
      agg_ci <- agg_ci[order(agg_ci[[on_x_axis]]),]
      ci_mat  <- agg_ci[[par_name]]
      colnames(agg_ci) <- c(on_x_axis, par_name)
      xs      <- agg_ci[[on_x_axis]]
      ci      <- agg_ci[[par_name]]
      trend   <- ci[, 2]
      ylim    <- range(ci)
    } else {
      # input_data is a p_vector, should be a p_matrix in this case
      # check if p_vector has the same number of rows as there's subjects in the dadm
      n_subj <- length(emc[[1]]$data)
      if(!isTRUE(nrow(input_data)==n_subj)) stop("For group aggregates, provide either posterior predictives or a parameter matrix with the same number of rows as subjects")

      updated <- do.call(rbind, lapply(1:nrow(input_data), function(i) get_pars_matrix_oo(input_data[i,,drop=FALSE], emc[[1]]$data[[i]], model = emc[[1]]$model())))

      # get dadm
      dadm_all <- do.call(rbind, emc[[1]]$data)
      if (!is.null(filter)) {
        idx <- filter(dadm_all)
        dadm_all <- dadm_all[idx, ]
        updated <- updated[idx, ]
      }

      grp_df   <- data.frame(x = dadm_all[[on_x_axis]], trend = updated[, par_name])
      grp_agg  <- aggregate(trend ~ x, grp_df, mean)  # across-subject mean
      grp_agg  <- grp_agg[order(grp_agg$x), ]
      xs       <- grp_agg$x
      trend    <- grp_agg$trend
      ci       <- NULL
      ylim     <- range(trend)
      agg_ci <- NULL
    }

    if (!is.null(true_trajectories)) {
      true_filt <- if (!is.null(filter)) true_trajectories[filter(true_trajectories), ] else true_trajectories
      true_grp  <- aggregate(stats::as.formula(paste0(par_name, " ~ ", on_x_axis)), true_filt, mean)
      true_grp <- true_grp[order(true_grp[[on_x_axis]]), ]
      ylim     <- range(ylim, true_grp[[par_name]])
    }

    if (!'xlab' %in% names(dots)) dots$xlab <- on_x_axis
    if (!'ylab' %in% names(dots)) dots$ylab <- par_name
    if (!'ylim' %in% names(dots)) dots$ylim <- ylim
    if (!'main' %in% names(dots)) dots$main <- "Group average"

    do.call(plot, c(list(x = xs, y = trend, type = 'l', col = ifelse(is.null(ci), 1, 'darkgreen'), lwd = 1.5), dots))

    if (!is.null(ci)) {
      if (pp_shaded) {
        polygon(c(xs, rev(xs)), c(ci[, 1], rev(ci[, 3])),
                col = adjustcolor(3, alpha.f = 0.3), border = NA)
      } else {
        for (draw in input_data_pp) {
          draw_filt <- if (!is.null(filter)) draw[filter(draw), ] else draw
          draw_agg  <- aggregate(
            stats::as.formula(paste0(par_name, " ~ ", on_x_axis)),
            data = draw_filt, FUN = mean
          )
          draw_agg <- draw_agg[order(draw_agg[[on_x_axis]]), ]
          lines(draw_agg[[on_x_axis]], draw_agg[[par_name]],
                lwd = 0.2, col = adjustcolor(3, alpha.f = 0.5))
        }
      }
    }

    if (!is.null(true_trajectories)) {
      lines(true_grp[[on_x_axis]], true_grp[[par_name]],
            col = 'black', lwd = 2, lty = 2)
      legend("topright", legend = c("95% CI", "True"),
             col = c(NA, "black"),
             lty = c(NA, 2), lwd = c(NA, 2),
             pch = c(22, NA),
             pt.bg = c(adjustcolor(3, alpha.f = 0.3), NA),
             pt.cex = 2, bty = "n")
      }

    return(invisible(agg_ci))
  }

  # --- Layout for multi-subject plot ---
  if (multi) {
    if (!isFALSE(dots$set.par)) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      if (!is.null(layout)) {
        par(mfrow = layout)
      } else {
        n     <- length(subjects_to_plot)
        ncols <- min(ceiling(sqrt(n)), 3L)
        nrows <- min(ceiling(n / ncols), 3L)
        par(mfrow = c(nrows, ncols))
      }
    }
    dots$set.par <- NULL
  }

  # --- Per-subject plot loop ---
  for (subj in subjects_to_plot) {
    do.call(.plot_trend_single,
            c(list(input_data        = input_data,
                   input_data_pp     = input_data_pp,
                   agg               = agg,
                   true_trajectories = true_trajectories,
                   emc               = emc,
                   par_name          = par_name,
                   subject           = subj,
                   filter            = filter,
                   on_x_axis         = on_x_axis,
                   pp_shaded         = pp_shaded),
              dots))
  }

  invisible(agg)
}


# Internal: plots trend for a single subject
.plot_trend_single <- function(input_data, input_data_pp = NULL, agg = NULL,
                               true_trajectories = NULL,
                               emc, par_name, subject,
                               filter    = NULL,
                               on_x_axis = 'trials',
                               pp_shaded = TRUE,
                               ...) {

  dadm <- emc[[1]]$data[[subject]]

  row_filter <- rep(TRUE, nrow(dadm))
  if (!is.null(filter)) {
    if (is.function(filter)) {
      row_filter <- filter(dadm)
    } else {
      stop("filter must be a function returning a logical vector")
    }
  }

  if (is.null(input_data_pp)) {
    # --- p_vector path ---
    updated <- get_pars_matrix_oo(p_vector = input_data,
                                  dadm     = dadm,
                                  model    = emc[[1]]$model())
    trend  <- updated[row_filter, par_name]
    x      <- dadm[row_filter, on_x_axis]
    ord    <- order(x)
    trend  <- trend[ord]
    x      <- x[ord]
    ylim   <- range(trend)
    ci     <- NULL
  } else {
    # --- posterior predictives path ---
    agg_sub <- agg[agg$subjects == subject, ]
    agg_sub <- agg_sub[order(agg_sub[[on_x_axis]]), ]
    ci      <- agg_sub[[par_name]]
    x       <- agg_sub[[on_x_axis]]
    trend   <- ci[, 2]
    ylim    <- if (pp_shaded) range(ci) else
      range(sapply(input_data_pp, function(d) {
        d_sub <- d[d$subjects == subject, ]
        if (!is.null(filter)) d_sub <- d_sub[filter(d_sub), ]
        range(d_sub[[par_name]])
      }))
  }

  # Expand ylim to include true trajectory if supplied
  if (!is.null(true_trajectories)) {
    true_sub <- true_trajectories[true_trajectories$subjects == subject, ]
    if (!is.null(filter)) true_sub <- true_sub[filter(true_sub), ]
    ylim <- range(ylim, true_sub[[par_name]])
  }

  dots <- list(...)
  if (!'xlab' %in% names(dots)) dots$xlab <- on_x_axis
  if (!'ylab' %in% names(dots)) dots$ylab <- par_name
  if (!'ylim' %in% names(dots)) dots$ylim <- ylim
  if (!'main' %in% names(dots)) dots$main <- as.character(subject)

  do.call(plot, c(list(x = x, y = trend, type = 'l', col = ifelse(is.null(ci), 1, 'darkgreen'), lwd = 1.5), dots))

  # --- Credible interval / individual lines ---
  if (!is.null(ci)) {
    if (pp_shaded) {
      polygon(c(x, rev(x)), c(ci[, 1], rev(ci[, 3])),
              col = adjustcolor(3, alpha.f = 0.3), border = NA)
    } else {
      for (draw in input_data_pp) {
        draw_sub <- draw[draw$subjects == subject, ]
        draw_sub <- draw_sub[order(draw_sub[[on_x_axis]]), ]
        if (!is.null(filter)) draw_sub <- draw_sub[filter(draw_sub), ]
        lines(draw_sub[[on_x_axis]], draw_sub[[par_name]],
              lwd = 0.2, col = adjustcolor(3, alpha.f = 0.5))
      }
    }
  }

  # --- True trajectory overlay ---
  if (!is.null(true_trajectories)) {
    true_sub <- true_trajectories[true_trajectories$subjects == subject, ]
    if (!is.null(filter)) true_sub <- true_sub[filter(true_sub), ]
    true_sub <- true_sub[order(true_sub[[on_x_axis]]), ]
    lines(true_sub[[on_x_axis]], true_sub[[par_name]],
          col = 'black', lwd = 2, lty = 2)
    legend("topright", legend = c("95% CI", "True"),
           col = c(NA, "black"),
           lty = c(NA, 2), lwd = c(NA, 2),
           pch = c(22, NA), pt.bg = c(adjustcolor(3, alpha.f = 0.3), NA),
           pt.cex = 2, bty = "n")
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
#'   Defaults: `list(spans=c(3, 5), detrend=FALSE, demean=TRUE, log=FALSE, taper=0)`.
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
  # Default from Wagenmakers et al., (2004), except detrend,
  # which we set to FALSE to capture the power of linear trends.
  # We don't get the log, we transform the log later.
  default_spec_args <- list(
    spans   = c(3, 5),
    detrend = FALSE,
    demean  = TRUE,
    log     = FALSE,
    taper   = 0,
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

  return(data.frame(freq = freq_grid, power = mean_power))
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
    if(!'xlab' %in% names(list(...))) plot_args$xlab = 'Period'
    do.call(plot, plot_args)

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


# ---------------------------------------------------------------------------
# credint plotting (ported from TC branch: plot_credint, plot_credint_map,
# plot_credints, plot_credints_map + helpers). PR-2 _plot.
# ---------------------------------------------------------------------------

#' Plot credible intervals from credint(map=TRUE)
#'
#' Plots the output of \code{\link{credint}(..., map=TRUE)} for one parameter
#' type, resolving the factorial structure embedded in the parameter names and
#' displaying lower, middle, and upper quantiles as points with capped vertical
#' interval bars.
#'
#' @param ci Length-1 list returned by \code{credint(..., map=TRUE)}.  The
#'   single element must be a matrix with row names encoding parameter types
#'   and factor levels (see Details) and column names that are quantile labels.
#' @param factors Character vector whose length equals the number of \code{_}
#'   characters in the selected row names (the number of factors).  Each
#'   element is a prefix that matches the start of the corresponding
#'   \code{_}-delimited segment; the remainder becomes the level label.  If
#'   \code{NULL} the filtered row names are printed and an informative error
#'   asks the user to supply \code{factors}.  Not required for zero-factor
#'   types.
#' @param type Character string selecting the parameter type to plot: the
#'   prefix of each row name before the first \code{_}.  For parameters with
#'   no \code{_} (e.g. \code{s}, \code{t0}) the full row name is the type and
#'   the parameter is plotted as a single point with an interval.  If
#'   \code{NULL} (default) the first type found in \code{ci[[1]]} is used and
#'   a message names it.
#' @param quants Integer vector of length 3, strictly increasing, indexing the
#'   quantile columns of \code{ci[[1]]} (default \code{1:3}).  Columns are
#'   interpreted as lower, middle (plotted point), and upper.  Errors if fewer
#'   than 3 columns exist or indices are out of range.
#' @param layout Controls the panel grid per page.  Default \code{NA}: choose a
#'   grid automatically -- with exactly 3 factors uses 1 row and
#'   \code{n_levels_f3} columns, with 4+ factors uses \code{c(2, 4)}; factors
#'   beyond the grid iterate across new pages.  An integer vector
#'   \code{c(rows, cols)} sets the grid explicitly.  \code{NULL} leaves the
#'   current \code{par(mfrow)}/\code{par(mfcol)} untouched, so panels are drawn
#'   into the caller's existing layout.
#' @param factor_labels Character vector of length \code{n_factors} \emph{or}
#'   a named list.  As a \strong{character vector}, overrides the factor names
#'   shown in axis labels, panel headings, and the legend title.  As a
#'   \strong{named list}, names must match entries in \code{factors} and each
#'   element is a character vector of replacement \emph{level} labels for that
#'   factor, merged with any \code{level_labels} also supplied
#'   (\code{factor_labels} takes precedence on conflict).
#' @param level_labels Unnamed or named list.  Each element is a character
#'   vector of replacement level names for the corresponding factor.  If
#'   \emph{unnamed}, must have length \code{n_factors} and elements are matched
#'   positionally.  If \emph{named}, names must match entries in \code{factors}
#'   and only those factors are relabelled; other factors keep their
#'   data-derived names.  An error identifies any name mismatch or wrong
#'   vector length.
#' @param level_order Unnamed or named list of integer vectors specifying the
#'   display order of levels for each factor.  Each element must be a
#'   permutation of \code{1:n_levels} for the corresponding factor.  Follows
#'   the same named/unnamed convention as \code{level_labels}: unnamed applies
#'   positionally across all factors; named targets only the specified factors.
#'   Applied after \code{level_labels}, so indices refer to the
#'   data-derived level order (before any relabelling).
#' @param factor_order Integer vector that is a permutation of \code{1:n_factors}.
#'   Reorders which factor plays which role: position 1 becomes the x axis,
#'   position 2 becomes the line groups, positions 3+ become panels.  For
#'   example, with three factors, \code{factor_order = c(2, 1, 3)} swaps the
#'   x axis and line factors.  Applied after \code{level_order}, so indices
#'   refer to the original factor positions as defined by \code{factors}.
#'   Default \code{NULL} leaves the factor order unchanged.
#' @param displace Controls horizontal offsetting of the line groups (levels of
#'   factor 2) so their points and intervals do not overlap at each x position.
#'   A single number (default \code{0.1}) spaces the groups evenly and
#'   symmetrically around each x location with that spacing (e.g. two groups at
#'   \eqn{\pm}\code{displace}/2; three at \code{-displace, 0, +displace}); set
#'   it to \code{0} to disable.  Alternatively a numeric vector of length equal
#'   to the number of line groups gives explicit per-group offsets.  With
#'   \code{credint(..., plot = TRUE)} this can be set globally via \code{...} or
#'   per parameter type via \code{plot_args[[type]]$displace}, the latter taking
#'   precedence.
#' @param cap_length Numeric.  Half-width of the interval end-caps as a
#'   fraction of the spacing between adjacent factor-1 x-positions
#'   (default \code{0.05}).
#' @param pt_cex Numeric. Character expansion factor for plotted points,
#'   passed as \code{cex} to \code{points()}.  Default \code{1}.  Does not
#'   affect text or interval line sizes.
#' @param legend_args Named list of arguments passed to
#'   \code{\link[graphics]{legend}()}, merged over these defaults:
#'   \code{x = "topleft"} (position), \code{bty = "n"} (no box),
#'   \code{title = <factor-2 name>}, \code{cex = 0.75}.  Any argument
#'   accepted by \code{legend()} can be supplied.  To move the legend use
#'   \code{list(x = "bottomright")}.  To suppress the legend pass
#'   \code{list(x = NA)}.  The data-driven arguments \code{legend} (label
#'   text), \code{lty}, \code{pch}, and \code{col} are always set by the
#'   function and cannot be overridden here.  Only drawn when there are two
#'   or more line groups.
#' @param ... Additional graphical arguments forwarded to \code{plot()}.  The
#'   following are intercepted with defaults before forwarding:
#'   \describe{
#'     \item{\code{xlab}}{X-axis label.  Default: name of factor 1 (after
#'       applying \code{factor_labels}), or \code{type} for zero-factor types.}
#'     \item{\code{ylab}}{Y-axis label.  Default: \code{type}.}
#'     \item{\code{main}}{Panel title(s).  Default: empty (\code{""}) for a
#'       single-panel type (fewer than three factors); factor-level combination
#'       labels (e.g. \code{"E: accuracy"}) for the multiple panels of a
#'       three-or-more-factor type.  A length-1
#'       \code{main} overrides the title of every panel with the same string;
#'       a character vector of length equal to the number of panels sets one
#'       title per panel (in panel-drawing order).  Any other length is an
#'       error.}
#'     \item{\code{ylim}}{Y-axis limits.  Default \code{NULL}: a single common
#'       range computed automatically to span every panel of this parameter
#'       type, so the panels share one y axis and are directly comparable.  A
#'       length-2 numeric vector \code{c(lower, upper)} applies the same limits
#'       to every panel; a list of such vectors (length equal to the number of
#'       panels) sets each panel's limits individually, in panel-drawing
#'       order.}
#'     \item{\code{lty}}{Line type(s), vector of length n_lines.
#'       Default: \code{1:n_lines}.}
#'     \item{\code{pch}}{Point character(s), vector of length n_lines.
#'       Default: \code{1:n_lines}.}
#'     \item{\code{cols}}{Colour(s), vector of length n_lines.
#'       Default: \code{"black"} for all groups.  Intervals use the same
#'       colour as their line group.}
#'   }
#'
#' @details
#' Row names of \code{ci[[1]]} follow the pattern
#' \code{type_f1level_f2level_...} where each \code{_}-delimited segment after
#' the type is a concatenation of a factor-name prefix and a level identifier
#' (e.g. \code{Sleft} = factor \code{S}, level \code{left}).  The
#' \code{factors} argument supplies these prefixes.
#'
#' Parameters with zero factors (no \code{_}) are plotted as a single point
#' with a vertical interval at x = 1.
#'
#' Panels are laid out with \code{par(mfrow = layout)}.  The function saves
#' and restores \code{par()} on exit.  Rows are ordered with factor 1 varying
#' slowest and the last factor varying fastest.
#'
#' @return Invisibly returns the plotting data frame (one column per factor
#'   plus \code{lower}, \code{mid}, \code{upper}).
#' @examples \donttest{
#' ci <- credint(samples_LNR, map = TRUE)
#' # type defaults to first type found; factors required for multi-factor types
#' plot_credint_map(ci, factors = "lM")
#' # explicit type, label overrides, named level_order to reverse display
#' plot_credint_map(ci, factors = "lM", type = "m",
#'              factor_labels = "Match",
#'              level_labels  = list(lM = c("Non-match", "Match")),
#'              level_order   = list(lM = c(2L, 1L)),
#'              ylab = "Mean log-RT", cols = c("blue", "red"))
#' }
#' @export
plot_credint_map <- function(ci, factors = NULL, type = NULL, quants = 1:3,
                         layout = NA, factor_labels = NULL,
                         level_labels = NULL, level_order = NULL,
                         factor_order = NULL, displace = 0.1,
                         cap_length = 0.05, pt_cex = 1,
                         legend_args = NULL, ...) {

  # ---- 1. Validate ci and quants ----------------------------------------
  if (!is.list(ci) || length(ci) != 1 || !is.matrix(ci[[1]]))
    stop("`ci` must be a length-1 list containing a matrix (output of credint(map=TRUE)).")
  arr <- ci[[1]]
  if (ncol(arr) < 3)
    stop("`ci[[1]]` must have at least 3 quantile columns; found ", ncol(arr), ".")
  if (length(quants) != 3 || !all(diff(quants) > 0))
    stop("`quants` must be a strictly increasing integer vector of length 3.")
  if (any(quants < 1L) || any(quants > ncol(arr)))
    stop("`quants` indices out of range 1:", ncol(arr), ".")

  # ---- 2. Resolve type (default to first) and filter rows ---------------
  rn        <- rownames(arr)
  row_types <- sub("_.*", "", rn)
  if (is.null(type)) {
    type <- unique(row_types)[1L]
    message("type not supplied; using '", type, "'.")
  }
  sel <- row_types == type
  if (!any(sel))
    stop("Type '", type, "' not found. Available: ",
         paste(sort(unique(row_types)), collapse = ", "), ".")
  sub_arr           <- arr[sel, quants, drop = FALSE]
  colnames(sub_arr) <- c("lower", "mid", "upper")
  sel_rn            <- rownames(sub_arr)

  # ---- 3. Count factors (number of _ in row names) ----------------------
  n_factors <- nchar(sel_rn[1L]) - nchar(gsub("_", "", sel_rn[1L], fixed = TRUE))
  nf_all    <- nchar(sel_rn) - nchar(gsub("_", "", sel_rn, fixed = TRUE))
  if (!all(nf_all == n_factors))
    stop("Selected rows have inconsistent numbers of '_' separators.")

  # ---- 4. Parse factorial structure (n_factors > 0) ---------------------
  if (n_factors > 0L) {
    remainders <- sub(paste0("^", type, "_"), "", sel_rn)
    parts      <- strsplit(remainders, "_", fixed = TRUE)

    if (is.null(factors)) {
      cat("Filtered row names for type '", type, "':\n", sep = "")
      print(sel_rn)
      stop("Supply `factors`: a character vector of length ", n_factors,
           " giving the name prefix for each of the ", n_factors, " factor(s).")
    }
    if (length(factors) != n_factors)
      stop("`factors` must have length ", n_factors, " (one per '_' separator), got ",
           length(factors), ".")

    raw_lvl_per_row <- lapply(seq_len(n_factors), function(f) {
      segs <- sapply(parts, `[`, f)
      bad  <- !startsWith(segs, factors[f])
      if (any(bad))
        stop("Factor prefix '", factors[f], "' (position ", f, ") does not match: ",
             paste(unique(segs[bad]), collapse = ", "), ".")
      sub(paste0("^", factors[f]), "", segs)
    })

    pdata <- as.data.frame(sub_arr, stringsAsFactors = FALSE)
    for (f in seq_len(n_factors))
      pdata[[paste0("f", f)]] <- raw_lvl_per_row[[f]]

    # Capture level order from row names before sorting
    factor_levels_list <- lapply(seq_len(n_factors), function(f)
      unique(raw_lvl_per_row[[f]]))

    # Sort pdata respecting row-name level order (not alphabetical)
    ord   <- do.call(order, lapply(seq_len(n_factors), function(f)
      match(pdata[[paste0("f", f)]], factor_levels_list[[f]])))
    pdata <- pdata[ord, , drop = FALSE]

  } else {
    pdata              <- as.data.frame(sub_arr, stringsAsFactors = FALSE)
    factor_levels_list <- list()
    factors            <- character(0)
  }

  # ---- 5. Helper: expand named/unnamed factor lists ---------------------
  # Converts a possibly-named list to a full-length list (NULL = unspecified).
  expand_flist <- function(x, arg) {
    if (is.null(x)) return(vector("list", n_factors))
    nms <- names(x)
    if (is.null(nms)) {
      if (length(x) != n_factors)
        stop("`", arg, "` (unnamed) must have length ", n_factors, ", got ", length(x), ".")
      return(x)
    }
    bad <- nms[!nms %in% factors]
    if (length(bad))
      stop("Names in `", arg, "` not matching any factor: ",
           paste(bad, collapse = ", "),
           ". Available: ", paste(factors, collapse = ", "), ".")
    out <- vector("list", n_factors)
    for (i in seq_along(nms)) out[[match(nms[i], factors)]] <- x[[i]]
    out
  }

  # ---- 6. Validate / convert factor_labels --------------------------------
  # Named-list form: treat entries as level_labels (factor_labels wins).
  if (is.list(factor_labels)) {
    for (nm in names(factor_labels))
      level_labels[[nm]] <- factor_labels[[nm]]
    factor_labels <- NULL
  }
  if (!is.null(factor_labels) && length(factor_labels) != n_factors)
    stop("`factor_labels` must have length ", n_factors, ", got ", length(factor_labels), ".")

  # ---- 7. Expand and validate level_labels ------------------------------
  ll_full <- expand_flist(level_labels, "level_labels")
  for (f in seq_len(n_factors)) {
    if (is.null(ll_full[[f]])) next
    exp <- length(factor_levels_list[[f]])
    got <- length(ll_full[[f]])
    if (got != exp) {
      nm <- if (!is.null(factor_labels)) factor_labels[f] else factors[f]
      stop("`level_labels` for factor '", nm, "' must have length ", exp,
           " matching levels: ", paste(factor_levels_list[[f]], collapse = ", "),
           ". Got length ", got, ".")
    }
  }

  display_factor_names <- if (!is.null(factor_labels)) factor_labels else factors
  display_levels_list  <- lapply(seq_len(n_factors), function(f)
    if (!is.null(ll_full[[f]])) ll_full[[f]] else factor_levels_list[[f]])

  # ---- 8. Expand, validate, and apply level_order -----------------------
  lo_full <- expand_flist(level_order, "level_order")
  for (f in seq_len(n_factors)) {
    if (is.null(lo_full[[f]])) next
    n_lev <- length(factor_levels_list[[f]])
    ord_f <- as.integer(lo_full[[f]])
    if (length(ord_f) != n_lev || !identical(sort(ord_f), seq_len(n_lev))) {
      nm <- if (!is.null(factor_labels)) factor_labels[f] else factors[f]
      stop("`level_order` for factor '", nm, "' must be a permutation of 1:",
           n_lev, " (", n_lev, " levels).")
    }
    factor_levels_list[[f]]  <- factor_levels_list[[f]][ord_f]
    display_levels_list[[f]] <- display_levels_list[[f]][ord_f]
  }

  # ---- 9. Apply factor_order --------------------------------------------
  if (!is.null(factor_order)) {
    if (n_factors == 0L)
      stop("`factor_order` cannot be used for zero-factor types.")
    fo <- as.integer(factor_order)
    if (length(fo) != n_factors || !identical(sort(fo), seq_len(n_factors)))
      stop("`factor_order` must be a permutation of 1:", n_factors,
           " (", n_factors, " factor(s)), got: ", paste(factor_order, collapse = ", "), ".")
    factor_levels_list   <- factor_levels_list[fo]
    display_levels_list  <- display_levels_list[fo]
    display_factor_names <- display_factor_names[fo]
    # Reorder factor columns in pdata and rename f1..fn accordingly
    pdata <- pdata[, c("lower", "mid", "upper", paste0("f", fo)), drop = FALSE]
    names(pdata)[4L:(3L + n_factors)] <- paste0("f", seq_len(n_factors))
    # Re-sort rows with new factor ordering (factor 1 slowest, last fastest)
    ord   <- do.call(order, lapply(seq_len(n_factors), function(f)
      match(pdata[[paste0("f", f)]], factor_levels_list[[f]])))
    pdata <- pdata[ord, , drop = FALSE]
  }

  raw_to_display <- lapply(seq_len(n_factors), function(f)
    setNames(display_levels_list[[f]], factor_levels_list[[f]]))

  # ---- 10. Extract and remove plot-control args from ... ----------------
  dots    <- list(...)
  n_lines <- if (n_factors >= 2L) length(factor_levels_list[[2L]]) else 1L

  lty  <- if (!is.null(dots$lty))  dots$lty  else seq_len(n_lines)
  pch  <- if (!is.null(dots$pch))  dots$pch  else seq_len(n_lines)
  cols <- if (!is.null(dots$cols)) dots$cols else rep("black", n_lines)

  xlab_default <- if (n_factors >= 1L) display_factor_names[1L] else type
  xlab         <- if (!is.null(dots$xlab)) dots$xlab else xlab_default
  ylab         <- if (!is.null(dots$ylab)) dots$ylab else type
  user_main    <- dots$main
  user_ylim    <- dots$ylim

  # par() settings intercepted from ...: applied after par(mfrow=) on each page
  # (par(mfrow=) resets mar/oma to defaults, so they must be re-applied after it)
  par_args <- c("mar","oma","mgp","tcl","las","cex","cex.axis","cex.lab","cex.main")
  par_settings <- dots[intersect(names(dots), par_args)]
  if (is.null(par_settings$mar)) par_settings$mar <- c(4.1, 4.1, 2.1, 0.1)
  dots[c("lty","pch","cols","xlab","ylab","main","ylim","xaxt", names(par_settings))] <- NULL

  # ---- 11. Displace ------------------------------------------------------
  # Scalar: even, symmetric spacing of the line groups around each x location
  #   (e.g. 2 lines -> +/- displace/2; 3 lines -> -displace, 0, +displace).
  # Length-n_lines vector: explicit per-line offsets.
  if (is.null(displace)) displace <- 0
  if (length(displace) == 1L) {
    if (!is.numeric(displace) || !is.finite(displace))
      stop("`displace` must be a single finite number (spacing) or a length-",
           n_lines, " vector of offsets.")
    displace <- (seq_len(n_lines) - (n_lines + 1) / 2) * displace
  } else if (length(displace) != n_lines) {
    stop("`displace` must be a single spacing value or have length ", n_lines,
         " (one per line group), got ", length(displace), ".")
  }

  # ---- 11. X-axis positions ---------------------------------------------
  if (n_factors >= 1L) {
    x_levels_raw  <- factor_levels_list[[1L]]
    x_levels_disp <- display_levels_list[[1L]]
    n_x           <- length(x_levels_raw)
  } else {
    x_levels_disp <- type
    n_x           <- 1L
  }
  x_base <- seq_len(n_x)
  cap_hw <- cap_length

  # ---- 12. Panel structure ----------------------------------------------
  if (n_factors <= 2L) {
    panel_combos <- NULL
    n_panels     <- 1L
    extra_f      <- integer(0)
  } else {
    extra_f      <- seq(3L, n_factors)
    panel_combos <- expand.grid(
      lapply(extra_f, function(f) factor_levels_list[[f]]),
      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    names(panel_combos) <- paste0("f", extra_f)
    n_panels <- nrow(panel_combos)
  }

  # `main` may be length 1 (shared by every panel, as before) or length n_panels
  # (one title per panel).
  if (!is.null(user_main) && length(user_main) != 1L &&
      length(user_main) != n_panels)
    stop("`main` must have length 1 or ", n_panels,
         " (one per panel), got ", length(user_main), ".")

  # `ylim` may be NULL (auto per panel), a length-2 vector (shared by every
  # panel), or a list of length-2 vectors (one per panel).
  if (!is.null(user_ylim)) {
    if (is.list(user_ylim)) {
      if (length(user_ylim) != n_panels)
        stop("`ylim` list must have length ", n_panels,
             " (one per panel), got ", length(user_ylim), ".")
      if (!all(vapply(user_ylim,
                      function(z) is.numeric(z) && length(z) == 2L, logical(1L))))
        stop("each element of a `ylim` list must be a numeric vector of length 2.")
    } else if (!is.numeric(user_ylim) || length(user_ylim) != 2L) {
      stop("`ylim` must be a length-2 numeric vector, or a list of such ",
           "vectors (one per panel).")
    }
  }

  # Default (no user `ylim`): a single common y range spanning the intervals of
  # every panel of this parameter type, so the panels are directly comparable.
  auto_ylim <- {
    yl  <- range(c(pdata$lower, pdata$upper), na.rm = TRUE)
    pad <- 0.05 * diff(yl); if (pad == 0) pad <- 0.5
    yl + c(-pad, pad)
  }

  # layout: NA = compute & apply default; NULL = don't touch mfrow; c(r,c) = use as-is
  set_layout <- !is.null(layout)
  if (isTRUE(is.na(layout))) {
    if (n_factors == 3L) {
      nr <- 1L
      nc <- length(factor_levels_list[[3L]])
    } else if (n_factors >= 4L) {
      nr <- 2L
      nc <- 4L
    } else {
      nr <- 1L; nc <- 1L
    }
    layout <- c(nr, nc)
  }
  # layout = NULL: don't touch mfrow; derive page size from the current device
  # so the page-boundary bookkeeping below still works.
  panels_per_page <- if (is.null(layout)) prod(par("mfrow")) else layout[1L] * layout[2L]

  # ---- 13. Draw helper --------------------------------------------------
  draw_group <- function(x_pos, grp_data, j) {
    if (n_factors >= 1L) {
      ord      <- match(grp_data$f1, x_levels_raw)
      grp_data <- grp_data[order(ord), , drop = FALSE]
    }
    lines(x_pos,  grp_data$mid, lty = lty[j], col = cols[j])
    points(x_pos, grp_data$mid, pch = pch[j], col = cols[j], cex = pt_cex)
    for (xi in seq_along(x_pos)) {
      x  <- x_pos[xi];  lo <- grp_data$lower[xi];  hi <- grp_data$upper[xi]
      segments(x,          lo, x,          hi, col = cols[j])
      segments(x - cap_hw, lo, x + cap_hw, lo, col = cols[j])
      segments(x - cap_hw, hi, x + cap_hw, hi, col = cols[j])
    }
  }

  # ---- 14. Plot loop ----------------------------------------------------
  # Save and restore only what this function explicitly changes:
  # mfrow (if set_layout) and any par_settings supplied by the user.
  # Saving the full par() snapshot would restore figure geometry (mfrow, fig,
  # fin, plt...) and reset the panel counter, breaking sequential multi-panel use.
  pars_to_save <- names(par_settings)
  if (set_layout) pars_to_save <- c("mfrow", pars_to_save)
  if (length(pars_to_save)) {
    old_par <- setNames(lapply(pars_to_save, par), pars_to_save)
    on.exit(par(old_par), add = TRUE)
  }

  for (panel_idx in seq_len(n_panels)) {

    if ((panel_idx - 1L) %% panels_per_page == 0L) {
      if (set_layout) par(mfrow = layout)
      if (length(par_settings)) do.call(par, par_settings)
    }

    # length-1 main is shared by all panels; length-n_panels gives one each
    um <- if (is.null(user_main)) NULL
          else if (length(user_main) == 1L) user_main else user_main[panel_idx]

    if (is.null(panel_combos)) {
      pdata_panel <- pdata
      panel_title <- if (!is.null(um)) um else ""
    } else {
      combo  <- panel_combos[panel_idx, , drop = FALSE]
      filter <- rep(TRUE, nrow(pdata))
      for (k in seq_along(extra_f)) {
        f      <- extra_f[k]
        filter <- filter & (pdata[[paste0("f", f)]] == combo[[k]])
      }
      pdata_panel <- pdata[filter, , drop = FALSE]
      panel_title <- if (!is.null(um)) um else
        paste(sapply(seq_along(extra_f), function(k) {
          f <- extra_f[k]
          paste0(display_factor_names[f], ": ", raw_to_display[[f]][combo[[k]]])
        }), collapse = ", ")
    }

    ylim <-
      if (is.null(user_ylim)) auto_ylim
      else if (is.list(user_ylim)) user_ylim[[panel_idx]]
      else user_ylim

    do.call(plot.default,
            c(list(x = numeric(0), y = numeric(0),
                   xlim = c(0.5, n_x + 0.5), ylim = ylim,
                   xaxt = "n", bty = "l",
                   xlab = xlab, ylab = ylab, main = panel_title),
              dots))
    axis(1, at = x_base, labels = x_levels_disp)

    if (n_factors == 0L) {
      segments(1,          pdata_panel$lower, 1,          pdata_panel$upper, col = cols[1])
      segments(1 - cap_hw, pdata_panel$lower, 1 + cap_hw, pdata_panel$lower, col = cols[1])
      segments(1 - cap_hw, pdata_panel$upper, 1 + cap_hw, pdata_panel$upper, col = cols[1])
      points(1, pdata_panel$mid, pch = pch[1], col = cols[1], cex = pt_cex)
    } else if (n_factors == 1L) {
      draw_group(x_base + displace[1L], pdata_panel, 1L)
    } else {
      for (j in seq_len(n_lines)) {
        raw_f2 <- factor_levels_list[[2L]][j]
        draw_group(x_base + displace[j],
                   pdata_panel[pdata_panel$f2 == raw_f2, , drop = FALSE], j)
      }
      if (n_lines > 1L) {
        leg_defaults <- list(x     = "topleft",
                             bty   = "n",
                             title = display_factor_names[2L],
                             cex   = 0.75)
        leg_call <- modifyList(leg_defaults,
                               if (!is.null(legend_args)) legend_args else list())
        if (!isTRUE(is.na(leg_call$x)))
          do.call(legend, c(list(legend = display_levels_list[[2L]],
                                 lty = lty, pch = pch, col = cols),
                            leg_call))
      }
    }
  }

  invisible(pdata)
}

#' Plot credible intervals from credint(map=FALSE)
#'
#' Displays posterior credible intervals for one parameter type from
#' \code{credint(..., map=FALSE)} on a two-scale plot: intercept parameters
#' on the left y-axis and effect parameters on the right y-axis, separated
#' by an extra gap on the shared x-axis.
#'
#' @param ci Length-1 list returned by \code{credint(..., map=FALSE)}.
#' @param type Character string selecting the parameter type: the prefix of
#'   each row name before the first \code{_}, or the full name if there is no
#'   \code{_}.  If \code{NULL} (default) the first type found is used and a
#'   message identifies it.
#' @param quants Integer vector of length 3, strictly increasing, indexing the
#'   quantile columns (default \code{1:3}: lower, middle, upper).
#' @param intercept Integer vector indicating which rows (within the selected
#'   type, in the order they appear in \code{ci[[1]]}) are intercept
#'   parameters.  Default \code{1} (first row only).  \code{NULL} or
#'   \code{integer(0)} produces an effects-only plot (no left axis).  If all
#'   rows are listed the plot is intercept-only (no right axis).  Ignored when
#'   \code{effects} is supplied.
#' @param effects Integer vector indicating which rows are \emph{effect}
#'   parameters; \code{intercept} is derived as the complement.  \code{NULL}
#'   (explicit) produces an intercept-only plot (no right axis).  If
#'   \code{effects} covers every row the plot has no intercepts and no left
#'   axis.  When not supplied, \code{intercept} determines the split.
#' @param yleft Named list controlling the left (intercept) y-axis.  All
#'   entries except \code{label} and \code{lim} are passed directly to
#'   \code{\link[graphics]{axis}(side = 2)}.  Recognised special keys:
#'   \describe{
#'     \item{\code{label}}{Y-axis title.  Takes precedence over
#'       \code{ylabs[1]} and the default \code{"Intercept"}.}
#'     \item{\code{lim}}{Numeric \code{c(lower, upper)} fixing the left axis
#'       limits.  Default \code{NULL} scales to fit the data with minimal
#'       padding.}
#'   }
#'   Other keys: \code{at}, \code{labels}, \code{cex.axis}, \code{las},
#'   \code{lwd}, \code{lwd.ticks}, \code{col}, \code{tcl}, etc.
#' @param yright Named list controlling the right (effect) y-axis, following
#'   the same structure as \code{yleft}.  \code{label} overrides
#'   \code{ylabs[2]} and the default \code{"Effect"}; \code{lim} fixes the
#'   right axis limits.
#' @param intercept_names Character vector of x-axis labels for the intercept
#'   parameters, in order.  Default \code{NULL} uses the full row name (e.g.
#'   \code{"m"}).
#' @param effect_names Character vector of x-axis labels for the effect
#'   parameters, in order.  Default \code{NULL} strips the type prefix (e.g.
#'   \code{"m_lMd"} becomes \code{"lMd"}).
#' @param layout Integer vector \code{c(rows, cols)} passed to
#'   \code{par(mfrow = layout)}.  Default \code{c(1, 1)} (single panel).
#'   Set to \code{NULL} to leave the current layout unchanged.
#' @param cap_length Numeric. Half-width of the interval end-caps as a
#'   fraction of the x spacing between adjacent parameters (default
#'   \code{0.05}).
#' @param pt_cex Numeric. Character expansion factor for plotted points,
#'   passed as \code{cex} to \code{points()}.  Default \code{1}.
#' @param effect_line Logical (default \code{TRUE}).  When \code{TRUE} and
#'   effect parameters are present, a horizontal grey dotted line is drawn at
#'   \code{y = 0} in the effect coordinate system, spanning the x range of the
#'   effect parameters.  Useful as a visual reference for null effects.
#' @param label_angle Numeric vector controlling x-axis label rotation in
#'   degrees.  Accepted lengths:
#'   \itemize{
#'     \item \strong{1} -- same angle for every label (default \code{90},
#'       vertical).
#'     \item \strong{2} -- \code{label_angle[1]} for all intercept labels,
#'       \code{label_angle[2]} for all effect labels.
#'     \item \strong{n_int + n_eff} -- one angle
#'       per label individually.
#'   }
#'   Text alignment (\code{adj}) is derived automatically from the angle:
#'   centred at 0 degrees, right-aligned at 90 degrees.  When non-default angles produce
#'   long labels, increase the bottom margin via \code{mar} in \code{...}.
#' @param xtick_stagger Logical.  If \code{TRUE}, every second x-axis tick label
#'   is dropped onto a lower line so that horizontal labels
#'   (\code{label_angle = 0}) do not overlap.  Default \code{FALSE}.  With many
#'   or long labels you may also need a larger bottom margin via \code{mar} in
#'   \code{...}.
#' @param ... Additional graphical arguments. The following are intercepted
#'   with defaults:
#'   \describe{
#'     \item{\code{main}}{Plot title. Default: \code{names(ci)[1]}, a colon,
#'       then the name of the first intercept parameter.}
#'     \item{\code{xlab}}{X-axis label. Default: \code{""}.}
#'     \item{\code{ylab}}{Left y-axis label. Default: \code{"Intercept"}.}
#'     \item{\code{ylab_right}}{Right y-axis label. Default: \code{"Effect"}.}
#'     \item{\code{ylim}}{Length-2 numeric limits for the left (intercept)
#'       y-axis.  Default: computed from the intercept intervals.  Equivalent
#'       to \code{yleft = list(lim = ...)}, which takes precedence if both are
#'       given.  The right (effect) axis is set via \code{yright}.}
#'     \item{\code{mar}}{Plot margins. Default: \code{c(4.1, 4.1, 2.1, 4.1)}
#'       (extra right margin for the right axis).}
#'     \item{\code{col_intercept}}{Colour for intercept points and intervals.
#'       Default: \code{"black"}.}
#'     \item{\code{col_effect}}{Colour for effect points and intervals.
#'       Default: \code{"black"}.}
#'     \item{\code{pch}}{Point character. Default: \code{16}.}
#'   }
#'
#' @details
#' The x-axis has \code{n_intercept} positions followed by a half-space gap
#' then \code{n_effect} positions.  The two groups of parameters are plotted in
#' their row-name order within \code{ci[[1]]}.  Points show the middle quantile;
#' vertical bars extend to the lower and upper quantiles with horizontal end-caps.
#'
#' @return Invisibly returns the filtered matrix of quantiles for the selected
#'   type.
#' @seealso \code{\link{plot_credint_map}} for \code{credint(map=TRUE)} output.
#' @examples \donttest{
#' ci <- credint(samples_LNR, map = FALSE)
#' plot_credint(ci, type = "m")
#' plot_credint(ci, type = "m",
#'              intercept_names = "Mean",
#'              effect_names    = "lM effect",
#'              col_intercept   = "steelblue",
#'              col_effect      = "tomato")
#' }
#' @export
plot_credint <- function(ci, type = NULL, quants = 1:3,
                         intercept = 1L, effects = NULL,
                         yleft = NULL, yright = NULL,
                         intercept_names = NULL, effect_names = NULL,
                         layout = NA,
                         cap_length = 0.05, pt_cex = 1,
                         effect_line = TRUE, label_angle = 90,
                         xtick_stagger = FALSE, ...) {

  effects_given <- !missing(effects)
  # layout: NA = c(1,1) for this single-panel function; NULL = don't touch mfrow
  if (isTRUE(is.na(layout))) layout <- c(1L, 1L)

  # ---- 1. Validate ci and quants ------------------------------------------
  if (!is.list(ci) || !all(vapply(ci, is.matrix, logical(1L))))
    stop("`ci` must be a list of matrices (output of credint()).")
  if (length(ci) > 1L)
    stop("For multi-element `ci` (e.g. correlations/covariances) use `plot_credints()`.")
  arr <- ci[[1]]
  if (ncol(arr) < 3)
    stop("`ci[[1]]` must have at least 3 quantile columns; found ", ncol(arr), ".")
  if (length(quants) != 3 || !all(diff(quants) > 0))
    stop("`quants` must be a strictly increasing integer vector of length 3.")
  if (any(quants < 1L) || any(quants > ncol(arr)))
    stop("`quants` indices out of range 1:", ncol(arr), ".")

  # ---- 2. Resolve type and filter rows ------------------------------------
  rn        <- rownames(arr)
  row_types <- sub("_.*", "", rn)
  if (is.null(type)) {
    type <- unique(row_types)[1L]
    message("type not supplied; using '", type, "'.")
  }

  # ---- Vector type: intercept-only multi-type plot -------------------------
  if (length(type) > 1L) {
    sel_rows <- vapply(type, function(t) {
      idx <- which(row_types == t)
      if (length(idx) == 0L)
        stop("Type '", t, "' not found. Available: ",
             paste(sort(unique(row_types)), collapse = ", "), ".")
      idx[1L]
    }, integer(1L))
    sub_arr           <- arr[sel_rows, quants, drop = FALSE]
    colnames(sub_arr) <- c("lower", "mid", "upper")
    n_pars            <- nrow(sub_arr)
    x_pos             <- seq_len(n_pars)
    x_range           <- c(0.5, n_pars + 0.5)

    x_labels <- if (!is.null(intercept_names)) {
      if (length(intercept_names) != n_pars)
        stop("`intercept_names` must have length ", n_pars,
             " (one per type), got ", length(intercept_names), ".")
      intercept_names
    } else rownames(sub_arr)

    if (length(label_angle) == 1L) {
      angles <- rep(label_angle, n_pars)
    } else if (length(label_angle) == n_pars) {
      angles <- label_angle
    } else {
      stop("`label_angle` must have length 1 or ", n_pars,
           " (one per type), got ", length(label_angle), ".")
    }
    adjs <- cbind(sin(angles * pi / 180) * 0.5 + 0.5, 1)

    dots          <- list(...)
    par_arg_names <- c("mar","oma","mgp","tcl","las","cex",
                       "cex.axis","cex.lab","cex.main")
    par_settings  <- dots[intersect(names(dots), par_arg_names)]
    if (is.null(par_settings$mar)) par_settings$mar <- c(4.1, 4.1, 2.1, 1.1)

    main <- if (!is.null(dots$main)) dots$main else ""
    xlab <- if (!is.null(dots$xlab)) dots$xlab else ""
    ylab <- if (!is.null(dots$ylab)) dots$ylab else "Intercept"
    col  <- if (!is.null(dots$col_intercept)) dots$col_intercept else "black"
    pch  <- if (!is.null(dots$pch)) dots$pch else 16L

    if (!is.null(yleft$label)) ylab <- yleft$label

    pad  <- function(r) { d <- diff(r); if (d == 0) d <- 0.5; r + c(-1,1)*0.05*d }
    ylim <- if (!is.null(yleft$lim)) yleft$lim else
            if (!is.null(dots$ylim)) dots$ylim else
              pad(range(sub_arr[, c("lower","upper")], na.rm = TRUE))

    pars_to_save <- names(par_settings)
    if (!is.null(layout)) pars_to_save <- c("mfrow", pars_to_save)
    if (length(pars_to_save)) {
      old_par <- setNames(lapply(pars_to_save, par), pars_to_save)
      on.exit(par(old_par), add = TRUE)
    }
    if (!is.null(layout)) par(mfrow = layout)
    if (length(par_settings)) do.call(par, par_settings)

    plot.default(NA, xlim = x_range, ylim = ylim,
                 axes = FALSE, xlab = xlab, ylab = ylab, main = main)
    axis(2L)
    cap_hw <- cap_length
    for (i in seq_len(n_pars)) {
      x  <- x_pos[i]
      lo <- sub_arr[i,"lower"]; hi <- sub_arr[i,"upper"]; md <- sub_arr[i,"mid"]
      segments(x, lo, x, hi, col = col)
      segments(x - cap_hw, lo, x + cap_hw, lo, col = col)
      segments(x - cap_hw, hi, x + cap_hw, hi, col = col)
      points(x, md, pch = pch, col = col, cex = pt_cex)
    }

    y0          <- par("usr")[3L]
    line_to_usr <- diff(par("usr")[3L:4L]) / par("pin")[2L] *
                   par("cin")[2L] * par("cex") * par("mgp")[2L]
    if (n_pars == 1L)
      segments(x_pos - 0.25, y0, x_pos + 0.25, y0, xpd = TRUE)
    else
      axis(1L, at = x_pos, labels = FALSE)

    cex_axis <- par("cex.axis")
    y_text   <- y0 - line_to_usr * cex_axis
    step     <- line_to_usr * cex_axis   # one text line, for staggered labels
    for (i in seq_len(n_pars)) {
      off <- if (xtick_stagger && i %% 2L == 0L) step else 0
      text(x_pos[i], y_text - off, labels = x_labels[i],
           srt = angles[i], xpd = TRUE, adj = adjs[i, ], cex = cex_axis)
    }

    return(invisible(sub_arr))
  }
  # --------------------------------------------------------------------------

  sel <- row_types == type
  if (!any(sel))
    stop("Type '", type, "' not found. Available: ",
         paste(sort(unique(row_types)), collapse = ", "), ".")
  sub_arr           <- arr[sel, quants, drop = FALSE]
  colnames(sub_arr) <- c("lower", "mid", "upper")
  n_rows            <- nrow(sub_arr)

  # ---- 3. Resolve intercept / effects -------------------------------------
  if (effects_given) {
    if (is.null(effects)) {
      intercept <- seq_len(n_rows)          # effects=NULL -> intercept-only
    } else {
      effects <- as.integer(effects)
      if (any(effects < 1L) || any(effects > n_rows))
        stop("`effects` indices out of range 1:", n_rows, ".")
      intercept <- as.integer(setdiff(seq_len(n_rows), effects))
    }
  } else {
    if (is.null(intercept)) intercept <- integer(0L)  # intercept=NULL -> no intercepts
    intercept <- as.integer(intercept)
    if (length(intercept) > 0L && (any(intercept < 1L) || any(intercept > n_rows)))
      stop("`intercept` indices out of range 1:", n_rows, ".")
  }
  effect_idx <- setdiff(seq_len(n_rows), intercept)
  int_arr       <- sub_arr[intercept,  , drop = FALSE]
  eff_arr       <- sub_arr[effect_idx, , drop = FALSE]

  # ---- 4. X positions and labels ------------------------------------------
  n_int   <- nrow(int_arr)
  n_eff   <- nrow(eff_arr)
  x_int   <- seq_len(n_int)
  # half-space gap only when both groups present; effects start at 1 when alone
  x_eff   <- if (n_int > 0L) seq_len(n_eff) + n_int else seq_len(n_eff)
  x_range <- c(0.5, max(c(x_int, x_eff), 0) + 0.5)

  # Default label for intercepts: full row name
  if (is.null(intercept_names))
    intercept_names <- rownames(int_arr)
  # Default label for effects: strip "type_" prefix
  if (is.null(effect_names))
    effect_names <- sub(paste0("^", type, "_?"), "", rownames(eff_arr))

  if (length(intercept_names) != n_int)
    stop("`intercept_names` must have length ", n_int, ".")
  if (length(effect_names) != n_eff)
    stop("`effect_names` must have length ", n_eff, ".")

  # ---- 5. Extract ... args ------------------------------------------------
  dots <- list(...)
  par_arg_names <- c("mar","oma","mgp","tcl","las","cex","cex.axis","cex.lab","cex.main")
  par_settings  <- dots[intersect(names(dots), par_arg_names)]
  if (is.null(par_settings$mar)) {
    mar <- c(4.1, 4.1, 2.1, 4.1)
    if (n_int == 0L) mar[2L] <- 1.1   # no left axis  -> shrink left margin
    if (n_eff == 0L) mar[4L] <- 1.1   # no right axis -> shrink right margin
    par_settings$mar <- mar
  }

  main       <- if (!is.null(dots$main))       dots$main else
                  if (n_int > 0L) paste0(names(ci)[1L], ": ", intercept_names[1L])
                  else paste0(names(ci)[1L], ": ", type)
  xlab       <- if (!is.null(dots$xlab))       dots$xlab       else ""
  ylab       <- if (!is.null(dots$ylab))       dots$ylab       else if (n_int > 0L) "Intercept" else ""
  ylab_right <- if (!is.null(dots$ylab_right)) dots$ylab_right else "Effect"
  if (!is.null(yleft$label))  ylab       <- yleft$label
  if (!is.null(yright$label)) ylab_right <- yright$label
  col_int    <- if (!is.null(dots$col_intercept)) dots$col_intercept else "black"
  col_eff    <- if (!is.null(dots$col_effect))    dots$col_effect    else "black"
  pch        <- if (!is.null(dots$pch)) dots$pch else 16L

  # ---- 6. Y limits --------------------------------------------------------
  pad <- function(r) { d <- diff(r); if (d == 0) d <- 0.5; r + c(-1,1)*0.05*d }
  intercept_ylim <- if (!is.null(yleft$lim))  yleft$lim else
                    if (!is.null(dots$ylim)) dots$ylim else
                    if (n_int > 0L) pad(range(int_arr[, c("lower","upper")], na.rm = TRUE))
  effect_ylim    <- if (!is.null(yright$lim)) yright$lim else
                    if (n_eff > 0L) pad(range(eff_arr[, c("lower","upper")], na.rm = TRUE))

  cap_hw <- cap_length

  # ---- 7. Draw helper (intervals + point) ---------------------------------
  draw_ci <- function(x_pos, dat, col) {
    for (i in seq_along(x_pos)) {
      x  <- x_pos[i]
      lo <- dat[i, "lower"];  hi <- dat[i, "upper"];  md <- dat[i, "mid"]
      segments(x, lo, x, hi, col = col)
      segments(x - cap_hw, lo, x + cap_hw, lo, col = col)
      segments(x - cap_hw, hi, x + cap_hw, hi, col = col)
      points(x, md, pch = pch, col = col, cex = pt_cex)
    }
  }

  # ---- 8. Resolve label_angle to per-label vector -------------------------
  x_all      <- c(x_int, x_eff)
  all_labels <- c(intercept_names, effect_names)
  n_labels   <- length(all_labels)

  if (length(label_angle) == 1L) {
    angles <- rep(label_angle, n_labels)
  } else if (length(label_angle) == 2L) {
    angles <- c(rep(label_angle[1L], n_int), rep(label_angle[2L], n_eff))
  } else if (length(label_angle) == n_labels) {
    angles <- label_angle
  } else {
    stop("`label_angle` must have length 1, 2, or ", n_labels,
         " (one per label), got ", length(label_angle), ".")
  }
  # adj[1]: 0.5 (centred) at 0 deg, 1 (right-aligned) at 90 deg
  adjs <- cbind(sin(angles * pi / 180) * 0.5 + 0.5, 1)

  # ---- 9. Plot ------------------------------------------------------------
  # Save and restore only what this function explicitly changes.
  pars_to_save <- names(par_settings)
  if (!is.null(layout)) pars_to_save <- c("mfrow", pars_to_save)
  if (length(pars_to_save)) {
    old_par <- setNames(lapply(pars_to_save, par), pars_to_save)
    on.exit(par(old_par), add = TRUE)
  }
  if (!is.null(layout)) par(mfrow = layout)
  if (length(par_settings)) do.call(par, par_settings)
  ylab_cex <- par("cex.lab")   # capture before overlay may reset it

  # Base plot -- use intercept scale when available, else effect scale
  base_ylim <- if (n_int > 0L) intercept_ylim else effect_ylim
  plot.default(NA, xlim = x_range, ylim = base_ylim,
               axes = FALSE, xlab = xlab, ylab = ylab, main = main)

  left_axis_args  <- yleft[!names(yleft)  %in% c("label", "lim")]
  right_axis_args <- yright[!names(yright) %in% c("label", "lim")]

  if (n_int > 0L) {
    do.call(axis, c(list(side = 2L), left_axis_args))
    draw_ci(x_int, int_arr, col_int)
  }

  # Right-axis overlay (effect scale) -- also acts as sole axis when n_int == 0
  if (n_eff > 0L) {
    if (n_int > 0L) par(new = TRUE)
    if (n_int > 0L)
      plot.default(NA, xlim = x_range, ylim = effect_ylim,
                   axes = FALSE, ann = FALSE, bty = "n")
    if (effect_line)
      segments(x_eff[1L] - 0.25, 0, x_eff[length(x_eff)] + 0.25, 0,
               col = "grey", lty = 3L)
    draw_ci(x_eff, eff_arr, col_eff)
    do.call(axis, c(list(side = 4L), right_axis_args))
    if (nchar(ylab_right) > 0L)
      # mtext() cex is absolute, but title() (left ylab) scales by par("cex"),
      # which shrinks in multi-panel layouts; match it so both labels agree.
      mtext(ylab_right, side = 4L, line = 3L, cex = ylab_cex * par("cex"))
  }

  # X-axis: full-width line (with break if both groups), then ticks only
  y0    <- par("usr")[3L]
  line_to_usr <- diff(par("usr")[3L:4L]) / par("pin")[2L] *
                 par("cin")[2L] * par("cex") * par("mgp")[2L]

  # Two disconnected axis segments, one per group.
  # Multi-value: axis() naturally spans first-to-last tick with tick marks.
  # Single-value: short line only (a tick with no extent looks odd).
  draw_xaxis_seg <- function(x_pos) {
    if (length(x_pos) == 1L)
      segments(x_pos - 0.25, y0, x_pos + 0.25, y0, xpd = TRUE)
    else
      axis(1L, at = x_pos, labels = FALSE)
  }
  if (n_int > 0L) draw_xaxis_seg(x_int)
  if (n_eff > 0L) draw_xaxis_seg(x_eff)

  cex_axis <- par("cex.axis")   # picks up cex.axis from par_settings if supplied
  y_text <- y0 - line_to_usr * cex_axis
  step   <- line_to_usr * cex_axis   # one text line, for staggered labels
  for (i in seq_along(x_all)) {
    off <- if (xtick_stagger && i %% 2L == 0L) step else 0
    text(x_all[i], y_text - off, labels = all_labels[i],
         srt = angles[i], xpd = TRUE,
         adj = adjs[i, ], cex = cex_axis)
  }

  invisible(sub_arr)
}

# ---- Shared helpers for plot_credints / plot_credints_map -------------------

.classify_types <- function(ci) {
  rn     <- rownames(ci[[1]])
  types  <- sub("_.*", "", rn)
  counts <- table(types)
  # preserve row-name order
  ord    <- match(names(counts), unique(types))
  counts <- counts[order(ord)]
  list(
    all_types    = names(counts),
    single_types = names(counts[counts == 1L]),
    multi_types  = names(counts[counts > 1L]),
    counts       = counts,
    rn           = rn,
    types        = types
  )
}

.resolve_order <- function(order, multi_types, drop = TRUE) {
  if (is.null(order)) return(multi_types)
  if (is.character(order)) {
    bad <- order[!order %in% multi_types]
    if (length(bad))
      stop("`order` contains types not among multi-value types: ",
           paste(bad, collapse = ", "), ".")
    return(order)
  }
  order <- as.integer(order)
  if (any(order < 1L) || any(order > length(multi_types)))
    stop("`order` indices out of range 1:", length(multi_types), ".")
  multi_types[order]
}

.panel_call <- function(fn, base_args, type_key, plot_args, dots,
                        main_override, sub_layout = NULL) {
  args <- modifyList(dots, if (!is.null(plot_args[[type_key]]))
                             plot_args[[type_key]] else list())
  args$main <- if (!is.null(plot_args[[type_key]]$main))
                 plot_args[[type_key]]$main else main_override
  # Resolve layout: plot_args wins, then sub_layout.
  # Must use c(..., list(layout=val)) rather than args$layout <- val because
  # $<- with NULL is a no-op in R and would lose layout=NULL.
  layout_val <- if (!is.null(plot_args[[type_key]]$layout))
                  plot_args[[type_key]]$layout else sub_layout
  args <- args[names(args) != "layout"]    # remove any stale layout key
  do.call(fn, c(base_args, list(layout = layout_val), args))
}

# Internal helper: single correlation/covariance panel (no @export -- not user-facing)
.plot_cor_panel <- function(mat, quants = 1:3, type = NULL,
                             ylab = "Correlation", main = "",
                             effect_line = TRUE, cap_length = 0.05,
                             pt_cex = 1, layout = NA, ...) {
  sub_mat           <- mat[, quants, drop = FALSE]
  colnames(sub_mat) <- c("lower", "mid", "upper")
  rn                <- rownames(sub_mat)
  if (!is.null(type)) {
    sel <- sub("_.*", "", rn) %in% type
    if (!any(sel)) stop("None of the specified types found in this element.")
    sub_mat <- sub_mat[sel, , drop = FALSE]
    rn      <- rownames(sub_mat)
  }
  ord     <- order(sub_mat[, "mid"], decreasing = TRUE)
  sub_mat <- sub_mat[ord, , drop = FALSE]
  rn      <- rownames(sub_mat)
  n       <- nrow(sub_mat)
  if (n == 0L) return(invisible(NULL))
  x_pos   <- seq_len(n)
  x_range <- c(0.5, n + 0.5)
  cap_hw  <- cap_length
  use_alt <- n > 3L

  dots          <- list(...)
  par_arg_names <- c("mar","oma","mgp","tcl","las","cex","cex.axis","cex.lab","cex.main")
  par_settings  <- dots[intersect(names(dots), par_arg_names)]

  xlab    <- if (!is.null(dots$xlab)) dots$xlab else ""
  col     <- if (!is.null(dots$col))  dots$col  else "black"
  pch     <- if (!is.null(dots$pch))  dots$pch  else 16L
  las_val <- if (!is.null(dots$las))  dots$las  else 0L

  if (is.null(par_settings$mar))
    par_settings$mar <- if (!use_alt) c(4.1, 4.1, 2.1, 1.1) else
                        if (las_val == 0L) c(6.1, 4.1, 6.1, 1.1) else c(4.1, 4.1, 4.1, 1.1)

  pad  <- function(r) { d <- diff(r); if (d == 0) d <- 0.5; r + c(-1,1)*0.05*d }
  ylim <- if (!is.null(dots$ylim)) dots$ylim else
            pad(range(sub_mat[, c("lower","upper")], na.rm = TRUE))

  if (isTRUE(is.na(layout))) layout <- c(1L, 1L)
  pars_to_save <- names(par_settings)
  if (!is.null(layout)) pars_to_save <- c("mfrow", pars_to_save)
  if (length(pars_to_save)) {
    old_par <- setNames(lapply(pars_to_save, par), pars_to_save)
    on.exit(par(old_par), add = TRUE)
  }
  if (!is.null(layout)) par(mfrow = layout)
  if (length(par_settings)) do.call(par, par_settings)

  plot.default(NA, xlim = x_range, ylim = ylim,
               axes = FALSE, xlab = xlab, ylab = ylab, main = main)
  axis(2L)

  if (effect_line)
    segments(x_range[1L], 0, x_range[2L], 0, col = "grey", lty = 3L)

  for (i in seq_len(n)) {
    lo <- sub_mat[i,"lower"]; hi <- sub_mat[i,"upper"]; md <- sub_mat[i,"mid"]
    segments(i, lo, i, hi, col = col)
    segments(i - cap_hw, lo, i + cap_hw, lo, col = col)
    segments(i - cap_hw, hi, i + cap_hw, hi, col = col)
    points(i, md, pch = pch, col = col, cex = pt_cex)
  }

  if (use_alt) {
    odd_idx  <- seq(1L, n, by = 2L)
    even_idx <- seq(2L, n, by = 2L)
    # Full-width axis lines
    axis(1L, at = c(1L, n), labels = FALSE, lwd = 1L, lwd.ticks = 0L)
    axis(3L, at = c(1L, n), labels = FALSE, lwd = 1L, lwd.ticks = 0L)
    # Ticks (all positions, both axes, no labels here)
    axis(1L, at = odd_idx,  labels = FALSE, lwd = 0L, lwd.ticks = 1L)
    axis(3L, at = even_idx, labels = FALSE, lwd = 0L, lwd.ticks = 1L)
    if (length(even_idx)) axis(1L, at = even_idx, labels = FALSE, lwd = 0L, lwd.ticks = 1L)
    if (length(odd_idx))  axis(3L, at = odd_idx,  labels = FALSE, lwd = 0L, lwd.ticks = 1L)
    # Labels: staggered mtext when las=0 (default), axis() otherwise
    if (las_val == 0L) {
      line_lo   <- 0.5
      line_hi   <- 2.0
      mtext_cex <- par("cex.axis")
      odd_lo    <- odd_idx[seq(1L, length(odd_idx),  by = 2L)]
      odd_hi    <- odd_idx[seq(2L, length(odd_idx),  by = 2L)]
      even_lo   <- even_idx[seq(1L, length(even_idx), by = 2L)]
      even_hi   <- even_idx[seq(2L, length(even_idx), by = 2L)]
      if (length(odd_lo))  mtext(rn[odd_lo],  side=1L, at=odd_lo,  line=line_lo, las=0L, cex=mtext_cex, col=col)
      if (length(odd_hi))  mtext(rn[odd_hi],  side=1L, at=odd_hi,  line=line_hi, las=0L, cex=mtext_cex, col=col)
      if (length(even_lo)) mtext(rn[even_lo], side=3L, at=even_lo, line=line_lo, las=0L, cex=mtext_cex, col=col)
      if (length(even_hi)) mtext(rn[even_hi], side=3L, at=even_hi, line=line_hi, las=0L, cex=mtext_cex, col=col)
    } else {
      axis(1L, at = odd_idx,  labels = rn[odd_idx],  las = las_val, lwd = 0L, lwd.ticks = 0L)
      axis(3L, at = even_idx, labels = rn[even_idx], las = las_val, lwd = 0L, lwd.ticks = 0L)
    }
  } else {
    axis(1L, at = x_pos, labels = rn, las = las_val)
  }
  invisible(sub_mat)
}

#' Plot all credint(map=FALSE) types in a multi-panel display
#'
#' Automatically classifies parameter types from \code{credint(..., map=FALSE)}
#' output, aggregates single-value types into one shared-scale panel, and
#' produces a standard \code{\link{plot_credint}} panel for each multi-value
#' type.
#'
#' @param ci Length-1 list returned by \code{credint(..., map=FALSE)}.
#' @param quants Integer vector of length 3 passed to each
#'   \code{\link{plot_credint}} call (default \code{1:3}).
#' @param effect_only Character vector of type names to be plotted as
#'   effects-only (no left axis).  Single-value types listed here get their
#'   own panel rather than being aggregated in panel 1.
#' @param intercept_only Character vector of type names to be plotted as
#'   intercept-only (no right axis).
#' @param order \code{NULL} (default, show all multi-value types in row-name
#'   order) or a character / integer vector selecting and ordering the
#'   multi-value type panels.  Unlisted types are dropped.
#' @param layout \code{NA} (default) computes a layout of at most 4 columns;
#'   \code{NULL} leaves the current \code{par(mfrow)} untouched; a
#'   \code{c(rows, cols)} vector applies that layout.
#' @param plot_args Named list of per-panel argument lists.  Keys are type
#'   names (as returned by the type-splitting logic) or \code{"intercepts"}
#'   for the aggregate single-value panel.  Each element is passed to the
#'   corresponding \code{plot_credint} call and overrides anything in
#'   \code{...}.
#' @param ... Additional arguments forwarded to every \code{plot_credint}
#'   call.
#'
#' @return Invisibly returns \code{NULL}.
#' @seealso \code{\link{plot_credint}}, \code{\link{plot_credints_map}}
#' @export
plot_credints <- function(ci, quants = 1:3,
                          effect_only = NULL, intercept_only = NULL,
                          order = NULL, layout = NA,
                          plot_args = NULL, ...) {

  if (!is.list(ci) || !all(vapply(ci, is.matrix, logical(1L))))
    stop("`ci` must be a list of matrices (output of credint()).")

  # ---- Multi-element branch (correlations / covariances) -------------------
  if (length(ci) > 1L) {
    elem_names <- names(ci)
    if (is.null(elem_names) || any(elem_names == ""))
      stop("All elements of a multi-element `ci` must be named.")
    dots <- list(...)
    if (isTRUE(is.na(layout))) {
      sub_layout <- NA
    } else if (is.null(layout)) {
      sub_layout <- NULL
    } else {
      par(mfrow = layout)
      sub_layout <- NULL
    }
    for (nm in elem_names) {
      args <- modifyList(dots, if (!is.null(plot_args[[nm]]))
                                 plot_args[[nm]] else list())
      args$main  <- if (!is.null(plot_args[[nm]]$main)) plot_args[[nm]]$main else ""
      layout_val <- if (!is.null(plot_args[[nm]]$layout)) plot_args[[nm]]$layout else sub_layout
      args <- args[names(args) != "layout"]
      if (is.null(args$ylab)) args$ylab <- nm
      do.call(.plot_cor_panel, c(list(mat=ci[[nm]], quants=quants, layout=layout_val), args))
    }
    return(invisible(NULL))
  }
  # ---- Single-element branch -----------------------------------------------

  cl <- .classify_types(ci)

  # validate effect_only / intercept_only
  bad_eo <- effect_only[!effect_only %in% cl$all_types]
  bad_io <- intercept_only[!intercept_only %in% cl$all_types]
  if (length(bad_eo)) stop("Types in `effect_only` not found: ",
                            paste(bad_eo, collapse = ", "), ".")
  if (length(bad_io)) stop("Types in `intercept_only` not found: ",
                            paste(bad_io, collapse = ", "), ".")
  both <- intersect(effect_only, intercept_only)
  if (length(both)) stop("Types in both `effect_only` and `intercept_only`: ",
                          paste(both, collapse = ", "), ".")

  ordered_multi <- .resolve_order(order, cl$multi_types)
  dots <- list(...)

  # Auto-detect effect_only: multi-value types where ALL rownames contain "_"
  # (no bare intercept row, e.g. "sv" with only sv_WBblack, sv_NGgun, ...).
  # Types where no rowname has "_" are left as intercepts and aggregate into
  # the combined panel rather than getting individual effect-only panels.
  auto_eff <- cl$multi_types[vapply(cl$multi_types, function(t) {
    all(grepl("_", cl$rn[cl$types == t], fixed = TRUE))
  }, logical(1L))]
  new_eff   <- auto_eff[!auto_eff %in% c(effect_only, intercept_only)]
  effect_only <- c(effect_only, new_eff)

  # single-value types: split into aggregate vs own panels
  sv_agg <- cl$single_types[!cl$single_types %in% c(effect_only, intercept_only)]
  sv_eff <- cl$single_types[cl$single_types %in% effect_only]
  sv_int <- cl$single_types[cl$single_types %in% intercept_only]

  # build panel list
  panels <- list()
  if (length(sv_agg))
    panels <- c(panels, list(list(kind="agg",  types=sv_agg,
                                  main=names(ci)[1L])))
  for (t in sv_eff)
    panels <- c(panels, list(list(kind="eff",  t=t, main=t)))
  for (t in sv_int)
    panels <- c(panels, list(list(kind="int",  t=t, main=t)))
  for (t in ordered_multi) {
    kind <- if (t %in% effect_only) "eff" else
            if (t %in% intercept_only) "int" else "norm"
    panels <- c(panels, list(list(kind=kind, t=t, main=t)))
  }

  n_panels <- length(panels)
  if (n_panels == 0L) stop("No panels to plot.")

  # layout=NA (default): each type gets its own page (sub-calls manage their
  #   own mfrow via layout=NA, which for plot_credint = c(1,1)).
  # layout=NULL: don't touch mfrow at all.
  # layout=c(r,c): put all types in one global grid.
  if (isTRUE(is.na(layout))) {
    sub_layout <- NA     # each sub-call starts a new page
  } else if (is.null(layout)) {
    sub_layout <- NULL   # user controls mfrow
  } else {
    par(mfrow = layout)  # explicit global grid
    sub_layout <- NULL   # sub-calls don't touch mfrow
  }

  for (p in panels) {
    key  <- if (p$kind == "agg") "intercepts" else p$t
    base <- list(ci = ci, quants = quants)
    if (p$kind == "agg") {
      base$type <- p$types
    } else {
      base$type <- p$t
      if (p$kind == "eff") base["intercept"] <- list(NULL)
      if (p$kind == "int") base["effects"]   <- list(NULL)
    }
    .panel_call(plot_credint, base, key, plot_args, dots, p$main, sub_layout)
  }

  invisible(NULL)
}

#' Plot all credint(map=TRUE) types in a multi-panel display
#'
#' Automatically classifies parameter types from \code{credint(..., map=TRUE)}
#' output, aggregates single-value (no factorial structure) types into one
#' shared-scale panel via \code{\link{plot_credint}}'s vector-type mode, and
#' produces a \code{\link{plot_credint_map}} panel for each multi-value type.
#' Factorial structure (\code{factors}) for each multi-value type must be
#' supplied via \code{plot_args}.
#'
#' @param ci Length-1 list returned by \code{credint(..., map=TRUE)}.
#' @param quants Integer vector of length 3 (default \code{1:3}).
#' @param order \code{NULL} (show all multi-value types in row-name order) or
#'   a character/integer vector selecting and ordering panels.  Unlisted types
#'   are dropped.
#' @param layout As in \code{\link{plot_credints}}.
#' @param plot_args Named list of per-panel argument lists keyed by type name
#'   or \code{"intercepts"} for the aggregate panel.  Supply \code{factors}
#'   for each multi-value type here, e.g.
#'   \code{plot_args = list(v = list(factors = c("S", "E")))}.
#' @param ... Additional arguments forwarded to every sub-call.
#'
#' @return Invisibly returns \code{NULL}.
#' @seealso \code{\link{plot_credint_map}}, \code{\link{plot_credints}}
#' @export
plot_credints_map <- function(ci, quants = 1:3,
                              order = NULL, layout = NA,
                              plot_args = NULL, ...) {

  if (!is.list(ci) || length(ci) != 1L || !is.matrix(ci[[1L]]))
    stop("`ci` must be a length-1 list containing a matrix.")

  cl   <- .classify_types(ci)
  dots <- list(...)

  ordered_multi <- .resolve_order(order, cl$multi_types)

  # friendly check: multi-value types with no factors in plot_args
  missing_f <- ordered_multi[vapply(ordered_multi, function(t)
    is.null(plot_args[[t]]$factors), logical(1L))]
  if (length(missing_f)) {
    lines <- vapply(missing_f, function(t) {
      rows_t   <- cl$rn[cl$types == t]
      suffixes <- sub(paste0("^", t, "_?"), "", rows_t)
      paste0("  ", t, ": row suffixes are [", paste(suffixes, collapse=", "),
             "] -- supply factors via plot_args$", t, "$factors")
    }, character(1L))
    stop("Multi-value type(s) need `factors` in `plot_args`:\n",
         paste(lines, collapse="\n"))
  }

  # build panels.  Default titles for map=TRUE: the intercepts (aggregate) panel
  # and single-panel parameter types are left blank; multi-panel types fall back
  # to plot_credint_map's per-panel factor-level labels (main = NULL here).
  panels <- list()
  if (length(cl$single_types))
    panels <- c(panels, list(list(kind="agg",
                                  types=cl$single_types,
                                  main="")))
  for (t in ordered_multi)
    panels <- c(panels, list(list(kind="map", t=t, main=NULL)))

  n_panels <- length(panels)
  if (n_panels == 0L) stop("No panels to plot.")

  if (isTRUE(is.na(layout))) {
    sub_layout <- NA
  } else if (is.null(layout)) {
    sub_layout <- NULL
  } else {
    par(mfrow = layout)
    sub_layout <- NULL
  }

  for (p in panels) {
    key  <- if (p$kind == "agg") "intercepts" else p$t
    base <- list(ci = ci, quants = quants)
    if (p$kind == "agg") {
      base$type <- p$types
      .panel_call(plot_credint, base, key, plot_args, dots, p$main, sub_layout)
    } else {
      base$type <- p$t
      .panel_call(plot_credint_map, base, key, plot_args, dots, p$main, sub_layout)
    }
  }

  invisible(NULL)
}
