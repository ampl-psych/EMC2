prior_samples_alpha <- function(theta_mu,theta_var,n=1e3)
  # samples from prior for alpha implied by hyper model
  # (mixture over hyper posterior samples)
{
  if(is.null(dim(theta_mu))){
    out <- mvtnorm::rmvnorm(n,theta_mu,theta_var)
  } else{
    out <- t(sapply(sample.int(dim(theta_mu)[2],n,TRUE),function(x){
      mvtnorm::rmvnorm(1,theta_mu[x,],theta_var[,,x])
    }))
  }
  return(out)
}

get_prior_samples <- function(samples,selection,filter,thin,subfilter,n_prior)
  # get matrix of prior samples for different parameter types
{
  if (inherits(samples, "pmwgs")) samps <- samples else samps <- samples[[1]]
  if (selection=="alpha") {
    if(!is.null(samples[[1]]$samples$theta_mu)){
      if (inherits(samples, "pmwgs")) {
        theta_mu <- as_Mcmc(samples,selection="mu",filter=filter,
                            thin=thin,subfilter=subfilter)
        theta_var <- get_sigma(samples,filter=filter,
                               thin=thin,subfilter=subfilter)
      } else {
        theta_mu <- do.call(rbind,as_mcmc.list(samples,selection="mu",
                                               filter=filter,thin=thin,subfilter=subfilter))
        theta_var <-   abind(lapply(samples,get_sigma,filter=filter,thin=thin,subfilter=subfilter))
      }
      psamples <- prior_samples_alpha(theta_mu,theta_var,n_prior)
      colnames(psamples) <- colnames(theta_mu)
    } else{
      theta_mu <- samples[[1]]$prior$theta_mu_mean
      theta_var <- samples[[1]]$prior$theta_mu_var
      psamples <- prior_samples_alpha(theta_mu,theta_var,n_prior)
      colnames(psamples) <- samples[[1]]$par_names
    }
    return(psamples)
  } else {
    variant_funs <- attr(samples[[1]], "variant_funs")
    design <- list(model = attr(samples[[1]]$data[[1]], "model"))
    attr(design, "p_vector") <- samples[[1]]$samples$alpha[,1,1] # just need the names in the right format
    psamples <- variant_funs$get_prior(design = design, N = n_prior, prior = samples[[1]]$prior, type = selection, map = FALSE)[[1]]
    return(psamples)
  }
}


#' plot_prior
#' Plots prior distributions by simulation
#' @param prior A list of specifying prior means and covariance, if NULL default
#' prior (mean=0, variance = 1, uncorrelated)
#' @param design Design corresponding to prior
#' @param plotp Names of parameters to plot (default NULL plots all)
#' @param type Type of prior (standard or single)
#' @param selection Select level prior to plot (default "alpha" for single and "mu" for standard)
#' @param mapped boolean for mapping mu or alpha back to the design cells on the natural scale.
#' @param data data frame, required when mapping involves priors
#' @param N Number of prior samples if data not provided
#' @param nrep Number of prior samples as multiple of numebr of rows in the data (minus any that violate "ok" constraint)
#' @param breaks Histogram breaks parameter
#' @param layout par(mfrow) setting (default c(3,3))
#' @param lower Lower quantile limit of values plotted. Default NULL all plotted,
#' numeric same for all parameters, parameter named list parameter specific
#' @param upper Upper quantile limit of values plotted. Default NULL all plotted,
#' numeric same for all parameters, parameter named list parameter specific
#' @param xlim List with parameter names of plot x limits or single pair same for all.
#' Any names not in list or if (default) NA xlim set as min and max.
#'
#' @return Invisible sampled columns of parameters as a data frame (covariate,
#' version first column = cell names) or otherwise a matrix
#' @export
plot_prior <- function(prior=NULL, design,plotp=NULL,
                       type = c("standard","single")[1],selection = NULL,
                       mapped=TRUE,data=NULL,
                       N=1e5, nrep=10,
                       breaks=50,layout=c(3,3),lower=NULL,upper=NULL,xlim=NA)
{
  if (is.null(selection)) {
    if (type=="single"){
      selection <- "alpha"
    } else{
      selection <- "mu"
    }
  }

  if(selection == "alpha" & type !="single"){
    lower = .05
    upper = .95
  }

  if(selection == "variance"){
    upper = .95
  }
  if(selection == "covariance"){
    upper = .95
    lower = .05
  }
  if (mapped & !(selection %in% c("alpha","mu"))){
    warning("For selections other than mu and alpha mapped must be FALSE, we set it to FALSE here")
    mapped <- F
  }


  if (!(type %in% c("standard","single", "diagonal")))
    stop("Only types standard, diagonal and single implemented")
  if (type=="single" & selection !="alpha")
    stop("Can only select alpha for single")
  if (!is.null(data) & is.null(design))
    stop("Must design when data provided")
  if (mapped & is.null(design))
    stop("Must provide design when mapped=TRUE")
  if (type=="standard") gp <- get_prior_standard
  if (type=="single") gp <- get_prior_single
  if (type=="diagonal") gp <- get_prior_diag
  if (mapped & !is.null(data)) { # Used for covariates
    message("Mapping prior based on data, use this option with covariates and Ttranform parameters")
    if(selection == "alpha" & type != "single"){
      pp_mu <- gp(prior, type = "mu", design = design,N=dim(data)[1])$mu
      pp_var <- gp(prior, type = "full_var", design = design,N=dim(data)[1])$full_var
      n_pars <- ncol(pp_mu)
      pp <- matrix(0, N, ncol(pp_mu))
      for(i in 1:dim(data)[1]){
        if(type == "diagonal"){
          pp_var_curr <- diag(pp_var[i,])
        } else{
          pp_var_curr <- matrix(pp_var[i,], n_pars, n_pars)
        }
        pp[i,] <- rmvnorm(1, pp_mu[i,], pp_var_curr)      }
      colnames(pp) <- colnames(pp_mu)
    } else{
      pp <- gp(prior, type = selection, design = design,N=dim(data)[1])[[selection]]
    }

    row.names(pp) <- 1:dim(data)[1]
    design$Ffactors$subjects <- row.names(pp)
    data$subjects <- factor(1:dim(data)[1])
    mp <- make_data(pp,design,data=data,mapped_p=TRUE)
    mpok <- mp[design$model()$rfun(pars=mp),]
    if (nrep>1) for (i in 2:nrep) {
      if(selection == "alpha" & type != "single"){
        pp_mu <- gp(prior, type = "mu", design = design,N=dim(data)[1])$mu
        pp_var <- gp(prior, type = "full_var", design = design,N=dim(data)[1])$full_var
        n_pars <- ncol(pp_mu)
        pp <- matrix(0, N, ncol(pp_mu))
        for(i in 1:dim(data)[1]){
          if(type == "diagonal"){
            pp_var_curr <- diag(pp_var[i,])
          } else{
            pp_var_curr <- matrix(pp_var[i,], n_pars, n_pars)
          }
          pp[i,] <- rmvnorm(1, pp_mu[i,], pp_var_curr)
        }
        colnames(pp) <- colnames(pp_mu)
      } else{
        pp <- gp(prior, type = selection, design = design,N=dim(data)[1])[[selection]]
      }

      row.names(pp) <- 1:dim(data)[1]
      mp <- make_data(pp,design,data=data,mapped_p=TRUE)
      mpok <- rbind(mpok,mp[design$model()$rfun(pars=mp),])
    }
    facs <- names(design$Ffactors)
    fnam <- apply(mpok[,facs[facs!="subjects"]],2,as.character)
    for (i in 1:dim(fnam)[2]) fnam[,i] <- paste0(dimnames(fnam)[[2]][i],fnam[,i])
    fnam <- apply(fnam,1,paste,collapse="_")
    par(mfrow=layout)
    par_names <- unique(c(design$model()$p_types,plotp))
    par_names <- par_names[!(par_names %in% names(design$constants))]
    if (!is.null(plotp)) {
      if (!all(plotp %in% par_names)) stop("plotp not in prior")
      par_names <- par_names[par_names %in% plotp]
    }
    lowers <- setNames(rep(0,length(par_names)),par_names)
    uppers <- setNames(rep(1,length(par_names)),par_names)
    if (!is.null(upper)) if (is.null(names(upper)))
      uppers[1:length(uppers)] <- upper else
        uppers[names(upper)] <- upper
    if (!is.null(lower)) if (is.null(names(lower)))
      lowers[1:length(lowers)] <- lower else
        lowers[names(lower)] <- lower
    xlims <- setNames(vector(mode="list",length=length(par_names)),par_names)
    if (is.list(xlim)) for (i in names(xlim)) xlims[[i]] <- xlim[[i]] else
      if (!all(is.na(xlim))) for (i in names(xlim)) xlims[[i]] <- xlim
    for (pnam in par_names) {
      lower <- quantile(mpok[,pnam], probs = lowers[[pnam]])
      upper <- quantile(mpok[,pnam], probs = uppers[[pnam]])
      if (is.null(xlims[[pnam]])) {
        tmp <- mpok[,pnam][(mpok[,pnam] >= lower) & (mpok[,pnam] <= upper)]
        xlims[[pnam]] <- c(min(tmp),max(tmp))
      }
      for (i in sort(unique(fnam))) {
        filtered <- mpok[fnam==i,pnam]
        filtered <- filtered[(filtered >= lower) & (filtered <= upper)]
        hist(filtered,prob = TRUE,
             main=i,xlab=pnam,xlim=xlims[[pnam]],breaks=breaks,
             cex.lab = 1.25, cex.main = 1.5, xlab = selection)
      }
    }
    invisible(cbind.data.frame(fnam,mpok[,par_names,drop=FALSE]))
  } else {

    if(selection == "alpha" & type != "single"){
      pp_mu <- gp(prior, type = "mu", design = design,N=N)$mu
      pp_var <- gp(prior, type = "full_var", design = design,N=N)$full_var
      n_pars <- ncol(pp_mu)
      samples <- matrix(0, N, ncol(pp_mu))
      for(i in 1:N){
        if(type == "diagonal"){
          pp_var_curr <- diag(pp_var[i,])
        } else{
          pp_var_curr <- matrix(pp_var[i,], n_pars, n_pars)
        }
        samples[i,] <- rmvnorm(1, pp_mu[i,], pp_var_curr)
      }
      colnames(samples) <- colnames(pp_mu)
    } else{
      samples <- gp(prior, type = selection, design = design,N=N)[[selection]]
    }

    if(selection %in% c("covariance", "correlation")){
      pnams <- names(attr(design, "p_vector"))
      lt <- lower.tri(diag(length(pnams)))
      pnams <- outer(pnams,pnams,paste,sep=".")[lt]
      colnames(samples) <- pnams
    }

    if (mapped)
      samples <- map_mcmc(samples,design,design$model,include_constants=FALSE)
    par(mfrow = layout)
    par_names <- colnames(samples)
    if (!is.null(plotp)) {
      if (!all(plotp %in% par_names)) stop("plotp not in prior")
      par_names <- par_names[par_names %in% plotp]
    }
    lowers <- setNames(rep(0,length(par_names)),par_names)
    uppers <- setNames(rep(1,length(par_names)),par_names)
    if (!is.null(upper)) if (is.null(names(upper)))
      uppers[1:length(uppers)] <- upper else
        uppers[names(upper)] <- upper
    if (!is.null(lower)) if (is.null(names(lower)))
      lowers[1:length(lowers)] <- lower else
        lowers[names(lower)] <- lower
    xlims <- setNames(vector(mode="list",length=length(par_names)),par_names)
    if (!is.null(xlim)) xlims[names(xlim)] <- xlim
    for(i in 1:ncol(samples)){
      lower <- quantile(quantile(samples[,i], probs = lowers[i]))
      upper <- quantile(quantile(samples[,i], probs = uppers[i]))
      filtered <- samples[,i][(samples[,i] >= lower) & (samples[,i] <= upper)]
      if (is.null(xlims[[i]]))
        hist(filtered, breaks = breaks, main = par_names[i], prob = TRUE,
             xlab = selection, cex.lab = 1.25, cex.main = 1.5) else
               hist(filtered, breaks = breaks, prob = TRUE,
                    xlab = selection, cex.lab = 1.25, cex.main = 1.5,xlim=xlims[[i]])
    }
    invisible(samples)
  }
}


#' make_prior
#'
#' Makes priors based on a design.
#'  For hierarchical models you specify prior on mean of the group mean (pmean) and standard deviations
#' (psd), and possibly these values taken from another prior (supplied by the
#' update argument), as long as those values have the same name as in the current
#' design and are not specified in pmean and psd. Where a value is not supplied
#' in arguments the user is prompted to enter numeric values (or functions that
#' evaluate to numbers).
#'
#' @param design Design for which a prior is constructed
#' @param pmean Named vector of prior means, or an unnamed scalar, which is then used for all.
#' @param psd Named vector of prior standard deviations, or an unnamed scalar, which is then used for all.
#' @param update Another prior from which to copy meand and sds.
#' @param verbose Boolean (default true) print values of prior to console (if
#' update only new values).
#' @param update_print When verbose print only new values (default TRUE)
#' @param type What type of model you plan on using, choice of standard, diagonal and single for this function
#' @param pscale For hierarchical models, the prior on the scale of the variances
#' @param df For hierarchical models, the prior on the degrees of freedom of the variances
#'
#' @return An EMC prior object
#' @export
make_prior <- function(design,pmean=NULL,psd=NULL,update=NULL,
                       verbose=TRUE,update_print=TRUE, type = "standard",
                       pscale = NULL, df = NULL)
{
  pm <- sampled_p_vector(design)
  pm <- ps <- psc <-  pm[1:length(pm)]
  if (!is.null(pmean) && is.null(names(pmean))) {
    pmean[1:length(pm)] <- pm
    pmean <- setNames(pmean,names(pm))
  }
  if (!is.null(psd) && is.null(names(psd))) {
    psd[1:length(ps)] <- psd
    psd <- setNames(psd,names(ps))
  }
  if (!is.null(pscale) && is.null(names(pscale))) {
    pscale[1:length(psc)] <- pscale
    pscale <- setNames(pscale,names(psc))
  }

  if ( !is.null(update) ) {
    pmu <- update$theta_mu_mean
    if(!type == "single"){
      pdfu <- update$v
      pscaleu <- update$A
      names(pscaleu) <- names(pmu)
      if(!is.null(pdfu) & is.null(df)){
        df <- pdfu
      }
    }
    psdu <- setNames(diag(update$theta_mu_var)^.5,names(pmu))
    isin <- names(pmu)[names(pmu) %in% names(pm)]
    if (length(isin)>0) {
      addn <- isin[!(isin %in% names(pmean))]
      pmean <- c(pmean,pmu[addn])
      addn <- isin[!(isin %in% names(psd))]
      psd <- c(psd,psdu[addn])
      if(!type == "single"){
        addn <- isin[!(isin %in% names(pscaleu))]
        psc <- c(psc,pscaleu[addn])
      }
    }
  }
  if (!all(names(pmean) %in% names(pm))) stop("pmean has names not in design")
  if (!is.null(pmean)) pm[names(pmean)] <- pmean
  todom <- !(names(pm) %in% names(pmean))
  if (any(todom)) {
    if(type == "single"){
      cat("Enter values for prior mean\n")
    } else{
      cat("Enter values for prior mean of group-level mean\n")
    }
    for (i in names(pm[todom])) {
      repeat {
        ans <- try(eval(parse(text=readline(paste0(i,": ")))),silent=TRUE)
        if (any(class(ans) %in% c("warning", "error", "try-error")) || is.na(ans) ) {
          cat("Must provide a numeric value\n")
        } else {
          pm[i] <- ans
          break
        }
      }
    }
  }
  if (!all(names(psd) %in% names(ps))) stop("psd has names not in design")
  if (!is.null(ps)) ps[names(psd)] <- psd
  todos <- !(names(ps) %in% names(psd))
  if (any(todos)) {
    if(type == "single"){
      cat("Enter values for prior standard deviation\n")
    } else{
      cat("Enter values for prior standard deviation of group-level mean\n")
    }
    for (i in names(pm[todos])) {
      repeat {
        ans <- try(eval(parse(text=readline(paste0(i,": ")))),silent=TRUE)
        if (any(class(ans) %in% c("warning", "error", "try-error")) || is.na(ans) ) {
          cat("Must provide a numeric value\n")
        } else {
          ps[i] <- ans
          break
        }
      }
    }
  } else if ( is.null(update) ) ps <- rep(1,length(ps))
  if(type != "single"){
    todoscale <- !(names(psc) %in% names(pscale))
    if (any(todoscale)) {
      cat("Enter values for prior scale of group-level variance, larger values leads to broader variance, default is 1 \n")
      for (i in names(pm[todos])) {
        repeat {
          ans <- try(eval(parse(text=readline(paste0(i,": ")))),silent=TRUE)
          if (any(class(ans) %in% c("warning", "error", "try-error")) || is.na(ans) ) {
            cat("Must provide a numeric value\n")
          } else {
            psc[i] <- ans
            break
          }
        }
      }
    } else if ( is.null(update) ) psc <- rep(1,length(psc))
    if(is.null(df)){
      if(type == "standard"){
        cat("Enter value for prior degrees of freedom for group-level variance, same for all parameters, 2 leads to uniform priors on correlations \n")

      } else{
        cat("Enter value for prior degrees of freedom for group-level variance, same for all parameters")
      }
      repeat {
        ans <- try(eval(parse(text=readline("df : "))),silent=TRUE)
        if (any(class(ans) %in% c("warning", "error", "try-error")) || is.na(ans) ) {
          cat("Must provide a numeric value\n")
        } else {
          df <- ans
          break
        }
      }
    }
  }


  if (verbose & (!update_print | (update_print & (!any(todom) | !any(todos))))) {
    if (!(update_print & !any(todom))) cat("Prior means\n")
    if (!update_print) print(pm) else if (any(todom)) print(pm[todom])
    if (!(update_print & !any(todos))) cat("\nPrior standard deviations\n")
    if (!update_print) print(ps) else if (any(todos)) print(ps[todos])
    if(type != "single"){
      if (!(update_print & !any(todos))) cat("\nPrior scale on the variances\n")
      if (!update_print) print(ps) else if (any(todoscale)) print(pscale[todoscale])
    }
  }
  if(type != "single"){
    return(  list(theta_mu_mean  = pm,theta_mu_var = diag(ps^2),
                  A = psc, v = df))
  } else{
    list(theta_mu_mean  = pm,theta_mu_var = diag(ps^2))
  }

}



