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


#' Prior plotting
#'
#' Plots the prior for a model by simulating from the user specified prior distributions.
#'
#' To get the default prior for a type, run: `get_prior_{type}(design = design, sample = F)`.
#' E.g.: `get_prior_diagonal(design = design, sample = F)`
#'
#' @param prior A list specifying the priors for the chosen type.
#' @param design Model design corresponding to the prior
#' @param plotp Names of parameters to plot (the default `NULL` plots all parameters)
#' @param type Type of prior (`standard`, `diagonal`, `single`, `blocked` or `factor`).
#' For `blocked`, the blocked parameters are still plotted in the covariance matrix.
#' @param selection Selects which parameter groups to plot (default `"alpha"`
#' for single and `"mu"` for standard)
#' @param mapped Boolean for mapping `mu` or `alpha` back to the design cells on the natural scale.
#' @param data Data frame, required when the mapping involves priors
#' @param N Number of prior samples if data not provided
#' @param nrep Number of prior samples as multiple of the number of rows in the data frame.
#' @param breaks breaks parameter for the histogram
#' @param layout `par(mfrow)` setting (the default is `c(3,3)`)
#' @param lower Defaults to `NULL`. Lower quantile limit of values plotted. If numeric, the same is used for all parameters,
#' alternatively, a (parameter) named list with individual quantile limits can be supplied.
#' @param upper Defaults to `NULL`. Upper quantile limit of values plotted. If numeric, the same is used for all parameters,
#' alternatively, a (parameter) named list with individual quantile limits can be supplied.
#' @param xlim A list with parameter names of x-axis limits or a single pair if the
#' same limits should be used for all parameters. Any names not in list or if (default) `NA`, `xlim` the minimum and maximum are set.
#' @param ... For additional arguments
#'
#' @return Invisible of the prior samples
#' @examples \dontrun{
#' # First define a design for the model
#' design_DDMaE <- make_design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then set up a prior using make_prior
#' p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#' t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' psd <- c(v_Sleft=1,v_Sright=1,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
#' t0=.4,Z=1,sv=.4,SZ=1)
#' # Here we left the variance prior at default
#' prior_DDMaE <- make_prior(design_DDMaE,pmean=p_vector,psd=psd)
#' # Now we can plot all sorts of (implied) priors
#' plot_prior(prior_DDMaE, design_DDMaE, selection = "mu")
#' plot_prior(prior_DDMaE, design_DDMaE, selection = "mu", mapped = FALSE)
#' plot_prior(prior_DDMaE, design_DDMaE, selection = "variance")
#' plot_prior(prior_DDMaE, design_DDMaE, selection = "covariance")
#' plot_prior(prior_DDMaE, design_DDMaE, selection = "correlation")
#' # We can also plot the implied prior on the participant level effects.
#' plot_prior(prior_DDMaE, design_DDMaE, selection = "alpha")
#' }
#' @export
plot_prior <- function(prior=NULL, design,plotp=NULL,
                       type = "standard",selection = NULL,
                       mapped=TRUE,data=NULL,
                       N=1e5, nrep=10,
                       breaks=50,layout=c(3,3),lower=NULL,upper=NULL,xlim=NA,
                       ...)
{
  if (is.null(selection)) {
    if (type=="single"){
      selection <- "alpha"
    } else{
      selection <- "mu"
    }
  }

  if(selection == "alpha" & type !="single"){
    if(is.null(lower))lower = .05
    if(is.null(upper)) upper = .95
  }

  if(selection == "variance"){
    if(is.null(upper)) upper = .95
  }
  if(selection == "covariance" || selection == "loadings"){
    if(is.null(lower))lower = .02
    if(is.null(upper)) upper = .98
  }
  if (mapped & !(selection %in% c("alpha","mu"))){
    # warning("For selections other than mu and alpha mapped must be FALSE, we set it to FALSE here")
    mapped <- FALSE
  }


  if (!(type %in% c("standard","single", "diagonal", "factor", "infnt_factor")))
    stop("Only types standard, diagonal, factor, infnt_factor and single implemented")
  if (type=="single" & selection !="alpha")
    stop("Can only select alpha for single")
  if (!is.null(data) & is.null(design))
    stop("Must design when data provided")
  if (mapped & is.null(design))
    stop("Must provide design when mapped=TRUE")
  if (type=="standard") gp <- get_prior_standard
  if (type=="single") gp <- get_prior_single
  if (type=="diagonal") gp <- get_prior_diag
  if (type=="infnt_factor") gp <- get_prior_infnt_factor
  if (type=="factor") gp <- get_prior_factor
  if (mapped & !is.null(data)) { # Used for covariates
    message("Mapping prior based on data, use this option with covariates and Ttranform parameters")
    if(selection == "alpha" & type != "single"){
      pp_mu <- gp(prior, type = "mu", design = design,N=dim(data)[1], ...)$mu
      pp_var <- gp(prior, type = "full_var", design = design,N=dim(data)[1], ...)$full_var
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
      pp <- gp(prior, type = selection, design = design,N=dim(data)[1], ...)[[selection]]
    }

    row.names(pp) <- 1:dim(data)[1]
    design$Ffactors$subjects <- row.names(pp)
    data$subjects <- factor(1:dim(data)[1])
    mp <- make_data(pp,design,data=data,mapped_p=TRUE)
    mpok <- mp[design$model()$rfun(pars=mp),]
    if (nrep>1) for (i in 2:nrep) {
      if(selection == "alpha" & type != "single"){
        pp_mu <- gp(prior, type = "mu", design = design,N=dim(data)[1], ...)$mu
        pp_var <- gp(prior, type = "full_var", design = design,N=dim(data)[1], ...)$full_var
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
        pp <- gp(prior, type = selection, design = design,N=dim(data)[1], ...)[[selection]]
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
             cex.lab = 1.25, cex.main = 1.5)
      }
    }
    invisible(cbind.data.frame(fnam,mpok[,par_names,drop=FALSE]))
  } else {

    if(selection == "alpha" & type != "single"){
      pp_mu <- gp(prior, type = "mu", design = design,N=N, ...)$mu
      pp_var <- gp(prior, type = "full_var", design = design,N=N, ...)$full_var
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
      samples <- gp(prior, type = selection, design = design,N=N, ...)[[selection]]
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

check_prior <- function(prior, sampled_p_names){
  if (is.null(names(prior$theta_mu_mean)))
    names(prior$theta_mu_mean) <- sampled_p_names
  if(length(prior$theta_mu_mean) != length(sampled_p_names))
    stop("prior mu should be same length as estimated parameters (p_vector)")
  if (!all(sort(names(prior$theta_mu_mean)) == sort(sampled_p_names)))
    stop("theta_mu_mean names not the same as sampled paramter names")
  # Make sure theta_mu_mean has same order as sampled parameters
  pnams <- names(prior$theta_mu_mean)
  prior$theta_mu_mean <- prior$theta_mu_mean[sampled_p_names]
  if(!is.matrix(prior$theta_mu_var)) {
    if(length(prior$theta_mu_var) != length(sampled_p_names))
      stop("prior theta should be same length as estimated parameters (p_vector)")
    # Make sure theta_mu_var has same order as sampled parameters
    names(prior$theta_mu_var) <- pnams
    prior$theta_mu_var <- prior$theta_mu_var[sampled_p_names]
  } else {
    if(nrow(prior$theta_mu_var) != length(sampled_p_names))
      stop("prior theta should have same number of rows as estimated parameters (p_vector)")
    if(ncol(prior$theta_mu_var) != length(sampled_p_names))
      stop("prior theta should have same number of columns as estimated parameters (p_vector)")
    # Make sure theta_mu_var has same order as sampled parameters
    dimnames(prior$theta_mu_var) <- list(pnams,pnams)
    prior$theta_mu_var <- prior$theta_mu_var[sampled_p_names,sampled_p_names]
  }
  return(prior)
}

merge_priors <- function(prior_list){
  out <- prior_list[[1]]
  if(length(prior_list) > 1){
    out_names <- names(out)
    for(j in 2:length(prior_list)){
      for(hyper in out_names){
        if(length(out[[hyper]]) > 1){
          if(length(dim(out[[hyper]])) == 2){
            out[[hyper]] <- adiag(out[[hyper]], prior_list[[j]][[hyper]])
          } else{
            out[[hyper]] <- c(out[[hyper]], prior_list[[j]][[hyper]])
          }
        }
      }
    }
  }
  return(out)
}


#' Prior specification
#'
#' Specify priors on mean of the group mean (`pmean`) and standard deviation (`psd`),
#' for the individual level for a `single` model and for the group level for
#' hierarchical models. These values are entered manually by default but can be
#' recycled from another prior (given in the `update` argument).
#'
#' Where a value is not supplied, the user is prompted to enter
#' numeric values (or functions that evaluate to numbers). For hierarchical models
#' the scale (`pscale`) and the  degrees of freedom (`df`) for the population variance
#' covariance matrix can also be specified. Set these arguments to `NULL` to
#'  enter values manually through prompts.
#'
#' To get the default prior for a type, run: `get_prior_{type}(design = design, sample = F)`
#' E.g.: `get_prior_diagonal(design = design, sample = F)`
#'
#' @param design Design for which a prior is constructed, typically the output of `make_design()`
#' @param pmean Named vector of prior means, or an unnamed scalar, which is then used for all.
#' @param psd Named vector of prior standard deviations, or an unnamed scalar,
#' which is then used for all. If `NA`, will assume 1 for all parameters
#' @param update Another prior from which to copy values
#' @param verbose Boolean, defaults to `TRUE`, print prior specification to console (if
#' `update`, only new values).
#' @param update_print When verbose, print only new values (defaults to `TRUE`)
#' @param type What type of model you plan on using, choice of `standard`, `diagonal`, `blocked`
#' and `single` for this function
#' @param pscale For hierarchical models, the prior on the scale of the variances.
#' The default `NA` will assume 1 for all parameters
#' @param df For hierarchical models, the prior on the degrees of freedom of the
#' variances. The default `NA` will assume 2
#'
#' @return An EMC prior list object
#' @examples \dontrun{
#' # First define a design for the model
#' design_DDMaE <- make_design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then set up a prior using make_prior
#' p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#' t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' psd <- c(v_Sleft=1,v_Sright=1,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
#' t0=.4,Z=1,sv=.4,SZ=1)
#' # Here we left the variance prior at default
#' prior_DDMaE <- make_prior(design_DDMaE,pmean=p_vector,psd=psd)
#' # Also add a group-level variance prior:
#' pscale <- c(v_Sleft=.6,v_Sright=.6,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
#' t0=.2,Z=.5,sv=.4,SZ=.3)
#' psd <- .4
#' prior_DDMaE <- make_prior(design_DDMaE,pmean=p_vector,psd=psd, pscale = pscale, psd = psd)
#' # If we specify a new design
#' design_DDMat0E <- make_design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~E, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # We can easily update the prior
#' prior_DDMat0E <- make_prior(design_DDMat0E, update = prior_DDMaE)
#' }
#' @export
make_prior <- function(design,pmean=NULL,psd=NULL,update=NULL,
                       verbose=TRUE,update_print=TRUE, type = "standard",
                       pscale = NA, df = NA)
{
  pm <- sampled_p_vector(design)
  pm <- ps <- psc <-  pm[1:length(pm)]

  # Make sure names are correct and fill in defaults if NA
  if (!is.null(pmean) && is.null(names(pmean))) {
    pmean[1:length(pm)] <- pm
    pmean <- setNames(pmean,names(pm))
  }
  if (!is.null(psd) && is.null(names(psd))) {
    if(any(is.na(psd))){
      psd[1:length(ps)] <- 1
      psd <- setNames(psd,names(ps))
    } else{
      psd[1:length(ps)] <- psd
      psd <- setNames(psd,names(ps))
    }
  }
  if (!is.null(pscale) && is.null(names(pscale))) {
    if(any(is.na(pscale))){
      pscale[1:length(psc)] <- 0.3
      pscale <- setNames(pscale,names(psc))
    } else{
      pscale[1:length(psc)] <- pscale
      pscale <- setNames(pscale,names(psc))
    }
  }
  if (is.na(df)) df <- 2

  if ( !is.null(update) ) {
    if(is(update[[1]], "pmwgs")) update <- update[[1]]$prior
    if (is.null(update$v)) type <- "single"
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
        addn <- isin[!(isin %in% names(pscale))]
        pscale <- c(pscale,pscaleu[addn])
      }
    }
  }
  if (!all(names(pmean) %in% names(pm))) stop("pmean has names not in design")
  if (!is.null(pmean)) pm[names(pmean)] <- pmean
  todom <- !(names(pm) %in% names(pmean))
  if ( any(todom) ) {
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
  if ( any(todos) ) {
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
  }
  if(type != "single"){
    if (!is.null(psc)) psc[names(pscale)] <- pscale
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
    }
    if ( is.null(df) ) {
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

make_prior_new <- function(design, type = "standard", update = NULL, ...){
  prior <- get_objects(design = design, type = type)
  args <- list(...)
  for(pri in names(prior$prior)){ # Loop over different objects in prior
    if(pri %in% names(prior$descriptions)){ # This excluded prior$theta_mu_invar
      if(pri %in% names(args)){ # Already specified in ellipsis, so fill in
        input <- args[[pri]]
        if(length(input) == nrow(prior$prior[[pri]])){
          input <- diag(input)
        }
        if(!is.null(dim(prior$prior[[pri]]))){
          if(!identical(dim(input),dim(prior$prior[[pri]]))){
            stop("make sure your prior matches the design and the type of model. See ?get_prior_{type}")
          }
        } else{
          if(length(input) != length(prior$prior[[pri]])){
            stop("make sure your prior matches the design and the type of model. See ?get_prior_{type}")
          }
        }
        prior$prior[[pri]] <- input
      } else{ # Ask the user to manually specify
        cat(paste0("Specify ", prior$descriptions[[pri]]," \n"))
        if(!is.null(dim(prior$prior[[pri]]))){
          tmp <- diag(prior$prior[[pri]])
        } else{
          tmp <- prior$prior[[pri]]
        }
        cat("Press enter to fill remaining with default value (", tmp[1], ")")
        filled <- F
        for(i in 1:length(tmp)){ # Loop over the length of the prior object
          if(!filled){
            name <- names(tmp)[i]
            if(is.null(name)){
              name <- pri
            }
            repeat {
              ans <- try(eval(parse(text=readline(paste0(name,": ")))),silent=TRUE)
              if(!is.null(ans)){
                if (any(class(ans) %in% c("warning", "error", "try-error")) || is.na(ans) || !is.numeric(ans)) {
                  cat("Must provide a numeric value\n")
                } else {
                  tmp[i] <- ans
                  break
                }
              } else{
                filled <- TRUE
                break
              }
            }
          }

        }
        if(!is.null(dim(prior$prior[[pri]]))){
          prior$prior[[pri]][,] <- diag(tmp)

        } else{
          prior$prior[[pri]] <- tmp
        }
      }
    }
  }
  # Fix that invar is the same as var..
  return(prior$prior)
}

