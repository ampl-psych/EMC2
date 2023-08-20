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
#' prior (mean=0, varaince = 1, uncorrelaterd)
#' @param design Design corresponding to prior
#' @param plotp Names of parameters to plot (default NULL plots all)
#' @param type Type of prior (standard or single)
#' @param selection Select level prior to plot (default "alpha" for single and "mu" for single)
#' @param mapped boolean for mapping mu or alpha back to the design cells on the natural scale.
#' @param data data frame, required when mapping involves priors
#' @param N Number of prior samples if data not provided
#' @param nrep Number of prior samples as multiple of numebr of rows in the data (minus any that violate "ok" constraint)
#' @param breaks Histogram breaks parameter
#' @param layout par(mfrow) setting (default c(3,3))
#' @param lower Lower quantile limit of values plotted. Default NULL all plotted,
#' integer same for all parameters, parameter named list parameter specific
#' @param lower Upper quantile limit of values plotted. Default NULL all plotted,
#' integer same for all parameters, parameter named list parameter specific
#' @param xlim List with parameter names of plot x limits or single pair same for all.
#' Any names not in list or if (defualt) NA xlim set as min and max.
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
    if (type=="standard") selection <- "mu"
    if (type=="single") selection <- "alpha"
  }
  if (mapped & !(selection %in% c("alpha","mu")))
    stop("For selections other than mu and alpha set mapped = FALSE")
  if (!(type %in% c("standard","single")))
    stop("Only types standard and single implemented")
  if (type=="single" & selection !="alpha")
    stop("Can only select alpha for single")
  if (!is.null(data) & is.null(design))
    stop("Must design when data provided")
  if (mapped & is.null(design))
    stop("Must provide design when mapped=TRUE")
  if (type=="standard") gp <- get_prior_standard
  if (type=="single") gp <- get_prior_single
  if (mapped & !is.null(data)) { # Used for covariates
    message("Mapping prior based on data, use this option with covariates and Ttranform parameters")
    pp <- gp(prior, design = design,N=dim(data)[1])[[selection]]
    row.names(pp) <- 1:dim(data)[1]
    design$Ffactors$subjects <- row.names(pp)
    data$subjects <- factor(1:dim(data)[1])
    mp <- make_data(pp,design,data=data,mapped_p=TRUE)
    mpok <- mp[design$model()$rfun(pars=mp),]
    if (nrep>1) for (i in 2:nrep) {
      pp <- gp(prior, design = design,N=dim(data)[1])[[selection]]
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
    samples <- gp(prior,design=design,N=N)[[selection]]
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
             xlab = type, cex.lab = 1.25, cex.main = 1.5) else
        hist(filtered, breaks = breaks, prob = TRUE,
             xlab = par_names[i], cex.lab = 1.25, cex.main = 1.5,xlim=xlims[[i]])
    }
    invisible(samples)
  }
}

# plot_prior <- function(prior, type = NULL,add_density=FALSE,adjust=1,breaks=50,
#                        layout=c(3,3),upper=NULL,xlim=NULL,mapped=TRUE,design=NULL){
#
#
#   if(is.null(type)) type <- names(prior)
#   for(typ in type) {
#     samples <- prior[[typ]]
#     if (mapped & (typ %in% c("alpha","mu"))) {
#       if (is.null(design)) stop("Must provide design when mapped=TRUE")
#       samples <- map_mcmc(samples,design,design$model,include_constants=FALSE)
#     }
#     par(mfrow = layout)
#     par_names <- colnames(samples)
#     uppers <- setNames(rep(.999,length(par_names)),par_names)
#     if (!is.null(upper)) if (is.null(names(upper)))
#       uppers[1:length(uppers)] <- upper else
#       uppers[names(upper)] <- upper
#     xlims <- setNames(vector(mode="list",length=length(par_names)),par_names)
#     if (!is.null(xlim)) xlims[names(xlim)] <- xlim
#     for(i in 1:ncol(samples)){
#       if(!any(samples[,i] < 0) || !any(samples[,i] > 0)){
#         quants <- quantile(abs(samples[,i]), probs = uppers[i])
#       } else{
#         quants <- quantile(abs(samples[,i]), probs = uppers[i])
#       }
#       filtered <- samples[,i][abs(samples[,i]) < quants]
#       if (is.null(xlims[[i]]))
#           hist(filtered, breaks = breaks, main = par_names[i], prob = TRUE,
#            xlab = type, cex.lab = 1.25, cex.main = 1.5) else
#           hist(filtered, breaks = breaks, main = par_names[i], prob = TRUE,
#            xlab = type, cex.lab = 1.25, cex.main = 1.5,xlim=xlims[[i]])
#       if (add_density) lines(density(filtered,adjust=adjust), col = "red")
#     }
#   }
# }
#
# prior_samples <- function(samps,type=c("mu","variance","covariance","correlation","sigma")[1],n=1e4)
#   # Samples from prior for standard at hyper levels for three types,
#   # or for single at subject level (in which case type="mu")
# {
#   if (type=="mu") {
#     out <- mvtnorm::rmvnorm(n, mean = samps$prior$theta_mu_mean,
#                             sigma = samps$prior$theta_mu_var)
#     dimnames(out)[[2]] <- samps$par_names
#     return(out)
#   } else {
#     var <- array(NA_real_, dim = c(samps$n_pars, samps$n_pars, n),
#                  dimnames = list(samps$par_names, samps$par_names, NULL))
#     # Can't specify n in riwish?
#     hyper <- attributes(samps)
#     for(i in 1:n){
#       a_half <- 1 / rgamma(n = samps$n_pars,shape = 1/2,
#                            rate = 1/(hyper$A_half^2))
#       var[,,i] <- riwish(hyper$v_half + samps$n_pars - 1, 2 * hyper$v_half * diag(1 / a_half))
#     }
#     if (type=="sigma") return(var)
#     if (type=="variance") return(t(apply(var,3,diag)))
#     if (type=="correlation")
#       var <- array(apply(var,3,cov2cor),dim=dim(var),dimnames=dimnames(var))
#     lt <- lower.tri(var[,,1])
#     return(t(apply(var,3,function(x){x[lt]})))
#   }
# }




