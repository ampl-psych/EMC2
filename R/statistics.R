std_error_IS2 <- function(IS_samples, n_bootstrap = 50000){
  log_marglik_boot= array(dim = n_bootstrap)
  for (i in 1:n_bootstrap){
    log_weight_boot = sample(IS_samples, length(IS_samples), replace = TRUE) #resample with replacement from the lw
    log_marglik_boot[i] <- median(log_weight_boot)
  }
  return(sd(log_marglik_boot))
}

dhalft <- function (x, scale = 25, nu = 1, log = FALSE)
{
  x <- as.vector(x)
  scale <- as.vector(scale)
  nu <- as.vector(nu)
  if (any(scale <= 0))
    stop("The scale parameter must be positive.")
  NN <- max(length(x), length(scale), length(nu))
  x <- rep(x, len = NN)
  scale <- rep(scale, len = NN)
  nu <- rep(nu, len = NN)
  dens <- log(2) - log(scale) + lgamma((nu + 1)/2) - lgamma(nu/2) -
    0.5 * log(pi * nu) - (nu + 1)/2 * log(1 + (1/nu) * (x/scale) *
                                            (x/scale))
  if (log == FALSE)
    dens <- exp(dens)
  return(dens)
}

logdinvGamma <- function(x, shape, rate){
  alpha <- shape
  beta <- 1/rate
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha +
                                                        1) * log(x) - (beta/x)
  return(pmax(log.density, -500)) #Roughly equal to 1e-22 on real scale
}

es_pmwg <- function(pmwg_mcmc,selection="alpha",summary_alpha=mean,
                    print_summary=TRUE,sort_print=TRUE,
                    filter="sample",thin=1,subfilter=NULL)
  # Effective size
{
  if (!(class(pmwg_mcmc[[1]]) %in% c("mcmc","mcmc.list"))) {
    if (is(pmwg_mcmc, "pmwgs")){
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter)
    }
    else{
      pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                thin=thin,subfilter=subfilter)
    }
  }
  if (attr(pmwg_mcmc,"selection")=="LL")
    stop("Effective size not sensible for LL\n")
  out <- do.call(rbind,lapply(pmwg_mcmc,effectiveSize))
  if (attr(pmwg_mcmc,"selection")=="alpha") {
    if (!is.null(summary_alpha)) out <- apply(out,2,summary_alpha)
    if (print_summary) if (sort_print) print(round(sort(out))) else
      print(round(out))
    invisible(out)
  } else {
    out <- apply(out,2,sum)
    if (print_summary) if (sort_print) print(round(sort(out))) else
      print(round(out))
    invisible(out)
  }
}


gd_pmwg <- function(pmwg_mcmc,return_summary=FALSE,print_summary=TRUE,
                    digits_print=2,sort_print=TRUE,autoburnin=FALSE,transform=TRUE, omit_mpsrf = TRUE,
                    selection="alpha",filter="sample",thin=1,subfilter=NULL,mapped=FALSE)
  # R hat, prints multivariate summary returns each participant unless +
  # multivariate as matrix unless !return_summary
{

  split_mcl <- function(mcl)
    # Doubles chains by splitting into first and secon half
  {
    if (!is.list(mcl)) mcl <- list(mcl)
    mcl2 <- mcl
    half <- floor(unlist(lapply(mcl,nrow))/2)
    for (i in 1:length(half)) {
      mcl2[[i]] <- coda::as.mcmc(mcl2[[i]][c((half[i]+1):(2*half[i])),])
      mcl[[i]] <- coda::as.mcmc(mcl[[i]][1:half[i],])
    }
    coda::as.mcmc.list(c(mcl,mcl2))
  }

  gelman_diag_robust <- function(mcl,autoburnin,transform, omit_mpsrf)
  {
    mcl <- split_mcl(mcl)
    gd <- try(gelman.diag(mcl,autoburnin=autoburnin,transform=transform, multivariate = !omit_mpsrf),silent=TRUE)
    if (is(gd, "try-error")){
      if(omit_mpsrf){
        return(list(psrf=matrix(Inf)))
      } else{
        return(list(psrf=matrix(Inf),mpsrf=Inf))
      }
    } else{
      return(gd)
    }
  }

  if ( selection=="LL" ) stop("Rhat not appropriate for LL")
  if (class(pmwg_mcmc[[1]]) %in% c("mcmc","mcmc.list")) {
    if (mapped) warning("Cannot transform to natural scale unless samples list provided")
  } else {
    if (is(pmwg_mcmc, "pmwgs")) {
      if (mapped) warning("Cannot transform to natural scale unless samples list provided")
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter)
    } else
      pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                thin=thin,subfilter=subfilter,mapped=mapped)
  }
  if (selection=="alpha" || selection == "random") {
    gd <- lapply(pmwg_mcmc,gelman_diag_robust,autoburnin = autoburnin, transform = transform, omit_mpsrf = omit_mpsrf)
    if(omit_mpsrf){
      out <- unlist(lapply(gd,function(x){max(x[[1]][,1])}))
    } else{
      out <- unlist(lapply(gd,function(x){x$mpsrf}))
    }
  } else {
    gd <- gelman_diag_robust(pmwg_mcmc,autoburnin = autoburnin, transform = transform, omit_mpsrf = omit_mpsrf)
    if(omit_mpsrf){
      out <- max(gd[[1]][,1])
    } else{
      out <- gd$mpsrf
    }
  }
  if (return_summary) return(out)
  if (sort_print) out <- sort(out)
  if (print_summary) print(round(out,digits_print))
  if (selection=="alpha" || selection == "random") invisible(
    cbind(do.call(rbind,lapply(gd,function(x){x[[1]][,1]})),
          mpsrf=unlist(lapply(gd,function(x){x[[2]]})))) else
            invisible(c(gd$psrf[,1],mpsrf=gd$mpsrf))
}


#' gd_summary
#'
#' Summarizes gelman_diag statistics for a samplers object, invisibly returning
#' a list of two lists containing univarite (psrf) and multivariate (mpsrf)
#' statistics.
#'
#' @param samplers Samples object with multiple chains
#' @param no_print Boolean for printing
#' @param digits Integer, number of digits for printing
#'
#' @return List of two lists names psrf and mpsrf.
#' @export
gd_summary <- function(samplers,no_print=TRUE,digits=2) {

  alpha <- gd_pmwg(samplers,selection="alpha",print_summary = FALSE)
  alphai <- alpha; alpha <- alpha[,"mpsrf"]; alphai <- alphai[,dimnames(alphai)[[2]]!="mpsrf"]
  hierarchical <- any(names(samplers[[1]]$samples)=="theta_mu")
  if (hierarchical) {
    mu <- gd_pmwg(samplers,selection="mu",print_summary = FALSE)
    variance <- gd_pmwg(samplers,selection="variance",print_summary = FALSE)
    correlation <- gd_pmwg(samplers,selection="correlation",print_summary = FALSE)
    mui <- mu; mu <- mu["mpsrf"]; mui <- mui[names(mui)!="mpsrf"]
    variancei <- variance; variance <- variance["mpsrf"]; variancei <- variancei[names(variancei)!="mpsrf"]
    correlationi <- correlation; correlation <- correlation["mpsrf"]; correlationi <- correlationi[names(correlationi)!="mpsrf"]
  }
  if (!(no_print)) {
    cat("ALPHA psrf\n")
    print(round(alphai,digits))
    cat("\nALPHA mpsrf\n")
    print(round(sort(alpha),digits))
    if (hierarchical) {
      cat("\nMU psrf\n")
      print(round(sort(mui),digits))
      cat("\nVARIANCE psrf\n")
      print(round(sort(variancei),digits))
      cat("\nCORRELATION psrf\n")
      print(round(sort(correlationi),digits))
      cat("\nHyper mpsrf\n")
      print(round(c(mu=mu,var=variance,corr=correlation),digits=digits))
    }
  }
  if (hierarchical)
    invisible(list(psrf=list(alpha=alphai,mu=mui,variance=variancei,correlation=correlationi),
                   mpsrf=list(alpha=alpha,mu=mu,variance=variance,correlation=correlation))) else
                     invisible(list(psrf=list(alpha=alphai), mpsrf=list(alpha=alpha)))
}


iat_pmwg <- function(pmwg_mcmc,
                     print_summary=TRUE,digits_print=2,sort_print=TRUE,summary_alpha=mean,
                     selection="alpha",filter="sample",thin=1,subfilter=NULL)
  # Integrated autocorrelation time, prints multivariate summary returns each participant unless +
  # multivariate as matrix unless !return_summary
{

  IAT <- function (x,verbose=FALSE)
    # From LaplacesDemon
  {
    dt <- x
    n <- length(x)
    mu <- mean(dt)
    s2 <- var(dt)
    maxlag <- max(3, floor(n/2))
    Ga <- rep(0, 2)
    Ga[1] <- s2
    lg <- 1
    Ga[1] <- Ga[1] + sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
    m <- 1
    lg <- 2 * m
    Ga[2] <- sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
    lg <- 2 * m + 1
    Ga[2] <- Ga[2] + sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
    IAT <- Ga[1]/s2
    while ((Ga[2] > 0) & (Ga[2] < Ga[1])) {
      m <- m + 1
      if (2 * m + 1 > maxlag) {
        if (verbose) cat("Not enough data, maxlag=", maxlag, "\n")
        break
      }
      Ga[1] <- Ga[2]
      lg <- 2 * m
      Ga[2] <- sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
      lg <- 2 * m + 1
      Ga[2] <- Ga[2] + sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
      IAT <- IAT + Ga[1]/s2
    }
    IAT <- -1 + 2 * IAT
    return(IAT)
  }

  get_IAT <- function(mcs) {
    if (!is(mcs, "mcmc.list")) apply(mcs,2,IAT) else
      apply(do.call(rbind,lapply(mcs,function(x){apply(x,2,IAT)})),2,mean)
  }

  if (!(class(pmwg_mcmc[[1]]) %in% c("mcmc","mcmc.list"))) {
    if (is(pmwg_mcmc, "pmwgs"))
      pmwg_mcmc <- as_Mcmc(pmwg_mcmc,selection=selection,filter=filter,
                           thin=thin,subfilter=subfilter) else
                             pmwg_mcmc <- as_mcmc.list(pmwg_mcmc,selection=selection,filter=filter,
                                                       thin=thin,subfilter=subfilter)
  }
  if ( selection=="LL" ) stop("IAT not appropriate for LL") else
    if (selection=="alpha") {
      out <- do.call(rbind,lapply(pmwg_mcmc,get_IAT) )
      if (!is.null(summary_alpha)) out <- apply(out,2,summary_alpha)
    } else out <- get_IAT(pmwg_mcmc)
  if (sort_print & !(selection == "alpha" & is.null(summary_alpha))) out <- sort(out)
  if (print_summary) print(round(out,digits_print))
  invisible(out)
}

#### Posterior parameter tests ----


# x_fun=NULL;y=NULL;y_name=x_name;y_fun=NULL;
# mapped=FALSE;x_subject=NULL;y_subject=NULL;
# mu=0;alternative = c("less", "greater")[1];
# probs = c(0.025,.5,.975);digits=2;p_digits=3;print_table=TRUE;
# filter="sample";selection="mu";subfilter=0
#
# x=sPNAS_a;x_name="t0"; mapped=TRUE

#' Posterior parameter tests
#'
#' Modeled after t.test. For a one sample test provide x and for two sample
#' also provide y.
#'
#' @param x A pmwgs object or list of these
#' @param x_name Name of the parameter to be tested
#' @param x_fun Function applied to each iteration's parameter vector to create
#' variable to be tested.
#' @param y A pmwgs object or list of these
#' @param y_name Name of the parameter to be tested
#' @param y_fun Function applied to each iteration's parameter vector to create
#' variable to be tested.
#' @param mapped Should samples be mapped before doing test.
#' @param x_subject Integer or name selecting a subject
#' @param y_subject Integer or name selecting a subject
#' @param mu Null value for single sample test (default 0)
#' @param alternative "less" or "greater" determining direction of test probability
#' @param probs Vector defining credible interval and central tendency quatiles
#' (default c(0.025,.5,.975))
#' @param digits Integer, significant digits for estimates in printed results
#' @param p_digits Integer, significant digits for probability in printed results
#' @param print_table Boolean (default TRUE) for printing results table
#' @param selection String designating parameter type (mu, variance, correlation, alpha = default)
#' @param filter A string. Specifies which stage you want to plot.
#' @param x_fun_name Name to give to quantity calculated by x_fun
#' @param y_fun_name Name to give to quantity calculated by y_fun
#' @param subfilter An integer or vector. If integer it will exclude up until
#'
#' @return Invisible results table with no rounding.
#' @export

p_test <- function(x,x_name=NULL,x_fun=NULL,x_fun_name="fun",
                   y=NULL,y_name=NULL,y_fun=NULL,y_fun_name="fun",
                   mapped=FALSE,
                   x_subject=NULL,y_subject=NULL,
                   mu=0,alternative = c("less", "greater")[1],
                   probs = c(0.025,.5,.975),digits=2,p_digits=3,print_table=TRUE,
                   filter="sample",selection="mu",subfilter=0)

{

  get_effect <- function(x,p_name=NULL,fun=NULL)
  {
    x <- do.call(rbind,x)
    if (!is.null(fun)) return(apply(x,1,fun))
    x[,p_name]
  }


  if (mapped & !(selection %in% c("mu","alpha")))
    stop("Can only analyze mapped mu or alpha parameters")
  if (is.null(x_name) & is.null(x_fun))
    stop("x_name or x_fun must be supplied")
  if (is.null(y_fun) && is.null(y_name)) y_name <- x_name
  if (is.null(y_fun) && !is.null(x_fun)) y_fun <- x_fun

  # Process x
  if (!is(x[[1]], "pmwgs")) stop("x must be a list of pmwgs objects")
  if (length(x[[1]]$data)==1) selection <- "alpha"
  x <- as_mcmc.list(x,selection=selection,filter=filter,
                    subfilter=subfilter,mapped=mapped)
  # Individual subject analysis
  if (selection != "alpha") x_subject <- NULL else
    if (is.null(x_subject)) x_subject <- names(x)[1] else
      if (is.numeric(x_subject)) x_subject <- names(x)[x_subject]
  if (!is.null(x_subject)) {
    if (!(x_subject %in% names(x))) stop("Subject x_subject not in x")
    message("Testing x subject ",x_subject)
    x <- x[[x_subject]]
  }
  # Check test is valid
  if (is.null(x_fun) && !all(x_name %in% dimnames(x[[1]])[[2]]) )
    stop("x_name not present in samples")
  if (length(x_name)>2)
    stop("x_name must be length 1 or 2")
  # Get x effect
  x <- get_effect(x,x_name,x_fun)
  if (length(x_name)>1) {
    x <- -apply(x,1,diff)
    x_name <- paste(x_name,collapse="-")
  }
  if (is.null(x_name)) x_name <- x_fun_name
  if (is.null(y)) {
    p <- mean(x<mu)
    if (alternative=="greater") p <- 1-p
    tab <- cbind(quantile(x,probs),c(NA,mu,NA))
    attr(tab,alternative) <- p
    dimnames(tab)[[2]] <- c(x_name,"mu")
  } else {
    if (!is(y[[1]], "pmwgs")) stop("y must be a list of pmwgs objects")
    y <- as_mcmc.list(y,selection=selection,filter=filter,
                      subfilter=subfilter,mapped=mapped)
    # Individual subject analysis
    if (selection != "alpha") y_subject <- NULL else
      if (is.null(y_subject)) y_subject <- names(y)[1] else
        if (is.numeric(y_subject)) y_subject <- names(y)[y_subject]
    if (!is.null(y_subject)) {
      if (!(y_subject %in% names(y))) stop("Subject y_subject not in y")
      message("Testing y subject ",y_subject)
      y <- y[[y_subject]]
    }
    if (is.null(y_fun) && !all(y_name %in% dimnames(y[[1]])[[2]]) )
      stop("y_name not present in samples")
    if (length(y_name)>2)
      stop("y_name must be length 1 or 2")
    y <- get_effect(y,y_name,y_fun)
    if (length(y_name)>1) {
      y <- -apply(y,1,diff)
      y_name <- paste(y_name,collapse="-")
    }
    if (length(x)>length(y)) x <- x[1:length(y)] else y <- y[1:length(x)]
    d <- x-y
    p <- mean(d<0)
    if (alternative=="greater") p <- 1-p
    tab <- cbind(quantile(x,probs),quantile(y,probs),quantile(d,probs))
    attr(tab,alternative) <- p
    if (is.null(y_name)) y_name <- y_fun_name
    if (x_name==y_name)
      dimnames(tab)[[2]] <- c(paste(x_name,c(x_subject,y_subject),sep="_"),
                              paste(x_subject,y_subject,sep="-")) else
                                dimnames(tab)[[2]] <- c(x_name,y_name,paste(x_name,y_name,sep="-"))
  }
  if (print_table) {
    ptab <- tab
    ptab <- round(ptab,digits)
    attr(ptab,alternative) <- round(attr(ptab,alternative),p_digits)
    print(ptab)
  }
  invisible(tab)
}

#' Calculate information criteria (DIC, BPIC), effective number of parameters and
#' constituent posterior deviance (D) summaries (meanD = mean of D, Dmean = D
#' for mean of posterior parameters and minD = minimum of D).
#'
#' @param samplers pmwgs object or list of these
#' @param filter A string. Specifies which stage you want to plot.
#' @param subfilter An integer or vector. If integer it will exclude up until
#' @param use_best_fit Boolean, default TRUE use best of minD and Dmean in
#' calculation otherwise always use Dmean
#' @param print_summary Boolean (default TRUE) print table of results
#' @param digits Integer, significant digits in printed table
#' @param subject Integer or string selecting a single subject, default NULL
#' returns sums over all subjects
#'
#' @return Table of DIC, BPIC, EffectiveN, meanD, Dmean, and minD
#' @export

IC <- function(samplers,filter="sample",subfilter=0,use_best_fit=TRUE,
               print_summary=TRUE,digits=0,subject=NULL)
  # Gets DIC, BPIC, effective parameters, mean deviance, and deviance of mean
{


  # Mean log-likelihood for each subject
  if (is(samplers, "pmwgs"))
    ll <- as_Mcmc(samplers,selection="LL",filter=filter,subfilter=subfilter) else
      ll <- as_mcmc.list(samplers,selection="LL",filter=filter,subfilter=subfilter)
    minDs <- -2*unlist(lapply(ll,function(x){max(unlist(x))}))
    mean_lls <- unlist(lapply(ll,function(x){mean(unlist(x))}))

    if (is(samplers, "pmwgs"))
      alpha <- as_Mcmc(samplers,selection="alpha",filter=filter,subfilter=subfilter) else
        alpha <- as_mcmc.list(samplers,selection="alpha",filter=filter,subfilter=subfilter)
    mean_pars <- lapply(alpha,function(x){apply(do.call(rbind,x),2,mean)})
    # log-likelihood for each subject using their mean parameter vector
    ll_func <- attr(samplers,"design_list")[[1]]$model()$log_likelihood
    data <- samplers[[1]]$data
    mean_pars_lls <- setNames(numeric(length(mean_pars)),names(mean_pars))
    for (sub in names(mean_pars))
      mean_pars_lls[sub] <- ll_func(mean_pars[[sub]],dadm = data[[sub]])
    Dmeans <- -2*mean_pars_lls
    if (use_best_fit) minDs <- pmin(minDs,Dmeans)

    if (!is.null(subject)) {
      Dmeans <- Dmeans[subject[1]]
      mean_lls <- mean_lls[subject[1]]
      minDs <- minDs[subject[1]]
    }

    # mean deviance(-2*ll of all data)
    mD <- sum(-2 * mean_lls)
    # Deviance of mean
    Dmean <- sum(Dmeans)
    # mimimum Deviance
    minD <- sum(minDs)

    # Use deviance of mean as best fit or use actual best fit
    if (!use_best_fit) Dm <- Dmean else Dm <- minD

    # effective number of parameters
    pD <- mD - Dm
    # DIC = mean deviance + effective number of parameters
    DIC <- mD + pD
    # BPIC = mean deviance + 2*effective number of parameters
    # Note this is the "easy" BPIC, instead of the complex 2007 one
    BPIC <- mD + 2*pD
    out <- c(DIC = DIC, BPIC = BPIC, EffectiveN = pD,meanD=mD,Dmean=Dmean,minD=minD)
    names(out) <- c("DIC","BPIC","EffectiveN","meanD","Dmean","minD")
    if (print_summary) print(round(out,digits))
    invisible(out)
}


#' Bayes Factors
#'
#' returns the Bayes Factor for two models
#'
#' @param MLL1 Numeric. Marginal likelihood of model 1. Obtained with `run_bridge_sampling`
#' @param MLL2 Numeric. Marginal likelihood of model 2.
#'
#' @return The BayesFactor for model 1 over model 2
#' @examples \dontrun{
#' # First get the marginal likelihood for two_models
#' # Here the full model is a list of samplers with the hypothesized effect
#' # The null model is a list of samplers without the hypothesized effect
#' MLL_full <- run_bridge_sampling(full_model, cores_per_prop = 3)
#' MLL_null <- run_bridge_sampling(null_model, cores_per_prop = 3)
#' # Now we can calculate their Bayes factor
#' get_BayesFactor(MLL_full, MLL_null)
#' }
#' @export
get_BayesFactor <- function(MLL1, MLL2){
  exp(MLL1 - MLL2)
}



#' Information criteria and marginal likelihoods
#'
#' Returns the BPIC/DIC or marginal deviance (-2*marginal likelihood) for a list of samples objects.
#'
#' @param sList List of samples objects
#' @param filter A string. Specifies which stage the samples are to be taken from `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`
#' @param subfilter An integer or vector. If an integer, it will exclude up until the integer within the filter.
#' @param use_best_fit Boolean, defaults to `TRUE`, uses the minimal or mean likelihood (whichever is better) in the
#' calculation, otherwise always uses the mean likelihood.
#' @param BayesFactor Boolean, defaults to `TRUE`. Include marginal likelihoods as estimated using WARP-III bridge sampling.
#' Usually takes a minute per model added to calculate
#' @param cores_for_props Integer, how many cores to use for the Bayes factor calculation, here 4 is the default for the 4 different proposal densities to evaluate, only 1, 2 and 4 are sensible.
#' @param cores_per_prop Integer, how many cores to use for the Bayes factor calculation if you have more than 4 cores available. Cores used will be cores_for_props * cores_per_prop. Best to prioritize cores_for_props being 4 or 2
#' @param print_summary Boolean (default `TRUE`), print table of results
#' @param digits Integer, significant digits in printed table except model weights
#' @param digits_p Integer, significant digits in printed table for model weights
#' @param ... Additional, optional arguments
#'
#' @return Matrix of effective number of parameters, mean deviance, deviance of
#' mean, DIC, BPIC, Marginal Deviance (if `BayesFactor=TRUE`) and associated weights.
#' @examples \dontrun{
#' # Define a samplers list of two (or more different models)
#' # Here the full model is a list of samplers with the hypothesized effect
#' # The null model is a list of samplers without the hypothesized effect
#' sList <- list(full_model, null_model)
#' # By default emc uses 4 cores to parallelize marginal likelihood estimation across proposals
#' # So cores_per_prop = 3 results in 12 cores used.
#' compare(sList, cores_per_prop = 3)
#' }
#' @export

compare <- function(sList,filter="sample",subfilter=0,use_best_fit=TRUE,
                    BayesFactor = TRUE, cores_for_props =4, cores_per_prop = 1,
                    print_summary=TRUE,digits=0,digits_p=3, ...) {
  if(is(sList, "emc")) sList <- list(sList)
  getp <- function(IC) {
    IC <- -(IC - min(IC))/2
    exp(IC)/sum(exp(IC))
  }
  if(!is.null(list(...)$subject)) subject <- list(...)$subject

  if (is.numeric(subfilter)) defaultsf <- subfilter[1] else defaultsf <- 0
  sflist <- as.list(setNames(rep(defaultsf,length(sList)),names(sList)))
  if (is.list(subfilter)) for (i in names(subfilter))
    if (i %in% names(sflist)) sflist[[i]] <- subfilter[[i]]

  ICs <- setNames(vector(mode="list",length=length(sList)),names(sList))
  for (i in 1:length(ICs)) ICs[[i]] <- IC(sList[[i]],filter=filter,
                                          subfilter=sflist[[i]],use_best_fit=use_best_fit,subject=NULL,print_summary=FALSE)
  ICs <- data.frame(do.call(rbind,ICs))
  DICp <- getp(ICs$DIC)
  BPICp <- getp(ICs$BPIC)
  out <- cbind.data.frame(DIC=ICs$DIC,wDIC=DICp,BPIC=ICs$BPIC,wBPIC=BPICp,ICs[,-c(1:2)])

  if(BayesFactor){
    MLLs <- numeric(length(sList))
    for(i in 1:length(MLLs)){
      MLLs[i] <- run_bridge_sampling(sList[[i]], filter = filter, subfilter = sflist[[i]], both_splits = FALSE,
                                     cores_for_props = cores_for_props, cores_per_prop = cores_per_prop)
    }
    MD <- -2*MLLs
    modelProbability <- getp(MD)
    out <- cbind.data.frame(MD = MD, wMD = modelProbability, out)
  }
  if (print_summary) {
    tmp <- out
    tmp$wDIC <- round(tmp$wDIC,digits_p)
    tmp$wBPIC <- round(tmp$wBPIC,digits_p)
    if(BayesFactor){
      tmp$wMD <- round(tmp$wMD, digits_p)
      tmp[,-c(2,4,6)] <- round(tmp[,-c(2,4,6)],digits=digits)
    } else{
      tmp[,-c(2,4)] <- round(tmp[,-c(2,4)],digits=digits)
    }
    print(tmp)
  }
  invisible(out)
}

#' The savage dickey ratio.
#'
#' Can be used to approximate the Bayes Factor for group level mean effects.
#' Note this is different to MD in `compare` since it only considers the group level
#' mean effect and not the whole model. For details see: Wagenmakers, Lodewyckx, Kuriyal, & Grasman (2010).
#'
#' @param samplers A list. A samplers object.
#' @param parameter A string. A parameter which you want to compare to H0. Will not be used if a FUN is specified.
#' @param H0 An integer. The H0 value which you want to compare to
#' @param filter A string. Specifies which stage the samples are to be taken from "preburn", "burn", "adapt", or "sample"
#' @param subfilter An integer or vector. If integer it will exclude up until that integer.
#' @param fun A function. Specifies an operation to be performed on the sampled or mapped parameters.
#' @param mapped A boolean. Whether the BF should be calculated for parameters mapped back to the real design, only works with selection = 'mu'.
#' @param selection A string. Default is mu. Whether to do the operation on the alpha, mu, covariance, variance, or correlations.
#' @param do_plot Boolean. Whether to include a plot of the prior and posterior density. With circles at H0.
#' @param xlim Vector, the xlimits for the plot.
#' @param subject Character. If type = "single" and multiple subjects were ran in one model, this is required.
#'
#' @return The BayesFactor for the hypothesis against H0.
#' @examples \dontrun{
#' # Here the samplers object has an effect parameter (e.g. B_Eneutral),
#' # that maps onto a certain hypothesis.
#' # The hypothesis here is that B_Eneutral is different from zero.
#' # We can test whether there's a group-level effect of Eneutral on B:
#' savage_dickey(samplers, parameter = "B_Eneutral")
#' # Alternatively we can also test whether two parameters from each other
#' Bdiff <- function(p)diff(p[c("B_Eneutral","B_Eaccuracy")])
#' savage_dickey(samplers,fun=Bdiff)
#' }
#' @export
savage_dickey <- function(samplers, parameter = NULL, H0 = 0, filter = "sample",
                          subfilter = 0, fun = NULL, mapped =F, selection = "mu",
                          do_plot = TRUE, xlim = NULL, subject = NULL){
  if(mapped & selection != "mu") stop("Mapped only works for mu")
  prior <- samplers[[1]]$prior
  type <- attr(samplers[[1]], "variant_funs")$type
  if(selection == "alpha" & type != "single") stop("For savage-dickey ratio, selection cannot be alpha")
  # type <- "standard"
  if(type == "standard") gp <- get_prior_standard
  if(type == "diagonal") gp <- get_prior_blocked
  if(type == "single") gp <- get_prior_single
  if(type == "blocked") gp <- get_prior_blocked
  psamples <- gp(prior = prior, type = selection)[[selection]]
  if(mapped){
    design <- attr(samplers, "design_list")[[1]]
    psamples <- map_mcmc(psamples, design, design$model,
                         include_constants = FALSE)
  }

  samples <- as_mcmc.list(samplers,selection=selection,filter=filter,
                          subfilter=subfilter,mapped=mapped)
  if(selection == "alpha"){
    if(length(samples) > 1 & is.null(subject)){
      stop("with non-hierarichal run with multiple subjects, you must specify which subject")
    } else if (length(samples) == 1){
      samples <- do.call(rbind, samples[[1]])
    } else{
      samples <- do.call(rbind, samples[[subject]])
    }
  } else{
    samples <- do.call(rbind, samples)
  }
  if(is.null(fun)){
    idx <- colnames(samples) == parameter
    samples <- samples[,idx]
    psamples <- psamples[,idx]
  } else{
    colnames(psamples) <- colnames(samples)
    samples <- apply(as.data.frame(samples), 1, fun)
    psamples <- apply(as.data.frame(psamples), 1, fun)
  }
  min_bound <- min(min(psamples), H0)
  max_bound <- max(max(psamples), H0)
  diff <- max_bound - min_bound
  pdensity <- density(psamples, bw = "sj", from = min_bound - diff/2, to = max_bound + diff/2)
  pdfun <-approxfun(pdensity)

  min_bound <- min(min(samples), H0)
  max_bound <- max(max(samples), H0)
  diff <- max_bound - min_bound
  post_density <- density(samples, bw = "sj", from = min_bound - diff/2, to = max_bound + diff/2)
  post_dfun <-approxfun(post_density)
  if(do_plot){
    if(is.null(xlim)){
      xmin <- min(quantile(samples, 0.025), quantile(psamples, 0.025))
      xmax <- max(quantile(samples, 0.975), quantile(psamples, 0.975))
      xlim <- c(xmin, xmax)
    }
    plot(post_density, xlim = xlim, lwd = 1.5, main = "Prior and posterior density")
    lines(pdensity, col = "red", lwd = 1.5)
    points(H0, post_dfun(H0), cex = 2)
    points(H0, pdfun(H0), col = "red", cex = 2)
    legend("topright", legend = c("posterior", "prior"), pch = c(1, 1), col = c("black", "red"))
  }
  return(pdfun(H0)/post_dfun(H0))
}

#' Information criteria for each participant
#'
#' Returns the BPIC/DIC based model weights for each participant in a list of samples objects
#'
#' @param sList List of samples objects
#' @param filter A string. Specifies which stage the samples are to be taken from "preburn", "burn", "adapt", or "sample"
#' @param subfilter An integer or vector. If integer it will exclude up until
#' @param use_best_fit Boolean, default TRUE use best of minimal likelihood and mean likelihood in
#' calculation otherwise always use mean likelihood.
#' @param print_summary Boolean (default TRUE) print table of results
#' @param digits Integer, significant digits in printed table
#'
#' @return List of matrices for each subject of effective number of parameters,
#' mean deviance, deviance of mean, DIC, BPIC and associated weights.
#' @examples \dontrun{
#' # Define a samplers list of two (or more different models)
#' # Here the full model is a list of samplers with the hypothesized effect
#' # The null model is a list of samplers without the hypothesized effect
#' sList <- list(full_model, null_model)
#' # By default emc uses 4 cores to parallelize marginal likelihood estimation across proposals
#' # So cores_per_prop = 3 results in 12 cores used.
#' compare_subject(sList, cores_per_prop = 3)
#' # prints a set of weights for each model for the different participants
#' # And returns the DIC and BPIC for each participant for each model.
#' }
#' @export
compare_subject <- function(sList,filter="sample",subfilter=0,use_best_fit=TRUE,
                            print_summary=TRUE,digits=3) {
  if(is(sList, "emc")) sList <- list(sList)
  subjects <- names(sList[[1]][[1]]$data)
  out <- setNames(vector(mode="list",length=length(subjects)),subjects)
  for (i in subjects) out[[i]] <- compare(sList,subject=i,BayesFactor=FALSE,
                                          filter=filter,subfilter=subfilter,use_best_fit=use_best_fit,print_summary=FALSE)
  if (print_summary) {
    wDIC <- lapply(out,function(x)x["wDIC"])
    wBPIC <- lapply(out,function(x)x["wBPIC"])
    pDIC <- do.call(rbind,lapply(wDIC,function(x){
      setNames(data.frame(t(x)),paste("wDIC",rownames(x),sep="_"))}))
    pBPIC <- do.call(rbind,lapply(wBPIC,function(x){
      setNames(data.frame(t(x)),paste("wBPIC",rownames(x),sep="_"))}))
    # if (BayesFactor) {
    #   pMD <- do.call(rbind,lapply(wMD,function(x){
    #   setNames(data.frame(t(x)),paste("wMD",rownames(x),sep="_"))}))
    #   print(round(cbind(pDIC,pBPIC,pMD),digits))
    #   mnams <- unlist(lapply(strsplit(dimnames(pDIC)[[2]],"_"),function(x){x[[2]]}))
    #   cat("\nWinners\n")
    #   print(rbind(DIC=table(mnams[apply(pDIC,1,which.max)]),
    #               BPIC=table(mnams[apply(pBPIC,1,which.max)]),
    #               MD=table(mnams[apply(pMD,1,which.max)])))
    #
    # } else {
    print(round(cbind(pDIC,pBPIC),digits))
    mnams <- unlist(lapply(strsplit(dimnames(pDIC)[[2]],"_"),function(x){x[[2]]}))
    cat("\nWinners\n")
    print(rbind(DIC=table(mnams[apply(pDIC,1,which.max)]),
                BPIC=table(mnams[apply(pBPIC,1,which.max)])))
    # }
  }
  invisible(out)
}


#' Calculate a table of model probabilities based for a list of samples objects
#' based on samples of marginal log-likelihood (MLL) added to these objects by
#' run_IS2. Probabilities estimated by a bootstrap ath picks a vector of MLLs,
#' one for each model in the list randomly with replacement nboot times,
#' calculates model probabilities and averages
#'
#'
#' @param mll List of samples objects with IS_samples attribute added by by run_IS2
#' @param nboot Integer number of bootstrap samples, the default (1e5) usually
#' gives stable results at 2 decimal places.
#' @param print_summary Boolean (default TRUE) print table of results
#' @param digits Integer, significant digits in printed table
#'
#' @return Vector of model probabilities with names from samples list.
#' @export

compare_MLL <- function(mll,nboot=1e5,digits=2,print_summary=TRUE)
  # mll is a list of vectors of marginal log-likelihoods for a set of models
  # picks a vector of mlls for each model in the list randomly with replacement
  # nboot times, calculates model probabilities and averages.
{
  pmp <- function(x)
    # posterior model probability for a vector of marginal log-likelihoods
  {
    x <- exp(x-max(x))
    x/sum(x)
  }

  out <- sort(apply(apply(do.call(rbind,lapply(mll,function(x){
    attr(x,"IS_samples")[sample(length(attr(x,"IS_samples")),nboot,replace=TRUE)]
  })),2,pmp),1,mean),decreasing=TRUE)
  print(round(out,digits))
  invisible(out)
}



condMVN <- function (mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE)
{
  if (missing(dependent.ind))
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given))
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(given.ind) == 0)
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(X.given) != length(given.ind))
    stop("lengths of `X.given' and `given.ind' must be same")
  if (check.sigma) {
    if (!isSymmetric(sigma))
      stop("sigma is not a symmetric matrix")
    eigenvalues <- eigen(sigma, only.values = TRUE)$values
    if (any(eigenvalues < 1e-08))
      stop("sigma is not positive-definite")
  }
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% chol2inv(chol(D))
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean = cMu, condVar = cVar)
}
