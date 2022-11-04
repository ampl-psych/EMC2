std_error_IS2 <- function(IS_samples, n_bootstrap = 50000){
  log_marglik_boot= array(dim = n_bootstrap)
  for (i in 1:n_bootstrap){
    log_weight_boot = sample(IS_samples, length(IS_samples), replace = TRUE) #resample with replacement from the lw
    log_marglik_boot[i] <- median(log_weight_boot)
  }
  return(sd(log_marglik_boot))
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


# return_summary=FALSE;print_summary=TRUE;digits_print=2;sort_print=TRUE;
# autoburnin=FALSE;transform=TRUE
# selection="alpha";filter="sample";thin=1;subfilter=NULL
# pmwg_mcmc=sVat0;selection="mu";filter="sample"
gd_pmwg <- function(pmwg_mcmc,return_summary=FALSE,print_summary=TRUE,
                    digits_print=2,sort_print=TRUE,autoburnin=FALSE,transform=TRUE,
                    selection="alpha",filter="sample",thin=1,subfilter=NULL,mapped=FALSE)
  # R hat, prints multivariate summary returns each participant unless +
  # multivariate as matrix unless !return_summary
{

  gelman_diag_robust <- function(mcl,autoburnin,transform)
  {
    gd <- try(gelman.diag(mcl,autoburnin=autoburnin,transform=transform),silent=TRUE)
    if (is(gd, "try-error")) list(psrf=matrix(Inf),mpsrf=Inf) else gd
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
  if (selection=="alpha") {
    gd <- lapply(pmwg_mcmc,gelman_diag_robust,autoburnin = autoburnin, transform = transform)
    out <- unlist(lapply(gd,function(x){x$mpsrf}))
  } else {
    gd <- gelman_diag_robust(pmwg_mcmc,autoburnin = autoburnin, transform = transform)
    out <- gd$mpsrf
  }
  if (return_summary) return(out)
  if (sort_print) out <- sort(out)
  if (print_summary) print(round(out,digits_print))
  if (selection=="alpha") invisible(
    cbind(do.call(rbind,lapply(gd,function(x){x[[1]][,1]})),
          mpsrf=unlist(lapply(gd,function(x){x[[2]]})))) else
            invisible(c(gd$psrf[,1],mpsrf=gd$mpsrf))
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
p_test <- function(x,x_name,x_fun=NULL,
                   y=NULL,y_name=x_name,y_fun=NULL,
                   mapped=FALSE,
                   x_subject=NULL,y_subject=NULL,
                   mu=0,alternative = c("less", "greater")[1],
                   probs = c(0.025,.5,.975),digits=2,p_digits=3,print_table=TRUE,
                   filter="sample",selection="mu",subfilter=0)

{

  get_effect <- function(x,p_name=NULL,fun=NULL)
    # Effect, must always be on mapped scale if lc_mat supplied.
  {
    x <- do.call(rbind,x)
    if (!is.null(fun)) return(apply(x,1,fun))
    x[,p_name]
  }


  if (mapped & !(selection %in% c("mu","alpha")))
    stop("Can only analyze mapped mu or alpha parameters")

  # Process x
  if (!is(x[[1]], "pmwgs")) stop("x must be a list of pmwgs objects")
  x <- as_mcmc.list(x,selection=selection,filter=filter,
                    subfilter=subfilter,mapped=mapped)
  # Individual subject analysis
  if (selection != "alpha") x_subject <- NULL else
    if (is.null(x_subject)) x_subject <- names(x)[1]
  if (!is.null(x_subject)) {
    if (!(x_subject %in% names(x))) stop("Subject x_subject not in x")
    message("Testing x subject ",x_subject)
    x <- x[[x_subject]]
  }
  # Check test is valid
  if (is.null(x_fun) && !(x_name %in% dimnames(x[[1]])[[2]]) )
    stop("x_name not present in samples")
  # Get x effect
  x <- get_effect(x,x_name,x_fun)
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
      if (is.null(y_subject)) y_subject <- names(y)[1]
    if (!is.null(y_subject)) {
      if (!(y_subject %in% names(y))) stop("Subject y_subject not in y")
      message("Testing y subject ",y_subject)
      y <- y[[y_subject]]
    }
    if (is.null(y_fun) && !(y_name %in% dimnames(y[[1]])[[2]]) )
      stop("y_name not present in samples")
    y <- get_effect(y,y_name,y_fun)
    if (length(x)>length(y)) x <- x[1:length(y)] else y <- y[1:length(x)]
    d <- x-y
    p <- mean(d<0)
    if (alternative=="greater") p <- 1-p
    tab <- cbind(quantile(x,probs),quantile(y,probs),quantile(d,probs))
    attr(tab,alternative) <- p
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


# ptype="mean"
# selection="mu";Flist=NULL;Clist=NULL
# p_tests <- function(samples,ptype,selection="mu",Flist=NULL,Clist=NULL)
#   # Performs multiple tests on mapped parameters
# {
#   if (!(selection %in% c("mu","alpha"))
#     stop("Can only analyze mapped mu or alpha parameters")
#
# }

# filter="sample";subfilter=0;use_best_fit=FALSE;print_summary=FALSE;digits=0
# filter="sample"
pmwg_IC <- function(samplers,filter="sample",subfilter=0,use_best_fit=TRUE,
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
    ll_func <- attr(samplers,"design_list")[[1]]$model$log_likelihood
    data <- samplers[[1]]$data
    mean_pars_lls <- setNames(numeric(length(mean_pars)),names(mean_pars))
    for (sub in names(mean_pars))
      mean_pars_lls[sub] <- ll_func(mean_pars[[sub]],dadm = data[[sub]])
    Dmeans <- -2*mean_pars_lls
    if (use_best_fit) minDs <- pmin(minDs,Dmeans)

    if (!is.null(subject)) {
      Dmeans <- Dmeans[subject]
      mean_lls <- mean_lls[subject]
      minDs <- minDs[subject]
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



compare_IC <- function(sList,filter="sample",subfilter=0,use_best_fit=TRUE,
                       print_summary=TRUE,digits=0,digits_p=3,subject=NULL) {

  getp <- function(IC) {
    IC <- -(IC - min(IC))/2
    exp(IC)/sum(exp(IC))
  }

  if (length(subfilter)==1) {
    tmp <- vector(mode="list",length(sList))
    for (i in 1:length(sList)) tmp[[i]] <- subfilter
    subfilter=tmp
  }
  if ( !is.list(subfilter) || length(subfilter)!=length(sList) )
    stop("If not a single digit, subfilter must be a list of the same length as sList")
  ICs <- setNames(vector(mode="list",length=length(sList)),names(sList))
  for (i in 1:length(ICs)) ICs[[i]] <- pmwg_IC(sList[[i]],filter=filter,
                                               subfilter=subfilter[[i]],use_best_fit=use_best_fit,subject=subject,print_summary=FALSE)
  ICs <- data.frame(do.call(rbind,ICs))
  DICp <- getp(ICs$DIC)
  BPICp <- getp(ICs$BPIC)
  out <- cbind.data.frame(DIC=ICs$DIC,wDIC=DICp,BPIC=ICs$BPIC,wBPIC=BPICp,ICs[,-c(1:2)])
  if (print_summary) {
    tmp <- out
    tmp$wDIC <- round(tmp$wDIC,digits_p)
    tmp$wBPIC <- round(tmp$wBPIC,digits_p)
    tmp[,-c(2,4)] <- round(tmp[,-c(2,4)])
    print(tmp)
  }
  invisible(out)
}

compare_ICs <- function(sList,filter="sample",subfilter=0,use_best_fit=TRUE,
                        print_summary=TRUE,digits=3,subject=NULL) {

  if (length(subfilter)==1) {
    tmp <- vector(mode="list",length(sList))
    for (i in 1:length(sList)) tmp[[i]] <- subfilter
    subfilter=tmp
  }
  if ( !is.list(subfilter) || length(subfilter)!=length(sList) )
    stop("If not a single digit, subfilter must be a list of the same length as sList")
  subjects <- names(sList[[1]][[1]]$data)
  out <- setNames(vector(mode="list",length=length(subjects)),subjects)
  for (i in subjects) out[[i]] <- compare_IC(sList,subject=i,
                                             filter=filter,subfilter=subfilter,use_best_fit=use_best_fit,print_summary=FALSE)
  if (print_summary) {
    wDIC <- lapply(out,function(x)x["wDIC"])
    wBPIC <- lapply(out,function(x)x["wBPIC"])
    pDIC <- do.call(rbind,lapply(wDIC,function(x){
      setNames(data.frame(t(x)),paste("wDIC",rownames(x),sep="_"))}))
    pBPIC <- do.call(rbind,lapply(wBPIC,function(x){
      setNames(data.frame(t(x)),paste("wBPIC",rownames(x),sep="_"))}))
    print(round(cbind(pDIC,pBPIC),digits))
    mnams <- unlist(lapply(strsplit(dimnames(pDIC)[[2]],"_"),function(x){x[[2]]}))
    cat("\nWinners\n")
    print(rbind(DIC=table(mnams[apply(pDIC,1,which.max)]),
                BPIC=table(mnams[apply(pBPIC,1,which.max)])))
  }
  invisible(out)
}


compare_MLL <- function(mll,nboot=100000,digits=2,print_summary=TRUE)
  # mll is a list of vectors of marginal log-likelihoods for a set of models
  # picks a vector of mlls for each model in the list randomly with replacement
  # nboot times, calculates model probabilities and averages, the default
  # nboot seems good for stable results at 2 decimal places.
{
  pmp <- function(x)
    # posterior model probability for a vector of marginal log-likelihoods
  {
    x <- exp(x-max(x))
    x/sum(x)
  }

  out <- sort(apply(apply(do.call(rbind,lapply(mll,function(x){
    x[sample(length(x),nboot,replace=TRUE)]})),2,pmp),1,mean),decreasing=TRUE)
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
