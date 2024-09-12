#' Information criteria and marginal likelihoods
#'
#' Returns the BPIC/DIC or marginal deviance (-2*marginal likelihood) for a list of samples objects.
#'
#' @param sList List of samples objects
#' @param stage A string. Specifies which stage the samples are to be taken from `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`
#' @param filter An integer or vector. If it's an integer, iterations up until the value set by `filter` will be excluded.
#' If a vector is supplied, only the iterations in the vector will be considered.
#' @param use_best_fit Boolean, defaults to `TRUE`, uses the minimal or mean likelihood (whichever is better) in the
#' calculation, otherwise always uses the mean likelihood.
#' @param BayesFactor Boolean, defaults to `TRUE`. Include marginal likelihoods as estimated using WARP-III bridge sampling.
#' Usually takes a minute per model added to calculate
#' @param cores_for_props Integer, how many cores to use for the Bayes factor calculation, here 4 is the default for the 4 different proposal densities to evaluate, only 1, 2 and 4 are sensible.
#' @param cores_per_prop Integer, how many cores to use for the Bayes factor calculation if you have more than 4 cores available. Cores used will be cores_for_props * cores_per_prop. Best to prioritize cores_for_props being 4 or 2
#' @param print_summary Boolean (default `TRUE`), print table of results
#' @param digits Integer, significant digits in printed table for information criteria
#' @param digits_p Integer, significant digits in printed table for model weights
#' @param ... Additional, optional arguments
#'
#' @return Matrix of effective number of parameters, mean deviance, deviance of
#' mean, DIC, BPIC, Marginal Deviance (if `BayesFactor=TRUE`) and associated weights.
#' @examples \dontrun{
#' # Define a list of two (or more different models)
#' # Here the full model is an emc object with the hypothesized effect
#' # The null model is an emc object without the hypothesized effect
#' design_full <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Now without a ~ E
#' design_null <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#'
#' full_model <- make_emc(forstmann, design_full)
#' full_model <- fit(full_model)
#'
#' null_model <- make_emc(forstmann, design_null)
#' null_model <- fit(null_model)
#' sList <- list(full_model, null_model)
#' # By default emc uses 4 cores to parallelize marginal likelihood estimation across proposals
#' # So cores_per_prop = 3 results in 12 cores used.
#' compare(sList, cores_per_prop = 3)
#' }
#' @export

compare <- function(sList,stage="sample",filter=NULL,use_best_fit=TRUE,
                        BayesFactor = TRUE, cores_for_props =4, cores_per_prop = 1,
                        print_summary=TRUE,digits=0,digits_p=3, ...) {
  if(is(sList, "emc")) sList <- list(sList)
  getp <- function(IC) {
    IC <- -(IC - min(IC))/2
    exp(IC)/sum(exp(IC))
  }
  if(!is.null(list(...)$subject)) subject <- list(...)$subject

  if (is.numeric(filter)) defaultsf <- filter[1] else defaultsf <- 0
  sflist <- as.list(setNames(rep(defaultsf,length(sList)),names(sList)))
  if (is.list(filter)) for (i in names(filter))
    if (i %in% names(sflist)) sflist[[i]] <- filter[[i]]
  dots <- add_defaults(list(...), group_only = FALSE)
  ICs <- setNames(vector(mode="list",length=length(sList)),names(sList))
  for (i in 1:length(ICs)) ICs[[i]] <- IC(sList[[i]],stage=stage,
                                          filter=sflist[[i]],use_best_fit=use_best_fit,subject=NULL,print_summary=FALSE,
                                          group_only = dots$group_only)
  ICs <- data.frame(do.call(rbind,ICs))
  DICp <- getp(ICs$DIC)
  BPICp <- getp(ICs$BPIC)
  out <- cbind.data.frame(DIC=ICs$DIC,wDIC=DICp,BPIC=ICs$BPIC,wBPIC=BPICp,ICs[,-c(1:2)])

  if(BayesFactor){
    MLLs <- numeric(length(sList))
    for(i in 1:length(MLLs)){
      MLLs[i] <- run_bridge_sampling(sList[[i]], stage = stage, filter = sflist[[i]], both_splits = FALSE,
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

std_error_IS2 <- function(IS_samples, n_bootstrap = 50000){
  log_marglik_boot= array(dim = n_bootstrap)
  for (i in 1:n_bootstrap){
    log_weight_boot = sample(IS_samples, length(IS_samples), replace = TRUE) #resample with replacement from the lw
    log_marglik_boot[i] <- median(log_weight_boot)
  }
  return(sd(log_marglik_boot))
}


# robust_diwish <- function (W, v, S) { #RJI_change: this function is to protect against weird proposals in the diwish function, where sometimes matrices weren't pos def
#   if (!is.matrix(S)) S <- matrix(S)
#   if (!is.matrix(W)) W <- matrix(W)
#   p <- nrow(S)
#   gammapart <- sum(lgamma((v + 1 - 1:p)/2))
#   ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
#   if (corpcor::is.positive.definite(W, tol=1e-8)){
#     cholW<-base::chol(W)
#   }else{
#     return(1e-10)
#   }
#   if (corpcor::is.positive.definite(S, tol=1e-8)){
#     cholS <- base::chol(S)
#   }else{
#     return(1e-10)
#   }
#   halflogdetS <- sum(log(diag(cholS)))
#   halflogdetW <- sum(log(diag(cholW)))
#   invW <- chol2inv(cholW)
#   exptrace <- sum(S * invW)
#   lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * exptrace
#   lpdf <- lnum - ldenom
#   out <- exp(lpdf)
#   if(!is.finite(out)) return(1e-100)
#   if(out < 1e-10) return(1e-100)
#   return(exp(lpdf))
# }

robust_diwish <- function (W, v, S) { #RJI_change: this function is to protect against weird proposals in the diwish function, where sometimes matrices weren't pos def
  if (!is.matrix(S)) S <- matrix(S)
  if (!is.matrix(W)) W <- matrix(W)
  p <- nrow(S)
  gammapart <- sum(lgamma((v + 1 - 1:p)/2))
  ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
  cholW <- base::chol(nearPD(W)$mat)
  cholS <- base::chol(nearPD(S)$mat)
  halflogdetS <- sum(log(diag(cholS)))
  halflogdetW <- sum(log(diag(cholW)))
  invW <- chol2inv(cholW)
  exptrace <- sum(S * invW)
  lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * exptrace
  lpdf <- lnum - ldenom
  out <- exp(lpdf)
  if(!is.finite(out)) return(1e-100)
  return(out)
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

rwish <- function(v, S){
  if (!is.matrix(S))
    S <- matrix(S)
  if (nrow(S) != ncol(S)) {
    stop(message = "S not square in rwish().\n")
  }
  if (v < nrow(S)) {
    stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * (p - 1)/2)
  }
  return(crossprod(Z %*% CC))
}


riwish <- function(v, S){
  return(solve(rwish(v, solve(S))))
}

logdinvGamma <- function(x, shape, rate){
  dgamma(1/x, shape, rate, log = TRUE) - 2 * log(x)
}

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

gelman_diag_robust <- function(mcl,autoburnin = FALSE,transform = TRUE, omit_mpsrf = TRUE)
{
  mcl <- split_mcl(mcl)
  gd <- try(gelman.diag(mcl,autoburnin=autoburnin,transform=transform, multivariate = !omit_mpsrf),silent=TRUE)
  gd_out <- gd[[1]][,1] # Remove CI
  if(!omit_mpsrf){
    gd_out <- c(gd_out, gd$mpsrf)
    names(gd_out)[length(gd_out)] <- "mpsrf"
  }

  if (is(gd, "try-error")){
    if(omit_mpsrf){
      return(list(psrf=matrix(Inf)))
    } else{
      return(list(psrf=matrix(Inf),mpsrf=Inf))
    }
  } else{
    return(gd_out)
  }
}

#' Calculate information criteria (DIC, BPIC), effective number of parameters and
#' constituent posterior deviance (D) summaries (meanD = mean of D, Dmean = D
#' for mean of posterior parameters and minD = minimum of D).
#'
#' @param emc emc object or list of these
#' @param stage A string. Specifies which stage you want to plot.
#' @param filter An integer or vector. If it's an integer, iterations up until the value set by `filter` will be excluded.
#' If a vector is supplied, only the iterations in the vector will be considered.
#' @param use_best_fit Boolean, default TRUE use best of minD and Dmean in
#' calculation otherwise always use Dmean
#' @param print_summary Boolean (default TRUE) print table of results
#' @param digits Integer, significant digits in printed table
#' @param subject Integer or string selecting a single subject, default NULL
#' returns sums over all subjects
#' @param group_only Boolean. If `TRUE` will calculate the IC for the group-level only
#'
#' @return Table of DIC, BPIC, EffectiveN, meanD, Dmean, and minD

IC <- function(emc,stage="sample",filter=0,use_best_fit=TRUE,
               print_summary=TRUE,digits=0,subject=NULL,
               group_only = FALSE)
  # Gets DIC, BPIC, effective parameters, mean deviance, and deviance of mean
{
  # Mean log-likelihood for each subject
  ll <- get_pars(emc, stage = stage, filter = filter, selection = "LL", merge_chains = TRUE)
  minDs <- -2*apply(ll[[1]][[1]], 2, min)
  mean_lls <- apply(ll[[1]][[1]], 2, mean)
  alpha <- get_pars(emc,selection="alpha",stage=stage,filter=filter, by_subject = TRUE, merge_chains = TRUE)
  mean_pars <- lapply(alpha,function(x){apply(do.call(rbind,x),2,mean)})
  # log-likelihood for each subject using their mean parameter vector
  ll_func <- attr(emc,"design_list")[[1]]$model()$log_likelihood
  data <- emc[[1]]$data
  mean_pars_lls <- setNames(numeric(length(mean_pars)),names(mean_pars))
  for (sub in names(mean_pars)){
    if(!is.data.frame(emc[[1]]$data[[1]])){
      mean_pars_lls[sub] <- log_likelihood_joint(t(mean_pars[[sub]]),dadms = data[[sub]])
    } else{
      mean_pars_lls[sub] <- ll_func(mean_pars[[sub]],dadm = data[[sub]])
    }
  }
  Dmeans <- -2*mean_pars_lls

  if (!is.null(subject)) {
    Dmeans <- Dmeans[subject[1]]
    mean_lls <- mean_lls[subject[1]]
    minDs <- minDs[subject[1]]
  } else{
    group_stats <- attr(emc[[1]], "variant_funs")$group_IC(emc, stage=stage,filter=filter)
    if(group_only){
      mean_lls <- group_stats$mean_ll
      minDs <- group_stats$minD
      Dmeans <- group_stats$Dmean
    } else{
      mean_lls <- c(mean_lls, group_stats$mean_ll)
      minDs <- c(minDs, group_stats$minD)
      Dmeans <- c(Dmeans, group_stats$Dmean)
    }
  }
  if (use_best_fit) minDs <- pmin(minDs,Dmeans)

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
#' @param MLL1 Numeric. Marginal likelihood of model 1. Obtained with `run_bridge_sampling()`
#' @param MLL2 Numeric. Marginal likelihood of model 2. Obtained with `run_bridge_sampling()`
#'
#' @return The BayesFactor for model 1 over model 2
#' @examples \dontrun{
#' # First get the marginal likelihood for two_models
#' # Here the full model is an emc object with the hypothesized effect
#' # The null model is an emc object without the hypothesized effect
#' MLL_full <- run_bridge_sampling(full_model, cores_per_prop = 3)
#' MLL_null <- run_bridge_sampling(null_model, cores_per_prop = 3)
#' # Now we can calculate their Bayes factor
#' get_BayesFactor(MLL_full, MLL_null)
#' }
#' @export
get_BayesFactor <- function(MLL1, MLL2){
  exp(MLL1 - MLL2)
}

#' Information criteria for each participant
#'
#' Returns the BPIC/DIC based model weights for each participant in a list of samples objects
#'
#' @param sList List of samples objects
#' @param stage A string. Specifies which stage the samples are to be taken from `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`
#' @param filter An integer or vector. If it's an integer, iterations up until the value set by `filter` will be excluded.
#' If a vector is supplied, only the iterations in the vector will be considered.
#' @param use_best_fit Boolean, defaults to `TRUE`, use minimal likelihood or mean likelihood
#' (whichever is better) in the calculation, otherwise always uses the mean likelihood.
#' @param print_summary Boolean (defaults to `TRUE`) print table of results
#' @param digits Integer, significant digits in printed table
#'
#' @return List of matrices for each subject of effective number of parameters,
#' mean deviance, deviance of mean, DIC, BPIC and associated weights.
#' @examples \dontrun{
#' # Define a list of two (or more different models)
#' # Here the full model is an emc object with the hypothesized effect
#' # The null model is an emc object without the hypothesized effect
#' design_full <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Now without a ~ E
#' design_null <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#'
#' full_model <- make_emc(forstmann, design_full)
#' full_model <- fit(full_model, cores_for_chains = 1)
#'
#' null_model <- make_emc(forstmann, design_null, cores_for_chains = 1)
#' null_model <- fit(null_model)
#' sList <- list(full_model, null_model)
#' compare_subject(sList)
#' # prints a set of weights for each model for the different participants
#' # And returns the DIC and BPIC for each participant for each model.
#' }
#' @export
compare_subject <- function(sList,stage="sample",filter=0,use_best_fit=TRUE,
                            print_summary=TRUE,digits=3) {
  if(is(sList, "emc")) sList <- list(sList)
  subjects <- names(sList[[1]][[1]]$data)
  out <- setNames(vector(mode="list",length=length(subjects)),subjects)
  for (i in subjects) out[[i]] <- compare(sList,subject=i,BayesFactor=FALSE,
                                          stage=stage,filter=filter,use_best_fit=use_best_fit,print_summary=FALSE)
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


make_nice_summary <- function(object, stat = "max", stat_only = FALSE, stat_name = NULL, ...){
  if(is.null(stat_name)) stat_name <- stat
  row_names <- names(object)
  col_names <- unique(unlist(lapply(object, names)))
  if(all(row_names %in% col_names)){
    col_names <- row_names
  }
  out_mat <- matrix(NA, nrow = length(row_names), ncol = length(col_names))
  for(i in 1:length(object)){
    idx <- col_names %in% names(object[[i]])
    out_mat[i,idx] <- object[[i]]
  }
  row_stat <- apply(out_mat, 1, FUN = get(stat), na.rm = T)
  out_mat <- cbind(out_mat, row_stat)

  if(nrow(out_mat) > 1){
    col_stat <- apply(out_mat, 2, FUN = get(stat), na.rm = T)
    col_stat[length(col_stat)] <- get(stat)(unlist(object))
    out_mat <- rbind(out_mat, c(col_stat))
    rownames(out_mat) <- c(row_names, stat_name)
  } else{
    rownames(out_mat) <- row_names
  }
  colnames(out_mat) <- c(col_names, stat_name)
  if(stat_only){
    out_mat <- out_mat[nrow(out_mat), ncol(out_mat)]
  }
  return(out_mat)
}


get_summary_stat <- function(emc, selection = "mu", fun, stat = NULL,
                             stat_only = FALSE, stat_name = NULL, digits = 3, ...){
  dots <- list(...)
  if(length(dots$subject) == 1 || emc[[1]]$n_subjects == 1) dots$by_subject <- TRUE
  MCMC_samples <- do.call(get_pars, c(list(emc = emc, selection = selection), fix_dots(dots, get_pars)))
  out <- vector("list", length = length(MCMC_samples))
  for(i in 1:length(MCMC_samples)){
    # cat("\n", names(MCMC_samples)[[i]], "\n")
    if(length(fun) > 1){
      outputs <- list()
      for(j in 1:length(fun)){
        outputs[[j]] <- do.call(fun[[j]], c(list(MCMC_samples[[i]]), fix_dots(dots, fun[[j]])))
      }
      out[[i]] <- do.call(cbind, outputs)
      if(!is.null(stat_name)){
        if(ncol(out[[i]]) != length(stat_name)) stop("make sure stat_name is the same length as function output")
        colnames(out[[i]]) <- stat_name
      }
    } else{
      out[[i]] <- do.call(fun, c(list(MCMC_samples[[i]]), fix_dots(dots, fun)))#fun(MCMC_samples[[i]], ...)
    }
  }
  names(out) <- names(MCMC_samples)
  if(length(fun) == 1 & !is.matrix(out[[i]]) & !is.null(stat)){
    out <- make_nice_summary(out, stat, stat_only, stat_name)
    out <- round(out, digits)
  } else{
    out <- lapply(out, round, digits)
  }
  return(out)
}


get_posterior_quantiles <- function(x, probs = c(0.025, .5, .975)){
  summ <- summary(x, probs)
  return(summ$quantiles)
}



