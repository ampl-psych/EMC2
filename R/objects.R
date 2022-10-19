thin_pmwg <- function(pmwg,thin=10)
  # removes samples from single pmwg object
{
  # mcmc creation functions allow thinning, but for cases where you want to
  # reduce the size, this will thin the sample component of a pmwg object
  if (class(pmwg) != "pmwgs") stop("Not a pmwgs object")
  if (thin >= pmwg$samples$idx) stop("Thin value to large\n")
  nmc_thin <- seq(thin,pmwg$samples$idx,by=thin)
  pmwg$samples <- base::rapply(pmwg$samples, f = function(x) filter_obj(x, nmc_thin), how = "replace")
  pmwg$samples$stage <- pmwg$samples$stage[nmc_thin]
  pmwg$samples$idx <- length(nmc_thin)
  return(pmwg)
}


filter_obj <- function(obj, idx){
  dims <- dim(obj)
  dim_names <- dimnames(obj)
  if(is.null(dims)) return(obj)
  if(length(dims) == 2){
    if(isSymmetric(round(obj, 3))) return(obj) #Don't extend priors and theta_mu_var_inv
  }
  obj <- obj[slice.index(obj, length(dims)) %in% idx]
  dims[length(dims)] <- length(idx)
  dim(obj) <- dims
  dimnames(obj) <- dim_names # Give back to the community
  return(obj)
}

get_stage <- function(pmwg,stage="burn")
  # Returns object with only stage samples
{
  if (class(pmwg) != "pmwgs") stop("Not a pmwgs object")
  if (!(stage %in% pmwg$samples$stage))
    stop("stage not available")
  idx <- which(pmwg$samples$stage == stage)
  pmwg$samples <- base::rapply(pmwg$samples, f = function(x) filter_obj(x, idx), how = "replace")
  pmwg$samples$stage <- pmwg$samples$stage[idx]
  pmwg$samples$idx <- length(idx)
  pmwg
}


remove_iterations <- function(pmwg,select,remove=TRUE,last_select=FALSE,filter=NULL)
  # Removes samples up to scalar remove or in remove vector
  # if remove=FALSE instead keeps these iterations
  # if last_select applies scalar select from last iteration
{
  if (class(pmwg) != "pmwgs") stop("Not a pmwgs object")
  if (length(select)==1) {
    if (is.null(filter)) {
      if (!last_select)
        select <- 1:select else
          select <- (pmwg$samples$idx-select+1):pmwg$samples$idx
    } else {
      n <- table(pmwg$samples$stage)
      if (filter=="burn") {
        if (select>n["burn"]) stop("Removing more than available in burn")
        if (!last_select)
          select <- 1:select else
            select <- (n["burn"]-select+1):n["burn"]
      } else if (filter=="adapt") {
        if (select>n["adapt"]) stop("Removing more than available in adapt")
        if (!last_select)
          select <- 1:select else
            select <- (n["adapt"]-select+1):n["adapt"]
          select <- select + n["burn"]
      } else if (filter=="sample") {
        if (select>n["sample"]) stop("Removing more than available in sample")
        if (!last_select)
          select <- 1:select else
            select <- (n["sample"]-select+1):n["sample"]
          select <- select + n["burn"] + n["adapt"]
      }
    }
  }
  if (any(select>pmwg$samples$idx))
    stop("select specifies iterations not present")
  ok <- 1:pmwg$samples$idx %in% select
  if (remove) ok <- !ok
  filter_idx <- which(ok)
  pmwg$samples <- base::rapply(pmwg$samples, f = function(x) filter_obj(x, filter_idx), how = "replace")
  pmwg$samples$stage <- pmwg$samples$stage[ok]
  pmwg$samples$idx <- sum(ok)
  pmwg
}

merge_samples <- function(samples){
  out_samples <- samples[[1]]
  # Only thing that differs between the chains is the samples$samples
  sampled_objects <- lapply(samples, FUN = function(x) return(x$samples))
  keys <- unique(unlist(lapply(sampled_objects, names)))
  sampled_objects <- setNames(do.call(mapply, c(abind, lapply(sampled_objects, '[', keys))), keys)
  sampled_objects$idx <- sum(sampled_objects$idx)
  out_samples$samples <- sampled_objects
  return(out_samples)
}


#### mcmc extraction ----

as_Mcmc <- function(sampler,filter=stages,thin=1,subfilter=0,
                    selection=c("alpha","mu","variance","covariance","correlation","LL","epsilon")[1])
  # replacement for pmwg package as_mccmc
  # allows for selection of subj_ll, and specifying selection as an integer
  # allows filter as single integer so can remove first filter samples
  # adds ability to thin and to subfilter (i.e., remove trials after removing
  # a stage), range of values or single integer (expanded to end)
{

  shorten <- function(x,r,d) {
    if (!is.numeric(r)) stop("subfilter must be numeric\n")
    if (length(r)==1) r <- (r+1):dim(x)[d]
    if (min(r)<1 || max(r)>dim(x)[d]) stop("subfilter out of range\n")
    if (d==2) x[,r,drop=FALSE] else x[,,r,drop=FALSE]
  }

  if (class(sampler) != "pmwgs") stop("sampler must be an pmwgs object\n")
  stages <- c("burn", "adapt", "sample")
  if (is.numeric(filter) && length(filter==1)) filter <- filter:sampler$samples$idx
  if (is.numeric(selection))
    selection <- selection=c("alpha","mu","variance","covariance","correlation","LL")[selection]
  if (selection %in% c("mu","variance","covariance","correlation") &&
      all(is.na((sampler$samples$theta_mu))))
    stop("Cannot select mu, variance, covariance or correlation for non-hierarchical model\n")
  if (all(filter %in% stages)) {
    filter <- which(sampler$samples$stage %in% filter)
  } else if (!all(filter %in% 1:sampler$samples$idx)) {
    stop("filter is not a vector of stage names, or integer vector of indices\n")
  }

  if (selection == "mu") {
    mu <- sampler$samples$theta_mu[, filter,drop=FALSE]
    if (is.null(subfilter)) mu <- t(mu) else
      mu <- t(shorten(mu,subfilter,2))
    if (thin > dim(mu)[1]) stop("Thin to large\n")
    out <- coda::mcmc(mu[seq(thin,dim(mu)[1],by=thin),,drop=FALSE])
    attr(out,"selection") <- selection
    return(out)
  } else if (selection %in% c("variance","covariance","correlation")) {
    tvar <- sampler$samples$theta_var[, , filter,drop=FALSE]
    if (selection=="correlation")
      tvar <- array(apply(tvar,3,stats::cov2cor),dim=dim(tvar),dimnames=dimnames(tvar))
    if (selection %in% c("covariance","correlation")) {
      lt <- lower.tri(tvar[,,1])
      pnams <- dimnames(tvar)[[1]]
      pnams <- outer(pnams,pnams,paste,sep=".")[lt]
      tvar <- apply(tvar,3,function(x){x[lt]})
      dimnames(tvar)[[1]] <- pnams
      if (all(tvar==0))
        stop("All covariances/correlations zero, model assumes independence.")
    } else tvar <- apply(tvar,3,diag)
    if (is.null(subfilter)) tvar <- t(tvar) else
      tvar <- t(shorten(tvar,subfilter,2))
    if (thin > dim(tvar)[1]) stop("Thin to large\n")
    out <- coda::mcmc(tvar[seq(thin,dim(tvar)[1],by=thin),,drop=FALSE])
    attr(out,"selection") <- selection
    return(out)
  } else if (selection == "alpha") {
    alpha <- sampler$samples$alpha[, , filter,drop=FALSE]
    if (!is.null(subfilter)) alpha <- shorten(alpha,subfilter,3)
    if (thin > dim(alpha)[3]) stop("Thin to large\n")
    alpha <- alpha[,,seq(thin,dim(alpha)[3],by=thin),drop=FALSE]
    out <- stats::setNames(lapply(
      seq(dim(alpha)[2]),
      function(x) {
        coda::mcmc(t(alpha[, x, ]))
      }
    ), names(sampler$data))
    attr(out,"selection") <- selection
    return(out)
  } else if (selection %in% c("LL","epsilon")) {
    if (selection == "epsilon")
      LL <- sampler$samples$epsilon[, filter,drop=FALSE] else
        LL <- sampler$samples$subj_ll[, filter,drop=FALSE]
      if (!is.null(subfilter)) LL <- shorten(LL,subfilter,2)
      if (thin > dim(LL)[2]) stop("Thin to large\n")
      out <- setNames(lapply(
        data.frame(t(LL[,seq(thin,dim(LL)[2],by=thin),drop=FALSE])),coda::mcmc),
        names(sampler$data))
      attr(out,"selection") <- selection
      return(out)
  }
  stop("Argument `selection` should be one of mu, variance, covariance, correlation, alpha, LL, or epsilon\n")
}




# samplers=pmwg_mcmc
as_mcmc.list <- function(samplers,
                         selection=c("alpha","mu","variance","covariance","correlation","LL","epsilon")[1],
                         filter="burn",thin=1,subfilter=0,mapped=FALSE,include_constants=FALSE)
  # Combines as_Mcmc of samplers and returns mcmc.list
  # mapped = TRUE map mu or alpha to model parameters (constants indicated by
  # attribute isConstant)
  # include_constants= TRUE, if not mapped add in constants
{

  stages <- c("burn", "adapt", "sample")
  if (length(filter)>1) filter <- filter[1]
  if (!class(samplers)=="list" || !all(unlist(lapply(samplers,class))=="pmwgs"))
    stop("samplers must be a list of pmwgs objects")
  if (!all(apply(do.call(rbind,lapply(samplers,function(x){
    x$par_names})),2,function(x){all(x[1]==x[-1])})))
    stop("samplers must have the same parameter names")
  if (!all(apply(do.call(rbind,lapply(samplers,function(x){
    x$subjects})),2,function(x){all(x[1]==x[-1])})))
    stop("samplers must have the same subjects")

  subjects <- names(samplers[[1]]$data)
  mcmcList <- lapply(samplers,as_Mcmc,selection=selection,filter=filter,
                     thin=thin,subfilter=subfilter)
  if (mapped) {
    if (selection == "alpha") {
      mcmcList <- lapply(mcmcList,function(x){lapply(x,map_mcmc, include_constants = include_constants,
                                                     design=attr(samplers,"design_list")[[1]],model=attr(samplers,"model_list")[[1]])})
    } else if (selection=="mu") {
      mcmcList <- lapply(mcmcList,map_mcmc, include_constants = include_constants,
                         design=attr(samplers,"design_list")[[1]],model=attr(samplers,"model_list")[[1]])
    } else warning("Can only map alpha or mu to model parameterization")
  }
  if (include_constants) {
    if (selection == "alpha") {
      mcmcList <- lapply(mcmcList,function(x){lapply(x,add_constants_mcmc,
                                                     constants=attr(samplers,"design_list")[[1]]$constants)})
    } else if (selection=="mu") {
      mcmcList <- lapply(mcmcList,add_constants_mcmc,
                         constants=attr(samplers,"design_list")[[1]]$constants)
    }
  }
  if (selection %in% c("LL","epsilon"))
    iter <- unlist(lapply(mcmcList,function(x){length(x[[1]])})) else
      if (selection == "alpha")
        iter <- unlist(lapply(mcmcList,function(x){dim(x[[1]])[1]})) else
          iter <- unlist(lapply(mcmcList,function(x){dim(x)[1]}))
  if (!all(iter[1]==iter[-1]))
    message("Chains have different numbers of samples, using first ",min(iter))
  iter <- min(iter)
  if (selection %in% c("alpha","LL","epsilon")) {
    nChains <- length(mcmcList)
    ns <- length(samplers[[1]]$subjects)
    out <- vector(mode="list",length=ns)
    lst <- vector(mode="list",length=nChains)
    for (i in 1:ns) {
      for (j in 1:nChains) if (selection %in% c("LL","epsilon"))
        lst[[j]] <- as.mcmc(mcmcList[[j]][[i]][1:iter]) else {
          isConstant <- attr(mcmcList[[j]][[i]],"isConstant")
          lst[[j]] <- as.mcmc(mcmcList[[j]][[i]][1:iter,])
          attr(lst[[j]],"isConstant") <- isConstant
        }
      out[[i]] <- mcmc.list(lst)
    }
    out <- setNames(out,subjects)
    attr(out,"selection") <- selection
    return(out)
  } else {
    out <- mcmc.list(mcmcList)
    attr(out,"selection") <- selection
    return(out)
  }
}




