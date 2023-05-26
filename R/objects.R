thin_pmwg <- function(pmwg,thin=10)
  # removes samples from single pmwg object
{
  # mcmc creation functions allow thinning, but for cases where you want to
  # reduce the size, this will thin the sample component of a pmwg object
  if (!is(pmwg, "pmwgs")) stop("Not a pmwgs object")
  if (thin >= pmwg$samples$idx) stop("Thin value to large\n")
  nmc_thin <- seq(thin,pmwg$samples$idx,by=thin)
  if(any(pmwg$nuisance)) pmwg$sampler_nuis$samples <- base::rapply(pmwg$sampler_nuis$samples, f = function(x) filter_obj(x, nmc_thin), how = "replace")
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
    if(nrow(obj) == ncol(obj)){
      if(nrow(obj) > 1){
        if(mean(abs(abs(rowSums(obj/max(obj))) - abs(colSums(obj/max(obj))))) < .01) return(obj)
      }
    }
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
  if (!is(pmwg, "pmwgs")) stop("Not a pmwgs object")
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
  if (!is(pmwg, "pmwgs")) stop("Not a pmwgs object")
  if (length(select)==1) {
    if (is.null(filter)) {
      if (!last_select)
        select <- 1:select else
          select <- (pmwg$samples$idx-select+1):pmwg$samples$idx
    } else {
      n <- table(pmwg$samples$stage)
      if (filter == "preburn"){
        if (select>n["preburn"]) stop("Removing more than available in preburn")
        if (!last_select)
          select <- 1:select else
            select <- (n["preburn"]-select+1):n["preburn"]
      }
      else if (filter=="burn") {
        if (select>n["burn"]) stop("Removing more than available in burn")
        if (!last_select)
          select <- 1:select else
            select <- (n["burn"]-select+1):n["burn"]
          select <- select +  n["preburn"]
      } else if (filter=="adapt") {
        if (select>n["adapt"]) stop("Removing more than available in adapt")
        if (!last_select)
          select <- 1:select else
            select <- (n["adapt"]-select+1):n["adapt"]
          select <- select + n["preburn"] + n["burn"]
      } else if (filter=="sample") {
        if (select>n["sample"]) stop("Removing more than available in sample")
        if (!last_select)
          select <- 1:select else
            select <- (n["sample"]-select+1):n["sample"]
          select <- select + n["preburn"] + n["burn"] + n["adapt"]
      }
    }
  }
  if (any(select>pmwg$samples$idx))
    stop("select specifies iterations not present")
  ok <- 1:pmwg$samples$idx %in% select
  if (remove) ok <- !ok
  filter_idx <- which(ok)
  if(any(pmwg$nuisance)) {
    pmwg$sampler_nuis$samples <- base::rapply(pmwg$sampler_nuis$samples, f = function(x) filter_obj(x, filter_idx), how = "replace")
    pmwg$sampler_nuis$samples$idx <- sum(ok)
  }
  pmwg$samples <- base::rapply(pmwg$samples, f = function(x) filter_obj(x, filter_idx), how = "replace")
  pmwg$samples$stage <- pmwg$samples$stage[ok]
  pmwg$samples$idx <- sum(ok)
  pmwg
}


#' Merge samples
#'
#' Merges samples from all chains together as one unlisted object
#'
#' @param samplers A list of samplers
#'
#' @return an unlisted sampler with all chains merged
#' @export
#'
merge_samples <- function(samplers){
  out_samples <- samplers[[1]]
  # Only thing that differs between the chains is the samples$samples
  sampled_objects <- lapply(samplers, FUN = function(x) return(x$samples))
  keys <- unique(unlist(lapply(sampled_objects, names)))
  sampled_objects <- setNames(do.call(mapply, c(abind, lapply(sampled_objects, '[', keys))), keys)
  sampled_objects$idx <- sum(sampled_objects$idx)
  out_samples$samples <- sampled_objects
  if(any(out_samples$nuisance)){
    sampled_objects <- lapply(samplers, FUN = function(x) return(x$sampler_nuis$samples))
    keys <- unique(unlist(lapply(sampled_objects, names)))
    sampled_objects <- setNames(do.call(mapply, c(abind, lapply(sampled_objects, '[', keys))), keys)
    out_samples$sampler_nuis$samples <- sampled_objects
  }
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

  if (!is(sampler, "pmwgs")) stop("sampler must be an pmwgs object\n")
  stages <- c("preburn", "burn", "adapt", "sample")
  if (is.numeric(filter) && length(filter==1)) filter <- filter:sampler$samples$idx
  if (is.numeric(selection)) selection <- c("alpha","mu","variance","covariance","correlation","LL")[selection]
  if (selection %in% c("mu","variance","covariance","correlation") && all(is.na((sampler$samples$theta_mu)))){
    stop("Cannot select mu, variance, covariance or correlation for non-hierarchical model\n")
  }
  if (all(filter %in% stages)) {
    filter <- which(sampler$samples$stage %in% filter)
  } else if (!all(filter %in% 1:sampler$samples$idx)) {
    stop("filter is not a vector of stage names, or integer vector of indices\n")
  }

  if (selection == "mu") {
    mu <- sampler$samples$theta_mu[, filter,drop=FALSE]
    if (is.null(subfilter)){
      mu <- t(mu)
    }  else{
      mu <- t(shorten(mu,subfilter,2))
    }
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
    } else {
      tvar <- apply(tvar,3,diag)
    }
    if (is.null(subfilter)){
      tvar <- t(tvar)
    } else{
      tvar <- t(shorten(tvar,subfilter,2))
    }
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
    if (selection == "epsilon"){
      LL <- sampler$samples$epsilon[, filter,drop=FALSE]
    } else{
      LL <- sampler$samples$subj_ll[, filter,drop=FALSE]
    }
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




# samplers=pmwg_mcmc; mapped=FALSE;include_constants=FALSE
as_mcmc.list <- function(samplers,
                         selection=c("alpha","mu","variance","covariance","correlation","LL","epsilon")[1],
                         filter="burn",thin=1,subfilter=0,mapped=FALSE,include_constants=FALSE)
  # Combines as_Mcmc of samplers and returns mcmc.list
  # mapped = TRUE map mu or alpha to model parameters (constants indicated by
  # attribute isConstant)
  # include_constants= TRUE, if not mapped add in constants
{

  stages <- c("preburn", "burn", "adapt", "sample")
  if(!all(filter %in% stages)) stop(paste0("must select filter from stages: ", stages))

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
    } else warning("Can only include constants in alpha or mu")
  }

  if (selection %in% c("alpha","LL")) {
    nChains <- length(mcmcList)
    ns <- length(samplers[[1]]$subjects)
    out <- vector(mode="list",length=ns)
    lst <- vector(mode="list",length=nChains)
    for (i in 1:ns) {
      for (j in 1:nChains){
          isConstant <- attr(mcmcList[[j]][[i]],"isConstant")
          lst[[j]] <- as.mcmc(mcmcList[[j]][[i]])
          attr(lst[[j]],"isConstant") <- isConstant
      }
      out[[i]] <- mcmc.list(lst)
    }
    out <- setNames(out,subjects)

  } else {
    out <- mcmc.list(mcmcList)
  }
  attr(out,"selection") <- selection
  return(out)
}

#' Returns the number of samples per chain per stage
#'
#' @param samplers A list of samplers, could be in any stage
#'
#' @return A table of iterations per stage per chain
#' @export

chain_n <- function(samplers)
{
  do.call(rbind,lapply(samplers, function(x){
    table(factor(x$samples$stage,levels=c("preburn", "burn","adapt","sample")))
  }))
}

extract_samples <- function(sampler, stage = c("adapt", "sample"), max_n_sample = NULL) {
  variant_funs <- attr(sampler, "variant_funs")
  samples <- sampler$samples

  if("sample" %in% stage & !is.null(max_n_sample)){
    sample_filter <- which(samples$stage %in% "sample" & seq_along(samples$stage) <= samples$idx)
    adapt_filter <- which(samples$stage %in% "adapt" & seq_along(samples$stage) <= samples$idx)
    if(length(sample_filter) > max_n_sample){
      sample_filter <- sample(sample_filter, max_n_sample)
    }
    full_filter <- c(adapt_filter, sample_filter)
  } else{
    full_filter <- which(samples$stage %in% stage & seq_along(samples$stage) <= samples$idx)
  }
  if(any(sampler$nuisance)){
    nuisance <- sampler$nuisance[!sampler$grouped]
    sampler$sampler_nuis$samples$alpha <- sampler$samples$alpha[nuisance,,]
    sampler$samples$alpha <- sampler$samples$alpha[!nuisance,,]
    out <- variant_funs$filtered_samples(sampler, full_filter)
    out$nuisance <- get_variant_funs("diagonal")$filtered_samples(sampler$sampler_nuis, full_filter)
  } else{
    out <- variant_funs$filtered_samples(sampler, full_filter)
  }
  return(out)
}

p_names <- function(samples,mapped=FALSE,design=FALSE)
  # parameter names of a pmwg object or list of pmwg objects or if mapped=TRUE
  # gets a list of names for mapped parameters of each type
{

  sp <- mp <- mapped_name_list(attr(samples,"design_list")[[1]],
                               attr(samples,"model_list")[[1]],design)
  if (mapped) return(mp)
  if (inherits(samples,"pmwg"))
    tmp <- dimnames(samples$samples$alpha)[[1]] else
      tmp <- dimnames(samples[[1]]$samples$alpha)[[1]]
    type <- unlist(lapply(strsplit(tmp,"_"),function(x){x[[1]]}))
    for (i in names(sp)) sp[[i]] <- tmp[type==i]
    sp
}

parameters_data_frame <- function(samples,filter="sample",thin=1,subfilter=0,
                                  mapped=FALSE,include_constants=FALSE,stat=NULL,
                                  selection=c("alpha","mu","variance","covariance","correlation","LL","epsilon")[2])
  # extracts and stacks chains into a matrix
{
  out <- as_mcmc.list(samples,selection=selection,filter=filter,thin=thin,
                      subfilter=subfilter,mapped=mapped,include_constants=include_constants)
  if (selection!="alpha") out <- as.data.frame(do.call(rbind,out)) else {
    out <- lapply(out,function(x){do.call(rbind,x)})
    subjects <- rep(names(out),lapply(out,function(x){dim(x)[1]}))
    out <- cbind.data.frame(subjects=factor(subjects,levels=names(out)),do.call(rbind,out))
  }
  if (!is.null(stat)) {
    if (selection=="alpha") {
      ss <- levels(out$subjects)
      outList <- setNames(vector(mode="list",length=length(ss)),ss)
      for (s in ss) outList[[s]] <- apply(out[out$subjects==s,-1],2,stat)
      out <- do.call(rbind,outList)
    } else out <- apply(out,2,stat)
  }
  out
}

get_sigma <- function(samps,filter="sample",thin=1,subfilter=0)
  # Get sigma matrix samples
{
  shorten <- function(x,r,d) {
    if (!is.numeric(r)) stop("subfilter must be numeric\n")
    if (length(r)==1) r <- (r+1):dim(x)[d]
    if (min(r)<1 || max(r)>dim(x)[d]) stop("subfilter out of range\n")
    if (d==2) x[,r,drop=FALSE] else x[,,r,drop=FALSE]
  }

  sigma <- samps$samples$theta_var[,,samps$samples$stage==filter]
  if (!is.null(subfilter)) sigma <- shorten(sigma,subfilter,3)
  if (thin >= samps$samples$idx) stop("Thin value to large\n")
  nmc_thin <- seq(thin,dim(sigma)[3],by=thin)
  sigma[,,nmc_thin]
}

subject_names <- function(samples)
  # Get names of subjects
{
  if (inherits(samples, "pmwg"))
    names(samples$data) else
      names(samples[[1]]$data)
}


#' Title
#'
#' @param samplers A list of samplers.
#' @param thin An integer. Indicates by how much to thin.
#'
#' @return A thinned list of samplers
#' @export

thin_chains <- function(samplers,thin=5)
  # thins a set of chains
{
  data_list <- attr(samplers,"data_list")
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  samplers <- lapply(samplers,thin_pmwg,thin=thin)
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  samplers
}
