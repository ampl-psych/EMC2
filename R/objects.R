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
