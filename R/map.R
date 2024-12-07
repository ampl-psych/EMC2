do_transform <- function(pars,transform)
  # pars is the parameter matrix, transform is a transform list
{
  isexp <- transform$func[colnames(pars)] == "exp"
  isprobit <- transform$func[colnames(pars)] == "pnorm"
  # Here instead of using isexp directly I use the column names in case v is replicated in pars (i.e. in map_mcmc)
  pars[,isexp] <- exp(sweep(pars[,isexp, drop = F], 2, transform$lower[colnames(pars)[isexp]], "-"))
  pars[,isprobit] <- pnorm(sweep(sweep(pars[,isprobit, drop = F], 2, transform$lower[colnames(pars)[isprobit]], "-")
                                 , 2, transform$upper[colnames(pars)[isprobit]]-transform$lower[colnames(pars)[isprobit]], "/"))
  pars
}

do_pre_transform <- function(p_vector,transform)
  # pars is the parameter matrix, transform is a transform list
{
  isexp <- transform$func[names(p_vector)] == "exp"
  isprobit <- transform$func[names(p_vector)] == "pnorm"
  # Here instead of using isexp directly I use the names in case v is replicated in p_vector (i.e. in map_mcmc)
  p_vector[isexp] <- exp(p_vector[isexp] - transform$lower[names(p_vector)[isexp]])
  p_vector[isprobit] <- pnorm((p_vector[isprobit] - transform$lower[names(p_vector)[isprobit]])/
                                (transform$upper[names(p_vector)[isprobit]]-transform$lower[names(p_vector)[isprobit]]))
  p_vector
}




get_pars_matrix <- function(p_vector,dadm) {
  # Order:
  # 1 pretransform
  # 2 add constants
  # 3 map
  # # - if premap trend:
  # #   First make trend pars matrix
  # #   Apply trends premap and remove trend pars from pars matrix
  # # - map
  # 4 if pretransform trend:
  # # - apply trends and remove trend pars from pars matrix
  # 5 transform
  # 6 if posttransform trend:
  # # - apply trends and remove trend pars from pars matrix
  # 7 trial-wise transform
  # 8 bound

  # Niek should constants be included in pre_transform? I think not?
  p_vector <- do_pre_transform(p_vector, attr(dadm, "model")()$pre_transform)
  # If there's any premap trends, they're done in map_p
  pars <- map_p(add_constants(p_vector,attr(dadm,"constants")),dadm)
  if(!is.null(attr(dadm, "model")()$trend) && attr(attr(dadm, "model")()$trend, "pretransform")){
    # This runs the trend and afterwards removes the trend parameters
    pars <- prep_trend(dadm, attr(dadm, "model")()$trend, pars)
  }
  pars <- do_transform(pars, attr(dadm,"model")()$transform)
  if(!is.null(attr(dadm, "model")()$trend) && attr(attr(dadm, "model")()$trend, "posttransform")){
    # This runs the trend and afterwards removes the trend parameters
    pars <- prep_trend(dadm, attr(dadm, "model")()$trend, pars)
  }
  pars <- attr(dadm,"model")()$Ttransform(pars, dadm)
  pars <- add_bound(pars, attr(dadm,"model")()$bound)
  return(pars)
}

make_pmat <- function(p_vector,design)
  # puts vector form of p_vector into matrix form
{
  ss <- design$Ffactors$subjects
  matrix(rep(p_vector,each=length(ss)),nrow=length(ss),
         dimnames=list(ss,names(p_vector)))
}

add_constants <- function(p,constants)
{
  if (is.null(constants)) return(p)
  if (is.matrix(p)) {
    nams <- c(dimnames(p)[[2]],names(constants))
    p <- cbind(p,matrix(rep(constants,each=dim(p)[1]),nrow=dim(p)[1]))
    dimnames(p)[[2]] <- nams
    return(p)
  } else{
    return(c(p,constants))
  }

}

#' Parameter mapping back to the design factors
#'
#' Maps a parameter vector that corresponds to sampled parameters
#' of the cognitive model back to the experimental design. The parameter vector
#' can be created using ``sampled_p_vector()``. The returned matrix shows whether/how parameters
#' differ across the experimental factors.
#'
#' @param p_vector A parameter vector. Must be in the form of ``sampled_p_vector(design)``
#' @param design A design list. Created by ``design``
#' @param model Optional model type (if not already specified in ``design``)
#' @param digits Integer. Will round the output parameter values to this many decimals
#' @param ... optional arguments
#' @param remove_subjects Boolean. Whether to include subjects as a factor in the design
#' @param covariates Covariates specified in the design can be included here.
#' @return Matrix with a column for each factor in the design and for each model parameter type (``p_type``).
#' @examples
#' # First define a design:
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then create a p_vector:
#' p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'           t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' # This will map the parameters of the p_vector back to the design
#' mapped_par(p_vector,design_DDMaE)
#'
#' @export

mapped_par <- function(p_vector,design,model=NULL,
                       digits=3,remove_subjects=TRUE,
                       covariates=NULL,...)
  # Show augmented data and corresponding mapped parameter
{
  remove_RACE <- TRUE
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if (is.null(covariates))
      Fcovariates <- design$Fcovariates else
      Fcovariates <- covariates
  if (is.null(model)) if (is.null(design$model))
  stop("Must specify model as not in design") else model <- design$model
  if (remove_subjects) design$Ffactors$subjects <- design$Ffactors$subjects[1]
  if (!is.matrix(p_vector)) p_vector <- make_pmat(p_vector,design)
  dadm <- design_model(make_data(p_vector,design,n_trials=1,Fcovariates=Fcovariates),
                       design,model,rt_check=FALSE,compress=FALSE, verbose = FALSE)
  ok <- !(names(dadm) %in% c("subjects","trials","R","rt","winner"))
  out <- cbind(dadm[,ok],round(get_pars_matrix(p_vector,dadm),digits))
  if (model()$type=="SDT")  out <- out[dadm$lR!=levels(dadm$lR)[length(levels(dadm$lR))],]
  if (model()$type=="DDM")  out <- out[,!(names(out) %in% c("lR","lM"))]
  if (any(names(out)=="RACE") && remove_RACE)
    out <- out[as.numeric(out$lR) <= as.numeric(as.character(out$RACE)),,drop=FALSE]
  return(out)
}

add_recalculated_pars <- function(pmat, model, cnams){
  modifiers <- unlist(lapply(strsplit(cnams,"_"),function(x){paste0(x[-1], collapse = "_")}))
  par_names <- colnames(pmat)
  unq_pars <- unique(par_names)
  par_table <- table(par_names)
  counts <- lapply(par_table, function(x) return(1:x))
  combn <- do.call(expand.grid, counts)
  colnames(combn) <- names(par_table)

  out <- list()
  modfs <- list()
  for(r in 1:nrow(combn)){
    pmat_in <- matrix(NA, nrow = nrow(pmat), ncol = length(unq_pars))
    colnames(pmat_in) <- unq_pars
    cur_modifiers <- setNames(numeric(length(unq_pars)), unq_pars)
    for(par in unq_pars){
      pmat_in[,par] <- pmat[,which(colnames(pmat) == par)[combn[r,par]]]
      cur_modifiers[par] <- modifiers[which(colnames(pmat) == par)[combn[r,par]]]
    }
    added <- model()$Ttransform(pmat_in)
    added <- added[, !(colnames(added) %in% colnames(pmat_in)), drop = F]
    attr(added, "ok") <- NULL
    out[[r]] <- added
    modfs[[r]] <- cur_modifiers
  }
  m_out <- matrix(0, nrow = nrow(pmat), ncol = 0)
  if(ncol(out[[1]]) == 0) return(NULL)
  for(i in 1:ncol(out[[1]])){
    cur_par <- lapply(out, function(x) x[,i, drop = F])
    not_dups <- !duplicated(cur_par)
    cur_combn <- combn[not_dups,]
    pars_vary <- colnames(cur_combn)[colMeans(cur_combn) != 1]
    to_add <- do.call(cbind, cur_par[not_dups])
    cur_modfs <- modfs[not_dups]
    fnams <- sapply(cur_modfs, function(x) paste0(unique(unlist(strsplit(x[pars_vary], split = "_"))), collapse = "_"))
    colnames(to_add) <- paste0(colnames(to_add)[1], "_", fnams)
    m_out <- cbind(m_out, to_add)
  }
  return(m_out)
}

map_mcmc <- function(mcmc,design,include_constants = TRUE, add_recalculated = FALSE, covariates = NULL)
  # Maps vector or matrix (usually mcmc object) of sampled parameters to native
  # model parameterization.
{
  if(is.null(design$model)){
    design <- design[[1]]
  }
  model <- design$model

  doMap <- function(mapi,pmat, covariates = NULL){
    if(!is.null(covariates)){
      if(nrow(covariates) != nrow(pmat)) covariates <- covariates[sample(1:nrow(covariates), nrow(pmat), replace = T),, drop = F]
      pmat <- cbind(pmat, covariates)
    }
    t(mapi %*% t(pmat[,dimnames(mapi)[[2]],drop=FALSE]))
  }

  get_p_types <- function(nams, reverse = FALSE){
    if(reverse){
      out <- unlist(lapply(strsplit(nams,"_"),function(x){x[[-1]]}))
    } else{
      out <- unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
    }
    return(out)
  }


  map <- attr(sampled_p_vector(design, add_da = TRUE, all_cells_dm = TRUE),"map")

  constants <- design$constants
  if (!is.matrix(mcmc) & !is.array(mcmc)) mcmc <- t(as.matrix(mcmc))
  if(length(dim(mcmc)) == 2){
    is_matrix <- TRUE
    mcmc_array <- array(mcmc, dim = c(nrow(mcmc), 1, ncol(mcmc)))
    rownames(mcmc_array) <- rownames(mcmc)
  } else{
    mcmc_array <- mcmc
    is_matrix <- FALSE
  }
  mp <- mapped_par(mcmc_array[,1,1],design,remove_RACE=FALSE, covariates = covariates)

  for(k in 1:ncol(mcmc_array)){
    mcmc <- t(mcmc_array[,k,])
    mcmc <- t(apply(mcmc, 1, do_pre_transform, design$model()$pre_transform))
    pmat <- add_constants(mcmc,constants)
    plist <- lapply(map,doMap,pmat=pmat, covariates)
    # Give mapped variables names and flag constant
    for (i in 1:length(plist)) {
      vars <- row.names(attr(terms(design$Flist[[i]]),"factors"))
      if(any(names(design$Fcovariates) %in% vars)){
        cov_vars <- vars[(vars %in% names(design$Fcovariates))]
        vars <- vars[!(vars %in% names(design$Fcovariates))]
        has_cov <- TRUE
      } else{
        has_cov <- FALSE
      }
      uniq <- !duplicated(apply(mp[,vars, drop = F],1,paste,collapse="_"))
      if (is.null(vars)) {
        colnames(plist[[i]]) <- names(plist)[i]
      }else {
        if(length(vars) == 1) colnames(plist[[i]]) <- vars
        else if(length(vars) == 2){
          colnames(plist[[i]]) <- paste(vars[1], apply(matrix(do.call(cbind, Map(paste0, vars[-1], mp[uniq, vars[-1]])))
                                                       ,1, paste0, collapse = "_"), sep = "_")
        } else{
          colnames(plist[[i]]) <- paste(vars[1], apply(do.call(cbind, Map(paste0, vars[-1], mp[uniq, vars[-1]]))
                                                       ,1, paste0, collapse = "_"), sep = "_")
        }
      }
      if(has_cov) colnames(plist[[i]]) <- paste(colnames(plist[[i]]), cov_vars, sep = "_")
      # if (is.matrix(plist[[i]])) isConstant <- c(isConstant, apply(plist[[i]],2,function(x){all(x[1]==x[-1])}))
    }
    pmat <- do.call(cbind,plist)
    cnams <- colnames(pmat)
    colnames(pmat) <- get_p_types(cnams)
    pmat[,] <- do_transform(pmat, model()$transform)[,1:length(cnams)]

    if(add_recalculated) extras <- add_recalculated_pars(pmat, model, cnams)
    colnames(pmat) <- cnams
    if(add_recalculated) pmat <- cbind(pmat, extras)
    is_constant <- apply(pmat,2,function(x){all(duplicated(x)[-1])})
    if(!include_constants){
      pmat <- pmat[,!is_constant, drop = F]
    }
    if(k == 1){
      out <- array(0, dim = c(ncol(pmat), ncol(mcmc_array), dim(mcmc_array)[3]))
      rownames(out) <- colnames(pmat)
      colnames(out) <- colnames(mcmc_array)
    }
    out[,k,] <- t(pmat)
  }
  if(is_matrix){
    out <- out[,1,]
  }
  attr(out,"isConstant") <- is_constant
  return(out)
}



fill_transform <- function(transform, model = NULL, p_vector = NULL,
                           supported=c("identity","exp","pnorm"),
                           has_lower=c("exp","pnorm"),has_upper=c("pnorm"),
                           is_pre = FALSE) {
  # Note that for filling pre_transform model is the sampled_p_vector
  if(!is.null(transform)){
    if (!all(transform$func %in% supported)) stop("Only ", paste(supported, collapse = ", "), " transforms supported")
    if(!is_pre){
      if (!all(names(transform$func) %in% names(model()$p_types)))stop("transform parameter not in the model")
    } else{
      if (!all(names(transform$func) %in% names(p_vector))) stop("pre_transform parameter not in the model")
    }
    if(!is.null(transform$lower) & !is.null(transform$upper) & is.null(transform$func)) stop("func must be provided if lower and/or upper are provided")
    if (!is.null(transform$lower)) {
      if (!all(names(transform$lower) %in% names(transform$func[transform$func %in% has_lower])))
        stop("lower can only apply to tranforms of type ",paste(has_lower,collapse=", "))
      if (!all(names(transform$upper) %in% names(transform$func[transform$func %in% has_upper])))
        stop("lower can only apply to tranforms of type ", paste(has_upper, collapse = ", "))
    }
  }
  if(!is_pre){
    filled_transform <- model()$transform$func
  } else{
    if(!is.null(model()$pre_transform$func)){
      filled_transform <- model()$pre_transform$func
      filled_transform[names(p_vector)[!names(p_vector) %in% names(filled_transform)]] <- "identity"
      filled_transform <- filled_transform[names(p_vector)]
    } else{
      filled_transform <- setNames(rep("identity", length(p_vector)), names(p_vector))
    }
  }
  filled_transform[names(transform$func)] <- transform$func
  filled_lower <- setNames(rep(-Inf,length(filled_transform)),names(filled_transform))
  filled_upper <- setNames(rep(Inf,length(filled_transform)),names(filled_transform))
  if(!is_pre){
    filled_lower[names(model()$transform$lower)] <- model()$transform$lower
    filled_upper[names(model()$transform$upper)] <- model()$transform$upper
  } else{
    filled_lower[names(model()$pre_transform$lower)] <- model()$pre_transform$lower
    filled_upper[names(model()$pre_transform$upper)] <- model()$pre_transform$upper
  }

  filled_lower[filled_transform %in% has_lower] <- 0
  filled_upper[filled_transform %in% has_upper] <- 1
  if (!is.null(transform$lower)) filled_lower[names(transform$lower)] <- transform$lower
  if (!is.null(transform$upper)) filled_lower[names(transform$upper)] <- transform$upper
  list(func=filled_transform,lower=filled_lower,upper=filled_upper)
}


fill_bound <- function(bound, model) {
  filled_bound <- model()$bound
  if (!is.null(bound)) {
    if (names(bound)[1] != "minmax")
      stop("first entry of bound must be named minmax")
    if (!all(colnames(bound$minmax) %in% names(model()$p_types)))
      stop("minmax column names must correspond to parameter types")
    if (!is.null(bound$exception) &&
        (!all(names(bound$exception) %in% names(model()$p_types))))
      stop("exception names must correspond to parameter types")
    filled_bound$minmax[,colnames(bound$minmax)] <- bound$minmax
    if (!is.null(bound$exception)) {
      filled_bound$exception <- c(bound$exception,filled_bound$exception)
      filled_bound$exception <- filled_bound$exception[!duplicated(names(filled_bound$exception))]
    }
  }
  filled_bound
}

# This form used in random number generation
do_bound <- function(pars,bound) {
  tpars <- t(pars[,colnames(bound$minmax),drop=FALSE])
  ok <- tpars > bound$minmax[1,] & tpars < bound$minmax[2,]
  if (!is.null(bound$exception)) ok[names(bound$exception),] <-
    ok[names(bound$exception),] |
    (tpars[names(bound$exception),] == bound$exception)
  apply(ok,2,all)
}

# This form used in get_pars
add_bound <- function(pars,bound) {
  attr(pars,"ok") <- do_bound(pars,bound)
  pars
}




