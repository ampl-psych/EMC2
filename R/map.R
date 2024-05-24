get_pars <- function(p_vector,dadm) {
  # Add constants, transform p_vector, map to design, transform mapped parameters
  # to the natural scale, and create trial-dependent parameters. Ordinal
  # parameters are first exponentiated.

  if (!is.null(attr(dadm,"ordinal")))
    if (is.matrix(p_vector))
       p_vector[,attr(dadm,"ordinal")] <- exp(p_vector[,attr(dadm,"ordinal")]) else
       p_vector[attr(dadm,"ordinal")] <- exp(p_vector[attr(dadm,"ordinal")])

  attr(dadm,"model")()$Ttransform(
    attr(dadm,"model")()$Ntransform(
      map_p(
        attr(dadm,"model")()$transform(add_constants(p_vector,attr(dadm,"constants"))),
        dadm)),
    dadm)
}

get_design <- function(samples)
  # prints out design from samples object
{
  design <- attr(samples,"design_list")[[1]]
  model <- attr(samples, "model_list")[[1]]
  design$Ffactors$subjects <- design$Ffactors$subjects[1]
  dadm <- design_model(make_data(sampled_p_vector(design,model),design,n_trials=1),design,model,
                       rt_check=FALSE,compress=FALSE)
  dadm[,!(names(dadm) %in% c("subjects","trials","R","rt","winner"))]
}

get_design_matrix <- function(samples){
  attr(sampled_p_vector(attr(samples,"design_list")[[1]]),"map")
}

get_map <- function(samples,add_design=TRUE) {
  out <- attr(attr(attr(samples,"design_list")[[1]],"p_vector"),"map")
  if (add_design) {
    mnl <- mapped_name_list(attr(samples,"design_list")[[1]],attr(samples,"model_list")[[1]],TRUE)
    for (i in names(out)) out[[i]] <- cbind(mnl[[i]],out[[i]])
  }
  out
}

pmat <- function(p_vector,design)
  # puts vector form of p_vector into matrix form
{
  ss <- design$Ffactors$subjects
  matrix(rep(p_vector,each=length(ss)),nrow=length(ss),
         dimnames=list(ss,names(p_vector)))
}



mapped_name_list <- function(design,model,save_design=FALSE)
  # makes a list, with entries for each parameter type, of names for mapped
  # parameters or with unique design columns
{
  doMap <- function(mapi,pmat) t(mapi %*% t(pmat[,dimnames(mapi)[[2]],drop=FALSE]))

  constants <- design$constants
  p_vector <- attr(design,"p_vector")
  mp <- mapped_par(p_vector,design)
  map <- attr(sampled_p_vector(design),"map")
  pmat <- model()$transform(add_constants(t(as.matrix(p_vector)),constants))
  plist <- lapply(map,doMap,pmat=pmat)
  if (model()$type=="SDT") {
    ht <- apply(map$threshold[,grepl("lR",dimnames(map$threshold)[[2]]),drop=FALSE],1,sum)
    plist$threshold <- plist$threshold[,ht!=max(ht),drop=FALSE]
  }
  # Give mapped variables names and remove constants
  for (i in 1:length(plist)) {
    vars <- row.names(attr(terms(design$Flist[[i]]),"factors"))
    if (is.null(vars)) dimnames(plist[[i]])[2] <- names(plist)[i] else {
      uniq <- !duplicated(apply(mp[,vars],1,paste,collapse="_"))
      if (save_design) plist[[i]] <- mp[uniq,vars[-1]] else
        dimnames(plist[[i]])[[2]] <-
          paste(vars[1],apply(mp[uniq,vars[-1],drop=FALSE],1,paste,collapse="_"),sep="_")
    }
  }
  if (save_design) plist else lapply(plist,function(x){dimnames(x)[[2]]})
}


make_pmat <- function(p_vector,design)
  # puts vector form of p_vector into matrix form
{
  ss <- design$Ffactors$subjects
  matrix(rep(p_vector,each=length(ss)),nrow=length(ss),
         dimnames=list(ss,names(p_vector)))
}

#' Augments parameter matrix or vector p with constant parameters (also used in data)
#'
#' @param p either a matrix or vector of parameters
#' @param constants a named vector of constants
#'
#' @return a matrix or vector, depending on input, with the varying parameters and constants combined.
#'
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

add_constants_mcmc <- function(p,constants){
  return(mcmc(add_constants(p,constants)))
}

#' Parameter mapping back to the design factors
#'
#' Maps a parameter vector that corresponds to sampled parameters
#' of the cognitive model back to the experimental design. The parameter vector
#' can be created using ``sampled_p_vector()``. The returned matrix shows whether/how parameters
#' differ across the experimental factors.
#'
#' @param p_vector A parameter vector. Must be in the form of ``sampled_p_vector(design)``
#' @param design A design list. Created ``by make_design``
#' @param model Optional model type (if not already specified in ``design``)
#' @param digits Integer. Will round the output parameter values to this many decimals
#' @param ... optional arguments
#' @param remove_subjects Boolean. Whether to include subjects as a factor in the design
#'
#' @return Matrix with a column for each factor in the design and for each model parameter type (``p_type``).
#' @examples
#' # First define a design:
#' design_DDMaE <- make_design(data = forstmann,model=DDM,
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
                       digits=3,remove_subjects=TRUE, ...)
  # Show augmented data and corresponding mapped parameter
{
  remove_RACE <- TRUE
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  Fcovariates <- design$Fcovariates
  if (is.null(model)) if (is.null(design$model))
  stop("Must specify model as not in design") else model <- design$model
  if (remove_subjects) design$Ffactors$subjects <- design$Ffactors$subjects[1]
  if (!is.matrix(p_vector)) p_vector <- make_pmat(p_vector,design)
  dadm <- design_model(make_data(p_vector,design,n_trials=1,Fcovariates=Fcovariates),
                       design,model,rt_check=FALSE,compress=FALSE)
  ok <- !(names(dadm) %in% c("subjects","trials","R","rt","winner"))
  out <- cbind(dadm[,ok],round(get_pars(p_vector,dadm),digits))
  if (model()$type=="SDT")  out <- out[dadm$lR!=levels(dadm$lR)[length(levels(dadm$lR))],]
  if (model()$type=="DDM")  out <- out[,!(names(out) %in% c("lR","lM"))]
  if (any(names(out)=="RACE") && remove_RACE)
    out <- out[as.numeric(out$lR) <= as.numeric(as.character(out$RACE)),,drop=FALSE]
  return(out)
}


# mcmc=mcmcList[[1]]; design=attr(samplers,"design_list")[[1]];model=attr(samplers,"model_list")[[1]]
map_mcmc <- function(mcmc,design,model, include_constants = TRUE)
  # Maps vector or matrix (usually mcmc object) of sampled parameters to native
  # model parameterization.
{
  doMap <- function(mapi,pmat) t(mapi %*% t(pmat[,dimnames(mapi)[[2]],drop=FALSE]))

  get_p_types <- function(nams)
    unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))

  if (!is.matrix(mcmc)) mcmc <- t(as.matrix(mcmc))

  if (!is.null(attr(design,"ordinal")))
    mcmc[,attr(design,"ordinal")] <- exp(mcmc[,attr(design,"ordinal")])

  map <- attr(sampled_p_vector(design, add_da = TRUE, all_cells_dm = TRUE),"map")
  constants <- design$constants
  mp <- mapped_par(mcmc[1,],design,remove_RACE=FALSE)
  pmat <- model()$transform(add_constants(mcmc,constants))
  plist <- lapply(map,doMap,pmat=pmat)
  if (model()$type=="SDT") {
    ht <- apply(map$threshold[,grepl("lR",dimnames(map$threshold)[[2]]),drop=FALSE],1,sum)
    plist$threshold <- plist$threshold[,ht!=max(ht),drop=FALSE]
  }
  # Give mapped variables names and flag constant
  isConstant <- NULL
  for (i in 1:length(plist)) {
    vars <- row.names(attr(terms(design$Flist[[i]]),"factors"))
    uniq <- !duplicated(apply(mp[,vars],1,paste,collapse="_"))
    if (is.null(vars)) dimnames(plist[[i]])[2] <- names(plist)[i] else {
      dimnames(plist[[i]])[[2]] <-
        paste(vars[1],apply(mp[uniq,vars[-1],drop=FALSE],1,paste,collapse="_"),sep="_")
    }
    if (dim(plist[[i]])[1]!=1) isConstant <- c(isConstant,
                                               apply(plist[[i]],2,function(x){all(x[1]==x[-1])}))
  }
  pmat <- do.call(cbind,plist)
  cnams <- dimnames(pmat)[[2]]
  dimnames(pmat)[[2]] <- get_p_types(cnams)
  out <- model()$Ntransform(pmat)[,1:length(cnams)]
  dimnames(out)[[2]] <- cnams
  if(!include_constants) out <- out[,!isConstant]
  out <- as.mcmc(out)
  attr(out,"isConstant") <- isConstant
  out
}

add_constants_mcmc <- function(p,constants){
  return(mcmc(add_constants(p,constants)))
}
