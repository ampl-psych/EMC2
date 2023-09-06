#' Binds together elements that make up a model design as a list
#'
#' This function returns a list with the model design.
#'
#' @param Flist A list. Contains the design formula of the type list(y ~ x).
#' @param Ffactors A named list. Participant factor and all factors used in design formulas that are defined as factors in the data frame to be fit.
#'
#' Must use "subjects" as name for participant factor. For others don't use "trial", "R", "rt", "lR", or "lM".
#' Also refrain from using "_" as factor names.
#'
#'
#' @param Rlevels A character vector. Contains the response factor levels.
#' @param model A function, Specifies the model type.
#' @param ddata A data frame that is used to determine Ffactors, Rlevels and Fcovariates,
#' in which case these arguments can be left NULL. Any numeric column except trials and
#' rt are made into covariates, R factor column makes Rlevels, and remaining factor
#' columns are used in Ffactors.
#' @param Clist A list. Contrast list.
#' @param matchfun A function. Specifies whether a response was correct or not.
#'
#' Example: `function(d)d$S==d$lR`
#'
#' @param constants A named vector. Sets constants. Any parameter named by `sampled_p_vector` can be set constant.
#' @param Fcovariates Covariate factors. Covariate measures which may be included in the model and/or mapped to factors.
#' @param Ffunctions Factor Functions. Functions to specify specific parameterisations, or functions of parameters.
#' @param adapt For future compatibility. Ignore.
#' @param report_p_vector TRUE (default) or FALSE. if TRUE, returns the vector of parameters to be estimated.
#' @param custom_p_vector A character vector. If specified, a custom likelihood function can be supplied.
#'
#' @return A design list.
#' @export
#'
#'
make_design <- function(Flist = NULL,Ffactors = NULL,Rlevels = NULL,model,ddata=NULL,
                        Clist=NULL,matchfun=NULL,constants=NULL,Fcovariates=NULL,Ffunctions=NULL,
                        adapt=NULL,report_p_vector=TRUE, custom_p_vector = NULL){
  if(!is.null(custom_p_vector)){
    design <- list(Flist = Flist, model = model, Ffactors = Ffactors)
    attr(design, "sampled_p_names") <-custom_p_vector
    attr(design, "custom_ll") <- TRUE
    return(design)
  }
  if (!is.null(ddata)) {
    facs <- lapply(ddata,levels)
    nfacs <- facs[unlist(lapply(facs,is.null))]
    facs <- facs[!unlist(lapply(facs,is.null))]
    Rlevels <- facs[["R"]]
    Ffactors <- facs[names(facs)!="R"]
    nfacs <- nfacs[!(names(nfacs) %in% c("trials","rt"))]
    if (length(nfacs)>0) Fcovariates <- nfacs
  }
  # # Frees up memory again by creating new enclosing environments, courtesy of Steven
  # if(!is.null(Ffunctions)){
  #   Ffunctions <- lapply(Ffunctions, function(f) {environment(f) <- new.env(parent=globalenv()); return(f)})
  # }
  # if(!is.null(Flist)) {
  #   Flist <- lapply(Flist, function(f) {environment(f) <- new.env(parent=globalenv()); return(f)})
  # }
  # if(!is.null(model)) environment(model) <- new.env(parent=globalenv())
  # if(!is.null(matchfun)) environment(matchfun) <- new.env(parent=globalenv())
  if (model()$type=="SDT") {
    Clist[["lR"]] <- contr.increasing(length(Rlevels),Rlevels)
  }

  design <- list(Flist=Flist,Ffactors=Ffactors,Rlevels=Rlevels,
                 Clist=Clist,matchfun=matchfun,constants=constants,
                 Fcovariates=Fcovariates,Ffunctions=Ffunctions,adapt=adapt,model=model)
  p_vector <- sampled_p_vector(design,design$model)

  if (model()$type=="SDT") {
    tnams <- dimnames(attr(p_vector,"map")$threshold)[[2]]
    max_threshold=paste0("lR",Rlevels[length(Rlevels)])
    tnams <- tnams[grepl(max_threshold,tnams)]
    if (!any(tnams %in% names(constants))) {
      design$constants <- stats::setNames(c(constants,rep(log(1e100),length(tnams))),
                                   c(names(constants),tnams))
      p_vector <- sampled_p_vector(design,design$model)
    }
  }
  attr(design,"p_vector") <- p_vector

  if (report_p_vector) {
    print(p_vector)
  }

  return(design)
}


#' contr.increasing()
#'
#' special contrast for SDT threshold factor
#' first = intercept, cumsum other (positive) levels to force non-decreasing
#'
#' @param n an integer. The number of items for which to create the contrast.
#' @param levels Character vector. the factor levels which will be the colnames of the returning matrix.
#'
#' @return a contrast matrix.
#' @export

contr.increasing <- function(n,levels=NULL)
{
  contr <- matrix(0,nrow=n,ncol=n-1,dimnames=list(NULL,2:n))
  contr[lower.tri(contr)] <- 1
  if (!is.null(levels)) dimnames(contr)[[2]] <- levels[-1]
  contr
}

#' contr.anova()
#'
#' orthogonal helmert contrast scaled to estimate differences between conditions. Use in make_design.
#'
#' @param n an integer. the number of items for which to create the contrast
#'
#' @return a contrast matrix.
#' @export

contr.anova <- function(n) {
  contr <- stats::contr.helmert(n)
  contr/rep(2*apply(abs(contr),2,max),each=dim(contr)[1])
}



#' sampled_p_vector()
#'
#' Makes an empty p_vector corresponding to model.
#' matchfun only needed in design if uses lM factor
#'
#' @param design a list of the design made with make_design.
#' @param model a model list. Default is the model specified in the design list.
#' @param doMap logical. If TRUE will
#'
#' @return Named vector with mapping attributes.
#' @export

sampled_p_vector <- function(design,model=NULL,doMap=TRUE)
  # Makes an empty p_vector corresponding to model.
{
  if (is.null(model)) model <- design$model
  if (is.null(model)) stop("Must supply model as not in design")

  Ffactors=c(design$Ffactors,list(R=design$Rlevels))
  data <- as.data.frame.table(array(dim=unlist(lapply(Ffactors,length)),
                                    dimnames=Ffactors))[,-(length(Ffactors)+1)]
  for (i in names(design$Ffactors))
    data[[i]] <- factor(data[[i]],levels=design$Ffactors[[i]])

  # if (!is.null(design$Ffunctions))
  #   data <- cbind.data.frame(data,data.frame(lapply(design$Ffunctions,function(f){f(data)})))

  if (!is.null(design$Fcovariates)) {
    covs <- matrix(0,nrow=dim(data)[1],ncol=length(design$Fcovariates),
                   dimnames=list(NULL,names(design$Fcovariates)))
    data <- cbind.data.frame(data,covs)
  }
  dadm <- design_model(
    add_accumulators(data,matchfun=design$matchfun,type=model()$type,Fcovariates=design$Fcovariates),
    design,model,add_acc=FALSE,verbose=FALSE,rt_check=FALSE,compress=FALSE)
  sampled_p_names <- attr(dadm,"sampled_p_names")
  out <- stats::setNames(numeric(length(sampled_p_names)),sampled_p_names)
  if (doMap) attr(out,"map") <-
    lapply(attributes(dadm)$designs,function(x){x[,,drop=FALSE]})
  out
}


add_accumulators <- function(data,matchfun=NULL,simulate=FALSE,type="RACE", Fcovariates=NULL) {
  if (!is.factor(data$R)) stop("data must have a factor R")
  factors <- names(data)[!names(data) %in% c("R","rt","trials",Fcovariates)]
  if (type %in% c("RACE","SDT")) {
    datar <- cbind(do.call(rbind,lapply(1:length(levels(data$R)),function(x){data})),
                   lR=factor(rep(levels(data$R),each=dim(data)[1]),levels=levels(data$R)))
    if (!is.null(matchfun)) {
      lM <- matchfun(datar)
      if (any(is.na(lM)) || !(is.logical(lM)))
        stop("matchfun not scoring properly")
      datar$lM <- factor(lM)
    }
  } else datar <- cbind(data,lR=factor(rep(levels(data$R)[1],dim(data)[1]),
                                       levels=levels(data$R)),lM=factor(rep(TRUE,dim(data)[1]))) # For DDM
  row.names(datar) <- NULL
  if (simulate) datar$rt <- NA else {
    datar$winner <- datar$lR==datar$R
    datar$winner[is.na(datar$winner)] <- FALSE
  }
  # sort cells together
  if ("trials" %in% names(data)){
    if(length(factors) > 1){
      datar[order(apply(datar[,c(factors)],1,paste,collapse="_"), as.numeric(datar$trials),as.numeric(datar$lR)),]
    } else{
      datar[order(datar[,c(factors)], as.numeric(datar$trials),as.numeric(datar$lR)),]
    }
  }
  else{
    if(length(factors) > 1){
      datar[order(apply(datar[,c(factors)],1,paste,collapse="_"), as.numeric(datar$lR)),]
    } else{
      datar[order(datar[,c(factors)], as.numeric(datar$lR)),]
    }
  }
}

design_model_custom_ll <- function(data, design, model){
  if (!is.factor(data$subjects)) {
    data$subjects <- factor(data$subjects)
    warning("subjects column was converted to a factor")
  }
  dadm <- data
  model_input <- model
  attr(dadm, "model") <- function(){
    return(list(log_likelihood = model_input))
  }
  attr(dadm,"sampled_p_names") <- attr(design, "sampled_p_names")
  attr(dadm, "custom_ll") <- TRUE
  return(dadm)

}


#' Combines a data frame with a design to create
#' a data augmented design model ("dadm") object. Augmentation
#' refers to replicating the data with one row for each accumulator.
#'
#' Usually called by make_samplers rather than directly by the user, except
#' where a dadm is needed for use with profile_pmwg.
#'
#' Performs a series to checks to make sure data frame and design match and
#' (by default) augments the data frame by adding accumulator factors and
#' compresses the result for efficient likelihood calculation.
#'
#'
#' @param data  data frame
#' @param design matching design
#' @param model if model not an attribute of design can be supplied here
#' @param add_acc default TRUE creates and 'augmented' data frame, which
#' replicates and stacks the supplied data frame for each accumulator and
#' adds factors to represent accumulators: lR = latent response, with a level
#' for each response (R) represented by the accumulator and lM = latent match,
#' a logical indicating if the accumulator represents the correct response as
#' specified by the design's matchfun.
#' @param rt_resolution maximum resolution of rt, NULL = no rounding
#' @param verbose if true reports compression outcome
#' @param compress default TRUE only keeps unique rows in terms of all
#' parameter design matrices R, lR and rt (at a given resolution)
#' @param rt_check checks if any truncation and censoring specified in the design
#' are respected.
#'
#' @return a (possibly) augmented and compressed data frame with attributes
#' specifying the design and how to decompress ready supporting likelihood
#' computation
#' @export

design_model <- function(data,design,model=NULL,
                         add_acc=TRUE,rt_resolution=0.02,verbose=TRUE,
                         compress=TRUE,rt_check=TRUE)
  # Flist is a list of formula objects, one for each p_type
  # da is augmented data (from add_accumulators), must have all of the factors
  #   and covariates that are used in formulas
  # Clist is a list of either a single unnamed contrast (as in the default)
  #   lists, one for each model()$p_type (allowed to have some p_types missing).
  # These elements of model define the paramaterization being used
  #   ptypes defines the parameter types for which designs must be specified
  #   transform if a function acting on p_vector before mapping
  #   Ntransform is a function acting on the output of map_p
{
  compress_dadm <- function(da,designs,Fcov,Ffun)
    # out keeps only unique rows in terms of all parameters design matrices
    # R, lR and rt (at given resolution) from full data set
  {
    nacc <- length(unique(da$lR))
    # contract output
    cells <- paste(
      apply(do.call(cbind,lapply(designs,function(x){
        apply(x[attr(x,"expand"),,drop=FALSE],1,paste,collapse="_")})
      ),1,paste,collapse="+"),da$subjects,da$R,da$lR,da$rt,sep="+")
    # Make sure that if row is included for a trial so are other rows
    if (nacc>1) cells <- paste0(rep(apply(matrix(cells,nrow=nacc),2,paste0,collapse="_"),
                      each=nacc),rep(1:nacc,times=length(cells)/nacc),sep="_")
    if (!is.null(Fcov))
      cells <- paste(cells,apply(da[,names(Fcov),drop=FALSE],1,paste,collapse="+"),sep="+")
    if (!is.null(Ffun))
      cells <- paste(cells,apply(da[,Ffun,drop=FALSE],1,paste,collapse="+"),sep="+")
    contract <- !duplicated(cells)
    out <- da[contract,,drop=FALSE]
    attr(out,"contract") <- contract
    attr(out,"expand") <- as.numeric(factor(cells,levels=unique(cells)))
    lR1 <- da$lR==levels(da$lR)[[1]]
    attr(out,"expand_winner") <- as.numeric(factor(cells[lR1],levels=unique(cells[lR1])))
    attr(out,"s_expand") <- da$subjects
    attr(out,"designs") <- lapply(designs,function(x){
      attr(x,"expand") <- attr(x,"expand")[contract]; x})

    # indices to use to contract further ignoring rt then expand back
    cells_nort <- paste(
      apply(do.call(cbind,lapply(designs,function(x){
        apply(x[attr(x,"expand"),,drop=FALSE],1,paste,collapse="_")})
      ),1,paste,collapse="+"),da$subjects,da$R,da$lR,sep="+")[contract]
    attr(out,"unique_nort") <- !duplicated(cells_nort)
    # Only first level WHY????
    cells <- cells[da$lR==levels(da$lR)[1]]
    cells_nort <- cells_nort[out$lR==levels(out$lR)[1]]
    attr(out,"expand_nort") <- as.numeric(factor(cells_nort,
       levels=unique(cells_nort)))[as.numeric(factor(cells,levels=unique(cells)))]

    # indices to use to contract ignoring rt and response (R), then expand back
    cells_nortR <- paste(apply(do.call(cbind,lapply(designs,function(x){
      apply(x[attr(x,"expand"),,drop=FALSE],1,paste,collapse="_")})),1,paste,collapse="+"),
      da$subjects,da$lR,sep="+")[contract]
    attr(out,"unique_nortR") <- !duplicated(cells_nortR)
    # Only first level WHY????
    cells_nortR <- cells_nortR[out$lR==levels(out$lR)[1]]
    attr(out,"expand_nortR") <- as.numeric(factor(cells_nortR,
       levels=unique(cells_nortR)))[as.numeric(factor(cells,levels=unique(cells)))]

    # Lower censor
    if (!any(is.na(out$rt))) { # Not an choice only model
      winner <- out$lR==levels(out$lR)[[1]]
      ok <- out$rt[winner]==-Inf
      if (any(ok)) {
        ok[ok] <- 1:sum(ok)
        attr(out,"expand_lc") <- ok[attr(out,"expand_winner")] + 1
      }
      # Upper censor
      ok <- out$rt[winner]==Inf
      if (any(ok)) {
        ok[ok] <- 1:sum(ok)
        attr(out,"expand_uc") <- ok[attr(out,"expand_winner")] + 1
      }
    }
    out
  }


  check_rt <- function(b,d,upper=TRUE)
    # Check bounds respected if present
  {
    if (!all(sort(levels(d$subjects))==sort(names(b))))
      stop("Bound vector must have same names as subjects")
    d <- d[!is.na(d$rt),]
    d <- d[is.finite(d$rt),]
    bound <- d$subjects
    levels(bound) <- unlist(b)[levels(bound)]
    if (upper)
      ok <- all(d$rt < as.numeric(as.character(bound))) else
        ok <- all(d$rt > as.numeric(as.character(bound)))
    if (!all(ok)) stop("Bound not respected in data")
  }

  if (is.null(model)) {
    if (is.null(design$model))
      stop("Model must be supplied if it has not been added to design")
    model <- design$model
  }
  if (model()$type=="SDT") rt_check <- FALSE
  if(model()$type == "MRI"){
    rt_check <- FALSE
    add_acc <- FALSE
    compress <- FALSE
  }
  if (!is.null(design$adapt)) compress=FALSE # RL models
  if (any(model()$p_types %in% names(data)))
    stop("Data cannot have columns with the same names as model parameters")
  if (!is.factor(data$subjects)) {
    data$subjects <- factor(data$subjects)
    warning("subjects column was converted to a factor")
  }

  if (!any(names(data)=="trials")) data$trials <- 1:dim(data)[1]

  if (rt_check) {
    # Truncation
    if (!is.null(attr(data,"UT"))) {
      if (length(attr(data,"UT"))==1 && is.null(names(attr(data,"UT"))))
        attr(data,"UT") <- stats::setNames(rep(attr(data,"UT"),length(levels(data$subjects))),
                                    levels(data$subjects))
      check_rt(attr(data,"UT"),data)
    }
    if (!is.null(attr(data,"LT"))) {
      if (length(attr(data,"LT"))==1 && is.null(names(attr(data,"LT"))))
        attr(data,"LT") <- stats::setNames(rep(attr(data,"LT"),length(levels(data$subjects))),
                                    levels(data$subjects))
      if (any(attr(data,"LT")<0)) stop("Lower truncation cannot be negative")
      check_rt(attr(data,"LT"),data,upper=FALSE)
    }
    if (!is.null(attr(data,"UT")) & !is.null(attr(data,"LT"))) {
      DT <- attr(data,"UT") - attr(data,"LT")
      if (!is.null(DT) && any(DT<0)) stop("UT must be greater than LT")
    }

    # Censoring
    if (!is.null(attr(data,"UC"))) {
      if (length(attr(data,"UC"))==1 && is.null(names(attr(data,"UC"))))
        attr(data,"UC") <- stats::setNames(rep(attr(data,"UC"),length(levels(data$subjects))),
                                    levels(data$subjects))
      check_rt(attr(data,"UC"),data)
      if (!is.null(attr(data,"UT")) && attr(data,"UT") < attr(data,"UC"))
        stop("Upper censor must be less than upper truncation")
    }
    if (!is.null(attr(data,"LC"))) {
      if (length(attr(data,"LC"))==1 && is.null(names(attr(data,"LC"))))
        attr(data,"LC") <- stats::setNames(rep(attr(data,"LC"),length(levels(data$subjects))),
                                    levels(data$subjects))
      if (any(attr(data,"LC")<0)) stop("Lower censor cannot be negative")
      check_rt(attr(data,"LC"),data,upper=FALSE)
      if (!is.null(attr(data,"LT")) && attr(data,"LT") > attr(data,"LC"))
        stop("Lower censor must be greater than lower truncation")
    }
    if (!any(data$rt[!is.na(data$rt)]==-Inf) & !is.null(attr(data,"LC")))
      attr(data,"LC") <- NULL else
        if (any(data$rt[!is.na(data$rt)]==-Inf) & is.null(attr(data,"LC")))
          stop("Data must have an LC attribute if any rt = -Inf")
    if (!any(data$rt[!is.na(data$rt)]==Inf) & !is.null(attr(data,"UC")))
      attr(data,"LC") <- NULL else
        if (any(data$rt[!is.na(data$rt)]==Inf) & is.null(attr(data,"UC")))
          stop("Data must have an UC attribute if any rt = Inf")
    if (!is.null(attr(data,"UC"))) check_rt(attr(data,"UC"),data)
    if (!is.null(attr(data,"LC"))) check_rt(attr(data,"LC"),data,upper=FALSE)
    if (!is.null(attr(data,"UC")) & !is.null(attr(data,"LC"))) {
      DC <- attr(data,"UC") - attr(data,"LC")
      if (!is.null(DC) && any(DC<0)) stop("UC must be greater than LC")
    }
  }

  if (!add_acc) da <- data else
    da <- add_accumulators(data,design$matchfun,type=model()$type,Fcovariates=design$Fcovariates)
  da <- da[order(da$subjects),] # fixes different sort in add_accumulators depending on subject type

  if (!is.null(design$Ffunctions)) for (i in names(design$Ffunctions)) {
    newF <- stats::setNames(data.frame(design$Ffunctions[[i]](da)),i)
    da <- cbind.data.frame(da,newF)
  }

  if (is.null(model()$p_types) | is.null(model()$transform) |
      is.null(model()$Ntransform) | is.null(model()$Ttransform))
    stop("p_types, transform and Ntransform must be supplied")
  if (!all(unlist(lapply(design$Flist,class))=="formula"))
    stop("Flist must contain formulas")
  nams <- unlist(lapply(design$Flist,function(x)as.character(stats::terms(x)[[2]])))
  names(design$Flist) <- nams
  if (!all(sort(model()$p_types)==sort(nams)) & model()$type != "MRI")
    stop("Flist must specify formulas for ",paste(model()$p_types,collapse = " "))
  if (is.null(design$Clist)) design$Clist=list(stats::contr.treatment)
  if (!is.list(design$Clist)) stop("Clist must be a list")
  if (!is.list(design$Clist[[1]])[1]) # same contrasts for all p_types
    design$Clist <- stats::setNames(lapply(1:length(model()$p_types),
                                    function(x)design$Clist),model()$p_types) else {
                                      missing_p_types <- model()$p_types[!(model()$p_types %in% names(design$Clist))]
                                      if (length(missing_p_types)>0) {
                                        nok <- length(design$Clist)
                                        for (i in 1:length(missing_p_types)) {
                                          design$Clist[[missing_p_types[i]]] <- list(stats::contr.treatment)
                                          names(design$Clist)[nok+i] <- missing_p_types[i]
                                        }
                                      }
                                    }
  if(model()$type != "MRI") for (i in model()$p_types) attr(design$Flist[[i]],"Clist") <- design$Clist[[i]]
  out <- lapply(design$Flist,make_dm,da=da,Fcovariates=design$Fcovariates)
  if (!is.null(rt_resolution) & !is.null(da$rt)) da$rt <- round(da$rt/rt_resolution)*rt_resolution
  if (compress) dadm <- compress_dadm(da,designs=out,
                                      Fcov=design$Fcovariates,Ffun=names(design$Ffunctions)) else {
                                        dadm <- da
                                        attr(dadm,"designs") <- out
                                        attr(dadm,"s_expand") <- da$subjects
                                        attr(dadm,"expand") <- 1:dim(dadm)[1]
                                      }
  if (verbose & compress) message("Likelihood speedup factor: ",
    round(dim(da)[1]/dim(dadm)[1],1)," (",dim(dadm[as.logical(dadm$lM),])[1]," unique trials)")
  p_names <-  unlist(lapply(out,function(x){dimnames(x)[[2]]}),use.names=FALSE)

  bad_constants <- names(design$constants)[!(names(design$constants) %in% p_names)]
  if (length(bad_constants) > 0)
    stop("Constant(s) ",paste(bad_constants,collapse=" ")," not in design")

  # Pick out constants
  sampled_p_names <- p_names[!(p_names %in% names(design$constants))]

  attr(dadm,"p_names") <- p_names
  attr(dadm,"model") <- model
  attr(dadm,"constants") <- design$constants
  attr(dadm,"sampled_p_names") <- sampled_p_names
  attr(dadm, "LT") <- attr(data,"LT")
  attr(dadm, "UT") <- attr(data,"UT")
  attr(dadm, "LC") <- attr(data,"LC")
  attr(dadm, "UC") <- attr(data,"UC")
  if (add_acc) {
    attr(dadm, "ok_dadm_winner") <- is.finite(dadm$rt) & dadm$winner
    attr(dadm, "ok_dadm_looser") <- is.finite(dadm$rt) & !dadm$winner
    attr(dadm, "ok_da_winner") <- attr(dadm, "ok_dadm_winner")[attr(dadm,"expand")]
    attr(dadm, "ok_da_looser") <- attr(dadm, "ok_dadm_looser")[attr(dadm,"expand")]
  }
  attr(dadm,"ok_trials") <- is.finite(data$rt)
  attr(dadm,"s_data") <- data$subjects
  # if (!is.null(design$adapt)) {
  #   attr(dadm,"adapt") <- stats::setNames(
  #     lapply(levels(dadm$subjects),augment,da=dadm,design=design),
  #     levels(dadm$subjects))
  #   attr(dadm,"adapt")$design <- design$adapt
  # }
  dadm
}


make_dm <- function(form,da,Clist=NULL,Fcovariates=NULL)
  # Makes a design matrix based on formula form from augmented data frame da
{

  compress_dm <- function(dm)
    # out keeps only unique rows, out[attr(out,"expand"),] gets back original.
  {
    cells <- apply(dm,1,paste,collapse="_")
    ass <- attr(dm,"assign")
    contr <- attr(dm,"contrasts")
    dups <- duplicated(cells)
    out <- dm[!dups,,drop=FALSE]
    attr(out,"expand") <- as.numeric(factor(cells,levels=unique(cells)))
    attr(out,"assign") <- ass
    attr(out,"contrasts") <- contr
    out
  }

  if (is.null(Clist)) Clist <- attr(form,"Clist")
  pnam <- stats::terms(form)[[2]]
  da[[pnam]] <- 1
  for (i in names(Clist)) if (i %in% names(da)) {
    if (!is.factor(da[[i]]))
      stop(i," must be a factor (design factors has a parameter name?)")
    levs <- levels(da[[i]])
    nl <- length(levs)
    if (class(Clist[[i]])[1]=="function")
      stats::contrasts(da[[i]]) <- do.call(Clist[[i]],list(n=levs)) else {
        if (!is.matrix(Clist[[i]]) || dim(Clist[[i]])[1]!=nl) {
          if (all(levs %in% row.names(Clist[[i]]))) # design with missing cells
            Clist[[i]] <- Clist[[i]][levs,] else
            stop("Clist for ",i," not a ",nl," row matrix")
        } else dimnames(Clist[[i]])[[1]] <- levs
        stats::contrasts(da[[i]],how.many=dim(Clist[[i]])[2]) <- Clist[[i]]
      }
  }
  out <- stats::model.matrix(form,da)
  if (dim(out)[2]==1) dimnames(out)[[2]] <- as.character(pnam) else {
    if (attr(stats::terms(form),"intercept")!=0) {
      cnams <- paste(pnam,dimnames(out)[[2]][-1],sep="_")
      dimnames(out)[[2]] <- c(pnam,cnams)
    } else dimnames(out)[[2]] <- paste(pnam,dimnames(out)[[2]],sep="_")
  }
  compress_dm(out)
}


#### Functions to look at parameters ----

map_p <- function(p,dadm)
  # Map p to dadm and returns matrix of mapped parameters
  # p is either a vector or a matrix (ncol = number of subjects) of p_vectors
{

  if ( is.matrix(p) ) {
    if (!all(sort(dimnames(p)[[2]])==sort(attr(dadm,"p_names"))))
      stop("p col.names must be: ",paste(attr(dadm,"p_names"),collapse=", "))
    if (!all(levels(dadm$subjects) %in% dimnames(p)[[1]]))
      stop("p must have rows named for every subject in dadm")
    p <- p[dadm$subjects,]
  } else if (!all(sort(names(p))==sort(attr(dadm,"p_names"))))
    stop("p names must be: ",paste(attr(dadm,"p_names"),collapse=", "))

  pars <- matrix(nrow=dim(dadm)[1],ncol=length(attr(dadm,"model")()$p_types),
                 dimnames=list(NULL,attr(dadm,"model")()$p_types))
  for (i in attr(dadm,"model")()$p_types) {
    if ( !is.matrix(p) ) {
      pm <- t(as.matrix(p[dimnames(attr(dadm,"designs")[[i]])[[2]]]))
      pm <- pm[rep(1,dim(pars)[1]),]
    } else pm <- p[,dimnames(attr(dadm,"designs")[[i]])[[2]],drop=FALSE]
    tmp <- pm*attr(dadm,"designs")[[i]][attr(attr(dadm,"designs")[[i]],"expand"),,drop=FALSE]
    tmp[is.nan(tmp)] <- 0 # 0 weight x Inf parameter fix
    pars[,i] <- apply(tmp,1,sum)
  }
  pars
}

# data generation

# Used in make_data and make_samplers
add_trials <- function(dat)
  # Add trials column, 1:n for each subject
{
  n <- table(dat$subjects)
  if (!any(names(dat)=="trials")) dat <- cbind.data.frame(dat,trials=NA)
  for (i in names(n)) dat$trials[dat$subjects==i] <- 1:n[i]
  dat
}

dm_list <- function(dadm)
  # Makes data model into subjects list for use by likelihood
  # Assumes each subject has the same design.
{

  sub_design <- function(designs,isin)
    lapply(designs,function(x) {
      attr(x,"expand") <- attr(x,"expand")[isin]
      x
    })


  model <- attr(dadm,"model")
  p_names <- attr(dadm,"p_names")
  sampled_p_names <- attr(dadm,"sampled_p_names")
  designs <- attr(dadm,"designs")
  expand <- attr(dadm,"expand")
  s_expand <- attr(dadm,"s_expand")
  unique_nort <- attr(dadm,"unique_nort")
  expand_nort <- attr(dadm,"expand_nort")
  unique_nortR <- attr(dadm,"unique_nortR")
  expand_nortR <- attr(dadm,"expand_nortR")
  ok_trials <- attr(dadm,"ok_trials")
  ok_dadm_winner <- attr(dadm,"ok_dadm_winner")
  ok_dadm_looser <- attr(dadm,"ok_dadm_looser")
  ok_da_winner <- attr(dadm,"ok_da_winner")
  ok_da_looser <- attr(dadm,"ok_da_looser")
  expand_uc <- attr(dadm,"expand_uc")
  expand_lc <- attr(dadm,"expand_lc")
  adapt <- attr(dadm,"adapt")

  # winner on expanded dadm
  expand_winner <- attr(dadm,"expand_winner")
  # subjects for first level of lR in expanded dadm
  slR1=dadm$subjects[expand][dadm$lR[expand]==levels(dadm$lR)[[1]]]

  dl <- stats::setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                        levels(dadm$subjects))
  for (i in levels(dadm$subjects)) {
    isin <- dadm$subjects==i         # dadm
    dl[[i]] <- dadm[isin,]
    dl[[i]]$subjects <- factor(as.character(dl[[i]]$subjects))
    if(is.null(attr(dadm, "custom_ll"))){
      isin1 <- s_expand==i             # da
      isin2 <- attr(dadm,"s_data")==i  # data


      attr(dl[[i]],"model") <- model
      attr(dl[[i]],"p_names") <- p_names
      attr(dl[[i]],"sampled_p_names") <- sampled_p_names
      attr(dl[[i]],"designs") <- sub_design(designs,isin)
      attr(dl[[i]],"expand") <- expand[isin1]-min(expand[isin1]) + 1
      attr(dl[[i]],"s_expand") <- NULL

      attr(dl[[i]],"ok_dadm_winner") <- ok_dadm_winner[isin]
      attr(dl[[i]],"ok_dadm_looser") <- ok_dadm_looser[isin]

      attr(dl[[i]],"ok_da_winner") <- ok_da_winner[isin1]
      attr(dl[[i]],"ok_da_looser") <- ok_da_looser[isin1]

      attr(dl[[i]],"unique_nort") <- unique_nort[isin]
      attr(dl[[i]],"unique_nortR") <- unique_nortR[isin]

      isinlR1 <- slR1==i
      if (!is.null(expand_nort)){
        attr(dl[[i]],"expand_nort") <-  expand_nort[isinlR1]- min( expand_nort[isinlR1]) + 1
      }

      if (!is.null(expand_nortR)){
        attr(dl[[i]],"expand_nortR") <- expand_nortR[isinlR1]-min(expand_nortR[isinlR1]) + 1
      }

      attr(dl[[i]],"ok_trials") <- ok_trials[isin2]
      if (!is.null(expand_winner)){
        attr(dl[[i]],"expand_winner") <- expand_winner[isin2]-min(expand_winner[isin2]) + 1
      }

      if (!is.null(attr(dadm,"expand_uc"))){
        attr(dl[[i]],"expand_uc") <- as.numeric(factor(expand_uc[isin2]))
      }
      if (!is.null(attr(dadm,"expand_lc"))){
        attr(dl[[i]],"expand_lc") <- as.numeric(factor(expand_lc[isin2]))
      }

      if (!is.null(attr(dadm,"LT"))){
        attr(dl[[i]],"LT") <- attr(dadm,"LT")[names(attr(dadm,"LT"))==i]
      }
      if (!is.null(attr(dadm,"UT"))){
        attr(dl[[i]],"UT") <- attr(dadm,"UT")[names(attr(dadm,"UT"))==i]
      }
      if (!is.null(attr(dadm,"LC"))){
        attr(dl[[i]],"LC") <- attr(dadm,"LC")[names(attr(dadm,"LC"))==i]
      }
      if (!is.null(attr(dadm,"UC"))){
        attr(dl[[i]],"UC") <- attr(dadm,"UC")[names(attr(dadm,"UC"))==i]
      }

      # adapt models
      if (!is.null(adapt)){
        attr(dl[[i]],"adapt") <- stats::setNames(list(adapt[[i]],adapt$design),c(i,"design"))
      }
    }
  }


  return(dl)
}



