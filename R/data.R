#' Simulates data
#'
#' Simulates data based on design and parameter vector arguments using the
#' random function for a model (usually a design attribute, but if not can be
#' supplied separately) by one of two methods
#' 1) creating a fully crossed and balanced design specified by the design
#' argument's formula list (Flist) and factor contrast list (Clist, if it is
#' null the data frame creation defaults are used) with number of trials per
#' cell specified by the trials argument ... or if the data argument is non-null
#' 2) using the design implicit in a data frame supplied, which allows creation
#' # of unbalanced and other irregular designs, and replacing previous data with
#' simulated data.
#' Accuracy scored by design matchfun argument.
#'
#' @param p_vector parameter vector
#' @param design design created by make_Design
#' @param model usually an attribute of design but if not must be given here.
#' @param trials number of trials per design cell
#' @param data A supplied that determines the design, with data sorted by subjects
#' and a trials = 1:n trials/subject factor added (or overwriting any existing)
#' @param expand replicates the design expand times
#' @param mapped_p if true instead returns a data frame with one row per design
#' cell and columns for each parameter specifying how they are mapped to the
#' design cells.
#' @param LT lower truncation bound below which data are removed
#' @param UT upper truncation bound above which data are removed
#' @param LC lower censoring bound below which data are turned into NA
#' @param UC upper censoring bound above which data are turned into NA
#' @param Fcovariates either a data frame of covariate valeus with the same
#' number of rows as the data or a list of functions specifying covariates for
#' each trial. Must have names specified in the design Fcovariates argument.
#' @param return_Ffunctions if false covariates are not returned
#'
#' @return
#' @export
#'
#' @examples
make_data <- function(p_vector,design,model=NULL,trials=NULL,data=NULL,expand=1,
                      mapped_p=FALSE,LT=NULL,UT=NULL,LC=NULL,UC=NULL,
                      Fcovariates=NULL,return_Ffunctions=FALSE)
{

  missingFilter <- function(data,LT,UT,LC,UC,Rmissing)
  {

    makeExclude <- function(data,L=NULL,U=NULL)
    {

      makeCrit <- function(data,bound) {
        exclude <- factor(data$subjects)
        levels(exclude) <- bound[levels(exclude)]
        as.numeric(as.character(exclude))
      }

      if (!is.null(L)) {
        if (length(L)==1) exclude <- data$rt < L else {
          if (!all(sort(levels(data$subjects))==names(L)))
            stop("Missing subjects")
          exclude <- data$rt < makeCrit(exclude)
        }
      } else {
        if (length(U)==1) exclude <- data$rt > U else {
          if (!all(sort(levels(data$subjects))==names(U)))
            stop("Missing subjects")
          exclude <- data$rt > makeCrit(exclude)
        }
      }
      exclude
    }

    if (!is.null(LT)) data <- data[!makeExclude(data,L=LT),]
    if (!is.null(UT)) data <- data[!makeExclude(data,U=UT),]
    if (!is.null(LT)) {
      exclude <- makeExclude(data,L=LC)
      data$rt[exclude] <- -Inf
      if (Rmissing) data$R[exclude] <- NA
    }
    if (!is.null(UT)) {
      exclude <- makeExclude(data,U=UC)
      data$rt[exclude] <- -Inf
      if (Rmissing) data$R[exclude] <- NA
    }
    data
  }

  if (is.null(model)) if (is.null(design$model))
    stop("Must specify model as not in design") else model <- design$model
    if (!is.matrix(p_vector)) p_vector <- make_pmat(p_vector,design)
    if ( is.null(data) ) {
      if (mapped_p) trials <- 1
      if ( is.null(trials) )
        stop("If data is not provided need to specify number of trials")
      Ffactors=c(design$Ffactors,list(trials=1:trials))
      data <- as.data.frame.table(array(dim=unlist(lapply(Ffactors,length)),
                                        dimnames=Ffactors))
      for (i in names(design$Ffactors))
        data[[i]] <- factor(data[[i]],levels=design$Ffactors[[i]])
      names(data)[dim(data)[2]] <- "R"
      data$R <- factor(data$R,levels=design$Rlevels)
      data$trials <- as.numeric(as.character(data$trials))
      # if (!is.null(design$Ffunctions)) data <-
      #   cbind.data.frame(data,data.frame(lapply(design$Ffunctions,function(f){f(data)})))
      LT <- UT <- LC <- UC <- Rmissing <- NULL
      # Add covariates
      if (!is.null(design$Fcovariates)) {
        if (!is.null(Fcovariates)) {
          if (!(all(names(Fcovariates)  %in% design$Fcovariates)))
            stop("All Fcovariates must be named in design$Fcovariates")
          if (!is.data.frame(Fcovariates)) {
            if (!all(unlist(lapply(Fcovariates,is.function))))
              stop("Fcovariates must be either a data frame or list of functions")
            nams <- names(Fcovariates)
            Fcovariates <- do.call(cbind.data.frame,lapply(Fcovariates,function(x){x(data)}))
            names(Fcovariates) <- nams
          }
          n <- dim(Fcovariates)[1]
          if (!(n==dim(data)[1])) stop("Fcovariates must specify ",dim(data)[1]," values per covariate")
          data <- cbind.data.frame(data,Fcovariates)
        }
        empty_covariates <- design$Fcovariates[!(design$Fcovariates %in% names(Fcovariates))]
        if (length(empty_covariates)>0) data[[empty_covariates]] <- 0
      }
    } else {
      LT <- attr(data,"LT"); UT <- attr(data,"UT")
      LC <- attr(data,"LC"); UC <- attr(data,"UC")
      Rmissing <- any(is.na(data$R))
      data <- add_trials(data[order(data$subjects),])
    }
    if (!is.factor(data$subjects)) data$subjects <- factor(data$subjects)
    if (!is.null(model)) {
      if (!is.function(model)) stop("model arguement must  be a function")
      if ( is.null(model()$p_types) ) stop("model()$p_types must be specified")
      if ( is.null(model()$transform) ) stop("model()$transform must be specified")
      if ( is.null(model()$Ntransform) ) stop("model()$Ntransform must be specified")
      if ( is.null(model()$Ttransform) ) stop("model()$Ttransform must be specified")
    }
    data <- design_model(
      add_accumulators(data,design$matchfun,simulate=TRUE,type=model()$type,Fcovariates=design$Fcovariates),
      design,model,add_acc=FALSE,compress=FALSE,verbose=FALSE,
      rt_check=FALSE)
    pars <- model()$Ttransform(model()$Ntransform(map_p(
      model()$transform(add_constants(p_vector,design$constants)),data
    )),data)
    if (!is.null(design$adapt)) {
      if (expand>1) {
        expand <- 1
        warning("Expand does not work with this type of model")
      }
      data <- adapt_data(data,design,model,pars,mapped_p=mapped_p,add_response = TRUE)
      if (mapped_p) return(data)
      adapt <- attr(data,"adapt")
      data <- data[data$lR==levels(data$lR)[1],!(names(data) %in% c("lR","lM"))]
      attr(data,"adapt") <- adapt
      return(data)
    }
    if (mapped_p) return(cbind(data[,!(names(data) %in% c("R","rt"))],pars))
    if (expand==1)
      Rrt <- model()$rfun(data$lR,pars) else
        Rrt <- model()$rfun(rep(data$lR,expand),apply(pars,2,rep,times=expand))
    if (expand>1) data <- cbind(rep=rep(1:expand,each=dim(data)[1]),
                                data.frame(lapply(data,rep,times=expand)))
    dropNames <- c("lR","lM")
    if (!return_Ffunctions && !is.null(design$Ffunctions))
      dropNames <- c(dropNames,names(design$Ffunctions) )
    data <- data[data$lR==levels(data$lR)[1],!(names(data) %in% dropNames)]
    for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[,i]
    missingFilter(data[,names(data)!="winner"],LT,UT,LC,UC,Rmissing)
}

add_Ffunctions <- function(data,design)
  # Adds columns created by Ffunctions (if not already there)
{
  Fdf <- data.frame(lapply(design$Ffunctions,function(f){f(data)}))
  ok <- !(names(Fdf) %in% names(data))
  if (!any(ok)) data else
    data <-  cbind.data.frame(data,Fdf[,ok,drop=FALSE])
}

post_predict <- function(samples,hyper=FALSE,n_post=100,expand=1,
                         filter="sample",subfilter=0,thin=1,n_cores=1,
                         use_par=c("random","mean","median")[1])
  # Post predictions for samples object, based on random samples or some
  # central tendency statistic.
  # n_post is number of parameter vectors used
  # expand=1 gives exact data design, larger values replicate whole design
  # filter/subfilter/thin as for as_mcmc.list
  # hyper=FALSE draws from alphas (participant level)
  # hyper=TRUE draws from hyper
{

  data <- attr(samples,"data_list")
  design <- attr(samples,"design_list")
  model <- attr(samples,"model_list")
  if(length(data) > 1){
    jointModel <- TRUE
    all_samples <- samples
  } else{
    jointModel <- FALSE
  }
  post_out <- vector("list", length = length(data))
  for(j in 1:length(data)){
    if(jointModel) samples <- single_out_joint(all_samples, j)
    subjects <- levels(data[[j]]$subjects)

    if (hyper) {
      pars <- vector(mode="list",length=n_post)
      for (i in 1:n_post) {
        pars[[i]] <- get_prior_samples(samples,selection="alpha",
                                       filter=filter,thin=thin,subfilter=subfilter,n_prior=length(subjects))
        row.names(pars[[i]]) <- subjects
      }
    } else {
      samps <- lapply(as_mcmc.list(samples,selection="alpha",
                                   filter=filter,subfilter=subfilter,thin=thin),function(x){do.call(rbind,x)})
      if (use_par != "random") {
        p <- do.call(rbind,lapply(samps,function(x){apply(x,2,use_par)}))
        row.names(p) <- subjects
      }
      pars <- vector(mode="list",length=n_post)
      for (i in 1:n_post) {
        if (use_par != "random") pars[[i]] <- p else {
          pars[[i]] <- do.call(rbind,lapply(samps,function(x){x[sample(1:dim(x)[1],1),]}))
          row.names(pars[[i]]) <- subjects
        }
      }
    }
    if (n_cores==1) {
      simDat <- vector(mode="list",length=n_post)
      for (i in 1:n_post) {
        cat(".")
        simDat[[i]] <- make_data(pars[[i]],design=design[[j]],model=model[[j]],data=data[[j]],expand=expand)
      }
      cat("\n")
    } else {
      simDat <- mclapply(1:n_post,function(i){
        make_data(pars[[i]],design=design[[j]],model=model[[j]],data=data[[j]],expand=expand)
      },mc.cores=n_cores)
    }
    if (!is.null(attr(simDat[[1]],"adapt"))) adapt <- attr(simDat[[1]],"adapt")
    out <- cbind(postn=rep(1:n_post,each=dim(simDat[[1]])[1]),do.call(rbind,simDat))
    if (!is.null(attr(simDat[[1]],"adapt"))) attr(out,"adapt") <- adapt
    if (n_post==1) pars <- pars[[1]]
    attr(out,"pars") <- pars
    post_out[[j]] <- out
  }
  if(!jointModel) post_out <- post_out[[1]]
  return(post_out)
}


