#' make_missing
#'
#' Truncate or censor data. is.na(rt) not truncated or censored.
#'
#' @param data Data frame with rt and R columns
#' @param LT lower truncation bound below which data are removed (scalar or subject named vector)
#' @param UT upper truncation bound above which data are removed (scalar or subject named vector)
#' @param LC lower censoring bound (scalar or subject named vector)
#' @param UC upper censoring bound (scalar or subject named vector)
#' @param LCresponse Boolean, default TRUE, if false set LC response to NA
#' @param UCresponse Boolean, default TRUE, if false set UC response to NA
#' @param LCdirection Boolean, default TRUE, set LC rt to 0, else to NA
#' @param UCdirection Boolean, default TRUE, set LC rt to Inf, else to NA
#'
#' @return Truncated and censored data frame

make_missing <- function(data,LT=0,UT=Inf,LC=0,UC=Inf,
    LCresponse=TRUE,UCresponse=TRUE,LCdirection=TRUE,UCdirection=TRUE)
{

  censor <- function(data,L=0,U=Inf,Ld=TRUE,Ud=TRUE,Lr=TRUE,Ur=TRUE)
  {
    if (Ld) Ld <- -Inf else Ld <- NA
    if (Ud) Ud <- Inf else Ud <- NA
    snams <- levels(data$subjects)
    if (length(L)==1) L <- setNames(rep(L,length(snams)),snams)
    if (length(U)==1) U <- setNames(rep(U,length(snams)),snams)
    for (i in snams) {
      pick <- data$subjects==i & data$rt < L[i]
      pick[is.na(pick)] <- FALSE
      data$rt[pick] <- Ld
      if (!Lr) data$R[pick] <- NA
      pick <- data$subjects==i & data$rt > U[i]
      pick[is.na(pick)] <- FALSE
      data$rt[pick] <- Ud
      if (!Ur) data$R[pick] <- NA
    }
    data
  }

  pick <- is.infinite(data$rt) | (data$rt>LT & data$rt<UT)
  pick[is.na(pick)] <- TRUE
  out <- censor(data[pick,],L=LC,U=UC,Lr=LCresponse,Ur=UCresponse,Ld=LCdirection,Ud=UCdirection)
  if (LC != 0) attr(out,"LC") <- LC
  if (UC != Inf) attr(out,"UC") <- UC
  if (LT != 0) attr(out,"LT") <- LT
  if (UT != Inf) attr(out,"UT") <- UT
  out
}


#' Simulates data
#'
#' Simulates data based on design and parameter vector arguments using the
#' random function for a model (usually a design attribute, but if not can be
#' supplied separately) by one of two methods
#' 1) creating a fully crossed and balanced design specified by the design
#' argument's formula list (formuula) and factor contrast list (contrasts), if it is
#' null the data frame creation defaults are used) with number of trials per
#' cell specified by the trials argument ... or if the data argument is non-null
#' 2) using the design implicit in a data frame supplied, which allows creation
#' of unbalanced and other irregular designs, and replacing previous data with
#' simulated data.
#' Accuracy scored by design matchfun argument.
#'
#' @param p_vector parameter vector used to simulate data, or a matrix with one
#' row per subject (with corresponding row names) or a lsit of tables produced
#' by plot_density (in which case posterior medians are used to simulate data)
#' @param design design created by make_Design
#' @param model usually an attribute of design but if not must be given here.
#' @param trials number of trials per design cell
#' @param data If supplied that determines the design, with data sorted by subjects
#' and a trials = 1:n trials/subject factor added (or overwriting any existing).
#' Truncation and censoring information also taken from data attributes and presence
#' of 0/Inf and NA in rt column and NA in R column.
#' @param expand replicates the design expand times
#' @param mapped_p if true instead returns a data frame with one row per design
#' cell and columns for each parameter specifying how they are mapped to the
#' design cells.
#' @param LT lower truncation bound below which data are removed (scalar or subject named vector)
#' @param UT upper truncation bound above which data are removed (scalar or subject named vector)
#' @param LC lower censoring bound (scalar or subject named vector)
#' @param UC upper censoring bound (scalar or subject named vector)
#' @param LCresponse Boolean, default TRUE, if false set LC response to NA
#' @param UCresponse Boolean, default TRUE, if false set UC response to NA
#' @param LCdirection Boolean, default TRUE, set LC rt to 0, else to NA
#' @param UCdirection Boolean, default TRUE, set LC rt to Inf, else to NA
#' @param force_direction Boolean, take direction from argument not data (default FALSE)
#' @param force_response Boolean, take response from argument not data (default FALSE)
#' @param rtContaminantNA Boolean, TRUE sets contaminant trial rt to NA, if FALSE
#' (the default) direction is taken from data or LCdirection or UCdirection (NB
#' if both of these are false an error occurs as then contamination is not identifiable).
#' @param Fcovariates either a data frame of covariate values with the same
#' number of rows as the data or a list of functions specifying covariates for
#' each trial. Must have names specified in the design Fcovariates argument.
#' @param return_Ffunctions if false covariates are not returned
#'
#' @return None
#' @export

make_data <- function(p_vector,design,model=NULL,trials=NULL,data=NULL,expand=1,
  mapped_p=FALSE,Fcovariates=NULL,return_Ffunctions=FALSE,LT=0,UT=Inf,
  LC=0,UC=Inf,LCresponse=TRUE,UCresponse=TRUE,LCdirection=TRUE,UCdirection=TRUE,
  force_direction=FALSE,force_response=FALSE,rtContaminantNA=FALSE)
{

  if (is.list(p_vector))
    p_vector <- do.call(rbind,lapply(p_vector,function(x)x[2,]))
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
    # Add covariates
    if (!is.null(design$Fcovariates)) {
      if (!is.null(Fcovariates)) {
        if (!(all(names(Fcovariates)  %in% names(design$Fcovariates))))
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
      empty_covariates <- names(design$Fcovariates)[!(names(design$Fcovariates) %in% names(data))]
      if (length(empty_covariates)>0) data[,empty_covariates] <- 0
    }
  } else {
    LT <- attr(data,"LT"); if (is.null(LT)) LT <- 0
    UT <- attr(data,"UT"); if (is.null(UT)) UT <- Inf
    LC <- attr(data,"LC"); if (is.null(LC)) LC <- 0
    UC <- attr(data,"UC"); if (is.null(UC)) UC <- Inf
    if (!force_direction) {
      ok <- data$rt==-Inf; ok[is.na(ok)] <- FALSE
      LCdirection <- any(ok)
      ok <- data$rt==Inf; ok[is.na(ok)] <- FALSE
      UCdirection=any(ok)
    }
    if (!force_response) {
      if (!any(is.infinite(data$rt)) & any(is.na(data$R))) {
        LCresponse <- UCresponse <- FALSE
      } else {
        ok <- data$rt==-Inf
        bad <- is.na(ok)
        LCresponse <- !any(ok[!bad] & is.na(data$R[!bad]))
        ok <- data$rt==Inf
        bad <- is.na(ok)
        UCresponse <- !any(ok[!bad] & is.na(data$R[!bad]))
      }
    }
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
  if (!is.null(attr(design,"ordinal")))
    p_vector[,attr(design,"ordinal")] <- exp(p_vector[,attr(design,"ordinal")])
  pars <- model()$Ttransform(model()$Ntransform(map_p(
    model()$transform(add_constants(p_vector,design$constants)),data
  )),data)
  if ( any(dimnames(pars)[[2]]=="pContaminant") && any(pars[,"pContaminant"]>0) )
    pc <- pars[data$lR==levels(data$lR)[1],"pContaminant"] else pc <- NULL
  if (!is.null(design$adapt)) {
    if (expand>1) {
      expand <- 1
      warning("Expand does not work with this type of model")
    }
    # data <- adapt_data(data,design,model,pars,mapped_p=mapped_p,add_response = TRUE)
    if (mapped_p) return(data)
    adapt <- attr(data,"adapt")
    data <- data[data$lR==levels(data$lR)[1],!(names(data) %in% c("lR","lM"))]
    if('Qvalues' %in% names(attributes(pars))) attr(data, 'Qvalues') <- attr(pars, 'Qvalues')
    if('predictionErrors' %in% names(attributes(pars))) attr(data, 'predictionErrors') <- attr(pars, 'predictionErrors')
    attr(data,"adapt") <- adapt
    return(data)
  }
  if (mapped_p) return(cbind(data[,!(names(data) %in% c("R","rt"))],pars))
  if (expand>1) {
    data <- cbind(rep=rep(1:expand,each=dim(data)[1]),
                  data.frame(lapply(data,rep,times=expand)))
    lR <- rep(data$lR,expand)
    pars <- apply(pars,2,rep,times=expand)
  } else lR <- data$lR
  if (any(names(data)=="RACE")) {
    Rrt <- matrix(ncol=2,nrow=dim(data)[1]/length(levels(data$lR)),
                  dimnames=list(NULL,c("R","rt")))
    RACE <- data[data$lR==levels(data$lR)[1],"RACE"]
    ok <- as.numeric(data$lR) <= as.numeric(as.character(data$RACE))
    for (i in levels(RACE)) {
      pick <- data$RACE==i
      lRi <- factor(data$lR[pick & ok])
      Rrti <- model()$rfun(lRi,pars[pick & ok,])
      Rrti$R <- as.numeric(Rrti$R)
      Rrt[RACE==i,] <- as.matrix(Rrti)
    }
    Rrt <- data.frame(Rrt)
    Rrt$R <- factor(Rrt$R,labels=levels(lR))
  } else Rrt <- model()$rfun(lR,pars)
  dropNames <- c("lR","lM","lSmagnitude")
  if (!return_Ffunctions && !is.null(design$Ffunctions))
    dropNames <- c(dropNames,names(design$Ffunctions) )
  data <- data[data$lR==levels(data$lR)[1],!(names(data) %in% dropNames)]
  for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[,i]
  data <- make_missing(data[,names(data)!="winner"],LT,UT,LC,UC,
    LCresponse,UCresponse,LCdirection,UCdirection)
  if ( !is.null(pc) ) {
    if (!any(is.infinite(data$rt)) & any(is.na(data$R)))
      stop("Cannot have contamination and censoring with no direction and response")
    contam <- runif(length(pc)) < pc
    data[contam,"R"] <- NA
    if ( LC!=0 | is.finite(UC) ) { # censoring
      if ( (LCdirection & UCdirection) &  !rtContaminantNA)
        stop("Cannot have contamination with a mixture of censor directions")
      if (rtContaminantNA & ((is.finite(LC) & !LCresponse & !LCdirection) |
                              (is.finite(UC) & !UCresponse & !UCdirection)))
        stop("Cannot have contamination and censoring with no direction and response")
      if (rtContaminantNA | (!LCdirection & !UCdirection)) data[contam,"rt"] <- NA else
        if (LCdirection) data[contam,"rt"] <- -Inf  else data[contam,"rt"] <- Inf
    } else data[contam,"rt"] <- NA
  }
  attr(data,"p_vector") <- p_vector;
  data
}



add_Ffunctions <- function(data,design)
  # Adds columns created by Ffunctions (if not already there)
{
  Fdf <- data.frame(lapply(design$Ffunctions,function(f){f(data)}))
  ok <- !(names(Fdf) %in% names(data))
  if (!any(ok)) data else
    data <-  cbind.data.frame(data,Fdf[,ok,drop=FALSE])
}

#' Generate posterior predictives.
#'
#' Generates ``n_post`` data sets based on posterior parameter estimates
#'
#' @param samplers A list of samples from which you want to generate posterior predictives.
#' @param hyper Boolean. Default is FALSE, if true simulates from the group level instead of subject-level parameters. Only makes sense when aggregating across subjects
#' @param n_post Integer. How many data sets do you want to generate from the posterior.
#' @param filter Character. Choice of the stages 'preburn', 'burn', 'adapt', 'sample', default is 'sample'. From which stage do you want to take samples to generate new data with.
#' @param subfilter Integer or numeric vector. If integer, will filter out the first x of samples, within your filter. If numeric vector will select those samples, within your filter.
#' @param thin Integer. By how much do you want to thin the chains before simulating from them. Will keep 1/thin samples.
#' @param n_cores Integer. Across how many cores do you want to parallellize.
#' @param stat Character. Can be mean, median or default random. Will take either random samples from the chain or use the mean or median of the parameter estimates.
#' @param ... Optional additional arguments
#' @return A list of simulated data sets of length n_post.
#' @examples \dontrun{
#' # based on a set of samplers ran by run_emc we can generate posterior predictives
#' post_predict(samplers, n_cores = 8)
#' }
#' @export


post_predict <- function(samplers,hyper=FALSE,n_post=100,
                         filter="sample",subfilter=0,thin=1,n_cores=1,
                         stat=c("random","mean","median")[1], ...)
  # Post predictions for samples object, based on random samples or some
  # central tendency statistic.
  # n_post is number of parameter vectors used
  # expand=1 gives exact data design, larger values replicate whole design
  # filter/subfilter/thin as for as_mcmc.list
  # hyper=FALSE draws from alphas (participant level)
  # hyper=TRUE draws from hyper
{
  # #' @param force_direction Boolean, take censor direction from argument not samples (default FALSE)
  # #' @param force_response Boolean, take censor response from argument not samples (default FALSE)
  # #' @param LCresponse Boolean, default TRUE, if false set LC response to NA
  # #' @param UCresponse Boolean, default TRUE, if false set UC response to NA
  # #' @param LCdirection Boolean, default TRUE, set LC rt to 0, else to NA
  # #' @param UCdirection Boolean, default TRUE, set LC rt to Inf, else to NA
  # #' @param expand Integer. Default is 1, exact same design for each subject. Larger values will replicate designs, so more trials per subject.

  LCresponse<-TRUE; UCresponse<-TRUE; LCdirection<-TRUE; UCdirection<-TRUE
  force_direction <- FALSE; force_response <- FALSE;expand <- 1
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  data <- attr(samplers,"data_list")
  design <- attr(samplers,"design_list")
  model <- attr(samplers,"model_list")
  if(length(data) > 1){
    jointModel <- TRUE
    all_samples <- samplers
  } else{
    jointModel <- FALSE
  }
  post_out <- vector("list", length = length(data))
  for(j in 1:length(data)){
    if(jointModel) samplers <- single_out_joint(all_samples, j)
    subjects <- levels(data[[j]]$subjects)

    if (hyper) {
      pars <- vector(mode="list",length=n_post)
      for (i in 1:n_post) {
        pars[[i]] <- get_prior_samples(samplers,selection="alpha",
                                       filter=filter,thin=thin,subfilter=subfilter,n_prior=length(subjects))
        row.names(pars[[i]]) <- subjects
      }
    } else {
      samps <- lapply(as_mcmc.list(samplers,selection="alpha",
                                   filter=filter,subfilter=subfilter,thin=thin),function(x){do.call(rbind,x)})
      if (stat != "random") {
        p <- do.call(rbind,lapply(samps,function(x){apply(x,2,stat)}))
        row.names(p) <- subjects
      }
      pars <- vector(mode="list",length=n_post)
      for (i in 1:n_post) {
        if (stat != "random") pars[[i]] <- p else {
          pars[[i]] <- do.call(rbind,lapply(samps,function(x){x[sample(1:dim(x)[1],1),]}))
          row.names(pars[[i]]) <- subjects
        }
      }
    }
    if (n_cores==1) {
      simDat <- vector(mode="list",length=n_post)
      for (i in 1:n_post) {
        cat(".")
        simDat[[i]] <- make_data(pars[[i]],design=design[[j]],model=model[[j]],data=data[[j]],expand=expand,
          force_direction = force_direction,force_response=force_response,
          LCresponse=LCresponse,UCresponse=UCresponse,LCdirection=LCdirection,UCdirection=UCdirection)
      }
      cat("\n")
    } else {
      simDat <- mclapply(1:n_post,function(i){
        make_data(pars[[i]],design=design[[j]],model=model[[j]],data=data[[j]],expand=expand,
          force_direction = force_direction,force_response=force_response,
          LCresponse=LCresponse,UCresponse=UCresponse,LCdirection=LCdirection,UCdirection=UCdirection)
      },mc.cores=n_cores)
    }
    if (!is.null(attr(simDat[[1]],"adapt"))) adapt <- attr(simDat[[1]],"adapt")
    out <- cbind(postn=rep(1:n_post,times=unlist(lapply(simDat,function(x)dim(x)[1]))),do.call(rbind,simDat))
    if (!is.null(attr(simDat[[1]],"adapt"))) attr(out,"adapt") <- adapt
    if (n_post==1) pars <- pars[[1]]
    attr(out,"pars") <- pars
    post_out[[j]] <- out
  }
  if(!jointModel) post_out <- post_out[[1]]
  return(post_out)
}

#' Make random effects.
#'
#' @param design A design list. The design as specified by `make_design`
#' @param group_means A numeric vector. The group level means for each parameter, in the same order as `sampled_p_vector(design)`
#' @param n_subj An integer. The number of random effects to create.
#' @param variance_proportion A double. Optional. If covariances aren't specified, the variances will be created by multiplying the means by this number. The covariances will be 0.
#' @param covariances A covariance matrix. Optional. Specify the intended covariance matrix.
#'
#' @return A matrix of random effects.
#' @export

make_random_effects <- function(design, group_means, n_subj, variance_proportion = .2, covariances = NULL){
  if(length(group_means) != length(sampled_p_vector(design))) stop("You must specify as many means as parameters in your design")
  if(is.null(covariances)) covariances <- diag(abs(group_means)*variance_proportion)
  random_effects <- mvtnorm::rmvnorm(n_subj,mean=group_means,sigma=covariances)
  colnames(random_effects) <- names(sampled_p_vector(design))
  rownames(random_effects) <- as.character(1:n_subj)
  return(random_effects)
}


