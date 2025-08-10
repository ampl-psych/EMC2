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

# Simulate data by looping over trials
make_data_unconditional <- function(data, pars, design, model, return_covariates, return_trialwise_parameters) {
  trialwise_parameters <- covariates <- NULL

  ## loop over trials, save covariates
  if(!is.null(model()$trend)) {
    covariate_names = covariates = NULL
    for(trend in model()$trend) {
      if(trend$kernel %in% c('delta', 'deltab')) covariate_names <- c(covariate_names, trend$trend_pnames[2])
      if(trend$kernel == 'deltab') data$rt <- 0  # ensure RT is not NA, so it won't be ignored because it's set to NA
    }
    if(return_covariates) {
      covariates <- matrix(NA, nrow=nrow(data), ncol=length(covariate_names))
      colnames(covariates) <- covariate_names
    }
  }

  includeColumns <- !colnames(data) %in% c('lR', 'lM', names(design$Ffunctions))
  if(!'lR' %in% colnames(data)) data$lR <- factor(rep(1, nrow(data)))  # for simulations of the DDM, assume all rows are lR==1
  all_trials <- sort(unique(data[,'trials']))
  for(trial_idx in 1:length(all_trials)) {
    this_pars <- pars

    # simulate *two* trials at the same time, to capture update of Q-value from previous trial
    this_trial <- all_trials[trial_idx]

    if(trial_idx > 1) {
      prev_trial <- all_trials[trial_idx-1]
      # Remove pars from subjects that have fewer than this_trial trials
      this_pars <- this_pars[as.character(unique(this_data$subjects)),]
      # set Q-values to previous one
      this_pars[,covariate_names] <- this_covariates
    } else {
      prev_trial <- NULL
    }

    this_data <- design_model(
      add_accumulators(data[data$trials%in%c(prev_trial, this_trial)&data$lR==levels(data$lR)[1],includeColumns],
                       design$matchfun,simulate=FALSE,type=model()$type,Fcovariates=design$Fcovariates),
      design,model,add_acc=FALSE,compress=FALSE,verbose=FALSE,
      rt_check=FALSE)

    # drop unused levels
    this_data$subjects <- droplevels(this_data$subjects)
    # Ensure trials are sorted
    this_data <- this_data[order(this_data$subjects,this_data$trials),]

    # Map single trial's parameters
    this_pars <- map_p(this_pars, this_data, model())

    if(!is.null(model()$trend) && attr(model()$trend, "pretransform")){
      # This runs the trend and afterwards removes the trend parameters
      this_pars <- prep_trend(this_data, model()$trend, this_pars)
    }
    this_pars <- do_transform(this_pars, model()$transform)
    if(!is.null(model()$trend) && attr(model()$trend, "posttransform")){
      # This runs the trend and afterwards removes the trend parameters
      this_pars <- prep_trend(this_data, model()$trend, this_pars)
    }
    this_pars <- model()$Ttransform(this_pars, this_data)

    # drop previous trial from pars, data
    this_pars <- this_pars[this_data$trials==this_trial,]
    this_data <- this_data[this_data$trials==this_trial,]

    ## save parameters if requested
    if(return_trialwise_parameters) {
      if(is.null(trialwise_parameters)) {
        trialwise_parameters <- matrix(NA, nrow=nrow(data), ncol=ncol(this_pars))
        colnames(trialwise_parameters) <- colnames(this_pars)
      }
      trialwise_parameters[data$trials==this_trial,] <- this_pars
    }

    if(!is.null(model()$trend)) {
      # save covariates if requested
      if(return_covariates) {
        covariates[data$trials==this_trial,covariate_names] <- this_pars[,covariate_names]
      }
      this_covariates <- this_pars[this_data$lR==levels(this_data$lR)[1],covariate_names,drop=FALSE]
      this_pars <- this_pars[,!colnames(this_pars) %in% covariate_names]
    }

    this_pars <- add_bound(this_pars, model()$bound, this_data$lR)

    Rrt <- model()$rfun(this_data,this_pars)

    # drop lR
    if(!is.null(this_data$lR)) this_data <- this_data[this_data$lR == levels(this_data$lR)[1],!names(this_data) %in% c('lR', 'lM', 'winner')]
    for (i in dimnames(Rrt)[[2]]) this_data[[i]] <- Rrt[,i]

    ## re-apply Ffunctions to new data
    if(!is.null(design$Ffunctions)) for(i in names(design$Ffunctions)) this_data[,i] <- design$Ffunctions[[i]](this_data)

    # add to data
    match_idx <- with(data, data$trials == this_trial & data$subjects %in% this_data$subjects)
    data[match_idx, colnames(this_data)] <- this_data[match(data$subjects[match_idx], this_data$subjects), ]
  }

  # remove lR, lM
  data <- data[data$lR == levels(data$lR)[1],!names(data) %in% c('lR', 'lM')]
  return(list(data=data, covariates=covariates, trialwise_parameters=trialwise_parameters))
}

#' Simulate Data
#'
#' Simulates data based on a model design and a parameter vector (`p_vector`) by one of two methods:
#' 1) Creating a fully crossed and balanced design specified by the design,
#' with number of trials per cell specified by the `n_trials` argument
#' 2) Using the design of a data frame supplied, which allows creation
#' of unbalanced and other irregular designs, and replacing previous data with
#' simulated data
#'
#' To create data for multiple subjects see ``?make_random_effects()``.
#'
#' @param parameters parameter vector used to simulate data.
#' Can also be a matrix with one row per subject (with corresponding row names)
#' or an emc object with sampled parameters
#' (in which case posterior medians of `alpha` are used to simulate data)
#' @param design Design list created by ``design()``
#' @param n_trials Integer. If ``data`` is not supplied, number of trials to create per design cell
#' @param data Data frame. If supplied, the factors are taken from the data. Determines the number of trials per level of the design factors and can thus allow for unbalanced designs
#' @param expand Integer. Replicates the ``data`` (if supplied) expand times to increase number of trials per cell.
#' @param staircase Default NULL, used with stop-signal paradigm simulation to specify a staircase
#' algorithm. If non-null and a list then passed through as is, if not it is assigned the
#' default list structure: list(p=.25,SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf)
#' @param functions List of functions you want to apply to the data generation.
#' @param ... Additional optional arguments
#' @return A data frame with simulated data
#' @examples
#' # First create a design
#' design_DDMaE <- design(factors = list(S = c("left", "right"),
#'                                            E = c("SPD", "ACC"),
#'                                            subjects = 1:30),
#'                             Rlevels = c("left", "right"), model = DDM,
#'                             formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                             constants=c(s=log(1)))
#' # Then create a p_vector:
#' parameters <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_EACC=log(2), t0=log(.2),
#'               Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#'
#' # Now we can simulate data
#' data <- make_data(parameters, design_DDMaE, n_trials = 30)
#'
#' # We can also simulate data based on a specific dataset
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                             formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                             constants=c(s=log(1)))
#' parameters <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'               t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#'
#' data <- make_data(parameters, design_DDMaE, data = forstmann)
#' @export

make_data <- function(parameters,design = NULL,n_trials=NULL,data=NULL,expand=1, staircase = NULL,
                      functions = NULL, ...)
{
  # #' @param LT lower truncation bound below which data are removed (scalar or subject named vector)
  # #' @param UT upper truncation bound above which data are removed (scalar or subject named vector)
  # #' @param LC lower censoring bound (scalar or subject named vector)
  # #' @param UC upper censoring bound (scalar or subject named vector)
  # #' @param LCresponse Boolean, default TRUE, if false set LC response to NA
  # #' @param UCresponse Boolean, default TRUE, if false set UC response to NA
  # #' @param LCdirection Boolean, default TRUE, set LC rt to 0, else to NA
  # #' @param UCdirection Boolean, default TRUE, set LC rt to Inf, else to NA
  # #' @param force_direction Boolean, take direction from argument not data (default FALSE)
  # #' @param force_response Boolean, take response from argument not data (default FALSE)
  # #' @param rtContaminantNA Boolean, TRUE sets contaminant trial rt to NA, if FALSE
  # #' (the default) direction is taken from data or LCdirection or UCdirection (NB
  # #' if both of these are false an error occurs as then contamination is not identifiable).
  # #' @param return_Ffunctions if false covariates are not returned

  if (!is.null(staircase)){
    staircase <- check_staircase(staircase)
  }
  # #' @param Fcovariates either a data frame of covariate values with the same
  # #' number of rows as the data or a list of functions specifying covariates for
  # #' each trial. Must have names specified in the design Fcovariates argument.
  check_bounds <- FALSE

  LT<-0
  UT<-Inf
  LC<-0
  UC<-Inf
  LCresponse<-TRUE
  UCresponse<-TRUE
  LCdirection<-TRUE
  UCdirection<-TRUE
  force_direction<-FALSE
  force_response<-FALSE
  rtContaminantNA<-FALSE
  return_Ffunctions <- FALSE
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if(is(parameters, "emc")){
    if(is.null(design)) design <- get_design(parameters)[[1]] # Currently not supported for multiple designs
    if(is.null(data)) data <- get_data(parameters)
    parameters <- do.call(rbind, credint(parameters, probs = 0.5, selection = "alpha", by_subject = TRUE))
  }

  # Make sure parameters are in the right format, either matrix or vector
  sampled_p_names <- names(sampled_pars(design))
  if(is.null(dim(parameters))){
    if(is.null(names(parameters))) names(parameters) <- sampled_p_names
  } else{
    if(length(rownames(parameters)) != length(design$Ffactors$subjects)){
      stop("input parameter matrix must have number of rows equal to number of subjects specified in design")
    }
    if(is.null(colnames(parameters))) colnames(parameters) <- sampled_p_names
    rownames(parameters) <- design$Ffactors$subjects
  }

  if(!is.null(attr(design, "custom_ll"))){
    data <- list()
    for(i in 1:nrow(parameters)){
      data[[i]] <- design$model()$rfun(parameters[i,], n_trials = n_trials, subject = i)
    }
    return(do.call(rbind, data))
  }

  model <- design$model

  if(grepl("MRI", model()$type)){
    return(make_data_wrapper_MRI(parameters, data, design))
  }
  if(is.data.frame(parameters)) parameters <- as.matrix(parameters)
  if (!is.matrix(parameters)) parameters <- make_pmat(parameters,design)
  if ( is.null(data) ) {
    design$Ffactors$subjects <- rownames(parameters)
    if ( is.null(n_trials) )
      stop("If data is not provided need to specify number of trials")
    design_in <- design
    design_in$Fcovariates <- design_in$Fcovariates[!design$Fcovariates %in% names(functions)]
    data <- minimal_design(design_in, covariates = list(...)$covariates,
                             drop_subjects = F, n_trials = n_trials, add_acc=F,
                           drop_R = F)
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
  if(!is.null(functions)){
    for(i in 1:length(functions)){
      data[[names(functions)[i]]] <- functions[[i]](data)
    }
  }
  if (!is.factor(data$subjects)) data$subjects <- factor(data$subjects)
  if (!is.null(model)) {
    if (!is.function(model)) stop("model argument must  be a function")
    if ( is.null(model()$p_types) ) stop("model()$p_types must be specified")
    if ( is.null(model()$Ttransform) ) stop("model()$Ttransform must be specified")
  }

  data <- design_model(
    add_accumulators(data,design$matchfun,simulate=TRUE,type=model()$type,Fcovariates=design$Fcovariates),
    design,model,add_acc=FALSE,compress=FALSE,verbose=FALSE,
    rt_check=FALSE)

  simulate_unconditional_on_data <- return_covariates <- return_trialwise_parameters <- FALSE
  if('conditional_on_data' %in% names(list(...))) {
    simulate_unconditional_on_data <- !list(...)$conditional_on_data
  }
  if('return_covariates' %in% names(list(...))) {
    return_covariates <- list(...)$return_covariates
  }
  if('return_trialwise_parameters' %in% names(list(...))) {
    return_trialwise_parameters <- list(...)$return_trialwise_parameters
    # if(return_trialwise_parameters&!simulate_unconditional_on_data) stop('Cannot return trialwise parameters when simulating conditional on data')
  }

  ## For both conditional and unconditional simulations...
  pars <- t(apply(parameters, 1, do_pre_transform, model()$pre_transform))
  pars <- add_constants(pars,design$constants)

  if(simulate_unconditional_on_data) {
    out <- make_data_unconditional(data=data, pars=pars, design=design, model=model,
                                   return_covariates=return_covariates,
                                   return_trialwise_parameters=return_trialwise_parameters)
    data <- out$data
    covariates <- out$covariates
    trialwise_parameters <- out$trialwise_parameters

    # add contamination
    if ( any(dimnames(pars)[[2]]=="pContaminant") && any(pars[,"pContaminant"]>0) )
      pc <- pars[data$lR==levels(data$lR)[1],"pContaminant"] else pc <- NULL

  } else {
    ## Default, vectorized way of simulating data
    pars <- map_p(pars, data, model())  # constants were already added above
    if(!is.null(model()$trend) && attr(model()$trend, "pretransform")){
      # This runs the trend and afterwards removes the trend parameters
      pars <- prep_trend(data, model()$trend, pars)
    }
    pars <- do_transform(pars, model()$transform)
    if(!is.null(model()$trend) && attr(model()$trend, "posttransform")){
      # This runs the trend and afterwards removes the trend parameters
      pars <- prep_trend(data, model()$trend, pars)
    }
    pars <- model()$Ttransform(pars, data)

    ## collect covariates and remove from pars
    if(!is.null(model()$trend)) {
      covariate_names = covariates = NULL
      for(trend in model()$trend) {
        if(trend$kernel == 'delta') covariate_names <- c(covariate_names, trend$trend_pnames[2])
      }
      if(return_covariates) {
        covariates <- pars[,covariate_names]
      }
      pars <- pars[,!colnames(pars) %in% covariate_names]
    }

    pars <- add_bound(pars, model()$bound, data$lR)

    ## SM: maybe of interest
    if(return_trialwise_parameters) trialwise_parameters <- pars
    pars_ok <- attr(pars, 'ok')
    if(mean(!pars_ok) > .1){
      warning("More than 10% of parameter values fall out of model bounds, see <model_name>$bounds()")
      return(FALSE)
    }
    if ( any(dimnames(pars)[[2]]=="pContaminant") && any(pars[,"pContaminant"]>0) )
      pc <- pars[data$lR==levels(data$lR)[1],"pContaminant"] else pc <- NULL
    if (expand>1) {
      data <- cbind(rep=rep(1:expand,each=dim(data)[1]),
                    data.frame(lapply(data,rep,times=expand)))
      pars <- apply(pars,2,rep,times=expand)
    }
    if (!is.null(staircase)) {
      attr(data, "staircase") <- staircase
    }
    if (any(names(data)=="RACE")) {
      Rrt <- RACE_rfun(data, pars, model)
    } else Rrt <- model()$rfun(data,pars)
    dropNames <- c("lR","lM","lSmagnitude")
    if (!return_Ffunctions && !is.null(design$Ffunctions))
      dropNames <- c(dropNames,names(design$Ffunctions))
    if(!is.null(data$lR)) data <- data[data$lR == levels(data$lR)[1],]
    data <- data[,!(names(data) %in% dropNames)]
    for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[,i]
  }

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
  attr(data,"p_vector") <- parameters;
  if(return_covariates) attr(data,'covariates') <- covariates
  if(return_trialwise_parameters) attr(data, 'trialwise_parameters') <- trialwise_parameters
  data
}

RACE_rfun <- function(data, pars, model){
  Rrt <- matrix(ncol=2,nrow=dim(data)[1]/length(levels(data$lR)),
         dimnames=list(NULL,c("R","rt")))
  RACE <- data[data$lR==levels(data$lR)[1],"RACE"]
  ok <- as.numeric(data$lR) <= as.numeric(as.character(data$RACE))
  for (i in levels(RACE)) {
    pick <- data$RACE==i
    data_in <- data[pick & ok,]
    data_in$lR <- factor(data$lR[pick & ok])
    tmp <- pars[pick & ok,]
    attr(tmp, "ok") <- rep(T, nrow(tmp))
    Rrti <- model()$rfun(data_in,tmp)
    Rrti$R <- as.numeric(Rrti$R)
    Rrt[RACE==i,] <- as.matrix(Rrti)
  }
  Rrt <- data.frame(Rrt)
  Rrt$R <- factor(Rrt$R, labels = levels(data$lR), levels = 1:length(levels(data$lR)))
  return(Rrt)
}

add_Ffunctions <- function(data,design)
  # Adds columns created by Ffunctions (if not already there)
{
  Fdf <- data.frame(lapply(design$Ffunctions,function(f){f(data)}))
  ok <- !(names(Fdf) %in% names(data))
  if (!any(ok)) data else
    data <-  cbind.data.frame(data,Fdf[,ok,drop=FALSE])
}

#' Generate Subject-Level Parameters
#'
#' Simulates subject-level parameters in the format required by ``make_data()``.
#'
#' @param design A design list. The design as specified by `design()`
#' @param group_means A numeric vector. The group level means for each parameter, in the same order as `sampled_pars(design)`
#' @param n_subj An integer. The number of subjects to generate parameters for. If `NULL` will be inferred from design
#' @param variance_proportion A double. Optional. If ``covariances`` are not specified, the variances will be created by multiplying the means by this number. The covariances will be 0.
#' @param covariances A covariance matrix. Optional. Specify the intended covariance matrix.
#'
#' @return A matrix of subject-level parameters.
#' @examples
#' # First create a design
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                             formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                             constants=c(s=log(1)))
#' # Then create a group-level means vector:
#' group_means =c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'                t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' # Now we can create subject-level parameters
#' subj_pars <- make_random_effects(design_DDMaE, group_means, n_subj = 19)
#'
#' # We can also define a covariance matrix to simulate from
#' subj_pars <- make_random_effects(design_DDMaE, group_means, n_subj = 19,
#'              covariances = diag(.1, length(group_means)))
#'
#' # The subject level parameters can be used to generate data
#' make_data(subj_pars, design_DDMaE, n_trials = 10)
#' @export

make_random_effects <- function(design, group_means, n_subj = NULL, variance_proportion = .2, covariances = NULL){
  if(is.null(n_subj)){
    n_subj <- length(design$Ffactors$subjects)
    subnames <- design$Ffactors$subjects
  } else{
    subnames <- as.character(1:n_subj)
  }
  if(length(group_means) != length(sampled_pars(design))) stop("You must specify as many means as parameters in your design")
  if(is.null(covariances)) covariances <- diag(abs(group_means)*variance_proportion)
  random_effects <- mvtnorm::rmvnorm(n_subj,mean=group_means,sigma=covariances)
  colnames(random_effects) <- names(sampled_pars(design))
  rownames(random_effects) <- subnames
  return(random_effects)
}


