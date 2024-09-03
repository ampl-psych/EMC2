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


#' Simulate data
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
#' @param mapped_p If `TRUE` instead returns a data frame with one row per design
#' cell and columns for each parameter specifying how they are mapped to the design cells.
#' @param hyper If `TRUE` the supplied parameters must be a set of samples, from which the group-level will be used to generate subject level parameters.
#' See also `make_random_effects` to generate subject-level parameters from a hyper distribution.
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

make_data <- function(parameters,design = NULL,n_trials=NULL,data=NULL,expand=1,
  mapped_p=FALSE, hyper = FALSE, ...)
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

  # #' @param Fcovariates either a data frame of covariate values with the same
  # #' number of rows as the data or a list of functions specifying covariates for
  # #' each trial. Must have names specified in the design Fcovariates argument.
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
  Fcovariates=NULL
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if(is(parameters, "emc")){
    if(is.null(design)) design <- attr(parameters, "design_list")[[1]]
    if(is.null(data)) data <- get_data(parameters)
    if(!hyper){
      parameters <- do.call(rbind, posterior_summary(parameters, probs = 0.5, selection = "alpha", by_subject = TRUE))
    } else{
      mu <- get_pars(parameters, selection = "mu", merge_chains = T, return_mcmc = F)
      Sigma <- get_pars(parameters, selection = "Sigma", merge_chains = T, return_mcmc = F)
      mu <- rowMeans(mu)
      Sigma <- apply(Sigma, 1:2, mean)
      parameters <- make_random_effects(design, group_means = mu, covariances = Sigma)
    }
  }

  sampled_p_names <- names(sampled_p_vector(design))
  if(is.null(dim(parameters))){
    if(is.null(names(parameters))) names(parameters) <- sampled_p_names
  } else{
    # design$Ffactors$subjects <- design$Ffactors$subjects[1:nrow(parameters)]
    # if(!is.null(data)){
    #   data<- data[data$subjects %in% design$Ffactors$subjects,]
    #   data$subjects <- factor(data$subjects)
    # }
    if(is.null(colnames(parameters))) colnames(parameters) <- sampled_p_names
    if(is.null(rownames(parameters))) rownames(parameters) <- design$Ffactors$subjects
  }

  if(!is.null(attr(design, "custom_ll"))){
    data <- list()
    for(i in 1:nrow(parameters)){
      data[[i]] <- attr(design, "rfun")(parameters[i,], n_trials = n_trials, subject = i)
    }
    return(do.call(rbind, data))
  }

  model <- design$model
  if (!is.matrix(parameters)) parameters <- make_pmat(parameters,design)
  if ( is.null(data) ) {
    design$Ffactors$subjects <- rownames(parameters)
    if (mapped_p) n_trials <- 1
    if ( is.null(n_trials) )
      stop("If data is not provided need to specify number of trials")
    Ffactors=c(design$Ffactors,list(trials=1:n_trials))
    data <- as.data.frame.table(array(dim=unlist(lapply(Ffactors,length)),
                                        dimnames=Ffactors))
    for (i in names(design$Ffactors))
      data[[i]] <- factor(data[[i]],levels=design$Ffactors[[i]])
    names(data)[dim(data)[2]] <- "R"
    data$R <- factor(data$R,levels=design$Rlevels)
    data$trials <- as.numeric(as.character(data$trials))
    # Add covariates
    if (!is.null(design$Fcovariates)) {
      if (!is.null(Fcovariates) & !all(unlist(lapply(Fcovariates,is.null)))) {
        if (!(all(names(Fcovariates)  %in% names(design$Fcovariates))))
          stop("All Fcovariates must be named in design$Fcovariates")
        if (!is.data.frame(Fcovariates) ) {
          if (!all(unlist(lapply(Fcovariates,is.function))))
            stop("Fcovariates must be either a data frame or list of functions")
          nams <- names(Fcovariates)
          Fcovariates <- do.call(cbind.data.frame,lapply(Fcovariates,function(x){x(data)}))
          names(Fcovariates) <- nams
        }
        n <- dim(Fcovariates)[1]
        if(n != nrow(data)) Fcovariates <- Fcovariates[sample(1:n, nrow(data), replace = TRUE),, drop = F]
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
    if (!is.function(model)) stop("model argument must  be a function")
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
    parameters[,attr(design,"ordinal")] <- exp(parameters[,attr(design,"ordinal")])
  pars <- model()$Ttransform(model()$Ntransform(map_p(
    model()$transform(add_constants(parameters,design$constants)),data
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
  attr(data,"p_vector") <- parameters;
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

#' Make random effects
#'
#' Simulates subject-level parameters in the format required by ``make_data()``.
#'
#' @param design A design list. The design as specified by `design()`
#' @param group_means A numeric vector. The group level means for each parameter, in the same order as `sampled_p_vector(design)`
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
#' subj_pars <- make_random_effects(design_DDMaE, group_means, n_subj = 5)
#'
#' # We can also define a covariance matrix to simulate from
#' subj_pars <- make_random_effects(design_DDMaE, group_means, n_subj = 5,
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
  if(length(group_means) != length(sampled_p_vector(design))) stop("You must specify as many means as parameters in your design")
  if(is.null(covariances)) covariances <- diag(abs(group_means)*variance_proportion)
  random_effects <- mvtnorm::rmvnorm(n_subj,mean=group_means,sigma=covariances)
  colnames(random_effects) <- names(sampled_p_vector(design))
  rownames(random_effects) <- subnames
  return(random_effects)
}


