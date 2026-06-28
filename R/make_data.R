get_missing <- function(supplied, data, bound_name, default,type) {

  subjectwise <- FALSE
  if (is.function(supplied)) supplied <- supplied(data)
  if (is.null(supplied)) {
    if (bound_name %in% colnames(data)) bound <- data[[bound_name]] else
      bound <- default
  } else {
    if (type=="logical" & !is.logical(supplied))
      stop(bound_name," must be logical")
    if (type=="numeric" & !is.numeric(supplied))
      stop(bound_name," must be numeric")
    subjectwise <- all(utils::hasName(supplied,levels(data$subjects)))
    if (!subjectwise) bound <- supplied else
      bound <- supplied[as.character(data$subjects)]
  }
  out <- rep(bound,length.out=nrow(data))
  attr(out,"subjectwise") <- subjectwise
  out
}

#' Add information about missing values to data and modify/filter accordingly.
#'
#' Columns corresponding to LC, UT, LC, UC, and pContaminant arguments are added
#' to the return, specifying, respectively, if a row is subject to lower or upper
#' truncation, lower or upper censoring, or contamination with the given probability.
#' Truncation removes rows. Contamination makes R = NA and the rt = Inf, and is
#' applied after truncation but before censoring. Lower/upper censoring sets
#'   R = NA if LCresponse/UCresponse = TRUE and
#'   rt = -Inf/Inf if LCdireciton/UCdireciton = TRUE (else NA).
#' All of these arguments can be: 1) NULL, in which case if the data frame has a column
#' of that name it is used, or if not the default is used (see argument definition),
#' 2) a scalar/logical (same value for every data row), 3) a subject-named  vector
#' (same value for each subject), 4) a vector with length matching data rows, or
#' 5) a function taking data as it argument that creates a column of appropriate
#' values. Note that if this function returns a subject named vector it will be
#' expanded as in case (3). For rows with response (column R) "nogo" the following
#' are enforced: UC=LC=LT=0, UT=Inf, UCdirection=UCresponse=TRUE
#'
#' @param data Data frame to be modified
#' @param LT Lower truncation bound below which data are removed, default 0.
#' @param UT Upper truncation bound above which data are removed, default Inf.
#' @param LC Lower censoring bound, default 0.
#' @param UC Upper censoring bound, default Inf.
#' @param LCresponse Logical. Default FALSE, set responses to NA on lower-censored trials.
#' @param UCresponse Logical. Default FALSE, set responses to NA on upper-censored trials.
#' @param LCdirection Logical. Default TRUE, lower-censored RTs are coded as -Inf; if FALSE, as NA.
#' @param UCdirection Logical. Default TRUE, upper-censored RTs are coded as Inf; if FALSE, as NA.
#' @param pContaminant Probability of contamination, default 0.
#' @param no_truncate Logical, default FALSE, for TRUE don't apply truncation to row (except if GO/NOGO).
#' @param no_censor Logical, default FALSE, for TRUE don't apply censor to row (except if GO/NOGO).
#' @param verbose Logical. Default FALSE, if TRUE report effects of filtering.
#' @param rt_resolution A double, see make_emc, specified here so binning of rt and LC/UC/LT/UT is consistent.
#'        The default is 1/60 as in make_emc, but when make_missing is called by make_data the default is to
#'        do nothing unless an explicit is value passed in the missing list.
#' @param digits Integer, number of decimal places used when rounding the
#'        percentage summaries reported when \code{verbose = TRUE}. Default 2.
#' @return A filtered and modified data frame with added/updated LC, UC, LT and UT columns
#' @examples
#' # First make some data
#' designRDM <- design(model = RDM,
#'   factors = list(subjects = 1:2, S = c("left", "right")),Rlevels = c("left", "right"),
#'   matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
#'   formula = list(B ~ 1, v ~ lM, A ~ 1, t0 ~ 1, s ~ lM),
#'   constants = c(s = log(1)))
#' p_vector <- log(c(B=2,A=.5,t0=0.2,v=1,v_lMTRUE=2,s_lMTRUE=.8))
#' dat <- make_data(p_vector, designRDM,n_trials = 10)
#'
#' # Filter data frame without LT/UC/LT/UT columns (as in most real data files)
#' data <- dat
#' mdata <- make_missing(dat,LT=.7,LC=.75,UC=1.5,UT=1.6,verbose=TRUE)
#' @export

make_missing <- function(data, LT = NULL, UT = NULL, LC = NULL, UC = NULL,
                         LCresponse = NULL, UCresponse = NULL,LCdirection = NULL, UCdirection = NULL,
                         no_truncate=FALSE,no_censor=FALSE,
                         pContaminant=NULL,verbose=FALSE,rt_resolution=1/60,digits = 2)
{

  no_truncate <- get_missing(no_truncate, data, "no_truncate",FALSE,"logical")
  LT <- get_missing(LT, data, "LT",0,"numeric")
  data$LT[!no_truncate] <- as.numeric(LT[!no_truncate])
  UT <- get_missing(UT, data, "UT",Inf,"numeric")
  data$UT[!no_truncate] <- as.numeric(UT[!no_truncate])

  # Go/no-go data are defined by explicit non-responses. Truncation is not a
  # coherent missing-data mechanism in this setting, so disable LT/UT globally.
  is_gng_data <- FALSE
  if ("R" %in% names(data)) {
    if (is.factor(data$R)) is_gng_data <- "nogo" %in% levels(data$R)
    if (!is_gng_data) is_gng_data <- any(data$R == "nogo", na.rm = TRUE)
  }
  if (is_gng_data) {
    requested_trunc <- any(((data$LT != 0) | is.finite(data$UT)) & !no_truncate, na.rm = TRUE)
    if (requested_trunc) {
      warning("Ignoring LT/UT truncation for go/no-go data (R includes 'nogo').")
    }
    no_truncate[] <- TRUE
    data$LT[] <- 0
    data$UT[] <- Inf
  }


  no_censor <- get_missing(no_censor, data, "no_censor",FALSE,"logical")
  LC <- get_missing(LC, data, "LC",0,"numeric")
  data$LC[!no_censor] <- as.numeric(LC[!no_censor])
  UC <- get_missing(UC, data, "UC",Inf,"numeric")
  data$UC[!no_censor] <- as.numeric(UC[!no_censor])

  if (!is.null(rt_resolution)) {
    data$rt <- .floor_to_rt_resolution(data$rt, rt_resolution)
    data$LC <- .floor_to_rt_resolution(data$LC, rt_resolution)
    data$UC <- .floor_to_rt_resolution(data$UC, rt_resolution)
    data$LT <- .floor_to_rt_resolution(data$LT, rt_resolution)
    data$UT <- .floor_to_rt_resolution(data$UT, rt_resolution)
  }

  LCresponse <- get_missing(LCresponse, data, "LCresponse",FALSE,"logical")
  UCresponse <- get_missing(UCresponse, data, "UCresponse",FALSE,"logical")
  LCdirection <- get_missing(LCdirection, data, "LCdirection",TRUE,"logical")
  UCdirection <- get_missing(UCdirection, data, "UCdirection",TRUE,"logical")

  LT_eff <- data$LT
  LC_eff <- data$LC
  UT_eff <- data$UT
  UC_eff <- data$UC

  isgng <- data$R == "nogo"
  isgng[is.na(isgng)] <- FALSE
  if (any(isgng)) {
    UC_eff[isgng] <- 0
    LC_eff[isgng] <- 0
    UT_eff[isgng] <- Inf
    LT_eff[isgng] <- 0
    UCdirection[isgng] <- TRUE
    UCresponse[isgng] <- FALSE
  }

  tol_l <- sqrt(.Machine$double.eps) * pmax(1, abs(LT_eff), abs(LC_eff))
  tol_u <- sqrt(.Machine$double.eps) * pmax(1, abs(UC_eff), abs(UT_eff))
  if (any((LT_eff - LC_eff) > tol_l & LC_eff != 0, na.rm = TRUE)) stop("LT > LC not allowed")
  if (any((UC_eff - UT_eff) > tol_u & UC_eff != Inf, na.rm = TRUE)) stop("UC > UT not allowed")

  # Only keep trials in LT-UT (inclusive) or infinite or NA
  cutL <- is.finite(data$rt) & (data$rt < LT_eff & is.finite(data$rt))
  cutL[is.na(cutL)] <- FALSE; cutL[no_truncate] <- FALSE
  cutU <- (data$rt > UT_eff & is.finite(data$rt))
  cutU[is.na(cutU)] <- FALSE; cutU[no_truncate] <- FALSE
  if (verbose) {
    if (!all(LT_eff==0)) {
      if (!attr(LT,"subjectwise")) stat <- mean(cutL) else
        stat <- tapply(cutL,data$subjects,mean)
      message("% lower truncation")
      print(round(100*stat,digits))
    }
    if (!all(UT_eff==Inf)) {
      if (!attr(UT,"subjectwise")) stat <- mean(cutU) else
        stat <- tapply(cutU,data$subjects,mean)
      message("% upper truncation")
      print(round(100*stat,digits))
    }
  }

  # Truncate
  data <- data[!cutL & !cutU, ]
  LT_eff <- LT_eff[!cutL & !cutU]
  LC_eff <- LC_eff[!cutL & !cutU]
  UT_eff <- UT_eff[!cutL & !cutU]
  UC_eff <- UC_eff[!cutL & !cutU]
  LCresponse <- LCresponse[!cutL & !cutU]
  UCresponse <- UCresponse[!cutL & !cutU]
  LCdirection <- LCdirection[!cutL & !cutU]
  UCdirection <- UCdirection[!cutL & !cutU]
  no_censor <- no_censor[!cutL & !cutU]

  pContaminant <- get_missing(pContaminant, data, "pContaminant",0,"numeric")
  if (!all(pContaminant==0)) {
    contam <- rbinom(nrow(data), 1, pContaminant) == 1
    data[contam, "rt"] <- Inf
    data[contam, "R"] <- NA
    if (verbose) {
      if (!attr(pContaminant,"subjectwise")) stat <- mean(contam) else
        stat <- tapply(contam,data$subjects,mean)
      message("% contaminated")
      print(round(100*stat,digits))
    }
  }

  # Censoring proportions (like truncation dont censor if equal to LC or UC)
  cutL <- is.finite(data$rt) & (data$rt < LC_eff)
  cutL[is.na(cutL)] <- TRUE; cutL[no_censor] <- FALSE
  cutU <- (data$rt > UC_eff)
  cutU[is.na(cutU)] <- TRUE; cutU[no_censor] <- FALSE
  if (verbose) {
    if (!all(LC_eff==0)) {
      if (!attr(LT,"subjectwise")) stat <- mean(cutL) else
        stat <- tapply(cutL,data$subjects,mean)
      message("% lower censoring (after truncation)")
      print(round(100*stat,digits))
    }
    if (!all(UC_eff==Inf)) {
      if (!attr(UT,"subjectwise")) stat <- mean(cutU) else
        stat <- tapply(cutU,data$subjects,mean)
      message("% upper censoring (after truncation)")
      print(round(100*stat,digits))
    }
  }

  # Lower censoring
  data$rt[cutL | cutU] <- NA   # RT wasn't recorded
  data$R[cutL & !LCresponse] <- NA
  data$R[cutU & !UCresponse] <- NA

  # Add missingness column
  data$missingness <- NA_integer_
  data$missingness[cutL & !cutU] <- 1L  # lower-censored only
  data$missingness[cutU & !cutL] <- 2L  # upper-censored only
  data$missingness[cutL &  cutU] <- 3L  # both

  data
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
  post_functions <- NULL
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
    if(!is.null(data)){
      if(nrow(parameters) == length(unique(data$subjects))){
        design$Ffactors$subjects <- unique(data$subjects)
      }
    }
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
                           do_functions = FALSE, ## SM: MIRRORS ZACH'S BRANCH, NOT SURE IF THIS IS WISE
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

  simulate_unconditional_on_data <- return_trialwise_parameters <- FALSE
  dots_local <- list(...)
  if (isFALSE(dots_local$conditional_on_data)) {
    simulate_unconditional_on_data <- TRUE
  } else if (!is.null(dots_local$conditional_on_data)) {
    simulate_unconditional_on_data <- !isTRUE(dots_local$conditional_on_data)
  }
  return_trialwise_parameters <- isTRUE(dots_local$return_trialwise_parameters)
  if('kernel_output_codes' %in% names(dots_local)) {
    kernel_output_codes <- dots_local$kernel_output_codes
  } else {
    kernel_output_codes <- c(1L)
  }

  ## For both conditional and unconditional simulations...
  pars <- t(apply(parameters, 1, do_pre_transform, model()$pre_transform))
  pars <- add_constants(pars,design$constants)
  if(simulate_unconditional_on_data) {
    res <- make_data_unconditional(data=data, pars=pars, design=design, model=model,
                                   return_trialwise_parameters, kernel_output_codes, optionals=optionals)

    data <- res$data
    trialwise_parameters <- res$trialwise_parameters
  } else {
    data <- design_model(
      add_accumulators(data,design$matchfun,simulate=TRUE,type=model()$type,Fcovariates=design$Fcovariates),
      design,model,add_acc=FALSE,compress=FALSE,verbose=FALSE,
      rt_check=FALSE)
    pars <- get_pars_oo(parameters, data, model())
    if(return_trialwise_parameters) {
      if(!is.null(model()$trend)) {
        trialwise_parameters <- cbind(subjects=data$subjects, trials=data$trials,
          pars,
          get_pars_oo(parameters, data, model(),
                      return_kernel_matrix = TRUE,
                      kernel_output_codes = kernel_output_codes)
        )
      } else {
        trialwise_parameters <- pars
      }
    }

    pars <- model()$Ttransform(pars, data)
    if (!is.null(optionals$nobound)) {
      attr(pars,"ok") <- rep(TRUE,nrow(pars))
    } else {
      pars <- fix_bound(pars, model()$bound, data$lR,fix=!is.null(optionals$shrink2bound))
    }
    pars_ok <- attr(pars, 'ok')
    if(mean(!pars_ok) > .1){
      warning("More than 10% of parameter values fall out of model bounds, see <model_name>$bounds()")
      return(FALSE)
    }
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
    dropNames <- c("lR","lM")
    if (!return_Ffunctions && !is.null(design$Ffunctions))
      dropNames <- c(dropNames,names(design$Ffunctions))
    if(!is.null(data$lR)) data <- data[data$lR == levels(data$lR)[1],]
    data <- data[,!(names(data) %in% dropNames)]
    for (i in dimnames(Rrt)[[2]]) data[[i]] <- Rrt[,i]
    if (is_choice_only_model_type(model()) &&
        "rt" %in% names(data) &&
        all(is.na(data$rt))) {
      data$rt <- NULL
    }
  }
  attr(data,"p_vector") <- parameters;
  if(!is.null(post_functions)){
    for(i in 1:length(post_functions)){
      data[[names(post_functions)[i]]] <- post_functions[[i]](data)
    }
  }
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
