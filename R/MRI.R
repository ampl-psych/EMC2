apply_contrasts <- function(events, contrast = NULL, cell_coding = FALSE, remove_intercept = TRUE) {
  factor_name <- events$factor[1]
  colnames(events)[colnames(events) == "event_type"] <- factor_name

  # If a contrast is provided, use it; otherwise let R default to its default contrasts.
  if(!is.null(contrast)){
    if(is.matrix(contrast)){
      if(!is.null(rownames(contrast))){
        events[[factor_name]] <- factor(events[[factor_name]], levels = rownames(contrast))
      } else{
        events[[factor_name]] <- factor(events[[factor_name]])
      }
      stats::contrasts(events[[factor_name]], how.many = ncol(contrast)) <- contrast
    } else {
      events[[factor_name]] <- factor(events[[factor_name]])
      stats::contrasts(events[[factor_name]]) <- do.call(contrast, list(n = length(unique(events[[factor_name]]))))
    }
  } else {
    events[[factor_name]] <- factor(events[[factor_name]])
    # R's default contrasts will be used.
  }

  if(length(unique(events[[factor_name]])) == 1){
    design <- matrix(1, nrow = nrow(events))
    colnames(design) <- factor_name
  } else if(cell_coding){
    design <- model.matrix(as.formula(paste0("~ 0 + ", factor_name)), events)
  } else {
    design <- model.matrix(as.formula(paste0("~ ", factor_name)), events)
    colnames(design)[1] <- paste0(factor_name, "0")
    if(remove_intercept) design <- design[, -1, drop = FALSE]
  }
  events$factor <- NULL
  events_design <- cbind(design, events)

  long_events <- reshape(events_design,
                         direction = "long",
                         varying = colnames(design),
                         v.names = "modulation",
                         timevar = "regressor",
                         times = colnames(design))
  rownames(long_events) <- NULL
  long_events <- long_events[long_events$modulation != 0, ]
  long_events <- long_events[, !(colnames(long_events) %in% c("id", factor_name))]
  long_events <- long_events[order(long_events$onset), ]
  return(long_events)
}


#' Reshape events data for fMRI analysis
#'
#' This function reshapes event data into a format suitable for fMRI analysis by
#' converting specified event_types into separate event types with appropriate modulation values.
#'
#' @param events A data frame containing event information with required columns 'subjects', 'run', and 'onset'
#' @param event_types A character vector of column names in the events data frame to be treated as event_types
#' @param duration Either a single numeric value (applied to all event_types), a list with named elements
#'        corresponding to event_types, or a function that takes the events data frame and returns durations
#' @param modulation Either a list with named elements corresponding to event_types, or a function that takes
#'          the events data frame and returns durations
#' @return A data frame with columns 'subjects', 'onset', 'run', 'modulation', 'duration', and 'event_type'
#' @export
#' @examples
#' # Create a simple events data frame
#' events <- data.frame(
#'   subjects = rep(1, 10),
#'   run = rep(1, 10),
#'   onset = seq(0, 90, by = 10),
#'   condition = rep(c("A", "B"), 5),
#'   rt = runif(10, 0.5, 1.5),
#'   accuracy = sample(0:1, 10, replace = TRUE)
#' )
#'
#' # Reshape with default duration
#' reshaped1 <- reshape_events(events, event_types = c("condition", "accuracy"))
#'
#' # Reshape with custom duration for each event_type
#' reshaped2 <- reshape_events(events,
#'                            event_types = c("condition", "accuracy", "rt"),
#'                            duration = list(condition = 0.5,
#'                                           accuracy = 0.2,
#'                                           rt = function(x) x$rt))
reshape_events <- function(events, event_types, duration = 0.001, modulation = NULL){
  if(!(all(c("onset", "run", "subjects") %in% colnames(events)))){
    stop("Expected columns: subjects, duration, onset, run")
  }


  # First check if only 1 numeric entry is present
  if(length(duration) == 1 && !is.list(duration)){
    duration <- rep(duration, length(event_types))
    duration <- lapply(duration, function(x) return(x)) # and make it into a list
  } else if(is.list(duration) & any(names(duration) %in% event_types)){
    duration_tmp <- replicate(length(event_types), list(0.001)) # Fill in the default spike function
    for(i in 1:length(event_types)){
      if(any(names(duration) %in% event_types[i])){
        duration_tmp[[i]] <- duration[[event_types[i]]]
      }
    }
    duration <- duration_tmp
  } else if(length(duration) != length(event_types)){
    stop("Length of duration must be 1 or equal to the number of event_types")
  }
  if(!is.list(duration)){
    duration <- lapply(duration, function(x) return(x))
  }
  out <- list()
  for(i in 1:length(event_types)){
    fact <- event_types[i]
    tmp <- events[,c('subjects', 'run', 'onset', fact)]
    if(is.character(tmp[,fact]) || is.factor(tmp[,fact])){
      tmp[,fact] <- paste0(fact, "_", tmp[,fact])
      colnames(tmp)[4] <- "event_type"
      if(is.null(modulation[[fact]])){
        tmp$modulation <- 1
      } else{
        if(is.function(modulation[[fact]])){
          tmp$modulation <- modulation[[fact]](events)
        } else{
          tmp$modulation <- modulation[[fact]]
        }
      }
    } else{
      if(!is.null(modulation[[fact]])){
        if(is.function(modulation[[fact]])){
          tmp[,4] <- modulation[[fact]](events)
        } else{
          tmp[,4] <- modulation[[fact]]
        }
      }
      colnames(tmp)[4] <- "modulation"
      tmp$event_type <- fact
    }

    if(is.function(duration[[i]])){
      tmp$duration <- duration[[i]](events)
    } else{
      tmp$duration <- duration[[i]]
    }

    out[[fact]] <- tmp
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out <- out[order(out$subjects, out$run, out$onset),]
  return(out)
}


#' Convolve Events with HRF to Construct Design Matrices
#'
#' This function convolves events with the HRF to construct design matrices for fMRI analysis.
#'
#' @param timeseries A data frame containing fMRI time series data with columns 'subjects', 'run', 'time', and at least one ROI column
#' @param events A data frame containing event information with required columns `subjects`, `run`, `onset`, `duration`, `event_type`, and `modulation`
#' @param factors A named list mapping factor names to event types
#' @param contrasts A named list of contrast matrices for each factor
#' @param covariates A character vector of event types to include as covariates
#' @param add_constant A boolean specifying whether a 1 should be included to the design matrix post convolution
#' @param hrf_model A character string specifying the HRF model to use ('glover', 'spm', 'glover + derivative', or 'spm + derivative')
#' @param cell_coding A character vector of factor names to use cell coding for
#' @param scale A boolean indicating whether to scale the design matrix.
#' @param high_pass Logical indicating whether to apply high-pass filtering.
#' Alternatively, specifying 'add' adds the regressors to the design matrix
#' @param high_pass_model Character indicating which type of high-pass filtering to apply ('cosine', 'poly')
#' @param cut_off A numeric value specifying the cutoff for the high-pass filter
#'
#' @return A list containing the design matrices
#' @export
#' @examples
#' # Generate a simple example timeseries
#' ts <- data.frame(
#'   subjects = rep(1, 100),
#'   run = rep(1, 100),
#'   time = seq(0, 99),
#'   ROI1 = rnorm(100)
#' )
#'
#' # Generate example events
#' events <- data.frame(
#'   subjects = rep(1, 4),
#'   run = rep(1, 4),
#'   onset = c(10, 30, 50, 70),
#'   duration = rep(0.5, 4),
#'   event_type = c("hard", "easy", "hard", "easy"),
#'   modulation = c(1, 1, 1, 1)
#' )
#'
#' # Build design matrices
#' design_matrices <-  convolve_design_matrix(
#'   timeseries = ts,
#'   events = events,
#'   factors = list(difficulty = c("hard", "easy")),
#'   contrasts = list(difficulty = matrix(c(-1, 1)))
#' )
convolve_design_matrix <- function(timeseries, events, factors = NULL, contrasts = NULL,
                                    covariates = NULL, add_constant = TRUE,
                                    hrf_model = 'glover', cell_coding = NULL,
                                    scale = TRUE, high_pass = TRUE,
                                    high_pass_model = "cosine", cut_off = 1e-12) {

  if(!(all(c("onset", "run", "subjects", "modulation", "duration") %in% colnames(events)))){
    stop("Expected columns in events: subjects, duration, onset, run, modulation, duration")
  }

  if(!(all(c("run", "subjects", "time") %in% colnames(timeseries)))){
    stop("Expected columns in frame_times: run, subjects, time")
  }

  if(!setequal(unique(timeseries$subjects), unique(events$subjects))){
    stop("please make sure timeseries and events have the same subjects")
  }

  if(!is.null(cell_coding) && !cell_coding %in% names(factors)) stop("Cell coded factors must have same name as factors argument")
  # Define double-gamma hyperparameters
  if(grepl("glover", hrf_model)){
    undershoot <- 12 # When does the negative time point occur
    dispersion <- .9 # Width of positive gamma
    u_dispersion <- .9 # Width of negative gamma
    ratio <- .48 # Relative size of undershoot compared to overshoot
  } else{ # Different settings for SPM models
    undershoot <- 16
    dispersion <- 1
    u_dispersion <- 1
    ratio <- 1/6
  }

  if(!is.data.frame(events)) events <- do.call(rbind, events)
  subjects <- unique(timeseries$subjects)
  # Holders for filtered design matrix
  all_dms <- list()
  for(subject in subjects){
    ev_sub <- events[events$subjects == subject, ]
    ts_sub <- timeseries[timeseries$subjects == subject, ]
    runs <- unique(ts_sub$run)
    # Define subject-wise new design matrix and timeseries
    dms_sub <- vector("list", length = length(runs))
    for(run in runs){
      ev_run <- ev_sub[ev_sub$run == run, ]
      ts_run <- ts_sub[ts_sub$run == run, ]
      ev_tmp <- data.frame()

      for(fact in names(factors)){
        idx <- ev_run$event_type %in% factors[[fact]]
        ev_run$factor[idx] <- fact
        tmp <- ev_run[idx, ]
        new_tmp <- apply_contrasts(tmp, contrast = contrasts[[fact]],
                                   cell_coding = fact %in% cell_coding)
        rownames(new_tmp) <- NULL
        new_tmp <- cbind(event_type = ev_run$event_type[idx], new_tmp)
        ev_tmp <- rbind(ev_tmp, new_tmp)
      }

      for(cov in covariates){
        idx <- ev_run$event_type %in% cov
        tmp <- ev_run[idx,c('event_type', 'subjects', 'run', 'onset', 'duration', 'modulation')]
        tmp$regressor <- cov
        ev_tmp <- rbind(ev_tmp, tmp)
      }

      ev_tmp <- ev_tmp[order(ev_tmp$onset), ]
      if((run == runs[1]) & (subject == subjects[1])){
        round_ev <- ev_tmp
        round_ev$duration <- round(scale(round_ev$duration))
        round_ev$modulation <- round(scale(round_ev$modulation))
        unq_idx <- !duplicated(round_ev[, !colnames(round_ev) %in% c("onset", "subjects")])
        print(ev_tmp[unq_idx,])
      }
      # Event_type was only included for printing
      ev_tmp <- ev_tmp[,colnames(ev_tmp) != "event_type"]
      dm <- construct_design_matrix(ts_run$time,
                                    events = ev_tmp,
                                    has_derivative = grepl("derivative", hrf_model),
                                    time_length = 32, # total hrf duration
                                    min_onset = -24, # Grid computation start time for oversampling
                                    oversampling = 50, # How many timepoints per tr
                                    onset = 0, # accounts for shifts in HRF
                                    delay = 6, # Time to peak of initial bump
                                    undershoot = undershoot,
                                    dispersion = dispersion,
                                    u_dispersion = u_dispersion,
                                    ratio = ratio,
                                    add_intercept = FALSE)
      if((run == runs[1]) & (subject == subjects[1]) & isTRUE(high_pass)){
        message("Filtering out high_pass noise, make sure you also use high_pass_filter(<timeseries>)")
      }
      if(!isFALSE(high_pass)){
        dm <- high_pass_filter(dm, high_pass_model, frame_times = ts_run$time, add=(high_pass == "add"))
      }
      if(add_constant) dm$constant <- 1
      dms_sub[[as.character(run)]] <- dm
    }
    dms_sub <- do.call(rbind, dms_sub)
    dms_sub[abs(dms_sub) < cut_off] <- 0
    rownames(dms_sub) <- NULL
    all_dms[[as.character(subject)]] <- dms_sub
  }
  if(scale){
    full_dm <- do.call(rbind, all_dms)
    maxs <- apply(full_dm, 2, max)
    all_dms <- lapply(all_dms, function(x){
      for(i in 1:ncol(x)){
        x[,i] <- x[,i]/maxs[i]
      }
      return(x)
    })
  }
  return(all_dms)
}

#' Split fMRI Timeseries Data by ROI Columns
#'
#' This function splits a timeseries data frame containing multiple ROI columns into a list
#' of data frames, where each data frame contains the common columns (subjects, run, time)
#' and one ROI column.
#'
#' @param timeseries A data frame containing fMRI timeseries data with required columns
#'   'subjects', 'run', and 'time', plus one or more ROI columns.
#' @param columns A character vector specifying which columns to split by. If NULL (default),
#'   all columns except 'subjects', 'run', and 'time' will be used.
#'
#' @return A named list of data frames, where each data frame contains the common columns
#'   (subjects, run, time) and one ROI column. The names of the list elements correspond
#'   to the ROI column names.
#'
#' @export
#'
#' @examples
#' # Create a simple example timeseries with multiple ROIs
#' set.seed(123)
#' n_frames <- 100
#'
#' # Create a data frame with multiple ROIs
#' timeseries <- data.frame(
#'   subjects = rep(1, n_frames),
#'   run = rep(1, n_frames),
#'   time = seq(0, n_frames-1),
#'   ROI1 = rnorm(n_frames),
#'   ROI2 = rnorm(n_frames),
#'   ROI3 = rnorm(n_frames)
#' )
#'
#' # Split the timeseries by all ROI columns
#' split_data <- split_timeseries(timeseries)
split_timeseries <- function(timeseries, columns = NULL){
  if(!(all(c("run", "subjects", "time") %in% colnames(timeseries)))){
    stop("Expected columns in timeseries: run, subjects, time")
  }
  if(is.null(columns)) columns <- colnames(timeseries)[!colnames(timeseries) %in% c('run', 'subjects', 'time')]
  out <- list()
  for(col in columns){
    if(!col %in% colnames(timeseries)) stop("Please ensure selected columns are in timeseries")
    out[[col]] <- timeseries[,c('subjects', 'run', 'time', col)]
  }
  return(out)
}


#' Apply High-Pass Filtering to fMRI Data
#'
#' This function applies high-pass filtering to fMRI data to remove low-frequency noise
#' and drift. It supports two filtering methods: cosine basis functions and polynomial
#' regressors.
#'
#' @param X A data frame or matrix containing the data to be filtered. If it contains
#'   columns 'subjects' and 'run', the function will apply filtering separately for
#'   each subject-run combination.
#' @param high_pass_model A character string specifying the high-pass filtering method.
#'   Options are 'cosine' (default) or 'poly' for polynomial regressors.
#' @param frame_times A numeric vector of time points for each frame. If NULL, the
#'   function will attempt to extract this from a 'time' column in X.
#' @param ... Additional arguments passed to the function.
#'
#' @return A data frame or matrix with the same structure as X, but with high-frequency
#'   components removed from the data columns.
#'
#' @export
#'
#' @examples
#' # Create a simple example data frame with drift
#' set.seed(123)
#' n_frames <- 100
#' time <- seq(0, 99)
#'
#' # Create a signal with low-frequency drift
#' drift <- 0.1 * time
#' signal <- sin(2 * pi * 0.1 * time) + drift
#' noise <- rnorm(n_frames, 0, 0.5)
#' data <- signal + noise
#'
#' # Create a data frame
#' df <- data.frame(
#'   time = time,
#'   signal = data
#' )
#'
#' # Apply high-pass filtering using cosine basis functions
#' filtered_df <- high_pass_filter(df, high_pass_model = "cosine")
high_pass_filter <- function(X, high_pass_model = 'cosine', frame_times = NULL, ...){
  if(is.null(frame_times)){
    if(!'time' %in% colnames(X)){
      stop("no column named 'time' for frame_times present, please separately provide")
    } else{
      frame_times <- X[,'time']
    }
  }
  if('subjects' %in% colnames(X) && is.null(list(...)$recursive)){
    out <- list()
    k <- 0
    message("Make sure you also high_pass_filter your events (set high_pass = TRUE in convolve_design_matrix)")
    for(sub in unique(X[,'subjects'])){
      tmp <- X[X[,'subjects'] == sub,]
      for(run in unique(tmp[,'run'])){
        k <- k + 1
        tmp_run <- tmp[tmp[,'run'] == run,]
        tmp_run <- high_pass_filter(tmp_run, recursive = TRUE)
        out[[k]] <- tmp_run
      }
    }
    return(do.call(rbind, out))
  }
  if(high_pass_model == "cosine"){
    nuisance <- cosine_drift(frame_times)
  } else if(high_pass_model == "poly"){
    nuisance <- poly_drift(frame_times)
  } else{
    stop("Only poly and cosine are supported as high_pass_model")
  }
  if(!is.null(list(...)$add)){
    if(list(...)$add){
      return(cbind(X, nuisance))
    }
  }

  gets_filter <- !colnames(X) %in% c('subjects', 'run', 'time')
  for(i in 1:ncol(X)){
    if(gets_filter[i]){
      fit <- lm(X[,i] ~ nuisance - 1)
      X[,i] <- residuals(fit)
    }
  }
  return(X)
}

poly_drift <- function(frame_times, order = 3) {
  # Ensure that 'order' is an integer
  order <- as.integer(order)
  n <- length(frame_times)
  # Compute maximum of frame_times (will be used to scale the time vector)
  tmax <- max(frame_times)

  # Create a matrix where column k corresponds to (frame_times/tmax)^k,
  # for k = 0, 1, ..., order; note that the 0th power yields a constant.
  pol <- sapply(0:order, function(k) (frame_times / tmax) ^ k)
  # 'pol' is now an n x (order+1) matrix

  # Orthogonalize the columns using QR decomposition.
  # The function qr.Q returns an orthonormal basis for the columns.
  pol_orth <- qr.Q(qr(pol))

  # Drop the constant
  result <- pol_orth[, -1, drop = FALSE]
  colnames(result) <- paste0("pol_", 1:ncol(result))
  return(result)
}


cosine_drift <- function(frame_times, high_pass = .01) {
  n_frames <- length(frame_times)
  n_times <- 0:(n_frames - 1)
  # Compute the time interval (dt) between frames.
  dt <- (frame_times[n_frames] - frame_times[1]) / (n_frames - 1)

  # Check if the product high_pass * dt is too high, issuing a warning if so.
  if (high_pass * dt >= 0.5) {
    warning(sprintf("High-pass filter will span all accessible frequencies and saturate the design matrix. You may want to reduce the high_pass value. The provided value is %.3f Hz", high_pass))
  }

  # Determine the number of cosine basis functions (excluding the constant).
  # This follows: order = min(n_frames - 1, floor(2 * n_frames * high_pass * dt))
  order <- min(n_frames - 1, floor(2 * n_frames * high_pass * dt))

  # Create a matrix to hold the cosine basis functions and a constant column.
  # The result will have (order + 1) columns.
  cosine_drift <- matrix(0, nrow = n_frames, ncol = order + 1)
  normalizer <- sqrt(2.0 / n_frames)

  # Fill the first 'order' columns with the cosine functions.
  # For each k = 1, 2, ..., order we compute:
  #   normalizer * cos( (pi/n_frames) * (n_times + 0.5) * k )
  for (k in 1:order) {
    cosine_drift[, k] <- normalizer * cos((pi / n_frames) * (n_times + 0.5) * k)
  }
  cosine_drift <- cosine_drift[,-ncol(cosine_drift)]

  # Set the last column to a constant of 1.
  colnames(cosine_drift) <- paste0("drift_", 1:ncol(cosine_drift))

  return(cosine_drift)
}


#' Create fMRI Design for EMC2 Sampling
#'
#' This function takes the output from convolve_design_matrix and transforms it into a design
#' suitable for sampling with EMC2. It properly configures parameter types, bounds, and transformations
#' for the specified model.
#'
#' @param design_matrix A list of design matrices,  the output from convolve_design_matrix
#' @param model A function that returns a model specification, options are MRI or MRI_AR1
#' @param ... Additional arguments passed to the model
#'
#' @return An object of class 'emc.design' suitable for EMC2 sampling
#' @export
#'
#' @examples
#' # Generate a simple example timeseries
#' ts <- data.frame(
#'  subjects = rep(1, 100),
#'  run = rep(1, 100),
#'  time = cumsum(rep(1.38, 100)),
#'  ROI1 = rnorm(100)
#' )
#'
#' # Generate example events
#' events <- data.frame(
#'  subjects = rep(1, 4),
#'  run = rep(1, 4),
#'  onset = c(10, 30, 50, 70),
#'  duration = rep(0.5, 4),
#'  event_type = c("A", "B", "A", "B"),
#'  modulation = c(1, 1, 1, 1)
#' )
#'
#' # Create convolved design matrix
#' design_matrix <- convolve_design_matrix(
#'  timeseries = ts,
#'  events = events,
#'  factors = list(condition = c("A", "B")),
#'  hrf_model = "glover"
#' )
#'
#' # Create fMRI design for EMC2
#' fmri_design <- design_fmri(design_matrix, model = MRI_AR1)

design_fmri <- function(design_matrix,
                             model = MRI_AR1, ...) {
  dots <- list(...)
  betas <- colnames(design_matrix[[1]])
  subjects <- names(design_matrix)
  model_list <- model()
  # Fill in new p_types
  p_not_beta <- model_list$p_types[names(model_list$p_types) != "beta"]
  model_list$p_types <- c(setNames(rep(model_list$p_types["beta"], length(betas)), betas), p_not_beta)
  # Fill in new bound
  new_mm <- do.call(cbind, rep(list(model_list$bound$minmax[,"beta"]), length(betas)))
  colnames(new_mm) <- betas
  model_list$bound$minmax <- cbind(model_list$bound$minmax, new_mm)
  model_list$bound$minmax <- model_list$bound$minmax[,colnames(model_list$bound$minmax) != "beta"]

  # Fill in new transforms
  new_t <- setNames(rep(model_list$transform$func["beta"], length(betas)), betas)
  model_list$transform$func <- c(model_list$transform$func, new_t)
  model_list$transform$func <- model_list$transform$func[names(model_list$transform$func) != "beta"]
  # Make pre transforms
  par_names <- c(betas, names(p_not_beta))
  model_list$pre_transform$func <- setNames(rep("identity", length(par_names)), par_names)
  model <- function() {return(model_list)}
  n_pars <- length(par_names)

  # Fill up final results
  model_list$transform <- fill_transform(dots$transform,model)
  model_list$bound <- fill_bound(dots$bound,model)
  model_list$pre_transform <- fill_transform(dots$pre_transform, model = model, p_vector = model_list$p_types, is_pre = TRUE)
  model <- function(){return(model_list)}



  Flist <- vector("list", n_pars)
  for(i in 1:n_pars){
    Flist[[i]] <- as.formula(paste0(par_names[i], "~1"))
  }

  design <- list(Flist = Flist, model = model, Ffactors = list(subjects = subjects))
  attr(design, "design_matrix") <- lapply(design_matrix, FUN=function(x) {
    y <- x[,colnames(x) != 'subjects']
    DM_tmp <- data.matrix(y)
    rownames(DM_tmp) <- NULL
    return(DM_tmp)
  })
  par_names <- setNames(numeric(length(par_names)), par_names)
  attr(design, "p_vector") <- par_names
  class(design) <- 'emc.design'
  return(design)
}


#' GLM model for fMRI data
#'
#' Creates a model specification for fMRI data using a normal distribution.
#' This model assumes that the observed BOLD signal follows a normal distribution
#' with a mean determined by the design matrix and betas, and a standard deviation
#' parameter for noise.
#'
#' @return A list containing model specification
#'
#' @details
#' The model uses a normal distribution to model fMRI BOLD signals.
#' Beta parameters represent the effect sizes for different conditions,
#' and the sd parameter represents the standard deviation of the noise.
#'
#' The log-likelihood function centers the predicted values by subtracting
#' the mean, which helps with model identifiability.
#'
#' @export
#' @examples
#' # Create a normal MRI model specification
#' model_spec <- MRI()
#'
#' # Access model parameters
#' model_spec$p_types
MRI <- function(){
  return(
    list(
      type="MRI",
      c_name = "MRI",
      p_types=c("beta" = 0, "sd" = log(1)),
      transform=list(func=c(beta = "identity", sd = "exp")),
      bound=list(minmax=cbind(beta=c(-Inf,Inf),sd=c(0.001,Inf))),
      Ttransform = function(pars, dadm) return(pars),
      rfun=function(pars){
        #   - Each row corresponds to an observation
        #   - All columns except the last are betas (already multiplied by the design matrix)
        #   - The last column is sigma (the noise standard deviation)
        # Extract sigma and betas
        sigma <- pars[, ncol(pars)]
        betas <- pars[, -ncol(pars)]
        # Compute the predicted mean for each observation as the sum of its betas
        y_hat <- rowSums(betas)
        # Generate simulated data: for each observation, add noise drawn from a normal distribution
        # with mean 0 and standard deviation sigma.
        y_sim <- y_hat + rnorm(n = length(y_hat), mean = 0, sd = sigma)
        return(y_sim)
      },
      log_likelihood=function(pars, dadm, model, min_ll=log(1e-10)){
        # Here pars is already multiplied by design matrix in map_p
        y <- as.matrix(dadm[,!colnames(dadm) %in% c("subjects", 'run', 'time', "trials")])
        # grab the right parameters
        sigma <- pars[,ncol(pars)]
        betas <- pars[,-ncol(pars)]
        y_hat <- rowSums(betas)
        ll <- sum(pmax(dnorm(as.matrix(y), mean = y_hat, sd = sigma, log = T), min_ll))
        return(ll)
      }
    )
  )
}

#' Create an AR(1) GLM model for fMRI data
#'
#' This function creates a model specification for MRI data with an AR(1) error structure.
#' The model includes beta parameters for the design matrix, a rho parameter for the
#' autocorrelation, and a standard deviation parameter for the noise.
#'
#' The AR(1) model accounts for temporal autocorrelation in the data, where each timepoint
#' is correlated with the previous timepoint according to the rho parameter.
#'
#' @return A list containing the model specifications
#'
#' @export
#' @examples
#' # Create an AR(1) GLM model for fMRI data
#' model_spec <- MRI_AR1()
#'
#' # Access model parameters
#' model_spec$p_types

MRI_AR1 <- function(){
  return(
    list(
      type="MRI_AR1",
      c_name = "MRI_AR1",
      p_types=c("beta" = 0, "rho" = pnorm(0.001), "sd" = log(1)),
      transform=list(func=c(beta = "identity",  rho = "pnorm", sd = "exp")),
      bound=list(minmax=cbind(beta=c(-Inf,Inf),sd=c(0.001,Inf), rho = c(0.0001, 1)),
                 exception=c(rho=0)),
      Ttransform = function(pars, dadm) return(pars),
      rfun=function(pars){
        n <- nrow(pars)
        m <- ncol(pars)

        # - betas: columns 1 to (m-2)
        # - rho: column (m-1)
        # - sigma: column m (stationary standard deviation)
        betas <- pars[,1:(m-2), drop = FALSE]
        rho   <- pars[,m-1]
        sigma <- pars[,m]

        # Compute the linear predictor (sum of beta contributions) and center it.
        y_hat <- rowSums(betas)
        # y_hat <- y_hat - mean(y_hat)

        # Allocate a vector for simulated data
        y_sim <- numeric(n)

        # Simulate the first observation
        y_sim[1] <- y_hat[1] + rnorm(1, mean = 0, sd = sigma[1])

        # Loop through time for the remaining observations
        for (t in 2:n) {
          # The conditional mean for observation t
          cond_mean <- y_hat[t] + rho[t] * (y_sim[t - 1] - y_hat[t - 1])
          # The conditional standard deviation for observation t
          cond_sd <- sigma[t] * sqrt(1 - rho[t]^2)
          # Simulate
          y_sim[t] <- cond_mean + rnorm(1, mean = 0, sd = cond_sd)
        }
        return(y_sim)
      },
      log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
        # Here pars is already multiplied by design matrix in map_p

        # Extract observed data (as a vector)
        y <- as.vector(as.matrix(dadm[, !colnames(dadm) %in% c("subjects", "run", "time", "trials")]))
        n <- length(y)

        m <- ncol(pars)
        betas <- pars[, 1:(m - 2), drop = FALSE]
        rho   <- pars[, m - 1]
        sigma <- pars[, m]

        y_hat <- rowSums(betas)
        # y_hat <- y_hat - mean(y_hat)

        # Log-likelihood for the first observation
        ll <- numeric(n)
        ll[1] <- dnorm(y[1], mean = y_hat[1], sd = sigma[1], log = TRUE)

        # For observations t = 2:n, compute conditional means
        cond_mean <- y_hat[-1] + rho[-1] * (y[-n] - y_hat[-n])
        cond_sd   <- sigma[-1] * sqrt(1 - rho[-1]^2)
        ll[-1] <- dnorm(y[-1], mean = cond_mean, sd = cond_sd, log = TRUE)

        ll <- pmax(ll, min_ll)
        return(sum(ll))
      }
    )
  )
}


# Plotting functions ------------------------------------------------------

#' Plot fMRI Design Matrix
#'
#' This function creates a visualization of an fMRI design matrix, showing the temporal
#' evolution of regressors over time. It can handle various input formats and provides
#' options to customize the visualization.
#'
#' @param design_matrix A design matrix for fMRI analysis. Can be a data frame, matrix,
#'   list of matrices, or an object of class 'emc.design'.
#' @param TRs The number of time points (TRs) to plot. Default is 100.
#' @param events A character vector specifying which regressors to plot. If NULL,
#'   all non-nuisance regressors will be plotted.
#' @param remove_nuisance Logical indicating whether to remove nuisance regressors
#'   (drift terms, polynomial terms, derivatives) from the plot. Default is TRUE.
#' @param subject The subject number to plot. Only applies for list of design matrices. Default is 1.
#' @param legend_pos Position of the legend. Default is "bottomleft".
#' @param ... Additional graphical parameters passed to matplot and legend.
#'
#' @return A plot showing the design matrix regressors over time.
#'
#' @export
#'
#' @examples
#' # Example time series
#' ts <- data.frame(
#'   subjects = rep(1, 100),
#'   run      = rep(1, 100),
#'   time     = seq(0, 99),
#'   ROI      = rnorm(100)
#' )
#' # Create a simple events data frame
#' events <- data.frame(
#'   subjects = rep(1, 10),
#'   run = rep(1, 10),
#'   onset = seq(0, 90, by = 10),
#'   condition = rep(c("A", "B"), 5),
#'   rt = runif(10, 0.5, 1.5),
#'   accuracy = sample(0:1, 10, replace = TRUE)
#' )
#' # Reshape with custom duration for each event_type
#' reshaped <- reshape_events(events,
#'                            event_types = c("condition", "accuracy", "rt"),
#'                            duration = list(condition = 0.5,
#'                                           accuracy = 0.2,
#'                                           rt = function(x) x$rt))
#' design_matrices <-  convolve_design_matrix(
#'                     timeseries = ts,
#'                     events = reshaped,
#'                     covariates = c('accuracy', 'rt'),
#'                     factors = list(cond = c("condition_A", "condition_B")),
#'                     contrasts = list(cond = matrix(c(-1, 1))))
#'
#' # Plot the design matrix
#' plot_design_fmri(design_matrices)
plot_design_fmri <- function(design_matrix, TRs = 100, events = NULL, remove_nuisance = TRUE, subject = 1,
                             legend_pos = "bottomleft", ...){
  if(is(design_matrix,"emc.design")){
    design_matrix <- attr(design_matrix, "design_matrix")
  }
  if(!is.data.frame(design_matrix) && !is.matrix(design_matrix)){
    if(is.list(design_matrix)) design_matrix <- design_matrix[[subject]]
  }
  enames <- colnames(design_matrix)
  if(remove_nuisance & is.null(events)){
    is_nuisance <- grepl("drift", enames) | grepl("poly", enames) | grepl("derivative", enames) | apply(design_matrix, 2, sd) == 0
    design_matrix <- design_matrix[,!is_nuisance, drop = F]
  }
  enames <- colnames(design_matrix)
  if(is.null(events)){
    events <- enames
  } else{
    if(any(!events %in% enames)){
      stop("events not in colnames design matrix")
    }
  }
  distinct_colors <- c(
    "#E6194B", "#3CB44B", "#0082C8", "#F58231", "#911EB4",
    "#46F0F0", "#F032E6", "#D2F53C", "#FABEBE", "#008080",
    "#E6BEFF", "#AA6E28", "#FFFAC8", "#800000", "#FFD8B1"
  )
  dots <- add_defaults(list(...), col = distinct_colors, lwd = 2, lty = 1, main = NULL, xlab = "TRs", ylab = "Amplitude")

  TRs <- min(c(TRs, nrow(design_matrix)))
  design_matrix <- design_matrix[1:TRs, events, drop = F]
  do.call(matplot, c(list(design_matrix, type = "l"), fix_dots_plot(dots)))
  do.call(legend, c(list(legend_pos, legend = events, bty = "n"), fix_dots(dots, legend)))
}



# -------------------------------------------------------------------------
# Filter events by event_type and bin/categorize the modulation
# -------------------------------------------------------------------------
prepare_event_groups <- function(events, event_type, n_bins = 4) {
  ev_sub <- events[events$event_type == event_type, ]
  if (nrow(ev_sub) == 0) {
    return(ev_sub)  # empty
  }
  if (!"modulation" %in% names(ev_sub)) {
    stop("The 'events' data frame must have a 'modulation' column.")
  }

  # Round to reduce floating precision issues
  ev_sub$modulation <- round(ev_sub$modulation, 6)
  mod_vals <- unique(ev_sub$modulation)
  n_unique <- length(mod_vals)

  if (n_unique > 6) {
    # Calculate quartile breakpoints for the 'modulation' data
    quartile_breaks <- quantile(ev_sub$modulation, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)

    # Bin the data into quartiles using these breakpoints
    ev_sub$mod_group <- cut(ev_sub$modulation, breaks = quartile_breaks, include.lowest = TRUE, labels = 1:n_bins)
    attr(ev_sub, "binned") <- TRUE

  } else {
    # treat as categorical
    ev_sub$mod_group <- factor(ev_sub$modulation, levels = sort(mod_vals))
    attr(ev_sub, "binned") <- FALSE
  }
  return(ev_sub)
}

# -----------------------------------------------------------------------------
# Fast FIR utilities -----------------------------------------------------------
# -----------------------------------------------------------------------------

# 1) Build an FIR design matrix for one set of onsets --------------------------
#    frame_times : numeric vector of acquisition times (s) **monotone**
#    onsets      : numeric vector of event onsets (s)  (single modulation group)
#    pre, post   : seconds before/after onset to model (positive numbers)
#
# Returns a dense matrix (length(frame_times) x K) where
#   K = floor((pre+post)/dt)+1  (dt = median diff(frame_times))
# Column k corresponds to lag t = -pre + (k-1)*dt.
# ---------------------------------------------------------------------------
# Build an FIR design with an optional externally‑specified TR  --------------
# ---------------------------------------------------------------------------
###############################################################################
## Fast FIR utilities + plot_fmri (robust version, 3 May 2025)               ##
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# 1)  FIR design – *global lag grid* ensured by `fixed_dt` argument
# ─────────────────────────────────────────────────────────────────────────────
build_fir_design <- function(frame_times, onsets, durations,
                             pre, post, fixed_dt, weights = NULL) {

  if (is.null(weights))   weights   <- rep(1, length(onsets))
  if (is.null(durations)) durations <- rep(0, length(onsets))  # impulse
  stopifnot(length(durations) == length(onsets))

  dt   <- fixed_dt
  lags <- seq(-pre, post, by = dt)
  K    <- length(lags)
  nTR  <- length(frame_times)

  X <- matrix(0, nrow = nTR, ncol = K,
              dimnames = list(NULL, sprintf("lag_%0.3f", lags)))

  ## pre‑compute frame centres for fast look‑up
  fc <- frame_times

  for (j in seq_along(onsets)) {
    onset <- onsets[j]
    dur   <- durations[j]
    w     <- weights[j]

    ## frames whose centre is within the event window
    in_evt <- which(fc >= onset & fc < onset + dur)
    if (!length(in_evt)) {              # shorter than first half‑TR → impulse
      in_evt <- which.min(abs(fc - onset))
    }

    ## update FIR cols for every lag
    for (k in seq_len(K)) {
      rows <- in_evt + round(lags[k] / dt)
      rows <- rows[rows >= 1 & rows <= nTR]
      if(length(rows))
        X[rows, k] <- X[rows, k] + w
    }
  }
  attr(X, "lag_seconds") <- lags
  X
}


# ─────────────────────────────────────────────────────────────────────────────
# 2)  AR(1) whitening
# ─────────────────────────────────────────────────────────────────────────────
estimate_rho <- function(y, clip = TRUE) {
  if (length(y) < 2L) return(0)

  y <- y - mean(y, na.rm = TRUE)        # 1 de‑mean
  rho <- sum(y[-1] * y[-length(y)]) /
    (sum(y[-length(y)]^2) + 1e-12)

  if (clip) rho <- max(-0.99, min(0.99, rho))  # 2 safe bounds
  rho
}

whiten_series <- function(mat_or_vec, rho) {
  if (abs(rho) < 1e-6) return(mat_or_vec)

  if (is.vector(mat_or_vec))
    return(c(mat_or_vec[1],
             mat_or_vec[-1] - rho * mat_or_vec[-length(mat_or_vec)]))

  rbind(mat_or_vec[1, , drop = FALSE],
        mat_or_vec[-1, , drop = FALSE] -
          rho * mat_or_vec[-nrow(mat_or_vec), , drop = FALSE])
}


# ─────────────────────────────────────────────────────────────────────────────
# 4)  Quick FIR GLM fit (OLS + AR(1) whitening, baseline correction)
# ─────────────────────────────────────────────────────────────────────────────
fit_fir_glm <- function(y, X, lags_sec, target) {
  rho <- estimate_rho(y)
  yw  <- whiten_series(y, rho)
  Xw  <- whiten_series(X, rho)

  beta <- as.vector(solve(crossprod(Xw), crossprod(Xw, yw)))
  beta <- beta[target]
  ## Baseline correction (pre‑stimulus lags)
  idx_pre <- lags_sec < 0
  if(any(idx_pre, na.rm = TRUE))
    beta <- beta - mean(beta[idx_pre], na.rm = TRUE)
  beta
}

# ─────────────────────────────────────────────────────────────────────────────
# 5)  Main extractor – returns long data.frame for plotting
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
# 5)  Main extractor – FIR for *all* events, returns long data‑frame
# ─────────────────────────────────────────────────────────────────────────────
get_fir_lines <- function(timeseries, events, event_type,
                          pre = 2, post = 18, n_bins = 4,
                          high_pass = TRUE, high_pass_model = "cosine") {

  ## basic checks
  stopifnot(all(c("run", "subjects", "time") %in% names(timeseries)))
  ROI_col <- setdiff(names(timeseries),
                     c("run", "subjects", "time", "postn"))
  if(length(ROI_col) != 1L) stop("Exactly one ROI column required")

  ## global lag grid
  global_dt   <- median(diff(sort(unique(timeseries$time))))
  global_lags <- seq(-pre, post, by = global_dt)
  K           <- length(global_lags)

  ## TARGET events (binned)
  ev_sub  <- prepare_event_groups(events, event_type, n_bins)
  if(nrow(ev_sub) == 0)
    stop("No events found for event_type = ", event_type)
  groups  <- sort(unique(ev_sub$mod_group))

  ## NUISANCE events  (everything else)
  ev_nuis <- events[events$event_type != event_type, ]
  nuis_types <- sort(unique(ev_nuis$event_type))    # might be empty

  betas <- setNames(vector("list", length(groups)), groups)

  ts_split <- split(timeseries, list(timeseries$subjects,
                                     timeseries$run), drop = TRUE)

  for(chunk in ts_split) {
    ft  <- chunk$time
    y   <- chunk[[ROI_col]]
    sid <- chunk$subjects[1]
    rid <- chunk$run[1]

    ## 1) build FIR blocks for *each nuisance type* in this run
    X_nuis <- NULL
    if(length(nuis_types)) {
      for(et in nuis_types) {
        ev_tmp <- ev_nuis[ev_nuis$event_type == et &
                            ev_nuis$subjects   == sid &
                            ev_nuis$run        == rid, ]
        if(!nrow(ev_tmp)) next

        if(sd(ev_tmp$modulation) == 0 & sd(ev_tmp$duration) == 0) next


        ## centre the modulator so ‘main effect’ and modulation are orthogonal
        if (sd(ev_tmp$modulation) > 0) {
          ev_tmp$modulation <- ev_tmp$modulation - mean(ev_tmp$modulation)
        } else{
          ev_tmp$modulation <- ev_tmp$duration - mean(ev_tmp$duration)
        }
        Xet <- build_fir_design(ft,
                                ev_tmp$onset,
                                ev_tmp$duration,
                                pre, post,
                                fixed_dt = global_dt,
                                weights  = ev_tmp$modulation)
        colnames(Xet) <- paste0(et, "_", colnames(Xet))
        X_nuis <- if(is.null(X_nuis)) Xet else cbind(X_nuis, Xet)
      }
    }

    ## 2) loop over modulation groups of the TARGET event
    for(g in groups) {
      ev_g <- ev_sub[ev_sub$mod_group == g &
                       ev_sub$subjects  == sid &
                       ev_sub$run       == rid, ]
      if(!nrow(ev_g)) next

      ## NOTE: weights = NULL → un‑scaled FIR for the curve we plot
      X_tar <- build_fir_design(ft,
                                ev_g$onset,
                                ev_g$duration,       # durations for the target
                                pre, post,
                                fixed_dt = global_dt,
                                weights  = NULL)     # no amplitude scaling

      X_full <- if(is.null(X_nuis)) X_tar else cbind(X_tar, X_nuis)

      if(!isFALSE(high_pass)){
        X_full <- high_pass_filter(X_full, high_pass_model, frame_times = ft, add=(high_pass == "add"))
      }
      b  <- fit_fir_glm(y, X_full, global_lags, target = seq_len(K))

      if(all(is.na(b))) next
      betas[[g]] <- rbind(betas[[g]],
                          cbind(t(b), nrow(ev_g)))
    }
  }

  ## 3) weight‑averaged β → long data.frame
  out <- lapply(names(betas), function(g) {
    mat <- betas[[g]]; if(is.null(mat)) return(NULL)
    B <- mat[, 1:K, drop = FALSE]; w <- mat[, K+1L]
    if(!is.matrix(B)) { B <- matrix(B, nrow = 1); w <- w[1] }

    beta_avg <- colSums(B * w) / sum(w)

    df <- data.frame(
      mod_group  = rep(g, K),
      time       = global_lags,
      avg_signal = beta_avg
    )
    if(attr(ev_sub, "binned")) {
      ev_all <- ev_sub[ev_sub$mod_group == g, ]
      df$min_mod <- rep(min(ev_all$modulation), K)
      df$max_mod <- rep(max(ev_all$modulation), K)
    }
    df
  })
  rownames_out <- NULL
  do.call(rbind, out)
}


#' Plot fMRI peri-stimulus time courses
#'
#' This function plots average BOLD response around specified events for a single ROI
#' by using FIR based event estimation, all event_types in events are taken into account in the FIR.
#' Posterior predictives can be overlaid via the `post_predict` argument.
#'
#' @param timeseries A data frame with columns 'subjects', 'run', 'time', and one ROI measurement column.
#' @param events A data frame with columns 'subjects', 'run', 'onset', 'duration', 'event_type', and 'modulation'.
#' @param event_type Character string specifying which `event_type` in `events` to plot.
#' @param high_pass Logical indicating whether to apply high-pass filtering.
#' Alternatively, specifying 'add' adds the regressors to the design matrix in the FIR.
#' The choice here should be the same as the choice for `convolve_design_matrix`
#' @param high_pass_model Character indicating which type of high-pass filtering to apply ('cosine', 'poly')
#' @param post_predict Optional posterior predictive samples data frame (not shown in examples).
#' @param posterior_args Named list of graphical parameters for posterior predictive lines.
#' @param legend_pos Position of the legend. Default: "topleft".
#' @param layout Panel layout matrix for multiple modulation groups. NULL leaves current layout
#' @param n_cores Number of cores to calculate FIR across subjects with.
#' @param ... Additional graphical parameters passed to plotting functions (e.g., col, lwd, lty).
#'
#' @return NULL. Produces plots as a side-effect.
#' @export
#'
#' @examples
#' ts <- data.frame(
#'   subjects = rep(1, 100),
#'   run      = rep(1, 100),
#'   time     = seq(0, 99),
#'   ROI      = rnorm(100)
#' )
#' events <- data.frame(
#'   subjects   = rep(1, 5),
#'   run        = rep(1, 5),
#'   onset      = c(10, 30, 50, 70, 90),
#'   event_type = rep("A", 5),
#'   modulation = rep(1, 5),
#'   duration   = rep(0.5, 5)
#' )
#' plot_fmri(ts, events = events, event_type = "A")
plot_fmri <- function(timeseries, post_predict = NULL, events, event_type,
                       high_pass = TRUE, high_pass_model = "cosine",
                       posterior_args = list(),
                       legend_pos = "topleft", layout = NA,
                       n_cores = 1, ...) {

  posterior_args <- add_defaults(posterior_args, col = "darkgreen", lwd = 2)
  plot_args <- add_defaults(list(...), col = "black", lwd = 2,
                            xlab = "time (s)", ylab = "BOLD response")

  ROI_col <- setdiff(names(timeseries),
                     c("run", "subjects", "time", "postn"))

  ts <- get_fir_lines(timeseries, events, event_type,
                      high_pass       = high_pass,
                      high_pass_model = high_pass_model)

  ## posterior predictive overlay
  if(!is.null(post_predict)) {
    pp_split <- split(post_predict, post_predict$postn)
    pp_lines <- auto_mclapply(pp_split, get_fir_lines, events, event_type,
                              high_pass       = high_pass,
                              high_pass_model = high_pass_model, mc.cores = n_cores)
    pp_df <- do.call(rbind, pp_lines)

    qfun  <- function(x) quantile(x, c(.025, .5, .975))
    agg   <- aggregate(avg_signal ~ mod_group + time,
                       data = pp_df, FUN = qfun)
    res_df <- data.frame(agg[, 1:2], agg$avg_signal)
    names(res_df)[3:5] <- c("p025", "p50", "p975")
  }

  ## layout
  n_plots <- length(unique(ts$mod_group))
  if(n_plots == 1){
    insert <- ""
  } else{
    insert <- "Q"
  }

  if(!is.null(layout)){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }
  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = n_plots, nplots = 1))
  } else {
    par(mfrow = layout)
  }

  ## loop over modulation groups
  for(mod in unique(ts$mod_group)){
    if(!is.null(post_predict)){
      ylim <- c(min(ts$avg_signal, res_df$p025), max(ts$avg_signal, res_df$p975))
    } else{
      ylim <- range(ts$avg_signal)
    }
    idx_ts <- ts$mod_group == mod
    tmp_plot_args <- add_defaults(plot_args, ylim = ylim,
                                  main = paste0(ROI_col, " - ", event_type, ifelse(n_plots == 1, "", paste0(": ", insert, mod))))
    do.call(plot, c(list(ts$time[idx_ts], ts$avg_signal[idx_ts], type = "l"), fix_dots_plot(tmp_plot_args)))
    abline(h = 0, lty = 2)
    if(!is.null(post_predict)){
      idx_pp <- res_df$mod_group == mod
      do.call(lines, c(list(res_df$time[idx_pp], res_df$p50[idx_pp]), fix_dots_plot(posterior_args)))
      polygon_args <- posterior_args
      polygon_args$col <- adjustcolor(polygon_args$col, alpha.f = .2)
      do.call(polygon, c(list(
        x = c(res_df$time[idx_pp], rev(res_df$time[idx_pp])),
        y = c(res_df$p025[idx_pp], rev(res_df$p975[idx_pp])), border = NA), fix_dots_plot(polygon_args)))
    }
    dots_legend <- plot_args
    dots_legend$col <- c(plot_args$col, posterior_args$col)
    dots_legend$lty <- c(plot_args$lty, posterior_args$lty)
    dots_legend$lwd <- c(plot_args$lwd, posterior_args$lwd)
    if(is.null(post_predict)){
      legend_in <- "data"
    } else{
      legend_in <- c("data", "posterior")
    }
    do.call(legend, c(list(legend_pos, legend = legend_in, bty = "n"), fix_dots(dots_legend, legend)))
  }
}



# Utility functions -------------------------------------------------------



make_mri_sampling_design <- function(design, sampled_p_names){
  out <- list()
  expand <- 1:nrow(design)
  rownames(design) <- NULL
  for(i in 1:ncol(design)){
    out[[i]] <- design[,i, drop = F]
    attr(out[[i]], "expand") <- expand
  }
  names(out) <- colnames(design)
  no_design <- sampled_p_names[!sampled_p_names %in% colnames(design)]
  for (nam in no_design){
    mat <- matrix(1, 1, 1)
    colnames(mat) <- nam
    out[[nam]] <- mat
    attr(out[[nam]], "expand") <- rep(1, nrow(design))
  }
  return(out)
}


make_data_wrapper_MRI <- function(parameters, data, design){
  data_list <- list()
  for(i in 1:nrow(parameters)){
    data_tmp <- data[data$subjects == unique(data$subjects)[i],]
    data_tmp$subjects <- factor(data_tmp$subjects)
    attr(data_tmp, "designs") <- design$fMRI_design[[i]]

    data_list[[i]] <- make_data_fMRI(parameters[i,,drop = F], design$model, data_tmp, design)
  }
  return(do.call(rbind, data_list))
}

make_data_fMRI <- function(parameters, model, data, design, ...){
  # if(is.null(attr(design, "design_matrix"))){
  #   stop("for fMRI simulation the original design needs to be passed to the simulation function")
  # }
  pars <- t(apply(parameters, 1, do_pre_transform, model()$pre_transform))
  pars <- map_p(add_constants(pars,design$constants),data, model())
  pars <- do_transform(pars, model()$transform)
  pars <- model()$Ttransform(pars, data)
  pars <- add_bound(pars, model()$bound)
  data[, !colnames(data) %in% c("subjects", "run", "time", "trials")] <- model()$rfun(pars)
  return(data)
}

add_design_fMRI_predict <- function(design, emc){
  design$fMRI_design <- lapply(emc[[1]]$data, function(x) return(attr(x,"designs")))
  return(design)
}

