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
#'
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
      if(names(duration) %in% event_types[i]){
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
                                    hrf_model = 'glover + derivative', cell_coding = NULL,
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
    ratio <- .35 # Relative size of undershoot compared to overshoot
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
        warning("Filtering out high_pass noise, make sure you also use high_pass_filter(<timeseries>)")
      }
      dm <- high_pass_filter(dm, high_pass_model, frame_times = ts_run$time, add=(high_pass == "add"))
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
#' @param model A function that returns a model specification, options are normal_mri or white_mri
#' @param ... Additional arguments passed to the model
#'
#' @return An object of class 'emc.design' suitable for EMC2 sampling
#' @export
#'
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
#'   event_type = c("condition", "condition", "condition", "condition"),
#'   modulation = c(1, 1, 1, 1)
#' )
#'
#' # Create convolved design matrix
#' design_matrix <- convolve_design_matrix(
#'   timeseries = ts,
#'   events = events,
#'   factors = list(condition = "condition"),
#'   hrf_model = "glover"
#' )
#'
#' # Create fMRI design for EMC2
#' fmri_design <- design_fmri(design_matrix, model = white_mri)

design_fmri <- function(design_matrix,
                             model = normal_mri, ...) {
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
#' model_spec <- normal_mri()
#'
#' # Access model parameters
#' model_spec$p_types
normal_mri <- function(){
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
#' model_spec <- white_mri()
#'
#' # Access model parameters
#' model_spec$p_types

white_mri <- function(){
  return(
    list(
      type="MRI_white",
      c_name = "MRI_white",
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
#' # Create a simple design matrix
#' set.seed(123)
#' n_trs <- 200
#'
#' # Create a design matrix with two conditions and drift terms
#' design <- data.frame(
#'   condition1 = c(rep(0, 20), rep(1, 10), rep(0, 170)),
#'   condition2 = c(rep(0, 100), rep(1, 10), rep(0, 90)),
#'   drift1 = seq(0, 1, length.out = n_trs),
#'   drift2 = seq(1, 0, length.out = n_trs),
#'   derivative = c(rep(0, 10), rep(1, 5), rep(0, 185))
#' )
#'
#' # Basic plot of the design matrix
#' plot_design_fmri(design)
#'
#' # Plot specific events only
#' plot_design_fmri(design, events = c("condition1", "condition2"))
#'
#' # Include nuisance regressors
#' plot_design_fmri(design, remove_nuisance = FALSE)
#'
#' # Customize the plot
#' plot_design_fmri(design,
#'                 TRs = 150,
#'                 legend_pos = "topright",
#'                 main = "fMRI Design Matrix",
#'                 lwd = 3)
#'
#' # Example with an emc.design object
#' # First create a design matrix using convolve_design_matrix
#' timeseries <- data.frame(
#'   subjects = rep(1, 100),
#'   run = rep(1, 100),
#'   time = seq(0, 99),
#'   ROI1 = rnorm(100)
#' )
#'
#' events <- data.frame(
#'   subjects = rep(1, 4),
#'   run = rep(1, 4),
#'   onset = c(10, 30, 50, 70),
#'   duration = rep(0.5, 4),
#'   event_type = c("A", "B", "A", "B"),
#'   modulation = c(1, 1, 1, 1)
#' )
#'
#' # Create convolved design matrix
#' design_matrix <- convolve_design_matrix(
#'   timeseries = timeseries,
#'   events = events,
#'   factors = list(condition = c("A", "B")),
#'   hrf_model = "glover"
#' )
#'
#' # Create fMRI design for EMC2
#' fmri_design <- design_fmri(design_matrix, model = normal_mri)
#'
#' # Plot the design matrix
#' plot_design_fmri(fmri_design)
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
    design_matrix <- design_matrix[,!is_nuisance]
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


  design_matrix <- design_matrix[1:TRs, events]
  do.call(matplot, c(list(design_matrix, type = "l"), fix_dots_plot(dots)))
  do.call(legend, c(list(legend_pos, legend = events, bty = "n"), fix_dots(dots, legend)))
}

get_peri_stim_lines <- function(
    timeseries,
    events,
    event_type,
    pre = 2,
    post = 18,
    n_bins = 4,
    uniform_tol = 1e-5
) {
  if(!(all(c("run", "subjects", "time") %in% colnames(timeseries)))){
    stop("Expected columns in timeseries: run, subjects, time")
  }
  timeseries <- timeseries[,colnames(timeseries) != "postn"]
  if(ncol(timeseries) > 4) stop("Can only supply one ROI")
  ROI_col <- colnames(timeseries)[!colnames(timeseries) %in% c("run", "subjects", "time")]
  # -------------------------------------------------------------------
  # 1) Prepare the events of interest (filter + bin/categorize modulation)
  # -------------------------------------------------------------------
  ev_sub <- prepare_event_groups(events, event_type, n_bins)
  if (nrow(ev_sub) == 0) {
    stop("No events found for event_type = ", event_type)
  }

  # -------------------------------------------------------------------
  # 2) Pre-split timeseries by subject-run, gather chunk metadata
  # -------------------------------------------------------------------
  chunk_info_list <- split_timeseries_chunks(timeseries, ROI_col, uniform_tol)

  # -------------------------------------------------------------------
  # 3) For each bin or unique mod_group, compute the average snippet
  # -------------------------------------------------------------------
  groups <- sort(unique(ev_sub$mod_group))

  # We'll accumulate the final results as a data frame:
  # columns => mod_group, time, avg_signal, optional min_mod, max_mod
  output_list <- list()

  for (g in groups) {
    # events for this group
    ev_g <- ev_sub[ev_sub$mod_group == g, ]
    if (nrow(ev_g) == 0) next

    # Compute the average snippet (time axis + average)
    out_snip <- compute_avg_snippet_for_group(
      ev_g,
      chunk_info_list,
      ROI_col,
      pre,
      post
    )
    if (is.null(out_snip)) {
      # no valid snippet
      next
    }

    # We'll store them in a data frame with columns:
    # mod_group, time, avg_signal, (optional min_mod / max_mod)
    # If we used binning, we can store the min/max modulation range
    if (attr(ev_sub, "binned") == TRUE) {
      # we used binning => figure out the min/max within this group
      mod_in_bin <- ev_g$modulation
      min_mod <- round(min(mod_in_bin, na.rm = TRUE), 3)
      max_mod <- round(max(mod_in_bin, na.rm = TRUE), 3)
      df_tmp <- data.frame(
        mod_group = g,
        time      = out_snip$time,
        avg_signal= out_snip$avg,
        min_mod   = min_mod,
        max_mod   = max_mod
      )
    } else {
      # we used discrete categories => store just mod_group
      df_tmp <- data.frame(
        mod_group = g,
        time      = out_snip$time,
        avg_signal= out_snip$avg
      )
    }

    output_list[[as.character(g)]] <- df_tmp
  }

  # Combine into one big data frame
  final_df <- do.call(rbind, output_list)
  final_df <- final_df[final_df$time >= -pre & final_df$time <= post,]
  rownames(final_df) <- NULL
  return(final_df)
}

plot_fmri <- function(timeseries, post_predict = NULL, events, event_type, posterior_args = list(),
                      legend_pos = "topright", layout = NA, ...){
  posterior_args <- add_defaults(posterior_args, col = "darkgreen", lwd = 2)
  plot_args <- add_defaults(list(...), col = "black", lwd = 2, xlab = "time (s)", ylab = "BOLD response")
  ROI_col <- colnames(timeseries)[!colnames(timeseries) %in% c("run", "subjects", "time")]
  ts <- get_peri_stim_lines(timeseries, events, event_type = event_type)
  if(!is.null(post_predict)){
    pps <- split(post_predict, post_predict$postn)
    pps <- lapply(pps, get_peri_stim_lines, events, event_type = event_type)
    df_all <- do.call(rbind, pps)
    # Split the data by these grouping columns
    split_list <- split(df_all, df_all[c("mod_group", "time")])

    # 4) For each group, compute the 2.5, 50, 97.5 percentiles of avg_signal
    res_list <- lapply(split_list, function(subdf) {
      # subdf is all rows with the same (mod_group, time)
      qs <- quantile(subdf$avg_signal, c(0.025, 0.5, 0.975), na.rm = TRUE)
      # Build a small output row with just the grouping columns
      out <- subdf[1, c("mod_group", "time"), drop = FALSE]
      # If min_mod/max_mod exist, keep them from the first row
      out$min_mod <- subdf$min_mod[1]
      out$max_mod <- subdf$max_mod[1]
      # Add the quantiles
      out$p025 <- qs[1]
      out$p50  <- qs[2]
      out$p975 <- qs[3]
      out
    })

    # 5) Reassemble into one data frame
    res_df <- do.call(rbind, res_list)
    rownames(res_df) <- NULL
  }
  par(mfrow = c(2,2))
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
  for(mod in unique(ts$mod_group)){
    idx_ts <- ts$mod_group == mod
    tmp_plot_args <- add_defaults(plot_args, ylim = c(min(ts$avg_signal, res_df$p025), max(ts$avg_signal, res_df$p975)),
                                  main = paste0(ROI_col, " - ", event_type, ifelse(n_plots == 1, "", paste0(": ", insert, mod))))
    do.call(plot, c(list(ts$time[idx_ts], ts$avg_signal[idx_ts], type = "l"), fix_dots_plot(tmp_plot_args)))
    idx_pp <- res_df$mod_group == mod
    abline(h = 0, lty = 2)
    if(!is.null(post_predict)){
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


# -------------------------------------------------------------------------
# Helper 1: Filter events by event_type and bin/categorize the modulation
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
    ev_sub$mod_group <- cut(ev_sub$modulation, breaks = quartile_breaks, include.lowest = TRUE, labels = 1:4)
    attr(ev_sub, "binned") <- TRUE

  } else {
    # treat as categorical
    ev_sub$mod_group <- factor(ev_sub$modulation, levels = sort(mod_vals))
    attr(ev_sub, "binned") <- FALSE
  }
  return(ev_sub)
}

# -------------------------------------------------------------------------
# Helper 2: Split timeseries by (subjects, run) and gather metadata
# -------------------------------------------------------------------------
split_timeseries_chunks <- function(timeseries, ROI_col, uniform_tol = 1e-5) {
  # Split
  ts_split <- split(timeseries, list(timeseries$subjects, timeseries$run), drop = TRUE)

  # For each chunk, store:
  #   df          => sorted by time
  #   time_vec    => numeric vector of times
  #   dt          => median diff
  #   is_uniform  => whether the data are within uniform_tol of dt
  # Return a named list keyed by "subj.run"
  out_list <- lapply(ts_split, function(df_sub) {
    df_sub <- df_sub[order(df_sub$time), ]
    time_vec <- df_sub$time
    if (length(time_vec) < 2) {
      return(list(
        df = df_sub,
        time_vec = time_vec,
        dt = NA_real_,
        is_uniform = FALSE
      ))
    }
    diffs <- diff(time_vec)
    dt_median <- median(diffs)
    max_dev <- max(abs(diffs - dt_median))
    is_uni <- (max_dev < uniform_tol)
    list(
      df = df_sub,
      time_vec = time_vec,
      dt = dt_median,
      is_uniform = is_uni
    )
  })
  return(out_list)
}

# -------------------------------------------------------------------------
# Helper 3: Extract a single snippet from a chunk
# -------------------------------------------------------------------------
extract_snippet <- function(ci, onset, pre, post, ROI_col) {
  # ci => chunk info, with $df, $time_vec, $dt, $is_uniform
  # returns baseline-corrected snippet or NULL if out-of-range

  df_sub   <- ci$df
  time_vec <- ci$time_vec
  dt       <- ci$dt
  n_pts    <- nrow(df_sub)

  if (n_pts < 1) return(NULL)

  start_t <- onset - pre
  end_t   <- onset + post

  if (ci$is_uniform && !is.na(dt) && dt > 0) {
    # Direct index-based approach
    t0 <- time_vec[1]
    start_idx <- round((start_t - t0)/dt) + 1
    end_idx   <- round((end_t   - t0)/dt) + 1

    if (start_idx < 1 || end_idx < 1 ||
        start_idx > n_pts || end_idx > n_pts ||
        start_idx > end_idx) {
      return(NULL)
    }
    snippet_vals <- df_sub[[ROI_col]][start_idx:end_idx]

    # Baseline correction: subtract mean in [-pre, 0)
    zero_idx   <- round((onset - t0)/dt) + 1
    pre_start  <- start_idx
    pre_end    <- max(zero_idx - 1, start_idx)
    if (pre_end < pre_start) {
      baseline_mean <- 0
    } else {
      baseline_vals <- df_sub[[ROI_col]][pre_start:pre_end]
      baseline_mean <- mean(baseline_vals, na.rm = TRUE)
    }
    snippet_vals_bc <- snippet_vals - baseline_mean
    return(snippet_vals_bc)

  } else {
    # fallback with findInterval
    start_idx <- findInterval(start_t, time_vec)
    end_idx   <- findInterval(end_t,   time_vec)
    if (start_idx < 1) start_idx <- 1
    if (end_idx   < 1) return(NULL)
    if (end_idx   > n_pts) end_idx <- n_pts
    if (start_idx > n_pts || start_idx > end_idx) return(NULL)

    snippet_rows <- df_sub[start_idx:end_idx, ]
    # baseline in [-pre, 0)
    baseline_rows <- snippet_rows[snippet_rows$time < onset, ]
    baseline_mean <- if (nrow(baseline_rows) > 0) {
      mean(baseline_rows[[ROI_col]], na.rm = TRUE)
    } else 0
    snippet_vals_bc <- snippet_rows[[ROI_col]] - baseline_mean
    return(snippet_vals_bc)
  }
}

# -------------------------------------------------------------------------
# Helper 4: Compute the average snippet for one group of events
# -------------------------------------------------------------------------
compute_avg_snippet_for_group <- function(ev_g, chunk_info_list, ROI_col, pre, post) {
  # ev_g => subset of events for a single mod_group
  # For each event, we do extract_snippet(...).
  # Then we combine them. If uniform, they'll have the same length.
  # We return a list with $time (the x-axis) and $avg (the averaged signal).

  if (nrow(ev_g) < 1) return(NULL)

  snippet_list <- vector("list", nrow(ev_g))
  for (i in seq_len(nrow(ev_g))) {
    subj_i  <- ev_g$subjects[i]
    run_i   <- ev_g$run[i]
    onset_i <- ev_g$onset[i]

    chunk_key <- paste0(subj_i, ".", run_i)
    ci <- chunk_info_list[[chunk_key]]
    if (is.null(ci)) {
      snippet_list[[i]] <- NULL
      next
    }
    snippet_vals <- extract_snippet(ci, onset_i, pre, post, ROI_col)
    snippet_list[[i]] <- snippet_vals
  }

  # Among all snippets, find max length
  lens <- sapply(snippet_list, length)
  max_len <- max(lens, na.rm = TRUE)
  if (max_len < 1) {
    return(NULL)
  }

  # Build a matrix of shape (max_len x #events), fill with NA
  mat_snips <- matrix(NA, nrow = max_len, ncol = length(snippet_list))
  for (col_i in seq_along(snippet_list)) {
    vals <- snippet_list[[col_i]]
    if (!is.null(vals) && length(vals) > 0) {
      mat_snips[1:length(vals), col_i] <- vals
    }
  }
  avg_snip <- rowMeans(mat_snips, na.rm = TRUE)

  # approximate dt from all used chunks
  # if truly uniform, they might differ across runs, but we do a single approach:
  dts <- sapply(chunk_info_list, `[[`, "dt")
  dt_global <- median(dts, na.rm = TRUE)
  # define a time axis from -pre.. with steps dt_global
  # length = max_len
  time_axis <- seq(-pre, by = dt_global, length.out = max_len)

  list(time = time_axis, avg = avg_snip)
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
