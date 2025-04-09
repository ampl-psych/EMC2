
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
#' design_matrices <-  convolved_design_matrix(
#'   timeseries = ts,
#'   events = events,
#'   factors = list(difficulty = c("hard", "easy")),
#'   contrasts = list(difficulty = matrix(c(-1, 1)))
#' )
convolved_design_matrix <- function(timeseries, events, factors = NULL, contrasts = NULL,
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


# High pass filtering -----------------------------------------------------
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
#' This function takes the output from convolved_design_matrix and transforms it into a design
#' suitable for sampling with EMC2. It properly configures parameter types, bounds, and transformations
#' for the specified model.
#'
#' @param design_matrix A list of design matrices,  the output from convolved_design_matrix
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
#' design_matrix <- convolved_design_matrix(
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
        # Center the predicted means (as in likelihood function)
        y_hat <- y_hat - mean(y_hat)
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
        y_hat <- y_hat - mean(y_hat)
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



plot_hrf_shape <- function(timeseries, events,
                     factors,
                     contrasts = NULL,
                     cell_coding = FALSE,
                     hrf_model = "glover + derivative",
                     min_onset = -24,
                     oversampling = 50,
                     factor_name = names(factors)[1],
                     aggregate = TRUE,
                     epoch_duration = 32) {
  # Extract the frame times from the timeseries.
  frame_times <- sort(unique(timeseries$time))

  # Subset the events to only those that belong to the chosen factor.
  ev_sub <- events[events$event_type %in% factors[[factor_name]], ]
  # Tag these events with the factor name and rename the event type column.
  ev_sub$factor <- factor_name
  colnames(ev_sub)[colnames(ev_sub) == "event_type"] <- factor_name

  # Apply contrasts if a contrast for this factor is provided.
  ev_proc <- apply_contrasts(ev_sub,
                             contrast = contrasts[[factor_name]],
                             cell_coding = factor_name %in% cell_coding)

  # ev_proc is in "long" format and should include a column 'regressor'
  reg_levels <- unique(ev_proc$regressor)

  # Set up a multi-panel plot (here we use 2 columns)
  n_panels <- length(reg_levels)
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(ceiling(n_panels/2), 2))

  # For each level (i.e. each regressor), compute and plot the convolved HRF.
  for (lev in reg_levels) {
    ev_lev <- ev_proc[ev_proc$regressor == lev, ]
    # Make sure we have the key columns: onset, duration, modulation.
    exp_condition <- as.matrix(ev_lev[, c("onset", "duration", "modulation")])

    # Compute the convolved regressor for this condition.
    reg_list <- compute_convolved_regressor(exp_condition,
                                            hrf_model,
                                            frame_times,
                                            con_id = as.character(lev),
                                            oversampling = oversampling,
                                            min_onset = min_onset)
    hrf_vector <- reg_list$computed_regressors[, 1]  # canonical HRF

    if (aggregate) {
      # Aggregate across events: extract epochs around each event onset and average.
      event_onsets <- ev_lev$onset
      # Determine the TR (assumes frame_times are equally spaced).
      tr <- diff(frame_times)[1]
      n_timepoints <- round(epoch_duration / tr)
      segments <- matrix(NA, nrow = length(event_onsets), ncol = n_timepoints)
      for (i in seq_along(event_onsets)) {
        onset_time <- event_onsets[i]
        # Find the index corresponding to the event onset.
        idx <- which.min(abs(frame_times - onset_time))
        # Make sure we have enough points after the event.
        if (idx + n_timepoints - 1 <= length(frame_times)) {
          segments[i, ] <- hrf_vector[idx:(idx + n_timepoints - 1)]
        }
      }
      # Compute the average HRF (ignoring any incomplete epochs).
      avg_hrf <- colMeans(segments, na.rm = TRUE)
      time_axis <- seq(0, epoch_duration, length.out = n_timepoints)

      # Plot the aggregated average HRF.
      plot(time_axis, avg_hrf, type = "l", lwd = 2, col = "blue",
           xlab = "Time (s) from event onset", ylab = "Average HRF amplitude",
           main = paste("Factor:", factor_name, "\nLevel:", lev, "\nAggregated HRF"))
    } else {
      # Plot the full computed regressor over the entire timeseries.
      plot(frame_times, hrf_vector, type = "l", lwd = 2, col = "blue",
           xlab = "Time (s)", ylab = "HRF amplitude",
           main = paste("Factor:", factor_name, "\nLevel:", lev))
    }
  }

  # Reset the plotting parameters.
  par(old_par)
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
  attr(data, "designs") <- design$fMRI_design[[data$subjects[1]]]
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
