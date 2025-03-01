
# Design matrix functions -------------------------------------------------

apply_contrasts <- function(events, contrast = NULL, cell_coding = FALSE, remove_intercept = FALSE, do_print = FALSE) {
  factor_name <- events$factor[1]
  colnames(events)[colnames(events) == "event_type"] <- factor_name

  # If a contrast is provided, use it; otherwise let R default to its default contrasts.
  if(!is.null(contrast)){
    if(is.matrix(contrast)){
      events[[factor_name]] <- factor(events[[factor_name]], levels = rownames(contrast))
      stats::contrasts(events[[factor_name]], how.many = ncol(contrast)) <- contrast
    } else {
      events[[factor_name]] <- factor(events[[factor_name]])
      stats::contrasts(events[[factor_name]]) <- do.call(contrast, list(n = length(unique(events[[factor_name]]))))
    }
  } else {
    events[[factor_name]] <- factor(events[[factor_name]])
    # R's default contrasts will be used.
  }

  if(cell_coding){
    design <- model.matrix(as.formula(paste0("~ 0 + ", factor_name)), events)
  } else {
    design <- model.matrix(as.formula(paste0("~ ", factor_name)), events)
    colnames(design)[1] <- paste0(factor_name, "0")
    if(remove_intercept) design <- design[, -1, drop = FALSE]
  }

  events_design <- cbind(events, design)
  if(do_print){
    print(unique(events_design[, !colnames(events_design) %in% c("onset", "subjects", "duration")]))
  }

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

# Build the full design matrices for all subjects and runs.
build_fmri_design_matrices <- function(timeseries, events, factors, contrasts,
                                       hrf_model = 'glover + derivative', cell_coding = NULL) {
  events$duration <- 0.01  # default duration
  if(!is.data.frame(events)) events <- do.call(rbind, events)
  subjects <- unique(timeseries$subjects)
  all_dms <- list()

  for(subject in subjects){
    ev_sub <- events[events$subjects == subject, ]
    ts_sub <- timeseries[timeseries$subjects == subject, ]
    runs <- unique(ts_sub$run)
    dms_sub <- vector("list", length = length(runs))

    for(run in runs){
      ev_run <- ev_sub[ev_sub$run == run, ]
      ts_run <- ts_sub[ts_sub$run == run, ]
      ev_tmp <- data.frame()

      for(fact in names(factors)){
        idx <- ev_run$event_type %in% factors[[fact]]
        ev_run$factor[idx] <- fact
        tmp <- ev_run[idx, ]
        ev_tmp <- rbind(ev_tmp, apply_contrasts(tmp, contrast = contrasts[[fact]],
                                                do_print = (run == runs[1]) & (subject == subjects[1]),
                                                cell_coding = fact %in% cell_coding))
      }

      ev_tmp <- ev_tmp[order(ev_tmp$onset), ]
      dms_sub[[as.character(run)]] <- construct_design_matrix(ts_run$time,
                                                              events = ev_tmp,
                                                              hrf_model = hrf_model,
                                                              min_onset = -24,
                                                              oversampling = 50,
                                                              add_intercept = FALSE)
    }

    dm_sub <- do.call(rbind, dms_sub)
    dm_sub$subjects <- subject
    all_dms[[as.character(subject)]] <- dm_sub
  }

  return(all_dms)
}


# Models ------------------------------------------------------------------



make_design_fmri <- function(data,
                             events,
                             model = normal_mri,
                             factors,
                             contrasts = NULL,
                             hrf_model='glover + derivative',
                             cell_coding = NULL,
                             ...) {
  dots <- list(...)
  if(!is.null(cell_coding) && !cell_coding %in% names(factors)) stop("Cell coded factors must have same name as factors argument")
  dots <- list(...)
  dots$add_intercept <- FALSE
  design_matrix <- build_fmri_design_matrices(data, events, factors, contrasts, hrf_model, cell_coding)
  data_names <- colnames(data)[!colnames(data) %in% c('subjects', 'run', 'time')]
  # First create the design matrix

  betas <- colnames(design_matrix[[1]])[colnames(design_matrix[[1]]) != "subjects"]


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

  design <- list(Flist = Flist, model = model, Ffactors = list(subjects = unique(data$subjects)))
  attr(design, "design_matrix") <- lapply(design_matrix, FUN=function(x) {
    y <- x[,colnames(x) != 'subjects']
    data.matrix(y)
  })
  par_names <- setNames(numeric(length(par_names)), par_names)
  attr(design, "p_vector") <- par_names
  return(design = design)
}

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
        y_hat <- y_hat - mean(y_hat)

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

        # Extract parameters:
        # - betas: columns 1 to (m-2)
        # - rho: column (m-1)
        # - sigma: column m (stationary standard deviation)
        betas <- pars[, 1:(m - 2), drop = FALSE]
        rho   <- pars[, m - 1]
        sigma <- pars[, m]

        # Compute predicted values and demean them
        y_hat <- rowSums(betas)
        y_hat <- y_hat - mean(y_hat)

        # Log-likelihood for the first observation
        ll <- numeric(n)
        ll[1] <- dnorm(y[1], mean = y_hat[1], sd = sigma[1], log = TRUE)

        # For observations t = 2:n, compute conditional means and sds vectorized:
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



plot_hrf <- function(timeseries, events,
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

