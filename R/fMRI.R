## fmri functions
quick_convolve <- function(regressor, modulator, hkernel, frame_times) {
  ## minimizes overhead for quick convolution; useful when called from within a likelihood function (speed-up ~20 ms per run)
  regressor$regressor[regressor$regressor==1] <- regressor$regressor[regressor$regressor==1]*modulator

  conv_reg = NULL
  for(i in 1:ncol(hkernel)) {      # reverse kernel -- numpy implementation differs from R
    conv_reg = cbind(conv_reg, convolve(regressor$regressor, rev(hkernel[,i]), type='o')[1:length(regressor$regressor)])
  }

  computed_regressors = resample_regressor_(conv_reg, regressor$hr_frame_times, frame_times)
  computed_regressors = orthogonalize_(computed_regressors)
  return(computed_regressors)
}

make_fmri_design_matrix_wrap <- function(timeseries, events, factors, contrasts){
  if(!is.data.frame(events)) events <- do.call(rbind, events)
  subject <- timeseries$subjects[1]
  runs <- unique(timeseries$run)
  events <- events[events$trial_type %in% unlist(factors),]
  events$factor <- events$trial_type
  events$modulation <- 1
  events$regressor <- events$trial_type
  out <- data.frame
  for(fact in names(factors)){
    idx <- events$trial_type %in% factors[[fact]]
    events$factor[idx] <- fact
    tmp <- events[idx,]
    dm <- contrasts[[fact]]
    tmp2 <- do.call(rbind, rep(list(tmp), ncol(dm)))
    tmp2 <- tmp2[order(tmp2$onset),]
    for(type in unique(tmp2$trial_type)){
      tmp2$modulation[tmp2$trial_type == type] <- dm[rownames(dm) == type]
      tmp2$regressor[tmp2$trial_type == type] <- colnames(dm)
    }
    tmp2 <- tmp2[!tmp2$modulation == 0,]
    events <- events[!idx,]
    events <- rbind(events, tmp2)
  }
  events <- events[order(events$onset),]
  runs <- unique(timeseries$run)
  dms <- vector(mode='list', length=length(runs))
  for(run in runs) {
    dms[[run]] <- make_fmri_design_matrix(timeseries[timeseries$run==run, 'time'],
                                          events=events[events$run == run,],
                                          hrf_model='glover + derivative', add_intercept=FALSE)
  }
  dm <- do.call(rbind, dms)
  dm$subjects <- subject
  return(dm)
}

make_design_fmri <- function(data,
                             events,
                             model,
                             factors,
                             contrasts,
                             hrf_model='glover + derivative',
                             add_intercept=FALSE,
                             ...) {

  design_matrix <- make_fmri_design_matrix_wrap(data, events, factors, contrasts)

  data_names <- colnames(data)[!colnames(data) %in% c('subjects', 'run', 'time')]

  # First create the design matrix

  # generate parameter names
  if(!is.null(design_matrix)) {
    par_names <- colnames(design_matrix)[colnames(design_matrix) != "subjects"]
  } else {
    par_names <- c()
  }
  if(!is.null(list(...)$par_regressors)) {
    if(hrf_model == 'glover + derivative') {
      par_names_from_events <- unique(events$trial_type)  # events contains all events
      par_names <- c(par_names, paste(rep(par_names_from_events, each=2), c('', '_derivative'), sep=''))
    } else if(hrf_model == 'glover') {
      par_names_from_events <- unique(events$trial_type)  # events contains all events
      par_names <- c(par_names, par_names_from_events)
    }
  }
  if(!('intercept' %in% par_names) & add_intercept) {
    par_names <- c(par_names, 'intercept')
  }

  if(!is.null(list(...)$par_regressors)) events_list <- split(list(...)$par_regressors, f=list(...)$par_regressors$subjects)

  dms_list <- split(design_matrix, f=design_matrix$subjects)
  model <- model()
  if(!is.null(list(...)$par_regressors)) {
    model$list(...)$par_regressors <- lapply(events_list, FUN = function(x){
      y <- x[,colnames(x) != "subjects"]
      y
    })
  }


  dots <- list(...)
  if('hkernel' %in% names(dots)) {
    model$hkernel <- dots$hkernel
  }
  if('regressors' %in% names(dots)) {
    model$regressors <- dots$regressors
    model$fast_convolve <- TRUE
  }


  df_par_names <- expand.grid(c(par_names, "sd"), data_names)
  par_names <- paste0(df_par_names[,2], "_", df_par_names[,1])

  n_pars <- length(par_names)
  Flist <- vector("list", n_pars)
  for(i in 1:n_pars){
    Flist[[i]] <- as.formula(paste0(par_names[i], "~1"))
  }
  model_function <- function() {return(model)}
  design <- list(Flist = Flist, model = model_function)
  attr(design, "design_matrix") <- lapply(dms_list, FUN=function(x) {
    y <- x[,colnames(x) != 'subjects']
    data.matrix(y)
  })
  attr(design, "p_vector") <- par_names
  return(design)
}



searchsorted <- function(x, insert) {
  ## SM, based on https://numpy.org/doc/stable/reference/generated/numpy.searchsorted.html
  return(findInterval(insert,x, left.open = TRUE)+1)
}

orthogonalize_ <- function(X) {
  if(is.null(dim(X))) {
    return(X)
  }
  if(dim(X)[2] == 1) {
    return(X)
  }

  for(i in 2:dim(X)[2]) {
    preds = (X[,i] %*% X[,1:(i-1)]) %*% ginv(as.matrix(X[,1:(i-1)]))
    X[,i] = X[,i]-preds
  }
  return(X)
}

regressor_names <- function(con_name, hrf_model, fir_delays=NULL) {
  names <- c(con_name)

  if(hrf_model %in% c('glover', 'spm')) {
    names = con_name
  } else if(hrf_model %in% c("glover + derivative", 'spm + derivative')) {
    names = c(con_name, paste0(con_name, "_derivative"))
  } else if(hrf_model %in% c('spm + derivative + dispersion',
                             'glover + derivative + dispersion')) {
    names = c(con_name, paste0(con_name, "_derivative"), paste0(con_name, "_dispersion"))
  }

  return(names)
}

resample_regressor_ <- function(hr_regressor, hr_frame_times, frame_times) {
  #  from scipy.interpolate import interp1d
  if(is.null(dim(hr_regressor))) {
    f = approxfun(hr_frame_times, hr_regressor)
    return(f(frame_times))
  } else {
    out <- NULL
    for(c in 1:ncol(hr_regressor)) {
      f = approxfun(hr_frame_times, hr_regressor[,c])
      out <- cbind(out, f(frame_times))
    }
    return(out)
  }
}

glover_hrf <- function(tr, oversampling=50, time_length=32., onset=0) {
  return(gamma_difference_hrf_(tr, oversampling, time_length, onset,
                               delay=6, undershoot=12., dispersion=.9,
                               u_dispersion=.9, ratio=.35))
}

gamma_difference_hrf_ <- function(tr, oversampling=50, time_length=32., onset=0.,
                                  delay=6, undershoot=16., dispersion=1.,
                                  u_dispersion=1., ratio=0.167) {
  dt = tr / oversampling
  time_stamps = seq(0, time_length,
                         length.out = round(as.double(time_length) / dt))
  time_stamps = time_stamps - onset

  # define peak and undershoot gamma functions
  peak_gamma = dgamma(x=(time_stamps-dt)/dispersion, delay/dispersion)
  # peak_gamma = gamma.pdf(
  #   time_stamps,
  #   delay / dispersion,
  #   loc=dt,
  #   scale=dispersion)
  undershoot_gamma = dgamma(x=(time_stamps-dt)/u_dispersion, undershoot/u_dispersion)
  # undershoot_gamma = gamma.pdf(
  #   time_stamps,
  #   undershoot / u_dispersion,
  #   loc=dt,
  #   scale=u_dispersion)

  # calculate the hrf
  hrf = peak_gamma - ratio * undershoot_gamma
  #  hrf = hrf / sum(hrf)
  hrf = hrf / max(hrf)

  return(hrf)
}

glover_time_derivative <- function(tr, oversampling=50, time_length=32., onset=0.) {
  do = .1
  dhrf = 1. / do * (
    glover_hrf(tr, oversampling, time_length, onset)- glover_hrf(tr, oversampling, time_length, onset + do)
  )
  return(dhrf)
}


hrf_kernel <- function(hrf_model, tr, oversampling=50, fir_delays=NULL) {
  acceptable_hrfs = c(
    #    'spm', 'spm + derivative', 'spm + derivative + dispersion',
    #    'fir',
    'glover', 'glover + derivative', 'glover + derivative + dispersion')

  if(!hrf_model %in% acceptable_hrfs) {
    stop('HRF model not supported.')
  }

  #  error_msg = ("Could not process custom HRF model provided.Please refer to the related documentation.")
  hrf_kernel <- NULL
  # if(hrf_model == 'glover') {
  #   hkernel = matrix(glover_hrf(tr, oversampling), ncol=1)
  # } else if(hrf_model == 'glover + derivative') {
  #   hkernel = cbind(glover_hrf(tr, oversampling),
  #                   glover_time_derivative(tr, oversampling))
  # } else if(hrf_model == 'glover + derivative + dispersion') {
  #   hkernel = cbind(glover_hrf(tr, oversampling),
  #                   glover_time_derivative(tr, oversampling),
  #                   glover_dispersion_derivative(tr, oversampling))
  # }
  return(hkernel)
}


sample_condition_ <- function(exp_condition, frame_times, oversampling=50,
                              min_onset=-24) {

  # Find the high-resolution frame_times
  n = length(frame_times)
  min_onset = as.double(min_onset)
  n_hr = ((n - 1) * 1. / (max(frame_times) - min(frame_times)) * (max(frame_times) * (1 + 1. / (n - 1)) - min(frame_times) - min_onset) * oversampling) + 1

  hr_frame_times = seq(min(frame_times) + min_onset,
                            max(frame_times) * (1 + 1. / (n - 1)),
                            length.out = round(n_hr))

  # Get the condition information
  #onsets, durations, values = tuple(map(np.asanyarray, exp_condition))
  onsets = exp_condition[,1]
  durations = exp_condition[,2]
  values = exp_condition[,3]
  #  if (onsets < frame_times[0] + min_onset).any():
  #    warnings.warn(('Some stimulus onsets are earlier than %s in the'
  #                   ' experiment and are thus not considered in the model'
  #                   % (frame_times[0] + min_onset)), UserWarning)

  # Set up the regressor timecourse
  tmax = length(hr_frame_times)
  regressor = matrix(0, ncol=length(hr_frame_times))

  t_onset = pmin(searchsorted(hr_frame_times, onsets), tmax - 1)
  for(i in 1:length(t_onset)) {
    t = t_onset[i]
    v = values[i]
    regressor[t] = regressor[t] + v
  }
  t_offset = pmin(searchsorted(hr_frame_times, onsets + durations), tmax - 1)

  # Handle the case where duration is 0 by offsetting at t + 1
  for(i in 1:length(t_offset)) {
    t <- t_offset[i]
    if((t < (tmax - 1)) & (t == t_onset[i])) {
      t_offset[i] = t_offset[i] + 1
    }
  }

  for(i in 1:length(t_offset)) {
    t = t_offset[i]
    v = values[i]
    regressor[t] = regressor[t]-v
  }
  regressor = cumsum(regressor)

  return(list('regressor'=regressor, 'hr_frame_times'=hr_frame_times))
}


compute_regressor <- function(exp_condition, hrf_model, frame_times, con_id='cond',
                              oversampling=50, fir_delays=NULL, min_onset=-24) {

  # fir_delays should be integers
  if(!is.null(fir_delays)) {
    for(i in 1:length(fir_delays)) {
      fir_delays[i] = as.integer(fir_delays[i])
    }
  }
  oversampling = as.integer(oversampling)

  # this is the average tr in this session, not necessarily the true tr
  tr = max(frame_times) / (length(frame_times) - 1)

  # 1. create the high temporal resolution regressor
  out = sample_condition_(
    exp_condition, frame_times, oversampling, min_onset)
  hr_regressor = out[['regressor']]
  hr_frame_times = out[['hr_frame_times']]

  # 2. create the  hrf model(s)
  hkernel = hrf_kernel(hrf_model, tr, oversampling, fir_delays)

  # 3. convolve the regressor and hrf, and downsample the regressor
  conv_reg = NULL
  # for(i in 1:ncol(hkernel)) {      # reverse kernel -- numpy implementation differs from R
  #   conv_reg = cbind(conv_reg, convolve(hr_regressor, rev(hkernel[,i]), type='o')[1:length(hr_regressor)])
  # }

  # 4. temporally resample the regressors
  #  if hrf_model == 'fir' and oversampling > 1:
  #   computed_regressors = _resample_regressor(
  #     conv_reg[:, oversampling - 1:],
  #     hr_frame_times[: 1 - oversampling],
  #     frame_times)
  # else:
  computed_regressors = resample_regressor_(conv_reg, hr_frame_times, frame_times)
  #
  # # 5. ortogonalize the regressors
  # if hrf_model != 'fir':
  computed_regressors = orthogonalize_(computed_regressors)
  #
  # # 6 generate regressor names
  reg_names = regressor_names(con_id, hrf_model, fir_delays=fir_delays)
  return(list(computed_regressors, reg_names))
}

check_events <- function(events) {
  if(!'modulation' %in% colnames(events)) {
    events['modulation'] = 1
  }
  return(events)
}

convolve_regressors_ <- function(events, hrf_model, frame_times, fir_delays=c(0),
                                 min_onset=-24, oversampling=50) {
  regressor_names = c()
  regressor_matrix = NULL

  events = check_events(events)
  regressor = events$regressor
  onset = events$onset
  duration = events$duration
  modulation = events$modulation

  for(condition in sort(unique(regressor))) {
    condition_mask = (regressor == condition)
    exp_condition = cbind(onset[condition_mask],
                          duration[condition_mask],
                          modulation[condition_mask])
    out = compute_regressor(
      exp_condition, hrf_model, frame_times, con_id=condition,
      fir_delays=fir_delays, oversampling=oversampling,
      min_onset=min_onset)
    reg = out[[1]]
    names = out[[2]]

    regressor_names = c(regressor_names, names)

    if(is.null(regressor_matrix)) {
      regressor_matrix = reg
    } else {
      regressor_matrix = cbind(regressor_matrix, reg)
    }
  }
  return(list(regressor_matrix, regressor_names))
}


# Function to create fMRI design matrix from a few basics
#' Function to create fMRI design matrix from a few basics
make_fmri_design_matrix <- function(frame_times, events=NULL, hrf_model='glover',
                                    drift_model='cosine',
                                    high_pass=0.01,
                                    drift_order=1,
                                    fir_delays=c(0),
                                    add_regs=NULL,
                                    add_reg_names=NULL,
                                    ROIs,
                                    min_onset=-24,
                                    oversampling=50, add_intercept=TRUE) {
  names_out <- c()
  matrix_out <- NULL

  # step 1: events-related regressors
  if(!is.null(events)) {
    # create the condition-related regressors
    if(is.character(hrf_model)) {
      hrf_model = tolower(hrf_model)
      out = convolve_regressors_(
        events, hrf_model, frame_times, fir_delays, min_onset,
        oversampling)
      matrix_out <- out[[1]]
      names_out <- out[[2]]
    }
  }
  df <- data.frame(matrix_out, row.names=frame_times)
  colnames(df) <- names_out
  if(add_intercept) df$intercept <- 1
  return(df)
}

normal <- function(){
  return(
    list(
      type="MRI",
      c_name = "MRI",
      p_types=c("sd" = log(1)), #This is a bit hacky for now
      # Transform to natural scale
      Ntransform=function(x) {

        if(is.null(dim(x))) {
          is_sd <- grepl("_sd", names(x))
          x[is_sd] <- exp(x[is_sd])+.001
        } else {
          is_sd <- grepl("_sd", dimnames(x)[[2]])
          x[,is_sd] <- exp(x[,is_sd])+.001
        }
        return(x)
      },
      # Trial dependent parameter transform
      Ttransform = function(pars,dadm=NULL)
      {
        pars
      },
      # p_vector transform
      transform = function(x) x,
      # Random function for racing accumulators
      rfun=function(lR,pars) rNORMAL(lR,pars),
      # Density function (PDF) for single accumulator
      dfun=function(rt,pars) dNORMAL(rt,pars),
      # Probability function (CDF) for single accumulator
      pfun=function(rt,pars) pNORMAL(rt,pars),
      # Race likelihood combining pfun and dfun
      log_likelihood=function(p_vector, dadm, predictionErrors=NULL, min_ll=log(1e-10), X2=NULL){
        # data
        y <- as.matrix(dadm[,!colnames(dadm) %in% c("subjects", 'run', 'time', "trials")])

        # first part of the design matrix is already generated, and fixed across parameters
        X <- attr(dadm, 'design_matrix_mri')

        #         if(!is.null(predictionErrors)) {
        #           events <- attr(dadm, "model")()$events[[as.character(dadm$subjects[1])]]
        #           events <- merge(events, predictionErrors, on='trial_nr')
        #           events$modulation <- (events$predictionErrors - mean(events$predictionErrors)) / sd(events$predictionErrors)
        #           runs <- unique(events$run)
        #           X2 <- as.matrix(do.call(rbind, lapply(runs, function(x) {
        #             make_fmri_design_matrix(frame_times=dadm[dadm$run==x,'time'],
        #                                     events=events[events$run==x,],
        #                                     hrf_model = 'glover + derivative',
        #                                     add_intercept=FALSE, oversampling = 5)
        #           })))
        #
        #           # combine 'fixed' and modulation parts of the DM
        #           X <- cbind(X, X2)
        #         }

        # transform parameters
        p_vector <- normal()$Ntransform(p_vector)
        X <- cbind(X, X2)
        # grab the right parameters
        is_sd <- grepl("sd", names(p_vector))
        sigma <- p_vector[is_sd]
        betas <- p_vector[!is_sd]
        betas <- matrix(betas, ncol = length(sigma))

        # get rid of intercept
        y_hat <- X %*% betas
        y_hat <- y_hat - mean(y_hat)
        ## SM style
        if(ncol(betas > 1)) {
          total_sum <- sum(pmax(dnorm(as.matrix(y), mean = y_hat, sd = rep(sigma, each=nrow(X)), log = T), min_ll))
        } else {
          total_sum <- sum(pmax(dnorm(y, mean = y_hat, sd = sigma, log = T), min_ll))
        }

        return(total_sum)
      }
    )
  )
}


# Script to run this ------------------------------------------------------
# frame_times = seq(0, 100, 1.38)
# events = data.frame('onset'=linspace(0, 50, 20), 'trial_type'=c(rep('cue',10),rep('stim',10)), 'duration'=1)
# events['modulation'] = rnorm(n=20)
# frame_times
#
#X <- make_fmri_design_matrix(frame_times=frame_times, events=events, hrf_model = 'glover + derivative')
# plot(X[,1], type='l')
# lines(X[,2], col=2)
# lines(X[,3], type='l')
# lines(X[,4], col=2)
#
