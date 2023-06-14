## fmri functions
quick_convolve <- function(regressor, modulator, hkernel, frame_times) {
  ## minimizes overhead for quick convolution; useful when called from within a likelihood function (speed-up ~20 ms per run)
  regressor$regressor[regressor$regressor==1] <- regressor$regressor[regressor$regressor==1]*modulator

  conv_reg = NULL
  for(i in 1:ncol(hkernel)) {      # reverse kernel -- numpy implementation differs from R
    conv_reg = cbind(conv_reg, convolve(regressor$regressor, rev(hkernel[,i]), type='o')[1:length(regressor$regressor)])
  }

  computed_regressors = EMC2:::resample_regressor_(conv_reg, regressor$hr_frame_times, frame_times)
  computed_regressors = EMC2:::orthogonalize_(computed_regressors)
  return(computed_regressors)
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
    preds = (X[,i] %*% X[,1:(i-1)]) %*% pinv(as.matrix(X[,1:(i-1)]))
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
  time_stamps = linspace(0, time_length,
                         round(as.double(time_length) / dt))
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

  if(hrf_model == 'glover') {
    hkernel = matrix(glover_hrf(tr, oversampling), ncol=1)
  } else if(hrf_model == 'glover + derivative') {
    hkernel = cbind(glover_hrf(tr, oversampling),
                    glover_time_derivative(tr, oversampling))
  } else if(hrf_model == 'glover + derivative + dispersion') {
    hkernel = cbind(glover_hrf(tr, oversampling),
                    glover_time_derivative(tr, oversampling),
                    glover_dispersion_derivative(tr, oversampling))
  }
  return(hkernel)
}


sample_condition_ <- function(exp_condition, frame_times, oversampling=50,
                              min_onset=-24) {

  # Find the high-resolution frame_times
  n = length(frame_times)
  min_onset = as.double(min_onset)
  n_hr = ((n - 1) * 1. / (max(frame_times) - min(frame_times)) * (max(frame_times) * (1 + 1. / (n - 1)) - min(frame_times) - min_onset) * oversampling) + 1

  hr_frame_times = linspace(min(frame_times) + min_onset,
                            max(frame_times) * (1 + 1. / (n - 1)),
                            round(n_hr))

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
  for(i in 1:ncol(hkernel)) {      # reverse kernel -- numpy implementation differs from R
    conv_reg = cbind(conv_reg, convolve(hr_regressor, rev(hkernel[,i]), type='o')[1:length(hr_regressor)])
  }

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
  trial_type = events$trial_type
  onset = events$onset
  duration = events$duration
  modulation = events$modulation

  for(condition in sort(unique(trial_type))) {
    condition_mask = (trial_type == condition)
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
#'
#' @return
#' @export
#'
#' @examples
make_fmri_design_matrix <- function(frame_times, events=NULL, hrf_model='glover',
                                    drift_model='cosine',
                                    high_pass=0.01,
                                    drift_order=1,
                                    fir_delays=c(0),
                                    add_regs=NULL,
                                    add_reg_names=NULL,
                                    min_onset=-24,
                                    oversampling=50, add_intercept=TRUE) {
  names = c()
  matrix = NULL

  # step 1: events-related regressors
  if(!is.null(events)) {
    # create the condition-related regressors
    if(is.character(hrf_model)) {
      hrf_model = tolower(hrf_model)
      out = convolve_regressors_(
        events, hrf_model, frame_times, fir_delays, min_onset,
        oversampling)
      matrix = out[[1]]
      names = out[[2]]
    }
  }
  df <- data.frame(matrix, row.names=frame_times)
  colnames(df) <- names
  if(add_intercept) df$intercept <- 1
  return(df)
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
