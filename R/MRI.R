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

  design <- list(Flist = Flist, model = model)
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
      rfun=function(lR,pars) return(pars),
      # Density function (PDF) for single accumulator
      log_likelihood=function(pars, dadm, model, min_ll=log(1e-10)){
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
      rfun=function(lR,pars) return(pars),
      # Density function (PDF) for single accumulator
      log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
        # Extract observed data (as a vector)
        y <- as.vector(as.matrix(dadm[, !colnames(dadm) %in% c("subjects", "run", "time", "trials")]))
        n <- length(y)

        # Total number of parameter columns
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

        # Replace any log-likelihood values below min_ll with min_ll
        ll <- pmax(ll, min_ll)

        sum(ll)
      }
    )
  )
}





