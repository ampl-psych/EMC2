#' @export
print.emc <- function(x, ...){
  n_chain <- chain_n(x)
  cat("Iterations: \n")
  print(chain_n(x))
  sub_names <- names(x[[1]][["data"]])
  cat("\n")
  cat("Subjects: \n")
  print(sub_names)
  par_names <- x[[1]][["par_names"]]
  cat("\n")
  cat("Parameters: \n")
  print(par_names)
  result <- list(
    iterations = n_chain, subjects = sub_names, parameters = par_names
  )
  return(invisible(result))
}

#' Summary Statistics for emc Objects
#'
#' Computes quantiles, `Rhat` and `ESS` for selected model parameters.
#'
#' Note that if `selection = alpha` and `by_subject = TRUE` (default)
#' is used, summary statistics are computed at the individual level.
#" Only the first subject's summary output is printed
#' to the console but summary statistics for all subjects are returned by the
#' function.
#'
#' If an emc object that has not been run with `fit` yet is supplied, summary of
#' the design will be returned.
#'
#' @param object An object of class `emc`
#' @param selection A character string indicating the parameter type
#' Defaults to `mu`, `sigma2`, and `alpha`. See below for more information.
#' @param probs The quantiles to be computed. Defaults to the the 2.5%, 50% and 97.5% quantiles.
#' @param digits An integer specifying rounding of output.
#' @param ... Optional arguments that can be passed to `get_pars`
#'
#' @return A list of summary output.
#' @export
summary.emc <- function(
    object,
    selection = c("mu", "sigma2", "alpha"),
    probs = c(0.025, .5, 0.975),
    digits = 3,
    ...
) {
  dots <- list(...)
  dots <- add_defaults(dots, by_subject = TRUE)
  if(object[[1]][["type"]] == "single") {
    selection <- "alpha"
  }
  if(!object[[1]][["init"]]) {
    warning("emc object has not been run with `fit` yet, design summary is returned.")
    summary(get_design(object))
    return(invisible(get_design(object)))
  }
  out_list <- list()
  for(select in selection) {
    stats <- do.call(
      get_summary_stat,
      c(
        list(
          object,
          fun = c(get_posterior_quantiles, r_hat, n_eff),
          probs = probs,
          stat_name = c(paste0(probs*100, "%"), "Rhat", "ESS"),
          selection = select
        ),
        dots
      )
    )
    for (i in seq_along(stats)) {
      stat <- round(stats[[i]], digits)
      stat[ , ncol(stat)] <- round(stat[ , ncol(stat)])
      if (select == "alpha" & dots[["by_subject"]]) {
        if (i == 1) {
          cat("\n", paste0(select, " ", names(stats)[i]), "\n")
          print(stat)
        }
      } else {
        if(length(stats) > 1) {
          cat("\n", paste0(select, " ", names(stats)[i]), "\n")
        } else {
          cat("\n", names(stats)[i], "\n")
        }
        print(stat)
      }
      out_list[[names(stats)[i]]] <- stat
    }
  }
  return(invisible(out_list))
}


#' Plot Function for emc Objects
#'
#' Makes trace plots for model parameters.
#'
#' If an emc object that has not been run with `fit` yet is supplied
#' prior plots will be returned.
#'
#' @param x An object of class `emc`
#' @param stage A character string indicating the sampling stage to be summarized.
#' Can be `preburn`, `burn`, `adapt`, or `sample`.
#' @param selection A character vector indicating the parameter group(s).
#' Defaults to `mu`, `sigma2`, and `alpha`.
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param ... Optional arguments that can be passed to `get_pars` or `plot.default` (see `par()`)
#' @return A trace/acf plot of the selected MCMC chains
#' @export
#'
#' @examples
#' plot(samples_LNR)
#' # Or trace autocorrelation for the second subject:
#' plot(samples_LNR, subject = 2, selection = "alpha")
#'
#' # Can also plot the trace of for example the group-level correlation:
#' plot(samples_LNR, selection = "correlation", col = c("green", "purple", "orange"), lwd = 2)
plot.emc <- function(
    x,
    stage = "sample",
    selection = c("mu", "sigma2", "alpha"),
    layout = NA,
    ...
) {
  if (!x[[1]][["init"]]) {
    warning("emc object has not been run with `fit` yet, prior plots are returned.")
    plot(get_prior(x))
    return(invisible(prior(x)))
  }
  if (x[[1]][["type"]] == "single" && selection[1] != "LL") {
    selection <- "alpha"
  }
  for (select in selection) {
    plot_mcmc_list(x, stage = stage, selection = select, layout = layout, ...)
  }
}



#' Generate Posterior/Prior Predictives
#'
#' Simulate ``n_post`` data sets using the posterior/prior parameter estimates
#'
#' @param object An emc or emc.prior object from which to generate predictives
#' @param hyper Boolean. Defaults to `FALSE`. If `TRUE`, simulates from the group-level (`hyper`)
#' parameters instead of the subject-level parameters.
#' @param n_post Integer. Number of generated datasets
#' @param n_cores Integer. Number of cores across which there should be parallellized
#' @param stat Character. Can be `mean`, `median` or `random` (i.e., the default).
#' Will take either random samples from the chain(s) or use the mean or median of the parameter estimates.
#' @param ... Optional additional arguments passed to `get_pars` or `make_data`
#' @return A list of simulated data sets of length `n_post`
#' @examples \donttest{
#' # based on an emc object ran by fit() we can generate posterior predictives
#' predict(samples_LNR, n_cores = 1, n_post = 10)
#' }
#' @export
predict.emc <- function(
    object,
    hyper = FALSE,
    n_post = 50,
    n_cores = 1,
    stat = c("random", "mean", "median"),
    ...
) {
  # #' @param force_direction Boolean, take censor direction from argument not samples (default FALSE)
  # #' @param force_response Boolean, take censor response from argument not samples (default FALSE)
  # #' @param LCresponse Boolean, default TRUE, if false set LC response to NA
  # #' @param UCresponse Boolean, default TRUE, if false set UC response to NA
  # #' @param LCdirection Boolean, default TRUE, set LC rt to 0, else to NA
  # #' @param UCdirection Boolean, default TRUE, set LC rt to Inf, else to NA
  # #' @param expand Integer. Default is 1, exact same design for each subject. Larger values will replicate designs, so more trials per subject.
  stat <- match.arg(stat)
  emc <- object
  dots <- list(...)
  data <- get_data(emc)
  design <- get_design(emc)
  return_trialwise_parameters <- isTRUE(dots[["return_trialwise_parameters"]])
  if (is.null(data[["subjects"]])) {
    jointModel <- TRUE
    all_samples <- emc
  } else {
    jointModel <- FALSE
    data <- list(data)
  }
  post_out <- vector("list", length = length(data))
  for (j in seq_along(data)) {
    if (jointModel) {
      emc <- single_out_joint(all_samples, j)
    }
    subjects <- levels(data[[j]][["subjects"]])
    if (grepl("MRI", design[[j]][["model"]]()[["type"]])) {
      design[[j]] <- add_design_fMRI_predict(design[[j]], emc)
      if (is.integer(design[[j]][["fMRI_design"]][[1]])) {
        design[[j]] <- design[[design[[j]][["fMRI_design"]][[1]]]]
      }
    }
    if (hyper) {
      mu <- do.call(
        get_pars,
        c(
          list(
            emc,
            selection = "mu",
            map = FALSE,
            return_mcmc = FALSE,
            merge_chains = TRUE,
            length.out = ceiling(n_post / length(emc))
          ),
          fix_dots(dots, get_pars)
        )
      )
      Sigma <- do.call(
        get_pars,
        c(
          list(
            emc,
            selection = "Sigma",
            map = FALSE,
            return_mcmc = FALSE,
            merge_chains = TRUE,
            remove_dup = FALSE,
            remove_constants = FALSE,
            length.out = ceiling(n_post / length(emc))
          ),
          fix_dots(dots, get_pars) # TODO
        )
      )
      pars <- get_alphas(mu, Sigma, subjects)
      pars <- pars[ , , 1:n_post] # With non-equally divisible n_post you get some remainder
      pars <- lapply(
        seq_len(dim(pars)[3]),
        function(i) {return(t(pars[ , , i]))}
      )
    } else {
      dots[["selection"]] <- "alpha"
      dots[["merge_chains"]] <- TRUE
      dots[["by_subject"]] <- TRUE
      samps <- do.call(
        get_pars,
        c(list(emc), fix_dots(dots, get_pars))
      )
      if (stat != "random") {
        p <- do.call(
          rbind,
          lapply(samps, function(x) {return(apply(x[[1]], 2, stat))})
        )
      }
      pars <- vector(mode = "list", length = n_post)
      for (i in seq_len(n_post)) {
        if (stat != "random") {
          pars[[i]] <- p
        } else {
          pars[[i]] <- do.call(
            rbind,
            lapply(samps, function(x) {return(x[[1]][sample(1:nrow(x[[1]]), 1), ])})
          )
        }
      }
    }
    simDat <- suppressWarnings(
      mclapply(
        seq_len(n_post),
        function(i) {
          do.call(
            make_data,
            c(
              list(pars[[i]], design = design[[j]], data = data[[j]]),
              fix_dots(dots, make_data)
            )
          )
        },
        mc.cores = n_cores
      )
    )
    in_bounds <- !sapply(simDat, is.logical)
    if (all(!in_bounds)) {
      stop("All samples fall outside of model bounds")
    }
    if (any(!in_bounds)) {
      good_post <- sample(seq_len(n_post), sum(!in_bounds))
      simDat[!in_bounds] <- suppressWarnings(
        mclapply(
          good_post,
          function(i) {
            do.call(
              make_data,
              c(
                list(pars[[i]], design = design[[j]], data = data[[j]], check_bounds = TRUE),
                fix_dots(dots, make_data)
              )
            )
          },
          mc.cores = n_cores
        )
      )
    }
    out <- cbind(
      postn = rep(
        seq_len(n_post),
        times = unlist(lapply(simDat, function(x) {return(dim(x)[1])}))
      ),
      do.call(rbind, simDat)
    )
    if (n_post == 1) {
      pars <- pars[[1]]
    }
    attr(out, "pars") <- pars
    if (return_trialwise_parameters) {
      attr(out, "trialwise_parameters") <- lapply(
        simDat,
        function(x) {return(attr(x, "trialwise_parameters"))}
      )
    }
    post_out[[j]] <- out
  }
  if(!jointModel) {
    post_out <- post_out[[1]]
  } else {
    joint_names <- get_joint_names(all_samples)
    names(post_out) <- joint_names
  }
  return(post_out)
}


# custom s3 funcs ---------------------------------------------------------
#' @rdname check
#' @export
check.emc <- function(
    emc,
    selection = c("mu", "sigma2", "alpha"),
    digits = 3,
    plot_worst = TRUE,
    version = c("old", "new"),
    ...
) {
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  version <- match.arg(version)
  dots <- list(...)
  out_list <- list()
  cat("Iterations:\n")
  print(chain_n(emc))
  if (emc[[1]][["type"]] == "single") {
    selection <- "alpha"
  }
  if (plot_worst) {
    mfrow <- coda_setmfrow(
      Nchains = length(emc),
      Nparms = length(selection),
      nplots = 1
    )
    par(mfrow = mfrow)
  }
  for (select in selection) {
    dots[["flatten"]] <- ifelse(select == "alpha", FALSE, TRUE)
    dots[["by_subject"]] <- TRUE
    ESS <- do.call(
      ess_summary,
      c(
        list(emc, selection = select, stat = NULL, version = version),
        fix_dots(dots, ess_summary)
      )
    )
    gds <- do.call(
      gd_summary,
      c(
        list(emc, selection = select, stat = NULL, version = version),
        fix_dots(dots, gd_summary)
      )
    )
    out <- list()
    max_gd <- -Inf
    for (name in names(ESS)) {
      combined <- rbind(
        round(gds[[name]], digits),
        round(ESS[[name]])
      )
      rownames(combined) <- c("Rhat", "ESS")
      out[[name]] <- combined
      if (max(gds[[name]]) > max_gd) {
        max_gd <- max(gds[[name]])
        cur_max <- name
        max_par <- names(gds[[name]])[which.max(gds[[name]])]
      }
    }
    if (length(ESS) > 1) {
      cat("\n", paste0(select, " highest Rhat : ", cur_max), "\n")
    } else {
      cat("\n", cur_max, "\n")
    }
    print(out[[cur_max]])
    if (plot_worst) {
      cur_dots <- dots
      if (select == "alpha") {
        cur_dots[["subject"]] <- cur_max
        cur_dots[["use_par"]] <- max_par
        cur_dots[["by_subject"]] <- TRUE
        MCMCs <- do.call(
          get_pars,
          c(
            list(emc, selection = select),
            fix_dots(cur_dots, get_pars)
          )
        )
        names(MCMCs) <- paste0("alpha : ", names(MCMCs))
      } else {
        cur_dots[["use_par"]] <- max_par
        MCMCs <- do.call(
          get_pars,
          c(
            list(emc, selection = select),
            fix_dots(cur_dots, get_pars)
          )
        )
      }
      cur_dots <- add_defaults(
        cur_dots,
        xlab = names(MCMCs)[1],
        ylab = "Highest Rhat parameter"
      )
      do.call(
        plot,
        c(
          list(MCMCs[[1]], auto.layout = FALSE, density = FALSE, ask = FALSE, smooth = FALSE),
          fix_dots_plot(cur_dots)
        )
      )
      legend("topleft", legend = paste0("Rhat : ", round(max_gd, digits)), bty = "n")
      legend("topright", legend = paste0("ESS : ", round(ESS[[cur_max]][max_par])), bty = "n")
    }
    out_list[[select]] <- out
  }
  return(invisible(out_list))
}

#' Convergence Checks for an emc Object
#'
#' Runs a series of MCMC diagnostic checks, prints those statistics to the
#' console, and makes traceplots of the worst converged parameter per selection.
#'
#' The potential scale reduction statistic \eqn{\hat{R}} is computed by [r_hat()],
#' and the effective sample size (ESS) is computed by [n_eff()]. In both cases,
#' by default implementations from the `coda` package are used
#' ([coda::gelman.diag()] and [coda::effectiveSize()], respectively), but more
#' up-to-date implementations from the `posterior` package are also supported,
#' by specifying `version = "new"` (see below).
#'
#' Full range of possible samples manipulations are described in `get_pars`.
#'
#' @param emc An emc object
#' @param selection A Character vector. Indicates which parameter types to check (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param digits Integer. How many digits to round the ESS and Rhat to in the plots
#' @param plot_worst Boolean. If `TRUE` also plots the chain plots for the worst parameter
#' @param version Character string, either `"old"` or `"new"`. `"old"` (default)
#'    calls [coda::gelman.diag()]) and [coda::effectiveSize()], while `"new"`
#'    uses more up-to-date implementations from the `posterior` package. See
#'    [r_hat()] and [n_eff()] for details.
#' @param ... Optional arguments that can be passed to `get_pars` or `plot.default` (see `par()`)
#'
#' @details The potential scale reduction statistic, \eqn{\hat{R}}, indicates
#' whether the chains are stationary and have converged to a common distribution
#' for a given parameter. If chains have not converged, \eqn{\hat{R}} will be
#' greater than 1; values above 1.1 are generally considered problematic.
#'
#' The effective sample size, ESS, indicates the number of independent draws from
#' a parameter's posterior distribution, which can be interpreted as the efficiency
#' of sampling.
#'
#' @return a list with the statistics for the worst converged parameter per selection
#' @examples
#' check(samples_LNR)
#' @export
check <- function(emc, ...){
  UseMethod("check")
}


#' @rdname parameters
#' @examples
#' # For posterior inference:
#' # Get 100 samples of the group-level mean (the default)
#' parameters(samples_LNR, N = 100)
#' # or from the individual-level parameters and mapped
#' parameters(samples_LNR, selection = "alpha", map = TRUE)
#' @export
parameters.emc <- function(x,selection = "mu", N = NULL, resample = FALSE, ...)
  # extracts and stacks chains into a matrix
{
  emc <- x

  if(is.null(N) || resample){
    nstage <- colSums(chain_n(emc))
    has_ran <- nstage[nstage != 0]
    if(!is.null(N)){
      N_res <- N
    } else{
      N_res <- nstage[names(has_ran)[length(has_ran)]]
    }
    N <- nstage[names(has_ran)[length(has_ran)]]
  }
  dots <- add_defaults(list(...), flatten = TRUE, length.out = N/length(emc))
  dots$merge_chains <- TRUE ; dots$return_mcmc <- FALSE
  out <- do.call(get_pars, c(list(emc,selection=selection), fix_dots(dots, get_pars)))
  if(selection == "alpha"){
    out <- aperm(out, c(3,1,2))
    if(resample){
      out <- apply(out, 2:3, sample, N_res, replace = T)
    }
    out <- apply(out, 3, identity, simplify=FALSE)
    for(i in 1:length(out)){
      out[[i]] <- as.data.frame(out[[i]])
      out[[i]] <- cbind(names(out)[i], out[[i]])
      colnames(out[[i]])[1] <- "subjects"
      if(!resample) out[[i]] <- out[[i]][1:N,]
    }
    out <- do.call(rbind, out)
    rownames(out) <- NULL
  } else{
    out <- as.data.frame(t(out))
    if(resample){
      out <- apply(out, 2, sample, N_res, replace = T)
    } else{
      out <- out[1:N,]
    }
  }
  out
}


#' Return Data Frame of Parameters
#'
#' @param x An emc or emc.prior object
#' @param selection String designating parameter type (e.g. mu, sigma2, correlation, alpha)
#' @param N Integer. How many samples to take from the posterior/prior. If `NULL` will return the full posterior
#' @param resample Boolean. If `TRUE` will sample N samples from the posterior with replacement
#' @param covariates For priors, possible covariates in the design
#' @param ... Optional arguments that can be passed to `get_pars`
#' @return A data frame with one row for each sample
#' (with a subjects column if selection = "alpha" and using draws from the posterior)
#' @export
parameters <- function(x, ...){
  UseMethod("parameters")
}


#' @rdname fit
#' @export

fit.emc <- function(emc, stage = NULL, iter = 1000, stop_criteria = NULL,
                    search_width = 1, step_size = 100, verbose = TRUE, fileName = NULL,
                    particle_factor=50, cores_per_chain = 1,
                    cores_for_chains = length(emc), max_tries = 20,
                    thin = FALSE,
                    ...){

  dots <- add_defaults(list(...), n_blocks = 1, verboseProgress = FALSE,
                       trim = TRUE, r_cores = 1)
  start_time <- Sys.time()
  stages_names <- c("preburn", "burn", "adapt", "sample")
  if(!is.null(stop_criteria) & !any(names(stop_criteria) %in% stages_names)){
    stop_criteria[["sample"]] <- stop_criteria
  }
  if(is.null(stop_criteria)){
    stop_criteria <- vector("list", length = 4)
    names(stop_criteria) <- stages_names
  }
  for(stage_name in stages_names){
    if(!stage_name %in% names(stop_criteria)) stop_criteria[stage_name] <- list(NULL)
  }

  stop_criteria <- stop_criteria[stages_names]
  stop_criteria <- mapply(get_stop_criteria, stages_names, stop_criteria, MoreArgs = list(type = emc[[1]]$type))
  if(is.null(stop_criteria[["sample"]]$iter)) stop_criteria[["sample"]]$iter <- iter
  names(stop_criteria) <- stages_names
  if (is.character(emc)) {
    emc <- fix_fileName(emc)
    if(is.null(fileName)) fileName <- emc
    emc <- loadRData(emc)
  }
  if(is.null(stage)){
    last_stage <- get_last_stage(emc)
  } else{
    last_stage <- stage
  }
  stages <- c('preburn', 'burn', 'adapt', 'sample')
  to_run <- stages[which(stages == last_stage):length(stages)]
  for(stage in to_run){
    emc <- run_emc(emc, stage = stage, stop_criteria[[stage]], cores_for_chains = cores_for_chains, search_width = search_width,
                   step_size = step_size,  verbose = verbose, verboseProgress = dots$verboseProgress,
                   fileName = fileName, particle_factor =  particle_factor, trim = dots$trim,
                   cores_per_chain = cores_per_chain, max_tries = max_tries, thin = thin, n_blocks = dots$n_blocks,
                   r_cores = dots$r_cores)
  }
  if (verbose) print(Sys.time()-start_time)
  return(emc)
}

#' Model Estimation in EMC2
#'
#' General purpose function to estimate models specified in EMC2.
#'
#' @details
#'
#' ``stop_criteria`` is either a list of lists with names of the stages,
#' or a single list in which case its assumed to be for the sample `stage` (see examples).
#' The potential stop criteria to be set are:
#'
#' ``selection`` (character vector): For which parameters the ``stop_criteria`` should hold
#'
#' ``mean_gd`` (numeric): The mean Gelman-Rubin diagnostic across all parameters in the selection
#'
#' ``max_gd`` (numeric): The max Gelman-Rubin diagnostic across all parameters in the selection
#'
#' ``min_unique`` (integer): The minimum number of unique samples in the MCMC chains across all parameters in the selection
#'
#' ``min_es`` (integer): The minimum number of effective samples across all parameters in the selection
#'
#' ``omit_mpsrf`` (Boolean): Whether to include the multivariate point-scale reduction factor in the Gelman-Rubin diagnostic. Default is ``FALSE``.
#'
#' ``iter`` (integer): The number of MCMC samples to collect.
#'
#' The estimation is performed using particle-metropolis within-Gibbs sampling.
#' For sampling details see:
#'
#' Gunawan, D., Hawkins, G. E., Tran, M.-N., Kohn, R., & Brown, S. (2020).
#' New estimation approaches for the hierarchical linear ballistic accumulator model.
#' *Journal of Mathematical Psychology* ,96, 102368. doi.org/10.1016/j.jmp.2020.102368
#'
#' Stevenson, N., Donzallaz, M. C., Innes, R. J., Forstmann, B., Matzke, D., & Heathcote, A. (2024).
#' EMC2: An R Package for cognitive models of choice. doi.org/10.31234/osf.io/2e4dq
#'
#'
#' @param emc An emc object created with ``make_emc``,
#' or a path to where the emc object is stored.
#' @param stage A string. Indicates which stage to start the run from, either ``preburn``, ``burn``, ``adapt`` or ``sample``.
#' If unspecified, it will run the subsequent stage (if there is one).
#' @param iter An integer. Indicates how many iterations to run in the sampling stage.
#' @param search_width A double. Tunes target acceptance probability of the MCMC process.
#' This fine-tunes the width of the search space to obtain the desired acceptance probability.
#' 1 is the default width, increases lead to broader search.
#' @param step_size An integer. After each step, the stopping requirements as specified
#' by ``stop_criteria`` are checked and proposal distributions are updated. Defaults to 100.
#' @param verbose Logical. Whether to print messages between each step with the current status regarding the ``stop_criteria``.
#' @param fileName A string. If specified, will auto-save emc object at this location on every iteration.
#' @param particle_factor An integer. ``particle_factor`` multiplied by the square
#' root of the number of sampled parameters determines the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across
#' participant calculations. Only available on Linux or Mac OS. For Windows, only
#' parallelization across chains (``cores_for_chains``) is available.
#' @param cores_for_chains An integer. How many cores to use across chains.
#' Defaults to the number of chains. The total number of cores used is equal to ``cores_per_chain`` * ``cores_for_chains``.
#' @param max_tries An integer. How many times should it try to meet the finish
#' conditions as specified by ``stop_criteria``? Defaults to 20. ``max_tries`` is
#' ignored if the required number of iterations has not been reached yet.
#' @param stop_criteria A list. Defines the stopping criteria and for which types
#' of parameters these should hold. See the details and examples section.
#' @param thin A boolean. If `TRUE` will automatically thin the MCMC samples, closely matched to the ESS.
#' Can also be set to a double, in which case 1/thin of the chain will be removed (does not have to be an integer).
#' @param ... Additional optional arguments
#' @return An emc object
#' @examples \donttest{
#' # Define a design first
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' # We also define a match function for lM
#' matchfun=function(d)d$S==d$lR
#'
#' # Drop most subjects for speed
#' dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
#' dat$subjects <- droplevels(dat$subjects)
#'
#' design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
#'                      formula=list(m~lM,s~1,t0~1),
#'                      contrasts=list(m=list(lM=ADmat)))
#' # Before fit can be called, we first need to make an emc object
#' LNR_s <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2)
#' # Run fit, here illustrating how to use stop_criteria (also for speed purposes)
#' LNR_s <- fit(LNR_s, cores_for_chains = 1, stop_criteria = list(
#'   preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
#'   sample = list(iter = 25, max_gd = 2)), verbose = FALSE, particle_factor = 30, step_size = 25)
#'}
#' @export
fit <- function(emc, ...){
  UseMethod("fit")
}

#' @rdname recovery
#' @export

recovery.emc <- function(emc, true_pars,
                          selection = "mu",
                          layout=NA,
                          do_CI = TRUE,
                          correlation = "pearson",
                          stat = "rmse",
                          digits = 3,
                          CI = .95, ci_plot_args = list(), ...)
{
  dots <- list(...)
  type <- emc[[1]]$type
  if(length(dots$subject) == 1 || emc[[1]]$n_subjects == 1) dots$by_subject <- TRUE
  dots$merge_chains <- TRUE
  MCMC_samples <- do.call(get_pars, c(list(emc, selection = selection), fix_dots(dots, get_pars)))
  true_MCMC_samples <- NULL
  if(!is(true_pars, "emc")){
    true_pars <- do.call(get_pars, c(list(emc, selection = selection, type = type, true_pars = true_pars),
                                        fix_dots(dots, get_pars, exclude = c("thin", "filter"))))
  } else{
    true_MCMC_samples <- do.call(get_pars, c(list(true_pars, selection = selection), fix_dots(dots, get_pars)))
    true_pars <- NULL
  }
  # pearson <- spearman <- rmse <- coverage <- setNames(numeric(length(MCMC_samples)),names(MCMC_samples))
  stats_list <- list()
  if(any(is.na(layout))){
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(MCMC_samples),
                                 nplots = 1))
  } else{par(mfrow=layout)}
  for(i in 1:length(MCMC_samples)){

    cur_name <- names(MCMC_samples)[i]
    stats <- get_recovery_stats(MCMC_samples[[i]], true_MCMC_samples[[i]],
                                true_pars[[i]], CI)
    if(do_CI){
      ylim <- range(c(stats$true, stats$recovered))
    } else{
      ylim <- range(c(stats$true, stats$recovered[,"50%"]))
    }
    main_name <- ifelse(length(MCMC_samples) == 1, cur_name, paste0(selection, ": ", cur_name))
    cur_dots <- add_defaults(dots, main = main_name, ylim = ylim, xlim = ylim, xlab = "Generated", ylab = "Estimated")
    do.call(plot, c(list(stats$true[,"50%"],stats$recovered[,"50%"]), fix_dots_plot(cur_dots)))
    abline(a=0,b=1,lty=3)
    if(do_CI){
      cur_ci_plot_args <- add_defaults(ci_plot_args, col = "grey", angle = 90, length = .05)
      do.call(arrows, c(list(stats$true[,"50%"],stats$recovered[,2],stats$true[,"50%"],stats$recovered[,3]),
                        fix_dots_plot(cur_ci_plot_args)))
      do.call(arrows, c(list(stats$true[,"50%"],stats$recovered[,2],stats$true[,"50%"],stats$recovered[,1]),
                        fix_dots_plot(cur_ci_plot_args)))
      if(is.null(true_pars)){
        do.call(arrows, c(list(stats$true[,"50%"],stats$recovered[,2],stats$true[,3],stats$recovered[,2]),
                          fix_dots_plot(cur_ci_plot_args)))
        do.call(arrows, c(list(stats$true[,"50%"],stats$recovered[,2],stats$true[,1],stats$recovered[,2]),
                          fix_dots_plot(cur_ci_plot_args)))
      }
    }
    if(correlation == "pearson"){
      legend("topleft",paste("r(pearson) = ",round(stats$pearson,digits)),bty="n")
    } else if(correlation == "spearman"){
      legend("topleft",paste("r(spearman) = ",round(stats$spearman,digits)),bty="n")
    }
    if (stat == "rmse") {
      legend("bottomright",paste("RMSE = ",round(stats$rmse,digits)),bty="n")
    } else if(stat == "coverage") {
      legend("bottomright",paste("95% Coverage = ",round(stats$coverage,digits)),bty="n")
    }
    stats <- make_recov_summary(stats)
    stats_list[[cur_name]] <- stats
  }
  return(invisible(stats_list))
}


#' Recovery Plots
#'
#' Plots recovery of data generating parameters/samples.
#' Full range of samples manipulations described in `get_pars`
#'
#' @param emc An emc object
#' @param true_pars A vector of data-generating parameters or an emc object with data-generating samples
#' @param selection A Character vector. Indicates which parameter types to plot (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param do_CI Boolean. If `TRUE` will also include bars representing the credible intervals
#' @param correlation Character. Which correlation to include in the plot. Options are either `pearson` or `spearman`
#' @param stat Character. Which statistic to include in the plot. Options are either `rmse` or `coverage`
#' @param digits Integer. How many digits to round the statistic and correlation in the plot to
#' @param CI Numeric. The size of the credible intervals. Default is .95 (95%).
#' @param ci_plot_args A list. Optional additional arguments to be passed to plot.default for the plotting of the credible intervals (see `par()`)
#' @param ... Optional arguments that can be passed to `get_pars` or `plot.default` (see `par()`)
#'
#' @return Invisible list with RMSE, coverage, and Pearson and Spearman correlations.
#' @export
#' @examples
#' # Make up some values that resemble posterior samples
#' # Normally this would be true values that were used to simulate the data
#' # Make up some values that resemble posterior samples
#' # Normally this would be true values that were used to simulate the data
#' pmat <- matrix(rnorm(12, mean = c(-1, -.6, -.4, -1.5), sd = .01), ncol = 4, byrow = TRUE)
#' # Conventionally this would be created before one makes data with true values
#' recovery(samples_LNR, pmat, correlation = "pearson", stat = "rmse", selection = "alpha")

#' # Similarly we can plot recovery of other parameters with a set of true samples
#' true_samples <- samples_LNR # Normally this would be data-generating samples
#' recovery(samples_LNR, true_samples, correlation = "pearson", stat = "rmse",
#'          selection = "correlation", cex = 1.5,
#'          ci_plot_args = list(lty = 3, length = .2, lwd = 2, col = "brown"))

recovery <- function(emc, ...){
  UseMethod("recovery")
}

#' @rdname hypothesis
#' @export
hypothesis.emc <- function(emc, parameter = NULL, H0 = 0, fun = NULL,selection = "mu",
                          do_plot = TRUE, use_prior_lim = TRUE,
                          N = 1e4, prior_args = list(), ...){
  dots <- add_defaults(list(...), flatten = TRUE)
  type <- emc[[1]]$type
  if (length(emc[[1]]$data)==1) selection <- "alpha"
  if(selection == "alpha" & type != "single") stop("For savage-dickey ratio, selection cannot be alpha for hierarchical models")
  prior <- get_prior(emc)


  psamples <-  get_objects(design = get_design(emc),
                           type = emc[[1]]$type, sample_prior = T, prior = prior,
                           selection = selection, N = N, sampler = emc)
  psamples <- do.call(get_pars, c(list(psamples, selection = selection, merge_chains = TRUE, return_mcmc = FALSE, by_subject = TRUE,
                                          type = emc[[1]]$type),
                                     fix_dots(dots, get_pars, exclude = c("thin", "filter"))))
  samples <- do.call(get_pars, c(list(emc, selection = selection, merge_chains = TRUE, return_mcmc = FALSE, by_subject = TRUE),
                                    fix_dots(dots, get_pars)))
  if(type == "single"){
    if(ncol(samples) > 1 & !is.null(dots$subject)){
      stop("with non-hierarichal run with multiple subjects, you must specify which subject")
    } else {
      samples <- samples[,1,]
      psamples <- psamples[,1,]
    }
  }
  if(is.null(fun)){
    idx <- rownames(samples) == parameter
    samples <- samples[idx,]
    idxp <- rownames(psamples) == parameter
    psamples <- psamples[idxp,]
  } else{
    samples <- apply(as.data.frame(samples), 2, fun)
    psamples <- apply(as.data.frame(psamples), 2, fun)
  }
  min_bound <- min(min(psamples), H0)
  max_bound <- max(max(psamples), H0)
  diff <- max_bound - min_bound
  pdensity <- density(psamples, from = min_bound - diff/2, to = max_bound + diff/2)
  pdfun <-approxfun(pdensity)

  min_bound <- min(min(samples), H0)
  max_bound <- max(max(samples), H0)
  diff <- max_bound - min_bound
  post_density <- density(samples, from = min_bound - diff/2, to = max_bound + diff/2)
  post_dfun <-approxfun(post_density)
  if(do_plot){
    if(is.null(dots$xlim)){
      if(use_prior_lim){
        dots$xlim <- range(c(quantile(samples, c(0.025, 0.975)),
                             quantile(psamples, c(0.025, 0.975)), H0 + .01, H0 - .01))
      } else{
        dots$xlim <- range(c(quantile(samples, c(0.025, 0.975)), H0 + .01, H0 - .01))
      }
    }
    prior_args <- add_defaults(prior_args, cex = 2, col = "red", lwd = 1.5)
    dots <- add_defaults(dots, cex = 2, col = "black", lwd = 1.5, main = "Prior and posterior density")
    do.call(plot, c(list(post_density), fix_dots_plot(dots)))
    do.call(lines, c(list(pdensity), fix_dots_plot(prior_args)))
    do.call(points, c(list(H0 , post_dfun(H0)), fix_dots_plot(dots)))
    do.call(points, c(list(H0, pdfun(H0)), fix_dots_plot(prior_args)))
    legend("topright", legend = c("posterior", "prior"), pch = c(1, 1), col = c(dots$col, prior_args$col))
  }
  return(pdfun(H0)/post_dfun(H0))
}

#' Within-Model Hypothesis Testing
#'
#' Approximates the Bayes factor for parameter effects using the savage-dickey ratio.
#'
#' Note this is different to the computation of the marginal deviance in `compare`
#' since it only considers the group level effect and not the whole model (i.e. subject-level parameters).
#' For details see: Wagenmakers, Lodewyckx, Kuriyal, & Grasman (2010).
#'
#' @param parameter A string. A parameter which you want to compare to H0. Will not be used if a FUN is specified.
#' @param H0 An integer. The H0 value which you want to compare to
#' @param fun A function. Specifies an operation to be performed on the sampled or mapped parameters.
#' @param do_plot Boolean. If `FALSE` will omit the prior-posterior plot and only return the savage-dickey ratio.
#' @inheritParams plot_pars
#' @return The Bayes factor for the hypothesis against H0.
#' @examples
#' # Here the emc object has an effect parameter (e.g. m),
#' # that maps onto a certain hypothesis.
#' # The hypothesis here is that m is different from zero.
#' # We can test whether there's a group-level effect on m:
#' hypothesis(samples_LNR, parameter = "m")
#' # Alternatively we can also test whether two parameters differ from each other
#' mdiff <- function(p)diff(p[c("m","m_lMd")])
#' hypothesis(samples_LNR,fun=mdiff)
#' @export
hypothesis <- function(emc, ...){
  UseMethod("hypothesis")
}




#' @rdname credible
#' @export
credible.emc <- function(x,x_name=NULL,x_fun=NULL,x_fun_name="fun", selection = "mu",
                   y=NULL,y_name=NULL,y_fun=NULL,y_fun_name="fun",
                   x_subject=NULL,y_subject=NULL,
                   mu=0,alternative = c("less", "greater")[1],
                   probs = c(0.025,.5,.975),digits=2,p_digits=3,print_table=TRUE,
                   ...)

{

  dots <- add_defaults(list(...), flatten = TRUE)
  get_effect <- function(x,p_name=NULL,fun=NULL)
  {
    x <- do.call(rbind,x)
    if (!is.null(fun)) return(apply(x,1,fun))
    x[,p_name]
  }
  if (is.null(x_name) & is.null(x_fun))
    stop("x_name or x_fun must be supplied")
  if (is.null(y_fun) && is.null(y_name)) y_name <- x_name
  if (is.null(y_fun) && !is.null(x_fun)) y_fun <- x_fun

  # Process x
  if (!is(x[[1]], "pmwgs")) stop("x must be a list of pmwgs objects")
  if (length(x[[1]]$data)==1) selection <- "alpha"
  if(is.numeric(x_subject)) x_subject <- names(x[[1]]$data)[x_subject]
  x <- do.call(get_pars, c(list(x,selection=selection, merge_chains = TRUE, by_subject = TRUE),
                              fix_dots(add_defaults(dots, subject = x_subject), get_pars)))[[1]]


  # Individual subject analysis
  # Check test is valid
  if (is.null(x_fun) && !all(x_name %in% dimnames(x[[1]])[[2]]) )
    stop("x_name not present in samples")
  if (length(x_name)>2)
    stop("x_name must be length 1 or 2")
  # Get x effect
  x <- get_effect(x,x_name,x_fun)
  if (length(x_name)>1) {
    x <- -apply(x,1,diff)
    x_name <- paste(x_name,collapse="-")
  }
  if (is.null(x_name)) x_name <- x_fun_name
  if (is.null(y)) {
    p <- mean(x<mu)
    if (alternative=="greater") p <- 1-p
    tab <- cbind(quantile(x,probs),c(NA,mu,NA))
    attr(tab,alternative) <- p
    dimnames(tab)[[2]] <- c(x_name,"mu")
  } else {
    if(is.numeric(y_subject)) y_subject <- names(y[[1]]$data)[y_subject]

    if (!is(y[[1]], "pmwgs")) stop("y must be a list of pmwgs objects")
    y <- do.call(get_pars, c(list(y,selection=selection, merge_chains = TRUE, by_subject = TRUE),
                                fix_dots(add_defaults(dots, subject = y_subject), get_pars)))[[1]]
    if (is.null(y_fun) && !all(y_name %in% dimnames(y[[1]])[[2]]) )
      stop("y_name not present in samples")
    if (length(y_name)>2)
      stop("y_name must be length 1 or 2")
    y <- get_effect(y,y_name,y_fun)
    if (length(y_name)>1) {
      y <- -apply(y,1,diff)
      y_name <- paste(y_name,collapse="-")
    }
    if (length(x)>length(y)) x <- x[1:length(y)] else y <- y[1:length(x)]
    d <- x-y
    p <- mean(d<0)
    if (alternative=="greater") p <- 1-p
    tab <- cbind(quantile(x,probs),quantile(y,probs),quantile(d,probs))
    attr(tab,alternative) <- p
    if (is.null(y_name)) y_name <- y_fun_name
    if (x_name==y_name){
      colnames(tab) <- c(paste(x_name,c(x_subject,y_subject),sep="_"),
                              paste(x_subject,y_subject,sep="-"))
    } else{
      colnames(tab) <- c(x_name,y_name,paste(x_name,y_name,sep="-"))
    }

  }
  if (print_table) {
    ptab <- tab
    ptab <- round(ptab,digits)
    attr(ptab,alternative) <- round(attr(ptab,alternative),p_digits)
    print(ptab)
  }
  invisible(tab)
}

#' Posterior Credible Interval Tests
#'
#' Modeled after `t.test`, returns the credible interval of the parameter or test
#' and what proportion of the posterior distribution (or the difference in posterior distributions
#' in case of a two sample test) overlaps with mu.
#' For a one sample test provide `x` and for two sample also provide `y`.
#' Note that for comparisons within one model, we recommend using `hypothesis()` if the priors
#' were well chosen.
#'
#' @param x An emc object
#' @param x_name A character string. Name of the parameter to be tested for `x`
#' @param x_fun Function applied to the MCMC chains to create
#' variable to be tested.
#' @param y A second emc object
#' @param y_name A character string. Name of the parameter to be tested for `y`
#' @param y_fun Function applied to the MCMC chains to create
#' variable to be tested.
#' @param x_subject Integer or name selecting a subject
#' @param y_subject Integer or name selecting a subject
#' @param mu Numeric. `NULL` value for single sample test if `y` is not supplied (default 0)
#' @param alternative `less` or `greater` determining direction of test probability
#' @param probs Vector defining quantiles to return.
#' @param digits Integer, significant digits for estimates in printed results
#' @param p_digits Integer, significant digits for probability in printed results
#' @param print_table Boolean (defaults to `TRUE`) for printing results table
#' @param selection A character string designating parameter type (e.g. `alpha` or `covariance`)
#' @param x_fun_name Name to give to quantity calculated by `x_fun`
#' @param y_fun_name Name to give to quantity calculated by `y_fun`
#' @param ... Additional optional arguments that can be passed to `get_pars`
#' @examples{
#' # Run a credible interval test (Bayesian ''t-test'')
#' credible(samples_LNR, x_name = "m")
#' # We can also compare between two sets of emc objects
#'
#' # # Now without a ~ E
#' # design_null <- design(data = forstmann,model=DDM,
#' #                            formula =list(v~0+S,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
#' #                            constants=c(s=log(1)))
#' #
#' # null_model <- make_emc(forstmann, design_null)
#' # null_model <- fit(null_model)
#' # credible(x = null_model, x_name = "a", y = full_model, y_name = "a")
#' #
#' # # Or provide custom functions:
#' # credible(x = full_model, x_fun = function(d) d["a_Eaccuracy"] - d["a_Eneutral"])
#' }
#' @return Invisible results table with no rounding.
#' @export
credible <- function(x, ...){
  UseMethod("credible")
}

#' Shorten an emc Object
#'
#' @inheritParams get_pars
#' @param ... additional optional arguments
#' @param x an emc object
#' @param keep_stages Boolean. If `TRUE`, will not remove samples from unselected stages.
#' @return A shortened emc object
#' @export
#'
#' @examples
#' subset(samples_LNR, length.out = 10)
subset.emc <- function(x, stage = "sample", filter = NULL, thin = 1, keep_stages = FALSE,
                       length.out = NULL, ...){
  x <- lapply(x, remove_samples, stage = stage, filter = filter,
                thin = thin, length.out = length.out, keep_stages = keep_stages)
  class(x) <- "emc"
  return(x)
}


#' @rdname gd_summary
#' @export
gd_summary.emc <- function(
    emc,
    selection = "mu",
    omit_mpsrf = TRUE,
    stat = "max",
    stat_only = FALSE,
    digits = 3,
    version = "old",
    ...
) {
  out <- get_summary_stat(
    emc, selection, r_hat,
    stat = stat, stat_only = stat_only, digits = digits,
    omit_mpsrf = omit_mpsrf, version = version, ...
  )
  return(out)
}

#' @rdname ess_summary
#' @export
ess_summary.emc <- function(
    emc,
    selection = "mu",
    stat = "min",
    stat_only = FALSE,
    digits = 1,
    version = "old",
    ...
) {
  out <- get_summary_stat(
    emc, selection, n_eff,
    stat = stat, stat_only = stat_only, digits = digits,
    version = version, ...
  )
  return(out)
}

#' @rdname credint
#' @export
credint.emc <- function(x, selection="mu", probs = c(0.025, .5, .975),
                                  digits = 3, ...){
  out <- get_summary_stat(x, selection, get_posterior_quantiles,
                          probs = probs, digits = digits, ...)
  return(out)
}

#' R-hat convergence diagnostic
#'
#' Computes the potential scale reduction factor (\eqn{\hat{R}}) for a selected
#' model parameter type in an `emc` object. Supports both the legacy Gelman–Rubin
#' diagnostic (Gelman & Rubin, 1992; Brooks & Gelman, 1998), and the more
#' up-to-date implementation proposed by Vehtari et al. (2021).
#'
#' @param emc An emc object.
#' @param selection A Character vector. Indicates which parameter types to check
#'    (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param omit_mpsrf Boolean. If `TRUE` (default) the multivariate point scale
#'    reduction factor (MPSRF) is *not* returned (see [coda::gelman.diag()] for
#'    details). Only relevant if `version = "old"` (see below).
#' @param stat A string. Should correspond to a function that can be applied to
#'    a vector, which will be performed on the vector/rows or columns of the
#'    matrix of the parameters
#' @param stat_only Boolean. If `TRUE` will only return the result of the
#'    applied `stat` function, otherwise returns both the stat result and the
#'    result of the function on all parameters.
#' @param digits Integer. How many digits to round the output to.
#' @param version Character string, either `"old"` or `"new"`. `"old"` (default)
#'    calls [coda::gelman.diag()] with split chains, while `"new"` uses a vendored
#'    implementation of [posterior::rhat()] based on Vehtari et al. (2021).
#'    See [r_hat()] for details.
#' @param ... Optional additional arguments that can be passed to `get_pars`
#'
#' @return A matrix or vector of R-hat values for the selected parameter type.
#'
#' @references
#' Gelman, A., & Rubin, D.B. (1992). Inference from iterative simulation using
#' multiple sequences. *Statistical Science*, *7*, 457–511.
#'
#' Brooks, S.P., & Gelman, A. (1998). General methods for monitoring convergence
#' of iterative simulations. *Journal of Computational and Graphical Statistics*,
#' *7*, 434–455.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: An improved \eqn{\hat{R}} for
#' assessing convergence of MCMC. *Bayesian Analysis*, *16*(2), 667–718.
#'
#' @examples
#' gd_summary(samples_LNR, selection = "correlation", stat = "mean", flatten = TRUE)
#'
#' @export
gd_summary <- function(emc, ...){
  UseMethod("gd_summary")
}

#' Effective Sample Size
#'
#' Computes the effective sample size (ESS) for a selected model parameter type
#' in an `emc` object. Supports both the legacy [coda::effectiveSize()] calculation
#' and a more modern implementation adapted from [posterior::ess_basic()] based
#' on Vehtari et al. (2021).
#'
#' @param emc An emc object.
#' @param selection A Character vector. Indicates which parameter types to check
#'    (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param stat A string. Should correspond to a function that can be applied to
#'    a vector, which will be performed on the vector/rows or columns of the
#'    matrix of the parameters
#' @param stat_only Boolean. If `TRUE` will only return the result of the
#'    applied `stat` function, otherwise returns both the stat result and the
#'    result of the function on all parameters.
#' @param digits Integer. How many digits to round the output to.
#' @param version Character string, either `"old"` or `"new"`. `"old"` (default)
#'    calls [coda::effectiveSize()], while `"new"` uses a vendored implementation
#'    of [posterior::ess_basic()] based on Vehtari et al. (2021).
#'    See [n_eff()] for details.
#' @param ... Optional additional arguments that can be passed to `get_pars`
#'
#' @return A matrix or vector of ESS values for the selected parameter type.
#'
#' @examples
#' ess_summary(samples_LNR, selection = "alpha")
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: An improved \eqn{\hat{R}} for
#' assessing convergence of MCMC. *Bayesian Analysis*, *16*(2), 667–718.
#'
#' @export
ess_summary <- function(emc, ...){
  UseMethod("ess_summary")
}

#' Posterior Quantiles
#'
#' Returns the quantiles of the selected parameter type.
#' Full range of possible samples manipulations described in `get_pars`.
#'
#' @inheritParams gd_summary.emc
#' @param x An emc or emc.prior object
#' @param probs A vector. Indicates which quantiles to return from the posterior.
#' @return A list of posterior quantiles for each parameter group in the selected parameter type.
#' @export
#'
#' @examples
#' credint(samples_LNR)
credint <- function(x, ...){
  UseMethod("credint")
}

#' @rdname get_data
#' @export
get_data.emc <- function(emc) {
  if(is.null(emc[[1]]$data)) return(NULL) # Prior samples
  if(is.null(emc[[1]]$data[[1]]$subjects)){ # Joint model
    dat <- vector("list", length(emc[[1]]$data[[1]]))
    for(i in 1:length(dat)){
      design <- get_design(emc)[[i]]
      tmp <- do.call(rbind,lapply(emc[[1]]$data,function(x){
        cur <- x[[i]]
        if(!is.null(cur$winner) && (length(unique(cur$lR)) > 1)){
          cur <- cur[cur$winner,]
        }
        expand <- attr(cur,"expand")
        if(is.null(expand)) expand <- 1:nrow(cur)
        return(cur[expand,])
      }))
      row.names(tmp) <- NULL
      tmp <- tmp[,!(colnames(tmp) %in% c("trials","lR","lM", "winner", "SlR", "RACE", names(design$Ffunctions)))]
      dat[[i]] <- tmp
    }
    names(dat) <- get_joint_names(emc)
  } else{
    design <- get_design(emc)[[1]]
    dat <- do.call(rbind,lapply(emc[[1]]$data,function(x){
      if(!is.null(x$winner) && (length(unique(x$lR)) > 1)){
        # Only expand winner for race models
        x <- x[x$winner,]
      }

      expand <- attr(x,"expand")
      if(is.null(expand)) expand <- 1:nrow(x)
      return(x[expand,])
    }))
    row.names(dat) <- NULL
    dat <- dat[,!(colnames(dat) %in% c("trials","lR","lM","winner", "SlR", "RACE", names(design$Ffunctions)))]
  }
  return(dat)
}
#' Get Data
#'
#' Extracts data from an emc object
#'
#' @details
#' emc adds columns and rows to a dataframe in order to facilitate efficient likelihood calculations.
#' This function will return the data as provided originally.
#' @param emc an emc object
#' @return A dataframe of the original data
#' @export
#' @examples
#' get_data(samples_LNR)
get_data <- function(emc){
  UseMethod("get_data")
}

#' @rdname get_prior
#' @export
get_prior.emc <- function(emc){
  prior <- emc[[1]]$prior
  attr(prior, "type") <- emc[[1]]$type
  class(prior) <- "emc.prior"
  return(prior)
}

#' Get Prior
#'
#' Extracts prior from an emc object
#'
#' @param emc an emc object
#' @return A prior with class emc.prior
#' @export
#' @examples
#' get_prior(samples_LNR)
get_prior <- function(emc){
  UseMethod("get_prior")
}

#' @rdname get_design
#' @export
get_design.emc <- function(x){
  # For backwards compatibility
  if(!is.null(attr(x, "design_list"))){
    emc_design <- attr(x, "design_list")
  } else{
    emc_design <- get_design(get_prior(x))
  }
  class(emc_design) <- "emc.design"
  return(emc_design)
}

#' @rdname get_group_design
#' @export
get_group_design.emc <- function(x){
  group_design <- get_group_design(get_prior(x))
  if(!is.null(group_design))  class(group_design) <- "emc.group_design"
  return(group_design)
}



#' Plot Design
#'
#' Makes design illustration by plotting simulated data based on the design
#'
#' @param x An `emc` or `emc.prior` object containing the design to plot
#' @param data Optional data to overlay on the design plot
#' @param factors Factors to use for varying parameters
#' @param plot_factor Optional. Make separate plots for each level of this factor
#' @param n_data_sim If data is provided, number of simulated datasets to generate for the plot. Default is 10.
#' @param p_vector Only needed when x is an `emc.design` object, which parameters to use for data generation.
#' @param functions A named list of functions that create additional columns in the data.
#' @param ... Additional arguments to pass to `make_design_plot`
#' @return No return value. Just plots the design
#' @export
plot_design <- function(x, data = NULL, factors = NULL, plot_factor = NULL, n_data_sim = 10, p_vector = NULL,
                        functions = NULL, ...){
  UseMethod("plot_design")
}

#' @rdname plot_design
#' @export
plot_design.emc <- function(x, data = NULL, factors = NULL, plot_factor = NULL, n_data_sim = 10, p_vector = NULL,
                            functions = NULL, ...){
  p_vector <- credint(x, probs = .5)[[1]]
  design <- get_design(x)[[1]]
  plot(design, p_vector, data = data, factors = factors, plot_factor = plot_factor, n_data_sim = n_data_sim,
       p_vector = p_vector, functions = functions, ...)
}

#' @rdname mapped_pars
#' @export
mapped_pars.emc <- function(x, p_vector = NULL, model = NULL, digits=3,remove_subjects=TRUE,
                                  covariates=NULL,...){
  if(is.null(p_vector)) p_vector <- credint(x, probs = .5)[[1]]
  design <- get_design(x)
  mapped_pars(design, p_vector, digits = digits, remove_subjects=remove_subjects,
              covariates=covariates,...)
}

#' Get Design
#'
#' Extracts design from an emc object
#'
#' @param x an `emc` or `emc.prior` object
#' @return A design with class emc.design
#' @export
#' @examples
#' get_design(samples_LNR)
get_design <- function(x){
  UseMethod("get_design")
}

#' Get Group Design
#'
#' Extracts group design from an emc object
#'
#' @param x an `emc` or `emc.prior` object
#' @return A design with class emc.group_design
#' @export
get_group_design <- function(x){
  UseMethod("get_group_design")
}



#' @rdname sampled_pars
#' @export
sampled_pars.emc <- function(x,group_design=NULL,doMap=FALSE, add_da = FALSE, all_cells_dm = FALSE, data = NULL){
  return(sampled_pars(get_design(x), group_design = group_design, doMap = doMap,
                          add_da = add_da, all_cells_dm = all_cells_dm, data = data))
}

#' @rdname auto_thin
#' @export
auto_thin.emc <- function(emc, stage = "sample", selection = c("alpha", "mu"), ...){
  ess <- 0
  for(select in selection){
    ess <- ess + ess_summary(emc, selection = select, stage = stage, stat_only = TRUE, stat = "mean")
  }
  ess <- ess/length(selection)
  return(subset(emc, stage = stage, length.out = min(ess/length(emc), sum(chain_n(emc)[1,stage]))))
}


#' Automatically Thin an emc Object
#'
#' Uses the effective sample size of `selection` to determine how much to optimally thin an emc object
#'
#' @inheritParams get_pars
#' @param selection Which parameter types (i.e. 'alpha' or 'mu' to consider when determining the effective sample size)
#' @param ... additional optional arguments
#'
#' @export
auto_thin <- function(emc, stage = "sample", selection = c("alpha", "mu"), ...){
  UseMethod("auto_thin")
}

