#' @export
print.emc <- function(x, ...){
  n_chain <- chain_n(x)
  cat("Iterations: \n")
  print(chain_n(x))
  sub_names <- names(x[[1]]$data)
  cat("\n")
  cat("Subjects: \n")
  print(sub_names)
  par_names <- x[[1]]$par_names
  cat("\n")
  cat("Parameters: \n")
  print(par_names)
  return(invisible(list(iterations = n_chain, subjects = sub_names, parameters = par_names)))
}

#' Summary statistics for emc objects
#'
#' Computes quantiles, `Rhat` and `ESS` for selected model parameters.
#'
#' Note that if `selection = alpha` and `by_subject = TRUE` (default)
#' is used, summary statistics are computed at the individual level.
#" Only the first subject's summary output is printed
#' to the console but summary statistics for all subjects are returned by the
#' function.
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
summary.emc <- function(object, selection = c("mu", "sigma2", "alpha"), probs = c(0.025, .5, 0.975),
                        digits = 3, ...){
  dots <- list(...)
  dots <- add_defaults(dots, by_subject = TRUE)
  if(attr(object[[1]], "variant_funs")$type == "single"){
    selection <- "alpha"
  }
  out_list <- list()
  for(select in selection){
    stats <- do.call(get_summary_stat, c(list(object, fun = c(get_posterior_quantiles, gelman_diag_robust, effectiveSize), probs = probs,
                     stat_name = c(paste0(probs*100, "%"), "Rhat", "ESS"), selection = select), dots))
    for(i in 1:length(stats)){
      stat <- round(stats[[i]], digits)
      stat[,ncol(stat)] <- round(stat[,ncol(stat)])
      if(select == "alpha" & dots$by_subject){
        if(i == 1){
          cat("\n", paste0(select, " ", names(stats)[i]), "\n")
          print(stat)
        }
      } else{
        if(length(stats) > 1){
          cat("\n", paste0(select, " ", names(stats)[i]), "\n")
        } else{
          cat("\n", names(stats)[i], "\n")
        }
        print(stat)
      }
      out_list[[names(stats)[i]]] <- stat
    }
  }
  return(invisible(out_list))
}


#' Plot function for emc objects
#'
#' Makes trace plots for model parameters.
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

plot.emc <- function(x, stage = "sample", selection = c("mu", "sigma2", "alpha"),
                     layout=NA, ...){
  if(attr(x[[1]], "variant_funs")$type == "single"){
    selection <- "alpha"
  }
  for(select in selection){
    plot_mcmc_list(x, stage = stage, selection = select, layout = layout, ...)
  }
}



#' Generate posterior predictives
#'
#' Simulate ``n_post`` data sets using the posterior parameter estimates
#'
#' @param object An emc object from which posterior predictives should
#' be generated
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
predict.emc <- function(object,hyper=FALSE,n_post=100,n_cores=1,
                        stat=c("random","mean","median")[1], ...)
{
  # #' @param force_direction Boolean, take censor direction from argument not samples (default FALSE)
  # #' @param force_response Boolean, take censor response from argument not samples (default FALSE)
  # #' @param LCresponse Boolean, default TRUE, if false set LC response to NA
  # #' @param UCresponse Boolean, default TRUE, if false set UC response to NA
  # #' @param LCdirection Boolean, default TRUE, set LC rt to 0, else to NA
  # #' @param UCdirection Boolean, default TRUE, set LC rt to Inf, else to NA
  # #' @param expand Integer. Default is 1, exact same design for each subject. Larger values will replicate designs, so more trials per subject.
  emc <- object
  dots <- list(...)
  data <- get_data(emc)
  design <- attr(emc,"design_list")
  if(is.null(data$subjects)){
    jointModel <- TRUE
    all_samples <- emc
  } else{
    jointModel <- FALSE
    data <- list(data)
  }
  post_out <- vector("list", length = length(data))
  for(j in 1:length(data)){
    if(jointModel) emc <- single_out_joint(all_samples, j)
    subjects <- levels(data[[j]]$subjects)

    if (hyper) {
      pars <- vector(mode="list",length=n_post)
      # for (i in 1:n_post) {
      #   pars[[i]] <- get_prior_samples(emc,selection="alpha",
      #                                  stage=stage,thin=thin,filter=filter,n_prior=length(subjects))
      #   row.names(pars[[i]]) <- subjects
      # }
    } else {
      dots$selection <- "alpha"; dots$merge_chains <- TRUE; dots$by_subject <- TRUE
      samps <- do.call(get_pars, c(list(emc), fix_dots(dots, get_pars)))
      if (stat != "random") {
        p <- do.call(rbind, lapply(samps, function(x) apply(x[[1]], 2, stat)))
      }
      pars <- vector(mode="list",length=n_post)
      for (i in 1:n_post) {
        if (stat != "random") pars[[i]] <- p else {
          pars[[i]] <- do.call(rbind,lapply(samps,function(x){x[[1]][sample(1:nrow(x[[1]]),1),]}))
        }
      }
    }
    if (n_cores==1) {
      simDat <- vector(mode="list",length=n_post)
      for (i in 1:n_post) {
        simDat[[i]] <- do.call(make_data, c(list(pars[[i]],design=design[[j]],data=data[[j]]), fix_dots(dots, make_data)))
      }
    } else {
      simDat <- mclapply(1:n_post,function(i){
        do.call(make_data, c(list(pars[[i]],design=design[[j]],data=data[[j]]), fix_dots(dots, make_data)))
      },mc.cores=n_cores)
    }
    if (!is.null(attr(simDat[[1]],"adapt"))) adapt <- attr(simDat[[1]],"adapt")
    out <- cbind(postn=rep(1:n_post,times=unlist(lapply(simDat,function(x)dim(x)[1]))),do.call(rbind,simDat))
    if (!is.null(attr(simDat[[1]],"adapt"))) attr(out,"adapt") <- adapt
    if (n_post==1) pars <- pars[[1]]
    attr(out,"pars") <- pars
    post_out[[j]] <- out
  }
  if(!jointModel) post_out <- post_out[[1]]
  return(post_out)
}


# custom s3 funcs ---------------------------------------------------------
#' @rdname check
#' @export
check.emc <- function(emc, selection = c('mu', 'sigma2', 'alpha'), digits = 3,
                      plot_worst = TRUE, ...){
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  dots <- list(...)
  out_list <- list()
  cat("Iterations:\n")
  print(chain_n(emc))
  if(attr(emc[[1]], "variant_funs")$type == "single") selection <- "alpha"
  if(plot_worst){
    mfrow <- coda_setmfrow(Nchains = length(emc), Nparms = length(selection),nplots = 1)
    par(mfrow = mfrow)
  }
  for(select in selection){
    dots$flatten <- ifelse(select == "alpha", FALSE, TRUE)
    dots$by_subject <- TRUE
    ESS <- do.call(ess_summary, c(list(emc, selection = select, stat= NULL), fix_dots(dots, ess_summary)))
    gds <- do.call(gd_summary, c(list(emc, selection = select, stat= NULL), fix_dots(dots, gd_summary)))
    out <- list()
    max_gd <- -Inf
    for(name in names(ESS)){
      combined <- rbind(round(gds[[name]], digits), round(ESS[[name]]))
      rownames(combined) <- c("Rhat", "ESS")
      out[[name]] <- combined
      if(max(gds[[name]]) > max_gd){
        max_gd <- max(gds[[name]])
        cur_max <- name
        max_par <- names(gds[[name]])[which.max(gds[[name]])]
      }
    }
    if(length(ESS) > 1){
      cat("\n", paste0(select, " highest Rhat : ", cur_max), "\n")
    } else{
      cat("\n", cur_max, "\n")
    }
    print(out[[cur_max]])
    if(plot_worst){
      cur_dots <- dots
      if(select == "alpha"){
        cur_dots$subject <- cur_max
        cur_dots$use_par <- max_par
        cur_dots$by_subject <- TRUE
        MCMCs <- do.call(get_pars, c(list(emc, selection = select), fix_dots(cur_dots, get_pars)))
        names(MCMCs) <- paste0("alpha : ", names(MCMCs))
      } else{
        cur_dots$use_par <- max_par
        MCMCs <- do.call(get_pars, c(list(emc, selection = select), fix_dots(cur_dots, get_pars)))
      }
      cur_dots <- add_defaults(cur_dots, xlab = names(MCMCs)[1], ylab = "Highest Rhat parameter")
      do.call(plot, c(list(MCMCs[[1]], auto.layout = FALSE, density = FALSE, ask = FALSE,smooth = FALSE),
                      fix_dots_plot(cur_dots)))
      legend("topleft",legend=paste0("Rhat : ",round(max_gd,digits)), bty = "n")
      legend("topright",legend=paste0("ESS : ", round(ESS[[cur_max]][max_par])), bty = "n")
    }
    out_list[[select]] <- out
  }
  return(invisible(out_list))
}

#' Convergence checks for an emc object
#'
#' Runs a series of convergence checks, prints statistics to the console, and
#' makes traceplots of the worst converged parameter per selection.
#'
#' Note that the `Rhat` is calculated by doubling the number of chains by
#' first splitting chains into first and second half, so it also a test of
#' stationarity.
#'
#' Efficiency of sampling is indicated by the effective
#' sample size (ESS) (from the `coda` R package).
#' Full range of possible samples manipulations described in `get_pars`.
#'
#' @param emc An emc object
#' @param selection A Character vector. Indicates which parameter types to check (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param digits Integer. How many digits to round the ESS and Rhat to in the plots
#' @param plot_worst Boolean. If `TRUE` also plots the chain plots for the worst parameter
#' @param ... Optional arguments that can be passed to `get_pars` or `plot.default` (see `par()`)
#'
#' @return a list with the statistics for the worst converged parameter per selection
#' @examples
#' check(samples_LNR)
#' @export
check <- function(emc, ...){
  UseMethod("check")
}


#' @rdname parameters
#' @export
parameters.emc <- function(emc,selection = "mu", N = NULL, resample = FALSE, ...)
  # extracts and stacks chains into a matrix
{
  dots <- list(...)

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
  dots$merge_chains <- TRUE ; dots$return_mcmc <- FALSE
  dots$flatten <- TRUE; dots$length.out <- N/length(emc)
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


#' Returns a parameter type from an emc object as a data frame.
#'
#' @param emc An emc object
#' @param selection String designating parameter type (e.g. mu, sigma2, correlation, alpha)
#' @param N Integer. How many samples to take from the posterior. If `NULL` will return the full posterior
#' @param resample Boolean. If `TRUE` will sample N samples from the posterior with replacement
#' @param ... Optional arguments that can be passed to `get_pars`
#' @return A data frame with one row for each sample (with a subjects column if selection = "alpha")
#' @export
parameters <- function(emc, ...){
  UseMethod("parameters")
}


#' @rdname fit
#' @export

fit.emc <- function(emc, stage = NULL, iter = 1000, stop_criteria = NULL,report_time=TRUE,
                    p_accept = .8, step_size = 100, verbose = TRUE, verboseProgress = FALSE, fileName = NULL,
                    particles = NULL, particle_factor=50, cores_per_chain = 1,
                    cores_for_chains = length(emc), max_tries = 20, n_blocks = 1,
                    ...){

  if (report_time) start_time <- Sys.time()
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
  stop_criteria <- mapply(get_stop_criteria, stages_names, stop_criteria, MoreArgs = list(type = attr(emc[[1]], "variant_funs")$type))
  if(is.null(stop_criteria[["sample"]]$iter)) stop_criteria[["sample"]]$iter <- iter
  names(stop_criteria) <- stages_names
  if (is.character(emc)) {
    emc <- fix_fileName(emc)
    if(is.null(fileName)) fileName <- emc
    emc <- loadRData(emc)
  }
  if(is.null(stage)){
    nstage <- colSums(chain_n(emc))
    if(all(nstage == 0)){
      stage <- "preburn"
    } else{
      has_ran <- nstage[nstage != 0]
      stage <- names(has_ran)[length(has_ran)]
    }
  }
  if(stage == "preburn"){
    emc <- run_emc(emc, stage = "preburn", stop_criteria[['preburn']], cores_for_chains = cores_for_chains, p_accept = p_accept,
                             step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                             fileName = fileName,
                             particles = particles, particle_factor =  particle_factor,
                             cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  }

  if(any(stage %in% c("preburn", "burn"))){
    emc <-  run_emc(emc, stage = "burn", stop_criteria[['burn']], cores_for_chains = cores_for_chains, p_accept = p_accept,
                              step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName,
                              particles = particles, particle_factor =  particle_factor,
                              cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  }
  if(any(stage %in% c("preburn", "burn", "adapt"))){
    emc <-  run_emc(emc, stage = "adapt", stop_criteria[['adapt']],cores_for_chains = cores_for_chains, p_accept = p_accept,
                              step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName,
                              particles = particles, particle_factor =  particle_factor,
                              cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  }
  if(any(stage %in% c("preburn", "burn", "adapt", "sample")) ){
    emc <-  run_emc(emc, stage = "sample",stop_criteria[['sample']],cores_for_chains = cores_for_chains, p_accept = p_accept,
                              step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName,
                              particles = particles, particle_factor = particle_factor,
                              cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  }
  if (report_time) print(Sys.time()-start_time)
  return(emc)
}

#' Model estimation in EMC2
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
#' @param p_accept A double. The target acceptance probability of the MCMC process.
#' This fine-tunes the width of the search space to obtain the desired acceptance probability. Defaults to .8
#' @param step_size An integer. After each step, the stopping requirements as specified
#' by ``stop_criteria`` are checked and proposal distributions are updated. Defaults to 100.
#' @param verbose Logical. Whether to print messages between each step with the current status regarding the ``stop_criteria``.
#' @param verboseProgress Logical. Whether to print a progress bar within each step or not.
#' Will print one progress bar for each chain and only if ``cores_for_chains = 1``.
#' @param fileName A string. If specified, will auto-save emc object at this location on every iteration.
#' @param particles An integer. How many particles to use, default is `NULL` and
#' ``particle_factor`` is used instead. If specified, ``particle_factor`` is overwritten.
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
#' @param n_blocks An integer. Number of blocks. Will block the parameter chains such that they are
#' updated in blocks. This can be helpful in extremely tough models with a large number of parameters.
#' @param stop_criteria A list. Defines the stopping criteria and for which types
#' of parameters these should hold. See the details and examples section.
#' @param report_time Boolean. If `TRUE`, the time taken to run the MCMC chains till completion of the `stop_criteria` will be printed.
#' @param ... Additional optional arguments
#' @return An emc object
#' @examples \dontrun{
#' # First define a design
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then make the emc object, we've omitted a prior here for brevity so default priors will be used.
#' emc_forstmann <- make_emc(forstmann, design)
#'
#' # With the emc object we can start sampling by simply calling fit
#' emc_forstmann <- fit(emc_forstmann, fileName = "intermediate_save_location.RData")
#'
#' # For particularly hard models it pays off to increase the ``particle_factor``
#' # and, although to a lesser extent, lower ``p_accept``.
#' emc_forstmann <- fit(emc_forstmann, particle_factor = 100, p_accept = .6)
#'
#' # Example of how to use the stop_criteria:
#' emc_forstmann <- fit(emc_forstmann, stop_criteria = list(mean_gd = 1.1, max_gd = 1.5,
#'             selection = c('alpha', 'sigma2'), omit_mpsrf = TRUE, min_es = 1000))
#' # In this case the stop_criteria are set for the sample stage, which will be
#' # run until the mean_gd < 1.1, the max_gd < 1.5 (omitting the multivariate psrf)
#' # and the effective sample size > 1000,
#' # for both the individual-subject parameters ("alpha")
#' # and the group-level variance parameters.
#'
#' # For the unspecified stages in the ``stop_criteria`` the default values
#' # are assumed which are found in Stevenson et al. 2024 <doi.org/10.31234/osf.io/2e4dq>
#'
#' # Alternatively, you can also specify the stop_criteria for specific stages by creating a
#' # nested list
#' emc_forstmann <- fit(emc_forstmann, stop_criteria = list("burn" = list(mean_gd = 1.1, max_gd = 1.5,
#'             selection = c('alpha')), "adapt" = list(min_unique = 100)))
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
  type <- attr(emc[[1]], "variant_funs")$type
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


#' Recovery plots
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
                          N = 1e4, prior_plot_args = list(), ...){
  dots <- add_defaults(list(...), flatten = TRUE)
  type <- attr(emc[[1]], "variant_funs")$type
  if (length(emc[[1]]$data)==1) selection <- "alpha"
  if(selection == "alpha" & type != "single") stop("For savage-dickey ratio, selection cannot be alpha for hierarchical models")
  prior <- emc[[1]]$prior


  psamples <-  get_objects(design = attr(emc,"design_list"),
                           type = attr(emc[[1]], "variant_funs")$type, sample_prior = T,
                           selection = selection, N = N, sampler = emc)
  psamples <- do.call(get_pars, c(list(psamples, selection = selection, merge_chains = TRUE, return_mcmc = FALSE, by_subject = TRUE,
                                          type = attr(emc[[1]], "variant_funs")$type),
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
    prior_plot_args <- add_defaults(prior_plot_args, cex = 2, col = "red", lwd = 1.5)
    dots <- add_defaults(dots, cex = 2, col = "black", lwd = 1.5, main = "Prior and posterior density")
    do.call(plot, c(list(post_density), fix_dots_plot(dots)))
    do.call(lines, c(list(pdensity), fix_dots_plot(prior_plot_args)))
    do.call(points, c(list(H0 , post_dfun(H0)), fix_dots_plot(dots)))
    do.call(points, c(list(H0, pdfun(H0)), fix_dots_plot(prior_plot_args)))
    legend("topright", legend = c("posterior", "prior"), pch = c(1, 1), col = c(dots$col, prior_plot_args$col))
  }
  return(pdfun(H0)/post_dfun(H0))
}

#' Within-model hypothesis testing
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

#' Posterior credible interval tests
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
#' @examples \dontrun{
#' # Run a credible interval test (Bayesian ''t-test'')
#' # Here the full model is an emc object with the hypothesized effect
#' design_full <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#'
#' full_model <- make_emc(forstmann, design_full)
#' full_model <- fit(full_model)
#' credible(full_model, x_name = "v")
#' # We can also compare between two sets of emc objects

#' # Now without a ~ E
#' design_null <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#'
#' null_model <- make_emc(forstmann, design_null)
#' null_model <- fit(null_model)
#' credible(x = null_model, x_name = "a", y = full_model, y_name = "a")
#'
#' # Or provide custom functions
#' credible(x = full_model, x_fun = function(d) d["a_Eaccuracy"] - d["a_Eneutral"])
#' }
#' @return Invisible results table with no rounding.
#' @export
credible <- function(x, ...){
  UseMethod("credible")
}

#' Shorten an emc object
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
  design_list <- attr(x, "design_list")
  x <- lapply(x, remove_samples, stage = stage, filter = filter,
                thin = thin, length.out = length.out, keep_stages = keep_stages)
  attr(x, "design_list") <- design_list
  class(x) <- "emc"
  return(x)
}


#' @rdname gd_summary
#' @export
gd_summary.emc <- function(emc,selection="mu", omit_mpsrf = TRUE,
                           stat = "max", stat_only = FALSE, digits = 3, ...){
  out <- get_summary_stat(emc, selection, gelman_diag_robust, stat = stat,
                          stat_only = stat_only, digits = digits, omit_mpsrf = omit_mpsrf, ...)
  return(out)
}

#' @rdname ess_summary
#' @export
ess_summary.emc <- function(emc,selection="mu", stat = "min", stat_only = FALSE,
                           digits = 1, ...){
  out <- get_summary_stat(emc, selection, effectiveSize,
                          stat = stat, stat_only = stat_only, digits = digits, ...)
  return(out)
}

#' @rdname posterior_summary
#' @export
posterior_summary.emc <- function(emc, selection="mu", probs = c(0.025, .5, .975),
                                  digits = 3, ...){
  out <- get_summary_stat(emc, selection, get_posterior_quantiles,
                          probs = probs, digits = digits, ...)
  return(out)
}

#' Gelman-Rubin statistic
#'
#' Returns the Gelman-Rubin diagnostics (otherwise known as the R-hat) of the selected parameter type;
#' i.e. the ratio of between to within MCMC chain variance.
#'
#' See: Gelman, A and Rubin, DB (1992)
#' Inference from iterative simulation using multiple sequences, *Statistical Science*, 7, 457-511.
#'
#' Full range of possible samples manipulations described in `get_pars`.
#'
#' @param emc An emc object
#' @param selection A Character vector. Indicates which parameter types to check (e.g., `alpha`, `mu`, `sigma2`, `correlation`).
#' @param omit_mpsrf Boolean. If `TRUE` also returns the multivariate point scale reduction factor (see `?coda::gelman.diag`).
#' @param stat A string. Should correspond to a function that can be applied to a vector,
#' which will be performed on the vector/rows or columns of the matrix of the parameters
#' @param stat_only Boolean. If `TRUE` will only return the result of the applied stat function,
#' otherwise returns both the stat result and the result of the function on all parameters.
#' @param digits Integer. How many digits to round the output to
#' @param ... Optional additional arguments that can be passed to `get_pars`
#'
#' @return A matrix or vector of R-hat values for the selected parameter type.
#' @export
#'
#' @examples
#' gd_summary(samples_LNR, selection = "correlation", stat = "mean", flatten = TRUE)
gd_summary <- function(emc, ...){
  UseMethod("gd_summary")
}

#' Effective sample size
#'
#' Returns the effective sample size (ESS) of the selected parameter type.
#' Full range of possible samples manipulations described in `get_pars`.
#'
#' @inheritParams gd_summary.emc
#' @return A matrix or vector of ESS values for the selected parameter type.
#' @export
#'
#' @examples
#' ess_summary(samples_LNR, selection = "alpha")
ess_summary <- function(emc, ...){
  UseMethod("ess_summary")
}

#' Posterior quantiles
#'
#' Returns the quantiles of the selected parameter type.
#' Full range of possible samples manipulations described in `get_pars`.
#'
#' @inheritParams gd_summary.emc
#' @param probs A vector. Indicates which quantiles to return from the posterior.
#' @return A list of posterior quantiles for each parameter group in the selected parameter type.
#' @export
#'
#' @examples
#' posterior_summary(samples_LNR)
posterior_summary <- function(emc, ...){
  UseMethod("posterior_summary")
}


#' @rdname get_data
#' @export
get_data.emc <- function(emc) {
  if(is.null(emc[[1]]$data[[1]]$subjects)){ # Joint model
    dat <- vector("list", length(emc[[1]]$data[[1]]))
    for(i in 1:length(dat)){
      design <- attr(emc, "design_list")[[i]]
      tmp <- lapply(emc[[1]]$data,\(x) x[[i]][attr(x[[i]],"expand"),])
      tmp <- do.call(rbind, tmp)
      row.names(tmp) <- NULL
      tmp <- tmp[tmp$lR == levels(tmp$lR)[1],]
      tmp <- tmp[,!(colnames(tmp) %in% c("trials","lR","lM", "winner", "SlR", "RACE", names(design$Ffunctions)))]
      dat[[i]] <- tmp
    }
  } else{
    design <- attr(emc, "design_list")[[1]]
    dat <- do.call(rbind,lapply(emc[[1]]$data,\(x) x[attr(x,"expand"),]))
    row.names(dat) <- NULL
    dat <- dat[dat$lR == levels(dat$lR)[1],]
    dat <- dat[,!(colnames(dat) %in% c("trials","lR","lM","winner", "SlR", "RACE", names(design$Ffunctions)))]
  }
  return(dat)
}
#' Get data
#'
#' Extracts data from an emc object
#'
#' @details
#' emc adds columns and rows to a dataframe in order to facilitate efficient likelihood calculations.
#' This function will return the data as provided originally.
#'
#' @param emc an emc object
#'
#' @return A dataframe of the original data
#'
#' @export
#'
#' @examples
#' get_data(samples_LNR)
get_data <- function(emc){
  UseMethod("get_data")
}
