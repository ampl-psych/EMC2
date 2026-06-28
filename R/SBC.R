#' Simulation-Based Calibration
#'
#' Runs SBC for an EMC2 model and associated design. Returns
#' normalized rank (between 0 and 1) and prior samples. For hierarchical models the group-level mean and
#' the (implied) group-level (co-)variance are returned.
#' For non-hierarchical models only the subject-level parameters rank is returned.
#'
#' @param design_in An emc design list. The design of the model to be used in SBC
#' @param prior_in An emc prior list. The prior for the design to be used in SBC
#' @param replicates Integer. The number of samples to draw from the prior
#' @param trials Integer. The number of trials of the simulated data (per subject)
#' @param n_subjects Integer. Only used for hierarchical models. The number of subjects to be used in data generation of each replicate
#' @param plot_data Boolean. Whether to plot the data simulated (aggregated across subjects)
#' @param verbose Verbose. Whether to print progress related messages
#' @param fileName Character. Highly recommended, saves temporary results to the fileName
#' @param ... A list of optional additional arguments that can be passed to `fit` and `make_emc`
#'
#' @return The ranks and prior samples. For hierarchical models also the prior-generated subject-level parameters.
#' @export
run_sbc <- function(design_in, prior_in, replicates = 250, trials = 100, n_subjects = 30,
                    plot_data = FALSE, verbose = TRUE,
                    fileName = NULL, ...){
  if(is.null(fileName)) message("Since SBC can take a while it's highly recommended to specify a fileName to save temporary results in case of crashes")
  type <- attr(prior_in, "type")
  if(type == "single"){
    out <- SBC_single(design_in, prior_in, replicates, trials,
                      plot_data, verbose, fileName, ...)
  } else{
    if(plot_data) warning("Hierarchical SBC does not support plotting")
    out <- SBC_hierarchical_parallel(design_in, prior_in, replicates, trials, n_subjects,
                                     verbose, fileName, ...)
  }
  invisible(out)
}


#' Recover Partial or Interrupted SBC Results
#'
#' Reassembles the object normally returned by `run_sbc()` from the intermediate
#' per-replicate files it writes, without running any further replicates. Use
#' this when a run was interrupted (e.g. a node reboot) or deliberately aborted
#' and you want to inspect the replicates completed so far (for example with
#' `plot_sbc_ecdf()` or `plot_sbc_hist()`).
#'
#' `run_sbc()` writes one file per completed replicate (in
#' `<fileName-without-extension>_temp/rep_<i>.rds`, and persistent
#' `<base>_rep_<i>.rds` files next to `fileName`) plus the prior draws
#' (`prior_samples.rds` in the temp directory, or `prior_mu`/`prior_var` /
#' `prior_alpha` in `fileName`). This function gathers whichever replicate files
#' survive and assembles them exactly as a completed run would.
#'
#' To instead *continue* an interrupted run to the target number of replicates,
#' simply re-issue the original `run_sbc()` call with the same `fileName`: it
#' resumes from the completed replicates and runs only the missing ones. Note
#' that faithful resumption requires the original prior draws (`prior_samples.rds`
#' or the prior saved in `fileName`); if only the `rep_<i>.rds` files survive,
#' the results can still be inspected here but the run cannot be resumed (a new
#' run would re-pair replicate indices with different true parameters).
#'
#' @param tempdir Character. The directory containing the SBC replicate files to
#'   recover. Either the `_temp` directory `run_sbc()` writes (`rep_<i>.rds`,
#'   plus `prior_samples.rds`), or a directory holding the persistent
#'   `<base>_rep_<i>.rds` files; the run's base name is auto-detected in the
#'   latter case. Must be supplied.
#' @param design The emc design used for the run. Used to name the recovered
#'   parameters (and, for hierarchical recovery, the group-mean `mu` columns).
#'   Must be supplied. May alternatively be a character vector of the parameter
#'   names directly.
#' @param fileName Character or NULL. Optional save location, behaving as in
#'   `run_sbc()`: the recovered object is written there in `run_sbc()`'s format
#'   (so the file is indistinguishable from a normal run that finished at the
#'   recovered number of replicates). When `NULL` (the default) nothing is
#'   written to disk; the result is only returned, so the recovered object is
#'   never placed in `tempdir` unless you explicitly point `fileName` there.
#' @param prior_in Optional emc prior list. Currently only used for validation;
#'   the prior draws are read from disk, not redrawn.
#' @param type Character. `"auto"` (default) detects single vs hierarchical from
#'   the rep files; otherwise force `"single"` or `"hierarchical"`.
#' @param verbose Logical. Whether to print progress/recovery messages.
#'
#' @return Invisibly, the same structure `run_sbc()` returns (for single fits an
#'   SBC list of rank/med/bias/coverage; for hierarchical fits a list with
#'   `rank`, `prior` and `rand_effects`), as if the run had finished at the
#'   recovered number of replicates. The integer replicate indices included are stored in
#'   the `"recovered_reps"` attribute. For hierarchical recovery when the prior
#'   draws are unavailable, `$prior` is `NULL` (ranks are still fully recovered).
#'   When `fileName` is supplied the object is also saved there: for single fits
#'   as `SBC` (plus `prior_alpha` when available), for hierarchical fits as
#'   `SBC_temp`, matching `run_sbc()`.
#' @export
recover_sbc <- function(tempdir, design, fileName = NULL, prior_in = NULL,
                        type = c("auto", "single", "hierarchical"),
                        verbose = TRUE) {
  type <- match.arg(type)
  if (missing(tempdir) || is.null(tempdir) || !dir.exists(tempdir))
    stop("`tempdir` must be an existing directory containing the SBC replicate ",
         "files; got: ", if (missing(tempdir)) "<missing>" else tempdir)
  if (missing(design) || is.null(design))
    stop("`design` (the emc design used for the run) must be supplied.")

  # --- gather replicate files from `tempdir`: rep_<i>.rds (temp-style) first,
  #     then auto-detected <base>_rep_<i>.rds (persistent style). ---
  reps      <- .sbc_read_rep_dir(tempdir)
  from_temp <- length(reps) > 0
  if (length(reps) == 0)
    reps <- .sbc_load_persist_reps(tempdir, persist_base = NULL)
  if (length(reps) == 0)
    stop("No replicate files found in `tempdir` (", tempdir,
         "); expected rep_<i>.rds or <base>_rep_<i>.rds files.")

  # --- prior draws: prior_samples.rds in tempdir, else (if given) from fileName ---
  prior_obj <- NULL
  ps_file   <- file.path(tempdir, "prior_samples.rds")
  if (file.exists(ps_file))
    prior_obj <- tryCatch(readRDS(ps_file), error = function(e) NULL)
  if (is.null(prior_obj) && !is.null(fileName) && file.exists(fileName)) {
    env <- new.env(parent = emptyenv())
    if (isTRUE(tryCatch({ load(fileName, envir = env); TRUE }, error = function(e) FALSE))) {
      if (all(c("prior_mu", "prior_var") %in% ls(env)))
        prior_obj <- list(prior_mu = env$prior_mu, prior_var = env$prior_var)
      else if ("prior_alpha" %in% ls(env))
        prior_obj <- env$prior_alpha
    }
  }

  # --- detect type from a rep file if not forced ---
  if (type == "auto") {
    first <- reps[[1]]
    type <- if (!is.null(first[["rank_mu_row"]])) "hierarchical"
            else if (!is.null(first[["rank"]]) || isTRUE(first[["failed"]])) "single"
            else stop("Could not auto-detect SBC type from rep files; pass type=")
  }

  par_names <- if (is.character(design)) design else names(sampled_pars(design))

  if (type == "single") {
    out    <- .sbc_assemble_single(reps, par_names = par_names)
    target <- if (is.matrix(prior_obj)) nrow(prior_obj) else NA_integer_
  } else {
    prior_mu <- prior_var <- NULL
    if (is.list(prior_obj) && !is.null(prior_obj[["prior_mu"]])) {
      prior_mu  <- prior_obj[["prior_mu"]]
      prior_var <- prior_obj[["prior_var"]]
    }
    out    <- .sbc_assemble_hierarchical(reps, prior_mu, prior_var, par_names)
    target <- if (!is.null(prior_mu)) ncol(prior_mu) else NA_integer_
    if (is.null(prior_mu))
      warning("Prior draws not found (no prior_samples.rds in '", tempdir,
              "'); $prior is NULL. Ranks are fully recovered, but the ",
              "run cannot be faithfully resumed.")
  }

  # --- save only if a fileName is given (run_sbc format) ---
  if (!is.null(fileName)) {
    out_dir <- dirname(fileName)
    if (nzchar(out_dir) && !dir.exists(out_dir))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (type == "single") {
      SBC <- out
      if (is.matrix(prior_obj)) {
        prior_alpha <- prior_obj
        save(SBC, prior_alpha, file = fileName)
      } else {
        save(SBC, file = fileName)
      }
    } else {
      SBC_temp <- out
      save(SBC_temp, file = fileName)
    }
  }

  used <- attr(out, "recovered_reps")
  if (verbose) {
    n_done <- length(used)
    msg <- paste0("Recovered ", n_done, " ", type, " replicate(s)")
    if (!is.na(target)) {
      msg <- paste0(msg, " of ", target, " requested")
      missing <- setdiff(seq_len(target), used)
      if (length(missing))
        msg <- paste0(msg, "; missing: ",
                      paste(utils::head(missing, 20), collapse = ","),
                      if (length(missing) > 20) ",..." else "")
    }
    saved <- if (!is.null(fileName)) paste0(" Saved to ", fileName, ".")
             else " Not saved (supply fileName to save)."
    resumable <- from_temp && file.exists(ps_file)
    message(msg, ".", saved, " ",
            if (resumable)
              "Re-run your original run_sbc() call with the same fileName to resume."
            else
              "Inspect-only: prior draws for resumption were not found in the temp directory.")
  }
  invisible(out)
}



# Internal helpers --------------------------------------------------------

.sbc_load_persist_reps <- function(persist_dir, persist_base = NULL) {
  if (is.null(persist_dir) || !dir.exists(persist_dir)) return(list())
  if (is.null(persist_base)) {
    # auto-detect the run's base name from <base>_rep_<i>.rds files present
    cand  <- list.files(persist_dir, pattern = "_rep_[0-9]+\\.rds$")
    bases <- unique(sub("_rep_[0-9]+\\.rds$", "", cand))
    if (length(bases) == 0L) return(list())
    if (length(bases) > 1L)
      stop("Multiple SBC runs found in '", persist_dir, "' (bases: ",
           paste(bases, collapse = ", "),
           "); point `tempdir` at a directory with a single run's rep files.")
    persist_base <- bases
  }
  pattern <- paste0("^", persist_base, "_rep_[0-9]+\\.rds$")
  files   <- list.files(persist_dir, pattern = pattern, full.names = TRUE)
  out <- list()
  for (f in files) {
    result <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(result)) {
      nm       <- sub(paste0(".*", persist_base, "_rep_([0-9]+)\\.rds$"), "\\1", f)
      out[[nm]] <- result
    }
  }
  out
}

.sbc_cleanup_persist_reps <- function(persist_dir, persist_base) {
  if (is.null(persist_dir) || is.null(persist_base)) return(invisible(NULL))
  pattern <- paste0("^", persist_base, "_rep_[0-9]+\\.rds$")
  files   <- list.files(persist_dir, pattern = pattern, full.names = TRUE)
  if (length(files)) file.remove(files)
  invisible(NULL)
}

# Read all rep_<i>.rds in a temp directory into a list keyed by the (character)
# replicate index. Unreadable files are skipped.
.sbc_read_rep_dir <- function(dir) {
  out <- list()
  if (is.null(dir) || !dir.exists(dir)) return(out)
  files <- list.files(dir, pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)
  for (f in files) {
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(r)) {
      i <- sub(".*rep_([0-9]+)\\.rds$", "\\1", basename(f))
      out[[i]] <- r
    }
  }
  out
}

# Assemble single (non-hierarchical) SBC replicate results into the object
# returned by SBC_single(). `reps` is a list keyed by character replicate index
# (as produced by the restart loaders or .sbc_read_rep_dir); NULL and failed
# replicates are skipped. The included replicate indices are recorded in the
# "recovered_reps" attribute. `par_names`, if supplied, overrides the (already
# named) per-replicate vectors.
.sbc_assemble_single <- function(reps, par_names = NULL) {
  idx  <- sort(as.integer(names(reps)))
  reps <- reps[as.character(idx)]
  ok   <- vapply(reps, function(r)
    !is.null(r) && !isTRUE(r$failed) && !is.null(r$rank), logical(1))
  if (!any(ok)) stop("No completed replicates found to assemble")
  used_idx <- idx[ok]
  SBC <- split_list_to_dfs(reps[ok])
  if (!is.null(par_names)) {
    for (comp in names(SBC)) colnames(SBC[[comp]][["alpha"]]) <- par_names
  }
  attr(SBC, "recovered_reps") <- used_idx
  SBC
}

# Assemble hierarchical SBC replicate results into the object returned by
# SBC_hierarchical_parallel(). NULL replicates are skipped; included indices are
# recorded in "recovered_reps". prior_mu/prior_var may be NULL (then $prior is
# NULL). rank_var column names come from the replicates' stored var_col_names;
# rank_mu column names come from par_names (or rownames(prior_mu) as a fallback).
.sbc_assemble_hierarchical <- function(reps, prior_mu = NULL, prior_var = NULL,
                                       par_names = NULL) {
  idx  <- sort(as.integer(names(reps)))
  reps <- reps[as.character(idx)]
  ok   <- vapply(reps, function(r)
    !is.null(r) && !is.null(r$rank_mu_row), logical(1))
  if (!any(ok)) stop("No completed replicates found to assemble")
  used_idx <- idx[ok]
  reps <- reps[ok]

  rank_mu          <- do.call(rbind, lapply(reps, `[[`, "rank_mu_row"))
  rank_var         <- do.call(rbind, lapply(reps, `[[`, "rank_var_row"))
  all_rand_effects <- lapply(reps, `[[`, "rand_effects")

  if (is.null(par_names) && !is.null(prior_mu)) par_names <- rownames(prior_mu)
  if (!is.null(par_names)) colnames(rank_mu) <- par_names
  colnames(rank_var) <- reps[[1]][["var_col_names"]]

  out <- list(rank         = list(mu = rank_mu, var = rank_var),
              prior        = list(mu = prior_mu, var = prior_var),
              rand_effects = all_rand_effects)
  attr(out, "recovered_reps") <- used_idx
  out
}


SBC_hierarchical <- function(design_in, prior_in, replicates = 250, trials = 100, n_subjects = 30,
                             plot_data = FALSE, verbose = TRUE,
                             fileName = NULL, ...){
  dots <- add_defaults(list(...), max_tries = 50, compress = FALSE, rt_resolution = 1e-12,
                       stop_criteria = list(min_es = 100, max_gd = 1.1,
                                            selection = c("alpha", "mu", "Sigma")))
  dots$verbose <- verbose
  type <- attr(prior_in, "type")
  if(type != "diagonal-gamma") warning("SBC in EMC2 is frequently biased for prior structures with heavy variance mass around 0 \n
                                       please consider using 'type = diagonal-gamma' to test your model. See also vignettes")
  # Draw prior samples
  # Should at a later point go to predict
  prior_mu <- plot(prior_in, design_in, do_plot = F, N = replicates, selection = "mu", return_mcmc = FALSE, map = FALSE)[[1]]
  prior_var <- plot(prior_in, design_in, do_plot = F, N = replicates, selection = "Sigma", return_mcmc = FALSE,
                          remove_constants = FALSE,map=FALSE)[[1]]
  mu_names <- rownames(prior_mu)
  var_names <- rownames(prior_var[, , 1])
  keep <- var_names %in% gsub(":", "_", mu_names)
  prior_var <- prior_var[keep, keep, , drop = FALSE]
  rank_mu <- data.frame()
  rank_var <- data.frame()
  par_names <- names(sampled_pars(design_in))
  if(!is.null(fileName)) save(prior_mu, prior_var, file = fileName)
  all_rand_effects <- vector("list", length = replicates)
  for(i in 1:replicates){
    if(verbose) print(paste0("Sample ", i, " out of ", replicates))
    rand_effects <- make_random_effects(design_in, prior_mu[,i], n_subj = n_subjects, covariances = prior_var[,,i])
    all_rand_effects[[i]] <- rand_effects
    data <- make_data(rand_effects,design_in, trials, model = design_in$model, ...)
    if(plot_data){
      plot_density(data, factors = names(design_in$Ffactors)[names(design_in$Ffactors) != "subjects"])
    }
    emc <-  do.call(make_emc, c(list(data = data, design = design_in, prior_list = prior_in, type = type), fix_dots(dots, make_emc)))
    emc <-  do.call(fit, c(list(emc = emc), fix_dots(dots, fit)))

    ESS_mu <- round(pmin(unlist(ess_summary(emc, selection = "mu", stat = NULL)), chain_n(emc)[1,"sample"]))
    mu_rec <- get_pars(emc, selection = "mu", return_mcmc = F, merge_chains = T,flatten = T)
    rank_mu <- rbind(rank_mu, mapply(get_ranks_ESS, split(mu_rec, row(mu_rec)), t(ESS_mu), prior_mu[,i]))

    ESS_var <- round(pmin(unlist(ess_summary(emc, selection = "Sigma", stat = NULL)), chain_n(emc)[1,"sample"]))
    var_rec <- get_pars(emc, selection = "Sigma", return_mcmc = F, merge_chains = T, flatten = T)
    prior_var_input <- get_pars(emc, selection = "Sigma", flatten = T, true_pars = prior_var[,,i])[[1]][[1]]
    rank_var <- rbind(rank_var, mapply(get_ranks_ESS, split(var_rec, row(var_rec)), t(ESS_var), prior_var_input))

    colnames(rank_mu) <- par_names
    colnames(rank_var) <- colnames(prior_var_input)
    if(!is.null(fileName)){
      SBC_temp <- list(rank = list(mu = rank_mu, var =  rank_var),
                     prior = list(mu = prior_mu, var = prior_var),
                     rand_effects = rand_effects, emc = emc)
      save(SBC_temp, file = fileName)
    }
  }
  return(list(rank = list(mu = rank_mu, var =  rank_var),
              prior = list(mu = prior_mu, var = prior_var),
              rand_effects = rand_effects))
}


sbc_running_counter <- function(i, temp_dir, offset = 0L) {
  if (!is.null(temp_dir)) {
    writeLines("", file.path(temp_dir, paste0("started_", i, ".flag")))
    length(list.files(temp_dir, pattern = "^started_.*\\.flag$")) + offset
  } else {
    i
  }
}

run_SBC_hierarchical_rep <- function(i, design_in, prior_mu, prior_var, trials, n_subjects,
                                     prior_in, type, dots, temp_dir, offset = 0L,
                                     persist_dir = NULL, persist_base = NULL) {
  dots[["cores_per_chain"]] <- 1L
  dots[["verbose"]] <- FALSE
  dots[["verboseProgress"]] <- FALSE
  message("Running data set ", sbc_running_counter(i, temp_dir, offset))
  rand_effects <- make_random_effects(design_in, prior_mu[, i], n_subj = n_subjects,
                                      covariances = prior_var[,, i])
  data <- do.call(make_data, c(list(rand_effects, design_in, trials, model = design_in$model),
                               fix_dots(dots, make_data)))
  emc <- suppressMessages(do.call(make_emc, c(list(data = data, design = design_in, prior_list = prior_in, type = type),
                              fix_dots(dots, make_emc))))
  emc <- do.call(fit, c(list(emc = emc), fix_dots(dots, fit)))

  ESS_mu  <- round(pmin(unlist(ess_summary(emc, selection = "mu",    stat = NULL)), chain_n(emc)[1, "sample"]))
  mu_rec  <- get_pars(emc, selection = "mu",    return_mcmc = FALSE, merge_chains = TRUE, flatten = TRUE)
  rank_mu_row <- unname(mapply(get_ranks_ESS, split(mu_rec, row(mu_rec)), t(ESS_mu), prior_mu[, i]))

  ESS_var         <- round(pmin(unlist(ess_summary(emc, selection = "Sigma", stat = NULL)), chain_n(emc)[1, "sample"]))
  var_rec         <- get_pars(emc, selection = "Sigma", return_mcmc = FALSE, merge_chains = TRUE, flatten = TRUE)
  prior_var_input <- get_pars(emc, selection = "Sigma", flatten = TRUE, true_pars = prior_var[,, i])[[1]][[1]]
  rank_var_row    <- unname(mapply(get_ranks_ESS, split(var_rec, row(var_rec)), t(ESS_var), prior_var_input))

  result <- list(rank_mu_row   = rank_mu_row,
                 rank_var_row  = rank_var_row,
                 var_col_names = colnames(prior_var_input),
                 rand_effects  = rand_effects)
  if (!is.null(temp_dir))
    saveRDS(result, file.path(temp_dir, paste0("rep_", i, ".rds")))
  # Also write to persistent storage (survives node-local temp cleanup)
  if (!is.null(persist_dir) && !is.null(persist_base))
    saveRDS(result, file.path(persist_dir, paste0(persist_base, "_rep_", i, ".rds")))
  result
}


SBC_hierarchical_parallel <- function(design_in, prior_in, replicates = 250, trials = 100,
                                      n_subjects = 30, verbose = TRUE,
                                      fileName = NULL, ...) {
  dots <- add_defaults(list(...), max_tries = 50, compress = FALSE, rt_resolution = 1e-12,
                       stop_criteria = list(min_es = 100, max_gd = 1.1,
                                            selection = c("alpha", "mu", "Sigma")),
                       cores_per_chain = 1)
  dots$verbose <- verbose
  type <- attr(prior_in, "type")
  if (type != "diagonal-gamma")
    warning("SBC in EMC2 is frequently biased for prior structures with heavy variance mass around 0 \n
             please consider using 'type = diagonal-gamma' to test your model. See also vignettes")

  temp_dir <- if (!is.null(fileName))
    paste0(tools::file_path_sans_ext(fileName), "_temp")
  else NULL

  # Persistent storage alongside fileName (survives node-local temp cleanup)
  persist_dir  <- if (!is.null(fileName)) dirname(normalizePath(fileName, mustWork = FALSE)) else NULL
  persist_base <- if (!is.null(fileName)) tools::file_path_sans_ext(basename(fileName)) else NULL

  # --- restart: recover completed replicates ---
  completed_results <- list()
  offset <- 0L
  if (!is.null(temp_dir) && dir.exists(temp_dir)) {
    # Normal restart: temp_dir survived
    prior_samples <- readRDS(file.path(temp_dir, "prior_samples.rds"))
    prior_mu  <- prior_samples$prior_mu
    prior_var <- prior_samples$prior_var
    existing  <- list.files(temp_dir, pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)
    for (f in existing) {
      result <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(result)) {
        file.remove(f)
        if (verbose) message("Deleted unreadable temp file: ", basename(f))
      } else {
        i <- as.integer(sub(".*rep_([0-9]+)\\.rds$", "\\1", f))
        completed_results[[as.character(i)]] <- result
      }
    }
    offset <- length(completed_results)
    file.remove(list.files(temp_dir, pattern = "^started_.*\\.flag$", full.names = TRUE))
    if (verbose)
      message("Restarting: ", offset, " of ", replicates,
              " replicates already complete, running remaining ",
              replicates - offset)
  } else {
    # Check for persistent rep files written alongside fileName (survive node cleanup)
    persist_reps <- .sbc_load_persist_reps(persist_dir, persist_base)
    if (length(persist_reps) > 0) {
      env <- new.env(parent = emptyenv())
      load(fileName, envir = env)
      if (!all(c("prior_mu", "prior_var") %in% ls(env)))
        stop("Cannot restart: prior_mu/prior_var not found in ", fileName)
      prior_mu <- env$prior_mu
      prior_var <- env$prior_var
      dir.create(temp_dir, showWarnings = FALSE)
      saveRDS(list(prior_mu = prior_mu, prior_var = prior_var),
              file.path(temp_dir, "prior_samples.rds"))
      for (nm in names(persist_reps)) {
        saveRDS(persist_reps[[nm]], file.path(temp_dir, paste0("rep_", nm, ".rds")))
        completed_results[[nm]] <- persist_reps[[nm]]
      }
      offset <- length(completed_results)
      if (verbose)
        message("Restarting from persistent rep files: ", offset, " of ", replicates,
                " replicates already complete, running remaining ",
                replicates - offset)
    } else {
      # Check for partial SBC_temp saved in fileName (e.g. from an older sequential run)
      partial_sbc <- NULL
      if (!is.null(fileName) && file.exists(fileName)) {
        env <- new.env(parent = emptyenv())
        loaded_ok <- isTRUE(tryCatch({ load(fileName, envir = env); TRUE }, error = function(e) FALSE))
        if (loaded_ok && "SBC_temp" %in% ls(env) &&
            !is.null(env$SBC_temp$rank$mu) && nrow(env$SBC_temp$rank$mu) > 0 &&
            !is.null(env$SBC_temp$prior$mu) && !is.null(env$SBC_temp$prior$var))
          partial_sbc <- env$SBC_temp
      }
      if (!is.null(partial_sbc)) {
        n_done    <- min(nrow(partial_sbc$rank$mu), replicates)
        prior_mu  <- partial_sbc$prior$mu
        prior_var <- partial_sbc$prior$var
        dir.create(temp_dir, showWarnings = FALSE)
        saveRDS(list(prior_mu = prior_mu, prior_var = prior_var),
                file.path(temp_dir, "prior_samples.rds"))
        var_col_names <- colnames(partial_sbc$rank$var)
        for (k in seq_len(n_done)) {
          result <- list(
            rank_mu_row   = unname(as.numeric(partial_sbc$rank$mu[k, ])),
            rank_var_row  = unname(as.numeric(partial_sbc$rank$var[k, ])),
            var_col_names = var_col_names,
            rand_effects  = NULL
          )
          saveRDS(result, file.path(temp_dir, paste0("rep_", k, ".rds")))
          completed_results[[as.character(k)]] <- result
          if (!is.null(persist_dir) && !is.null(persist_base))
            saveRDS(result, file.path(persist_dir, paste0(persist_base, "_rep_", k, ".rds")))
        }
        offset <- n_done
        # Overwrite fileName with prior_mu/prior_var so the persist_reps branch works on crash-restart
        if (!is.null(fileName)) save(prior_mu, prior_var, file = fileName)
        if (verbose)
          message("Restarting from partial SBC_temp in ", fileName, ": ",
                  n_done, " of ", replicates, " replicates complete, running remaining ",
                  replicates - n_done)
      } else {
        # Fresh run — draw prior samples and set up temp_dir
        prior_mu  <- plot(prior_in, design_in, do_plot = FALSE, N = replicates, selection = "mu",
                          return_mcmc = FALSE, map = FALSE)[[1]]
        prior_var <- plot(prior_in, design_in, do_plot = FALSE, N = replicates, selection = "Sigma",
                          return_mcmc = FALSE, remove_constants = FALSE, map = FALSE)[[1]]
        if (!is.null(temp_dir)) {
          dir.create(temp_dir, showWarnings = FALSE)
          saveRDS(list(prior_mu = prior_mu, prior_var = prior_var),
                  file.path(temp_dir, "prior_samples.rds"))
        }
        if (!is.null(fileName)) save(prior_mu, prior_var, file = fileName)
      }
    }
  }

  par_names    <- names(sampled_pars(design_in))
  missing_reps <- setdiff(1:replicates, as.integer(names(completed_results)))

  if (length(missing_reps) > 0) {
    if (verbose)
      message("Processing ", dots[["cores_per_chain"]], " data sets in parallel")
    res_new <- auto_mclapply(
      X   = missing_reps,
      FUN = run_SBC_hierarchical_rep,
      design_in, prior_mu, prior_var, trials, n_subjects, prior_in, type, dots, temp_dir, offset,
      persist_dir, persist_base,
      mc.cores = dots[["cores_per_chain"]]
    )
    # recover from disk for any workers that returned NULL
    for (idx in seq_along(missing_reps)) {
      i <- missing_reps[[idx]]
      if (is.null(res_new[[idx]])) {
        f_temp <- if (!is.null(temp_dir)) file.path(temp_dir, paste0("rep_", i, ".rds")) else NULL
        f_pers <- if (!is.null(persist_dir) && !is.null(persist_base))
          file.path(persist_dir, paste0(persist_base, "_rep_", i, ".rds")) else NULL
        if (!is.null(f_temp) && file.exists(f_temp))
          res_new[[idx]] <- readRDS(f_temp)
        else if (!is.null(f_pers) && file.exists(f_pers))
          res_new[[idx]] <- readRDS(f_pers)
      }
      if (!is.null(res_new[[idx]]))
        completed_results[[as.character(i)]] <- res_new[[idx]]
    }
  }

  # Assemble in replicate order (shared with recover_sbc)
  out <- .sbc_assemble_hierarchical(completed_results, prior_mu, prior_var, par_names)
  attr(out, "recovered_reps") <- NULL
  if (!is.null(fileName)) {
    SBC_temp <- out
    save(SBC_temp, file = fileName)
    unlink(temp_dir, recursive = TRUE)
    # Clean up persistent rep files now that the main file is complete
    .sbc_cleanup_persist_reps(persist_dir, persist_base)
  }
  return(out)
}


run_SBC_subject <- function(rep, design_in, prior_alpha, trials, prior_in, dots, temp_dir, offset = 0L,
                            persist_dir = NULL, persist_base = NULL){
  dots[["verbose"]] <- FALSE
  dots[["verboseProgress"]] <- FALSE
  message("Running data set ", sbc_running_counter(rep, temp_dir, offset))
  p_vector <- prior_alpha[rep,]
  data <- do.call(make_data, c(list(parameters = p_vector, design = design_in, n_trials = trials), fix_dots(dots, make_data)))
  emc <- suppressMessages(do.call(make_emc, c(list(data = data, design = design_in, prior_list = prior_in, type = "single"), fix_dots(dots, make_emc))))

  p_vector_dir <- if (!is.null(temp_dir)) temp_dir else "."
  fit_result <- tryCatch({
    do.call(fit, c(list(emc = emc), fix_dots(dots, fit)))
  }, warning = function(w) {
    filename <- file.path(p_vector_dir, paste0("p_vector_rep", rep, ".Rdata"))
    save(p_vector, file = filename)
    warning("A warning occurred during fitting for replication ", rep,
            ". The input parameters have been saved as ", filename,
            ". Warning message: ", w$message)
    warning(w)
    return(NULL)
  }, error = function(e) {
    filename <- file.path(p_vector_dir, paste0("p_vector_rep", rep, ".Rdata"))
    save(p_vector, file = filename)
    warning("An error occurred during fitting for replication ", rep,
            ". The input parameters have been saved as ", filename,
            ". Error message: ", e$message)
    return(NULL)
  })

  if (is.null(fit_result))
    return(list(rank = NULL, med = NULL, bias = NULL, coverage = NULL, failed = TRUE))

  emc <- fit_result

  ESS       <- ess_summary(emc, stat = NULL)[[1]]
  alpha_rec <- get_pars(emc, selection = "alpha", return_mcmc = F, merge_chains = T, flatten = T)[,1,]
  rank      <- mapply(get_ranks_ESS, split(alpha_rec, row(alpha_rec)), t(ESS), p_vector)
  names(rank) <- names(sampled_pars(design_in))
  CI       <- credint(emc)[[1]]
  med      <- CI[, 2]
  bias     <- med - p_vector
  coverage <- p_vector > CI[, 1] & p_vector < CI[, 3]

  result <- list(rank = rank, med = med, bias = bias, coverage = coverage)
  if (!is.null(temp_dir))
    saveRDS(result, file.path(temp_dir, paste0("rep_", rep, ".rds")))
  if (!is.null(persist_dir) && !is.null(persist_base))
    saveRDS(result, file.path(persist_dir, paste0(persist_base, "_rep_", rep, ".rds")))
  result
}

split_list_to_dfs <- function(lst, type = "alpha") {
  comps <- names(lst[[1]])   # get component names from the first element
  out <- lapply(comps, function(nm) {
    res <- list()
    res[[type]] <- do.call(rbind, lapply(lst, `[[`, nm))
    return(res)
  })
  names(out) <- comps
  out
}



SBC_single <- function(
  design_in,
  prior_in,
  replicates = 250,
  trials = 100,
  plot_data = FALSE,
  verbose = TRUE,
  fileName = NULL,
  ...
) {
  if (attr(prior_in, "type") != "single") {
    stop("can only use `type = single`")
  }
  dots <- add_defaults(
    list(...),
    max_tries = 50,
    compress = FALSE, rt_resolution = 1e-12,
    stop_criteria = list(
      min_es = 100, max_gd = 1.1, selection = c("alpha", "mu", "Sigma")
    ),
    cores_per_chain = 1
  )
  dots[["verbose"]] <- verbose

  temp_dir <- if (!is.null(fileName))
    paste0(tools::file_path_sans_ext(fileName), "_temp")
  else NULL

  persist_dir  <- if (!is.null(fileName)) dirname(normalizePath(fileName, mustWork = FALSE)) else NULL
  persist_base <- if (!is.null(fileName)) tools::file_path_sans_ext(basename(fileName)) else NULL

  # --- restart: recover completed replicates ---
  completed_results <- list()
  offset <- 0L
  if (!is.null(temp_dir) && dir.exists(temp_dir)) {
    prior_alpha <- readRDS(file.path(temp_dir, "prior_samples.rds"))
    existing    <- list.files(temp_dir, pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)
    for (f in existing) {
      result <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(result)) {
        file.remove(f)
        if (verbose) message("Deleted unreadable temp file: ", basename(f))
      } else {
        i <- as.integer(sub(".*rep_([0-9]+)\\.rds$", "\\1", f))
        completed_results[[as.character(i)]] <- result
      }
    }
    offset <- length(completed_results)
    file.remove(list.files(temp_dir, pattern = "^started_.*\\.flag$", full.names = TRUE))
    if (verbose)
      message("Restarting: ", offset, " of ", replicates,
              " replicates already complete, running remaining ",
              replicates - offset)
  } else {
    persist_reps <- .sbc_load_persist_reps(persist_dir, persist_base)
    if (length(persist_reps) > 0) {
      env <- new.env(parent = emptyenv())
      load(fileName, envir = env)
      if (!"prior_alpha" %in% ls(env))
        stop("Cannot restart: prior_alpha not found in ", fileName)
      prior_alpha <- env$prior_alpha
      dir.create(temp_dir, showWarnings = FALSE)
      saveRDS(prior_alpha, file.path(temp_dir, "prior_samples.rds"))
      for (nm in names(persist_reps)) {
        saveRDS(persist_reps[[nm]], file.path(temp_dir, paste0("rep_", nm, ".rds")))
        completed_results[[nm]] <- persist_reps[[nm]]
      }
      offset <- length(completed_results)
      if (verbose)
        message("Restarting from persistent rep files: ", offset, " of ", replicates,
                " replicates already complete, running remaining ",
                replicates - offset)
    } else {
      prior_alpha <- parameters(prior_in, N = replicates, selection = "alpha")
      if (!is.null(temp_dir)) {
        dir.create(temp_dir, showWarnings = FALSE)
        saveRDS(prior_alpha, file.path(temp_dir, "prior_samples.rds"))
      }
      if (!is.null(fileName)) save(prior_alpha, file = fileName)
    }
  }

  missing_reps <- setdiff(1:replicates, as.integer(names(completed_results)))

  if (length(missing_reps) > 0) {
    if (verbose)
      message("Processing ", dots[["cores_per_chain"]], " data sets in parallel")
    res_new <- auto_mclapply(
      X = missing_reps,
      FUN = run_SBC_subject,
      design_in, prior_alpha, trials, prior_in, dots, temp_dir, offset,
      persist_dir, persist_base,
      mc.cores = dots[["cores_per_chain"]]
    )
    for (idx in seq_along(missing_reps)) {
      i <- missing_reps[[idx]]
      if (is.null(res_new[[idx]])) {
        f_temp <- if (!is.null(temp_dir)) file.path(temp_dir, paste0("rep_", i, ".rds")) else NULL
        f_pers <- if (!is.null(persist_dir) && !is.null(persist_base))
          file.path(persist_dir, paste0(persist_base, "_rep_", i, ".rds")) else NULL
        if (!is.null(f_temp) && file.exists(f_temp))
          res_new[[idx]] <- readRDS(f_temp)
        else if (!is.null(f_pers) && file.exists(f_pers))
          res_new[[idx]] <- readRDS(f_pers)
      }
      if (!is.null(res_new[[idx]]))
        completed_results[[as.character(i)]] <- res_new[[idx]]
    }
  }

  # Assemble in replicate order (shared with recover_sbc)
  SBC <- .sbc_assemble_single(completed_results)
  attr(SBC, "recovered_reps") <- NULL
  if (!is.null(fileName)) {
    save(SBC, prior_alpha, file = fileName)
    unlink(temp_dir, recursive = TRUE)
    .sbc_cleanup_persist_reps(persist_dir, persist_base)
  }
  return(SBC)
}

calc_sbc_stats <- function(stats){
  if(is.null(stats$coverage)) return(NULL)
  out <- list()
  out_names <- names(stats[[1]])
  for(i in 1:length(stats[[1]])){
    out[[out_names[i]]] <- list(
      coverage = colMeans(stats$coverage[[i]]),
      # precision: SD across replicates of the posterior median (lower = more precise)
      precision = apply(stats$med[[i]], 2, sd),
      bias = colMeans(stats$bias[[i]])
    )
  }
  return(out)
}


# Validate the add_to_main argument of the SBC plot functions. add_to_main is
# NULL (no prefix), a single string (the same prefix for every panel), or a
# character vector of length n_panels (one prefix per panel).
.sbc_check_add_to_main <- function(add_to_main, n_panels) {
  if (is.null(add_to_main)) return(NULL)
  if (!is.character(add_to_main))
    stop("add_to_main must be NULL or a character vector")
  if (length(add_to_main) != 1L && length(add_to_main) != n_panels)
    stop("add_to_main must be a single string or a character vector of length ",
         n_panels, " (the number of panels being plotted), not ", length(add_to_main))
  add_to_main
}

# Per-panel prefix for add_to_main (see .sbc_check_add_to_main).
.sbc_main_prefix <- function(add_to_main, panel) {
  if (is.null(add_to_main)) "" else
    if (length(add_to_main) == 1L) add_to_main else add_to_main[panel]
}

# Validate the `main` argument of the SBC plot functions: NULL (use the
# auto-generated per-panel title), a single string (the same title on every
# panel), or a character vector of length n_panels (one title per panel).
.sbc_check_main <- function(main, n_panels) {
  if (is.null(main)) return(NULL)
  if (!is.character(main))
    stop("main must be NULL or a character vector")
  if (length(main) != 1L && length(main) != n_panels)
    stop("main must be a single string or a character vector of length ",
         n_panels, " (the number of panels being plotted), not ", length(main))
  main
}

# Per-panel title from `main` (see .sbc_check_main): NULL when `main` is NULL so
# the caller can fall back to its auto-generated title.
.sbc_main_at <- function(main, panel) {
  if (is.null(main)) NULL else
    if (length(main) == 1L) main else main[panel]
}

# Valid legend() position keywords for add_stats.
.sbc_legend_keywords <- c("bottomright", "bottom", "bottomleft", "left", "topleft",
                          "top", "topright", "right", "center")

# Resolve the add_stats argument of the SBC plot functions into a normalized
# spec, or NULL for none. add_stats is TRUE (all three at defaults), FALSE
# (none), or a named list keyed by a subset of c("coverage","bias","precision").
# Each statistic's entry is either a single legend() position keyword (same for
# every panel), FALSE/NA/NULL to suppress it, or itself a named list/vector keyed
# by parameter giving a per-panel position (parameters not named use the stat's
# default; a per-parameter value of FALSE/NA suppresses that one panel). Omitted
# statistics keep their default position on every panel. The result is a named
# list (over the shown statistics); each element is list(mode="all", pos=) or
# list(mode="per", per=<named chr, NA = suppressed>, default=). Query it per
# panel with .sbc_panel_pos().
.sbc_resolve_stat_spec <- function(add_stats, all_params) {
  defaults <- c(coverage = "topleft", bias = "topright", precision = "bottomleft",
                pvalue = "bottomright")
  if (is.logical(add_stats) && length(add_stats) == 1L) {
    if (isTRUE(add_stats))
      return(lapply(defaults, function(p) list(mode = "all", pos = unname(p))))
    return(NULL)
  }
  if (!is.list(add_stats))
    stop("add_stats must be TRUE/FALSE or a named list")
  bad <- setdiff(names(add_stats), names(defaults))
  if (is.null(names(add_stats)) || length(bad))
    stop("add_stats list names must be a subset of 'coverage', 'bias', 'precision', 'pvalue'; got: ",
         paste(bad, collapse = ", "))
  is_suppress <- function(v) is.null(v) || (length(v) == 1L && (isFALSE(v) || isTRUE(is.na(v))))
  chk_kw <- function(v, ctx) {
    if (!(is.character(v) && length(v) == 1L && v %in% .sbc_legend_keywords))
      stop(ctx, " must be a single legend position keyword (one of: ",
           paste(.sbc_legend_keywords, collapse = ", "), "), or FALSE/NA to suppress")
    v
  }
  spec <- list()
  for (stat in names(defaults)) {
    if (!(stat %in% names(add_stats))) {                       # not mentioned: default everywhere
      spec[[stat]] <- list(mode = "all", pos = unname(defaults[[stat]])); next
    }
    v <- add_stats[[stat]]
    if (is_suppress(v)) next                                   # suppressed entirely
    if (is.list(v) || !is.null(names(v))) {                    # per-parameter
      if (is.null(names(v)) || any(names(v) == ""))
        stop("add_stats$", stat, " per-parameter entries must all be named by parameter")
      pbad <- setdiff(names(v), all_params)
      if (length(pbad))
        stop("add_stats$", stat, " names unknown parameter(s): ", paste(pbad, collapse = ", "))
      per <- stats::setNames(rep(NA_character_, length(v)), names(v))
      for (pn in names(v)) {
        vv <- if (is.list(v)) v[[pn]] else v[[pn]]
        per[pn] <- if (is_suppress(vv)) NA_character_ else chk_kw(vv, paste0("add_stats$", stat, "$", pn))
      }
      spec[[stat]] <- list(mode = "per", per = per, default = unname(defaults[[stat]]))
    } else {                                                   # single keyword for all panels
      spec[[stat]] <- list(mode = "all", pos = chk_kw(v, paste0("add_stats$", stat)))
    }
  }
  if (length(spec) == 0L) NULL else spec
}

# Positions of the statistics to draw on the panel for parameter `param`, given a
# spec from .sbc_resolve_stat_spec(). Returns a named character vector (subset of
# coverage/bias/precision) of legend keywords; suppressed statistics are omitted.
.sbc_panel_pos <- function(spec, param) {
  out <- character(0)
  for (stat in names(spec)) {
    s <- spec[[stat]]
    p <- if (identical(s$mode, "all")) s$pos
         else if (param %in% names(s$per)) s$per[[param]]
         else s$default
    if (!is.na(p)) out[stat] <- p
  }
  out
}

# The three recovery statistics that require single-subject SBC output (and so
# may be unavailable). The calibration p-value ("pvalue") is computed from the
# ranks directly and is always available, so it is excluded here.
.sbc_recovery_stats <- c("coverage", "bias", "precision")

# Format a calibration p-value for display in a panel corner.
.sbc_fmt_p <- function(p) {
  if (length(p) != 1L || is.na(p)) return("NA")
  if (p < 0.001) "<0.001" else formatC(p, format = "f", digits = 3)
}

# Kolmogorov-Smirnov p-value for uniformity of normalized SBC ranks (used by the
# ECDF-difference plot, the curve it matches). Deterministic asymptotic p; the
# ranks are discrete, so this is mildly conservative -- the simultaneous bands on
# the plot carry the exact visual verdict. High p = consistent with calibration.
.sbc_ks_p <- function(r_norm) {
  r_norm <- r_norm[is.finite(r_norm)]
  if (length(r_norm) < 2L) return(NA_real_)
  suppressWarnings(stats::ks.test(r_norm, "punif")$p.value)
}

# Chi-square goodness-of-fit-to-uniform p-value over the displayed histogram bins
# (used by the rank-histogram plot, so it is the exact test of the bars drawn).
# `counts` are the per-bin counts; df = number of bins - 1. High p = calibrated.
.sbc_chisq_p <- function(counts) {
  if (length(counts) < 2L || sum(counts) == 0) return(NA_real_)
  suppressWarnings(stats::chisq.test(counts)$p.value)
}

#' Plot the Histogram of the Observed Rank Statistics of SBC
#'
#' Note that this plot is dependent on the number of bins, and a more general
#' visualization is to use `plot_sbc_ecdf`
#'
#' @param ranks A list of named dataframes of the rank statistic
#' @param bins An integer specifying the number of bins to use when plotting the histogram
#' @param layout Optional. A numeric vector specifying the layout using `par(mfrow = layout)`.
#'   `NA` (the default) auto-computes a grid. `NULL` leaves the current device
#'   layout (`mfrow`/`mfcol`) untouched, drawing the panels into the existing
#'   grid with no new page between them.
#' @param add_stats Logical. `TRUE` (the default) appends a chi-square
#'   goodness-of-fit-to-uniform p-value, computed over the displayed bins (the
#'   exact test of the bars drawn; high p = consistent with calibration), as a
#'   second row of each panel's title. `FALSE` omits it.
#' @param main Optional. `NULL` (the default) uses the auto-generated per-panel
#'   title (the first title row). A single string is used as the title on every
#'   panel. A character vector sets one title per panel and must have the same
#'   length as the number of panels being plotted. The `add_stats` p-value, when
#'   shown, is always added as a second title row beneath this.
#' @param add_to_main Optional. `NULL` (the default) does nothing. A single
#'   string is prepended to the start of the (auto-generated or `main`) title of
#'   every panel. A character vector is prepended panel-by-panel and must have
#'   the same length as the number of panels being plotted (no separator is
#'   inserted).
#' @return No returns
#' @export
plot_sbc_hist <- function(ranks, bins = 10, layout = NA, add_stats = TRUE,
                          main = NULL, add_to_main = NULL){
  if (!is.null(ranks[["rank"]])) ranks <- ranks[["rank"]]
  selects <- names(ranks)
  # add_stats: TRUE (default) appends the chi-square goodness-of-fit-to-uniform
  # p-value (over the displayed bins) as a second title row; FALSE omits it.
  show_p <- !isFALSE(add_stats)
  # layout = NULL means: do not touch the device layout (mfrow/mfcol); draw the
  # panels into whatever grid is already set up, with no new page between them.
  # (layout = NA auto-computes a grid; a numeric vector sets it via par(mfrow).)
  manage_layout <- !is.null(layout)
  if (manage_layout) {
    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar)) # code line i + 1
  }

  n_sample <- nrow(ranks[[1]])
  low <- qbinom(0.025, n_sample, 1/bins)
  mid <- qbinom(0.5, n_sample, 1/bins)
  high <- qbinom(0.975, n_sample, 1/bins)
  par_names <- colnames(ranks[[1]])
  n_panels <- sum(vapply(ranks, ncol, integer(1)))
  add_to_main <- .sbc_check_add_to_main(add_to_main, n_panels)
  main <- .sbc_check_main(main, n_panels)
  panel <- 0L
  for (j in seq_along(ranks)) {
    if (manage_layout) {
      if (any(is.na(layout))) {
        par(mfrow = coda_setmfrow(Nparms = ncol(ranks[[1]])))
      } else {
        par(mfrow = layout)
      }
    }
    rank <- ranks[[j]]
    for (i in 1:ncol(rank)) {
      panel <- panel + 1L
      base_main <- if (is.null(main)) paste0(selects[j], " - ", par_names[i])
                   else .sbc_main_at(main, panel)
      title_main <- paste0(.sbc_main_prefix(add_to_main, panel), base_main)
      h <- hist(rank[ , i], breaks = bins, plot = FALSE)
      # y-limit must cover both the uniform reference band (high) and the tallest
      # bar -- bars can exceed the band when miscalibrated -- plus some headroom.
      ylim_hi <- max(high, h$counts) + 2
      plot(h, main = "", xlab = "rank", ylim = c(0, ylim_hi))
      # Second title row: chi-square calibration p-value over the drawn bins.
      if (show_p)
        title_main <- paste0(title_main, "\n",
                             "p(chisq) : ", .sbc_fmt_p(.sbc_chisq_p(h$counts)))
      title(main = title_main)
      abline(h = low, lty = 2)
      abline(h = mid, lty = 2)
      abline(h = high, lty = 2)
    }
  }
}



get_gamma <- function (N, K, conf_level = 0.95)
{
  p_interior <- function (p_int, x1, x2, z1, z2, gamma, N)
  {
    z_tilde <- (z2 - z1)/(1 - z1)
    N_tilde <- rep(N - x1, each = length(x2))
    p_int <- rep(p_int, each = length(x2))
    x_diff <- outer(x2, x1, "-")
    p_x2_int <- p_int * dbinom(x_diff, N_tilde, z_tilde)
    list(p_int = rowSums(p_x2_int), x1 = x2)
  }
  target <- function(gamma, conf_level, N, K) {
    z <- 1:(K - 1)/K
    z1 <- c(0, z)
    z2 <- c(z, 1)
    x2_lower <- qbinom(gamma/2, N, z2)
    x2_upper <- c(N - rev(x2_lower)[2:K], N)
    x1 <- 0
    p_int <- 1
    for (i in seq_along(z1)) {
      tmp <- p_interior(p_int, x1 = x1, x2 = x2_lower[i]:x2_upper[i],
                        z1 = z1[i], z2 = z2[i], gamma = gamma, N = N)
      x1 <- tmp$x1
      p_int <- tmp$p_int
    }
    abs(conf_level - sum(p_int))
  }
  optimize(target, c(0, 1 - conf_level), conf_level, N = N,
           K = K)$minimum
}

get_lims <- function (N, K, gamma)
{
  lims <- list()
  z <- seq(0, 1, length.out = K)
  lims$lower <- qbinom(gamma/2, N, z)/N - z
  lims$upper <- qbinom(1 - gamma/2, N, z)/N - z
  lims$z <- z
  lims
}

make_smooth <- function(x, y, N = 1000){
  lo <- smooth.spline(x, y, spar=0.5)
  xl <- seq(0, 1, 1/N)
  return(predict(lo,xl)$y)
}

#' Plot the ECDF Difference in SBC Ranks
#'
#' Plots  F_hat(z) - z  where F_hat is the
#' empirical CDF of the (normalized) ranks and z is the Uniform(0,1) CDF.
#' The shaded band is the simultaneous (1 - gamma) envelope for F_hat(z) - z.
#'
#' @param ranks A list of named dataframes of the rank statistic (raw or normalized)
#' @param layout Optional. A numeric vector specifying the layout using `par(mfrow = layout)`.
#'   `NA` (the default) auto-computes a grid. `NULL` leaves the current device
#'   layout (`mfrow`/`mfcol`) untouched, drawing the panels into the existing
#'   grid with no new page between them.
#' @param add_stats Controls the statistics overlaid on each panel: the three
#'   recovery statistics (coverage, bias, precision) and `pvalue`, a
#'   Kolmogorov-Smirnov p-value for uniformity of the ranks (the test matching
#'   the ECDF-difference curve; high p = consistent with calibration). `TRUE`
#'   (the default) shows all four at default legend positions (coverage
#'   `"topleft"`, bias `"topright"`, precision `"bottomleft"`, pvalue
#'   `"bottomright"`); `FALSE` shows none. Alternatively a named list keyed by a
#'   subset of `"coverage"`, `"bias"`, `"precision"`, `"pvalue"` may be given:
#'   each value is a `legend()` position keyword (e.g. `"topright"`,
#'   `"bottomleft"`, `"center"`) placing that statistic, an entry of `FALSE` or
#'   `NA` suppresses it, and omitted names keep their default position. A
#'   statistic's value may instead be a named list/vector keyed by parameter (a
#'   subset; parameters not named use that statistic's default position) to set
#'   its position per panel, with `FALSE`/`NA` suppressing it on a given panel.
#'   The three recovery statistics are only available for single-subject SBC
#'   (omitted, with a warning, for hierarchical SBC output); `pvalue` is computed
#'   from the ranks and is always shown.
#' @param main Optional. `NULL` (the default) uses the auto-generated per-panel
#'   title. A single string is used as the title on every panel. A character
#'   vector sets one title per panel and must have the same length as the number
#'   of panels being plotted.
#' @param add_to_main Optional. `NULL` (the default) does nothing. A single
#'   string is prepended to the start of the (auto-generated or `main`) title of
#'   every panel. A character vector is prepended panel-by-panel and must have
#'   the same length as the number of panels being plotted (no separator is
#'   inserted).
#' @param K Optional. Effective sample size of the MCMC that produced the ranks.
#' @return No returns
#' @export
plot_sbc_ecdf <- function(ranks, layout = NA, add_stats = TRUE, main = NULL,
                          add_to_main = NULL, K = 500) {

  # stats extraction (keep existing behaviour)
  stats <- NULL
  if (!is.null(ranks[["rank"]])) {
    stats <- calc_sbc_stats(ranks)
    ranks <- ranks[["rank"]]
  }

  selects <- names(ranks)
  n_panels <- sum(vapply(ranks, ncol, integer(1)))
  add_to_main <- .sbc_check_add_to_main(add_to_main, n_panels)
  main <- .sbc_check_main(main, n_panels)
  panel <- 0L
  stat_spec <- .sbc_resolve_stat_spec(add_stats, unique(unlist(lapply(ranks, colnames))))
  if (is.null(stats) && !is.null(stat_spec) && any(names(stat_spec) %in% .sbc_recovery_stats))
    warning("add_stats was requested but no recovery statistics are available; ",
            "coverage/bias/precision are only produced for single-subject SBC. Skipping ",
            "them (the p(KS) calibration value is still shown).")

  # layout = NULL means: do not touch the device layout (mfrow/mfcol); draw the
  # panels into whatever grid is already set up, with no new page between them.
  # (layout = NA auto-computes a grid; a numeric vector sets it via par(mfrow).)
  manage_layout <- !is.null(layout)
  if (manage_layout) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  # N = number of SBC simulations (rows)
  N <- nrow(ranks[[1]])

  # Envelope for F_hat(z) - z
  gamma <- get_gamma(N, K)
  res <- get_lims(N, K, gamma)  # res$z, res$lower, res$upper

  for (j in seq_along(ranks)) {

    if (manage_layout) {
      if (any(is.na(layout))) {
        par(mfrow = coda_setmfrow(Nparms = ncol(ranks[[j]]), nplots = length(ranks)))
      } else {
        par(mfrow = layout)
      }
    }

    rank <- ranks[[j]]
    stat <- if (!is.null(stats)) stats[[j]] else NULL
    par_names <- colnames(rank)

    for (i in seq_len(ncol(rank))) {

      r <- rank[, i]
      r <- r[is.finite(r)]

      # Convert to normalized ranks in [0,1] if needed
      r_norm <- r
      r_min <- suppressWarnings(min(r_norm, na.rm = TRUE))
      r_max <- suppressWarnings(max(r_norm, na.rm = TRUE))

      if (!(r_max <= 1 && r_min >= 0)) {
        # If looks 1-based integer ranks: 1..K
        if (r_min >= 1 && r_max <= K) {
          r_norm <- (r_norm - 1) / K
        } else if (r_min >= 0 && r_max <= (K - 1)) {
          # If looks 0-based: 0..(K-1)
          r_norm <- r_norm / K
        } else {
          # Fallback: scale by K (keeps it in a reasonable range if ranks are 0..ESS)
          r_norm <- r_norm / K
        }
      }

      # ECDF difference: F_hat(z) - z
      Fn <- ecdf(r_norm)
      z <- res[["z"]]
      y <- Fn(z) - z

      panel <- panel + 1L
      base_main <- if (is.null(main)) paste0(selects[j], " - ", par_names[i])
                   else .sbc_main_at(main, panel)
      main_i <- paste0(.sbc_main_prefix(add_to_main, panel), base_main)

      ylim_lo <- min(res[["lower"]], y, na.rm = TRUE) - 0.01
      ylim_hi <- max(res[["upper"]], y, na.rm = TRUE) + 0.03

      plot(
        x = z,
        y = y,
        type = "s",
        lwd = 2,
        ylim = c(ylim_lo, ylim_hi),
        xlim = c(0, 1),
        ylab = "ECDF Difference  F(z) - z",
        xlab = "Normalized Rank Statistic (z)",
        main = main_i
      )

      polygon(
        x = c(z, rev(z)),
        y = c(res[["lower"]], rev(res[["upper"]])),
        col = adjustcolor("cornflowerblue", 0.2),
        border = NA
      )

      # redraw ECDF diff on top of ribbon
      lines(z, y, type = "s", lwd = 2)

      abline(h = 0, lty = 2, col = "gray50")

      if (!is.null(stat_spec)) {
        pos <- .sbc_panel_pos(stat_spec, par_names[i])
        prints <- character(0)
        if (!is.null(stat)) {
          prints["coverage"]  <- paste0("coverage : ",  round(stat[["coverage"]][i], 2))
          prints["bias"]      <- paste0("bias : ",       round(stat[["bias"]][i], 3))
          prints["precision"] <- paste0("precision : ",  round(stat[["precision"]][i], 3))
        }
        if ("pvalue" %in% names(pos))
          prints["pvalue"] <- paste0("p(KS) : ", .sbc_fmt_p(.sbc_ks_p(r_norm)))
        for (s in names(pos))
          if (s %in% names(prints))
            legend(x = pos[[s]], legend = prints[[s]], bty = "n")
      }
    }
  }

  invisible(NULL)
}

get_ranks_ESS <- function(posterior, ESS, prior){
  posterior <- posterior[seq(1, length(posterior), length.out = ESS)]
  return(pmin(rank(c(prior, posterior))[1]/(ESS+1),1))
}
