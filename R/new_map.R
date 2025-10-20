minimal_design <- function(design, covariates = NULL, drop_subjects = TRUE,
                           n_trials = 1, add_acc = TRUE, drop_R = TRUE,
                           drop_R_levels = TRUE, do_functions = TRUE, ...) {
  dots <- add_defaults(list(...), verbose = TRUE)
  if(!is.null(design$Ffactors)) design <- list(design)
  out <- list()
  for(i in 1:length(design)){
    cur_des <- design[[i]]
    ## 1.  Check and construct the factorial backbone
    if (is.null(cur_des$Ffactors) || !is.list(cur_des$Ffactors)) {
      stop("`cur_des` must contain a list element called `Ffactors`.")
    }
    # Add Rlevels
    cur_des$Ffactors$R <- cur_des$Rlevels
    # Remove subjects
    if(drop_subjects){
      cur_des$Ffactors <- cur_des$Ffactors[names(cur_des$Ffactors) != "subjects"]
    }
    if(n_trials > 1){
      cur_des$Ffactors <- c("trials" = list(1:n_trials), cur_des$Ffactors)
    }
    fac_df <- do.call(expand.grid, c(cur_des$Ffactors))
    n <- nrow(fac_df)
    fac_df[] <- lapply(fac_df, as.factor)
    if(n_trials > 1 & !drop_subjects){
      fac_df$trials <- as.numeric(ave(as.character(fac_df$subjects), as.character(fac_df$subjects), FUN = seq_along))
    } else if(n_trials > 1){
      fac_df$trials <- 1:nrow(fac_df)
    }

    ## 2.  Add covariates (if requested)
    if (!is.null(cur_des$Fcovariates)) {
      for (cv in cur_des$Fcovariates) {
        if(!is.null(list(...)$emc)){
          dat <- list(...)$emc
          if(is.data.frame(dat)){
            covariates[[cv]] <- dat[,cv]
          } else if(!(is(dat, "emc") & is.null(dat[[1]]$subjects))){
            dat <- get_data(dat)
            covariates[[cv]] <- dat[,cv]
          }
        }
        if(!cv %in% names(covariates)){
          if(dots$verbose) message(paste0("Imputing ", cv, " with random values"))
          vec <- rnorm(n)
        } else{
          vec <- covariates[[cv]]
          if(length(vec) < nrow(fac_df)){
            vec <- rep(vec, ceiling(nrow(fac_df)/length(vec)))
          }
        }
        fac_df[[cv]] <- vec[seq_len(n)]
      }
    }

    ## 3. Add accumulators
    if(add_acc){
      fac_df <- add_accumulators(fac_df, matchfun = cur_des$matchfun, type = cur_des$model()$type)
    }
    if(!is.null(fac_df$R) & drop_R_levels){
      fac_df <- fac_df[fac_df$R == unique(fac_df$R)[1],]
    }
    if(!is.null(fac_df$R) & drop_R){
      fac_df <- fac_df[,!colnames(fac_df) %in% c("R", "winner")]
    }

    ## 4.  Derive additional columns via Ffunctions (if any)
    if (!is.null(cur_des$Ffunctions) & do_functions) {
      funs <- cur_des$Ffunctions
      if (!is.list(funs)) funs <- list(funs)

      for (j in seq_along(funs)) {
        res <- funs[[j]](fac_df)
        nm <- names(funs)[j]
        fac_df[[nm]] <- res
      }
    }

    out[[i]] <- fac_df
  }
  if(length(out) == 1) out <- out[[1]]
  return(out)
}

do_map <- function(draws, design,
                   n_trials = 30, data = NULL, functions = NULL,
                   add_recalculated = FALSE,...){
  par_idx <- 0
  all_out <- list()
  for(i in 1:length(design)){
    cur_des <- design[[i]]
    # See which design we have to get out
    cur_idx <- par_idx + seq_len(length(sampled_pars(cur_des)))
    par_idx <- max(cur_idx)
    cur_draws <- draws[,cur_idx]
    if(grepl("MRI", cur_des$model()$type)){
      # MRI mapping is rudimentary for now
      all_out[[i]] <- map_MRI(cur_draws)
      next
    }
    joint <- FALSE
    if(any(grepl("|", colnames(cur_draws), fixed =  T))){
      # Get out joint
      joint <- TRUE
      prefix <- unique(gsub("[|].*", "", colnames(cur_draws)))
      colnames(cur_draws) <- sub("^[^|]*\\|", "", colnames(cur_draws))
    }

    combined <- par_data_map(cur_draws, cur_des, n_trials = n_trials,
                             data = data, functions = functions,
                             add_recalculated = add_recalculated, ...)
    # Now I have the data, and the parameters mapped to them.
    # Next we loop over the columns of the parameters and calculate their mappings
    # An open question is how to deal with added recalculated pars.
    # For the other parameters I can simply identify if they come up in a formula
    # and if so include that column in the mapping (otherwise don't)
    # That's less trivial for recalculated parameters
    if(joint){
      colnames(out) <- paste0(prefix, "|", colnames(out))
    }
    all_out[[i]] <- out
  }
}



par_data_map <- function(par_mcmc, design, n_trials = NULL, data = NULL,
                         functions = NULL, add_recalculated = FALSE, ...){
  design <- design[[1]]
  model <- design$model
  if ( is.null(data) ) {
    design$Ffactors$subjects <- rownames(parameters)
    if ( is.null(n_trials) )
      stop("If data is not provided need to specify number of trials")
    design$Fcovariates <- design$Fcovariates[!design$Fcovariates %in% names(functions)]
    data <- minimal_design(design, covariates = list(...)$covariates,
                           drop_subjects = F, n_trials = n_trials, add_acc=F,
                           drop_R = F)
  } else {
    data <- add_trials(data[order(data$subjects),])
  }
  if(!is.null(functions)){
    for(i in 1:length(functions)){
      data[[names(functions)[i]]] <- functions[[i]](data)
    }
  }
  if (!is.factor(data$subjects)) data$subjects <- factor(data$subjects)
  if (!is.null(model)) {
    if (!is.function(model)) stop("model argument must  be a function")
    if ( is.null(model()$p_types) ) stop("model()$p_types must be specified")
    if ( is.null(model()$Ttransform) ) stop("model()$Ttransform must be specified")
  }
  model <- design$model
  data <- design_model(
    add_accumulators(data,design$matchfun,simulate=TRUE,type=model()$type,Fcovariates=design$Fcovariates),
    design,model,add_acc=FALSE,compress=FALSE,verbose=FALSE,
    rt_check=FALSE)

  n_mcmc <- dim(par_mcmc)[3]
  n_pars <- length(model()$p_types)
  n_subs <- ncol(par_mcmc)
  for(i in 1:n_mcmc){
    parameters <- t(as.matrix(par_mcmc[,,i], nrow = n_pars, ncol = n_subs))
    pars <- t(apply(parameters, 1, do_pre_transform, model()$pre_transform))
    pars <- map_p(add_constants(pars,design$constants),data, model())
    if(!is.null(model()$trend) && attr(model()$trend, "pretransform")){
      # This runs the trend and afterwards removes the trend parameters
      pars <- prep_trend(data, model()$trend, pars)
    }
    pars <- do_transform(pars, model()$transform)
    if(!is.null(model()$trend) && attr(model()$trend, "posttransform")){
      # This runs the trend and afterwards removes the trend parameters
      pars <- prep_trend(data, model()$trend, pars)
    }
    if(add_recalculated) pars <- model()$Ttransform(pars, data)
    if(i == 1){
      # Ttransform could add unwanted friends, so safest to just
      # figure out the dimensions here
      out <- array(NA, dim = c(nrow(pars), n_mcmc, ncol(pars)),
                   dimnames = list(NULL, NULL, colnames(pars)))
    }
    out[,i,] <- pars
  }
  return(list(data = data, pars = out))
}


do_map_old <- function(draws, design, add_recalculated = FALSE, ...) {
  if(!is.matrix(draws)){
    draws <- draws[,1,]
    is_array <- TRUE
  } else{
    is_array <- FALSE
  }
  draws <- t(draws)
  n_draws <- nrow(draws)

  if(!is.null(design$Ffactors)) design <- list(design)
  all_out <- list()
  par_idx <- 0
  for(i in 1:length(design)){
    cur_des <- design[[i]]

    cur_idx <- par_idx + seq_len(length(sampled_pars(cur_des)))
    par_idx <- max(cur_idx)
    cur_draws <- draws[,cur_idx]
    if(grepl("MRI", cur_des$model()$type)){
      all_out[[i]] <- map_MRI(cur_draws)
      next
    }

    if(any(grepl("|", colnames(cur_draws), fixed =  T))){
      joint <- T
      prefix <- unique(gsub("[|].*", "", colnames(cur_draws)))
      colnames(cur_draws) <- sub("^[^|]*\\|", "", colnames(cur_draws))
    } else joint <- FALSE
    ## 0.  coefficient draw matrix & constants ------------------------------
    if (!is.null(cur_des$constants)) {
      cur_draws <- cur_draws[, !colnames(cur_draws) %in% names(cur_des$constants), drop = FALSE]
    }
    cur_draws <- t(apply(cur_draws, 1, do_pre_transform, cur_des$model()$pre_transform))
    if(!is.null(cur_des$constants)){
      constants <- matrix(rep(cur_des$constants, each = nrow(cur_draws)), nrow = nrow(cur_draws), dimnames = list(NULL, names(cur_des$constants)))
      cur_draws <- cbind(cur_draws, constants)
    }

    ## 1.  formula list ------------------------------------------------------
    fmls <- if (!is.null(cur_des$Fformula)) cur_des$Fformula else cur_des$Flist
    if (is.null(fmls))
      stop("`design` must contain `Fformula` or `Flist` with parameter formulas.")
    if (is.null(names(fmls)) || any(!nzchar(names(fmls))))
      names(fmls) <- vapply(fmls, function(f) deparse(f[[2L]]), character(1))

    ## 2.  minimal experimental design --------------------------------------
    design_df <- minimal_design(cur_des,
                                drop_subjects = TRUE,
                                n_trials      = 100, # sample 100 here so that covariate is semi-accurately represented
                                add_acc       = TRUE,
                                ...)

    factor_cols    <- unique(c(names(cur_des$Ffactors), names(cur_des$Ffunctions), 'lM', 'lR'))
    covariate_cols <- if (!is.null(cur_des$Fcovariates)) cur_des$Fcovariates else character(0)

    base_cells <- unique(design_df[, colnames(design_df) %in% factor_cols, drop = FALSE])

    ## covariate source rows for sampling
    cov_source <- if (length(covariate_cols))
      design_df[, covariate_cols, drop = FALSE] else NULL
    n_cov_rows <- if (is.null(cov_source)) 1L else nrow(cov_source)

    ## 3.  loop over parameters ---------------------------------------------
    mapped_list <- vector("list", length(fmls))
    names(mapped_list) <- names(fmls)

    for (par in names(fmls)) {
      ## covariates that this formula actually uses -------------------------
      par_covs <- intersect(covariate_cols, all.vars(fmls[[par]]))

      ## model‑matrix template (1 dummy row just to grab column names) -------
      X <- make_full_dm(fmls[[par]], cur_des$Clist, design_df)
      if (!ncol(X)) next                                  # nothing to map
      beta <- cur_draws[, colnames(X), drop = FALSE]

      ## -- build factor‑only matrix to detect unique cells -----------------
      ##  -- build factor-only matrix & find unique cells ------------------------
      zero_df <- base_cells

      if (length(par_covs)) {
        # add zero-filled placeholders for any covariates that appear in the formula
        zero_df[par_covs] <- 0
      }

      # model matrix that ignores covariates by construction
      mm_factor <- make_full_dm(fmls[[par]], cur_des$Clist, zero_df)

      dup  <- duplicated(mm_factor)   # TRUE for repeated rows
      u_df <- zero_df[!dup, , drop = FALSE]

      mm_order <- colnames(mm_factor) # == colnames(X)
      n_cell   <- nrow(u_df)

      ## allocate result matrix --------------------------------------------
      map_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_cell)

      ## --- mapping --------------------------------------------------------
      if (length(par_covs)) {
        ## covariates ARE in this formula → sample per draw
        for (d in seq_len(n_draws)) {
          z <- cov_source[sample.int(n_cov_rows, 1L), par_covs, drop = FALSE]

          this_df <- u_df
          this_df[, par_covs] <- z[rep(1L, n_cell), , drop = FALSE]  # overwrite 0s

          mm_full <- make_full_dm(fmls[[par]], cur_des$Clist, this_df)
          colnames(mm_full) <- colnames(X)
          map_mat[d, ] <- beta[d, ] %*% t(mm_full[, mm_order, drop = FALSE])
        }
      } else {
        ## no covariates in formula → vectorised multiply
        mm_full <- make_full_dm(fmls[[par]], cur_des$Clist, u_df)
        colnames(mm_full) <- colnames(X)
        map_mat[,] <- beta %*% t(mm_full[, mm_order, drop = FALSE])
      }

      ## label columns ------------------------------------------------------
      vary_facs <- intersect(all.vars(fmls[[par]][[3]]), factor_cols)

      ## (1) factor part of the label
      if (length(vary_facs)) {
        fac_lbl <- apply(u_df[, vary_facs, drop = FALSE], 1,
                         function(z) paste(paste0(names(z), z), collapse = "_"))
      } else {
        fac_lbl <- rep("", n_cell)
      }

      ## (2) covariate part of the label (same for every cell if present)
      cov_lbl <- if (length(par_covs)) paste(par_covs, collapse = "_") else ""

      ## (3) combine: par[_fac][_cov]
      colnames(map_mat) <- trimws(paste(par,
                                        ifelse(fac_lbl == "", "", paste0("_", fac_lbl)),
                                        if (cov_lbl != "") paste0("_", cov_lbl) else "",
                                        sep = ""), which = "right")
      tmp_mat <- map_mat
      colnames(tmp_mat) <- rep(par, ncol(tmp_mat))
      map_mat[] <- do_transform(tmp_mat, cur_des$model()$transform)
      mapped_list[[par]] <- map_mat
    }

    ## 4.  combine & drop duplicates ----------------------------------------
    out <- do.call(cbind, mapped_list)
    if (ncol(out) > 1){
      out <- out[, !duplicated(colnames(out)), drop = FALSE]
    }
    if(add_recalculated){
      map_names <- colnames(out)
      tmp <- out
      colnames(tmp) <- get_p_types(map_names)
      tmp <- add_constants(tmp, cur_des$constants)
      tmp[,colnames(tmp) %in% names(cur_des$constants)] <-
        do_transform(tmp[,colnames(tmp) %in% names(cur_des$constants), drop = F], cur_des$model()$transform)
      extra <- add_recalculated_pars(tmp, cur_des$model, map_names)
      out <- cbind(out, extra)
    }
    if(joint){
      colnames(out) <- paste0(prefix, "|", colnames(out))
    }
    all_out[[i]] <- out
  }
  out <- t(do.call(cbind, all_out))
  if(is_array){
    out <- array(out, dim = c(nrow(out), 1, ncol(out)), dimnames = list(rownames(out), "alpha", colnames(out)))
  }
  return(out)
}
