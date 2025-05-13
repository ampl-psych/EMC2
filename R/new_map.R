minimal_design <- function(design, covariates = NULL, drop_subjects = TRUE,
                           n_trials = 1, add_acc = TRUE, ...) {
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
      cur_des$Ffactors$trials <- 1:n_trials
    }
    fac_df <- do.call(expand.grid, c(cur_des$Ffactors, stringsAsFactors = FALSE))
    if(n_trials > 1 & !drop_subjects){
      fac_df$trials <- ave(fac_df$subjects, fac_df$subjects, FUN = seq_along)
    }
    n <- nrow(fac_df)
    fac_df[] <- lapply(fac_df, as.factor)

    ## 2.  Add covariates (if requested)
    if (!is.null(cur_des$Fcovariates)) {
      for (cv in cur_des$Fcovariates) {
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

    ## 3.  Derive additional columns via Ffunctions (if any)
    if (!is.null(cur_des$Ffunctions)) {
      funs <- cur_des$Ffunctions
      if (!is.list(funs)) funs <- list(funs)

      for (i in seq_along(funs)) {
        res <- funs[[i]](fac_df)
        nm <- names(funs)[i]
        fac_df[[nm]] <- res
      }
    }

    ## 4. Add accumulators
    if(add_acc){
      fac_df <- add_accumulators(fac_df, matchfun = cur_des$matchfun, type = cur_des$model()$type)
    }
    if(!is.null(fac_df$R)){
      fac_df <- fac_df[fac_df$R == unique(fac_df$R)[1],]
      fac_df <- fac_df[,!colnames(fac_df) %in% c("R", "winner")]
    }
    out[[i]] <- fac_df
  }
  if(length(out) == 1) out <- out[[1]]
  return(out)
}

add_map <- function(draws, design, add_recalculated = FALSE, ...) {

  ## 0.  coefficient draw matrix & constants ------------------------------
  draws <- as.matrix(draws)
  n_draws <- nrow(draws)
  if (!is.null(design$constants)) {
    draws <- draws[, !colnames(draws) %in% names(design$constants), drop = FALSE]
  }
  draws <- t(apply(draws, 1, do_pre_transform, design$model()$pre_transform))
  ## 1.  formula list ------------------------------------------------------
  fmls <- if (!is.null(design$Fformula)) design$Fformula else design$Flist
  if (is.null(fmls))
    stop("`design` must contain `Fformula` or `Flist` with parameter formulas.")
  if (is.null(names(fmls)) || any(!nzchar(names(fmls))))
    names(fmls) <- vapply(fmls, function(f) deparse(f[[2L]]), character(1))

  ## 2.  minimal experimental design --------------------------------------
  design_df <- minimal_design(design,
                              drop_subjects = TRUE,
                              n_trials      = 100,
                              add_acc       = FALSE,
                              ...)

  factor_cols    <- unique(c(names(design$Ffactors), names(design$Ffunctions)))
  covariate_cols <- if (!is.null(design$Fcovariates)) design$Fcovariates else character(0)

  base_cells <- unique(design_df[, colnames(design_df) %in% factor_cols, drop = FALSE])

  ## covariate source rows for sampling
  cov_source <- if (length(covariate_cols))
    design_df[, covariate_cols, drop = FALSE] else NULL
  n_cov_rows <- if (is.null(cov_source)) 1L else nrow(cov_source)

  ## 3.  loop over parameters ---------------------------------------------
  mapped_list <- vector("list", length(fmls))
  names(mapped_list) <- names(fmls)

  for (par in names(fmls)) {

    rhs_only <- stats::update(fmls[[par]], NULL ~ .)

    ## covariates that this formula actually uses -------------------------
    par_covs <- intersect(covariate_cols, all.vars(rhs_only))

    ## model‑matrix template (1 dummy row just to grab column names) -------
    X <- model.matrix(rhs_only,
                      data           = design_df[1, , drop = FALSE],
                      contrasts.arg  = design$Fcontrasts)

    ## rename coefficient columns so they match the draws -----------------
    if (ncol(X) == 1) {
      colnames(X) <- par
    } else if (attr(stats::terms(fmls[[par]]), "intercept") != 0) {
      colnames(X) <- c(par, paste(par, colnames(X)[-1], sep = "_"))
    } else {
      colnames(X) <- paste(par, colnames(X), sep = "_")
    }
    X <- X[, !colnames(X) %in% names(design$constants), drop = FALSE]
    if (!ncol(X)) next                                  # nothing to map
    beta <- draws[, colnames(X), drop = FALSE]

    ## -- build factor‑only matrix to detect unique cells -----------------
    zero_df <- base_cells
    if (length(par_covs))
      for (cv in par_covs) zero_df[[cv]] <- 0

    mm_factor <- model.matrix(rhs_only, zero_df,
                              contrasts.arg = design$Fcontrasts)
    colnames(mm_factor) <- colnames(X)

    keep_cell     <- !duplicated(t(mm_factor))
    base_cells_u  <- zero_df[keep_cell, , drop = FALSE]
    mm_order      <- colnames(X)
    n_cell        <- nrow(base_cells_u)

    ## allocate result matrix --------------------------------------------
    map_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_cell)

    ## --- mapping --------------------------------------------------------
    if (length(par_covs)) {
      ## covariates ARE in this formula → sample per draw
      for (d in seq_len(n_draws)) {
        z <- cov_source[sample.int(n_cov_rows, 1L), par_covs, drop = FALSE]

        this_df <- base_cells_u
        this_df[, par_covs] <- z[rep(1L, n_cell), , drop = FALSE]  # overwrite 0s

        mm_full <- model.matrix(rhs_only, this_df,
                                contrasts.arg = design$Fcontrasts)
        colnames(mm_full) <- colnames(X)
        map_mat[d, ] <- beta[d, ] %*% t(mm_full[, mm_order, drop = FALSE])
      }
    } else {
      ## no covariates in formula → vectorised multiply
      mm_full <- model.matrix(rhs_only, base_cells_u,
                              contrasts.arg = design$Fcontrasts)
      colnames(mm_full) <- colnames(X)
      map_mat[,] <- beta %*% t(mm_full[, mm_order, drop = FALSE])
    }

    ## label columns ------------------------------------------------------
    vary_facs <- intersect(attr(stats::terms(rhs_only), "term.labels"), factor_cols)

    ## (1) factor part of the label
    if (length(vary_facs)) {
      fac_lbl <- apply(base_cells_u[, vary_facs, drop = FALSE], 1,
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
    map_mat[] <- do_transform(tmp_mat, design$model()$transform)
    mapped_list[[par]] <- map_mat
  }

  ## 4.  combine & drop duplicates ----------------------------------------
  out <- do.call(cbind, mapped_list)
  if (ncol(out) > 1){
    out <- out[, !duplicated(colnames(out)), drop = FALSE]
  }
  if(add_recalculated){
    browser()
    map_names <- colnames(out)
    tmp <- out
    colnames(tmp) <- get_p_types(map_names)
    tmp <- add_constants(tmp, design$constants)
    tmp[,colnames(tmp) %in% names(design$constants)] <-
      do_transform(tmp[,colnames(tmp) %in% names(design$constants), drop = F], design$model()$transform)
    extra <- add_recalculated_pars(tmp, design$model, map_names)
    out <- cbind(out, extra)
  }
  out
}
