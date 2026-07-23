# Apply shared renaming to a character vector of pnames
.apply_shared_pnames <- function(pnames, shared) {
  for (new_name in names(shared)) {
    pnames[pnames %in% shared[[new_name]]] <- new_name
  }
  unique(pnames)
}

# Apply shared renaming to a named character vector of transforms
.apply_shared_transforms <- function(transforms, shared) {
  if (is.null(transforms) || !length(transforms)) return(transforms)
  nms <- names(transforms)
  for (new_name in names(shared)) {
    nms[nms %in% shared[[new_name]]] <- new_name
  }
  transforms <- stats::setNames(transforms, nms)
  # keep first transform per name (should be identical in practice)
  transforms[!duplicated(names(transforms))]
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Create a kernel specification
#'
#' @param cov_names Character vector of covariate column names to apply the
#'   kernel to.
#' @param type Character string specifying the kernel type (e.g. `"delta"`,
#'   `"exp_incr"`). See [trend_help()] for all options.
#' @param par_input Optional character vector of parameter names to use as
#'   additional inputs to the kernel.
#' @param kernel_args Optional named list of kernel-specific arguments.
#'   Currently supported:
#'   \describe{
#'     \item{`q_reset_column`}{For delta-family kernels (`"delta"`,
#'       `"delta2kernel"`, `"delta2lr"`) only. Name of a logical or integer
#'       column in `data` indicating trials on which the Q-value should be
#'       reset to `q0` before the prediction error is computed.
#'       `TRUE`/`1` triggers a reset; `FALSE`/`0` does not.}
#'     \item{`belief_reset_column`}{For DBM-family kernels (`"beta_binomial"`,
#'       `"beta_binomial_decay"`, `"beta_binomial_window"`, `"DBM"`, `"TPM"`)
#'       only. Name of a logical or integer column in `data` indicating trials
#'       on which the predictive distribution should be re-initialised, that is,
#'       all preceding observations should be forgotten.
#'       `TRUE`/`1` triggers a reset; `FALSE`/`0` does not.}#'
#'   }
#' @param custom_kernel A custom kernel registered with [register_kernel()].
#'   Required when `kernel = "custom"`.
#' @param at If `NULL`, the kernel is applied to every row in `dadm`. If a
#'   factor column name (e.g. `"lR"`), the kernel is applied only to entries
#'   corresponding to the first level of that factor and fed forward to the
#'   other levels. Defaults to `"lR"`. For DDMs, `at` should be set to `NULL`.
#'
#' @return An object of class `emc2_kernel`.
#' @seealso [make_base()], [make_trend()], [trend_help()]
#' @export
make_kernel <- function(cov_names,
                        type,
                        par_input          = NULL,
                        kernel_args        = NULL,
                        custom_kernel      = NULL,
                        at                 = "lR") {

  # ---- validate kernel type ----
  known <- names(trend_help(type, return_types = TRUE)$kernels)
  if (!type %in% known)
    stop("Kernel type '", type, "' not recognised. See trend_help().")

  if (identical(type, "custom") && is.null(custom_kernel))
    stop("custom_kernel must be provided when type = 'custom'.")

  # ---- normalise cov_names / par_input ----
  if (is.null(cov_names)) cov_names <- character(0)
  if (is.null(par_input)) par_input <- character(0)

  # ---- validate / normalise kernel_args ----
  delta_kernels <- .sequential_kernels()
  if (!is.null(kernel_args)) {
    if (!is.list(kernel_args))
      stop("kernel_args must be NULL or a named list.")
    if (!is.null(kernel_args$q_reset_column) && !type %in% delta_kernels) {
      warning("kernel_args$q_reset_column is only meaningful for delta-family ",
              "kernels; ignored for kernel '", type, "'.")
      kernel_args$q_reset_column <- NULL
    }
    if (!is.null(kernel_args$belief_reset_column) && !type %in% delta_kernels) {
      warning("kernel_args$belief_reset_column is only meaningful for DBM-family ",
              "kernels; ignored for kernel '", type, "'.")
      kernel_args$belief_reset_column <- NULL
    }
  }

  # ---- generic parameter names and transforms (no prefix yet) ----
  if (identical(type, "custom")) {
    ck <- custom_kernel
    if (!inherits(ck, "emc2_custom_kernel")) stop("custom_kenrel must be created by register_kernel().")
    generic_pnames     <- ck$kernel_pnames
    # custom transforms: named character vector or NULL
    ckf <- ck$transforms
    generic_transforms <- if (!is.null(ckf)) {
      stats::setNames(as.character(ckf), ck$kernel_pnames)
    } else {
      stats::setNames(rep("identity", length(generic_pnames)), generic_pnames)
    }
  } else {
    kinfo              <- trend_help(type, do_return = TRUE, show_experimental=TRUE)
    generic_pnames     <- kinfo$default_pars %||% character(0)
    raw_tf             <- kinfo$transforms$func %||% list()
    generic_transforms <- stats::setNames(
      as.character(unlist(raw_tf)),
      names(raw_tf)
    )
  }

  tmp_id <- paste0("k_tmp_", sample.int(1e9, 1))
  structure(
    list(
      kernel_id          = tmp_id,
      type        = type,
      cov_names          = cov_names,
      par_input          = par_input,
      kernel_args        = kernel_args,
      kernel_pointer     = if (!is.null(custom_kernel)) custom_kernel$kernel_pointer else NULL,
      at                 = at,
      sequential         = type %in% .sequential_kernels(),
      # generic (unprefixed) — finalised to prefixed in make_trend()
      generic_pnames     = generic_pnames,
      generic_transforms = generic_transforms
    ),
    class = "emc2_kernel"
  )
}


#' Create a base specification linking a kernel output to a model parameter
#'
#' @param target_parameter Character string naming the model parameter to which
#'   the trend is applied (e.g. `"v"`, `"B"`).
#' @param type Character string specifying the base function. One of
#'   `"lin"`, `"centered"`, `"add"`, `"identity"`. See [trend_help()] for
#'   details.
#' @param kernel An `emc2_kernel` object created by [make_kernel()].
#' @param kernel_output Integer selecting which output of the kernel to use.
#'   Non-sequential kernels always have one output (`1L`). Sequential kernels
#'   (delta family) may have more; use [trend_help()] to see available outputs.
#'   Defaults to `1L`.
#' @param coding Optional named list of functions that translate trial-level
#'   covariates into accumulator-specific predictors. Each function must accept
#'   `(dadm, cov_names)` and return a matrix of size
#'   `(nrow(dadm), length(cov_names))`, coding how each covariate contributes
#'   to each accumulator on each trial.
#' @param phase Character string; one of `"premap"`, `"pretransform"`,
#'   `"posttransform"`. Controls when the trend is applied in the parameter
#'   computation pipeline. Defaults to `"premap"`.
#'
#' @return An object of class `emc2_base`.
#'
#' @details
#' ## Coding
#'
#' The `coding` argument allows you to specify how trial-by-trial covariates
#' influence model parameters for each accumulator. This is useful when the
#' same covariate enters different accumulators with different signs or weights
#' (e.g. an advantage coding for RL-DDM or RL-ARD models).
#'
#' A coding function must have the signature `function(dadm, cov_names)` and
#' return a numeric matrix with `nrow(dadm)` rows and `length(cov_names)`
#' columns. The matrix codes how each covariate contributes to the current
#' accumulator on each trial.
#'
#' A minimal example — an *advantage* map that codes the chosen accumulator
#' as +1 and the unchosen accumulator as -1:
#'
#' ```r
#' advantage_map <- function(dadm, cov_names) {
#'   lS      <- paste0("cov", ifelse(dadm$lR == "left",  dadm$cov_left,  dadm$cov_right))
#'   lSother <- paste0("cov", ifelse(dadm$lR == "right", dadm$cov_left,  dadm$cov_right))
#'   plus  <- sapply(cov_names, function(x) ifelse(lS      == x,  1, 0))
#'   minus <- sapply(cov_names, function(x) ifelse(lSother == x, -1, 0))
#'   plus + minus
#' }
#' ```
#'
#' @examples
#' # --- Single covariate map (advantage coding) ---
#'
#' data <- data.frame(
#'   subjects  = rep(1, 2),
#'   S         = c("left", "right"),
#'   cov_left  = c("1", "4"),
#'   cov_right = c("3", "2"),
#'   rt        = c(1.2, 0.8),
#'   R         = factor(c("left", "right")),
#'   cov1 = c(1, NA), cov2 = c(NA, 1),
#'   cov3 = c(NA, 1), cov4 = c(1,  1)
#' )
#'
#' advantage_map <- function(dadm, cov_names) {
#'   lS      <- paste0("cov", ifelse(dadm$lR == "left",  dadm$cov_left,  dadm$cov_right))
#'   lSother <- paste0("cov", ifelse(dadm$lR == "right", dadm$cov_left,  dadm$cov_right))
#'   map_plus  <- sapply(cov_names, function(col) ifelse(lS      == col,  1, 0))
#'   map_minus <- sapply(cov_names, function(col) ifelse(lSother == col, -1, 0))
#'   map_plus + map_minus
#' }
#'
#' k <- make_kernel(cov_names=c("cov1", "cov2", "cov3", "cov4"), type="delta", at="lR")
#' b <- make_base(target_parameter="v", type="lin", kernel=k,
#'                coding = list(differences = advantage_map))
#' trend <- make_trend(b)
#'
#' design_RDM <- design(
#'   model   = RDM,
#'   data    = data,
#'   formula = list(B ~ 1, v ~ 1, t0 ~ 1),
#'   trend   = trend
#' )
#'
#' emc <- make_emc(data, design_RDM, type = "single")
#' attr(emc[[1]]$data[[1]], "covariate_coding")
#'
#'
#' # --- Multiple covariate coding (advantage + sum, for RL-ARD models) ---
#'
#' sum_map <- function(dadm, cov_names) {
#'   lS      <- paste0("cov", ifelse(dadm$lR == "left",  dadm$cov_left,  dadm$cov_right))
#'   lSother <- paste0("cov", ifelse(dadm$lR == "right", dadm$cov_left,  dadm$cov_right))
#'   map_this  <- sapply(cov_names, function(col) ifelse(lS      == col, 1, 0))
#'   map_other <- sapply(cov_names, function(col) ifelse(lSother == col, 1, 0))
#'   map_this + map_other
#' }
#'
#' b2 <- make_base(target_parameter = "v", type="lin", kernel = k,
#'                 coding = list(sums = sum_map))
#' trend2 <- make_trend(b, b2)
#'
#' # The model now includes two base parameters (e.g. v.w_differences, v.w_sums)
#' get_trend_pnames(trend2)
#'
#' design_RDM2 <- design(
#'   model   = RDM,
#'   data    = data,
#'   formula = list(B ~ 1, v ~ 1, t0 ~ 1),
#'   trend   = trend2
#' )
#'
#' emc2 <- make_emc(data, design_RDM2, type = "single")
#' attr(emc2[[1]]$data[[1]], "covariate_coding")
#'
#' @seealso [make_kernel()], [make_trend()], [trend_help()]
#' @export
make_base <- function(target_parameter,
                      type,
                      kernel,
                      kernel_output = 1L,
                      coding        = NULL,
                      phase         = "premap") {

  # ---- validate ----
  if (!inherits(kernel, "emc2_kernel"))
    stop("kernel must be an emc2_kernel object created by make_kernel().")

  known_bases <- names(trend_help(return_types = TRUE)$bases)
  if (!type %in% known_bases)
    stop("base type '", type, "' not recognised. See trend_help().")

  if (!phase %in% c("premap", "pretransform", "posttransform"))
    stop("phase must be one of 'premap', 'pretransform', 'posttransform'.")

  kernel_output <- as.integer(kernel_output)

  # ---- validate coding ----
  if (!is.null(coding)) {
    if (!is.list(coding) || !all(vapply(coding, is.function, logical(1))))
      stop("coding must be a named list of functions.")
    if (is.null(names(coding)) || any(names(coding) == ""))
      stop("All entries in coding must be named.")
    if (length(coding) > 1)
      stop("coding must contain exactly one entry. To use multiple coding schemes, ",
           "create one make_base() per scheme and combine them in make_trend().")
  }

  # ---- generic base parameter names and transforms (no prefix yet) ----
  binfo              <- trend_help(base = type, do_return = TRUE, coding = coding)
  generic_pnames     <- binfo$default_pars %||% character(0)
  raw_tf             <- binfo$transforms$func %||% list()
  generic_transforms <- stats::setNames(
    as.character(unlist(raw_tf)),
    names(raw_tf)
  )

  structure(
    list(
      type               = type,
      target_parameter   = target_parameter,
      kernel_id          = kernel$kernel_id,
      kernel_output      = kernel_output,
      coding             = coding,
      phase              = phase,
      # generic (unprefixed) — finalised to prefixed in make_trend()
      generic_pnames     = generic_pnames,
      generic_transforms = generic_transforms,
      # carry kernel object so make_trend() can retrieve it
      .kernel_obj        = kernel
    ),
    class = "emc2_base"
  )
}


#' Create a trend specification
#'
#' Combines one or more [make_base()] specifications (each of which references
#' a [make_kernel()]) into a complete trend object ready for use in [design()].
#'
#' @param ... One or more `emc2_base` objects created by [make_base()].
#' @param shared Optional named list for parameter sharing across bases or
#'   kernels. Each entry maps a new shared name to a character vector of
#'   existing prefixed parameter names that should all resolve to that single
#'   name. Example: `list(alpha = c("v.alpha", "B.alpha"))`. Validation is
#'   performed: all referenced names must exist in the assembled parameter set.
#'
#' @return An object of class `emc2_trend` with components:
#'   \describe{
#'     \item{`kernels`}{Named list of `emc2_kernel` objects keyed by
#'       `kernel_id`. Each has finalised `pnames` (prefixed parameter names)
#'       and `transforms` slots.}
#'     \item{`bases`}{List of `emc2_base` objects, each with finalised
#'       `pnames` and `transforms` slots.}
#'     \item{`sequential`}{Logical; `TRUE` if any kernel is from the delta
#'       family and requires trial-by-trial sequential updating.}
#'   }
#'
#' @details
#' ## Parameter naming
#'
#' Parameter names are assembled in `make_trend()` by combining the
#' `target_parameter` of the first base that references a given kernel with
#' the kernel's generic parameter names (e.g. `q0`, `alpha` → `v.q0`,
#' `v.alpha`). Base parameters are prefixed with their own `target_parameter`
#' (e.g. `v.w`).
#'
#' @examples
#' # Simple trend on v with a delta kernel
#' k <- make_kernel(cov_names = "reward", type="delta", )
#' b <- make_base(target_parameter = "v", type="lin", kernel = k)
#' trend <- make_trend(b)
#' get_trend_pnames(trend)
#'
#' # Reusing one kernel for two parameters
#' k    <- make_kernel(cov_names = "reward", type="delta")
#' b_v  <- make_base(target_parameter = "v", type="lin", kernel = k, kernel_output = 1L)
#' b_B  <- make_base(target_parameter = "B", type="add", kernel = k, kernel_output = 2L)
#' trend <- make_trend(b_v, b_B)
#' get_trend_pnames(trend)
#'
#' @seealso [make_kernel()], [make_base()], [trend_help()]
#' @export
make_trend <- function(..., shared = NULL) {

  bases <- list(...)

  if (length(bases) == 0)
    stop("Provide at least one emc2_base object.")

  not_base <- !vapply(bases, inherits, logical(1), what = "emc2_base")
  if (any(not_base))
    stop("All arguments to make_trend() must be emc2_base objects created by make_base().")

  # ------------------------------------------------------------------
  # Step 1: collect unique kernels by temporary id (order of first appearance).
  # Must happen before any demotion so .kernel_obj is still accessible.
  # ------------------------------------------------------------------
  kernels_list <- list()
  first_target <- list()

  for (b in bases) {
    tmp_kid <- b$kernel_id
    if (is.null(kernels_list[[tmp_kid]])) {
      kernels_list[[tmp_kid]] <- b$.kernel_obj
      first_target[[tmp_kid]] <- b$target_parameter
    } else if (!identical(first_target[[tmp_kid]], b$target_parameter)) {
      warning(
        "A kernel is referenced by bases targeting '",
        first_target[[tmp_kid]], "' and '", b$target_parameter, "'. ",
        "Its parameter names will be prefixed with '", first_target[[tmp_kid]], "'. ",
        "Consider using `shared` to give them a meaningful shared name."
      )
    }
  }

  # ------------------------------------------------------------------
  # Step 2: replace temporary ids with clean local ids (k1, k2, ...)
  # ------------------------------------------------------------------
  tmp_ids   <- names(kernels_list)
  local_ids <- paste0("k", seq_along(kernels_list))
  old_to_new <- stats::setNames(local_ids, tmp_ids)

  names(kernels_list) <- local_ids
  for (kid in local_ids)
    kernels_list[[kid]]$kernel_id <- kid

  # update kernel_id references on bases
  for (i in seq_along(bases))
    bases[[i]]$kernel_id <- old_to_new[[bases[[i]]$kernel_id]]

  # also update first_target keys
  names(first_target) <- local_ids

  # ------------------------------------------------------------------
  # Step 3: finalise kernels — prefix pnames + transforms,
  # ------------------------------------------------------------------
  for (kid in local_ids) {
    k      <- kernels_list[[kid]]
    prefix <- first_target[[kid]]
    pnames <- if (length(k$generic_pnames) > 0)
      paste0(prefix, ".", k$generic_pnames)
    else
      character(0)
    #pnames <- paste0(prefix, ".", k$generic_pnames)

    transforms <- if (length(k$generic_transforms)) {
      stats::setNames(as.character(k$generic_transforms), pnames)
    } else {
      NULL
    }

    # attr(kernels_list[[kid]], "generic_pnames")     <- k$generic_pnames
    # attr(kernels_list[[kid]], "generic_transforms") <- k$generic_transforms
    kernels_list[[kid]]$generic_pnames     <- NULL
    kernels_list[[kid]]$generic_transforms <- NULL

    kernels_list[[kid]]$pnames     <- pnames
    kernels_list[[kid]]$transforms <- transforms
  }

  # ------------------------------------------------------------------
  # Step 3b: detect colliding base pnames before finalising
  # ------------------------------------------------------------------
  all_candidate_pnames <- lapply(bases, function(b) {
    if (length(b$generic_pnames) > 0)
      paste0(b$target_parameter, ".", b$generic_pnames)
    else
      character(0)
  })

  # find pnames that appear in more than one base
  flat <- unlist(all_candidate_pnames, use.names = FALSE)
  colliding <- unique(flat[duplicated(flat)])

  # ------------------------------------------------------------------
  # Step 4: finalise bases — prefix pnames + transforms,
  #         demote generic fields and .kernel_obj to attributes
  # ------------------------------------------------------------------
  seen_base_pnames <- character(0)

  for (i in seq_along(bases)) {
    b      <- bases[[i]]
    prefix <- b$target_parameter
    pnames <- if (length(b$generic_pnames) > 0)
      paste0(prefix, ".", b$generic_pnames)
    else
      character(0)

    pnames <- vapply(pnames, function(p) {
      if (!p %in% colliding) return(p)
      p_typed <- paste0(p, "_", kernels_list[[b$kernel_id]]$type)
      if (!p_typed %in% seen_base_pnames) return(p_typed)
      k <- 2L
      candidate <- paste0(p_typed, "_", k)
      while (candidate %in% seen_base_pnames) {
        k <- k + 1L
        candidate <- paste0(p_typed, "_", k)
      }
      candidate
    }, character(1))

    seen_base_pnames <- c(seen_base_pnames, pnames)

    transforms <- if (length(b$generic_transforms)) {
      stats::setNames(as.character(b$generic_transforms), pnames)
    } else {
      NULL
    }

    bases[[i]]$generic_pnames     <- NULL
    bases[[i]]$generic_transforms <- NULL
    bases[[i]]$.kernel_obj        <- NULL

    bases[[i]]$pnames     <- pnames
    bases[[i]]$transforms <- transforms
  }

  # ------------------------------------------------------------------
  # Step 5: validate and apply shared renaming
  # ------------------------------------------------------------------
  if (!is.null(shared)) {
    all_pnames <- unique(c(
      unlist(lapply(kernels_list, `[[`, "pnames"), use.names = FALSE),
      unlist(lapply(bases,        `[[`, "pnames"), use.names = FALSE)
    ))
    for (new_name in names(shared)) {
      bad <- shared[[new_name]][!shared[[new_name]] %in% all_pnames]
      if (length(bad))
        stop("shared references unknown parameter name(s): ",
             paste(bad, collapse = ", "))
    }
    for (kid in local_ids) {
      kernels_list[[kid]]$pnames     <- .apply_shared_pnames(kernels_list[[kid]]$pnames, shared)
      kernels_list[[kid]]$transforms <- .apply_shared_transforms(kernels_list[[kid]]$transforms, shared)
    }
    for (i in seq_along(bases)) {
      bases[[i]]$pnames     <- .apply_shared_pnames(bases[[i]]$pnames, shared)
      bases[[i]]$transforms <- .apply_shared_transforms(bases[[i]]$transforms, shared)
    }
  }

  ## Add phase info for kernels as well - easier to use downstream in C++
  # Compute earliest phase for each kernel across all referencing bases
  phase_order <- c(premap = 1L, pretransform = 2L, posttransform = 3L)

  for (kid in names(kernels_list)) {
    referencing_phases <- vapply(
      Filter(function(b) b$kernel_id == kid, bases),
      function(b) b$phase,
      character(1)
    )
    earliest <- names(which.min(phase_order[referencing_phases]))
    kernels_list[[kid]]$phase <- earliest
  }

  structure(
    list(kernels = kernels_list, bases = bases),
    class = "emc2_trend"
  )
}




#' Get help information for trend kernels and bases
#'
#' @param kernel Character string specifying the kernel type to get information about
#' @param base Character string specifying the base type to get information about
#' @param show_experimental Boolean, if TRUE, will also show information on experimental kernels. These could produce unexpected results or crash. Use at your own risk.
#' @param ... Additional arguments
#' @return Formatted trend information
#' @export
#' @examples
#' # Get information about exponential increasing kernel
#' trend_help(kernel = "exp_incr")
#'
#' # Get information about linear base
#' trend_help(base = "lin")
#'
#' # Return available kernel and base types
#' trend_help()
trend_help <- function(kernel = NULL, base = NULL, show_experimental=FALSE, ...){
  dots <- add_defaults(list(...), do_return = FALSE, return_types = FALSE)
  bases <- get_bases()
  n_pars <- sapply(lapply(bases, '[[', 'default_pars'),length)
  base_2p <- names(n_pars)[n_pars==1] #names(bases)[1:3]
  base_1p <- names(n_pars)[n_pars==0] #names(bases)[4:5]
  kernels <- get_kernels()
  if(dots$return_types){
    return(list(kernels = kernels, bases = bases))
  }

  # don't print experimental kernels
  if (!show_experimental) {
    kernels <- kernels[!sapply(kernels, function(k) isTRUE(k$experimental))]
  }

  if (is.null(kernel) && is.null(base)) {
    cat("Available kernels:\n")
    for (k in names(kernels)) {
      cat(paste0("  ", k, ": ", kernels[[k]]$description, "\n"))
    }
    cat("\nAvailable base types:\n")
    for (m in names(bases)) {
      cat(paste0("  ", m, ": ", bases[[m]]$description, "\n"))
    }
    cat("\nPhase options:\n")
    cat("  premap: Trend is applied before parameter mapping. This means the trend parameters\n")
    cat("          are mapped first, then used to transform cognitive model parameters before \n")
    cat("          their mapping.\n")
    cat("  pretransform: Trend is applied after parameter mapping but before transformations.\n")
    cat("                Cognitive model parameters are mapped first, then trend is applied, \n")
    cat("                followed by transformations.\n")
    cat("  posttransform: Trend is applied after both mapping and transformations.\n")
    cat("                 Cognitive model parameters are mapped and transformed first, \n")
    cat("                 then trend is applied.\n")
  } else {
    if (!is.null(kernel)) {
      if(dots$do_return) return(kernels[[kernel]])
      if (kernel %in% names(kernels)) {
        cat("Description: \n")
        cat(kernels[[kernel]]$description, "\n \n")
        if(!is.null(kernels[[kernel]]$transforms)){
          cat("Default transformations (in order): \n")
          cat(paste0(kernels[[kernel]]$transforms, "\n \n"))
        }
        cat("Available bases, first is the default: \n")
        cat(paste0(kernels[[kernel]]$bases, collapse = ", "), "\n \n")
      } else {
        cat(paste0("Kernel '", kernel, "' is not recognized.\n"))
        cat("Run trend_help() for all available kernels.\n")
      }
    }
    if (!is.null(base)) {
      if(dots$do_return) {
        if(!is.null(dots$coding) & base != 'add') {
          # Update default parameters
          old_base_pars <- bases[[base]]$default_pars
          suffixes <- names(dots$coding)
          new_base_pars <- paste0(old_base_pars, "_", suffixes)

          bases[[base]]$default_pars <- new_base_pars

          # Update transforms
          transforms <- bases[[base]]$transforms
          new_transforms <- list()

          updated <- FALSE
          for (arg in c("func", "lower", "upper")) {

            # Skip if this transform component doesn't exist
            if (!arg %in% names(transforms)) next

            tr_arg <- transforms[[arg]]
            new_tr_arg <- list()

            # Only keep old parameters that actually exist in the transform list
            valid_pars <- intersect(old_base_pars, names(tr_arg))

            # Replicate each old transform value across all new parameter names
            for (p in valid_pars) {
              for (p_new in new_base_pars) {
                new_tr_arg[[p_new]] <- tr_arg[[p]]
                if(!updated) updated <- TRUE
              }
            }

            new_transforms[[arg]] <- new_tr_arg
          }
          if(updated) bases[[base]]$transforms <- new_transforms
        }
        return(bases[[base]])
      }
      if (base %in% names(bases)) {
        cat("Description: \n")
        cat(bases[[base]]$description, "\n \n")
        if(!is.null(bases[[base]]$transforms)){
          cat("Default transformations: \n")
          cat(paste0(bases[[base]]$transforms, "\n \n"))
        }
      } else {
        cat(paste0("base '", base, "' is not recognized.\n"))
        cat("Run trend_help() for all available base types.\n")
      }
    }
  }
}




#' Check and update formula list for trend parameters
#'
#' @param trend An `emc2_trend` object created by [make_trend()].
#' @param covariates Character vector of covariate column names in the data.
#' @param model A model function, or NULL.
#' @param formula List of formulas, or NULL.
#' @param parameter_design A parameter_design list, or NULL.
#' @return Updated formula list with intercept formulas added for missing trend
#'   parameters.
check_trend <- function(trend, covariates = NULL, model = NULL,
                        formula = NULL, parameter_design = NULL) {

  # ---- non-premap bases must target existing model parameters ----
  if (!is.null(model)) {
    model_pnames <- names(model()$p_types)
    for (b in trend$bases) {
      if (!identical(b$phase, "premap") && !b$target_parameter %in% model_pnames)
        stop("pretransform/posttransform base targets '", b$target_parameter,
             "' which is not a model parameter.")
    }
  }

  if (is.null(covariates))
    stop("must specify covariates when using trend")

  # ---- auto-add intercept formulas for missing trend pnames ----
  trend_pnames <- get_trend_pnames(trend)

  if (!is.null(formula)) {
    pd_targets <- NULL
    if (!is.null(parameter_design)) {
      expand_over <- parameter_design$expand_over
      pd_targets  <- if (!is.null(expand_over)) {
        as.vector(outer(rownames(parameter_design$weights),
                        expand_over, paste, sep = "."))
      } else {
        rownames(parameter_design$weights)
      }
    }

    formula_lhs <- unlist(lapply(formula, function(x) all.vars(x)[1]))
    isin <- trend_pnames %in% formula_lhs
    if (!is.null(pd_targets))
      isin <- isin | (trend_pnames %in% pd_targets)

    if (any(!isin)) {
      formula <- c(formula,
                   lapply(trend_pnames[!isin],
                          function(x) stats::as.formula(paste(x, "~ 1"))))
      message("Intercept formula added for trend_pars: ",
              paste(trend_pnames[!isin], collapse = ", "))
    }
  }

  formula
}


#' Update a model function with trend parameter types and transforms
#'
#' @param trend An `emc2_trend` object created by [make_trend()].
#' @param model A model function.
#' @return An updated model function.
update_model_trend <- function(trend, model) {
  model_list <- model()

  seen_kernel_ids <- character(0)

  for (b in trend$bases) {
    # Add base parameters and transforms
    if (!is.null(b$transforms))
      model_list$transform$func <- c(model_list$transform$func, b$transforms)
    model_list$p_types <- c(
      model_list$p_types,
      stats::setNames(numeric(length(b$pnames)), b$pnames)
    )

    # Add corresponding kernel parameters and transforms (once per kernel)
    kid <- b$kernel_id
    if (!kid %in% seen_kernel_ids) {
      k <- trend$kernels[[kid]]
      if (!is.null(k$transforms))
        model_list$transform$func <- c(model_list$transform$func, k$transforms)
      model_list$p_types <- c(
        model_list$p_types,
        stats::setNames(numeric(length(k$pnames)), k$pnames)
      )
      seen_kernel_ids <- c(seen_kernel_ids, kid)
    }
  }

  # Avoid duplicated names (may be caused by shared pnames)
  model_list$transform$func <- model_list$transform$func[!duplicated(names(model_list$transform$func))]
  model_list$p_types        <- model_list$p_types[!duplicated(names(model_list$p_types))]

  model_list$trend <- trend
  model <- function() model_list
  model
  #   model_list <- model()
  #
  #   for (k in trend$kernels) {
  #     if (!is.null(k$transforms))
  #       model_list$transform$func <- c(model_list$transform$func, k$transforms)
  #     model_list$p_types <- c(
  #       model_list$p_types,
  #       stats::setNames(numeric(length(k$pnames)), k$pnames)
  #     )
  #   }
  #
  #   for (b in trend$bases) {
  #     if (!is.null(b$transforms))
  #       model_list$transform$func <- c(model_list$transform$func, b$transforms)
  #     model_list$p_types <- c(
  #       model_list$p_types,
  #       stats::setNames(numeric(length(b$pnames)), b$pnames)
  #     )
  #   }
  #
  #   model_list$trend <- trend
  #   model <- function() model_list
  #   model
}

##' Register a custom C++ trend kernel
##'
##' Compiles and registers a user-provided C++ function that maps per-trial
##' kernel parameters and inputs to a numeric vector. The C++ function must have
##' signature:
##'   NumericVector f(NumericMatrix trend_pars, NumericMatrix input)
##' and provide an exported pointer creator using EMC2_MAKE_PTR.
##'
##' @param kernel_parameters Character vector of kernel parameter names (in order).
##' @param file Path to the C++ file implementing the custom kernel. The file
##'   should include EMC2/userfun.hpp and define a pointer creator (via
##'   EMC2_MAKE_PTR) that is exported to R.
##' @param transforms Optional named character vector or list mapping each custom
##'   kernel parameter name to a transform name (e.g., "identity", "exp", "pnorm").
##'   Length must match `kernel_parameters`. If unnamed but the correct length,
##'   the order is assumed to match `kernel_parameters`.
##' @return An object to pass to `make_base(custom_kernel=...)`, carrying the
##'   pointer, parameter names, default base, and optional transform mapping.
##' @export
register_kernel <- function(kernel_parameters, file, transforms = NULL){
  if (!is.character(kernel_parameters) || length(kernel_parameters) == 0)
    stop("kernel_parameters must be a non-empty character vector")
  if (!file.exists(file)) stop("C++ file not found: ", file)

  # Normalize transforms to a character vector in the order of kernel_parameters
  trf_vec <- NULL
  if (!is.null(transforms)) {
    if (is.list(transforms)) transforms <- unlist(transforms, recursive = FALSE, use.names = TRUE)
    transforms <- unlist(transforms, use.names = TRUE)
    if (length(transforms) != length(kernel_parameters)) {
      stop("length(transforms) must match length(kernel_parameters)")
    }
    if (is.null(names(transforms))) {
      trf_vec <- as.character(transforms)
      names(trf_vec) <- kernel_parameters
    } else {
      # Reorder to match kernel_parameters
      if (!all(sort(names(transforms)) == sort(kernel_parameters))) {
        stop("names(transforms) must match kernel_parameters if names are supplied")
      }
      trf_vec <- as.character(transforms[kernel_parameters])
    }
  }

  # Compile and load user code
  Rcpp::sourceCpp(file)
  maker <- "EMC2_make_custom_trend_ptr"
  if (!exists(maker, mode = "function", inherits = TRUE)) {
    # Try to derive maker name from macro usage in the file
    lines <- tryCatch(readLines(file), error = function(e) character(0))
    m <- regmatches(lines, regexpr("EMC2_MAKE_PTR\\(([^)]+)\\)", lines))
    if (length(m) == 1) {
      sym <- sub(".*EMC2_MAKE_PTR\\(([^)]+)\\).*", "\\1", m)
      cand <- paste0("EMC2_make_", sym, "_ptr")
      if (exists(cand, mode = "function", inherits = TRUE)) maker <- cand
    }
  }
  if (!exists(maker, mode = "function", inherits = TRUE)) {
    stop("Could not find exported pointer maker function. Expected '", maker,
         "' or a function generated by EMC2_MAKE_PTR(name).")
  }
  ptr <- do.call(maker, list())

  obj <- list(kernel_pnames = as.character(kernel_parameters),
              maker = maker,
              file = normalizePath(file),
              kernel_pointer = ptr,
              transforms = trf_vec)
  class(obj) <- "emc2_custom_kernel"
  obj
}

.sequential_kernels <- function() {
  names(Filter(function(k) isTRUE(k$sequential), get_kernels()))
}

#' Get parameter names from a trend object
#'
#' @param trend An `emc2_trend` object created by [make_trend()].
#' @return A character vector of unique parameter names used in the trend,
#'   covering both kernel parameters (e.g. `v.q0`, `v.alpha`) and base
#'   parameters (e.g. `v.w`). Shared parameter renaming is already baked into
#'   the names at `make_trend()` time, so no post-hoc substitution is needed.
#' @export
#' @examples
#' k <- make_kernel("exp_incr", cov_names = "trial")
#' b <- make_base("lin", target_parameter = "v", kernel = k)
#' trend <- make_trend(b)
#' get_trend_pnames(trend)
get_trend_pnames <- function(trend) {
  kernel_pnames <- unlist(lapply(trend$kernels, `[[`, "pnames"), use.names = FALSE)
  base_pnames   <- unlist(lapply(trend$bases,   `[[`, "pnames"), use.names = FALSE)
  unique(c(kernel_pnames, base_pnames))

  # BELOW CODE IS ONLY USEFUL FOR CHECKING AGAINST OLD UNIT TESTS, WHICH ASSUMED A DIFFERENT PARAMETER ORDER
  # seen_kernel_ids <- character(0)
  # pnames <- character(0)
  #
  # for (b in trend$bases) {
  #   pnames <- c(pnames, b$pnames)
  #
  #   kid <- b$kernel_id
  #   if (!kid %in% seen_kernel_ids) {
  #     pnames <- c(pnames, trend$kernels[[kid]]$pnames)
  #     seen_kernel_ids <- c(seen_kernel_ids, kid)
  #   }
  # }
  #
  # shared_names <- names(trend$shared)
  # non_shared   <- unique(pnames[!pnames %in% shared_names])
  # c(non_shared, shared_names)
}

has_trend_coding <- function(model) {
  if (is.null(model) || is.function(model) || is.null(model$trend)) return(FALSE)
  any(vapply(model$trend$bases, function(x) !is.null(x$coding) && length(x$coding) > 0, logical(1)))
}

has_conditional_covariates <- function(design) {
  # find define covariates that depend on behavior -- these are rt, R, or any of the outputs of the functions provided.
  # they can be lRfiltered or not
  function_output_columns <- names(design$Ffunctions)
  behavioral_covariates <- c('rt', 'R', function_output_columns)

  # find actual covariates, look for a match
  trend <- design$model()$trend
  any(vapply(trend$kernels, function(k) any(k$cov_names %in% behavioral_covariates), logical(1)))
}

get_bases <- function() {
  bases <- list(
    lin = list(description = "Linear base: parameter + w * k",
               transforms = list(func = list("w" = "identity")),
               default_pars = "w"),
    centered = list(description = "Centered mapping: parameter + w*(k - 0.5)",
                    transforms = list(func = list("w" = "identity")),
                    default_pars = "w"),
    add = list(description = "Additive base: parameter + k",
               transforms = NULL,
               default_pars = character(0)),
    identity = list(description = "Identity base: k",
                    transforms = NULL,
                    default_pars = character(0))
  )
  bases
}

get_kernels <- function() {
  bases <- get_bases()
  n_pars <- sapply(lapply(bases, '[[', 'default_pars'),length)
  base_2p <- names(n_pars)[n_pars==1]
  base_1p <- names(n_pars)[n_pars==0]

  kernels <- list(
    custom = list(description = "Custom C++ kernel: provided via register_trend().",
                  transforms = NULL,
                  default_pars = character(0),
                  bases = names(bases),
                  sequential = TRUE,          # User decides
                  n_outputs  = NA_integer_,   # determined by the user's C++ implementation
                  NA_allowed = TRUE
                  ),
    lin_decr = list(description = "Decreasing linear kernel: k = -c",
                    transforms = NULL,
                    default_pars = character(0),
                    bases = base_2p,
                    sequential   = FALSE,
                    n_outputs    = 1L,
                    NA_allowed   = FALSE),
    lin_incr = list(description = "Increasing linear kernel: k = c",
                    transforms = NULL,
                    default_pars = character(0),
                    bases = base_2p,
                    sequential   = FALSE,
                    n_outputs    = 1L,
                    NA_allowed   = FALSE),
    exp_decr = list(description = "Decreasing exponential kernel: k = exp(-d_ed * c)",
                    transforms = list(func =list("d_ed" = "exp")),
                    default_pars = "d_ed",
                    bases = base_2p,
                    sequential   = FALSE,
                    n_outputs    = 1L,
                    NA_allowed   = FALSE),
    exp_incr = list(description = "Increasing exponential kernel: k = 1 - exp(-d_ei * c)",
                    transforms = list(func =list("d_ei" = "exp")),
                    default_pars = "d_ei",
                    bases = base_2p,
                    sequential   = FALSE,
                    n_outputs    = 1L,
                    NA_allowed   = FALSE),
    pow_decr = list(description = "Decreasing power kernel: k = (1 + c)^(-d_pd)",
                    transforms = list(func =list("d_pd" = "exp")),
                    default_pars = "d_pd",
                    bases = base_2p,
                    sequential   = FALSE,
                    n_outputs    = 1L,
                    NA_allowed   = FALSE),
    pow_incr = list(description = "Increasing power kernel: k = 1 - (1 + c)^(-d_pi)",
                    transforms = list(func =list("d_pi" = "exp")),
                    default_pars = "d_pi",
                    bases = base_2p,
                    sequential   = FALSE,
                    n_outputs    = 1L,
                    NA_allowed   = FALSE),
    poly2 = list(description = "Quadratic polynomial: k = d1 * c + d2 * c^2",
                 transforms = list(func = list("d1" = "identity", "d2" = "identity")),
                 default_pars = c("d1", "d2"),
                 bases = base_1p,
                 sequential   = FALSE,
                 n_outputs    = 1L,
                 NA_allowed   = FALSE),
    poly3 = list(description = "Cubic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3",
                 transforms = list(func = list("d1" = "identity", "d2" = "identity", "d3" = "identity")),
                 default_pars = c("d1", "d2", "d3"),
                 bases = base_1p,
                 sequential   = FALSE,
                 n_outputs    = 1L,
                 NA_allowed   = FALSE),
    poly4 = list(description = "Quartic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3 + d4 * c^4",
                 transforms = list(func = list("d1" = "identity", "d2" = "identity", "d3" = "identity", "d4" = "identity")),
                 default_pars = c("d1", "d2", "d3", "d4"),
                 bases = base_1p,
                 sequential   = FALSE,
                 n_outputs    = 1L,
                 NA_allowed   = FALSE),
    delta = list(description = paste(
                 "Standard delta rule kernel: k = q[i].\n",
                 "        Updates q[i] = q[i-1] + alpha * (c[i-1] - q[i-1]).\n",
                 "        Parameters: q0 (initial value), alpha (learning rate)."
                 ),
                 default_pars = c("q0", "alpha"),
                 transforms = list(func = list("q0" = "identity", "alpha" = "pnorm")),
                 bases = base_2p,
                 sequential   = TRUE,
                 n_outputs    = 2L,
                 NA_allowed   = TRUE),
    delta2lr = list(description = paste(
                "Dual learning rate delta rule: k = q[i].\n",
                "         Like the standard delta rule, but with separate\n",
                "         learning rates for positive and negative prediction errors.\n",
                "         Parameters: q0 (initial value), alphaPos (learning rate for positive PEs),\n",
                "         alphaNeg (learning rate for negative PEs)."
               ),
              default_pars = c("q0", "alphaPos", "alphaNeg"),
              transforms = list(func = list("q0" = "identity",
                                            "alphaPos" = "pnorm",
                                            "alphaNeg" = "pnorm")),
              bases = base_2p,
              sequential   = TRUE,
              n_outputs    = 2L,
              NA_allowed   = TRUE),
    delta_decoupled = list(description = paste(
              "Delta rule: k = q[i].\n",
              "         Standard delta rule:\n",
              "         Q_{t+1} = Q_{t} + alpha*(outcome - Q_{t})\n",
              "                 = (1-alpha)*Q_{t} + alpha*outcome\n",
              "         Decoupled delta rule:\n",
              "         Q_{t+1} = Q_{t} + alpha*outcome - lambda*Q_{t}\n",
              "                = (1-lambda)*Q_{t} + alpha*outcome\n",
              "         In the standard rule, lambda = alpha (coupled).\n",
              "         Here, lambda and alpha are free to vary independently."
            ),
    default_pars = c("q0", "alpha", "lambda"),
    transforms = list(func = list("q0" = "identity",
                                  "alpha" = "pnorm",
                                  "lambda" = "pnorm")),
    bases = base_2p,
    sequential=TRUE,
    experimental=TRUE,
    n_outputs = 2L,
    NA_allowed=TRUE),
  rescorlawagner = list(description = paste(
                "Rescorla Wagner delta rule: k = q[i].\n",
                "         Like the standard delta rule, but with compound prediction errors:\n",
                "         PE = (r - sum(Q)) - with sum over all active covariates on a trial.\n",
                "         Parameters: q0 (initial value), alpha (learning rate)"
              ),
              default_pars = c("q0", "alpha"),
              transforms = list(func = list("q0" = "identity",
                                            "alpha" = "pnorm")),
              bases = base_2p,
              sequential   = TRUE,
              n_outputs    = 2L,
              experimental = TRUE,
              NA_allowed=TRUE),
  delta2kernel = list(description = paste(
                "Dual kernel delta rule: k = q[i].\n",
                "         Combines fast and slow learning rates\n",
                "         and switches between them based on dSwitch.\n",
                "         Parameters: q0 (initial value), alphaFast (fast learning rate),\n",
                "         propSlow (alphaSlow = propSlow * alphaFast), dSwitch (switch threshold)."
              ),
              default_pars = c("q0", "alphaFast", "propSlow", "dSwitch"),
              transforms = list(func = list("q0" = "identity", "alphaFast" = "pnorm",
                                            "propSlow" = "pnorm", "dSwitch" = "pnorm")),
              bases = base_2p,
              experimental=TRUE,
              NA_allowed=TRUE),
  beta_binomial = list(
    description = paste(
      "Beta-Binomial learning kernel:\n",
      "        k = predicted probability of the current trial's observation\n",
      "        being X as opposed to Y.\n",
      "        Assumes a binary (Bernoulli) sequence of observations.\n",
      "        Parameters: a0, b0 (shape parameters of the Beta prior)."
    ),
    default_pars = c("a0", "b0"),
    transforms   = list(func = list("a0" = "exp", "b0" = "exp")),
    bases        = base_2p,
    sequential   = TRUE,
    n_outputs    = 4L,
    experimental = TRUE,
    NA_allowed=TRUE),
  beta_binomial_decay = list(
    description = paste(
      "Beta-Binomial learning kernel with leaky integration:\n",
      "        k = predicted probability of the current trial's observation\n",
      "        being X as opposed to Y.\n",
      "        Assumes a binary (Bernoulli) sequence of observations.\n",
      "        Applies exponential decay ('leaky integration') to the count\n",
      "        of past observations by multiplying past accumulated counts by exp(-1/decay).\n",
      "        Parameters: a0, b0 (shape parameters of the Beta prior), decay."
    ),
    default_pars = c("a0", "b0", "decay"),
    transforms   = list(func = list("a0" = "exp", "b0" = "exp", "decay" = "exp")),
    bases        = base_2p,
    sequential   = TRUE,
    n_outputs    = 4L,
    experimental = TRUE,
    NA_allowed=TRUE),
  beta_binomial_window = list(
    description = paste(
      "Beta-Binomial learning kernel with sliding window on memory:\n",
      "        k = predicted probability of the current trial's observation\n",
      "        being X as opposed to Y.\n",
      "        Assumes a binary (Bernoulli) sequence of observations.\n",
      "        Limits memory of past events to a fixed window.\n",
      "        Parameters: a0, b0 (shape parameters of the Beta prior), window."
    ),
    default_pars = c("a0", "b0", "window"),
    transforms   = list(func = list("a0" = "exp", "b0" = "exp", "window" = "exp")),
    bases        = base_2p,
    sequential   = TRUE,
    n_outputs    = 4L,
    experimental = TRUE,
    NA_allowed=TRUE),
  dbm = list(
    description = paste(
      "Dynamic Belief Model (DBM) kernel:\n",
      "        k = predicted probability of the current trial's observation\n",
      "        being X as opposed to Y.\n",
      "        Assumes a binary (Bernoulli) sequence of observations.\n",
      "        Assumes that the trial-wise latent probability that the observation is X\n",
      "        is a mixture of the previous trial's posterior belief about that probability\n",
      "        and a fixed prior belief.",
      "        The parameter controlling this mixture can be interpreted as the assumed\n",
      "        probability of a change-point in the environment (i.e., volatility).\n",
      "        Parameters: cp (change-point probability),\n",
      "        mu0 (mean of Beta prior), s0 (scale of Beta prior)."
    ),
    default_pars = c("cp", "mu0", "s0"),
    transforms   = list(func = list("cp" = "pnorm", "mu0" = "pnorm", "s0" = "exp")),
    bases        = base_2p,
    sequential   = TRUE,
    n_outputs    = 4L,
    experimental = TRUE,
    NA_allowed=TRUE),
  tpm = list(
    description = paste(
      "Transition Probability Model (TPM) kernel:\n",
      "        k = predicted probability of the current trial's observation\n",
      "        being X as opposed to Y, conditional on the previous trial's observation.\n",
      "        Assumes a binary (Bernoulli) sequence of observations and\n",
      "        estimates first-order transition probabilities.\n",
      "        Incorporates a trial-wise change-point probability cp[i] controlling\n",
      "        belief volatility.\n",
      "        Note that k(t) is a one-step-ahead prediction: It uses only observations\n",
      "        up to and including trial t-1, and is computed before trial t's observation\n",
      "        is known. This is in contrast to the original Matlab toolbox implementation\n",
      "        (Meyniel et al., 2016), where the output for trial t is the model's forecast\n",
      "        for trial t+1, computed from observations through trial t.",
      "        Parameters: cp (change-point probability),\n",
      "        a0, b0 (shared shape parameters of the two - assumed independent -\n",
      "        Beta priors over transition probabilities p(X|X) and p(X|Y))."
    ),
    default_pars = c("cp", "a0", "b0"),
    transforms   = list(func = list("cp" = "pnorm", "a0" = "exp", "b0" = "exp")),
    bases        = base_2p,
    sequential   = TRUE,
    n_outputs    = 4L,
    experimental = TRUE,
    NA_allowed=TRUE)
  # delta2kernel2 = list(description = paste(
  #               "Steven fucking around with the delta2kernel. You shouldn't see this! Dual kernel delta rule: k = q[i].\n",
  #               "         Combines fast and slow learning rates\n",
  #               "         and switches between them based on dSwitch.\n",
  #               "         Parameters: q0 (initial value), alphaFast (fast learning rate),\n",
  #               "         propSlow (alphaSlow = propSlow * alphaFast), dSwitch (switch threshold)."
  #             ),
  #             default_pars = c("q0", "alphaFast", "propSlow", "dSwitch"),
  #             transforms = list(func = list("q0" = "identity", "alphaFast" = "pnorm",
  #                                           "propSlow" = "pnorm", "dSwitch" = "pnorm")),
  #             bases = base_2p),
              )
  kernels
}


format_kernel <- function(kernel, kernel_pars=NULL) {
  kernels <- get_kernels()
  eq_string <- kernels[[kernel]]$description
  eq_string <- strsplit(eq_string, ': k = ')[[1]][[2]]
  if(kernel %in% .sequential_kernels()) eq_string <- strsplit(eq_string, '\\.')[[1]][[1]]
  if(kernel %in% c('exp_incr', 'pow_incr', 'poly1', 'poly2', 'poly3', 'poly4')) eq_string <- paste0('(', eq_string, ')')

  # add placeholders
  if(!is.null(kernel_pars)) {
    for(kernel_par_n in 1:length(kernel_pars)) {
      old_name <- kernels[[kernel]]$default_pars[kernel_par_n]
      eq_string <- gsub(old_name, kernel_pars[kernel_par_n], eq_string)
    }
  }

  eq_string
}

format_base <- function(base) {
  bases <- get_bases()
  eq_string <- bases[[base]]$description
  eq_string <- strsplit(eq_string, ': ')[[1]][[2]]
  eq_string
}


verbal_trend <- function(design_matrix, trend) {
  dm_cn <- colnames(design_matrix)
  bases <- trend$bases
  trend_str <- c()
  trend_par_names <- c()

  for(base in bases) {
    base_type <- base$type
    trend_par_name <- base$target_parameter
    kernel <- trend$kernels[[base$kernel_id]]
    kernel_type <- kernel$type
    covariate <- kernel$cov_names
    base_pars <- base$pnames
    kernel_pars <- if (length(kernel$pnames) > 0) kernel$pnames else NULL

    # Only emit a string for bases whose target parameter appears in the design matrix
    if (!trend_par_name %in% dm_cn) next

    # format kernel and base
    kernel_formatted <- format_kernel(kernel_type, kernel_pars=kernel_pars)
    base_formatted <- format_base(base_type)
    plain_par_name <- sub("_t$", "", base$target_parameter)  # always strip, safe if no _t
    display_name <- paste0(plain_par_name, "_t")              # LHS always gets _t
    base_formatted <- paste0(display_name, " = ", base_formatted)

    # replace all in one go, use placeholders to prevent cascading replacements
    w_par <- if (length(base_pars) > 0) base_pars[1] else NULL
    replacements <- c('k'=gsub('c', covariate[1], kernel_formatted), 'w'=w_par, 'parameter'=plain_par_name)
    if(base$phase != 'premap') replacements['parameter + '] <- ''

    patterns <- names(replacements)
    placeholders <- paste0("___PLACEHOLDER", seq_along(patterns), "___")
    for (i in seq_along(patterns)) {
      base_formatted <- gsub(patterns[i], placeholders[i], base_formatted, fixed = TRUE)
    }
    for (i in seq_along(placeholders)) {
      base_formatted <- gsub(placeholders[i], replacements[i], base_formatted, fixed = TRUE)
    }

    # Add additional covariates
    if(length(covariate) > 1) {
      for(cov_n in 2:length(covariate)) {
        kernel_formatted <- format_kernel(kernel_type, kernel_pars=kernel_pars)
        additional_base <- format_base(base_type)

        w_par <- if (length(base_pars) > 0) base_pars[1] else NULL
        replacements <- c('k'=gsub('c', covariate[cov_n], kernel_formatted), 'w'=w_par, 'parameter + ' = '')
        patterns <- names(replacements)
        placeholders <- paste0("___PLACEHOLDER", seq_along(patterns), "___")
        for (i in seq_along(patterns)) {
          additional_base <- gsub(patterns[i], placeholders[i], additional_base, fixed = TRUE)
        }
        for (i in seq_along(placeholders)) {
          additional_base <- gsub(placeholders[i], replacements[i], additional_base, fixed = TRUE)
        }

        base_formatted <- paste0(base_formatted, ' + ', additional_base)
      }
    }
    trend_str <- c(trend_str, base_formatted)
    trend_par_names <- c(trend_par_names, trend_par_name)

  }
  if(length(trend_str) > 0) {
    trend_str <- setNames(trend_str, trend_par_names)
  }
  trend_str
}






## Functions for working with custom kernels

##' Extract pointers of custom C++ trend kernels from trend list or emc object
##'
##' Extract the pointers so they can be re-added to an emc object after loading it from disk.
##'
##' @param input_data Either an emc object or an `emc2_trend` object.
##' @return A list-of-lists of custom pointers, one sub-list per model, keyed
##'   by `kernel_id`. Kernels without a custom pointer return `NULL`. Returns
##'   `NULL` if no custom pointers are found anywhere.
get_custom_kernel_pointers <- function(input_data) {
  if(inherits(input_data, 'emc')) {
    if(is.list(input_data[[1]]$model)) {
      # joint model: list of model functions
      trend <- lapply(input_data[[1]]$model, function(x) x()$trend)
    } else {
      # single model function — wrap for uniform structure
      trend <- list(input_data[[1]]$model()$trend)
    }
  } else {
    # raw trend list — wrap for uniform structure
    trend <- list(input_data)
  }

  if(all(sapply(trend, is.null))) return(NULL)

  # make list
  # ptr_list <- lapply(trend, function(x) lapply(x, function(y) attr(y, 'custom_ptr')))
  ptr_list <- lapply(trend, function(t) {
    if (!inherits(t, "emc2_trend")) return(NULL)
    lapply(t$kernels, function(k) k$kernel_pointer)
  })
  if(all(sapply(ptr_list, function(ptrs) all(sapply(ptrs, is.null))))) return(NULL)

  return(ptr_list)
}


##' (Re-)Set pointers of custom C++ trend kernels to an emc object
##'
##' When an emc object is loaded from disk, or returned by forked processes, the pointers to custom kernels
##' need to be re-created. This is a convenience function to do this.
##'
##' @param emc An emc object
##' @param ptrs A list of pointer lists as generated by
##'   `get_custom_kernel_pointers()`, keyed by `kernel_id`.
##' @return An emc object with the custom pointers re-instated.
set_custom_kernel_pointers <- function(emc, ptrs) {
  if(is.null(ptrs)) return(emc)   # nothing to set

  # validation
  if(length(ptrs) != length(emc[[1]]$model)) {
    stop('List of potential pointers not equal to number of models')
  }
  if(!is.list(emc[[1]]$model) && is.null(emc[[1]]$model()$trend)) {
    stop('emc object has no trends, nothing to set...')
  }

  set_ptrs_on_model <- function(model_fn, ptrs_for_model) {
    if(length(ptrs_for_model) == 0) return(model_fn)
    model_list <- model_fn()
    trend <- model_list$trend
    # Iterate by kernel_id name so order doesn't matter
    for (kid in names(ptrs_for_model)) {
      if (!is.null(ptrs_for_model[[kid]]) && !is.null(trend$kernels[[kid]]))
        trend$kernels[[kid]]$kernel_pointer <- ptrs_for_model[[kid]]
    }
    model_list$trend <- trend
    # Modify the existing closure environment in-place rather than creating a new closure
    environment(model_fn)$model_list <- model_list
    return(model_fn)
  }

  # fix models
  for(chain_ in seq_along(emc)) {
    if('model' %in% names(emc[[chain_]])) {
      if(is.list(emc[[chain_]]$model)) {
        for(j in seq_along(emc[[chain_]]$model)) {
          emc[[chain_]]$model[[j]] <- set_ptrs_on_model(emc[[chain_]]$model[[j]], ptrs[[j]])
        }
      } else {
        emc[[chain_]]$model <- set_ptrs_on_model(emc[[chain_]]$model, ptrs[[1]])
      }
    }
    # fix models hidden in priors
    if('prior' %in% names(emc[[chain_]])) {
      for(design_n in seq_along(attr(emc[[chain_]]$prior, 'design'))) {
        if('model' %in% names(attr(emc[[chain_]]$prior, 'design')[[design_n]])) {
          # only 1 model here by definition - priors are per submodel
          attr(emc[[chain_]]$prior, 'design')[[design_n]]$model <- set_ptrs_on_model(attr(emc[[chain_]]$prior, 'design')[[design_n]]$model, ptrs[[design_n]])
        }
      }
    }
  }
  return(emc)
}


# pointer_reset_wrapper <- function(sub_emc, emc){
#   # TO FIX, make custom kernel pointers work with joint models
#   # for now just return no updates
#   if(is.list(emc[[1]]$model)){ # Joint model!!
#     return(sub_emc)
#   } else{
#     return(set_custom_kernel_pointers(sub_emc, get_custom_kernel_pointers(emc)))
#   }
# }

##' Reset pointers of custom C++ trend kernels to an emc object
##'
##' When an emc object is loaded from disk, or returned by forked processes, the pointers to custom kernels
##' need to be re-created. This is a convenience function to do this.
##'
##' @param emc A target emc object with missing pointers
##' @param pointer_source Either an `emc2_trend` object with valid pointers, or
##'   an emc object with valid pointers.
##' @param model_number If `pointer_source` is a raw trend list and `emc` is a joint model, specifies
##'   which model the trend belongs to. Defaults to 1.
##' @return An emc object with the custom pointers re-instated.
##' @export
fix_custom_kernel_pointers <- function(emc, pointer_source, model_number = 1) {
  ptrs <- get_custom_kernel_pointers(pointer_source)
  if(!inherits(pointer_source, 'emc') && is.list(emc[[1]]$model)) {
    # raw trend list + joint model: slot pointers into the correct position
    n_models <- length(emc[[1]]$model)
    ptrs_joint <- vector('list', n_models)
    ptrs_joint[[model_number]] <- ptrs[[1]]
    ptrs <- ptrs_joint
  }
  return(set_custom_kernel_pointers(emc, ptrs))
}




# Simulation --------------------------------------------------------------
make_data_unconditional <- function(data, pars, design, model,
                                    return_trialwise_parameters,
                                    kernel_output_codes = c(1L),
                                    optionals = NULL,
                                    n_context_trials = 1L) {
  model_fun  <- model
  model_list <- model()
  includeColumns <- colnames(data)

  # -----------------------------------------------------------------------
  # Step 1: Build the full dadm ONCE for all subjects and trials.
  # -----------------------------------------------------------------------
  dadm_full <- design_model(
    add_accumulators(data, design$matchfun, simulate = FALSE,
                     type = model_list$type, Fcovariates = design$Fcovariates),
    design, model_fun, add_acc = FALSE, compress = FALSE,
    verbose = FALSE, rt_check = FALSE, compress_dms = FALSE
  )
  if (!"R"  %in% names(dadm_full)) dadm_full$R  <- NA
  if (!"rt" %in% names(dadm_full)) dadm_full$rt <- NA

  # Name Flist by the LHS of each formula if not already named
  if (is.null(names(design$Flist))) {
    names(design$Flist) <- sapply(design$Flist, function(f)
      as.character(stats::terms(f)[[2]])
    )
  }

  # Number of accumulators (rows per trial)
  n_acc <- sum(dadm_full$trials == dadm_full$trials[1] &
                 dadm_full$subjects == dadm_full$subjects[1])

  # -----------------------------------------------------------------------
  # Step 2: Set up design cache.
  # -----------------------------------------------------------------------
  factor_cols <- setdiff(names(design$Ffactors), "subjects")
  ffun_cols   <- names(design$Ffunctions)
  p_types     <- names(design$Flist)

  pnames <- names(design$Flist)
  if (!is.list(design$Clist[[1]])) {
    design$Clist <- stats::setNames(
      lapply(seq_along(pnames), function(x) design$Clist),
      pnames
    )
  } else {
    missing_p_types <- pnames[!(pnames %in% names(design$Clist))]
    if (length(missing_p_types) > 0) {
      nok <- length(design$Clist)
      for (i in seq_along(missing_p_types)) {
        design$Clist[[missing_p_types[i]]] <- list(stats::contr.treatment)
        names(design$Clist)[nok + i] <- missing_p_types[i]
      }
    }
  }
  for (i in pnames) attr(design$Flist[[i]], "Clist") <- design$Clist[[i]]

  uses_ffun <- sapply(p_types, function(x) {
    any(all.vars(design$Flist[[x]]) %in% ffun_cols)
  })
  cached_pars   <- p_types[!uses_ffun]
  uncached_pars <- p_types[ uses_ffun]

  # Split parameter_design output rows into cached vs uncached
  # (uncached if any non-zero weight column is a function)
  cached_pd_pars <- if (!is.null(design$parameter_design)) rownames(design$parameter_design$weights) else character(0)
  uncached_pd_pars <- character(0)

  make_designs_cached <- local({
    cache <- list()
    function(dadm_slice, key) {
      if (is.null(cache[[key]])) {
        regular <- lapply(
          stats::setNames(cached_pars, cached_pars),
          function(x) make_dm(design$Flist[[x]], da = dadm_slice,
                              Fcovariates = design$Fcovariates,
                              compress_dms = FALSE)
        )
        pd_cached <- if (length(cached_pd_pars) > 0) {
          expand_parameter_design(
            list(weights = design$parameter_design$weights[cached_pd_pars, , drop = FALSE]),
            dadm_slice, compress_dms = FALSE
          )
        } else list()

        cache[[key]] <<- c(regular, pd_cached)
      }

      fresh_regular <- if (length(uncached_pars) > 0) {
        lapply(
          stats::setNames(uncached_pars, uncached_pars),
          function(x) make_dm(design$Flist[[x]], da = dadm_slice,
                              Fcovariates = design$Fcovariates,
                              compress_dms = FALSE)
        )
      } else list()

      all_designs <- c(cache[[key]], fresh_regular)
      pd_names    <- cached_pd_pars
      all_designs[c(p_types, pd_names[pd_names %in% names(all_designs)])]
      # message("p_types: ", paste(p_types, collapse=", "))
      # message("names(all_designs): ", paste(names(all_designs), collapse=", "))
      # message("ffun_cols: ", paste(ffun_cols, collapse=", "))
      # message("cached_pars: ", paste(cached_pars, collapse=", "))
      # message("uncached_pars: ", paste(uncached_pars, collapse=", "))
    }
  })


  # Identify whether any trend has covariate coding
  trend        <- model_list$trend  # list(kernels = list(...), bases = list(...))
  has_trend    <- !is.null(trend)

  # Step 5: covariate coding — collect all bases that have coding schemes
  bases_with_coding <- if (has_trend) {
    Filter(function(b) !is.null(b$coding), trend$bases)
  } else list()
  has_covariate_coding <- length(bases_with_coding) > 0

  # Step 9: feedback — collect all kernels that have feedback functions
  kernels_with_feedback <- if (has_trend) {
    Filter(function(k) !is.null(k$feedback), trend$kernels)
  } else list()
  has_feedback <- length(kernels_with_feedback) > 0

  has_ffunctions <- !is.null(design$Ffunctions)

  # -----------------------------------------------------------------------
  # Step 3: Per-subject, per-trial loop.
  # -----------------------------------------------------------------------
  trialwise_parameters <- NULL
  subj_levels <- levels(dadm_full$subjects)
  constants   <- attr(dadm_full, "constants")
  if (is.null(constants)) constants <- NA

  for (subj in subj_levels) {
    sub_trialwise_parameters <- NULL
    subj_mask <- dadm_full$subjects == subj
    if (!any(subj_mask)) next
    subj_rows <- which(subj_mask)

    dadm_subj   <- dadm_full[subj_rows, , drop = FALSE]
    trial_vals  <- sort(unique(dadm_subj$trials))
    n_rows_subj <- nrow(dadm_subj)

    idx_by_trial <- split(seq_len(n_rows_subj), dadm_subj$trials)

    # Helper: returns indices for up to n_context_trials previous trials
    # concatenated with the current trial's indices.
    get_context_idx <- function(j) {
      if (j == 1L || n_context_trials == 0L) return(idx_by_trial[[j]])
      lookback  <- seq(max(1L, j - n_context_trials), j - 1L)
      ctx_idx   <- unlist(idx_by_trial[lookback], use.names = FALSE)
      c(ctx_idx, idx_by_trial[[j]])
    }

    particle_matrix <- matrix(
      as.numeric(pars[which(subj == subj_levels), , drop = FALSE]),
      nrow = 1
    )
    colnames(particle_matrix) <- colnames(pars)

    # Pre-allocate designs_prefix from dadm_full's designs, zeroed out
    # Includes parameter_design matrices already (appended in design_model())
    designs_prefix <- lapply(attr(dadm_full, "designs"), function(m) {
      out <- m[subj_rows, , drop = FALSE]
      attr(out, "parameter_design") <- attr(m, "parameter_design")
      out[] <- 0
      out
    })

    # Pre-allocate covariate_coding_prefix: bootstrap structure from trial 1
    if (has_covariate_coding) {
      idx_t1  <- idx_by_trial[[1]]
      dadm_t1 <- dadm_subj[idx_t1, , drop = FALSE]

      covariate_coding_prefix <- list()
      for (base in bases_with_coding) {
        kernel    <- trend$kernels[[base$kernel_id]]
        cov_names <- kernel$cov_names
        for (map_name in names(base$coding)) {
          trial1_result <- base$coding[[map_name]](dadm = dadm_t1, cov_names)
          if (is.null(dim(trial1_result))) {
            trial1_result <- matrix(trial1_result, nrow = 1,
                                    dimnames = list(NULL, names(trial1_result)))
          }
          covariate_coding_prefix[[map_name]] <- matrix(
            0,
            nrow = n_rows_subj,
            ncol = ncol(trial1_result),
            dimnames = list(NULL, colnames(trial1_result))
          )
        }
      }
    }

    dadm_subj_df <- as.list(dadm_subj)
    class(dadm_subj_df) <- "data.frame"
    attr(dadm_subj_df, "row.names") <- .set_row_names(n_rows_subj)

    R_col  <- match("R",  names(dadm_subj_df))
    rt_col <- match("rt", names(dadm_subj_df))

    for (j in seq_along(trial_vals)) {
      current_trial        <- trial_vals[j]
      idx_curr             <- idx_by_trial[[as.character(current_trial)]]
      idx_ctx              <- get_context_idx(j)
      is_last_trial        <- j == length(trial_vals)
      tmp_return_trialwise <- is_last_trial && return_trialwise_parameters

      # 1. Materialize current trial slice
      dadm_current <- lapply(dadm_subj_df, `[`, idx_curr)
      class(dadm_current) <- "data.frame"
      attr(dadm_current, "row.names") <- .set_row_names(length(idx_curr))

      # 2. Ffunction pass 1: runs before get_pars_c_wrapper_oo.
      #    Handles Ffunctions that affect condition/design (e.g. depend on
      #    previous trial's R/rt/feedback, or purely on design factors).
      #    Updates dadm_current and dadm_subj_df so key + make_designs_cached
      #    see the correct values.
      # if (has_ffunctions) {
      #   for (i in names(design$Ffunctions)) {
      #     result <- design$Ffunctions[[i]](dadm_current)
      #     dadm_current[[i]]           <- result
      #     dadm_subj_df[[i]][idx_curr] <- result
      #   }
      # }
      if (has_ffunctions) {
        dadm_ctx <- lapply(dadm_subj_df, `[`, idx_ctx)
        class(dadm_ctx) <- "data.frame"
        attr(dadm_ctx, "row.names") <- .set_row_names(length(idx_ctx))

        for (i in names(design$Ffunctions)) {
          result_full             <- design$Ffunctions[[i]](dadm_ctx)
          result_curr             <- utils::tail(result_full, length(idx_curr))
          dadm_current[[i]]       <- result_curr
          dadm_subj_df[[i]][idx_curr] <- result_curr
        }
      }


      # 3. Compute condition key from updated dadm_subj_df
      key <- paste(vapply(factor_cols, function(fc)
        as.integer(dadm_subj_df[[fc]][idx_curr[1]]),
        integer(1)), collapse = "_")
      if (nchar(key) == 0) key <- "intercept_only"
      # key <- paste(vapply(factor_cols, function(fc)
      #   as.integer(dadm_subj_df[[fc]][idx_curr[1]]),
      #   integer(1)), collapse = "_")

      # 4. Get current-trial designs (cached + fresh) and write into prefix
      designs_current <- make_designs_cached(dadm_current, key)
      for (nm in names(designs_current)) {
        designs_prefix[[nm]][idx_curr, ] <- designs_current[[nm]]
      }

      # 5. Compute covariate coding for current trial only, write into buffer
      if (has_covariate_coding) {
        for (base in bases_with_coding) {
          kernel    <- trend$kernels[[base$kernel_id]]
          cov_names <- kernel$cov_names
          for (scheme_name in names(base$coding)) {
            result <- base$coding[[scheme_name]](dadm = dadm_current, cov_names)
            covariate_coding_prefix[[scheme_name]][idx_curr, ] <- as.matrix(result)
          }
        }
        attr(dadm_subj_df, "covariate_coding") <- covariate_coding_prefix
      }

      # 6. Get parameter matrix for full subject buffer
      pm <- get_pars_c_wrapper(
        particle_matrix      = particle_matrix,
        data                 = dadm_subj_df,
        constants            = constants,
        designs              = designs_prefix,
        bounds               = model_list$bound,
        transforms           = model_list$transform,
        pretransforms        = model_list$pre_transform,
        trend                = model_list$trend,
        return_kernel_matrix = FALSE,
        return_all_pars      = TRUE
      )

      if (tmp_return_trialwise && !is.null(model_list$trend)) {
        covariates <- get_pars_c_wrapper(
          particle_matrix      = particle_matrix,
          data                 = dadm_subj_df,
          constants            = constants,
          designs              = designs_prefix,
          bounds               = model_list$bound,
          transforms           = model_list$transform,
          pretransforms        = model_list$pre_transform,
          trend                = model_list$trend,
          return_kernel_matrix = TRUE,
          kernel_output_codes  = kernel_output_codes,
          return_all_pars      = TRUE
        )
        attr(pm, "trialwise_parameters") <- covariates
      }

      # 7. Ttransform + bounds on current-trial rows only
      pr <- model_list$Ttransform(pm[idx_curr, , drop = FALSE], dadm_current)

      if (!is.null(optionals$nobound)) {
        attr(pr, "ok") <- rep(TRUE, nrow(pr))
      } else {
        pr <- fix_bound(pr, model_list$bound, dadm_current$lR,
                        fix = !is.null(optionals$shrink2bound))
      }

      # 8. Simulate R and rt
      if (any(names(dadm_current) == "RACE")) {
        Rrt <- RACE_rfun(dadm_current, pr, model_fun)
      } else {
        Rrt <- model_list$rfun(dadm_current, pr)
      }

      dadm_subj_df[[R_col]][idx_curr]  <- Rrt[, "R"]
      dadm_subj_df[[rt_col]][idx_curr] <- Rrt[, "rt"]

      # 9. Feedback functions (trend)
      if (has_feedback) {
        dadm_current <- lapply(dadm_subj_df, `[`, idx_curr)
        class(dadm_current) <- "data.frame"
        attr(dadm_current, "row.names") <- .set_row_names(length(idx_curr))

        for (kernel in kernels_with_feedback) {
          for (col_name in names(kernel$feedback)) {
            dadm_subj_df[[col_name]][idx_curr] <- kernel$feedback[[col_name]](dadm_current)
          }
        }
      }

      # 10. Ffunction pass 2: runs after simulation and feedback.
      #     Handles Ffunctions that depend on the current trial's R, rt, or
      #     feedback values. Updates dadm_subj_df only (dadm_current is no
      #     longer needed this trial); results are available to pass 1 of
      #     the next trial.
      #     Now context-aware
      if (has_ffunctions) {
        dadm_ctx <- lapply(dadm_subj_df, `[`, idx_ctx)
        class(dadm_ctx) <- "data.frame"
        attr(dadm_ctx, "row.names") <- .set_row_names(length(idx_ctx))

        for (i in names(design$Ffunctions)) {
          result_full             <- design$Ffunctions[[i]](dadm_ctx)
          dadm_subj_df[[i]][idx_curr] <- utils::tail(result_full, length(idx_curr))
        }
      }

      # 11. Collect trialwise parameters on last trial
      if (tmp_return_trialwise) {
        sub_trialwise_parameters <- as.data.frame(cbind(pm, attr(pm, "trialwise_parameters")))
        sub_trialwise_parameters$subject <- subj
        sub_trialwise_parameters$trial   <- rep(trial_vals, each = n_acc)
      }
    }

    if (return_trialwise_parameters) {
      trialwise_parameters <- rbind(trialwise_parameters, sub_trialwise_parameters)
    }

    # Write subject results back into dadm_full
    missing_in_full <- setdiff(names(dadm_subj_df), names(dadm_full))
    if (length(missing_in_full)) {
      for (nm in missing_in_full) dadm_full[[nm]] <- NA
    }
    dadm_full[subj_rows, names(dadm_subj_df)] <- dadm_subj_df
  }

  # -----------------------------------------------------------------------
  # Step 4: Final pass, trim output columns.
  # -----------------------------------------------------------------------
  if (n_acc > 1) {
    first_lR <- levels(dadm_full$lR)[1]
    dadm_full <- dadm_full[dadm_full$lR == first_lR, , drop = FALSE]
  }
  dadm_full <- dadm_full[, unique(c(includeColumns, "R", "rt")), drop = FALSE]
  dadm_full <- dadm_full[, !colnames(dadm_full) %in% c("lR", "lM"), drop = FALSE]

  list(data = dadm_full, trialwise_parameters = trialwise_parameters)
}

has_delta_rules <- function(model) {
  trend <- model()$trend
  if(is.null(trend)) return(FALSE)

  for(kernel_id in seq_along(trend$kernels)) {
    if(trend$kernels[[kernel_id]]$type %in% .sequential_kernels()) return(TRUE)
  }
  return(FALSE)
}


#' Apply a kernel implied in an emc object
#'
#' Applies the trend-specific kernel associated with an `emc` model to the
#' a subject's data and returns the resulting kernel matrix.
#'
#' @description
#' This function extracts the appropriate trend, and applies its implied kernel
#' using either the pure R implementation or the Rcpp version.
#' When `mode = "compare"`, the function checks whether the
#' two implementations produce identical output.
#'
#' @param kernel_pars A named vector of kernel parameters on the natural scale.
#'   Use `NULL` for kernels that do not require parameters.
#'
#' @param emc An `emc`.
#'
#' @param subject Subject index for which to apply the kernel. Defaults to `1`.
#'
#' @param input_pars Optional parameter matrix containing externally supplied
#'   parameter values (e.g., trend parameters). Only needed for custom kernels.
#'
#' @return
#' Returns a kernel matrix produced by the corresponding implementation.
#' @export
apply_kernel <- function(kernel_pars, emc, subject=1, input_pars=NULL) {
  dadm <- emc[[1]]$data[[subject]]
  model <- emc[[1]]$model()
  trend <- model$trend
  if(length(trend$kernels) > 1) {
    stop("Multiple kernels provided. This function only supports a single trend")
  }

  # Build parameter matrix
  p_vector <- sampled_pars(emc)
  if(!is.null(kernel_pars)) {
    p_vector[names(p_vector) %in% names(kernel_pars)] <- kernel_pars
  }
  if(!is.null(input_pars)) {
    p_vector[names(p_vector) %in% colnames(input_pars)] <- input_pars[1]
  }
  p_mat <- t(as.matrix(p_vector))
  colnames(p_mat) <- names(p_vector)
  out <- get_pars_oo(p_mat, dadm, model, return_kernel_matrix = TRUE)
  out <- out[, grepl(paste0("^", names(trend$kernels)[1], "\\."), colnames(out)), drop = FALSE]
  colnames(out) <- trend$kernels[[1]]$covariate
  out
}
