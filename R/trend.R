# Custom kernel: operate on all input columns at once; compress by at; exclude rows with any NA; expand back
run_kernel_custom <- function(trend_pars = NULL, input, funptr, at_factor = NULL) {
  if (!is.matrix(input)) input <- matrix(input, ncol = 1)
  n <- nrow(input)
  if (is.null(trend_pars)) trend_pars <- matrix(nrow = n, ncol = 0)

  # Compress to first-level rows when at provided
  if (!is.null(at_factor)) {
    if (!is.factor(at_factor)) stop("'at' column must be a factor")
    first_level <- at_factor == levels(at_factor)[1]
    expand_idx <- make_expand_idx(first_level)
    input_comp <- input[first_level, , drop = FALSE]
    tpars_comp <- if (nrow(trend_pars)) trend_pars[first_level, , drop = FALSE] else matrix(nrow = sum(first_level), ncol = 0)
  } else {
    expand_idx <- seq_len(n)
    input_comp <- input
    tpars_comp <- trend_pars
  }

  # Exclude any rows with at least one NA across columns
  good <- rowSums(is.na(input_comp)) == 0
  comp_out <- numeric(nrow(input_comp))
  if (any(good)) {
    in_good <- input_comp[good, , drop = FALSE]
    tp_good <- if (ncol(tpars_comp)) tpars_comp[good, , drop = FALSE] else matrix(nrow = sum(good), ncol = 0)
    contrib <- EMC2_call_custom_trend(tp_good, in_good, funptr)
    contrib[is.na(contrib)] <- 0
    comp_out[good] <- contrib
  }
  # Expand back to full rows, return as single-column matrix
  matrix(comp_out[expand_idx], ncol = 1)
}

#' Create a trend specification for model parameters
#'
#' @param par_names Character vector specifying which parameters to apply trend to
#' @param cov_names Character vector specifying which covariates to use for each trend
#' @param kernels Character vector specifying which kernel function to use for each trend
#' @param bases Optional character vector specifying which base function to use for each trend
#' @param shared Named list with entries the parameter names to be shared and the names the new names of the shared parameter.
#' @param trend_pnames Optional character vector specifying custom parameter names
#' @param phase Character vector (length 1 or `length(par_names)`) specifying the phase for each trend entry;
#'        one of "premap", "pretransform", or "posttransform". Defaults to "premap".
#' @param par_input Optional character vector(s) of parameter names to use as additional inputs for the trend
#' @param at If NULL (default), trend is applied everywhere. If a factor name (e.g., "lR"), trend is applied only to entries corresponding to the first level of that factor.
#' @param custom_trend A trend registered with `register_trend`
#'
#' @return A list containing the trend specifications for each parameter
#' @export
#'
#' @examples
#' # Put trend on B and v parameters
#' trend <- make_trend(
#'   par_names = c("B", "v"),
#'   cov_names = "strial",
#'   kernels = c("exp_incr", "poly3"),
#'   phase = "premap",
#'   shared = list(shrd = list("B.B0", "v.d1"))
#' )
#' get_trend_pnames(trend)
#'
make_trend <- function(par_names, cov_names = NULL, kernels, bases = NULL,
                       shared = NULL, trend_pnames = NULL,
                       phase = "premap",
                       par_input = NULL, at = NULL,
                       custom_trend = NULL){
  if(!(length(par_names) == length(kernels))){
    stop("Make sure that par_names and kernels have the same length")
  }
  if (is.null(cov_names)) {
    cov_names <- rep(list(character(0)), length(par_names))
  } else if(length(cov_names) != length(par_names)){
    if(length(cov_names) == 1){
      cov_names <- rep(cov_names, length(par_names))
    } else{
      stop("Make sure that cov_names and par_names have the same length")
    }
  }
  # normalize par_input to align with par_names (each entry vector of names or character(0))
  if (is.null(par_input)) {
    par_input <- rep(list(character(0)), length(par_names))
  } else if (length(par_input) != length(par_names)) {
    if (length(par_input) == 1) {
      par_input <- rep(par_input, length(par_names))
    } else {
      stop("Make sure that par_input and par_names have the same length or par_input is NULL/length 1")
    }
  }
  # normalize phase
  if (length(phase) != length(par_names)) phase <- rep(phase, length(par_names))
  if (!all(phase %in% c("premap","pretransform","posttransform"))) {
    stop("phase must be one of 'premap', 'pretransform', 'posttransform'")
  }
  if(!is.null(bases)){
    if(length(kernels) != length(bases)){
      stop("If bases not is NULL make sure you specify as many bases as kernels")
    }
  }
  # Normalize custom_trend to a per-parameter list if provided
  if (!is.null(custom_trend)) {
    if (inherits(custom_trend, "emc2_custom_trend")) {
      custom_trend <- rep(list(custom_trend), length(par_names))
    } else if (is.list(custom_trend)) {
      if (length(custom_trend) == 1) custom_trend <- rep(custom_trend, length(par_names))
      if (length(custom_trend) != length(par_names))
        stop("custom_trend must be a single registered trend or a list aligned with par_names")
      if (!all(vapply(custom_trend, inherits, logical(1), what = "emc2_custom_trend")))
        stop("Items in custom_trend must be created by register_trend()")
    } else {
      stop("custom_trend must be NULL, a single registered trend, or a list of them")
    }
  }

  trends_out <- list()
  for(i in 1:length(par_names)){
    trend <- list()
    # Kernel
    if(!kernels[i] %in% names(trend_help(kernels[i], return_types = TRUE)$kernels)){
      stop("Kernel type not support see `trend_help()`")
    } else  {
      trend$kernel <- kernels[i]
      trend$at <- at
    }

    # base
    if (is.null(bases)) {
      if (identical(kernels[i], "custom")) {
        if (is.null(custom_trend)) stop("custom_trend must be provided when using kernel='custom'")
        ct <- custom_trend[[i]]
        trend$base <- ct$base
      } else {
        trend$base <- trend_help(kernels[i], do_return = TRUE)$bases[1]
      }
    } else {
      if (identical(kernels[i], "custom")) {
        # For custom kernels, accept any of the standard bases the user specifies.
        base_ok <- c("lin","exp_lin","centered","add","identity")
        if (!(bases[i] %in% base_ok)) stop("Unknown base '", bases[i], "' for custom kernel. Pick one of ", paste(base_ok, collapse = ", "))
        trend$base <- bases[i]
      } else {
        if(bases[i] %in% names(trend_help(kernels[i], do_return = TRUE)$bases)){
          stop("base type not supported with kernel, see `trend_help(<kernel>)`")
        }
        trend$base <- bases[i]
      }
    }
    # add par names
    trend_pnames <- trend_help(base = trend$base, do_return = TRUE)$default_pars
    # Kernel parameter names:
    if (identical(kernels[i], "custom")) {
      if (is.null(custom_trend)) stop("custom_trend must be provided when using kernel='custom'")
      ct <- custom_trend[[i]]
      trend_pnames <- c(trend_pnames, ct$trend_pnames)
      # Attach the external pointer and optional transforms to this trend entry
      attr(trend, "custom_ptr") <- attr(ct, "custom_ptr")
      if (!is.null(attr(ct, "custom_transforms"))) {
        # ensure order matches ct$trend_pnames
        ctf <- attr(ct, "custom_transforms")
        if (!is.null(names(ctf))) ctf <- ctf[ct$trend_pnames]
        attr(trend, "custom_transforms") <- unname(ctf)
      }
    } else {
      trend_pnames <- c(trend_pnames, trend_help(kernel = kernels[i], do_return = TRUE)$default_pars)
    }
    trend$trend_pnames <- paste0(par_names[i], ".", trend_pnames)
    trend$covariate <- unlist(cov_names[i])
    trend$par_input <- unlist(par_input[[i]])
    trend$phase <- phase[i]
    trends_out[[par_names[i]]] <- trend
  }
  if(!is.null(shared)){
    # For each group of shared parameters
    for (main_par in names(shared)) {
      # Get the parameters to be replaced
      to_replace <- shared[[main_par]]
      # Loop through all trends
      for (trend_name in names(trends_out)) {
        # Get current trend parameter names
        curr_pnames <- trends_out[[trend_name]]$trend_pnames

        # Check if any of the parameters to be replaced exist
        curr_pnames[curr_pnames %in% to_replace] <- main_par
        trends_out[[trend_name]]$trend_pnames <- curr_pnames
      }
    }
  }

  attr(trends_out, "sequential") <- any(kernels %in% c("delta", "delta2"))
  return(trends_out)
}


#' Get parameter types from trend object
#'
#' @param trend A trend object created by make_trend()
#' @return A character vector of parameter names used in the trend
#' @export
#' @examples
#' trend <- make_trend(par_names = "v", cov_names = "trial", kernels = "exp_incr")
#' get_trend_pnames(trend)
#'
get_trend_pnames <- function(trend){
  out <- unlist(lapply(trend, function(x) x$trend_pnames))
  names(out) <- NULL
  return(unique(out))
}

#' Get help information for trend kernels and bases
#'
#' @param kernel Character string specifying the kernel type to get information about
#' @param base Character string specifying the base type to get information about
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
trend_help <- function(kernel = NULL, base = NULL, ...){
  dots <- add_defaults(list(...), do_return = FALSE, return_types = FALSE)
  bases <- get_bases()
  base_2p <- names(bases)[1:3]
  base_1p <- names(bases)[4:5]
  kernels <- get_kernels()
  if(dots$return_types){
    return(list(kernels = kernels, bases = bases))
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
    cat("\nTrend options:\n")
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
      if(dots$do_return) return(bases[[base]])
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

# Helper to compute expand index from first-level mask
make_expand_idx <- function(first_level) {
  idx <- cumsum(first_level)
  if (any(idx == 0)) stop("Found rows before first 'at' level within subject. Cannot anchor expansion.")
  idx
}


run_kernel <- function(trend_pars = NULL, kernel, input, funptr = NULL, at_factor = NULL) {
  # input: vector or matrix; apply per column and sum contributions; handle NA by zeroing; optional at compression/expansion
  if (!is.matrix(input)) input <- matrix(input, ncol = 1)
  n <- nrow(input)
  if (is.null(trend_pars)) trend_pars <- matrix(nrow = n, ncol = 0)
  out <- rep(0.0, n)

  # Custom kernels: operate on full matrix at once; returns n x 1 matrix
  if (identical(kernel, "custom")) {
    if (is.null(funptr)) stop("Missing function pointer for custom kernel. Pass 'funptr'.")
    return(run_kernel_custom(trend_pars, input, funptr, at_factor))
  }

  # Precompute at compression/expansion and compressed trend parameters
  if (!is.null(at_factor)) {
    if (!is.factor(at_factor)) stop("'at' column must be a factor")
    first_level <- at_factor == levels(at_factor)[1]
    expand_idx <- make_expand_idx(first_level)
    tpars_comp <- if (nrow(trend_pars)) trend_pars[first_level, , drop = FALSE] else matrix(nrow = sum(first_level), ncol = 0)
    use_at <- TRUE
  } else {
    first_level <- rep(TRUE, n)
    expand_idx <- seq_len(n)
    tpars_comp <- trend_pars
    use_at <- FALSE
  }
  # Per-column contribution, then return matrix with one column per input
  cols <- ncol(input)
  out_mat <- matrix(0, nrow = n, ncol = cols)
  for (j in seq_len(ncol(input))) {
    covariate_full <- input[, j]
    # 1) Compress to first-level rows if at_factor provided
    covariate_comp <- covariate_full[first_level]

    # 2) Initialize compressed output with zeros
    comp_len <- length(covariate_comp)
    comp_out <- numeric(comp_len)

    # 3) Exclude NAs
    good <- !is.na(covariate_comp)
    if (any(good)) {
      # 4) Run kernel on good subset only
      if (kernel == "custom") {
        if (is.null(funptr)) stop("Missing function pointer for custom kernel. Pass 'funptr'.")
        # Build 1-col input matrix for custom kernel
        in_good <- matrix(covariate_comp[good], ncol = 1)
        tp_good <- if (ncol(tpars_comp)) tpars_comp[good, , drop = FALSE] else matrix(nrow = sum(good), ncol = 0)
        contrib <- EMC2_call_custom_trend(tp_good, in_good, funptr)
        contrib[is.na(contrib)] <- 0
        comp_out[good] <- contrib
      } else {
        # Built-in kernels (use only rows in 'good')
        # Access parameters by column index as before
        if (kernel == "lin_decr") {
          comp_out[good] <- -covariate_comp[good]
        } else if (kernel == "lin_incr") {
          comp_out[good] <- covariate_comp[good]
        } else if (kernel == "exp_decr") {
          comp_out[good] <- exp(-tpars_comp[good, 1] * covariate_comp[good])
        } else if (kernel == "exp_incr") {
          comp_out[good] <- 1 - exp(-tpars_comp[good, 1] * covariate_comp[good])
        } else if (kernel == "pow_decr") {
          comp_out[good] <- (1 + covariate_comp[good])^(-tpars_comp[good, 1])
        } else if (kernel == "pow_incr") {
          comp_out[good] <- 1 - (1 + covariate_comp[good])^(-tpars_comp[good, 1])
        } else if (kernel == "poly2") {
          comp_out[good] <- tpars_comp[good, 1] * covariate_comp[good] + tpars_comp[good, 2] * covariate_comp[good]^2
        } else if (kernel == "poly3") {
          comp_out[good] <- tpars_comp[good, 1] * covariate_comp[good] + tpars_comp[good, 2] * covariate_comp[good]^2 + tpars_comp[good, 3] * covariate_comp[good]^3
        } else if (kernel == "poly4") {
          comp_out[good] <- tpars_comp[good, 1] * covariate_comp[good] + tpars_comp[good, 2] * covariate_comp[good]^2 + tpars_comp[good, 3] * covariate_comp[good]^3 + tpars_comp[good, 4] * covariate_comp[good]^4
        } else if (kernel == "delta") {
          tmp <- run_delta(tpars_comp[good, 1], tpars_comp[good, 2], covariate_comp[good])
          comp_out[good] <- tmp
        } else if (kernel == "delta2") {
          tmp <- run_delta2(tpars_comp[good, 1], tpars_comp[good, 2], tpars_comp[good, 3], tpars_comp[good, 4], covariate_comp[good])
          comp_out[good] <- tmp
        } else {
          stop("Unknown kernel type")
        }
      }
    }

    # 5) Expand back to full subject rows and store into output matrix column
    out_mat[, j] <- comp_out[expand_idx]
  }
  out_mat
}

# Helper: Apply forward-fill to covariates when using 'at' filtering
apply_forward_fill <- function(values, dadm,at) {
  idx <- dadm[,at] == levels(dadm[,at])[1] # assumes first level occurs first within each subject
  values[!idx] <- NA
  # Forward-fill within each subject separately
  filled <- values
  subs <- levels(dadm$subjects)
  for (s in subs) {
    m <- dadm$subjects == s
    if (!any(m)) next
    filled[m] <- na_locf(filled[m], na.rm = FALSE)
  }
  if (any(is.na(filled))) {
    stop("Found NA after forward-fill. This should not happen.")
  }
  return(filled)
}

prep_trend_phase <- function(dadm, trend, pars, phase, return_trialwise_parameters = FALSE){
  # Apply only trends in the requested phase, sequentially
  tnames <- names(trend)
  all_remove <- character(0)
  if(return_trialwise_parameters) tpars <- list()
  for (idx in seq_along(trend)){
    cur_trend <- trend[[idx]]
    if (!identical(cur_trend$phase, phase)) next
    par <- tnames[idx]
    all_remove <- c(all_remove, cur_trend$trend_pnames)
    updated <- run_trend(dadm, cur_trend, pars[, par], pars[, cur_trend$trend_pnames, drop = FALSE], pars,
                             return_trialwise_parameters = return_trialwise_parameters)
    if(return_trialwise_parameters){
      trialwise_parameters <- attr(updated, "trialwise_parameters")
      colnames(trialwise_parameters) <- paste0(par, "_", c(cur_trend$covariate, cur_trend$par_input))
      tpars[[par]] <- trialwise_parameters
    }

    pars[,par] <- updated

  }
  if (length(all_remove)) pars <- pars[, !(colnames(pars) %in% unique(all_remove)), drop = FALSE]
  if(return_trialwise_parameters) attr(pars, "trialwise_parameters") <- do.call(cbind, tpars)
  return(pars)
}

# Probably no need to loop and idx subjects
run_trend <- function(dadm, trend, param, trend_pars, pars_full = NULL,
                      return_trialwise_parameters = FALSE){
  n_base_pars <- switch(trend$base,
                        lin = 1,
                        exp_lin = 1,
                        centered = 1,
                        add = 0,
                        identity = 0)

  # Fix dimension for single-column trend_pars
  if(is.null(dim(trend_pars))) trend_pars <- t(t(trend_pars))

  out <- numeric(nrow(dadm))
  cov_cols <- trend$covariate

  # Check if this is a delta-rule kernel requiring special handling
  is_delta_kernel <- trend$kernel %in% c('delta', 'delta2')
  use_at_filter <- !is.null(trend$at)

  # Build par_input columns if needed
  par_in_cols <- if (!is.null(trend$par_input)) trend$par_input else character(0)
  par_input_matrix <- NULL
  if (length(par_in_cols) > 0) {
    par_input_matrix <- matrix(NA_real_, nrow(dadm), length(par_in_cols))
    for (j in seq_along(par_in_cols)) {
      par_input_matrix[, j] <- pars_full[, par_in_cols[j]]
    }
  }

  # Build a single input matrix across covariates and par_input (match C++ behavior)
  cov_mat <- NULL
  if (length(cov_cols) > 0) {
    cov_mat <- matrix(NA_real_, nrow(dadm), length(cov_cols))
    for (j in seq_along(cov_cols)) cov_mat[, j] <- dadm[, cov_cols[j]]
  }
  input_matrix <- cov_mat
  if (!is.null(par_input_matrix)) {
    input_matrix <- if (is.null(input_matrix)) par_input_matrix else cbind(input_matrix, par_input_matrix)
  }

  # Extract kernel parameters (excluding base parameters)
  if (ncol(trend_pars) > n_base_pars) {
    kernel_pars <- trend_pars[, seq.int(n_base_pars + 1, ncol(trend_pars)), drop = FALSE]
  } else {
    kernel_pars <- matrix(nrow = nrow(trend_pars), ncol = 0)
  }
  funptr <- if (identical(trend$kernel, "custom")) attr(trend, "custom_ptr") else NULL

  if(return_trialwise_parameters) tlist <- list()
  for(s in 1:length(unique(dadm$subjects))){
    s_idx <- dadm$subjects == unique(dadm$subjects)[s]
    dat <- dadm[s_idx,]
    if (is.null(input_matrix)) {
      k_sum <- rep(0, sum(s_idx))
    } else {
      subset_input <- input_matrix[s_idx,, drop = FALSE]
      at_fac <- if (use_at_filter) dat[, trend$at] else NULL
      kern_mat <- run_kernel(kernel_pars[s_idx,,drop = FALSE], trend$kernel, subset_input,
                             funptr = funptr, at_factor = at_fac)
      if(return_trialwise_parameters){
        tlist[[s]] <- kern_mat
      }
      # Optional per-column map weighting prior to summing
      if (!is.null(trend$map)) {
        map_mat <- trend$map[s_idx,, drop = FALSE]
        kern_mat <- kern_mat * map_mat
      }
      # Sum across columns
      if (ncol(kern_mat) == 0) {
        k_sum <- rep(0, nrow(kern_mat))
      } else {
        k_sum <- rowSums(kern_mat)
      }
    }
    out[s_idx] <- out[s_idx] + k_sum
  }

  # Do the mapping
  out <- switch(trend$base,
                lin = param + trend_pars[,1]*out,
                exp_lin = exp(param) + trend_pars[,1]*out,
                centered = param + trend_pars[,1]*(out-.5),
                add = param + out,
                identity = out
  )
  if(return_trialwise_parameters) attr(out, "trialwise_parameters") <- do.call(rbind, tlist)
  return(out)
}

check_trend <- function(trend, covariates = NULL, model = NULL, formula = NULL) {
  if(!is.null(model)){
    # non-premap trend targets must be model parameters
    tnames <- names(trend)
    ok <- vapply(seq_along(trend), function(i){
      if (identical(trend[[i]]$phase, "premap")) return(TRUE)
      tnames[i] %in% names(model()$p_types)
    }, logical(1))
    if (!all(ok)) stop("pretransform/posttransform trend has a parameter name not in the model")
  }
  if (is.null(covariates)) stop("must specify covariates when using trend")
  covnames <- unlist(lapply(trend,function(x)x$covariate))
  if (!all(covnames %in% covariates)){
    stop("trend has covnames not in covariates")
  }
  # Premap + par_input: allowed. Scalars will be replicated to vector length in C++ mapping
  trend_pnames <- get_trend_pnames(trend)
  if (!is.null(formula)) {
    isin <-  trend_pnames %in% unlist(lapply(formula,function(x)all.vars(x)[1]))
    if(any(!isin)){
      # Add missing trend parameters to formula with intercept-only model
      formula <- c(formula, lapply(trend_pnames[!isin], function(x) as.formula(paste(x, "~ 1"))))
      message("Intercept formula added for trend_pars: ", paste(trend_pnames[!isin],collapse=", "))
    }
  }
  return(formula)
}


update_model_trend <- function(trend, model) {
  # Get model list to modify
  model_list <- model()

  # For each parameter in the trend
  tnames <- names(trend)
  for (i in seq_along(trend)) {
    par <- tnames[i]
    cur_trend <- trend[[i]]

    # Get default transforms from base and kernel
    base_transforms <- trend_help(base = cur_trend$base, do_return = TRUE)$transforms
    if (identical(cur_trend$kernel, "custom")) {
      ctf <- attr(cur_trend, "custom_transforms")
      kernel_transforms <- if (is.null(ctf)) NULL else list(func = ctf)
    } else {
      kernel_transforms <- trend_help(cur_trend$kernel, do_return = TRUE)$transforms
    }

    # Combine transforms
    if (!is.null(kernel_transforms) || !is.null(base_transforms)) {
      tmp <- c(base_transforms$func, kernel_transforms$func)
      names(tmp) <- cur_trend$trend_pnames
      # Update the appropriate transform list based on premap
      model_list$transform$func <- c(model_list$transform$func, unlist(tmp))
    }
    model_list$p_types <- c(model_list$p_types, setNames(numeric(length(cur_trend$trend_pnames)), cur_trend$trend_pnames))
  }
  # Ensure that shared parameters are removed
  model_list$p_types <- model_list$p_types[unique(names(model_list$p_types))]
  model_list$trend <- trend
  # Return updated model function
  model <- function() { return(model_list) }
  return(model)
}

run_delta <- function(q0,alpha,covariate) {
  q <- pe <- numeric(length(covariate))
  q[1] <- q0[1]
  if(length(q) == 1) return(q)  # only 1 trial, cannot be updated

  for (i in 2:length(q)) {
    pe[i-1] <- covariate[i-1]-q[i-1]
    q[i] <- q[i-1] + alpha[i-1]*pe[i-1]
  }
  return(q)
}

run_delta2 <- function(q0,alphaFast,propSlow,dSwitch,covariate) {
  q <- qFast <- qSlow <- peFast <- peSlow <- numeric(length(covariate))
  q[1] <- qFast[1] <- qSlow[1] <- q0[1]
  alphaSlow <- propSlow*alphaFast
  for (i in 2:length(covariate)) {
    peFast[i-1] <- covariate[i-1]-qFast[i-1]
    peSlow[i-1] <- covariate[i-1]-qSlow[i-1]
    qFast[i] <- qFast[i-1] + alphaFast[i-1]*peFast[i-1]
    qSlow[i] <- qSlow[i-1] + alphaSlow[i-1]*peSlow[i-1]
    if (abs(qFast[i]-qSlow[i])>dSwitch[i]){
      q[i] <- qFast[i]
    } else{
      q[i] <- qSlow[i]
    }
  }
  return(q)
}

##' Register a custom C++ trend kernel
##'
##' Compiles and registers a user-provided C++ function that maps per-trial
##' kernel parameters and inputs to a numeric vector. The C++ function must have
##' signature:
##'   NumericVector f(NumericMatrix trend_pars, NumericMatrix input)
##' and provide an exported pointer creator using EMC2_MAKE_PTR.
##'
##' @param trend_parameters Character vector of kernel parameter names (in order).
##' @param file Path to the C++ file implementing the custom kernel. The file
##'   should include EMC2/userfun.hpp and define a pointer creator (via
##'   EMC2_MAKE_PTR) that is exported to R.
##' @param transforms Optional named character vector or list mapping each custom
##'   kernel parameter name to a transform name (e.g., "identity", "exp", "pnorm").
##'   Length must match `trend_parameters`. If unnamed but the correct length,
##'   the order is assumed to match `trend_parameters`.
##' @param base Default base to use when creating trends with this custom kernel
##'   if no `bases` argument is supplied to `make_trend`. One of
##'   c("lin","exp_lin","centered","add","identity"). Default "add".
##' @return An object to pass to `make_trend(custom_trend=...)`, carrying the
##'   pointer, parameter names, default base, and optional transform mapping.
##' @export
register_trend <- function(trend_parameters, file, transforms = NULL, base = "add"){
  if (!is.character(trend_parameters) || length(trend_parameters) == 0)
    stop("trend_parameters must be a non-empty character vector")
  if (!file.exists(file)) stop("C++ file not found: ", file)
  base_ok <- c("lin","exp_lin","centered","add","identity")
  if (!is.character(base) || length(base) != 1L || !(base %in% base_ok))
    stop("base must be one of ", paste(base_ok, collapse = ", "))

  # Normalize transforms to a character vector in the order of trend_parameters
  trf_vec <- NULL
  if (!is.null(transforms)) {
    if (is.list(transforms)) transforms <- unlist(transforms, recursive = FALSE, use.names = TRUE)
    transforms <- unlist(transforms, use.names = TRUE)
    if (length(transforms) != length(trend_parameters)) {
      stop("length(transforms) must match length(trend_parameters)")
    }
    if (is.null(names(transforms))) {
      names(transforms) <- trend_parameters
      trf_vec <- as.character(transforms)
    } else {
      # Reorder to match trend_parameters
      if (!all(sort(names(transforms)) == sort(trend_parameters))) {
        stop("names(transforms) must match trend_parameters if names are supplied")
      }
      trf_vec <- as.character(transforms[trend_parameters])
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

  obj <- list(trend_pnames = as.character(trend_parameters),
              base = base,
              maker = maker,
              file = normalizePath(file))
  class(obj) <- "emc2_custom_trend"
  attr(obj, "custom_ptr") <- ptr
  if (!is.null(trf_vec)) attr(obj, "custom_transforms") <- trf_vec
  obj
}

# apply_lR_filter <- function(d, cov_name) {
#   if(!lR %in% colnames(d)) {
#     d[levels(d$lR)!=levels(d$lR)[1],cov_name] <- NA
#   }
#   d
# }

has_delta_rules <- function(model) {
  trend <- model()$trend
  if(is.null(trend)) return(FALSE)

  for(trend_n in 1:length(trend)) {
    if(trend[[trend_n]]$kernel %in% c('delta', 'delta2')) return(TRUE)
  }
  return(FALSE)
}

has_conditional_covariates <- function(design) {
  # find define covariates that depend on behavior -- these are rt, R, or any of the outputs of the functions provided.
  # they can be lRfiltered or not
  function_output_columns <- names(design$Ffunctions)
  behavioral_covariates <- c('rt', 'R', function_output_columns)

  # find actual covariates, look for a match
  trend <- design$model()$trend
  for(trend_n in 1:length(trend)) {
    for(cov in trend[[trend_n]]$covariate) {
      if(cov  %in% behavioral_covariates) return(TRUE)
    }
  }
  return(FALSE)
}

get_bases <- function() {
  bases <- list(
    lin = list(description = "Linear base: parameter + w * k",
               transforms = list(func = list("w" = "identity")),
               default_pars = "w"),
    exp_lin = list(description = "Exponential linear base: exp(parameter) + exp(w) * k",
                   transforms = list(func = list("w" = "exp")),
                   default_pars = "w"),
    centered = list(description = "Centered mapping: parameter + w*(k - 0.5)",
                    transforms = list(func = list("w" = "identity")),
                    default_pars = "w"),
    add = list(description = "Additive base: parameter + k",
               transforms = NULL),
    identity = list(description = "Identity base: k",
                    transforms = NULL)
  )
  bases
}

get_kernels <- function() {
  bases <- get_bases()
  base_2p <- names(bases)[1:3]
  base_1p <- names(bases)[4:5]

  kernels <- list(
    custom = list(description = "Custom C++ kernel: provided via register_trend().",
                  transforms = NULL,
                  default_pars = NULL,
                  bases = names(bases)),
    lin_decr = list(description = "Decreasing linear kernel: k = -c",
                    transforms = NULL,
                    default_pars = NULL,
                    bases = base_2p),
    lin_incr = list(description = "Increasing linear kernel: k = c",
                    transforms = NULL,
                    default_pars = NULL,
                    bases = base_2p),
    exp_decr = list(description = "Decreasing exponential kernel: k = exp(-d_ed * c)",
                    transforms = list(func =list("d_ed" = "exp")),
                    default_pars = "d_ed",
                    bases = base_2p),
    exp_incr = list(description = "Increasing exponential kernel: k = 1 - exp(-d_ei * c)",
                    transforms = list(func =list("d_ei" = "exp")),
                    default_pars = "d_ei",
                    bases = base_2p),
    pow_decr = list(description = "Decreasing power kernel: k = (1 + c)^(-d_pd)",
                    transforms = list(func =list("d_pd" = "exp")),
                    default_pars = "d_pd",
                    bases = base_2p),
    pow_incr = list(description = "Increasing power kernel: k = 1 - (1 + c)^(-d_pi)",
                    transforms = list(func =list("d_pi" = "exp")),
                    default_pars = "d_pi",
                    bases = base_2p),
    poly2 = list(description = "Quadratic polynomial: k = d1 * c + d2 * c^2",
                 transforms = list(func = list("d1" = "identity", "d2" = "identity")),
                 default_pars = c("d1", "d2"),
                 bases = base_1p),
    poly3 = list(description = "Cubic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3",
                 transforms = list(func = list("d1" = "identity", "d2" = "identity", "d3" = "identity")),
                 default_pars = c("d1", "d2", "d3"),
                 bases = base_1p),
    poly4 = list(description = "Quartic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3 + d4 * c^4",
                 transforms = list(func = list("d1" = "identity", "d2" = "identity", "d3" = "identity", "d4" = "identity")),
                 default_pars = c("d1", "d2", "d3", "d4"),
                 bases = base_1p),
    delta = list(description = paste(
      "Standard delta rule kernel: k = q[i].\n",
      "        Updates q[i] = q[i-1] + alpha * (c[i-1] - q[i-1]).\n",
      "        Parameters: q0 (initial value), alpha (learning rate)."
    ),
    default_pars = c("q0", "alpha"),
    transforms = list(func = list("q0" = "identity", "alpha" = "pnorm")),
    bases = base_2p),
    delta2 = list(description = paste(
      "Dual kernel delta rule: k = q[i].\n",
      "         Combines fast and slow learning rates\n",
      "         and switches between them based on dSwitch.\n",
      "         Parameters: q0 (initial value), alphaFast (fast learning rate),\n",
      "         propSlow (alphaSlow = propSlow * alphaFast), dSwitch (switch threshold)."
    ),
    default_pars = c("q0", "alphaFast", "propSlow", "dSwitch"),
    transforms = list(func = list("q0" = "identity", "alphaFast" = "pnorm",
                                  "propSlow" = "pnorm", "dSwitch" = "pnorm")),
    bases = base_2p)
  )
  kernels
}


format_kernel <- function(kernel, kernel_pars=NULL) {
  kernels <- get_kernels()
  eq_string <- kernels[[kernel]]$description
  eq_string <- strsplit(eq_string, ': k = ')[[1]][[2]]
  if(kernel %in% c('delta', 'delta2')) eq_string <- strsplit(eq_string, '\\.')[[1]][[1]]
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
  trend_par_names <- dm_cn[dm_cn %in% names(trend)]
  trend_str <- c()
  for(trend_par_name in trend_par_names) {
    base <- trend[[trend_par_name]]$base
    kernel <- trend[[trend_par_name]]$kernel
    covariate <- trend[[trend_par_name]]$covariate
    trend_pnames <- trend[[trend_par_name]]$trend_pnames
    n_base_pars <- switch(base,
                          lin = 1,
                          exp_lin = 1,
                          centered = 1,
                          add = 0,
                          identity = 0)
    if(length(trend_pnames) > n_base_pars) {
      kernel_pars <- trend_pnames[(n_base_pars+1):length(trend_pnames)]
    } else {
      kernel_pars <- NULL
    }
    base_pars=NULL
    if(n_base_pars > 0) base_pars <- trend_pnames[1:n_base_pars]

    # format kernel and base
    kernel_formatted <- format_kernel(kernel, kernel_pars=kernel_pars)
    base_formatted <- format_base(base)
    if(attr(trend, 'premap')) {
      trend_par_name <- gsub('_t', '', trend_par_name)
    }
    base_formatted <- paste0(trend_par_name, '_t = ', base_formatted)

    # replace all in one go, use placeholders to prevent cascading replacements
    replacements <- c('k'=gsub('c', covariate[1], kernel_formatted), 'w'=base_pars[1], 'parameter'=trend_par_name)
    if(!attr(trend, 'premap')) {
      replacements <- c('k'=gsub('c', covariate[1], kernel_formatted), 'w'=base_pars[1], 'parameter + ' = '')
    }
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
        kernel_formatted <- format_kernel(kernel, kernel_pars=kernel_pars)
        additional_base <- format_base(base)

        replacements <- c('k'=gsub('c', covariate[cov_n], kernel_formatted), 'w'=base_pars[1], 'parameter + ' = '')
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
  }
  if(length(trend_str) > 0) {
    trend_str <- setNames(trend_str, trend_par_names)
  }
  trend_str
}


make_data_unconditional <- function(data, pars, design, model, return_trialwise_parameters) {
  model_fun <- model
  model_list <- model()
  includeColumns <- colnames(data)
  # Initial scaffolding (attributes and factor setup)
  data <- design_model(
    add_accumulators(data,design$matchfun,simulate=FALSE,type=model_list$type,Fcovariates=design$Fcovariates),
    design,model_fun,add_acc=FALSE,compress=FALSE,verbose=FALSE,
    rt_check=FALSE)
  trialwise_parameters <- NULL
  # Iterate per subject, then per trial
  subj_levels <- levels(data$subjects)
  for (subj in subj_levels) {
    idx_subj_all <- which(data$subjects == subj)
    if (!length(idx_subj_all)) next
    trials_subj <- data$trials[idx_subj_all]
    trial_vals <- sort(unique(trials_subj))

    for (j in seq_along(trial_vals)) {
      tmp_return_trialwise <- ifelse(j == length(trial_vals) & return_trialwise_parameters, TRUE, FALSE)

      current_trial <- trial_vals[j]
      prefix_rows <- idx_subj_all[trials_subj %in% trial_vals[seq_len(j)]]
      current_rows <- idx_subj_all[trials_subj == current_trial]

      # Rebuild design for the current prefix so that map_p uses updated designs
      dm <- design_model(data[prefix_rows, ],design, model_fun, add_acc = FALSE, compress = FALSE, verbose = FALSE, rt_check = FALSE)

      mask_current <- dm$subjects == subj & dm$trials == current_trial
      if (!any(mask_current)) next

      # Standard mapping + trends + transforms on the prefix
      pm <- map_p(pars, dm, model_list, tmp_return_trialwise)
      tr <- model_list$trend
      if (!is.null(tr)) {
        phases <- vapply(tr, function(x) x$phase, character(1))
        if (any(phases == "pretransform")) pm <- prep_trend_phase(dm, tr, pm, "pretransform", tmp_return_trialwise)
      }
      pm <- do_transform(pm, model_list$transform)
      if (!is.null(tr) && any(phases == "posttransform")) pm <- prep_trend_phase(dm, tr, pm, "posttransform", tmp_return_trialwise)
      cur_dm <- dm[mask_current, , drop = FALSE]
      pr <- model_list$Ttransform(pm[mask_current, , drop = FALSE], cur_dm)
      pr <- add_bound(pr, model_list$bound, cur_dm$lR)

      # Identify current-trial rows inside the prefix design


      # Simulate current trial rows
      if (any(names(dm) == "RACE")) {
        Rrt <- RACE_rfun(cur_dm, pr, model_fun)
      } else {
        Rrt <- model_list$rfun(cur_dm, pr)
      }
      # Write outputs back to original data rows for the current trial
      target_rows <- prefix_rows[mask_current]
      for (nm in dimnames(Rrt)[[2]]) data[target_rows, nm] <- Rrt[, nm]

      # NS I don't actually think this is necessary couldn't this be specified
      # As a standard function in the design?
      # # Optional per-trend feedback → next trial for this subject
      # if(!is.null(tr) && !is.null(tr$feedback_fun)){
      #   nams <- names(tr$feedback_fun)
      #   # Build the window: current prefix plus next-trial rows (if any)
      #   if (j < length(trial_vals)) {
      #     next_rows <- idx_subj_all[trials_subj == trial_vals[j+1]]
      #     window_rows <- c(prefix_rows, next_rows)
      #   } else {
      #     window_rows <- prefix_rows
      #   }
      #   for(i in 1:length(nams)){
      #     fb_vec <- tr$feedback_fun[[i]](data[window_rows,,drop=FALSE])
      #     data[window_rows, nams[i]] <- fb_vec
      #   }
      # }

      # Store trialwise parameters if requested
      # if (!is.null(trialwise_parameters)) trialwise_parameters[target_rows, ] <- pr
      if(tmp_return_trialwise){
        trialwise_parameters <- cbind(pm, attr(pm, "trialwise_parameters"))
      }
    }
  }
  if(is.null(data$lR)) data$lR <- 1
  data <- data[data$lR == unique(data$lR)[1], unique(c(includeColumns, "R", "rt"))]
  data <- data[,!colnames(data) %in% c('lR', 'lM')]
  return(list(data = data, trialwise_parameters = trialwise_parameters))
}


## Functions for working with custom kernels

##' Extract pointers of custom C++ trend kernels from trend list or emc object
##'
##' Extract the pointers so they can be re-added to an emc object after loading it from disk.
##'
##' @param input_data Either an emc object or a trend list.
##' @return A list of custom pointers. The list is of the same size as total number of trends; trends without a custom kernel return NULL.
get_custom_kernel_pointers <- function(input_data) {
  if(inherits(input_data, 'emc')) {
    trend <- input_data[[1]]$model()$trend
  } else {
    trend <- input_data
  }
  if(is.null(trend)) return(NULL)

  ptrs <- lapply(trend, function(x) attr(x, 'custom_ptr'))
  if(all(sapply(ptrs, is.null))) return(NULL)

  return(ptrs)
}

##' (Re-)Set pointers of custom C++ trend kernels to an emc object
##'
##' When an emc object is loaded from disk, or returned by forked processes, the pointers to custom kernels
##' need to be re-created. This is a convenience function to do this.
##'
##' @param emc An emc object
##' @param ptrs A list of pointers, as generated by `get_custom_kernel_pointers()`
##' @return An emc object with the custom pointers re-instated.
set_custom_kernel_pointers <- function(emc, ptrs) {
  if(is.null(ptrs)) return(emc)   # nothing to set
  if(is.null(emc[[1]]$model()$trend)) stop('emc object has no trends, nothing to set...')
  if(length(ptrs) != length(emc[[1]]$model()$trend)) {
    stop('List of potential pointers not equal to number of trends')
  }

  for(chain_ in 1:length(emc)) {
    # update model() function in emc
    if('model' %in% names(emc[[chain_]])) {
      model_list <- emc[[chain_]]$model()
      trend <- model_list$trend
      for(i in 1:length(trend)) {
        if(!is.null(ptrs[[i]])) attr(trend[[i]], 'custom_ptr') <- ptrs[[i]]
      }
      model_list$trend <- trend
      emc[[chain_]]$model <- function() return(model_list)
    }
    ## Update model() function hidden in the design hidden in the prior
    if('prior' %in% names(emc[[chain_]])) {
      for(design_n in 1:length(emc[[chain_]]$prior, 'design')) {
        if('model' %in% names(attr(emc[[chain_]]$prior, 'design')[[design_n]])) {
          model_list <- attr(emc[[chain_]]$prior, 'design')[[design_n]]$model()
          trend <- model_list$trend
          for(i in 1:length(trend)) {
            if(!is.null(ptrs[[i]])) attr(trend[[i]], 'custom_ptr') <- ptrs[[i]]
          }
          model_list$trend <- trend
          attr(emc[[chain_]]$prior, 'design')[[design_n]]$model <- function() return(model_list)
        }
      }
    }
  }
  return(emc)
}


##' Reset pointers of custom C++ trend kernels to an emc object
##'
##' When an emc object is loaded from disk, or returned by forked processes, the pointers to custom kernels
##' need to be re-created. This is a convenience function to do this.
##'
##' @param emc A target emc object with missing pointers
##' @param pointer_source Either a trend object with correct pointers or another emc object with correct pointers
##' @return An emc object with the custom pointers re-instated.
##' @export
fix_custom_kernel_pointers <- function(emc, pointer_source) {
  return(set_custom_kernel_pointers(emc, get_custom_kernel_pointers(pointer_source)))
}

