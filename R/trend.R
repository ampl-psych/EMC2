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
                       par_input = NULL,
                       phase = "premap",
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
    } else{
      trend$kernel <- kernels[i]
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
  bases <- list(
    lin = list(description = "Linear base: parameter + b * k",
               transforms = list(func = list("B0" = "identity")),
               default_pars = "B0"),
    exp_lin = list(description = "Exponential linear base: exp(parameter) + exp(b) * k",
                   transforms = list(func = list("B0" = "exp")),
                   default_pars = "B0"),
    centered = list(description = "Centered mapping: parameter + b*(k - 0.5)",
                    transforms = list(func = list("B0" = "identity")),
                    default_pars = "B0"),
    add = list(description = "Additive base: parameter + k",
               transforms = NULL),
    identity = list(description = "Identity base: k",
                    transforms = NULL)
  )
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
    exp_decr = list(description = "Decreasing exponential kernel: k = exp(-d * c)",
                    transforms = list(func =list("d_ed" = "exp")),
                    default_pars = "d_ed",
                    bases = base_2p),
    exp_incr = list(description = "Increasing exponential kernel: k = 1 - exp(-d * c)",
                    transforms = list(func =list("d_ei" = "exp")),
                    default_pars = "d_ei",
                    bases = base_2p),
    pow_decr = list(description = "Decreasing power kernel: k = (1 + c)^(-d)",
                    transforms = list(func =list("d_pd" = "exp")),
                    default_pars = "d_pd",
                    bases = base_2p),
    pow_incr = list(description = "Increasing power kernel: k = 1 - (1 + c)^(-d)",
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
      "Standard delta rule kernel:",
      "Updates q[i] = q[i-1] + alpha * (c[i-1] - q[i-1]).",
      "Parameters: q0 (initial value), alpha (learning rate)."
    ),
    default_pars = c("q0", "alpha"),
    transforms = list(func = list("q0" = "identity", "alpha" = "pnorm")),
    bases = base_2p),
    delta2 = list(description = paste(
      "Dual kernel delta rule:",
      "Combines fast and slow learning rates",
      "and switches between them based on dSwitch.",
      "Parameters: q0 (initial value), alphaFast (fast learning rate),",
      "propSlow (alphaSlow = propSlow * alphaFast), dSwitch (switch threshold)."
    ),
    default_pars = c("q0", "alphaFast", "propSlow", "dSwitch"),
    transforms = list(func = list("q0" = "identity", "alphaFast" = "pnorm",
                                  "propSlow" = "pnorm", "dSwitch" = "pnorm")),
    bases = base_2p)
  )
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
    cat("          are mapped first, then used to transform cognitive model parameters before their mapping.\n")
    cat("  pretransform: Trend is applied after parameter mapping but before transformations.\n")
    cat("                Cognitive model parameters are mapped first, then trend is applied, followed by transformations.\n")
    cat("  posttransform: Trend is applied after both mapping and transformations.\n")
    cat("                 Cognitive model parameters are mapped and transformed first, then trend is applied.\n")
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

run_kernel <- function(trend_pars = NULL, kernel, input, funptr = NULL) {
  # input: vector or matrix; apply per column and sum contributions; ignore NA contributions
  if (!is.matrix(input)) input <- matrix(input, ncol = 1)
  n <- nrow(input)
  out <- rep(0.0, n)
  if (identical(kernel, "custom")) {
    # trend_pars provided here should already exclude base parameter columns
    if (is.null(funptr)) stop("Missing function pointer for custom kernel. Pass 'funptr'.")
    contrib <- EMC2_call_custom_trend(trend_pars, input, funptr)
    contrib[is.na(contrib)] <- 0
    return(contrib)
  }
  for (j in seq_len(ncol(input))) {
    covariate <- input[, j]
    contrib <- switch(kernel,
      lin_decr = -covariate,
      lin_incr = covariate,
      exp_decr = exp(-trend_pars[,1]*covariate),
      exp_incr = 1-exp(-trend_pars[,1]*covariate),
      pow_decr = (1+covariate)^(-trend_pars[,1]),
      pow_incr = (1-(1+covariate)^(-trend_pars[,1])),
      poly2 =  trend_pars[,1]*covariate + trend_pars[,2]*covariate^2,
      poly3 =  trend_pars[,1]*covariate + trend_pars[,2]*covariate^2 + trend_pars[,3]*covariate^3,
      poly4 =  trend_pars[,1]*covariate + trend_pars[,2]*covariate^2 + trend_pars[,3]*covariate^3 + trend_pars[,4]*covariate^4,
      delta = {
        good <- !is.na(covariate)
        tmp <- run_delta(trend_pars[good,1], trend_pars[good,2], covariate[good])
        z <- rep(0.0, n); z[good] <- tmp; z
      },
      delta2 = {
        good <- !is.na(covariate)
        tmp <- run_delta2(trend_pars[good,1], trend_pars[good,2], trend_pars[good,3], trend_pars[good,4], covariate[good])
        z <- rep(0.0, n); z[good] <- tmp; z
      }
    )
    contrib[is.na(contrib)] <- 0
    out <- out + contrib
  }
  out
}

prep_trend_phase <- function(dadm, trend, pars, phase){
  # Apply only trends in the requested phase, sequentially
  tnames <- names(trend)
  all_remove <- character(0)
  for (idx in seq_along(trend)){
    cur_trend <- trend[[idx]]
    if (!identical(cur_trend$phase, phase)) next
    par <- tnames[idx]
    all_remove <- c(all_remove, cur_trend$trend_pnames)
    pars[, par] <- run_trend(dadm, cur_trend, pars[, par], pars[, cur_trend$trend_pnames, drop = FALSE], pars)
  }
  if (length(all_remove)) pars <- pars[, !(colnames(pars) %in% unique(all_remove)), drop = FALSE]
  return(pars)
}

run_trend <- function(dadm, trend, param, trend_pars, pars_full = NULL){
  n_base_pars <- switch(trend$base,
                        lin = 1,
                        exp_lin = 1,
                        centered = 1,
                        add = 0,
                        identity = 0)
  out <- numeric(nrow(dadm))
  # Build combined input matrix: covariates + par_input
  cov_cols <- trend$covariate
  par_in_cols <- if (!is.null(trend$par_input)) trend$par_input else character(0)
  n_inputs <- length(cov_cols) + length(par_in_cols)
  if (n_inputs > 0) {
    input_all <- matrix(NA_real_, nrow(dadm), n_inputs)
    if (length(cov_cols) > 0) {
      for (i in seq_along(cov_cols)) input_all[, i] <- dadm[, cov_cols[i]]
    }
    if (length(par_in_cols) > 0) {
      for (j in seq_along(par_in_cols)) input_all[, length(cov_cols) + j] <- pars_full[, par_in_cols[j]]
    }
    if (ncol(trend_pars) > n_base_pars) {
      kernel_pars <- trend_pars[, seq.int(n_base_pars + 1, ncol(trend_pars)), drop = FALSE]
    } else {
      kernel_pars <- trend_pars[, 0, drop = FALSE]
    }
    funptr <- if (identical(trend$kernel, "custom")) attr(trend, "custom_ptr") else NULL
    out <- out + run_kernel(kernel_pars, trend$kernel, input_all, funptr = funptr)
  }
  # Do the mapping
  out <- switch(trend$base,
                lin = param + trend_pars[,1]*out,
                exp_lin = exp(param) + trend_pars[,1]*out,
                centered = param + trend_pars[,1]*(out-.5),
                add = param + out,
                identity = out
  )
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
  # Premap + par_input not supported (no full parameter matrix yet)
  for (i in seq_along(trend)) {
    if (identical(trend[[i]]$phase, "premap") && length(trend[[i]]$par_input) > 0) {
      warning("par_input is ignored for premap trends in current implementation")
      break
    }
  }
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



