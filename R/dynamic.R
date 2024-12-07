#' Create a trend specification for dynamic parameters
#'
#' @param par_names Character vector specifying which parameters should be made dynamic
#' @param cov_names Character vector specifying which covariates to use for each trend
#' @param kernels Character vector specifying which kernel function to use for each trend
#' @param bases Optional character vector specifying which base function to use for each trend
#' @param shared Logical indicating if parameters should be shared across trends on different parameters
#' @param trend_pnames Optional character vector specifying custom parameter names
#' @param premap Logical indicating if trend should be applied before or after parameter mapping
#' @param pretransform If !premap, logical indicating if trend should be applied before or after parameter transformation
#' @return A list containing the trend specifications for each parameter
#' @export
#'
#' @examples
#' # Make B and v parameters dynamic using exponential kernels
#' dynamic <- make_trend(
#'   par_names = c("B", "v"),
#'   cov_names = "strial",
#'   kernels = c("exp_incr", "poly3"),
#'   premap = TRUE
#' )
#'
make_trend <- function(par_names, cov_names, kernels, bases = NULL,
                       shared = FALSE, trend_pnames = NULL, premap = TRUE,
                       pretransform = FALSE){
  if(pretransform & premap){
    warning("Setting pretransform has no effect if premap = TRUE")
  }
  if(!(length(par_names) == length(kernels))){
    stop("Make sure that par_names and kernels have the same length")
  }
  if(length(cov_names) != length(par_names) && length(cov_names) == 1){
    cov_names <- rep(cov_names, length(par_names))
  } else{
    stop("Make sure that cov_names and par_names have the same length")
  }
  if(!is.null(bases)){
    if(length(kernels) != length(bases)){
      stop("If bases not is NULL make sure you specify as many bases as kernels")
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
    if(is.null(bases)){
      trend$base <- trend_help(kernels[i], do_return = TRUE)$bases[1]
    } else{
      if(bases[i] %in% names(trend_help(kernels[i], do_return = TRUE)$bases)){
        stop("base type not supported with kernel, see `trend_help(<kernel>)`")
      }
      trend$base <- bases[i]
    }
    # add par names
    trend_pnames <- trend_help(base = trend$base, do_return = TRUE)$default_pars
    trend_pnames <- c(trend_pnames, trend_help(kernel = kernels[i], do_return = TRUE)$default_pars)
    trend$trend_pnames <- paste0(par_names[i], ".", trend_pnames)
    trend$covariate <- cov_names[i]
    trend$shared <- shared
    trends_out[[par_names[i]]] <- trend
  }
  attr(trends_out, "premap") <- premap
  attr(trends_out, "pretransform") <- pretransform & !premap
  attr(trends_out, "posttransform") <- !pretransform & !premap
  return(trends_out)
}


#' Get parameter names from trend object
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
  return(out)
}

trend_help <- function(kernel = NULL, base = NULL, ...){
  dots <- add_defaults(list(...), do_return = FALSE, return_types = FALSE)
  bases <- list(
    lin = list(description = "Linear base: parameter + b * k",
               transforms = list(func = list("B0" = "identity")),
               default_pars = "B0"),
    exp_lin = list(description = "Exponential linear base: exp(parameter) + exp(b) * k",
                   transforms = list(func = list("B0" = "exp")),
                   default_pars = "B0"),
    add = list(description = "Additive base: parameter + k",
               transforms = NULL),
    identity = list(description = "Identity base: k",
                    transforms = NULL)
  )
  base_2p <- names(bases)[1:2]
  base_1p <- names(bases)[3:4]
  kernels <- list(
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
      transforms = list(func = list("q0" = "identity", "alpha" = "exp")),
      bases = base_2p),
    delta2 = list(description = paste(
      "Dual kernel delta rule:",
      "Combines fast and slow learning rates",
      "and switches between them based on dSwitch.",
      "Parameters: q0 (initial value), alphaFast (fast learning rate),",
      "propSlow (alphaSlow = propSlow * alphaFast), dSwitch (switch threshold)."
    ),
    default_pars = c("q0", "alphaFast", "propSlow", "dSwitch"),
    transforms = list(func = list("q0" = "identity", "alphaFast" = "exp",
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

#' Run a kernel function for trend analysis
#'
#' @param param Parameter value to be modified by the trend
#' @param trend_pars Matrix or vector of trend parameters
#' @param kernel String specifying which kernel function to use. Options include:
#'   \itemize{
#'     \item lin_decr - Linear decreasing
#'     \item lin_incr - Linear increasing
#'     \item exp_decr - Exponential decreasing
#'     \item exp_incr - Exponential increasing
#'     \item pow_decr - Power decreasing
#'     \item pow_incr - Power increasing
#'     \item poly2 - Second order polynomial
#'     \item poly3 - Third order polynomial
#'     \item poly4 - Fourth order polynomial
#'     \item delta - Single kernel learning rule
#'     \item delta2 - Dual kernel learning rule
#'   }
#' @param covariate Vector of covariate values
#'
#' @return Vector of transformed values based on the specified kernel function
run_kernel <- function(param, trend_pars = NULL, kernel, covariate) {
  # Covariate
  out <- switch(kernel, # kernel
                  # Linear decreasing, linear base type
                  lin_decr = -covariate,
                  # Linear increasing, linear base type
                  lin_incr = covariate,
                  # Exponential decreasing, linear base type
                  exp_decr = exp(trend_pars[,1]*covariate),
                  # Exponential increasing - help
                  exp_incr = (1-exp(-trend_pars[1]*covariate)),
                  # Decreasing power
                  pow_decr = (1+covariate)^(-trend_pars[1]),
                  # Increasing power
                  pow_incr = (1-(1+covariate)^(-trend_pars[1])),
                  # Second to fourth order polynomials, for these base is default none? Or linear?
                  poly2 =  trend_pars[,1]*covariate + trend_pars[,2]*covariate^2,
                  poly3 = trend_pars[,1]*covariate + trend_pars[,2]*covariate^2 + trend_pars[,3]*covariate^3,
                  poly4 = trend_pars[,1]*covariate + trend_pars[,2]*covariate^2 + trend_pars[,3]*covariate^3 + trend_pars[,4]*covariate^4,
                  # Single and dual kernel learning rule
                  delta = run_delta(trend_pars[,1],trend_pars[,2],
                                covariate),
                  delta2 = run_delta2(trend_pars[,1],trend_pars[,2], trend_pars[,3], trend_pars[,4],
                                  covariate),
  )
  return(out)
}

prep_trend <- function(dadm, trend, pars){
  for(par in names(trend)){
    pars[,par] <- run_trend(dadm, trend, pars[,par], pars[,trend[[par]]$trend_pnames, drop = FALSE])
  }
  pars[,!colnames(pars) %in% get_trend_pnames(trend)]
  return(pars)
}

run_trend <- function(dadm, trend, param, trend_pars){
  # For now assume no sequential process

  n_base_pars <- switch(trend$base,
                        lin = 1,
                        exp_lin = 1,
                        add = 0,
                        identity = 0)
  out <- numeric(nrow(dadm))
  for(cov_name in trend$covariate){
    # Loop over covariates, first filter out NAs
    NA_idx <- is.na(dadm[,cov_name])
    cov_tmp <- dadm[!NA_idx, cov_name]
    param_tmp <- param[!NA_idx]
    trend_pars_tmp <- trend_pars[!NA_idx,, drop = FALSE]
    if(trend$kernel %in% c("delta", "delta2")){
      # These trends have a sequential nature, don't filter duplicate entries
      filter <- rep(T, length(cov_tmp))
    } else{
      # Else filter out duplicated entries
      together <- cbind(cov_tmp, trend_pars_tmp)
      filter <- !duplicated(together)
    }
    # Prep trend parameters to filter out base parameters
    kernel_pars <- trend_pars_tmp[filter,(n_base_pars+1):ncol(trend_pars), drop = FALSE]
    # Keep an expansion index
    unq_idx <- cumsum(filter)
    # Run trend
    output <- run_kernel(param_tmp[filter], kernel_pars, trend$kernel, cov_tmp[filter])
    # Decompress output and map back based on non-NA
    out[!NA_idx] <- out[!NA_idx] + output[unq_idx]
  }
  # Do the mapping
  out <- switch(trend$base,
                        lin = param + trend_pars[,1]*out,
                        exp_lin = exp(param) + trend_pars[,1]*out,
                        add = param + out,
                        identity = out
  )
  return(out)
}

check_trend <- function(trend, covariates = NULL, model = NULL, formula = NULL) {
  if(!is.null(model)){
    if(!attr(trend, "premap")){
      param_names <- names(trend)
      if (!all(param_names %in% names(model()$p_types))){
        stop("trend premap has a parameter type name not in the model")
      }
    }
  }
  if (is.null(covariates)) stop("must specify covariates when using trend")
  covnames <- unlist(lapply(trend,function(x)x$covnames))
  if (!all(covnames %in% names(covariates))){
    stop("trend has covnames not in covariates")
  }
  if (any(duplicated(names(trend)))) stop("Duplicate parameter names in trend")
  for(par in names(trend)){
    cur_dyn <- trend[[par]]
    if (!is.null(formula)) {
      isin <-  cur_dyn$trend_pnames %in% unlist(lapply(formula,function(x)all.vars(x)[1]))
      if(any(!isin)){
      # Add missing trend parameters to formula with intercept-only model
        formula <- c(formula, lapply(cur_dyn$trend_pnames[!isin], function(x) as.formula(paste(x, "~ 1"))))
        message("Intercept formula added for trend_pars: ", paste(cur_dyn$trend_pnames[!isin],collapse=", "))
      }
    }
    trend[[par]] <- cur_dyn
  }
  return(formula)
}


update_model_trend <- function(trend, model) {
  # Get model list to modify
  model_list <- model()

  # For each parameter in the trend
  for (par in names(trend)) {
    cur_trend <- trend[[par]]

    # Get default transforms from kernel and base
    kernel_transforms <- trend_help(cur_trend$kernel, do_return = TRUE)$transforms
    base_transforms <- trend_help(base = cur_trend$base, do_return = TRUE)$transforms

    # Combine transforms
    if (!is.null(kernel_transforms) || !is.null(base_transforms)) {
      tmp <- c(kernel_transforms$func, base_transforms$func)
      names(tmp) <- cur_trend$trend_pnames
      # Update the appropriate transform list based on premap
      model_list$transform$func <- c(model_list$transform$func, unlist(tmp))
    }
    model_list$p_types <- c(model_list$p_types, setNames(numeric(length(cur_trend$trend_pnames)), cur_trend$trend_pnames))
  }
  model_list$trend <- trend
  # Return updated model function
  model <- function() { return(model_list) }
  return(model)
}

run_delta <- function(q0,alpha,covariate) {
  q <- pe <- numeric(length(target))
  q[1] <- q0[1]
  for (i in 2:length(q)) {
    pe[i-1] <- covariate[i-1]-q[i-1]
    q[i] <- q[i-1] + alpha[i-1]*pe[i-1]
  }
  return(q)
}

run_delta2_i <- function(q0,alphaFast,propSlow,dSwitch,covariate) {
  q <- qFast <- qSlow <- peFast <- peSlow <- numeric(length(target))
  q[1] <- qFast[1] <- qSlow[1] <- q0[1]
  alphaSlow <- propSlow*alphaFast
  for (i in 2:length(target)) {
    peFast[i-1] <- target[i-1]-qFast[i-1]
    peSlow[i-1] <- target[i-1]-qSlow[i-1]
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



