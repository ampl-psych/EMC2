
# Rewrite of the factor.switching package ---------------------------------

# sp_new <- function(iter, lambda_varimax, q, p, dim_all_c, all_c, lambda_hat, st,
#                    cost.matrix, perm){
#   lambda <- matrix(lambda_varimax[,,iter], nrow = p, ncol = q)
#   costs <- numeric(dim_all_c)
#   for (i in 1:dim_all_c) {
#     c_vec <- all_c[i, ]
#     lambda_switch <- matrix(c_vec, ncol = q, nrow = p,
#                             byrow = T) * lambda_hat
#     for (j in 1:q) {
#       temp <- (lambda - lambda_switch[, j])^2
#       cost.matrix[j, ] <- colSums(temp)
#     }
#     matr <- lpSolve::lp.assign(cost.matrix)$solution
#     perm[i,] <- order(which(matr > 0, arr.ind = T)[,1])
#     costs[i] <- sum(cost.matrix * matr)
#   }
#   minIndex <- order(costs)[1]
#   switchedMatrix <- matrix(all_c[minIndex, ], ncol = q, nrow = p, byrow = T) * lambda[,perm[minIndex, ]]
#   return(list(min_perm = perm[minIndex, ], min_c = all_c[minIndex, ], min_cost =  costs[minIndex],
#               switchedMatrix = switchedMatrix))
# }


#' Reorder MCMC Samples of Factor Loadings
#'
#' This function reorders MCMC samples of factor loadings to address the label switching problem
#' in Bayesian factor analysis. It implements a parallelized version of the code and
#' algorithm proposed by Papastamoulis and Ntzoufras (2022)
#' @references Papastamoulis, P., & Ntzoufras, I. (2022). On the identifiability of Bayesian factor
#' analytic models. *Statistical Computing*, 32(2), 1-29. doi: 10.1007/s11222-022-10084-4
#'
#' @param emc an 'emc' object of type `infnt_factor`.
#' @param lambda Needs to be supplied if emc is not supplied.
#' Array of factor loadings with dimensions p (variables) x q (factors) x n (MCMC iterations)
#' @param verbose Logical; whether to print progress information
#' @param n_cores Number of cores for parallel processing
#' @param rotate_fun A function that returns an orthogonally rotated factor loadings matrix. If NULL uses `varimax`
#' @return A list containing:
#'   \item{lambda_reordered}{Array of reordered loadings}
#'   \item{lambda_reordered_mcmc}{Array of reordered loadings as MCMC object}
#'   \item{lambda_hat}{Matrix of mean loadings after reordering}
#'   \item{v_vectors}{Matrix of permutation vectors}
#'   \item{c_vectors}{Matrix of sign-switching vectors}
#'
#' @examples
#' # This function works natively with emc objects, but also factor arrays:
#' # Simulate a small example with 5 variables, 2 factors, and 10 MCMC iterations
#' set.seed(123)
#' p <- 5  # Number of variables
#' q <- 2  # Number of factors
#' n <- 10 # Number of MCMC iterations
#'
#' # Create random factor loadings with label switching
#' lambda <- array(0, dim = c(p, q, n))
#' for (i in 1:n) {
#'   # Generate base loadings
#'   base_loadings <- matrix(rnorm(p*q, 0, 0.5), p, q)
#'   base_loadings[1:3, 1] <- abs(base_loadings[1:3, 1]) + 0.5  # Strong loadings on factor 1
#'   base_loadings[4:5, 2] <- abs(base_loadings[4:5, 2]) + 0.5  # Strong loadings on factor 2
#'
#'   # Randomly switch labels and signs
#'   if (runif(1) > 0.5) {
#'     # Switch factor order
#'     base_loadings <- base_loadings[, c(2, 1)]
#'   }
#'   if (runif(1) > 0.5) {
#'     # Switch sign of factor 1
#'     base_loadings[, 1] <- -base_loadings[, 1]
#'   }
#'   if (runif(1) > 0.5) {
#'     # Switch sign of factor 2
#'     base_loadings[, 2] <- -base_loadings[, 2]
#'   }
#'
#'   lambda[,,i] <- base_loadings
#' }
#'
#' # Align the loadings
#' result <- align_loadings(lambda = lambda, verbose = TRUE, n_cores = 1)
#'
#' # Examine the aligned loadings
#' print(result)
#'
#' @export
align_loadings <- function (emc = NULL, lambda = NULL, n_cores = 1, verbose = TRUE, rotate_fun = NULL)
{
  maxIter <- 100; threshold <- 1e-06;
  rotate <- TRUE; printIter <- 1000;
  if(is.null(lambda) & is.null(emc)) stop("Need to supply either emc or lambda")
  if(is.null(lambda)){
    lambda <- get_pars(emc, selection = "loadings", merge_chains = TRUE, return_mcmc = FALSE)
  }
  energy <- apply(lambda^2, 2, mean)          # length q
  lambda <- lambda[,energy / max(energy) > 0.0075,, drop = F]
  q <- ncol(lambda)
  p <- nrow(lambda)
  mcmcIterations <- dim(lambda)[3]
  threshold <- threshold * mcmcIterations * p * q
  lambda_varimax <- lambda
  if(is.null(rotate_fun)){
    rotate_fun <- function(loadings){
      varimax(loadings, normalize = F)$loadings
    }
  }
  if (rotate) {
    if (q > 1) {
      for (iter in 1:mcmcIterations) {
        lambda_varimax[,,iter] <- rotate_fun(lambda[,,iter])
      }
    }
  }
  all_c <- array(data = 1, dim = c(2^q, q))
  l <- 1
  if (q > 1) {
    for (i in 1:(q - 1)) {
      pos <- combn(1:q, i)
      j <- dim(pos)[2]
      for (k in 1:j) {
        l <- l + 1
        all_c[l, pos[, k]] <- rep(-1, i)
      }
    }
  }
  all_c[l + 1, ] <- rep(-1, q)
  lambda_hat <- apply(lambda_varimax, 1:2, mean)
  lambda_hat_zero <- lambda_hat
  st <- 1:q
  dim_all_c <- 2^q
  dim_all_v <- factorial(q)
  perm <- matrix(1:q, ncol = q, nrow = dim_all_c, byrow = TRUE)
  costs <- numeric(dim_all_c)
  f <- numeric(maxIter)
  lambda_hat_values <- array(data = NA, dim = c(p, q, maxIter))
  totalIterations <- 0
  criterion = TRUE
  cost.matrix <- matrix(numeric(q * q), nrow = q, ncol = q)
  t1 <- numeric(maxIter)
  start_time <- Sys.time()
  while ((criterion == TRUE) & (totalIterations < maxIter)) {
    totalIterations <- totalIterations + 1
    if(verbose) cat(paste0("* iteration: ", totalIterations), "\n")
    lambda_hat_new <- 0 * lambda_hat
    sp_res <- parallel::mclapply(X = 1:mcmcIterations, FUN = sp_new, lambda_varimax, q, p, dim_all_c, all_c, lambda_hat,
                                 st, cost.matrix, perm, mc.cores = n_cores)
    min_costs <- sapply(sp_res, FUN = function(x) return(x$min_cost))
    v_vectors <- do.call(rbind, lapply(sp_res, FUN = function(x) return(x$min_perm)))
    c_vectors <- do.call(rbind, lapply(sp_res, FUN = function(x) return(x$min_c)))
    objective <- sum(min_costs)
    lambda_hat_new <- Reduce('+', lapply(sp_res, FUN = function(x) return(x$switchedMatrix)))
    f[totalIterations] <- objective
    if (totalIterations > 1) {
      if (f[totalIterations - 1] - f[totalIterations] <
          threshold) {
        criterion = FALSE
      }
    }
    if (verbose) {
      cat(paste0("   -----  objective function = ", round(objective,
                                                          3), " -----"), "\n")
      cat("\n")
    }
    lambda_hat_new <- lambda_hat_new/mcmcIterations
    lambda_hat <- lambda_hat_new
    lambda_hat_values[, , totalIterations] <- lambda_hat
    end_time <- Sys.time()
    t1[totalIterations] <- as.numeric(difftime(end_time,
                                               start_time, units = "min"))
  }
  c_vec <- rep(1, q)
  v_vec <- 1:q
  f_zero <- 0
  for (i in 1:mcmcIterations) {
    lambda <- matrix(lambda_varimax[,,i], nrow = p, ncol = q)
    switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow = T) *
      lambda[, v_vec]
    f_zero <- f_zero + sum((switchedMatrix - lambda_hat)^2)
  }
  t_exact <- c(0, t1[1:totalIterations])
  f_exact <- c(f_zero, f[1:totalIterations])
  objective_function <- data.frame(time = t_exact, value = f_exact)
  lambda_reordered <- lambda_varimax
  for (i in 1:mcmcIterations) {
    lambda <- matrix(lambda_varimax[,,i], nrow = p, ncol = q)
    switchedMatrix <- matrix(c_vectors[i, ], ncol = q, nrow = p,
                             byrow = T) * lambda[, v_vectors[i, ]]
    lambda_reordered[,,i] <- switchedMatrix
  }
  lambda_reordered <- rearrange_loadings(lambda_reordered)
  if(!is.null(emc)){
    emc <- subset(emc, stage = get_last_stage(emc))
    n_iter <- chain_n(emc)[1,get_last_stage(emc)]
    for(i in 1:length(emc)){
      emc[[i]]$samples$lambda <- lambda_reordered[,,((i-1)*(n_iter)+ 1):(i*n_iter), drop = FALSE]
    }
    class(emc) <- "emc"
    return(emc)
  }
  return(lambda_reordered)
}



# Functions useful for descriptives of the loadings -----------------------
rearrange_loadings <- function(lambda, metric = c("ssq", "absmean"))
{
  metric <- match.arg(metric)
  # ----- collapse MCMC draws to a single p × q summary -----------------
  if (length(dim(lambda)) == 3L) {
    lambda_hat <- apply(lambda, 1:2, mean)           # p × q
  } else {
    lambda_hat <- lambda
    lambda      <- array(lambda, dim = c(dim(lambda), 1L))  # unify shape
  }

  # ----- dominance score for each factor --------------------------------
  score <- switch(metric,
                  ssq     = colSums(lambda_hat^2),
                  absmean = colMeans(abs(lambda_hat)))

  ord <- order(score, decreasing = TRUE)             # strongest → first

  # ----- sign-flip so column mean ≥ 0 -----------------------------------
  sign_flip <- sign(colMeans(lambda_hat[, ord, drop = F]))
  sign_flip[sign_flip == 0] <- 1                     # treat exact zeros as +

  # ----- apply permutation + sign flips over *all* iterations -----------
  lambda_re <- lambda[, ord, , drop = FALSE]         # permute columns
  for (j in seq_along(sign_flip)) {
    lambda_re[, j, ] <- sign_flip[j] * lambda_re[, j, ]
  }
  return(lambda_re)
}


standardize_loadings_SEM <- function(loadings, variances, psi_inv){
  new_loadings <- loadings
  n_pars <- nrow(loadings)
  for(i in 1:dim(loadings)[3]){
    Dinv <- 1/sqrt(diag(variances[,,i]))
    psi <- solve(psi_inv[,,i])
    Deta   <-  sqrt(diag(psi))     # length‑k vector
    new_loadings[,,i] <- (loadings[,,i]*Dinv)*rep(Deta, each = n_pars)
  }
  return(new_loadings)
}


standardize_loadings <- function(loadings, residuals){
  stdize_set <- function(samples = NULL, idx = NULL, loadings = NULL, residuals = NULL){
    if(is.null(loadings)) loadings <- samples$samples$lambda[,,idx, drop = F]
    if(is.null(residuals)) residuals <- samples$samples$epsilon_inv[,idx]
    new_loadings <- loadings
    for(i in 1:dim(loadings)[3]){
      sigma_err <- residuals[,i]
      vars <- loadings[,,i] %*% t(loadings[,,i]) + diag(sigma_err)
      Dinv <- 1/sqrt(diag(vars))
      new_loadings[,,i] <- loadings[,,i]*Dinv
    }
    return(new_loadings)
  }
  out <- stdize_set(loadings = loadings, residuals = residuals)
  if(is.null(colnames(out))) colnames(out) <- paste0("F", 1:ncol(out))
  return(out)
}

#' Cut Factors Based on Credible Loadings
#'
#' This function removes factors that do not have more than one credible loading
#' based on the specified confidence interval.
#'
#' @param emc An 'emc' object containing factor analysis results
#' @param CI Numeric. Confidence interval percentage (default is 95)
#'
#' @return An 'emc' object with factors that don't meet the credibility criterion removed
#'
#' @export
cut_factors <- function(emc, CI = 95){
  lower <- (100-CI)/2
  upper <- (100-CI)/2 + CI
  credints <- credint(emc, selection = "std_loadings", probs = c(lower/100, .5, upper/100))
  is_credible <- get_credible_factors(credints)
  if(!any(is_credible)) stop("no factors with more than 1 credible loading")
  for(i in 1:length(emc)){
    emc[[i]]$samples$lambda <- emc[[i]]$samples$lambda[,is_credible,]
    emc[[i]]$samples$eta <- emc[[i]]$samples$eta[,is_credible,]
  }
  class(emc) <- "emc"
  return(emc)
}

get_credible_factors <- function(credints){
  sapply(credints, function(x) (sum(x[,1] > 0) + sum(x[,3] < 0)) > 1)
}

#' Plot Group-Level Relations
#'
#' An adjusted version of the `corrplot` package function `corrplot()` tailored
#' to `EMC2` and the plotting of estimated correlations.
#'
#' @param emc An EMC2 object, commonly the output of `run_emc()`.
#' @param stage Character. The stage from which to take the samples, defaults to
#' the sampling stage `sample`.
#' @param plot_cred Boolean. Whether to plot the 95 percent credible intervals or not
#' @param plot_means Boolean. Whether to plot the means or not
#' @param only_cred Boolean. Whether to only plot credible values
#' @param nice_names Character string. Alternative names to give the parameters
#' @param selection Character. Whether to plot correlations or loadings
#' @param use_par Character. Which parameters to include. If null, includes all.
#' @param ... Optional additional arguments
#' @return No return value, creates a plot of group-level relations
#' @examples
#' # For a given set of hierarchical model samples we can make a
#' # correlation matrix plot.
#' plot_relations(samples_LNR, only_cred = TRUE, plot_cred = TRUE)
#' # We can also only plot the correlations where the credible interval does not include zero
#' plot_relations(samples_LNR, plot_means = TRUE, only_cred = TRUE)
#'
#' @export

plot_relations <- function(emc = NULL, stage = "sample",  plot_cred = FALSE,
                           plot_means = TRUE, only_cred = TRUE, nice_names = NULL,
                           selection = "correlation", use_par = NULL, ...){

  dots <- add_defaults(list(...), probs = c(0.025, .975))
  if(!selection %in% c("correlation", "loadings", "std_loadings")){
    stop("Selection must be either 'correlation', 'std_loadings', or 'loadings'")
  }
  do_corr <- FALSE
  if(selection == "correlation"){
    do_corr <- TRUE
  }
  addCoef.col <- "black"
  values <- get_pars(emc, selection = selection, return_mcmc = F, merge_chains = TRUE, remove_constants = F)
  if(!is.null(use_par)){
    if(any(!use_par %in% rownames(values))) stop("some parameter names in use_par not in estimated parameters")
    if(selection == "correlation"){
      values <- values[use_par, use_par,, drop = F]
    } else{
      values <- values[use_par,,, drop = F]
    }
  }
  means <- apply(values, 1:2, mean)



  # You might have to play around with the x and y limits of the legend.

  # Only add this if you want credible intervals
  if(plot_cred || only_cred){
    cred <- aperm(apply(values, 1:2, quantile, probs = dots$probs))
    conf <- paste0("[", format(cred[,,1,drop = F], digits=1), ":", format(cred[,,2,drop = F], digits=1), "]")
    if(only_cred){
      is_cred <- unique(c(which(cred[,,1] > 0), which(cred[,,2] < 0)))
      conf[!1:length(conf) %in% is_cred] <- ""
      is_cred <- unique(c(which(t(cred[,,1] > 0)), which(t(cred[,,2] < 0))))
      means[!1:length(means) %in% is_cred] <- NA
    }
    cred <- round(cred, 2)
  }
  if(do_corr){
    if(!is.null(nice_names)) colnames(means) <- rownames(means) <- nice_names

    if(only_cred){
      p.mat <- means
      p.mat[is.na(means)] <- 1
      p.mat[!is.na(means)] <- 0
      means[is.na(means)] <- 0
      corrplot(means, addCoef.col =addCoef.col, number.cex = .75, tl.col = "black",
               p.mat = p.mat, insig = "blank", sig.level = 0.05)
    } else{
      corrplot(means, addCoef.col =addCoef.col, number.cex = .75, tl.col = "black")
    }

  } else{
    cols <- diverging_hcl(200, palette = "Red-Green")
    if(!is.null(nice_names)) rownames(means) <- nice_names

    # You might have to play around with the x and y limits of the legend.
    if(only_cred){
      p.mat <- means
      p.mat[is.na(means)] <- 1
      p.mat[!is.na(means)] <- 0
      means[is.na(means)] <- 0
      colnames(p.mat) <- colnames(means) <- 1:ncol(means)
      p <- corrplot(means, is.corr = F, col = cols, col.lim = c(-1, 1),
                    cl.pos = "n", addCoef.col =addCoef.col, number.cex = .75, tl.col = "black",
                    p.mat = p.mat, insig = "blank", sig.level = 0.05)
    } else{
      p <- corrplot(means, is.corr = F, col = cols, col.lim = c(-1, 1),
                    cl.pos = "n", addCoef.col =addCoef.col, number.cex = .75, tl.col = "black")
    }
    max_x <- max(p$corrPos$x)
    max_y <- max(p$corrPos$y)
    colorlegend(xlim=c(max_x + 1, max_x + 3), ylim=c(max_y/2 - max_y/5,max_y/2 + max_y/5), cols, c(seq(-1,1,.25)), align="l", vertical=TRUE, addlabels=TRUE)
  }

  if(plot_cred){
    xs <- row(matrix(cred[,,1], nrow = ncol(means)))
    ys <- (ncol(matrix(cred[,,1], nrow = ncol(means)))+1) - col(matrix(cred[,,2], nrow = ncol(means))) - 0.05
    conf[conf == "[ 1.0:1.0]"] <- ""
    conf[conf == "[ 1.00:1.0]"] <- ""
    if(plot_means){
      ys <- ys - 0.1
      text(xs, ys, conf, pos=1, cex=0.55, font = 2)
    } else{
      text(xs, ys, conf, pos=1, cex=0.7, font = 2)
    }
  }
}

#'
#' Factor diagram plot
#' #Makes a factor diagram plot. Heavily based on the fa.diagram function of the `psych` package.
#'
#' @param emc An emc object
#' @param stage Character. The stage from which to take the samples
#' @param loadings An array of loadings. Can be alternatively supplied if emc is not supplied
#' @param standardize Boolean. Whether to standardize the loadings
#' @param simple Boolean. Whether the factor diagram should be simplified for visual clarity.
#' @param only_cred Boolean. Whether to only plot the credible loadings
#' @param cut Numeric. Mean loadings beneath this number will be excluded.
#' @param nice_names Character vector. Alternative names to give the parameters
#' @param factor_names Character vector. Names to give the different factors
#' @param sort Boolean. Whether to sort the paramaters before plotting for visual clarity.
#' @param adj Integer. Adjust to adjust loading values positions in the diagram if illegible.
#' @param main Character vector. Title of the plot
#' @param cex Integer. Font size
#' @return NULL

factor_diagram <- function(emc = NULL, stage = "sample",
                                loadings = NULL, standardize = TRUE,
                                simple = FALSE, only_cred = TRUE,
                                cut = 0, nice_names = NULL,
                                factor_names = NULL, sort = TRUE,
                                adj = 1, main = NULL, cex = NULL){
  if(standardize){
    loadings <- get_pars(emc, selection = "std_loadings", return_mcmc = F, merge_chains = TRUE)
  } else{
    loadings <- get_pars(emc, selection = "loadings", return_mcmc = F, merge_chains = TRUE)
  }
  means <- apply(loadings, 1:2, mean)
  if(only_cred){
    cred <- aperm(apply(loadings, 1:2, quantile, probs = c(0.025, 0.975)))
    is_cred <- unique(c(which(t(cred[,,1] > 0)), which(t(cred[,,2] < 0))))
    means[!1:length(means) %in% is_cred] <- 0
    if(cut == 0) cut <- 0.0001
  }
  if(!is.null(nice_names)){
    rownames(means) <- nice_names
  }
  if(!is.null(factor_names)){
    colnames(means) <- factor_names
  }
  fa.diagram(means, cut = cut, simple = simple, sort = sort, adj = adj,
             main = main, cex = cex)
}


is_double_decimals <- function(x) {
  decimals <- abs(x - round(x, 2))
  return(decimals > .Machine$double.eps^0.5)
}


.extract_group_design_info <- function(group_design, par_names) {
  random_blocks <- stats::setNames(vector("list", length(par_names)), par_names)

  for (par in par_names) {
    random <- NULL
    map <- NULL

    if (!is.null(group_design) && !is.null(group_design[[par]])) {
      obj <- group_design[[par]]
      if (is.list(obj) && !is.data.frame(obj)) {
        random <- obj$random
        map <- obj$map
      }
    }

    blocks <- list()
    if (is.list(random) && length(random) > 0) {
      for (j in seq_along(random)) {
        z <- random[[j]]
        if (is.null(z) || is.null(dim(z))) {
          next
        }
        gname <- names(random)[j]
        if ((is.null(gname) || gname == "") && !is.null(map) && length(map) >= j) {
          gname <- map[[j]]$gname
        }
        if (is.null(gname) || gname == "") {
          gname <- paste0("random", j)
        }

        cols_random <- colnames(z)
        if (is.null(cols_random) || length(cols_random) == 0) {
          cols_random <- paste0(gname, "_", seq_len(ncol(z)))
        }

        blocks[[length(blocks) + 1]] <- list(
          name = paste0(par, "_", gname),
          gname = gname,
          n_cols = ncol(z),
          cols = cols_random
        )
      }
    }
    random_blocks[[par]] <- blocks
  }

  return(random_blocks)
}


.dot_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  return(x)
}


.edge_line <- function(from, to, label = NULL) {
  if (is.null(label) || !nzchar(label)) {
    return(paste0("  ", from, " -> ", to, ";"))
  }
  return(paste0("  ", from, " -> ", to, " [label=\"", .dot_escape(label), "\"];"))
}


#' Graphical Model
#'
#' Draws a probabilistic graphical model (PGM) using circles, squares, arrows,
#' and plates for the standard EMC hierarchical model.
#'
#' @param emc an emc object of type `standard`, `blocked`, or `diagonal`
#'
#' @returns Invisibly returns a DiagrammeR `grViz` graph object
#' @export
graphical_model <- function(emc) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package 'DiagrammeR' is required for this function. Please install it.")
  }
  if (is.null(emc) || !is.list(emc) || length(emc) == 0 || is.null(emc[[1]])) {
    stop("`emc` must be a non-empty EMC object.")
  }
  model_type <- emc[[1]]$type
  if (is.null(model_type) || !(model_type %in% c("standard", "blocked", "diagonal"))) {
    stop("`graphical_model()` currently supports only standard/blocked/diagonal EMC objects.")
  }

  par_names <- emc[[1]]$par_names
  if (is.null(par_names) || length(par_names) == 0) {
    alpha_samples <- emc[[1]]$samples$alpha
    if (!is.null(alpha_samples) && !is.null(dimnames(alpha_samples)[[1]])) {
      par_names <- dimnames(alpha_samples)[[1]]
    }
  }
  if (is.null(par_names) || length(par_names) == 0) {
    stop("Could not determine parameter names from `emc`.")
  }

  n_subjects <- emc[[1]]$n_subjects
  if (is.null(n_subjects)) {
    n_subjects <- length(emc[[1]]$subjects)
  }
  if (is.null(n_subjects) || n_subjects < 1) {
    n_subjects <- 1
  }

  raw_group_design <- emc[[1]]$group_designs
  has_group_design <- !is.null(raw_group_design) && length(raw_group_design) > 0
  group_design <- if (has_group_design) raw_group_design else NULL

  group_info <- NULL
  if (has_group_design) {
    group_info <- .extract_group_design_info(group_design, par_names)
  }

  has_random <- FALSE
  if (has_group_design) {
    has_random <- any(vapply(group_info, length, integer(1)) > 0)
    has_random <- has_random || tryCatch(get_n_random(par_names, group_design) > 0, error = function(e) FALSE)
  }

  prior <- emc[[1]]$prior
  has_prior_mu <- !is.null(prior$theta_mu_mean) || !is.null(prior$theta_mu_var)
  has_prior_sigma <- !is.null(prior$A) || !is.null(prior$v)
  has_prior_s <- has_random && (!is.null(prior$A_z) || !is.null(prior$v_Z))
  trial_plate_label <- "trials l = 1...L"

  all_random_blocks <- list()
  if (has_group_design) {
    all_random_blocks <- unlist(group_info, recursive = FALSE, use.names = FALSE)
  }

  # Safe compact layout knobs: rankdir, nodesep, ranksep, pad, fontsize, arrowsize.
  lines <- c(
    "digraph emc_pgm {",
    "  graph [rankdir=LR, splines=true, nodesep=0.05, ranksep=0.1, pad=0.02];",
    "  node [fontname=\"Helvetica\", fontsize=14, color=\"#2F2F2F\", penwidth=1.0];",
    "  edge [fontname=\"Helvetica\", fontsize=13, color=\"#4F4F4F\", arrowsize=0.7, penwidth=1.0];",
    "",
    "  Sigma [label=\"Sigma\", shape=circle, style=\"filled\", fillcolor=\"white\"];"
  )
  if (has_group_design) {
    lines <- c(lines, "  beta [label=\"beta\", shape=circle, style=\"filled\", fillcolor=\"white\"];")
  } else {
    lines <- c(lines, "  mu [label=\"mu\", shape=circle, style=\"filled\", fillcolor=\"white\"];")
  }

  if (has_random) {
    lines <- c(lines, "  s_g [label=<s<SUB>k,g</SUB>>, shape=circle, style=\"filled\", fillcolor=\"white\"];")
  }
  if (has_prior_mu) {
    lines <- c(lines, "  prior_mu [label=<mu<SUB>0</SUB>, Sigma<SUB>0</SUB>>, shape=box, style=\"filled\", fillcolor=\"white\"];")
  }
  if (has_prior_sigma) {
    lines <- c(lines, "  prior_sigma [label=\"A, v\", shape=box, style=\"filled\", fillcolor=\"white\"];")
  }
  if (has_prior_s) {
    lines <- c(lines, "  prior_s [label=<A<SUB>z</SUB>, V<SUB>z</SUB>>, shape=box, style=\"filled\", fillcolor=\"white\"];")
  }

  if (has_random) {
    lines <- c(
      lines,
      "",
      "  subgraph cluster_groups {",
      "    label=\"random groups h = 1...H_g\";",
      "    color=\"#9A9A9A\"; style=\"rounded,dashed\";",
      "    u_g [label=<u<SUB>k,g,h</SUB>>, shape=circle, style=\"filled\", fillcolor=\"white\"];",
      "  }"
    )
  }

  lines <- c(
    lines,
    "",
    "  subgraph cluster_subjects {",
    paste0("    label=\"subjects i = 1..", n_subjects, "\";"),
    "    color=\"#9A9A9A\"; style=\"rounded,dashed\";"
  )
  if (has_group_design) {
    lines <- c(
      lines,
      "    W_i [label=<W<SUB>i</SUB>>, shape=box, style=\"filled\", fillcolor=\"#E6E6E6\"];"
    )
  }

  if (has_random) {
    lines <- c(lines, "    Z_i [label=<Z<SUB>i</SUB>>, shape=box, style=\"filled\", fillcolor=\"#E6E6E6\"];")
  }

  lines <- c(
    lines,
    "",
    "    subgraph cluster_coefficients {",
    "      label=\"coefficients k = 1..K\";",
    "      color=\"#C8C8C8\"; style=\"rounded,dashed\";",
    "      alpha_i [label=<alpha<SUB>i</SUB>>, shape=circle, style=\"filled\", fillcolor=\"white\"];"
  )
  if (has_group_design) {
    lines <- c(lines, "      mu_i [label=<mu<SUB>i</SUB>>, shape=doublecircle, style=\"filled\", fillcolor=\"white\"];")
  }
  lines <- c(lines, "    }")

  lines <- c(
    lines,
    "",
    "    subgraph cluster_trials {",
    paste0("      label=\"", .dot_escape(trial_plate_label), "\";"),
    "      color=\"#BDBDBD\"; style=\"rounded,dashed\";",
    "      X_il [label=<X<SUB>i,l</SUB>>, shape=box, style=\"filled\", fillcolor=\"#E6E6E6\"];",
    "      theta_il [label=<theta<SUB>i,l</SUB>>, shape=doublecircle, style=\"filled\", fillcolor=\"white\"];",
    "      R_il [label=<R<SUB>i,l</SUB>>, shape=circle, style=\"filled\", fillcolor=\"#E6E6E6\"];",
    "      rt_il [label=<rt<SUB>i,l</SUB>>, shape=circle, style=\"filled\", fillcolor=\"#E6E6E6\"];",
    "    }",
    "  }",
    ""
  )

  if (has_prior_mu) {
    lines <- c(lines, .edge_line("prior_mu", if (has_group_design) "beta" else "mu"))
  }
  if (has_prior_sigma) {
    lines <- c(lines, .edge_line("prior_sigma", "Sigma"))
  }
  if (has_prior_s) {
    lines <- c(lines, .edge_line("prior_s", "s_g"))
  }
  if (has_random) {
    lines <- c(lines, .edge_line("s_g", "u_g"))
  }

  if (has_group_design) {
    lines <- c(lines, .edge_line("beta", "mu_i"))
    lines <- c(lines, .edge_line("W_i", "mu_i"))
    if (has_random) {
      lines <- c(lines, .edge_line("u_g", "mu_i"))
      lines <- c(lines, .edge_line("Z_i", "mu_i"))
    }
    lines <- c(lines, .edge_line("mu_i", "alpha_i"))
  } else {
    lines <- c(lines, .edge_line("mu", "alpha_i"))
  }
  lines <- c(lines, .edge_line("Sigma", "alpha_i"))
  lines <- c(lines, .edge_line("alpha_i", "theta_il"))
  lines <- c(lines, .edge_line("X_il", "theta_il"))
  lines <- c(lines, .edge_line("theta_il", "R_il"))
  lines <- c(lines, .edge_line("theta_il", "rt_il"))

  lines <- c(lines, "}")
  dot <- paste(lines, collapse = "\n")
  graph <- DiagrammeR::grViz(dot)
  print(graph)
  return(invisible(graph))
}


#' Make SEM Diagram
#'
#' @param emc an emc object
#' @param plot_values whether to plot the values or just the nodes/edges
#' @param cred_only whether to only plot credible values
#' @param par_names optional, if specified will overwrite the parameter names with
#' user-defined names
#' @param cut optional. A numeric value factor loadings smaller than this value will be omitted.
#' @param ... optional additional arguments passed to render_graph from DiagrammeR
#'
#' @returns Invisibly returns a DiagrammeR graph object
#' @export

make_SEM_diagram <- function(emc,
                          plot_values = TRUE,
                          cred_only = FALSE,
                          par_names = NULL,
                          cut = NULL,
                          ...){
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package 'DiagrammeR' is required for this function. Please install it.")
  }
  sem_settings <- emc[[1]]$sem_settings
  Lambda_mat <- sem_settings$Lambda_mat
  B_mat <- sem_settings$B_mat
  K_mat <- sem_settings$K_mat
  G_mat <- sem_settings$G_mat

  if(is.null(Lambda_mat)){
    if(!is.null(emc[[1]]$samples$lambda)){
      Lambda_mat <- emc[[1]]$samples$lambda[,,2]
      Lambda_mat[is_double_decimals(Lambda_mat)] <- Inf
    } else{
      stop("No loading matrix defined in your emc object")
    }
  }
  if(is.null(colnames(Lambda_mat))){
    colnames(Lambda_mat) <- paste0("F", 1:ncol(Lambda_mat))
  }

  n_pars <- nrow(Lambda_mat)
  n_factors <- ncol(Lambda_mat)

  if(is.null(B_mat)){
    B_mat <- matrix(0, nrow = n_factors, ncol = n_factors)
  }
  if(is.null(K_mat)){
    K_mat <- matrix(0, nrow = n_pars, ncol = 0)
  }
  if(is.null(G_mat)){
    G_mat <- matrix(0, nrow = n_factors, ncol = 0)
  }

  if(!is.null(par_names)){
    rownames(Lambda_mat) <- par_names
    rownames(K_mat) <- par_names
  }

  if(all(chain_n(emc) == 0)) plot_values <- FALSE

  n_cov <- max(ncol(K_mat), ncol(G_mat))
  covnames <- unique(c(colnames(K_mat), colnames(G_mat)))
  all_names <- c(rownames(Lambda_mat), colnames(Lambda_mat), covnames)
  all_shapes <- c(rep("circle", nrow(Lambda_mat[rowSums(abs(Lambda_mat) != 0),]) + ncol(Lambda_mat)),
                  rep("square", n_cov))
  n <- length(all_shapes)
  ndf <- DiagrammeR::create_node_df(n = n, shape = all_shapes, fontsize = 14, fixedsize = F,
                        label = all_names,color = "black", fillcolor = "white",
                        style = "filled")

  from <- c()
  to <- c()
  label <- c()
  if(plot_values){
    L_vals <- credint(emc, selection = "std_loadings", remove_constants = FALSE, digits = 2)
    if(any(B_mat != 0)) B_vals <- credint(emc, selection = "structural_regressors", remove_constants = FALSE, digits = 2)
    if(ncol(K_mat) > 1) K_vals <- credint(emc, selection = "regressors", remove_constants = FALSE, digits = 2)
    if(ncol(G_mat) > 1) G_vals <- credint(emc, selection = "factor_regressors", remove_constants = FALSE, digits = 2)
  }

  for(i in 1:ncol(Lambda_mat)){
    is_free <- Lambda_mat[,i] != 0
    if(cred_only){
      is_free <- is_free & apply(apply(L_vals[[i]], 1, sign), 2, FUN = function(x) return(length(unique(x)))) == 1
    }
    if(!is.null(cut)){
      is_free <- is_free & (abs(L_vals[[i]][,2]) > cut)
    }
    from <- c(from, rep(n_pars + i, sum(is_free)))
    to <- c(to, which(is_free))
    if(plot_values){
      label <- c(label, L_vals[[i]][,2][is_free])
    } else{
      label_tmp <- rep("*", sum(is_free))
      label_tmp[which(Lambda_mat[is_free,i] != Inf)] <- Lambda_mat[is_free & Lambda_mat[,i] != Inf,i]
      label <- c(label, label_tmp)
    }
  }

  if(any(B_mat != 0)){
    for(i in 1:ncol(B_mat)){
      is_free <- B_mat[,i] != 0
      if(cred_only) is_free <- is_free &
          apply(apply(B_vals[[i]], 1, sign), 2, FUN = function(x) return(length(unique(x)))) == 1
      from <- c(from, rep(n_pars + i, sum(is_free)))
      to <- c(to, n_pars + which(is_free))
      if(plot_values){
        label <- c(label, B_vals[[i]][,2][is_free])
      } else{
        label_tmp <- rep("*", sum(is_free))
        label_tmp[which(B_mat[is_free,i] != Inf)] <- B_mat[is_free & B_mat[,i] != Inf,i]
        label <- c(label, label_tmp)
      }
    }
  }

  if(any(K_mat != 0)){
    for(i in 1:ncol(K_mat)){
      is_free <- K_mat[,i] != 0
      if(cred_only) is_free <- is_free &
          apply(apply(K_vals[[i]], 1, sign), 2, FUN = function(x) return(length(unique(x)))) == 1
      from <- c(from, rep(n_pars + n_factors + i, sum(is_free)))
      to <- c(to, which(is_free))
      if(plot_values){
        label <- c(label, K_vals[[i]][,2][is_free])
      } else{
        label_tmp <- rep("*", sum(is_free))
        label_tmp[which(K_mat[is_free,i] != Inf)] <- K_mat[is_free & K_mat[,i] != Inf,i]
        label <- c(label, label_tmp)
      }
    }
  }

  if(any(G_mat) != 0){
    for(i in 1:ncol(G_mat)){
      is_free <- G_mat[,i] != 0
      if(cred_only) is_free <- is_free &
          apply(apply(G_vals[[i]], 1, sign), 2, FUN = function(x) return(length(unique(x)))) == 1
      from <- c(from, rep(n_pars + n_factors + i, sum(is_free)))
      to <- c(to, n_pars + which(is_free))
      if(plot_values){
        label <- c(label, G_vals[[i]][,2][is_free])
      } else{
        label_tmp <- rep("*", sum(is_free))
        label_tmp[which(G_mat[is_free,i] != Inf)] <- G_mat[is_free & G_mat[,i] != Inf,i]
        label <- c(label, label_tmp)
      }
    }
  }

  edf <- DiagrammeR::create_edge_df(from = from, to = to,
                        color = "#9B9999", fontsize = 14, label = label)

  all_incl <- unique(c(edf$from, edf$to))
  ndf$style[!ndf$id %in% all_incl] <- "invisible"
  graph <-DiagrammeR::create_graph(
      nodes_df = ndf,
      edges_df = edf,
      attr_theme = "tb")

  print(DiagrammeR::render_graph(graph, ...))
  return(invisible(graph))
}


rotater <- function(L_array, rot_fun, sign_convention = TRUE) {
  stopifnot(length(dim(L_array)) == 3)
  p <- dim(L_array)[1]; m <- dim(L_array)[2]; iters <- dim(L_array)[3]

  # Posterior mean loadings
  Abar <- apply(L_array, c(1, 2), mean)

  # Fit rotation once on the center matrix
  fit <- rot_fun(Abar)

  # Get rotated mean and the orthogonal transform T
  Abar_rot <- if (!is.null(fit$loadings)) unclass(fit$loadings) else
    if (!is.null(fit$L))        unclass(fit$L)        else
      stop("Rotation result missing $loadings/$L.")
  Tmat <- if (!is.null(fit$Th)) fit$Th else
    if (!is.null(fit$Tmat)) fit$Tmat else
      if (!is.null(fit$rotmat)) fit$rotmat else
        stop("Rotation result missing $Th/$Tmat/$rotmat (orthogonal transform).")

  # Optional: deterministic column sign (largest |loading| positive)
  if (sign_convention) {
    s <- rep(1, m)
    for (j in seq_len(m)) {
      idx <- which.max(abs(Abar_rot[, j]))
      s[j] <- if (Abar_rot[idx, j] < 0) -1 else 1
    }
    S <- diag(s, m, m)
    Tmat     <- Tmat %*% S
    Abar_rot <- Abar_rot %*% S
  }

  # Apply the SAME orthogonal T to every draw
  L_rot <- array(NA_real_, dim = c(p, m, iters), dimnames = dimnames(L_array))
  for (k in seq_len(iters)) {
    L_rot[, , k] <- L_array[, , k] %*% Tmat
  }
  return(L_rot)
}


#' Rotate loadings based on posterior median
#'
#' This function rotates factor loadings using a rotation function
#' based on the posterior median.
#'
#' @param emc An 'emc' object containing factor analysis results
#' @param rot_fun a rotation function for factor loadings, see also `GPArotation`W
#'
#' @return An 'emc' object with rotated factor loadings
#'
#' @export
rotate_loadings <- function(emc, rot_fun) {
  L_rot <- rotater(do.call(abind, lapply(emc, function(x) x$samples$lambda)), rot_fun)
  iters <- dim(L_rot)[3]
  start <- 0
  for(i in 1:length(emc)){
    end <-i*(iters/length(emc))
    emc[[i]]$samples$lambda <- L_rot[,,(start + 1):end]
    start <- end
  }
  class(emc) <- "emc"
  return(emc)
}
