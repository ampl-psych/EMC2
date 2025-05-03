
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
#'
#' @return A list containing:
#'   \item{lambda_reordered}{Array of reordered loadings}
#'   \item{lambda_reordered_mcmc}{Array of reordered loadings as MCMC object}
#'   \item{lambda_hat}{Matrix of mean loadings after reordering}
#'   \item{v_vectors}{Matrix of permutation vectors}
#'   \item{c_vectors}{Matrix of sign-switching vectors}
#'
#' @examples
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
#' result <- align_loadings(lambda, verbose = TRUE, n_cores = 1)
#'
#' # Examine the aligned loadings
#' print(result$lambda_hat)
#'
#' @export
align_loadings <- function (emc = NULL, lambda = NULL, n_cores = 1, verbose = TRUE)
{
  maxIter <- 100; threshold <- 1e-06;
  rotate <- TRUE; printIter <- 1000;
  if(is.null(lambda) & is.null(emc)) stop("Need to supply either emc or lambda")
  if(is.null(lambda)){
    lambda <- get_pars(emc, selection = "loadings", merge_chains = TRUE, return_mcmc = FALSE)
  }
  energy <- apply(lambda^2, 2, mean)          # length q
  lambda <- lambda[,energy / max(energy) > 0.01,, drop = F]
  q <- ncol(lambda)
  p <- nrow(lambda)
  mcmcIterations <- dim(lambda)[3]
  threshold <- threshold * mcmcIterations * p * q
  lambda_varimax <- lambda
  if (rotate) {
    if (q > 1) {
      for (iter in 1:mcmcIterations) {
        lambda_varimax[,,iter] <- varimax(lambda[,,iter], normalize = F)$loadings
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
    n_iter <- emc[[1]]$samples$idx
    for(i in 1:length(emc)){
      emc[[i]]$samples$lambda <- lambda_reordered[,,((i-1)*(n_iter)+ 1):(i*n_iter), drop = FALSE]
    }
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

standardize_loadings <- function(loadings, residuals){
  stdize_set <- function(samples = NULL, idx = NULL, loadings = NULL, residuals = NULL){
    if(is.null(loadings)) loadings <- samples$samples$theta_lambda[,,idx, drop = F]
    if(is.null(residuals)) residuals <- samples$samples$theta_residuals[,idx]
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
  return(out)
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
#' @param ... Optional additional arguments
#'
#' @return No return value, creates a plot of group-level relations
#' @examples
#' # For a given set of hierarchical model samples we can make a
#' # correlation matrix plot.
#' plot_relations(samples_LNR, only_cred = TRUE, plot_cred = TRUE)
#' # We can also only plot the correlations where the credible interval does not include zero
#' plot_relations(samples_LNR, plot_means = TRUE, only_cred = TRUE)
#'
#' @export

plot_relations <- function(emc = NULL, stage = "sample",  plot_cred = TRUE,
                           plot_means = TRUE, only_cred = FALSE, nice_names = NULL,
                           selection = "correlation", ...){

  if(!selection %in% c("correlation", "loadings")){
    stop("Selection must be either 'correlation' or 'loadings'")
  }

  # for future factor model compatibility
  loadings <- NULL
  standardize <- TRUE
  do_corr <- TRUE
  corrs <- NULL

  # overwrite the optionals that were supplied
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if(!is.null(loadings)) do_corr <- F
  addCoef.col <- "black"
  if(!plot_means) addCoef.col <- NULL
  if(!is.null(emc)) sampled <- merge_chains(emc)
  if(do_corr || !is.null(corrs)){
    if(is.null(corrs)){
      values <- sampled$samples$theta_var[,,sampled$samples$stage == stage, drop = F]
    } else{
      values <- corrs
    }
    means <- cov2cor(apply(values, 1:2, mean))
  } else{
    # Now we assume loadings
    if(is.null(loadings)){
      if(standardize){
        loadings <- standardize_loadings(emc, stage = stage)
      } else{
        samples <- merge_chains(emc)
        loadings <- samples$samples$theta_lambda[,,samples$samples$stage == stage]
      }
    }
    values <- loadings
    means <- apply(values, 1:2, mean)
  }



  # You might have to play around with the x and y limits of the legend.

  # Only add this if you want confidence intervals
  if(plot_cred || only_cred){
    if(do_corr){
      for(i in 1:dim(values)[3]){
        values[,,i] <- cov2cor(values[,,i])
      }
    }

    cred <- aperm(apply(values, 1:2, quantile, probs = c(0.025, 0.975)))
    conf <- paste0("[", format(cred[,,1,drop = F], digits=1), ":", format(cred[,,2,drop = F], digits=1), "]")
    if(only_cred){
      is_cred <- unique(c(which(cred[,,1] > 0), which(cred[,,2] < 0)))
      conf[!1:length(conf) %in% is_cred] <- ""
      is_cred <- unique(c(which(t(cred[,,1] > 0)), which(t(cred[,,2] < 0))))
      means[!1:length(means) %in% is_cred] <- NA
    }
    cred <- round(cred, 2)


  }
  if(do_corr || !is.null(corrs)){
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
    if(plot_means){
      ys <- ys - 0.1
      text(xs, ys, conf, pos=1, cex=0.55, font = 2)
    } else{
      text(xs, ys, conf, pos=1, cex=0.7, font = 2)
    }
  }
}

# #'
# #' Factor diagram plot
# #' #Makes a factor diagram plot. Heavily based on the fa.diagram function of the `psych` package.
# #'
# #' @param emc An emc object
# #' @param stage Character. The stage from which to take the samples
# #' @param loadings An array of loadings. Can be alternatively supplied if emc is not supplied
# #' @param standardize Boolean. Whether to standardize the loadings
# #' @param simple Boolean. Whether the factor diagram should be simplified for visual clarity.
# #' @param only_cred Boolean. Whether to only plot the credible loadings
# #' @param cut Numeric. Mean loadings beneath this number will be excluded.
# #' @param nice_names Character vector. Alternative names to give the parameters
# #' @param factor_names Character vector. Names to give the different factors
# #' @param sort Boolean. Whether to sort the paramaters before plotting for visual clarity.
# #' @param adj Integer. Adjust to adjust loading values positions in the diagram if illegible.
# #' @param main Character vector. Title of the plot
# #' @param cex Integer. Font size
# #' @return NULL
# #' @examples \donttest{
# #' # For a given set of hierarchical factor model samples we can make a factor diagram
# #' make_factor_diagram(emc, only_cred = T)
# #' # We can also specify nice names and adjust the loading positions
# #' make_factor_diagram(emc, nice_names = paste0("V", 1:10), adj = 2)
# #' }

make_factor_diagram <- function(emc = NULL, stage = "sample",
                                loadings = NULL, standardize = TRUE,
                                simple = FALSE, only_cred = FALSE,
                                cut = 0, nice_names = NULL,
                                factor_names = NULL, sort = TRUE,
                                adj = 1, main = NULL, cex = NULL){
  if(is.null(loadings)){
    if(standardize){
      loadings <- standardize_loadings(emc, stage = stage)
    } else{
      samples <- merge_chains(emc)
      loadings <- samples$samples$theta_lambda[,,samples$samples$stage == stage]
    }
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


