# donothing <- function(i) repeat{}
#
# pdonothing <- function(i)
#   parallel::mclapply(1:3,donothing,mc.cores=3)
#
# test_parallel <- function() {
#   ncores=3
#   parallel::mclapply(1:ncores,pdonothing,mc.cores=ncores)
# }


rmvn <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                     method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE,
                     checkSymmetry = TRUE, rnorm = stats::rnorm)
# Alternative to rmvnorm not using %*%
{

  mymm <- function(y,z) t(apply(y,1,\(x) apply(x*z,2,sum)))


  if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                                    check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma))
    stop("mean and sigma have non-conforming size")
  R <- chol(sigma, pivot = TRUE)
  R <- R[, order(attr(R, "pivot"))]
  retval <- mvmult(matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994),R)
  # retval <- mymm(matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994),R)
  # retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}


dmvn <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE,
                   checkSymmetry = TRUE)
{
  if (is.vector(x)) x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean)))
      dim(mean) <- NULL
    if (length(mean) != p)
      stop("x and mean have non-conforming size")
  }
  if (!missing(sigma)) {
    if (p != ncol(sigma))
      stop("x and sigma have non-conforming size")
    if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                                      check.attributes = FALSE))
      stop("sigma must be a symmetric matrix")
  }
  dec <- tryCatch(base::chol(sigma), error = function(e) e)
  if (inherits(dec, "error")) {
    x.is.mu <- colSums_cpp(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf
  } else {
    tmp <- rcpp_forwardsolve(t(dec), t(x) - mean)
    rss <- colSums_cpp(tmp^2)
    logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if (log) logretval else exp(logretval)
}


# library(tmvtnorm)

mvmult <- function(a,b) {
  if (!is.matrix(a)) a <- matrix(a,nrow=1)
  if (!is.matrix(b)) b <- matrix(b,ncol=1)
  mat_mult(a,b)
}

checkSymmetricPositiveDefinite <- function (x, name = "sigma")
    {
    if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
        stop(sprintf("%s must be a symmetric matrix", name))
    }
    if (NROW(x) != NCOL(x)) {
        stop(sprintf("%s must be a square matrix", name))
    }
    if (any(diag(x) <= 0)) {
        stop(sprintf("%s all diagonal elements must be positive",
            name))
    }
    if (det(x) <= 0) {
        stop(sprintf("%s must be positive definite", name))
    }
}

checkTmvArgs <- function (mean, sigma, lower, upper)
{
    if (is.null(lower) || any(is.na(lower)))
        stop(sQuote("lower"), " not specified or contains NA")
    if (is.null(upper) || any(is.na(upper)))
        stop(sQuote("upper"), " not specified or contains NA")
    if (!is.numeric(mean) || !is.vector(mean))
        stop(sQuote("mean"), " is not a numeric vector")
    if (is.null(sigma) || any(is.na(sigma)))
        stop(sQuote("sigma"), " not specified or contains NA")
    if (!is.matrix(sigma)) {
        sigma <- as.matrix(sigma)
    }
    if (NCOL(lower) != NCOL(upper)) {
        stop("lower and upper have non-conforming size")
    }
    checkSymmetricPositiveDefinite(sigma)
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    if (length(lower) != length(mean) || length(upper) != length(mean)) {
        stop("mean, lower and upper must have the same length")
    }
    if (any(lower >= upper)) {
        stop("lower bound should be strictly less than the upper bound (lower<upper)")
    }
    cargs <- list(mean = mean, sigma = sigma, lower = lower,
        upper = upper)
    return(cargs)
}


rtmvnorm.rejection <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
    lower = rep(-Inf, length = length(mean)), upper = rep(Inf,
        length = length(mean)), D = diag(length(mean)))
{
    k <- length(mean)
    Y <- matrix(NA, n, k)
    numSamples <- n
    numAcceptedSamplesTotal <- 0
    r <- length(lower)
    d <- length(mean)
    if (r == d & identical(D, diag(d))) {
        alpha <- pmvnorm(lower = lower, upper = upper, mean = mean,
            sigma = sigma)
        if (alpha <= 0.01)
            warning(sprintf("Acceptance rate is very low (%s) and rejection sampling becomes inefficient. Consider using Gibbs sampling.",
                alpha))
        estimatedAlpha <- TRUE
    }
    else {
        alpha <- 1
        estimatedAlpha <- FALSE
    }
    while (numSamples > 0) {
        nproposals <- ifelse(numSamples/alpha > 1e+06, numSamples,
            ceiling(max(numSamples/alpha, 10)))
        X <- rmvn(nproposals, mean = mean, sigma = sigma)
        X2 <- mvmult(X,t(D))
        ind <- logical(nproposals)
        for (i in 1:nproposals) {
            ind[i] <- all(X2[i, ] >= lower & X2[i, ] <= upper)
        }
        numAcceptedSamples <- length(ind[ind == TRUE])
        if (length(numAcceptedSamples) == 0 || numAcceptedSamples ==
            0)
            next
        if (!estimatedAlpha) {
            alpha <- numAcceptedSamples/nproposals
            if (alpha <= 0.01)
                warning(sprintf("Acceptance rate is very low (%s) and rejection sampling becomes inefficient. Consider using Gibbs sampling.",
                  alpha))
        }
        numNeededSamples <- min(numAcceptedSamples, numSamples)
        Y[(numAcceptedSamplesTotal + 1):(numAcceptedSamplesTotal +
            numNeededSamples), ] <- X[which(ind)[1:numNeededSamples],
            ]
        numAcceptedSamplesTotal <- numAcceptedSamplesTotal +
            numAcceptedSamples
        numSamples <- numSamples - numAcceptedSamples
    }
    Y
}

rtmvnorm.linear.constraints <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
    H = NULL, lower = rep(-Inf, length = length(mean)), upper = rep(Inf,
        length = length(mean)), D = diag(length(mean)), algorithm,
    ...)
{
    d <- length(mean)
    if (!is.matrix(D) || det(D) == 0) {
        stop("D must be a (n x n) matrix with full rank n!")
    }
    alpha <- as.vector(lower - mvmult(D,mean))
    beta <- as.vector(upper - mvmult(D,mean))
    Dinv <- solve(D)
    if (!is.null(H)) {
        Tinv <- mvmult(mvmult(t(Dinv),H),Dinv)
        Z <- rtmvnorm(n, mean = rep(0, d), sigma = diag(d), H = Tinv,
            lower = alpha, upper = beta, algorithm = algorithm,
            ...)
    }
    else {
        T <- mvmult(mvmult(D,sigma),t(D))
        Z <- rtmvnorm(n, mean = rep(0, d), sigma = T, H = NULL,
            lower = alpha, upper = beta, algorithm = algorithm,
            ...)
    }
    X <- sweep(mvmult(Z,t(Dinv)), 2, FUN = "+", mean)
    return(X)
}


rtmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
    lower = rep(-Inf, length = length(mean)), upper = rep(Inf, length = length(mean)),
    D = diag(length(mean)), H = NULL, algorithm = c("rejection", "gibbs", "gibbsR"), ...)
{
    algorithm <- match.arg(algorithm)
    if (is.null(mean) && (is.null(sigma) || is.null(H))) {
        stop("Invalid arguments for ", sQuote("mean"), " and ",
            sQuote("sigma"), "/", sQuote("H"), ". Need at least mean vector and covariance or precision matrix.")
    }
    cargs <- checkTmvArgs(mean, sigma, lower, upper)
    mean <- cargs$mean
    sigma <- cargs$sigma
    lower <- cargs$lower
    upper <- cargs$upper
    if (!is.null(H) && !identical(sigma, diag(length(mean)))) {
        stop("Cannot give both covariance matrix sigma and precision matrix H arguments at the same time")
    }
    else if (!is.null(H) && !inherits(H, "sparseMatrix")) {
        checkSymmetricPositiveDefinite(H, name = "H")
    }
    if (n < 1 || !is.numeric(n) || n != as.integer(n) || length(n) >
        1) {
        stop("n must be a integer scalar > 0")
    }
    if (!is.matrix(D) || det(D) == 0) {
        stop("D must be a (n x n) matrix with full rank n!")
    }
    if (!identical(D, diag(length(mean)))) {
        retval <- rtmvnorm.linear.constraints(n = n, mean = mean,
            sigma = sigma, H = H, lower = lower, upper = upper,
            D = D, algorithm = algorithm, ...)
        return(retval)
    }
    else {
        if (algorithm == "rejection") {
            if (!is.null(H)) {
                retval <- rtmvnorm.rejection(n, mean, sigma = solve(H),
                  lower, upper, ...)
            }
            else {
                retval <- rtmvnorm.rejection(n, mean, sigma,
                  lower, upper, ...)
            }
        }
    }
    return(retval)
}


checkmvArgs <- function (lower, upper, mean, corr, sigma)
{
    if (!is.numeric(lower) || !is.vector(lower))
        stop(sQuote("lower"), " is not a numeric vector")
    if (!is.numeric(upper) || !is.vector(upper))
        stop(sQuote("upper"), " is not a numeric vector")
    if (!is.numeric(mean) || !is.vector(mean))
        stop(sQuote("mean"), " is not a numeric vector")
    if (is.null(lower) || any(is.na(lower)))
        stop(sQuote("lower"), " not specified or contains NA")
    if (is.null(upper) || any(is.na(upper)))
        stop(sQuote("upper"), " not specified or contains NA")
    rec <- cbind(lower, upper, mean)
    lower <- rec[, "lower"]
    upper <- rec[, "upper"]
    if (!all(lower <= upper))
        stop("at least one element of ", sQuote("lower"), " is larger than ",
            sQuote("upper"))
    mean <- rec[, "mean"]
    if (any(is.na(mean)))
        stop("mean contains NA")
    if (is.null(corr) && is.null(sigma)) {
        corr <- diag(length(lower))
    }
    else if (!is.null(corr) && !is.null(sigma)) {
        sigma <- NULL
        warning("both ", sQuote("corr"), " and ", sQuote("sigma"),
            " specified: ignoring ", sQuote("sigma"))
    }
    UNI <- FALSE
    if (!is.null(corr)) {
        if (!is.numeric(corr))
            stop(sQuote("corr"), " is not numeric")
        if (!is.matrix(corr)) {
            if (length(corr) == 1)
                UNI <- TRUE
            if (length(corr) != length(lower))
                stop(sQuote("diag(corr)"), " and ", sQuote("lower"),
                  " are of different length")
        }
        else {
            if (length(corr) == 1) {
                UNI <- TRUE
                corr <- corr[1, 1]
                if (length(lower) != 1)
                  stop(sQuote("corr"), " and ", sQuote("lower"),
                    " are of different length")
            }
            else {
                if (length(diag(corr)) != length(lower))
                  stop(sQuote("diag(corr)"), " and ", sQuote("lower"),
                    " are of different length")
                if (!chkcorr(corr))
                  stop(sQuote("corr"), " is not a correlation matrix")
            }
        }
    }
    if (!is.null(sigma)) {
        if (!is.numeric(sigma))
            stop(sQuote("sigma"), " is not numeric")
        if (!is.matrix(sigma)) {
            if (length(sigma) == 1)
                UNI <- TRUE
            if (length(sigma) != length(lower))
                stop(sQuote("diag(sigma)"), " and ", sQuote("lower"),
                  " are of different length")
        }
        else {
            if (length(sigma) == 1) {
                UNI <- TRUE
                sigma <- sigma[1, 1]
                if (length(lower) != 1)
                  stop(sQuote("sigma"), " and ", sQuote("lower"),
                    " are of different length")
            }
            else {
                if (length(diag(sigma)) != length(lower))
                  stop(sQuote("diag(sigma)"), " and ", sQuote("lower"),
                    " are of different length")
                if (!isTRUE(all.equal(sigma, t(sigma))) || any(diag(sigma) <
                  0))
                  stop(sQuote("sigma"), " is not a covariance matrix")
            }
        }
    }
    list(lower = lower, upper = upper, mean = mean, corr = corr,
        sigma = sigma, uni = UNI)
}

# mvt <- function (lower, upper, df, corr, delta, algorithm = GenzBretz(),
#     ...)
# {
#     if (...length() > 0)
#         algorithm <- GenzBretz(...)
#     else if (is.function(algorithm) || is.character(algorithm))
#         algorithm <- do.call(algorithm, list())
#     if (any(abs(lower - upper) < sqrt(.Machine$double.eps) *
#         (abs(lower) + abs(upper)) | lower == upper))
#         return(list(value = 0, error = 0, msg = "lower == upper"))
#     n <- ncol(corr)
#     if (is.null(n) || n < 2)
#         stop("dimension less then n = 2")
#     if (length(lower) != n)
#         stop("wrong dimensions")
#     if (length(upper) != n)
#         stop("wrong dimensions")
#     if (n > 1000)
#         stop("only dimensions 1 <= n <= 1000 allowed")
#     infin <- rep(2, n)
#     infin[isInf(upper)] <- 1
#     infin[isNInf(lower)] <- 0
#     infin[isNInf(lower) & isInf(upper)] <- -1
#     if (n >= 2 && ((isMiwa <- inherits(algorithm, "Miwa")) ||
#         inherits(algorithm, "TVPACK"))) {
#         if (any(infin == -1)) {
#             WhereBothInfIs <- which(infin == -1)
#             n <- n - length(WhereBothInfIs)
#             corr <- corr[-WhereBothInfIs, -WhereBothInfIs]
#             upper <- upper[-WhereBothInfIs]
#             lower <- lower[-WhereBothInfIs]
#             if (!missing(delta))
#                 delta <- delta[-WhereBothInfIs]
#             if (n <= 1) {
#                 if (n && !missing(delta)) {
#                   upper <- upper - delta
#                   lower <- lower - delta
#                 }
#                 return(list(value = if (n == 0) 1 else pnorm(upper) -
#                   pnorm(lower), error = 0, msg = "Normal Complettion (dim reduced to 1)"))
#             }
#             infin <- infin[-WhereBothInfIs]
#         }
#         if (isMiwa && any(infin == 0)) {
#             WhereNegativInfIs <- which(infin == 0)
#             inversecorr <- rep(1, n)
#             inversecorr[WhereNegativInfIs] <- -1
#             corr <- mvmult(mvmult(diag(inversecorr),corr),diag(inversecorr))
#             infin[WhereNegativInfIs] <- 1
#             tempsaveupper <- upper[WhereNegativInfIs]
#             upper[WhereNegativInfIs] <- -lower[WhereNegativInfIs]
#             lower[WhereNegativInfIs] <- -tempsaveupper
#         }
#     }
#     if (all(infin < 0))
#         return(list(value = 1, error = 0, msg = "Normal Completion"))
#     if (inherits(algorithm, "GenzBretz") && n > 1) {
#         corr <- matrix(as.vector(corr), ncol = n, byrow = TRUE)
#         corr <- corr[upper.tri(corr)]
#     }
#     ret <- probval(algorithm, n, df, lower, upper, infin, corr,
#         delta)
#     inform <- ret$inform
#     msg <- if (inform == 0)
#         "Normal Completion"
#     else if (inform == 1)
#         "Completion with error > abseps"
#     else if (inform == 2)
#         "Dimension greater 1000 or dimension < 1"
#     else if (inform == 3)
#         "Covariance matrix not positive semidefinite"
#     else inform
#     list(value = ret$value, error = ret$error, msg = msg, algo = class(algorithm))
# }



pmvnorm <- function (lower = -Inf, upper = Inf, mean = rep(0, length(lower)),
    corr = NULL, sigma = NULL, algorithm = GenzBretz(), keepAttr = TRUE,
    seed = NULL, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    carg <- checkmvArgs(lower = lower, upper = upper, mean = mean,
        corr = corr, sigma = sigma)
    if (!is.null(carg$corr)) {
        corr <- carg$corr
        if (carg$uni) {
            stop(sQuote("sigma"), " not specified: cannot compute pnorm")
        }
        else {
            lower <- carg$lower - carg$mean
            upper <- carg$upper - carg$mean
            mean <- rep(0, length(lower))
            RET <- mvtnorm:::mvt(lower = lower, upper = upper, df = 0,
                corr = corr, delta = mean, algorithm = algorithm,
                ...)
        }
    }
    else {
        if (carg$uni) {
            RET <- list(value = pnorm(carg$upper, mean = carg$mean,
                sd = sqrt(carg$sigma)) - pnorm(carg$lower, mean = carg$mean,
                sd = sqrt(carg$sigma)), error = 0, msg = "univariate: using pnorm")
        }
        else {
            lower <- (carg$lower - carg$mean)/sqrt(diag(carg$sigma))
            upper <- (carg$upper - carg$mean)/sqrt(diag(carg$sigma))
            mean <- rep(0, length(lower))
            corr <- cov2cor(carg$sigma)
            RET <- mvtnorm:::mvt(lower = lower, upper = upper, df = 0,
                corr = corr, delta = mean, algorithm = algorithm,
                ...)
        }
    }
    if (keepAttr)
        structure(RET$value, error = RET$error, msg = RET$msg,
            algorithm = RET$algorithm)
    else RET$value
}
