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
  retval <- mat_mult(matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994),R)
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
