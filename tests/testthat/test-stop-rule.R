test_that("split Rhat matches its definition", {
  set.seed(101)
  n <- 400; h <- n %/% 2
  X <- array(rnorm(n * 3 * 3), c(n, 3, 3))
  X[, 2, 2] <- X[, 2, 2] + 1
  reference <- sapply(1:3, function(j) {
    Y <- cbind(X[1:h, j, ], X[(h + 1):n, j, ])
    W <- mean(apply(Y, 2, var))
    sqrt(((h - 1) / h * W + var(colMeans(Y))) / W)
  })
  expect_equal(unname(split_rhat(X)), reference, tolerance = 1e-10)
  expect_equal(unname(split_rhat(X + 1e6)), reference, tolerance = 1e-6)
})
