# split_rhat() and the stop rule that reads it (check_gd()).

ar1 <- function(n, rho) as.numeric(stats::filter(rnorm(n, 0, sqrt(1 - rho^2)), rho, "recursive"))

test_that("split_rhat is the split potential scale reduction factor, with no correction and no transform", {
  set.seed(101)
  n <- 400; k <- 5; C <- 3
  X <- array(sapply(seq_len(k * C), function(i) ar1(n, .6)), c(n, k, C), dimnames = list(NULL, paste0("p", 1:k), NULL))
  X[, 2, ] <- exp(X[, 2, ])                                # all positive: coda would take logs
  X[150:260, 3, 2] <- X[150:260, 3, 2] - 3                 # one chain visits a shoulder
  r <- split_rhat(X)
  h <- n %/% 2
  direct <- sapply(seq_len(k), function(j){
    Y <- cbind(X[1:h, j, ], X[(h + 1):(2 * h), j, ])
    W <- mean(apply(Y, 2, var)); B <- h * var(colMeans(Y))
    sqrt(((h - 1) / h * W + B / h) / W)
  })
  expect_equal(unname(r), direct, tolerance = 1e-10)
  expect_identical(names(r), paste0("p", 1:k))
  # coda's point estimate carries a degrees-of-freedom correction, largest where
  # the half-chains' variances differ (the shoulder), and a log transform
  mcl <- coda::as.mcmc.list(lapply(seq_len(C), function(c) coda::as.mcmc(X[, , c])))
  cd <- coda::gelman.diag(split_mcl(mcl), autoburnin = FALSE, transform = TRUE, multivariate = FALSE)[[1]][, 1]
  expect_gt(cd[3], r[3])
  expect_equal(unname(gelman_diag_robust(mcl)), unname(r), tolerance = 1e-10)
  expect_identical(names(gelman_diag_robust(mcl, omit_mpsrf = FALSE))[k + 1], "mpsrf")
  # a variable that does not vary has no Rhat
  X[, 4, ] <- 1
  expect_true(is.nan(split_rhat(X)[4]))
  # large offsets do not cost precision
  expect_equal(unname(split_rhat(X[, 1, , drop = FALSE] + 1e6)), unname(r[1]), tolerance = 1e-6)
})

test_that("the stop rule's one-pass subject-level Rhats are gd_summary()'s", {
  g <- stage_gds(samples_LNR, c("alpha", "mu"), "sample")
  ref <- gd_summary(samples_LNR, selection = "alpha", stat = NULL, digits = 6)
  expect_equal(g$alpha, unname(unlist(ref)), tolerance = 1e-5)       # parameters, subjects within
  expect_equal(unname(g$other), unname(unlist(gd_summary(samples_LNR, selection = "mu", stat = NULL, digits = 6))), tolerance = 1e-5)
  expect_length(g$gd, length(g$alpha) + length(g$other))
  # without the first 16 draws of the stage: what subset() would keep, for both selections
  g2 <- stage_gds(samples_LNR, c("alpha", "mu"), "sample", filter = 16)
  short <- subset(samples_LNR, stage = "sample", filter = 16)
  expect_equal(unname(chain_n(short)[1, "sample"]), 34)
  expect_equal(g2$alpha, unname(unlist(gd_summary(short, selection = "alpha", stat = NULL, digits = 6))), tolerance = 1e-5)
  expect_equal(unname(g2$other), unname(unlist(gd_summary(short, selection = "mu", stat = NULL, digits = 6))), tolerance = 1e-5)
})

# samples_LNR with its sample-stage subject and group-mean draws replaced
fake_emc <- function(n, drift = 0, shoulder = FALSE, seed = 1){
  set.seed(seed)
  emc <- samples_LNR
  p <- dim(emc[[1]]$samples$alpha)[1]; ns <- dim(emc[[1]]$samples$alpha)[2]
  for(i in seq_along(emc)){
    S <- emc[[i]]$samples
    grow <- function(a) { d <- dim(a); k <- length(d); idx <- rep(seq_len(d[k]), length.out = n)
      if(k == 1) a[idx] else if(k == 2) a[, idx, drop = FALSE] else a[, , idx, drop = FALSE] }
    S <- rapply(S, function(a) if(is.null(dim(a))) { if(length(a) == S$idx) a[rep(seq_len(S$idx), length.out = n)] else a } else if(dim(a)[length(dim(a))] == S$idx) grow(a) else a, how = "replace")
    tr <- drift * exp(-seq_len(n) / (n / 12))                 # a transient every chain shares
    for(j in seq_len(p)) for(s in seq_len(ns)) S$alpha[j, s, ] <- .2 * ar1(n, .5) + tr
    for(j in seq_len(p)) S$theta_mu[j, ] <- .1 * ar1(n, .5) + tr
    if(shoulder && i == 2) S$alpha[1, 1, round(.55 * n):round(.75 * n)] <- S$alpha[1, 1, round(.55 * n):round(.75 * n)] - 1.2
    S$stage <- rep("sample", n); S$idx <- n
    emc[[i]]$samples <- S
  }
  emc
}

test_that("check_gd passes stationary chains and drops the start of chains that are still arriving", {
  emc <- fake_emc(600)
  out <- check_gd(emc, "sample", max_gd = 1.1, mean_gd = NULL, omit_mpsrf = TRUE, trys = 1, verbose = FALSE, selection = c("alpha", "mu"), iter = 600)
  expect_true(out$gd_done)
  expect_equal(unname(chain_n(out$emc)[1, "sample"]), 600)                    # nothing dropped
  emc <- fake_emc(600, drift = 3)
  out <- check_gd(emc, "sample", max_gd = 1.1, mean_gd = NULL, omit_mpsrf = TRUE, trys = 1, verbose = FALSE, selection = c("alpha", "mu"), iter = 600)
  expect_equal(unname(chain_n(out$emc)[1, "sample"]), 400)                    # the first third went
  expect_lt(max(out$gd), max(stage_gds(emc, c("alpha", "mu"), "sample")$gd))
  # the mean criterion of the burn stage reads the same statistic
  out <- check_gd(fake_emc(300), "sample", max_gd = NULL, mean_gd = 1.1, omit_mpsrf = TRUE, trys = 1, verbose = FALSE, selection = c("alpha", "mu"), iter = 300)
  expect_true(out$gd_done)
})

test_that("gd_quantile judges the subject level by a quantile and the rest by the largest", {
  emc <- fake_emc(1000, shoulder = TRUE, seed = 3)
  g <- stage_gds(emc, c("alpha", "mu"), "sample")
  expect_gt(max(g$alpha), 1.1)                                                # one subject's parameter, one chain
  expect_lt(gd_top(g, .9), 1.1)
  expect_identical(gd_top(g), max(g$gd))
  args <- list(emc, "sample", max_gd = 1.1, mean_gd = NULL, omit_mpsrf = TRUE, trys = 1, verbose = FALSE, selection = c("alpha", "mu"), iter = 1000)
  expect_false(do.call(check_gd, args)$gd_done)
  expect_true(do.call(check_gd, c(args, list(gd_quantile = .9)))$gd_done)
  # a group-level parameter above the bound fails whatever the quantile
  g$other[1] <- 1.3
  expect_gt(gd_top(g, .9), 1.1)
  expect_error(get_stop_criteria("sample", list(max_gd = 1.1, gd_quantile = 1.5), "standard"), "gd_quantile")
  expect_error(get_stop_criteria("burn", list(mean_gd = 1.1, gd_quantile = .99), "standard"), "max_gd")
  expect_identical(get_stop_criteria("sample", list(max_gd = 1.1, gd_quantile = .99), "standard")$gd_quantile, .99)
})
