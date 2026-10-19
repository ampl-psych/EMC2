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

test_that("the progress line reports the rule's value, and the plain maximum when it differs", {
  expect_equal(EMC2:::rhat_message(1.0523, c(1.01, 1.0523)), "Rhat=1.052")
  expect_equal(EMC2:::rhat_message(1.08, c(1.01, 1.9)), "Rhat=1.080 (rule; max 1.900)")
  expect_null(EMC2:::rhat_message(NULL, numeric(0)))
  g <- EMC2:::check_gd(samples_LNR, "sample", max_gd = 1.1, mean_gd = NULL, omit_mpsrf = TRUE,
                       selection = c("alpha", "mu"), gd_quantile = .5)
  expect_equal(g$rule, max(stats::quantile(EMC2:::stage_gds(g$emc, c("alpha", "mu"), "sample")$alpha, .5),
                           EMC2:::stage_gds(g$emc, c("alpha", "mu"), "sample")$other))
})

test_that("check_progress counts the draws kept after a discard, and max_tries binds in sample", {
  crit <- list(max_gd = 1.0001, selection = "mu", omit_mpsrf = TRUE)     # never met: the first third may be dropped
  p <- suppressWarnings(EMC2:::check_progress(samples_LNR, "sample", iter = 1e6, stop_criteria = crit, max_tries = 1,
                                             step_size = 10, n_cores = 1, verbose = FALSE, progress = list(trys = 0, iters_total = 0)))
  expect_equal(p$total_iters_stage, unname(chain_n(p$emc)[1, "sample"]))
  expect_true(p$done)                                                    # iteration target unmet, but max_tries reached
  q <- EMC2:::check_progress(samples_LNR, "sample", iter = 1e6, stop_criteria = crit, max_tries = 5,
                             step_size = 10, n_cores = 1, verbose = FALSE, progress = list(trys = 0, iters_total = 0))
  expect_false(q$done)
})
