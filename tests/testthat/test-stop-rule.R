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

test_that("the opt-in sample criteria are validated", {
  sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
  expect_equal(sc("sample", list(min_ess = 100, max_ll_gap = 100, max_subj_spread = 15))$max_subj_spread, 15)
  expect_error(sc("burn", list(min_ess = 100)), "sample stage")
  expect_error(sc("sample", list(max_ll_gap = -1)), "non-negative")
  expect_error(EMC2:::get_stop_criteria("sample", list(min_ess = 100), type = "single"), "single")
})

test_that("the sample stop statistics see a chain in another region and a subject apart", {
  st <- EMC2:::sample_stop_stats(samples_LNR)
  expect_true(all(c("ess_min", "ll_gap", "subj_spread", "rhat_b", "rhat_ident") %in% names(st)))
  expect_gt(st$ess_min, 0); expect_lt(st$ll_gap, 100)
  far <- samples_LNR; far[[2]]$samples$subj_ll <- far[[2]]$samples$subj_ll - 300   # every subject 300 down: sum 300 x n
  sf <- EMC2:::sample_stop_stats(far)
  expect_gt(sf$ll_gap, 300); expect_equal(sf$ll_worst_chain, 2)
  one <- samples_LNR; one[[3]]$samples$subj_ll[1, ] <- one[[3]]$samples$subj_ll[1, ] - 40
  so <- EMC2:::sample_stop_stats(one)
  expect_gt(so$subj_spread, 30); expect_equal(so$subj_spread_subject, rownames(one[[1]]$samples$subj_ll)[1])
  # the veto holds a run the Rhat rule alone would stop, and the reason is in the statistics
  crit <- list(max_gd = 10, selection = "mu", omit_mpsrf = TRUE)
  p0 <- EMC2:::check_progress(far, "sample", iter = 10, stop_criteria = crit, max_tries = 20, step_size = 10, n_cores = 1,
                              verbose = FALSE, progress = list(trys = 0, iters_total = 0))
  p1 <- EMC2:::check_progress(far, "sample", iter = 10, stop_criteria = c(crit, max_ll_gap = 100), max_tries = 20, step_size = 10,
                              n_cores = 1, verbose = FALSE, progress = list(trys = 0, iters_total = 0))
  expect_true(p0$done); expect_false(p1$done); expect_gt(p1$stop_stats$ll_gap, 100)
  p2 <- EMC2:::check_progress(samples_LNR, "sample", iter = 10, stop_criteria = c(crit, min_ess = 1e6), max_tries = 20, step_size = 10,
                              n_cores = 1, verbose = FALSE, progress = list(trys = 0, iters_total = 0))
  expect_false(p2$done)
  expect_null(p2$curr_min_es)                 # min_ess is not read as the older min_es by partial matching
  expect_match(EMC2:::stop_stats_message(st), "min ESS=.*LL gap=.*subject spread=")
})

test_that("projected Rhats: the disagreement axis sees chains apart, the identified axis does not", {
  set.seed(3)
  # two group means on a ridge: chains apart along (1, -1), the sum (1, 1) pinned
  MU <- lapply(1:4, function(i) { u <- rnorm(400, i * 2, 1); rbind(a = u + rnorm(400, 0, .05), b = -u + rnorm(400, 0, .05)) })
  r <- EMC2:::projected_rhats(MU)
  expect_gt(r$rhat_b, 1.5); expect_lt(r$rhat_ident, 1.05)
})

test_that("the sample stop statistics survive group variances that are not sampled (SEM: theta_var all NA)", {
  sem <- samples_LNR
  for(i in seq_along(sem)) sem[[i]]$samples$theta_var[] <- NA_real_
  st <- EMC2:::sample_stop_stats(sem)
  expect_true(is.finite(st$ess_min))                      # the group means still count
  expect_false(grepl("^lsd:", st$ess_min_par))
  expect_true(is.finite(st$ll_gap)); expect_true(is.finite(st$subj_spread))
  expect_true(is.na(EMC2:::ess_series(c(1, NA, 3, 4, 5)))); expect_true(is.na(EMC2:::ess_series(rep(2, 10))))
  none <- sem; for(i in seq_along(none)) none[[i]]$samples$theta_mu[] <- 0     # no finite non-constant series at all
  sn <- EMC2:::sample_stop_stats(none)
  expect_true(is.na(sn$ess_min)); expect_true(is.na(sn$ess_min_par))
  expect_false(grepl("min ESS", EMC2:::stop_stats_message(sn)))
  p <- EMC2:::check_progress(sem, "sample", iter = 10, stop_criteria = list(max_gd = 10, selection = "mu", omit_mpsrf = TRUE),
                             max_tries = 20, step_size = 10, n_cores = 1, verbose = FALSE, progress = list(trys = 0, iters_total = 0))
  expect_true(p$done)
})
