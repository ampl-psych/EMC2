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

test_that("the chain-count default is 6 for hierarchical types, 3 for single, and a given count is kept", {
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]; dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
  expect_length(suppressMessages(make_emc(dat, des, rt_resolution = 0.05)), 6)
  expect_length(suppressMessages(make_emc(dat, des, rt_resolution = 0.05, type = "single")), 3)
  expect_length(suppressMessages(make_emc(dat, des, rt_resolution = 0.05, n_chains = 4)), 4)
})

test_that("the sample ETA counts the run's own target and at least one more check until done", {
  e <- EMC2:::estimate_remaining_total_time("sample", tries_done = 12, elapsed_dt = as.difftime(60, units = "secs"), max_tries = 40,
                                            current_iters = 5000, target_iters = 5100, step_size = 100)
  expect_gt(as.numeric(e$min_time), 0); expect_gt(as.numeric(e$max_time), as.numeric(e$min_time))
  met <- EMC2:::estimate_remaining_total_time("sample", tries_done = 12, elapsed_dt = as.difftime(60, units = "secs"), max_tries = 40,
                                              current_iters = 6000, target_iters = 5100, step_size = 100)
  expect_equal(as.numeric(met$min_time), 60)                          # target met, not done: one more check
  fin <- EMC2:::estimate_remaining_total_time("sample", tries_done = 12, elapsed_dt = as.difftime(60, units = "secs"), max_tries = 40,
                                              current_iters = 6000, target_iters = 5100, step_size = 100, done = TRUE)
  expect_equal(as.numeric(fin$max_time), 0)
})

test_that("the first third is dropped only when the stop is then met, and the drop is recorded", {
  crit <- list(max_gd = 1.3, selection = "mu", omit_mpsrf = TRUE)     # all draws: about 1.7; without the first third: about 1.0
  n <- unname(chain_n(samples_LNR)[1, "sample"])
  it <- which(samples_LNR[[1]]$samples$stage == "sample")
  bad <- samples_LNR; k <- it[seq_len(round(n / 3))]
  bad[[1]]$samples$theta_mu[, k] <- bad[[1]]$samples$theta_mu[, k] + 20     # chain 1 far off in the first third only
  args <- list(stage = "sample", iter = 10, max_tries = 20, step_size = 10, n_cores = 1, verbose = FALSE,
               progress = list(trys = 0, iters_total = 0))
  # the Rhat rule fails on all draws and passes without the first third: dropped, stage stops, drop recorded
  p <- do.call(EMC2:::check_progress, c(list(bad, stop_criteria = crit), args))
  expect_true(p$done); expect_equal(p$total_iters_stage, n - round(n / 3))
  expect_length(p$discards, 1); expect_equal(p$discards[[1]]$removed, round(n / 3)); expect_equal(p$gd$discarded, round(n / 3))
  # the same, but another criterion would still fail without the first third: all draws kept, nothing recorded
  q <- do.call(EMC2:::check_progress, c(list(bad, stop_criteria = c(crit, min_ess = 1e6)), args))
  expect_false(q$done); expect_equal(q$total_iters_stage, n); expect_length(q$discards, 0); expect_equal(q$gd$discarded, 0)
  # check_gd() itself never subsets: it offers the shorter draws only when they pass the rule
  g <- EMC2:::check_gd(bad, "sample", max_gd = 1.3, mean_gd = NULL, omit_mpsrf = TRUE, selection = "mu")
  expect_false(g$gd_done); expect_equal(unname(chain_n(g$emc)[1, "sample"]), n); expect_true(g$short$gd_done)
  expect_null(EMC2:::check_gd(bad, "sample", max_gd = 1.0001, mean_gd = NULL, omit_mpsrf = TRUE, selection = "mu")$short)
})

test_that("a failing sample stage keeps all its draws, check after check", {
  skip_on_cran()
  old <- options(emc.adapt_converge = FALSE); on.exit(options(old))
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:3], ]; dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
  sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
  a <- list(cores_for_chains = 1, step_size = 50, verbose = FALSE, particle_factor = 30)
  set.seed(4); e <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 3)
  e <- do.call(run_emc, c(list(e, "preburn", sc("preburn", list(iter = 30))), a))
  e <- do.call(run_emc, c(list(e, "burn", sc("burn", list(iter = 100))), a))
  e <- do.call(run_emc, c(list(e, "adapt", sc("adapt", list(min_unique = 50))), a))
  # an unreachable Rhat target: 4 checks of 50 must leave all 200 draws (the old discard held the stage near 2/3)
  r <- suppressWarnings(do.call(run_emc, c(list(e, "sample", sc("sample", list(iter = 1000, max_gd = 1.0001))),
                                           c(a, max_tries = 4))))
  expect_equal(unname(chain_n(r)[1, "sample"]), 200)
  expect_length(r[[1]]$stop_record$discards, 0); expect_false(r[[1]]$stop_record$met)
  expect_null(r[[1]]$stop_record$aborted)
  # the same with the early stop at check 2: the stage ends there, recorded as aborted, not met
  expect_warning(s <- do.call(run_emc, c(list(e, "sample", sc("sample", list(iter = 1000, max_gd = 1.0001, abort_gd = 1.0001,
                                                                             abort_try = 2))), c(a, max_tries = 4))),
                 "stopped early at check 2")
  expect_equal(unname(chain_n(s)[1, "sample"]), 100)
  expect_false(s[[1]]$stop_record$met); expect_equal(s[[1]]$stop_record$aborted$try, 2)
  expect_gt(s[[1]]$stop_record$aborted$rule, 1.0001); expect_equal(s[[1]]$stop_record$criteria$abort_try, 2)
})

test_that("the early stop: validated, and taken only at its check with the rule's Rhat above the threshold", {
  sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
  expect_equal(sc("sample", list(max_gd = 1.1, abort_gd = 1.5))$abort_try, 20)
  expect_error(sc("burn", list(abort_gd = 1.5)), "sample stage")
  expect_error(sc("sample", list(abort_try = 20)), "needs `abort_gd`")
  expect_error(sc("sample", list(abort_gd = 1)), "above 1")
  expect_error(sc("sample", list(abort_gd = 1.5, abort_try = 2.5)), "positive integer")
  bad <- samples_LNR; it <- which(bad[[1]]$samples$stage == "sample")
  bad[[1]]$samples$theta_mu[, it] <- bad[[1]]$samples$theta_mu[, it] + 20      # chain 1 apart throughout: rule Rhat far above 1.5
  crit <- list(max_gd = 1.1, selection = "mu", omit_mpsrf = TRUE, abort_gd = 1.5, abort_try = 3)
  cp <- function(x, crit, trys) EMC2:::check_progress(x, "sample", iter = 10, stop_criteria = crit, max_tries = 20, step_size = 10,
                                                      n_cores = 1, verbose = FALSE, progress = list(trys = trys, iters_total = 0))
  expect_warning(p <- cp(bad, crit, 2), "check 3: the rule's Rhat")
  expect_true(p$done); expect_false(p$met); expect_equal(p$aborted$try, 3); expect_gt(p$aborted$rule, 1.5)
  q <- cp(bad, crit, 1)                                                          # check 2: not yet
  expect_false(q$done); expect_null(q$aborted)
  q <- cp(bad, crit, 3)                                                          # check 4: past it, only at the check itself
  expect_false(q$done); expect_null(q$aborted)
  ok <- cp(samples_LNR, c(crit[names(crit) != "max_gd"], max_gd = 1.0001), 2)    # unmet but below 1.5: carries on
  expect_false(ok$done); expect_null(ok$aborted)
  expect_null(cp(bad, crit[!names(crit) %in% c("abort_gd", "abort_try")], 2)$aborted)
})
