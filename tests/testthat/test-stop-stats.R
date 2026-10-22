# The sample-stage stop statistics hardfit() reads next to the Rhat rule, and
# the one-check decision (sample_check()) on stored draws.

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
  expect_match(EMC2:::stop_stats_message(st), "min ESS=.*LL gap=.*subject spread=")
  # a non-finite subject log-likelihood: the likelihood statistics are NA, the ESS block survives
  inf <- samples_LNR; inf[[1]]$samples$subj_ll[1, inf[[1]]$samples$idx] <- -Inf
  si <- EMC2:::sample_stop_stats(inf)
  expect_true(is.na(si$ll_gap)); expect_true(is.na(si$subj_spread)); expect_true(is.finite(si$ess_min))
})

test_that("the vetoes hold a check the Rhat rule alone would pass, and NULL turns a criterion off", {
  n <- unname(chain_n(samples_LNR)[1, "sample"])
  far <- samples_LNR; far[[2]]$samples$subj_ll <- far[[2]]$samples$subj_ll - 300
  chk <- function(x, target = n, ...) EMC2:::sample_check(x, target = target, max_gd = 10, gd_quantile = .99, selection = "mu", ...)
  expect_true(chk(samples_LNR, min_ess = NULL, max_ll_gap = NULL, max_subj_spread = NULL)$met)
  expect_true(chk(far, min_ess = NULL, max_ll_gap = NULL, max_subj_spread = NULL)$met)
  r <- chk(far, min_ess = NULL, max_ll_gap = 100, max_subj_spread = NULL)
  expect_false(r$met); expect_gt(r$st$ll_gap, 100)
  expect_false(chk(samples_LNR, min_ess = 1e6, max_ll_gap = NULL, max_subj_spread = NULL)$met)
  one <- samples_LNR; one[[3]]$samples$subj_ll[1, ] <- one[[3]]$samples$subj_ll[1, ] - 40
  expect_false(chk(one, min_ess = NULL, max_ll_gap = NULL, max_subj_spread = 15)$met)
  expect_false(chk(samples_LNR, min_ess = NULL, max_ll_gap = NULL, max_subj_spread = NULL, target = n + 1)$met)
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
})

test_that("without a group level (type single) the subject parameters are the series", {
  single <- samples_LNR
  for(i in seq_along(single)){ single[[i]]$samples$theta_mu <- NULL; single[[i]]$samples$theta_var <- NULL }
  single[[1]]$type <- "single"
  st <- EMC2:::sample_stop_stats(single)
  expect_true(is.finite(st$ess_min)); expect_match(st$ess_min_par, ":")
  expect_true(is.finite(st$rhat_b)); expect_true(is.finite(st$ll_gap))
  expect_equal(EMC2:::rule_defaults("single")[c("max_gd", "gd_quantile")], list(max_gd = 1.05, gd_quantile = 1))
  expect_equal(EMC2:::rule_defaults("standard", max_gd = 1.2)[c("max_gd", "gd_quantile")], list(max_gd = 1.2, gd_quantile = .99))
})

test_that("the first third is dropped only when every criterion is then met", {
  n <- unname(chain_n(samples_LNR)[1, "sample"])
  it <- which(samples_LNR[[1]]$samples$stage == "sample")
  bad <- samples_LNR; k <- it[seq_len(round(n / 3))]
  bad[[1]]$samples$theta_mu[, k] <- bad[[1]]$samples$theta_mu[, k] + 20     # chain 1 far off in the first third only
  chk <- function(x, max_gd = 1.3) EMC2:::sample_check(x, target = 10, max_gd = max_gd, gd_quantile = .99, selection = "mu",
                                                       min_ess = NULL, max_ll_gap = NULL, max_subj_spread = NULL)
  # the Rhat rule fails on all draws and passes without the first third: dropped and met
  r <- chk(bad)
  expect_true(r$met); expect_equal(r$dropped, round(n / 3)); expect_equal(r$draws, n - round(n / 3))
  expect_equal(unname(chain_n(r$emc)[1, "sample"]), n - round(n / 3))
  # another criterion still fails without the first third: all draws kept, not met
  q <- EMC2:::sample_check(bad, target = 10, max_gd = 1.3, gd_quantile = .99, selection = "mu",
                           min_ess = 1e6, max_ll_gap = NULL, max_subj_spread = NULL)
  expect_false(q$met); expect_equal(q$dropped, 0L); expect_equal(unname(chain_n(q$emc)[1, "sample"]), n)
  # nothing to drop when the rule passes outright
  expect_equal(chk(samples_LNR, max_gd = 10)$dropped, 0L)
})
