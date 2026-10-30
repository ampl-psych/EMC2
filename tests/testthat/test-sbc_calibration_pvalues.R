test_that(".sbc_ks_p flags non-uniform ranks, high p when uniform", {
  u <- (seq_len(500) - 0.5) / 500          # near-perfectly uniform
  expect_gt(EMC2:::.sbc_ks_p(u), 0.5)
  expect_lt(EMC2:::.sbc_ks_p(u^3), 0.001)  # piled toward 0
  expect_true(is.na(EMC2:::.sbc_ks_p(numeric(0))))
})

test_that(".sbc_chisq_p flags non-uniform bin counts, high p when flat", {
  expect_gt(EMC2:::.sbc_chisq_p(rep(50, 10)), 0.99)        # perfectly flat
  expect_lt(EMC2:::.sbc_chisq_p(c(200, rep(5, 9))), 0.001) # piled in one bin
  expect_true(is.na(EMC2:::.sbc_chisq_p(c(0, 0))))         # degenerate
})

test_that("plot_sbc_ecdf warns when requested recovery statistics are unavailable", {
  set.seed(1)
  ranks <- list(alpha = data.frame(v = runif(200), a = runif(200)))
  pdf(NULL)
  on.exit(dev.off())
  # bare ranks have no recovery stats -> warns about those, still draws p(KS)
  expect_warning(plot_sbc_ecdf(ranks, layout = NULL), "recovery statistics")
  # suppressing the three recovery stats leaves only p(KS): no warning
  expect_error(plot_sbc_ecdf(ranks, layout = NULL,
                             add_stats = list(coverage = FALSE, bias = FALSE,
                                              precision = FALSE)), NA)
})

test_that("plot_sbc_hist calculates p-values from the plotted bins", {
  set.seed(1)
  ranks <- list(alpha = data.frame(v = runif(200), a = runif(200)))
  pdf(NULL)
  on.exit(dev.off())
  bin_counts <- list()
  chisq_p <- EMC2:::.sbc_chisq_p
  local_mocked_bindings(.sbc_chisq_p = function(counts) {
    bin_counts[[length(bin_counts) + 1L]] <<- counts
    chisq_p(counts)
  }, .package = "EMC2")

  expected <- lapply(ranks$alpha, function(x) hist(x, breaks = 10, plot = FALSE)$counts)
  expect_silent(plot_sbc_hist(ranks, layout = NULL))
  expect_equal(bin_counts, unname(expected))
  expect_equal(vapply(bin_counts, chisq_p, numeric(1)),
               vapply(expected, function(x) chisq.test(x)$p.value, numeric(1)),
               ignore_attr = TRUE)

  expect_silent(plot_sbc_hist(ranks, layout = NULL, add_stats = FALSE))
})
