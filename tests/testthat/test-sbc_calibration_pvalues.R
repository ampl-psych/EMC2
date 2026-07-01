# Tests for the calibration p-value added to the SBC plots: the p(KS) /
# p(chisq) corner annotation and its add_stats wiring. The drawing helpers are
# internal, so they are reached with EMC2:::.

test_that(".sbc_fmt_p formats with a <0.001 floor", {
  expect_equal(EMC2:::.sbc_fmt_p(0.5), "0.500")
  expect_equal(EMC2:::.sbc_fmt_p(0.0001), "<0.001")
  expect_equal(EMC2:::.sbc_fmt_p(NA_real_), "NA")
})

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

test_that("add_stats resolves 'pvalue' into the free (bottomright) corner", {
  spec <- EMC2:::.sbc_resolve_stat_spec(TRUE, c("v", "a"))
  expect_equal(spec$pvalue$pos, "bottomright")
  expect_null(EMC2:::.sbc_resolve_stat_spec(list(pvalue = FALSE), c("v", "a"))$pvalue)
  expect_equal(EMC2:::.sbc_resolve_stat_spec(list(pvalue = "topleft"),
                                             c("v", "a"))$pvalue$pos, "topleft")
  expect_error(EMC2:::.sbc_resolve_stat_spec(list(foo = "top"), c("v", "a")),
               "subset of")
})

test_that("plot_sbc_ecdf keeps the corner stats (warns without recovery stats)", {
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

test_that("plot_sbc_hist shows only p(chisq) and never warns about recovery stats", {
  set.seed(1)
  ranks <- list(alpha = data.frame(v = runif(200), a = runif(200)))
  pdf(NULL)
  on.exit(dev.off())
  # no recovery-stat machinery anymore: silent, p(chisq) goes in the title
  expect_silent(plot_sbc_hist(ranks, layout = NULL))
  expect_silent(plot_sbc_hist(ranks, layout = NULL, add_stats = FALSE))
})
