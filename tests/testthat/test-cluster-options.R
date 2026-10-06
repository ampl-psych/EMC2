test_that("socket workers (the Windows path of auto_mclapply) see the session's emc.* options", {
  skip_on_cran()
  withr::local_options(emc.sampler = "legacy")
  out <- EMC2:::cluster_lapply(1:2, function(i) getOption("emc.sampler"), cores = 1)
  expect_equal(unlist(out), c("legacy", "legacy"))
})
