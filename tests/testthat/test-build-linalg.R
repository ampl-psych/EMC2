# On macOS configure links the Accelerate framework for BLAS / LAPACK. The link
# order decides which library a symbol binds to (two-level namespace: the first
# library on the line that exports it), so Accelerate must precede R's own
# -lRblas / -lRlapack, or every Armadillo product in the package (the sampler's
# linear algebra, the neural-likelihood flows) runs on R's reference BLAS: on an
# Apple M2 Max a flow's log density costs 68 us per trial instead of 4.5 us.
# Checked on the installed shared object with otool.

test_that("macOS: Accelerate precedes R's BLAS and LAPACK in EMC2's shared object", {
  skip_if_not(Sys.info()[["sysname"]] == "Darwin", "not macOS")
  skip_if_not(nzchar(Sys.which("otool")), "otool not found")
  so <- getLoadedDLLs()[["EMC2"]][["path"]]
  libs <- system2("otool", c("-L", shQuote(so)), stdout = TRUE)[-1]
  libs <- trimws(sub(" \\(.*$", "", libs))
  accelerate <- grep("Accelerate\\.framework", libs)
  rblas <- grep("libRblas|libRlapack", libs)
  skip_if(!length(accelerate), "Accelerate not linked (configure did not add it)")
  order <- sort(c(accelerate, rblas))
  expect_true(all(accelerate[1] < rblas),
              label = sprintf("link order: %s", paste(basename(libs[order]), collapse = " < ")))
})
