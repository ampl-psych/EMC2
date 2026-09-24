.onLoad <- function(libname, pkgname) {
  # On macOS EMC2's compiled code uses Accelerate's BLAS (configure links it).
  # Accelerate runs large matrix products (e.g. the neural-likelihood
  # evaluators' GEMM/GEMV on per-trial inputs) on several threads through
  # libdispatch, which crashes in forked processes: parallel chains
  # (mclapply), run_sbc() and the nn_* checks. EMC2's parallelism comes from
  # chains, replicates and OpenMP, so keep Accelerate single-threaded unless
  # the user has set VECLIB_MAXIMUM_THREADS. Accelerate reads it at its first
  # BLAS call, so this must happen before any EMC2 computation.
  if (identical(Sys.info()[["sysname"]], "Darwin") && !nzchar(Sys.getenv("VECLIB_MAXIMUM_THREADS")))
    Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
}
