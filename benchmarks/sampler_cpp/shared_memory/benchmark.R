# DDM5 sample-stage benchmark of the installed sampler.
args <- commandArgs(TRUE)
root <- Sys.getenv('EMC2_BENCHMARK_ROOT', '/home/nsteven/2026/StanCompare/emc2_shared')
.libPaths(c(Sys.getenv('EMC2_BENCHMARK_LIB', file.path(root, 'lib')), .libPaths()))
library(EMC2)
cores <- as.integer(args[1])
prefix <- args[2]
x <- readRDS(file.path(root, 'ddm5_start.rds'))
draws <- 100L
batches <- 10L
chains <- length(x)

last <- function(chain) EMC2:::remove_samples(chain, stage = 'sample', filter = chain$samples$idx - 1L)

setup <- system.time({
  context <- EMC2:::sampler_contexts(x, cores)
}, gcFirst = FALSE)[['elapsed']]
kept <- vector('list', batches)
invisible(gc())
elapsed <- system.time({
  for (batch in seq_len(batches)) {
    out <- EMC2:::run_native_chains(x, 'sample', draws, n_cores = cores,
      cores_for_chains = chains, context = context)
    kept[[batch]] <- lapply(out, `[[`, 'samples')
    x <- lapply(out, last)
  }
}, gcFirst = FALSE)[['elapsed']]
saveRDS(kept, paste0(prefix, '.rds'))
write.csv(data.frame(method = 'shared', cores_per_chain = cores, chains, draws = draws * batches,
                     setup_seconds = setup, sample_seconds = elapsed,
                     total_seconds = setup + elapsed), paste0(prefix, '.csv'), row.names = FALSE)
cat(cores, 'cores per chain; setup', setup, 'sample', elapsed, '\n')
