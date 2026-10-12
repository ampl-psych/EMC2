root <- Sys.getenv('EMC2_BENCHMARK_ROOT', '/home/nsteven/2026/StanCompare/emc2_shared')
results <- file.path(root, 'results')
files <- list.files(results, pattern = '^(fork|socket|socket_samples|shared)-(2|15)-[0-3]\\.rds$', full.names = TRUE)
reference <- readRDS(files[1])
for (file in files[-1]) stopifnot(identical(reference, readRDS(file)))
cat('All methods and core counts produced identical samples and RNG states.\n')
m <- read.csv(file.path(results, 'measurements.csv'))
m$memory_run <- as.logical(m$memory_run)
timing <- aggregate(cbind(setup_seconds, sample_seconds, total_seconds) ~ method + cores_per_chain,
                    m[!m$memory_run, ], median)
ram <- m[m$memory_run, c('method', 'cores_per_chain', grep('^peak_', names(m), value = TRUE))]
out <- merge(timing, ram, by = c('method', 'cores_per_chain'))
write.csv(out, file.path(results, 'summary.csv'), row.names = FALSE)
print(out, row.names = FALSE)
