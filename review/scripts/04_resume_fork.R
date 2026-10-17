# PROVES (run with NOT_CRAN=true and PR202_TESTS=<dir holding PR #202's test-resume.R and test-omp-fork.R>, copies in ../pr202-tests): resume at every stage boundary incl. more adapt after the tail and adapt from a sample object; a threaded run then threaded fits in forked children (libgomp) return.
# PRINTED 8b257e3d: test-resume "more adapt after the sample-kernel tail" 6/6 pass; "particle settings are resized" 5/6 (fails only on on_exhausted); test-omp-fork second test 2/2 pass, first errors (omp_runtime not found); plain mclapply check: children returned after 0 s with idx 11 as expected.
# PRINTED 3f24a807: all pass.
suppressPackageStartupMessages({library(testthat); library(EMC2)})
cat("EMC2 from:", find.package("EMC2"), "\n")
cat("\n==== test-resume.R (from PR #202) ====\n")
print(test_file(file.path(Sys.getenv("PR202_TESTS", "../pr202-tests"), "test-resume.R"), reporter = "summary"))
cat("\n==== test-omp-fork.R (from PR #202) ====\n")
print(test_file(file.path(Sys.getenv("PR202_TESTS", "../pr202-tests"), "test-omp-fork.R"), reporter = "summary"))
cat("\n==== fork hang check without the helper functions: threaded run, then mclapply fits ====\n")
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]; dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
emc <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 2)
sc <- function(iter) EMC2:::get_stop_criteria("preburn", list(iter = iter), type = "standard")
emc <- run_emc(emc, "preburn", sc(5), cores_for_chains = 2, cores_per_chain = 2, verbose = FALSE)
job <- parallel::mcparallel(parallel::mclapply(1:2, function(i) run_emc(emc, "preburn", sc(10), cores_for_chains = 2, cores_per_chain = 2, verbose = FALSE)[[1]]$samples$idx, mc.cores = 2))
t0 <- Sys.time(); res <- parallel::mccollect(job, wait = FALSE, timeout = 240)
if (is.null(res)) { system(paste("pkill -9 -P", job$pid)); tools::pskill(job$pid, tools::SIGKILL); parallel::mccollect(job, wait = FALSE) }
cat("forked children returned:", !is.null(res), " after", round(as.numeric(Sys.time() - t0, units = "secs")), "s; result:", paste(unlist(res), collapse = ","), " expected:", emc[[1]]$samples$idx + 5, "\n")
