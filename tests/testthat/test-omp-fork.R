# Fork safety with GNU libgomp: a process that has opened an OpenMP parallel
# region and then forks leaves its children waiting on its thread pool for
# ever. The sampler releases the pool when a threaded entry point returns
# (omp_release_pool), so fit() in forked children (run_sbc(), mclapply()) after
# a fit in the session does not hang. The hang needs libgomp, so only a Linux
# gcc build (the ubuntu CI runner) exercises it; elsewhere the test just runs.

test_that("the OpenMP runtime is reported and the pool release is libgomp-only", {
  runtime <- EMC2:::omp_runtime()
  expect_type(runtime, "character")
  release <- EMC2:::omp_release_pool()
  expect_true(release %in% c(0L, -1L))
  if (!grepl("gomp", runtime)) expect_identical(release, -1L)
})

test_that("a threaded run, then threaded runs in forked children", {
  skip_on_os("windows")
  skip_on_cran()
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR,
                formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
  emc <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 2)
  sc <- function(iter) EMC2:::get_stop_criteria("preburn", list(iter = iter), type = "standard")
  emc <- run_emc(emc, "preburn", sc(5), cores_for_chains = 2, cores_per_chain = 2, verbose = FALSE)
  job <- parallel::mcparallel(parallel::mclapply(1:2, function(i)
    run_emc(emc, "preburn", sc(10), cores_for_chains = 2, cores_per_chain = 2,
            verbose = FALSE)[[1]]$samples$idx, mc.cores = 2))
  res <- parallel::mccollect(job, wait = FALSE, timeout = 300)
  if (is.null(res)) {
    system(paste("pkill -9 -P", job$pid)); tools::pskill(job$pid, tools::SIGKILL)
    parallel::mccollect(job, wait = FALSE)
  }
  expect_false(is.null(res), label = "forked children returned (no hang)")
  expect_equal(unlist(res[[1]]), c(10, 10))
})
