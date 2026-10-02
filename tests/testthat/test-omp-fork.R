# Fork safety of the OpenMP likelihood backend (omp_state, R/sampling.R).
# The hang itself (GNU libgomp: a process forked from one with a live thread
# pool waits for ever in its first parallel region) cannot be a unit test;
# what is tested is the bookkeeping around it, with the platform, the pids and
# the pool release faked.

omp_dat <- local({
  des <- design(factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
                matchfun = function(d) d$S == d$lR, model = LBA,
                formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1), constants = c(sv = log(1)),
                report_p_vector = FALSE)
  p <- c(v = 1, v_lMTRUE = 1.5, B = log(1), A = log(.5), t0 = log(.2))
  set.seed(3)
  emc <- make_emc(make_data(p, des, n_trials = 20), des, type = "single", n_chains = 2, verbose = FALSE)
  props <- matrix(p, 20, length(p), byrow = TRUE, dimnames = list(NULL, names(p))) + rnorm(100, 0, .05)
  list(dadm = emc[[1]]$data[[1]], model = emc[[1]]$model, props = props)
})
omp_ll <- function() EMC2:::calc_ll_manager(omp_dat$props, omp_dat$dadm, omp_dat$model)

# set the package's fork-safety state for one test and put it back afterwards
local_omp_state <- function(..., env = parent.frame()) {
  st <- EMC2:::omp_state
  old <- as.list(st, all.names = TRUE)
  withr::defer({
    rm(list = ls(st, all.names = TRUE), envir = st)
    list2env(old, envir = st)
  }, envir = env)
  new <- list(...)
  for (nm in names(new)) assign(nm, new[[nm]], envir = st)
  st$warned <- FALSE
  invisible(st)
}
threaded <- function(env = parent.frame())
  withr::local_options(list(emc.ll_backend = "multithreaded", emc.n_threads = 2), .local_envir = env)

test_that("the package records the process it was loaded in", {
  st <- EMC2:::omp_state
  expect_identical(st$top_pid, Sys.getpid())
  expect_identical(st$hazard, EMC2:::omp_fork_hazard())
  if (Sys.info()[["sysname"]] %in% c("Darwin", "Windows")) expect_false(st$hazard)
  expect_null(st$pool_pid)
})

test_that("the fork hazard is told from the OpenMP runtime", {
  skip_on_os("windows")
  hz <- EMC2:::omp_fork_hazard
  expect_true(hz("/lib/x86_64-linux-gnu/libgomp.so.1", "Linux"))
  expect_true(hz("", "Linux"))                       # not recognised: assume the worst
  expect_false(hz("/apps/intel/lib/libiomp5.so", "Linux"))
  expect_false(hz("/usr/lib/llvm-14/lib/libomp.so.5", "Linux"))
  expect_false(hz("none", "Linux"))                  # built without OpenMP
  expect_false(hz("/lib/libgomp.so.1", "Darwin"))
  expect_type(EMC2:::omp_runtime(), "character")
  # the release is only ever done with libgomp
  if (!grepl("gomp", EMC2:::omp_runtime())) expect_identical(EMC2:::omp_release_pool(), -1L)
})

test_that("top-level process: the pool is released after a threaded likelihood", {
  serial <- omp_ll()
  st <- local_omp_state(hazard = TRUE, top_pid = Sys.getpid(), pool_pid = NULL)
  threaded()
  n <- 0
  local_mocked_bindings(omp_release_pool = function() { n <<- n + 1; 0L }, .package = "EMC2")
  expect_no_warning(l <- omp_ll())
  expect_equal(l, serial)
  expect_no_warning(omp_ll())
  expect_equal(n, 2)
  expect_null(st$pool_pid)       # nothing to warn its children about
})

test_that("forked process: keeps its pool, no release, and says so to its own children", {
  st <- local_omp_state(hazard = TRUE, top_pid = Sys.getpid() + 1L, pool_pid = NULL)
  threaded()
  local_mocked_bindings(omp_release_pool = function() stop("release called in a forked process"), .package = "EMC2")
  expect_no_warning(omp_ll())
  expect_identical(st$pool_pid, Sys.getpid())
  expect_no_warning(omp_ll())    # still threaded in this process
})

test_that("release not available: the session keeps threading, warns once and records its pid", {
  serial <- omp_ll()
  st <- local_omp_state(hazard = TRUE, top_pid = Sys.getpid(), pool_pid = NULL)
  threaded()
  local_mocked_bindings(omp_release_pool = function() -1L, .package = "EMC2")
  expect_warning(l <- omp_ll(), "could not be released")
  expect_equal(l, serial)
  expect_identical(st$pool_pid, Sys.getpid())
  expect_no_warning(omp_ll())
})

test_that("process forked from one with a live pool: serial likelihood, one warning", {
  serial <- omp_ll()
  st <- local_omp_state(hazard = TRUE, top_pid = Sys.getpid() + 1L, pool_pid = Sys.getpid() + 1L)
  threaded()
  local_mocked_bindings(
    calc_ll_multithreaded = function(...) stop("threaded likelihood called"),
    omp_release_pool = function() stop("release called"), .package = "EMC2")
  expect_warning(l <- omp_ll(), "serial likelihood")
  expect_identical(l, serial)
  expect_no_warning(l2 <- omp_ll())
  expect_identical(l2, serial)
  expect_identical(st$pool_pid, Sys.getpid() + 1L)
})

test_that("no hazard (macOS, Windows): nothing is released, recorded or replaced", {
  st <- local_omp_state(hazard = FALSE, top_pid = Sys.getpid(), pool_pid = Sys.getpid() + 1L)
  threaded()
  n <- 0
  local_mocked_bindings(
    calc_ll_multithreaded = function(...) { n <<- n + 1; matrix(0, 1, 1) },
    omp_release_pool = function() stop("release called"), .package = "EMC2")
  expect_no_warning(omp_ll())
  expect_equal(n, 1)
  expect_identical(st$pool_pid, Sys.getpid() + 1L)
})

test_that("omp_release_pool reports success or 'not available'", {
  expect_true(EMC2:::omp_release_pool() %in% c(0L, -1L))
})

test_that("threaded likelihood, then forked threaded children, on this platform", {
  skip_on_os("windows")
  skip_on_cran()
  serial <- omp_ll()
  threaded()
  expect_equal(omp_ll(), serial)
  job <- parallel::mcparallel(parallel::mclapply(1:2, function(i) omp_ll(), mc.cores = 2))
  res <- parallel::mccollect(job, wait = FALSE, timeout = 60)
  if (is.null(res)) {
    system(paste("pkill -9 -P", job$pid)); tools::pskill(job$pid, tools::SIGKILL)
    parallel::mccollect(job, wait = FALSE)
  }
  expect_false(is.null(res), label = "forked children returned (no hang)")
  for (l in res[[1]]) expect_equal(l, serial)
})
