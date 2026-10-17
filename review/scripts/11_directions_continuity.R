# PRINTED 8b257e3d: directions after 3 burn steps identical to step-1: FALSE (rebuilt every step; loc restarts); burn->adapt rebuilds TRUE; tail rebuilds TRUE; sample keeps TRUE; resume paths ran.
# PRINTED fix/joint-move-directions: identical TRUE (uses 25/75, loc kept adapting 2.3 1.3 1.7 1.6 -> 2.1 2.1 2.3 3.9; pre-tail adapt loc 5.2 5.8 3.5 3.5 after 400 uses); burn->adapt rebuilds TRUE; tail rebuilds TRUE; sample keeps TRUE; more adapt after the tail and adapt from a sample object ran.
# PROVES: the location move's directions are built once per stage (kept across the steps of burn), rebuilt at the burn->adapt
# boundary and at the sample-kernel tail, and kept in sample; the resume paths still run. Every run_emc(adapt) call ends with the tail,
# so the pre-tail state is read back from fileName (saved after each step, before the tail).
suppressPackageStartupMessages(library(EMC2)); cat("EMC2 from:", find.package("EMC2"), "\n")
old <- options(emc.adapt_converge = FALSE)
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]; dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
f <- tempfile(fileext = ".RData")
args <- list(cores_for_chains = 2, step_size = 25, verbose = FALSE, fileName = f)
sm <- function(e) attr(e[[1]]$samples, "scale_move"); fmt <- function(x) paste(round(x, 2), collapse = " ")
set.seed(7); e <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 2)
e <- do.call(run_emc, c(list(e, "preburn", sc("preburn", list(iter = 10))), args))
set.seed(1); b1 <- do.call(run_emc, c(list(e, "burn", sc("burn", list(iter = 25))), args)); s1 <- sm(b1)   # one step
set.seed(1); b3 <- do.call(run_emc, c(list(e, "burn", sc("burn", list(iter = 75))), args)); s3 <- sm(b3)   # three steps, same start
cat("burn: directions after 3 steps identical to the step-1 directions:", identical(s1$directions, s3$directions),
    "| uses", s1$uses, "/", s3$uses, "| loc after step 1:", fmt(s1$loc), "| after step 3:", fmt(s3$loc), "\n")
set.seed(2); a1 <- do.call(run_emc, c(list(b3, "adapt", sc("adapt", list(iter = 50))), args)); t1 <- sm(a1)   # two adapt steps, then the tail
load(f); pre <- sm(emc)                                                                                      # state saved after the last adapt step, before the tail
cat("burn -> adapt rebuilds directions:", !identical(s3$directions, pre$directions), "| adapt loc before the tail:", fmt(pre$loc), "uses", pre$uses, "\n")
cat("tail rebuilds directions:", !identical(pre$directions, t1$directions), "| frozen loc after the tail:", fmt(t1$loc), "| uses", t1$uses, "\n")
set.seed(3); s <- do.call(run_emc, c(list(a1, "sample", sc("sample", list(iter = 30))), args)); u <- sm(s)
cat("sample keeps step/loc/directions:", identical(u[c("step", "loc", "directions")], t1[c("step", "loc", "directions")]), "\n")
set.seed(4); more <- try(do.call(run_emc, c(list(a1, "adapt", sc("adapt", list(iter = 100))), args)))
cat("more adapt after the tail:", if (inherits(more, "try-error")) "ERROR" else "ran", "\n")
set.seed(5); back <- try(do.call(run_emc, c(list(s, "adapt", sc("adapt", list(iter = 100))), args)))
cat("adapt from a sample object:", if (inherits(back, "try-error")) "ERROR" else "ran", "\n")
options(old)
