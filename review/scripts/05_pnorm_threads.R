# PROVES: a DDM with the probit-scaled Z runs on subject threads and gives bit-identical draws across core settings (8 subjects, 20 preburn iterations); LNR likewise. Set OPENBLAS_NUM_THREADS=1 for clean timings.
# PRINTED 8b257e3d (machine also running R CMD check): DDM 1x1 3.5 s, 1x4 2.8 s, 2x2 2.1 s, all bit-identical to 1x1; LNR user/elapsed 0.99 (1x1), 1.75 (1x4), 1.93 (2x2).
suppressPackageStartupMessages(library(EMC2))
cat("EMC2 from:", find.package("EMC2"), "\n")
sub8 <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:8], ]; sub8$subjects <- droplevels(sub8$subjects)
sc <- function(stage, crit, type = "standard") EMC2:::get_stop_criteria(stage, crit, type = type)
cat("\n==== DDM with probit Z: are subject threads used? (8 subjects, 20 preburn iterations) ====\n")
des_ddm <- design(data = sub8, model = DDM, formula = list(v ~ 0 + S, a ~ E, t0 ~ 1, Z ~ 1, s ~ 1), constants = c(s = log(1)), report_p_vector = FALSE)
cat("DDM transforms:", paste(names(des_ddm$model()$transform$func), des_ddm$model()$transform$func, collapse = ", "), "\n")
set.seed(1); e <- make_emc(sub8, des_ddm, n_chains = 2, rt_resolution = 0.02)
for (cfg in list(c(1, 1), c(1, 4), c(2, 2))) {
  set.seed(1); t <- system.time(r <- run_emc(e, "preburn", sc("preburn", list(iter = 20)), cores_for_chains = cfg[1], cores_per_chain = cfg[2], verbose = FALSE))
  cat(sprintf("cores_for_chains=%d cores_per_chain=%d: elapsed %.1fs user %.1fs (user/elapsed = %.2f)\n", cfg[1], cfg[2], t[["elapsed"]], t[["user.self"]], t[["user.self"]] / t[["elapsed"]]))
  if (cfg[1] == 1 && cfg[2] == 1) ref <- r
  cat("  bit-identical to 1x1:", identical(lapply(r, function(x) x$samples$alpha), lapply(ref, function(x) x$samples$alpha)), "\n")
}
cat("\n==== LNR, 8 subjects: thread use ====\n")
des <- design(data = sub8, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
set.seed(1); e <- make_emc(sub8, des, n_chains = 2, rt_resolution = 0.02)
for (cfg in list(c(1, 1), c(1, 4), c(2, 2))) {
  set.seed(1); t <- system.time(r <- run_emc(e, "preburn", sc("preburn", list(iter = 50)), cores_for_chains = cfg[1], cores_per_chain = cfg[2], verbose = FALSE))
  cat(sprintf("cores_for_chains=%d cores_per_chain=%d: elapsed %.1fs user %.1fs (user/elapsed = %.2f)\n", cfg[1], cfg[2], t[["elapsed"]], t[["user.self"]], t[["user.self"]] / t[["elapsed"]]))
}
cat("done\n")
