# Timing of the neural-likelihood evaluators (src/flow_ddm.cpp, src/flow_race.cpp)
# on the shipped artefacts: shared theta (one cell), per-trial theta (every row
# distinct), and a 5-member ensemble (the cap256w_c4 member replicated), for
# n = 1, 1000, 10000 trials. Median milliseconds per call.
#
#   Rscript $(Rscript -e 'cat(system.file("scripts", "nle_bench_eval.R", package = "EMC2"))') [out.rds]
#
# The GEMMs go to whatever BLAS EMC2.so was linked against: Accelerate on macOS
# (configure adds -framework Accelerate), R's BLAS elsewhere. Report that with
# the numbers; on reference BLAS the batched products are ~10x slower.

suppressPackageStartupMessages(library(EMC2))
E <- asNamespace("EMC2")
out <- commandArgs(TRUE)[1]
cat(sprintf("EMC2 %s @ %s\nR %s, %s; R's BLAS: %s\n", as.character(utils::packageVersion("EMC2")),
            find.package("EMC2"), getRversion(), R.version$platform, utils::sessionInfo()$BLAS))
emc2_build_info()

tm <- function(f, budget = 1.5, max_reps = 2000) {
  f(); t0 <- proc.time()[[3]]; f(); t1 <- proc.time()[[3]] - t0
  reps <- max(3L, min(max_reps, as.integer(budget / max(t1, 1e-6))))
  x <- numeric(reps)
  for (i in seq_len(reps)) { s <- Sys.time(); f(); x[i] <- as.numeric(Sys.time() - s, units = "secs") }
  stats::median(x) * 1e3
}
draw <- function(fl, n) {
  l <- fl$bounds_sampled$lower; h <- fl$bounds_sampled$upper
  t(replicate(n, l + (h - l) * stats::runif(length(l), .05, .95)))
}

set.seed(1)
fl <- readRDS(file.path(E$nle_dir(), "ddm_cap256w_c4.rds"))$members[[1]]
p1 <- E$ddm_build(fl)
p5 <- E$ddm_build_ensemble(rep(list(fl), 5))
r <- readRDS(file.path(E$nle_dir(), "rdm_small.rds"))
pr <- E$flow_build(r)
rows <- list()
for (n in c(1, 1000, 10000)) {
  rt <- exp(stats::runif(n, log(.2), log(3))); R <- sample(1:2, n, TRUE)
  Tsh <- draw(fl, 1)[rep(1, n), , drop = FALSE]; Tpt <- draw(fl, n)
  Rsh <- draw(r, 1)[rep(1, n), , drop = FALSE]; Rpt <- draw(r, n)
  add <- function(model, case, ms)
    rows[[length(rows) + 1]] <<- data.frame(model = model, case = case, n = n, ms = ms)
  add("DDM cap256w_c4", "shared theta", tm(function() E$ddm_eval_trials_cpp(p1, Tsh, rt, R)))
  add("DDM cap256w_c4", "per-trial theta", tm(function() E$ddm_eval_trials_cpp(p1, Tpt, rt, R)))
  add("DDM cap256w_c4", "ens5 shared theta", tm(function() E$ddm_ens_eval_trials_cpp(p5, Tsh, rt, R)))
  add("DDM cap256w_c4", "ens5 per-trial theta", tm(function() E$ddm_ens_eval_trials_cpp(p5, Tpt, rt, R)))
  add("RDM rdm_small", "shared theta", tm(function() E$flow_eval_trials_cpp(pr, Rsh, rt)))
  add("RDM rdm_small", "per-trial theta", tm(function() E$flow_eval_trials_cpp(pr, Rpt, rt)))
}
res <- do.call(rbind, rows)
res <- res[order(res$model, res$case, res$n), ]
print(res, digits = 3, row.names = FALSE)
if (!is.na(out)) saveRDS(res, out)
