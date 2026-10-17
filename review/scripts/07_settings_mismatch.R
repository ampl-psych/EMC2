# PRINTED on both builds: (a) mix length after the adapt tail: 4; more adapt after the tail: ran; (b) short counts vector: ran without complaint (the R side resets and resizes the settings before the C++ sampler sees them, so PR #202's C++ length check is not reached in practice).
# Proves: what happens when a saved object's particle settings do not match the kernel it is run with.
# (a) the realistic case, more adapt after adapt's tail tuned the 4-component sample kernel: handled by the R-side truncation on both builds;
# (b) a hand-edited mismatch that the R side does not repair: PR #202 stops with an EMC2 message naming the subject and lengths,
#     8b257e3d reaches the C++ adaptation arithmetic.
suppressPackageStartupMessages(library(EMC2)); cat("EMC2 from:", find.package("EMC2"), "\n")
old <- options(emc.adapt_converge = FALSE)
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]; dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
args <- list(cores_for_chains = 1, particle_factor = 30, step_size = 25, verbose = FALSE)
set.seed(1); emc <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 2)
emc <- do.call(run_emc, c(list(emc, "preburn", sc("preburn", list(iter = 10))), args))
emc <- do.call(run_emc, c(list(emc, "burn", sc("burn", list(iter = 25))), args))
emc <- do.call(run_emc, c(list(emc, "adapt", sc("adapt", list(min_unique = 20))), args))
cat("(a) mix length after the adapt tail:", length(attr(emc[[1]]$samples, "pm_settings")[[1]][[1]]$mix), "\n")
r <- try(do.call(run_emc, c(list(emc, "adapt", sc("adapt", list(iter = 25))), args)))
cat("(a) more adapt after the tail:", if (inherits(r, "try-error")) paste("ERROR:", conditionMessage(attr(r, "condition"))) else "ran", "\n")
# (b) hand-edited: proposal_counts and acc_counts one entry short of the kernel
bad <- emc
for (ch in seq_along(bad)) {
  pm <- attr(bad[[ch]]$samples, "pm_settings")
  pm[[1]][[1]]$proposal_counts <- pm[[1]][[1]]$proposal_counts[-1]; pm[[1]][[1]]$acc_counts <- pm[[1]][[1]]$acc_counts[-1]
  attr(bad[[ch]]$samples, "pm_settings") <- pm
}
r <- try(do.call(run_emc, c(list(bad, "adapt", sc("adapt", list(iter = 25))), args)))
cat("(b) short counts vector:", if (inherits(r, "try-error")) paste("ERROR:", conditionMessage(attr(r, "condition"))) else "ran without complaint", "\n")
options(old)
