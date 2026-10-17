# PRINTED 8b257e3d: ERROR: chol(): decomposition failed.   PRINTED fix/joint-move-directions: ran; directions 4x4 finite TRUE, loc 4.06 2.01 1.65 1.65.
# PROVES: a non-positive-definite prior covariance of the group means no longer stops the stage with an Armadillo message.
suppressPackageStartupMessages(library(EMC2)); cat("EMC2 from:", find.package("EMC2"), "\n")
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]; dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
set.seed(3); e <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 2)
for (i in seq_along(e)) e[[i]]$prior$theta_mu_var[1, 1] <- 0     # a zero prior variance on the first group mean
e <- run_emc(e, "preburn", sc("preburn", list(iter = 10)), cores_for_chains = 1, verbose = FALSE)
r <- try(run_emc(e, "burn", sc("burn", list(iter = 20)), cores_for_chains = 1, verbose = FALSE))
if (inherits(r, "try-error")) cat("burn with a zero prior variance: ERROR:", conditionMessage(attr(r, "condition")), "\n") else {
  sm <- attr(r[[1]]$samples, "scale_move"); cat("burn with a zero prior variance: ran; directions", paste(dim(sm$directions), collapse = "x"), "finite:", all(is.finite(sm$directions)), "loc:", paste(round(sm$loc, 2), collapse = " "), "\n") }
