# PROVES: posterior agreement with the joint moves on (default) and off (options(emc.scale_move = FALSE)) on a 4-subject LNR, 3 chains, fit stages, 4000 sample iterations; saves the two objects for 09.
# PRINTED 8b257e3d: mu |z| max 1.95 (0 above 2), sd ratio 0.93-0.98; sigma2 |z| max 1.76; alpha (subject 1) |z| max 2.59; frozen loc 2.55 2.70 2.65 4.52.
suppressPackageStartupMessages(library(EMC2)); cat("EMC2 from:", find.package("EMC2"), "\n")
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:4], ]; dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
one <- function(on, seed) {
  options(emc.scale_move = if (on) NULL else FALSE)
  set.seed(seed); e <- make_emc(dat, des, n_chains = 3, rt_resolution = 0.02)
  e <- run_emc(e, "preburn", sc("preburn", list(iter = 150)), cores_for_chains = 3, verbose = FALSE)
  e <- run_emc(e, "burn", sc("burn", list(mean_gd = 1.1)), cores_for_chains = 3, verbose = FALSE)
  e <- run_emc(e, "adapt", sc("adapt", list(min_unique = 600)), cores_for_chains = 3, verbose = FALSE)
  e <- run_emc(e, "sample", sc("sample", list(iter = 4000)), cores_for_chains = 3, verbose = FALSE)
  sm <- attr(e[[1]]$samples, "scale_move"); if (!is.null(sm)) cat("  frozen loc:", paste(round(sm$loc, 3), collapse = " "), " step:", paste(round(sm$step, 3), collapse = " "), "\n")
  e
}
stats <- function(e, selection) {
  m <- as.matrix(do.call(rbind, get_pars(e, selection = selection, stage = "sample", return_mcmc = TRUE, merge_chains = FALSE, flatten = TRUE)[[1]]))
  ess <- unlist(ess_summary(e, selection = selection, stage = "sample", stat = NULL, flatten = TRUE))
  data.frame(mean = colMeans(m), sd = apply(m, 2, sd), ess = as.numeric(ess)[seq_len(ncol(m))], row.names = colnames(m))
}
cmp <- function(a, b, what) {
  z <- (a$mean - b$mean) / sqrt(a$sd^2 / a$ess + b$sd^2 / b$ess); r <- a$sd / b$sd
  cat(sprintf("%-7s %3d parameters | mean difference |z|: max %.2f, n > 2: %d | sd ratio on/off: %.2f-%.2f | median ESS on/off: %.0f / %.0f\n",
              what, length(z), max(abs(z)), sum(abs(z) > 2), min(r), max(r), median(a$ess), median(b$ess)))
}
t0 <- Sys.time(); on <- one(TRUE, 11); off <- one(FALSE, 11); cat(sprintf("two fits: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
saveRDS(list(on = on, off = off), "posterior_onoff_emc.rds")
for (sel in c("mu", "sigma2", "alpha")) { a <- stats(on, sel); b <- stats(off, sel); cmp(a, b, sel); if (sel == "mu") print(round(cbind(on_mean = a$mean, off_mean = b$mean, on_sd = a$sd, off_sd = b$sd, on_ess = a$ess, off_ess = b$ess), 3)) }
