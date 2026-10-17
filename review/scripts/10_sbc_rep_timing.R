# PROVES: the cost of one hierarchical SBC replicate per candidate cell through the package's own worker (run_SBC_hierarchical_rep), fit() defaults, 20 subjects, 100 trials per cell, 3 chains x 1 core.
# PRINTED on a 4-vCPU Linux box: LBA (vignette design/prior) 2.1 min; DDM with probit Z (vignette design/prior, diagonal-gamma) 52.1 min.
# Time one hierarchical SBC replicate per candidate cell, through the package's own per-replicate worker.
suppressPackageStartupMessages(library(EMC2)); cat("EMC2 from:", find.package("EMC2"), "\n")
n_subjects <- 20; trials <- 100
matchfun <- function(d) d$S == d$lR
design_LBA <- design(factors = list(subjects = 1:n_subjects, S = c("left", "right")), Rlevels = c("left", "right"), matchfun = matchfun,
                     formula = list(v ~ lM, B ~ 1, t0 ~ 1, sv ~ 1, A ~ 1), constants = c(sv = log(1)), model = LBA)
prior_LBA <- prior(design_LBA, type = "diagonal-gamma", mu_mean = c(1.3, .7, log(.8), log(.2), log(.3)), mu_sd = c(.2, .1, .1, .05, .05), shape = 10, rate = c(.2, .2, .2, .1, .1))
design_DDM <- design(factors = list(subjects = 1:n_subjects, S = c("left", "right")), Rlevels = c("left", "right"),
                     formula = list(v ~ 1, a ~ 1, t0 ~ 1, s ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1, st0 ~ 1), constants = c(s = log(1)), model = DDM)
prior_DDM <- prior(design_DDM, type = "diagonal-gamma", pmean = c(1.2, log(.8), log(.3), qnorm(.5), log(.1), qnorm(.05), log(.05)),
                   psd = c(.15, .15, .1, .05, .1, .1, .15), shape = 10, rate = c(.2, .2, .2, .1, .1, .1, .1))
one <- function(label, design_in, prior_in, seed = 1) {
  dots <- EMC2:::add_defaults(list(cores_for_chains = 3), max_tries = 50, compress = FALSE, rt_resolution = 1e-12,
                              stop_criteria = list(min_es = 100, max_gd = 1.1, selection = c("alpha", "mu", "Sigma")), cores_per_chain = 1)
  set.seed(seed)
  prior_mu  <- plot(prior_in, design_in, do_plot = FALSE, N = 2, selection = "mu", return_mcmc = FALSE, remove_constants = FALSE, map = FALSE)[[1]]
  prior_var <- plot(prior_in, design_in, do_plot = FALSE, N = 2, selection = "Sigma", return_mcmc = FALSE, remove_constants = FALSE, map = FALSE)[[1]]
  td <- tempfile(paste0("sbc_", label)); dir.create(td)
  t0 <- Sys.time()
  res <- EMC2:::run_SBC_hierarchical_rep(1, design_in, prior_mu, prior_var, trials, n_subjects, prior_in, attr(prior_in, "type"), dots, td)
  mins <- as.numeric(Sys.time() - t0, units = "mins")
  cat(sprintf("%s: one replicate (%d subjects, %d trials/cell, 3 chains x 1 core) took %.1f min; rank_mu: %s\n", label, n_subjects, trials, mins, paste(res$rank_mu_row, collapse = " ")))
  f <- list.files(td, pattern = "rep_1.rds", full.names = TRUE); cat("  wrote:", basename(f), format(file.size(f), big.mark = ","), "bytes\n")
}
one("LBA_vignette", design_LBA, prior_LBA)
one("DDM_Z_vignette", design_DDM, prior_DDM)
