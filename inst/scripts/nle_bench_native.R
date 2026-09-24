# Likelihood cost of the neural-likelihood models on the R path (dfun/pfun per
# particle, c_name = NULL) vs the native branch of the compiled pipeline
# (calc_ll type "NN", src/model_NN.h): milliseconds per particle, the
# multithreaded backend's scaling, and whole fits (identical seeds) with the
# share of fit time spent in the likelihood.
#
#   Rscript $(Rscript -e 'cat(system.file("scripts", "nle_bench_native.R", package = "EMC2"))') [out.rds]
#
# DDMnn = cap256w_c4 with all 8 DDM parameters sampled; RDMnn = rdm_small, two
# accumulators. N = trials. The GEMMs go to whatever BLAS EMC2.so was linked
# against (Accelerate on macOS); report that with the numbers.

suppressPackageStartupMessages(library(EMC2))
E <- asNamespace("EMC2")
out <- commandArgs(TRUE)[1]
cat(sprintf("EMC2 %s @ %s\nR %s, %s; R's BLAS: %s\n", as.character(utils::packageVersion("EMC2")),
            find.package("EMC2"), getRversion(), R.version$platform, utils::sessionInfo()$BLAS))
emc2_build_info()
E$nle_artefact_info("ddm_cap256w_c4")
E$nle_artefact_info("rdm_small")

r_path <- function(model) { ml <- model(); ml$c_name <- NULL; function() ml }
tm <- function(f, reps = 5) { f(); stats::median(vapply(seq_len(reps), function(i) system.time(f())[["elapsed"]], 0)) }

setups <- list(
  DDMnn = list(model = DDMnn, Rlevels = c("a", "b"),
               formula = list(v ~ S, a ~ 1, t0 ~ 1, Z ~ 1, SZ ~ 1, sv ~ 1, st0 ~ 1), constants = c(s = 0),
               p = c(v = 1, v_Sb = -.5, a = log(1.2), t0 = log(.25), Z = 0, SZ = qnorm(.3),
                     sv = log(.5), st0 = log(.1))),
  RDMnn = list(model = RDMnn, Rlevels = c("a", "b"),
               formula = list(v ~ lR, B ~ 1, t0 ~ 1, A ~ 1), constants = c(s = 0),
               p = c(v = log(1), v_lRb = log(1.8), B = log(1), t0 = log(.25), A = log(.3))))

make_case <- function(st, N, n_chains = 3) {
  des <- suppressMessages(design(factors = list(subjects = 1, S = c("a", "b")), Rlevels = st$Rlevels,
                                 model = st$model, formula = st$formula, constants = st$constants,
                                 report_p_vector = FALSE))
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- st$p[names(p)]
  set.seed(1)
  dat <- make_data(p, des, n_trials = N / 2)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = n_chains, compress = TRUE,
                                   rt_resolution = .001, verbose = FALSE))
  list(p = p, emc = emc, dadm = emc[[1]]$data[[1]], model = emc[[1]]$model)
}

res <- list(per_particle = NULL, threads = NULL, fits = NULL)
n_prop <- 300
for (nm in names(setups)) for (N in c(400, 2000)) {
  cs <- make_case(setups[[nm]], N)
  set.seed(2)
  P <- matrix(cs$p, n_prop, length(cs$p), byrow = TRUE, dimnames = list(NULL, names(cs$p)))
  P <- P + matrix(stats::rnorm(length(P), sd = .1), n_prop)
  llN <- E$calc_ll_manager(P, cs$dadm, cs$model)
  llR <- E$calc_ll_manager(P, cs$dadm, r_path(cs$model))
  tR <- tm(function() E$calc_ll_manager(P, cs$dadm, r_path(cs$model)), 3) / n_prop * 1e3
  tN <- tm(function() E$calc_ll_manager(P, cs$dadm, cs$model)) / n_prop * 1e3
  res$per_particle <- rbind(res$per_particle, data.frame(
    model = nm, N = N, dadm_rows = nrow(cs$dadm), R_path_ms = tR, native_ms = tN, speedup = tR / tN,
    max_rel_diff = max(abs(llN - llR) / abs(llR))))
  old <- options(emc.ll_backend = "multithreaded")
  for (th in c(1, 2, 4, 8)) {
    options(emc.n_threads = th)
    tT <- tm(function() E$calc_ll_manager(P, cs$dadm, cs$model)) / n_prop * 1e3
    res$threads <- rbind(res$threads, data.frame(model = nm, N = N, threads = th, native_ms = tT))
  }
  options(old)
}
cat("\nPer particle (single thread; 300 particles around the generating values):\n")
print(res$per_particle, digits = 3, row.names = FALSE)
cat("\nMultithreaded backend (options(emc.ll_backend = \"multithreaded\", emc.n_threads = k)):\n")
thr <- res$threads
thr$speedup_vs_1 <- ave(thr$native_ms, thr$model, thr$N, FUN = function(x) x[1] / x)
print(thr, digits = 3, row.names = FALSE)

# Whole fits, 3 chains x 250 iterations, one core; identical seeds
sc <- list(preburn = list(iter = 50), burn = list(iter = 50), adapt = list(iter = 50), sample = list(iter = 100))
fit_with <- function(emc, model, prof = FALSE) {
  e <- emc
  for (k in seq_along(e)) e[[k]]$model <- model
  set.seed(99)
  if (prof) { pf <- tempfile(); utils::Rprof(pf, interval = .005) }
  t0 <- Sys.time()
  r <- suppressWarnings(suppressMessages(fit(e, stop_criteria = sc, cores_for_chains = 1, verbose = FALSE,
                                             step_size = 50)))
  t <- as.numeric(Sys.time() - t0, units = "secs")
  share <- NA
  if (prof) {
    utils::Rprof(NULL)
    sp <- utils::summaryRprof(pf)$by.total
    i <- grep("\"calc_ll_manager\"", rownames(sp), fixed = TRUE)[1]
    share <- if (is.na(i)) NA else sp$total.pct[i]
  }
  list(t = t, share = share, r = r)
}
for (nm in names(setups)) for (N in c(400, 2000)) {
  cs <- make_case(setups[[nm]], N)
  fR <- fit_with(cs$emc, r_path(cs$model), prof = TRUE)
  fN <- fit_with(cs$emc, cs$model, prof = TRUE)
  aR <- get_pars(fR$r, stage = "sample", merge_chains = TRUE, return_mcmc = FALSE)
  aN <- get_pars(fN$r, stage = "sample", merge_chains = TRUE, return_mcmc = FALSE)
  res$fits <- rbind(res$fits, data.frame(
    model = nm, N = N, R_path_s = fR$t, native_s = fN$t, speedup = fR$t / fN$t,
    ll_share_R = fR$share, ll_share_native = fN$share, max_draw_diff = max(abs(aR - aN))))
}
cat("\nFits (3 chains x 250 iterations, 1 core, same seed; ll_share = % of fit time in calc_ll_manager):\n")
print(res$fits, digits = 3, row.names = FALSE)
if (!is.na(out)) saveRDS(res, out)
