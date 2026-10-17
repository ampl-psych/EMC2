# PROVES (run as R_LIBS=<lib> Rscript 02_behaviour.R): A set.seed() on a fresh object; B (invalid: 5 preburn iterations were already done, nothing ran; see 03); C cores 1x1 == 2x2; D DDM probit Z bit-identical across cores; E CDM serial warning; F on_singular$on_exhausted; G the joint moves fire (scale_move attribute after burn/adapt/sample, draws differ on vs off).
# PRINTED 8b257e3d: A TRUE/FALSE as expected; C TRUE; D 1x1==1x4==2x2 TRUE; E warnings: <none>; F error "unknown on_singular field(s): on_exhausted"; G after burn uses=60, directions 4x4, loc 2.5-3.0; sample stage keeps step/loc/directions TRUE; on-vs-off draws differ TRUE.
# PRINTED 3f24a807 (PR #202): E warns "Chains and subjects run one at a time ... the CDM likelihood draws from R's generator"; F accepted with a warning; G: no scale_move attribute after burn (moves start in adapt there), so the script stops at G.
# Behavioural checks of the C++ sampler. Run as: R_LIBS=<lib> Rscript behaviour.R
suppressPackageStartupMessages(library(EMC2))
cat("EMC2 from:", find.package("EMC2"), "\n"); print(packageDescription("EMC2")$RemoteSha); try(EMC2::emc2_build_info())
section <- function(x) cat("\n==== ", x, " ====\n", sep = "")
matchfun <- function(d) d$S == d$lR
sub2 <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]; sub2$subjects <- droplevels(sub2$subjects)
sub4 <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:4], ]; sub4$subjects <- droplevels(sub4$subjects)
des <- design(data = sub2, model = LNR, matchfun = matchfun, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
run <- function(emc, stage, crit, ...) run_emc(emc, stage, sc(stage, crit), verbose = FALSE, ...)
al <- function(e) lapply(e, function(x) x$samples$alpha)
mk <- function(d = sub2, ds = des, ...) make_emc(d, ds, rt_resolution = 0.05, n_chains = 2, ...)

section("A. set.seed() before a run on a fresh (uninitialised) object")
set.seed(1); r1 <- run(mk(), "preburn", list(iter = 10), cores_for_chains = 1)
set.seed(2); r2 <- run(mk(), "preburn", list(iter = 10), cores_for_chains = 1)
set.seed(1); r3 <- run(mk(), "preburn", list(iter = 10), cores_for_chains = 1)
cat("seed 1 == seed 1 (reproduces):", identical(al(r1), al(r3)), "\n")
cat("seed 1 == seed 2 (should be FALSE):", identical(al(r1), al(r2)), "\n")

section("B. set.seed() before a run on an initialised object (resume)")
set.seed(1); base <- run(mk(), "preburn", list(iter = 5), cores_for_chains = 1)
set.seed(1); b1 <- run(base, "preburn", list(iter = 5), cores_for_chains = 1)
set.seed(2); b2 <- run(base, "preburn", list(iter = 5), cores_for_chains = 1)
set.seed(1); b3 <- run(base, "preburn", list(iter = 5), cores_for_chains = 1)
cat("same seed reproduces:", identical(al(b1), al(b3)), "\n")
cat("different seed changes the draws (should be TRUE):", !identical(al(b1), al(b2)), "\n")

section("C. cores 1x1 vs 2x2 bit-identical (LNR)")
set.seed(1); c11 <- run(mk(), "preburn", list(iter = 10), cores_for_chains = 1, cores_per_chain = 1)
set.seed(1); c22 <- run(mk(), "preburn", list(iter = 10), cores_for_chains = 2, cores_per_chain = 2)
cat("1x1 == 2x2:", identical(al(c11), al(c22)), "\n")

section("D. probit-scaled Z (DDM) threaded and bit-identical across core settings")
des_ddm <- design(data = sub4, model = DDM, formula = list(v ~ 0 + S, a ~ E, t0 ~ 1, Z ~ 1, s ~ 1), constants = c(s = log(1)), report_p_vector = FALSE)
set.seed(1); e <- make_emc(sub4, des_ddm, n_chains = 2, rt_resolution = 0.02)
set.seed(1); t1 <- system.time(d11 <- run(e, "preburn", list(iter = 10), cores_for_chains = 1, cores_per_chain = 1))[["elapsed"]]
set.seed(1); t4 <- system.time(d14 <- run(e, "preburn", list(iter = 10), cores_for_chains = 1, cores_per_chain = 4))[["elapsed"]]
set.seed(1); t22 <- system.time(d22 <- run(e, "preburn", list(iter = 10), cores_for_chains = 2, cores_per_chain = 2))[["elapsed"]]
cat(sprintf("elapsed 1x1 %.1fs, 1x4 %.1fs, 2x2 %.1fs (threaded if 1x4 and 2x2 are clearly faster)\n", t1, t4, t22))
cat("1x1 == 1x4:", identical(al(d11), al(d14)), " 1x1 == 2x2:", identical(al(d11), al(d22)), "\n")
cat("pretransform of Z:", paste(names(des_ddm$model()$pre_transform$func), des_ddm$model()$pre_transform$func, collapse = ", "), "\n")

section("E. CDM: warns once that it runs serially?")
cdm_design <- design(model = CDM, factors = list(subjects = 1), Rlevels = "resp",
  formula = lapply(c("v", "theta", "a", "t0", "s", "sv"), function(p) as.formula(paste0(p, "~ 1"))), report_p_vector = FALSE)
cdm_pars <- c(v = log(1.1), theta = qnorm(0.55), a = log(1.0), t0 = log(0.2), s = log(1.0), sv = log(0.1))
set.seed(1); cdat <- make_data(cdm_pars, cdm_design, n_trials = 40)
cemc <- make_emc(cdat, cdm_design, type = "single", compress = FALSE, n_chains = 2, verbose = FALSE)
w <- NULL
res <- withCallingHandlers(try(run_emc(cemc, "preburn", EMC2:::get_stop_criteria("preburn", list(iter = 5), type = "single"), cores_for_chains = 2, cores_per_chain = 2, verbose = FALSE)),
  warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
cat("warnings:", if (length(w)) paste(w, collapse = " | ") else "<none>", "\n")

section("F. on_singular$on_exhausted from earlier scripts")
print(try(EMC2:::resolve_on_singular(list(max_retries = 2, on_exhausted = "carry_forward"))))

section("G. the joint moves fire: stages, directions, loc adaptation")
set.seed(3); g <- mk()
g <- run(g, "preburn", list(iter = 20), cores_for_chains = 2)
cat("after preburn scale_move attr:", !is.null(attr(g[[1]]$samples, "scale_move")), "\n")
g <- run(g, "burn", list(iter = 60), cores_for_chains = 2)
sm <- attr(g[[1]]$samples, "scale_move")
cat("after burn: uses =", sm$uses, " loc =", paste(round(sm$loc, 3), collapse = " "), " directions:", if (is.null(sm$directions)) "NULL" else paste(dim(sm$directions), collapse = "x"), "\n")
cat("lik_prec present after burn:", !is.null(g[[1]]$lik_prec), "\n")
old <- options(emc.adapt_converge = FALSE)
g <- run(g, "adapt", list(min_unique = 20), cores_for_chains = 2)
sm <- attr(g[[1]]$samples, "scale_move")
cat("after adapt+tail: uses =", sm$uses, " loc =", paste(round(sm$loc, 3), collapse = " "), " step =", paste(round(sm$step, 3), collapse = " "), "\n")
if (!is.null(sm$directions)) { cat("directions (rows = parameters, cols = moves):\n"); print(round(sm$directions, 3)) }
g2 <- run(g, "sample", list(iter = 30), cores_for_chains = 2)
sm2 <- attr(g2[[1]]$samples, "scale_move")
cat("sample stage keeps directions/loc/step fixed:", identical(sm2[c("step", "loc", "directions")], sm[c("step", "loc", "directions")]), "\n")
options(old)
# Does the location move change the draws at all? Compare a sample stage with the moves on and off from the same object and seed.
set.seed(5); on <- run(g, "sample", list(iter = 10), cores_for_chains = 1)
options(emc.scale_move = FALSE); set.seed(5); off <- run(g, "sample", list(iter = 10), cores_for_chains = 1); options(emc.scale_move = NULL)
cat("draws differ with the joint moves on vs off:", !identical(al(on), al(off)), "\n")

section("H. type with a group design: location move off, scale move on")
cat("done\n")
