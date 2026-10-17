# PROVES: with a group-level design the location move is off (scale move on, loc untouched, no directions); for types diagonal, diagonal-gamma and blocked the directions are built and loc adapts.
# PRINTED 8b257e3d: group design: step adapted TRUE, loc unchanged (.2) TRUE, directions NULL; diagonal/diagonal-gamma/blocked: directions 4x4, loc 1.2-2.7 after 30 burn iterations.
suppressPackageStartupMessages(library(EMC2)); cat("EMC2 from:", find.package("EMC2"), "\n")
sc <- function(stage, crit, type = "standard") EMC2:::get_stop_criteria(stage, crit, type = type)
sub4 <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:4], ]; sub4$subjects <- droplevels(sub4$subjects)
sub4$x <- c(-1, -.3, .3, 1)[as.integer(sub4$subjects)]
des4 <- design(data = sub4, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
cat("\n==== group design: location move must be off, scale move on ====\n")
gd <- group_design(formula = list(m ~ x), data = sub4, subject_design = des4)
pri <- prior(des4, group_design = gd)
set.seed(2); g <- make_emc(sub4, des4, group_design = gd, prior_list = pri, n_chains = 2, rt_resolution = 0.05)
g <- run_emc(g, "preburn", sc("preburn", list(iter = 10)), cores_for_chains = 2, verbose = FALSE)
g <- run_emc(g, "burn", sc("burn", list(iter = 30)), cores_for_chains = 2, verbose = FALSE)
sm <- attr(g[[1]]$samples, "scale_move")
if (is.null(sm)) cat("no scale_move attribute after burn\n") else
cat("scale_move after burn with a group design: step adapted:", !all(sm$step == .5), " loc unchanged (.2):", all(sm$loc == .2), " directions:", if (is.null(sm$directions)) "NULL (location move off)" else "PRESENT", "\n")
for (type in c("diagonal", "diagonal-gamma", "blocked")) {
  cat("\n==== type", type, " ====\n")
  set.seed(2); d <- try(make_emc(sub4, des4, type = type, n_chains = 2, rt_resolution = 0.05, par_groups = if (type == "blocked") c(1, 1, 2, 2) else NULL))
  if (inherits(d, "try-error")) next
  d <- run_emc(d, "preburn", sc("preburn", list(iter = 10), type), cores_for_chains = 2, verbose = FALSE)
  d <- run_emc(d, "burn", sc("burn", list(iter = 30), type), cores_for_chains = 2, verbose = FALSE)
  sm <- attr(d[[1]]$samples, "scale_move")
  if (is.null(sm)) cat("no scale_move attribute\n") else cat("directions", if (is.null(sm$directions)) "NULL" else paste(dim(sm$directions), collapse = "x"), " loc:", paste(round(sm$loc, 2), collapse = " "), " step:", paste(round(sm$step, 2), collapse = " "), "\n")
}
cat("done\n")
