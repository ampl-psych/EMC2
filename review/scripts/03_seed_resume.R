# PROVES: whether set.seed() before a run on an already initialised object changes the chain's own stream (5 more preburn iterations through run_emc, and native_run_chains alone with identical prepared arguments).
# PRINTED 8b257e3d: seed1==seed2 TRUE for both (set.seed() ignored after init).  PRINTED 3f24a807: seed1==seed2 FALSE for both (Random::mix reseeds from R).
suppressPackageStartupMessages(library(EMC2)); cat("EMC2 from:", find.package("EMC2"), "\n")
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]; dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR, formula = list(m ~ lM, s ~ 1, t0 ~ 1), report_p_vector = FALSE)
sc <- function(stage, crit) EMC2:::get_stop_criteria(stage, crit, type = "standard")
al <- function(e) lapply(e, function(x) x$samples$alpha)
set.seed(1); base <- run_emc(make_emc(dat, des, rt_resolution = 0.05, n_chains = 2), "preburn", sc("preburn", list(iter = 5)), cores_for_chains = 1, verbose = FALSE)
# preburn again, 5 more iterations (iter is the stage total)
set.seed(1); b1 <- run_emc(base, "preburn", sc("preburn", list(iter = 10)), cores_for_chains = 1, verbose = FALSE)
set.seed(2); b2 <- run_emc(base, "preburn", sc("preburn", list(iter = 10)), cores_for_chains = 1, verbose = FALSE)
set.seed(1); b3 <- run_emc(base, "preburn", sc("preburn", list(iter = 10)), cores_for_chains = 1, verbose = FALSE)
cat(sprintf("5 more preburn iterations from an initialised object: seed1==seed1 %s; seed1==seed2 %s (TRUE means set.seed() is ignored)\n",
    identical(al(b1), al(b3)), identical(al(b1), al(b2))))
# the native entry point alone, same prepared arguments
base <- EMC2:::restore_duplicates(base)
ctx <- EMC2:::sampler_contexts(base, 1)
args <- lapply(base, function(x) EMC2:::prepare_native_stage(x, "preburn", 5L, 100L, list(search_width = 1, kernel = NULL), NULL))
set.seed(1); n1 <- EMC2:::native_run_chains(args, ctx, "preburn", 5L, 1L, FALSE)
set.seed(2); n2 <- EMC2:::native_run_chains(args, ctx, "preburn", 5L, 1L, FALSE)
set.seed(1); n3 <- EMC2:::native_run_chains(args, ctx, "preburn", 5L, 1L, FALSE)
g <- function(x) lapply(x, function(ch) ch$samples$alpha)
cat("native_run_chains: seed1==seed1", identical(g(n1), g(n3)), "; seed1==seed2", identical(g(n1), g(n2)), "(TRUE means the chain stream ignores R's seed)\n")
