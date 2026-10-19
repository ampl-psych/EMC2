# The validation kit (R/nn_validate.R) on the shipped neural likelihoods, in
# the cells the NLE project used for them, printed next to the project's own
# recorded numbers (Wuth, 2026; handover logs sbc/ref_sb.log,
# sbc/st0zero_sbc.log, sbc/rdm_small_screen.log, sbc/postshift_*.log,
# sbc/compare_lan.log). Their estimates came from other draws and data sets
# (and a WienR 1e-8 control), so agreement is to Monte Carlo error.
#
#   Rscript $(Rscript -e 'cat(system.file("scripts", "nle_validation_kit.R", package = "EMC2"))') \
#     [cores] [lan_card.rds] [out.rds]
#
# lan_card.rds: the converted HSSM ddm_uniform_st LAN (not shipped; in the
# source tree at tests/testthat/golden/lan_ddm_uniform_st.rds). Skipped if absent.

suppressPackageStartupMessages(library(EMC2))
args <- commandArgs(TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 1L
lan_path <- if (length(args) >= 2) args[2] else ""
out <- if (length(args) >= 3) args[3] else ""
cat(sprintf("EMC2 %s @ %s; cores %d; BLAS %s; VECLIB_MAXIMUM_THREADS=%s\n",
            as.character(utils::packageVersion("EMC2")), find.package("EMC2"), cores,
            utils::sessionInfo()$BLAS, Sys.getenv("VECLIB_MAXIMUM_THREADS")))
res <- list()
t_start <- proc.time()[["elapsed"]]
lap <- function(what) cat(sprintf("[%s: %.1f min]\n", what, (proc.time()[["elapsed"]] - t_start) / 60))
recorded <- function(x) { cat("  recorded: "); cat(x, sep = "\n            "); cat("\n") }

# --- cells (sbc_specs.R / SBCddm7.R), by name -----------------------------------------
ddm_mean <- c(v = 1, a = log(.8), t0 = log(.3), Z = qnorm(.5), SZ = qnorm(.2), sv = log(.1), st0 = log(.05))
ddm_sd <- c(v = .15, a = .15, t0 = .1, Z = .05, SZ = .15, sv = .15, st0 = .15)
cells <- list(
  cap256w_c4 = nn_cell(DDMnn, mean = ddm_mean, sd = ddm_sd),
  st0zero = nn_cell(function() DDMnn("ddm_st0zero"), mean = ddm_mean[-7], sd = ddm_sd[-7]),
  rdm_small = nn_cell(RDMnn, factors = list(S = c("left", "right")), Rlevels = c("left", "right"),
                      matchfun = function(d) d$S == d$lR,
                      formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
                      mean = c(v = log(1.5), v_lMTRUE = .4, B = log(1), A = log(.3), t0 = log(.2)), sd = .1))
for (nm in names(cells)) EMC2:::nle_artefact_info(cells[[nm]]$model)

# --- LAN, with its analytic twin (EMC2's DDM under the documented mapping) -------------
lan_twin <- function() list(
  type = "DDM",
  p_types = c(v = 0, a = log(1), z = qnorm(.5), t = log(.5), st = log(.05)),
  transform = list(func = c(v = "identity", a = "exp", z = "pnorm", t = "exp", st = "exp")),
  # the LAN's box for v, t, st (t - st >= 0 keeps the analytic t0 valid); a and
  # z loose, because make_data() checks bounds after Ttransform, where a and z
  # are the DDM's (2 a, Z * a)
  bound = list(minmax = cbind(v = c(-3, 3), a = c(0, Inf), z = c(0, Inf), t = c(.25, 2.25), st = c(.001, .25))),
  Ttransform = function(pars, dadm) {
    d <- cbind(v = pars[, "v"], a = 2 * pars[, "a"], sv = 0, t0 = pars[, "t"] - pars[, "st"],
               st0 = 2 * pars[, "st"], s = 1, Z = pars[, "z"], SZ = 0)
    d <- EMC2::DDM()$Ttransform(d, dadm)
    attr(d, "ok") <- attr(pars, "ok")
    d
  },
  rfun = function(data = NULL, pars) EMC2::DDM()$rfun(data, pars),
  dfun = function(rt, R, pars) EMC2::DDM()$dfun(rt, R, pars),
  pfun = function(rt, R, pars) EMC2::DDM()$pfun(rt, R, pars),
  log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
    EMC2:::log_likelihood_ddm(pars = pars, dadm = dadm, model = model, min_ll = min_ll))
lan_support <- function(p) p[, "t"] - p[, "st"]
if (nzchar(lan_path) && file.exists(lan_path)) {
  LAN <- register_nn_model(lan_path, twin = lan_twin)
  # score_bias_lan.R's draw region, in LAN units: v U(-2,2), a U(.5,2), z U(.35,.65),
  # t U(.3,1), st U(.02,.2); here a normal cell inside it
  cells$lan <- nn_cell(LAN, mean = c(v = 0, a = log(1), z = 0, t = log(.55), st = log(.07)),
                       sd = c(v = .45, a = .17, z = .09, t = .14, st = .3))
}
print(lapply(cells, function(c) c$box[!c$box$inside, ]))

# --- total mass -------------------------------------------------------------------------
cat("\n==== Total mass (50 prior draws) ====\n")
res$mass <- lapply(names(cells), function(nm) {
  tm <- nn_total_mass(cells[[nm]], n = 50, cores = cores,
                      support = if (nm == "lan") lan_support)
  print(tm); tm$summary
})
names(res$mass) <- names(cells)
lap("total mass")

# --- score bias ---------------------------------------------------------------------------
cat("\n==== Score bias (h = 0.01) ====\n")
sb <- function(nm, n, ...) { s <- nn_score_bias(cells[[nm]], n_draws = n, cores = cores, ...); print(s); s$summary }
res$score_bias <- list(
  cap256w_c4 = sb("cap256w_c4", 200),
  st0zero = sb("st0zero", 200),
  rdm_small = sb("rdm_small", 200))
recorded(c("cap256w_c4 (200 SBC draws): v +0.00036(.0001) a +0.01182(.0014) Z -0.00630(.0005)",
           "  sv -0.00009 SZ -0.00171(.0001) st0 +0.00392(.0003); control a +0.00007",
           "st0zero: v -0.00158(.0001) a +0.03356(.0012) Z +0.00142(.0004) sv +0.00002 SZ -0.00168(.0001)",
           "rdm_small: v -0.00350 v_lMTRUE -0.00184 B +0.00300 A +0.00080; control v -0.00012 B -0.00021"))
if (!is.null(cells$lan)) {
  # slow parameter vectors (a_DDM up to ~4, v ~ 0): the grid must reach 60 s
  # (as in score_bias_lan.R), or the control's g on a is -0.005, not 0
  res$score_bias$lan <- sb("lan", 80, pars = c("v", "a", "z"), support = lan_support,
                           rt_max = 60, n_grid = 240)
  recorded("LAN (compare_lan.log, 80 draws, another region): v +0.0038(.0024) log_a -0.1221(.0063) probit_w +0.0025(.0043)")
}
lap("score bias")

# --- posterior shift ----------------------------------------------------------------------
cat("\n==== Posterior shift (60 data sets x 400 trials, 2000 draws) ====\n")
ps <- function(nm) { s <- nn_posterior_shift(cells[[nm]], n_datasets = 60, n_trials = 400, n_draws = 2000,
                                             cores = cores); print(s); s[c("summary", "ess", "minutes")] }
res$posterior_shift <- list(cap256w_c4 = ps("cap256w_c4"))
recorded(c("cap256w_c4 (grid simulator) control/network/shift: v +.022/+.053/+.036 a -.137/-.186/-.044",
           "  t0 +.020/+.280/+.226 Z +.267/+.234/-.040 SZ +.010/+.001/-.010 sv -.067/-.059/+.006",
           "  st0 -.144/-.225/-.074; SBC N400 (500 reps) network t0 +.319, control t0 -.038"))
res$posterior_shift$st0zero <- ps("st0zero")
recorded(c("st0zero (exact simulator): v +.091/+.032/-.066 a +.127/+.304/+.160 t0 -.057/+.116/+.198",
           "  Z -.001/+.111/+.111 SZ +.018/+.059/+.042 sv -.061/-.053/+.008; SBC N400 t0 +.276"))
res$posterior_shift$rdm_small <- ps("rdm_small")
recorded("rdm_small SBC N400 (500 reps) network std bias: v +.034 v_lMTRUE -.026 B +.041 A +.090 t0 -.001")
if (!is.null(cells$lan)) res$posterior_shift$lan <- ps("lan")
lap("posterior shift")

if (nzchar(out)) saveRDS(res, out)
