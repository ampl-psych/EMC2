# Acceptance check for the four NLE networks registered through register_nn_model():
# log densities of the registered models (dfun, the same path the likelihood uses)
# against the NLE project's own reference evaluators, at 1e-10.
#
#   Rscript nle_acceptance.R <port dir> <reg_s402 card.rds> <reg_s402 regress.json> <DDMreg.R> [n]
#
# <port dir>   the R port (flow_ddm.R, flow_race.R): tests/testthat/port in the source tree
# reg_s402     the card made by inst/scripts/regress_to_card.R, its regress.json, and the
#              NLE project's evaluation/DDMreg.R (the reference evaluator of that net)
#
# The three shipped artefacts are compared with the R port (per row, per-trial theta;
# flows), reg_s402 with DDMreg.R's reg_log_pdf. A log-density difference d bounds the
# score-bias error at d / (2 h) (h = 0.01): 1e-10 -> 5e-9, inside the 1e-6 the
# migration must reproduce.
suppressPackageStartupMessages(library(EMC2))
args <- commandArgs(TRUE)
stopifnot(length(args) >= 4)
port_dir <- args[1]; card_path <- args[2]; json_path <- args[3]; reg_r <- args[4]
n <- if (length(args) >= 5) as.integer(args[5]) else 400L
nle <- function(f) getFromNamespace(f, "EMC2")
cat(sprintf("EMC2 %s @ %s; BLAS %s\n", as.character(utils::packageVersion("EMC2")), find.package("EMC2"),
            utils::sessionInfo()$BLAS))
port <- new.env()
sys.source(file.path(port_dir, "flow_race.R"), envir = port)
sys.source(file.path(port_dir, "flow_ddm.R"), envir = port)
port_fix <- function(m) {
  fix <- function(mlp) { mlp$layers <- lapply(mlp$layers, function(l)
    list(W = as.matrix(l$W), b = as.numeric(l$b))); mlp }
  if (!is.null(m$mlp)) m$mlp <- fix(m$mlp)
  if (!is.null(m$flow_mlp)) { m$flow_mlp <- fix(m$flow_mlp); m$classifier_mlp <- fix(m$classifier_mlp) }
  m
}
set.seed(6001)
res <- list()

# in-box rows on the natural scale, columns in the net's order, and the same on the sampled scale
draw <- function(reg, n) {
  m <- t(replicate(n, reg$lower + (reg$upper - reg$lower) * runif(length(reg$lower), .02, .98)))
  colnames(m) <- reg$pars
  th <- m
  for (j in seq_along(reg$pars)) th[, j] <- switch(reg$transforms[j], identity = m[, j], log = log(m[, j]),
                                                   probit = qnorm(m[, j]))
  list(pars = m, theta = th)
}
rt <- exp(runif(n, log(.06), log(3))); R <- sample(1:2, n, TRUE)

flows <- c(cap256w_c4 = "ddm_cap256w_c4", st0zero = "ddm_st0zero", rdm_small = "rdm_small")
for (nm in names(flows)) {
  m <- register_nn_model(flows[[nm]])()
  reg <- m$nn
  nle("nle_artefact_info")(flows[[nm]])
  x <- readRDS(file.path(nle("nle_dir")(), paste0(flows[[nm]], ".rds")))
  fl <- port_fix(if (!is.null(x$members)) x$members[[1]] else x)
  d <- draw(reg, n)
  if (reg$kind == "flow_joint") {
    ours <- log(m$dfun(rt, R, d$pars))
    ref <- vapply(seq_len(n), function(i) port$ddm_eval(fl, d$theta[i, ], rt[i], R[i])$log_pdf, 0)
  } else {
    ours <- log(m$dfun(rt, d$pars))
    ref <- vapply(seq_len(n), function(i) port$flow_eval(fl, d$theta[i, ], rt[i])$log_pdf, 0)
  }
  res[[nm]] <- max(abs(ours - ref))
  cat(sprintf("%-11s %s: max |log density - R port| over %d per-trial rows = %.2e\n", nm, reg$kind, n, res[[nm]]))
}

# reg_s402: registered by path; reference = the project's DDMreg.R evaluator
m <- register_nn_model(card_path)()
reg <- m$nn
ref_env <- new.env()
Sys.setenv(REG_JSON = json_path)
sys.source(reg_r, envir = ref_env)
d <- draw(reg, n)
stopifnot(identical(colnames(d$pars), c("v", "a", "t0", "s", "Z", "SZ", "sv", "st0")))
ours <- log(m$dfun(rt, R, d$pars))
ref <- ref_env$reg_log_pdf(rt, R, d$pars)
# dfun returns a density: rows whose log density is below -600 underflow to 0 there (a net
# has no floor; rt << t0 rows sit near exp(-4000)); compare the rest on the log scale
ok <- ref > -600
stopifnot(all(is.finite(ours[ok])), all(ours[!ok] == -Inf | ours[!ok] < -600))
res$reg_s402 <- max(abs(ours[ok] - ref[ok]))
cat(sprintf("reg_s402    %s (sha256 %s): max |log density - DDMreg.R| over %d rows (%d not underflowed) = %.2e\n",
            reg$kind, substr(reg$sha256, 1, 12), n, sum(ok), res$reg_s402))
cat(sprintf("\nlargest difference %.2e (limit 1e-10); implied score-bias error <= %.1e (limit 1e-6)\n",
            max(unlist(res)), max(unlist(res)) / 0.02))
stopifnot(max(unlist(res)) < 1e-10)
cat("ACCEPTED\n")
