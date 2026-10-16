# The one verification command of the neural-likelihood layer: installs
# nothing and confirms that the installed EMC2 evaluates each card as its
# golden vectors say it should.
#
#   Rscript nle_verify.R [<card> ...]
#   Rscript -e 'source(system.file("scripts", "nle_verify.R", package = "EMC2"))'
#
# Without arguments the shipped artefacts (ddm_cap256w_c4, ddm_st0zero,
# rdm_small) are checked against the golden files shipped with them. Each
# further argument is a card file (<tag>_weights.json or .rds) checked against
# <tag>_golden.csv beside it; its twin is the card's `model` (DDMsc -> DDM).
# For a scaled DDM card the compiled likelihood is also compared with the R
# path over prior draws, and the R path with the analytic unit-scale DDM in
# place of the network (the oracle) where the NLE R reference is at hand.
#
# Prints one line per check, PASS or FAIL, and exits with status 1 on any
# FAIL. Run it after installing EMC2 and whenever a card is added.
suppressPackageStartupMessages(library(EMC2))
nle <- function(f) getFromNamespace(f, "EMC2")
args <- commandArgs(trailingOnly = TRUE)
fails <- 0L
say <- function(ok, ...) {
  cat(if (ok) "PASS " else "FAIL ", ..., "\n", sep = "")
  if (!ok) fails <<- fails + 1L
  invisible(ok)
}
try_check <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) { say(FALSE, label, ": ", conditionMessage(r)); return(invisible(NULL)) }
  r
}

# --- the environment -------------------------------------------------------------
lib <- find.package("EMC2")
sha_file <- file.path(dirname(lib), "EMC2_SHA.txt")
cat(sprintf("%s; EMC2 %s @ %s%s\n", R.version.string, as.character(utils::packageVersion("EMC2")), lib,
            if (file.exists(sha_file)) paste0(" (", readLines(sha_file, n = 1), ")") else ""))
algo <- nle("nle_hash_algo")()
cat(sprintf("artefact hashes: %s%s\n", algo,
            if (algo == "md5") " (this R has no tools::sha256sum; the MANIFEST's md5 lines pin the artefacts)" else ""))

# --- shipped artefacts -------------------------------------------------------------
for (a in c("ddm_cap256w_c4", "ddm_st0zero", "rdm_small")) {
  r <- try_check(a, nn_golden_check(a, flow_only = (a == "ddm_cap256w_c4"), error = FALSE))
  if (!is.null(r)) say(r$pass, sub("^(PASS|FAIL) ", "", format(r)))
}

# --- cards given on the command line -------------------------------------------------
scaled_checks <- function(card, m) {
  reg <- m()$nn
  # compiled vs R path through a design, 50 prior draws, N = 100, the full model
  PM <- c(v = 1, a = log(.8), t0 = log(.3), Z = qnorm(.5), sv = log(.1), SZ = qnorm(.2), st0 = log(.05))
  PS <- c(v = .15, a = .15, t0 = .1, Z = .05, sv = .15, SZ = .15, st0 = .15)
  cell <- suppressWarnings(nn_cell(m, mean = PM, sd = PS,
                                   formula = lapply(names(PM), function(p) stats::as.formula(paste(p, "~ 1"))),
                                   factors = list(S = c("left", "right")), Rlevels = c("left", "right")))
  set.seed(7)
  dat <- suppressMessages(make_data(cell$mean, cell$control_design, n_trials = 50))
  X <- nle("nn_cell_draws")(cell, 50)
  dadm <- nle("nn_dadm")(dat, cell$design)
  prepared <- cell$design$model()        # the design's prepared model, compiled likelihood off
  r_path <- function() { ml <- prepared; ml$c_name <- NULL; ml }
  native <- as.numeric(nle("calc_ll_manager")(X, dadm, cell$design$model))
  rp <- as.numeric(nle("calc_ll_manager")(X, dadm, r_path))
  d <- max(abs(native - rp))
  say(all(is.finite(native)) && d <= 1e-11,
      sprintf("compiled vs R path of '%s' over %d prior draws (N = 100): max |diff| %.2e against 1e-11",
              reg$label, nrow(X), d))
  # the oracle: the network replaced by the analytic unit-scale DDM (needs the
  # NLE reference ddm_sc_oracle; shipped with the tests, or given by EMC2_NLE_PORT)
  port_dir <- Sys.getenv("EMC2_NLE_PORT", "")
  port_file <- if (nzchar(port_dir)) file.path(port_dir, "ddm_scaled.R") else
    system.file("tests", "testthat", "port", "ddm_scaled.R", package = "EMC2")
  if (nzchar(port_file) && file.exists(port_file)) {
    port <- new.env(); sys.source(port_file, envir = port)
    set.seed(20261013); worst <- 0
    for (k in 1:60) {
      Z <- runif(1, .2, .8)
      p <- c(v = runif(1, -3, 3), a = runif(1, .4, 2.5), sv = if (k %% 3 == 0) 0 else runif(1, .05, 2),
             t0 = runif(1, .1, .5), st0 = 0, s = runif(1, .3, 2), Z = Z, SZ = if (k %% 4 == 0) 0 else runif(1, .05, .9))
      if (abs(p[["v"]] * p[["a"]] / p[["s"]]^2) > 9.5 || p[["sv"]] * p[["a"]] / p[["s"]]^2 > 7.5) next
      rt <- p[["t0"]] + (p[["a"]] / p[["s"]])^2 * rexp(8, 1 / .3)
      R <- factor(sample(1:2, length(rt), TRUE), levels = 1:2)
      P <- matrix(p, length(rt), length(p), byrow = TRUE, dimnames = list(NULL, names(p)))
      Pw <- P; Pw[, "SZ"] <- nn_sz_convert(P[, "SZ"], P[, "Z"], "proportion", "width")
      a <- nle("dDDM")(rt, R, Pw, precision = 1e-8)
      o <- nle("nn_eval_scaled")(reg, rt, as.integer(R), P, "pdf", evaluator = port$ddm_sc_oracle)
      big <- a > 1; if (any(big)) worst <- max(worst, abs(log(o[big]) - log(a[big])))
    }
    say(worst < 1e-9, sprintf("oracle (network replaced by the analytic unit-scale DDM) of '%s': max |dlog| %.2e against 1e-9",
                              reg$label, worst))
  } else cat("skip oracle check: the NLE reference ddm_scaled.R is not at hand (set EMC2_NLE_PORT)\n")
}

for (card in args) {
  if (!file.exists(card)) { say(FALSE, "no such card ", card); next }
  m <- try_check(card, register_nn_model(card))
  if (is.null(m)) next
  reg <- m()$nn
  cat(sprintf("card %s: kind %s, card_version %s, %s=%s\n", basename(card), reg$kind,
              format(reg$card_version), algo, reg$sha256))
  r <- try_check(card, nn_golden_check(m, error = FALSE))
  if (!is.null(r)) say(r$pass, sub("^(PASS|FAIL) ", "", format(r)))
  if (reg$kind == "flow_scaled") try_check(paste(card, "scaled checks"), scaled_checks(card, m))
}

cat(if (fails) sprintf("%d check(s) FAILED\n", fails) else "all checks passed\n")
quit(save = "no", status = if (fails) 1L else 0L)
