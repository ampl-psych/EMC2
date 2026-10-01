# Stage M7 fixture (e): BETE's likelihood with confidence guessing (log_likelihood_mt, branch
# origin/BETE 9eab6fd0, R/likelihood.R lines 393-415) on the Stage M0 grid at rho = 0.
#   Run from this directory, after make_bete_density.R:  Rscript make_bete_guess_density.R
#
# BETE's guessing is a mixture over TRIALS: with probability gp both accumulators race to
# their guess thresholds gb (a plain positive-drift LBA race, under either rule), and the
# rating is drawn from the pg probabilities; otherwise the model without guessing applies:
#   L = (1 - gp) * ll + gp * LBAll * pR
#   ll    = the response density without guessing (fixture (b), bete_density.rds)
#   LBAll = dLBA(rt; winner, b = gb) * (1 - pLBA(rt; loser, b = gb)), winner = the choice
#   pR    = get_pgm(pg) (last rating = the rest), rating from the response by
#           rate <- c(K:1, 1:K)  (even number of responses) or c(K:1, 2:K) (odd, don't know)
#   don't know (odd): ll and LBAll are the AVERAGE of the two middle terms.
# The guess race is evaluated here with the don't-know paper's dlba.norm / plba.norm (verbatim
# extract; the standard positive-drift LBA formulas, as EMC2's LBA() that BETE calls), and the
# rating probabilities with BETE's get_pgm (verbatim extract). gb = A + gB * B (BETE's *PLBA
# models). BETE's guess race ignores the drift correlation r, so only rho = 0 is a parity target.
source("mtlba_ref_lib.R")
source("ref_extract/dk_LBA-MTR_dists_extract.R")
source("ref_extract/bete_model_MTR_extract.R")
base <- readRDS("bete_density.rds")$cells
base <- base[base$rho == 0, ]

# guess settings: gp, gB (per accumulator), pg1 .. pg3 (the first K - 1 are used)
configs <- list(
  mid   = list(gp = .1, gB = c(.5, .5), pg = c(.5, .3, .4)),
  asym  = list(gp = .2, gB = c(.3, .8), pg = c(.2, .6, .1)),
  low   = list(gp = .15, gB = c(.6, .4), pg = c(1, 0, 0)))   # BETE's default: always guess low

lba_race <- function(dt, p, w, gb) {
  l <- 3 - w
  dlba.norm(dt, A = p$A[w], b = gb[w], mean_v = p$v[w], sd_v = p$sv[w], posdrift = TRUE) *
    (1 - plba.norm(dt, A = p$A[l], b = gb[l], mean_v = p$v[l], sd_v = p$sv[l], posdrift = TRUE))
}

cells <- list()
for (set in names(mtlba_sets)) for (K in 1:4) for (rule in c("BE", "TC"))
  for (cf in names(configs)) for (dk in c(FALSE, if (K > 1) TRUE)) {
  p <- mtlba_sets[[set]]; g <- configs[[cf]]
  gb <- p$A + g$gB * p$B
  b0 <- base[base$set == set & base$K == K & base$rule == rule, ]
  pm <- if (K > 1) get_pgm(matrix(g$pg[seq_len(K - 1)], nrow = 1)) else matrix(0, 1, 0)
  pm <- cbind(pm, 1 - sum(pm))
  nr <- if (dk) 2 * K - 1 else 2 * K
  rate <- if (dk) c(K:1, 2:K) else c(K:1, 1:K)
  for (R in seq_len(nr)) for (dt in mtlba_dt) {
    std <- function(Reven) b0$dens[b0$R == Reven & b0$dt == dt]
    if (dk && R == K) {            # the middle (don't-know) response: average of the two
      ll <- (std(K) + std(K + 1)) / 2
      lba <- lba_race(dt, p, 1, gb) / 2 + lba_race(dt, p, 2, gb) / 2
    } else {
      Reven <- if (dk && R > K) R + 1 else R
      ll <- std(Reven)
      lba <- lba_race(dt, p, if (Reven <= K) 1 else 2, gb)
    }
    pR <- pm[1, rate[R]]
    cells[[length(cells) + 1]] <- data.frame(
      set = set, K = K, rule = rule, config = cf, dk = dk, R = R, dt = dt, rt = p$t0 + dt,
      std = ll, guess = lba, pR = pR, dens = (1 - g$gp) * ll + g$gp * lba * pR)
  }
}
cells <- do.call(rbind, cells)
stopifnot(!anyNA(cells$dens), all(cells$dens >= 0))

saveRDS(list(
  cells = cells, sets = mtlba_sets, configs = configs,
  provenance = list(
    what = "BETE log_likelihood_mt with guessing, per response and time: (1 - gp) * ll + gp * LBAll * pR at rho = 0",
    source = "EMC2 branch origin/BETE (9eab6fd0), R/likelihood.R lines 393-415 and R/model_MTR.R (get_pgm); copy in rating-work/mtlba/ref/bete",
    notes = c(
      "trial-level mixture: MTLBA(guess = 'trial'); the manuscript's per-accumulator trigger model is NOT this",
      "R = BETE's folded response 1 .. 2K (dk = FALSE) or 1 .. 2K - 1 (dk = TRUE, R = K the don't-know response)",
      "don't know: BETE averages the two middle terms; MTLBA sums them, so MTLBA = 2 * dens on those rows",
      "gb = A + gB * B per accumulator; pg stick-breaking via get_pgm, pg1 = the LOWEST rating",
      "BETE reads gp and the pg from the SECOND accumulator's row only (likelihood and simulator), so gp ~ lR does not give a trigger probability per accumulator there; here they are equal for both accumulators",
      "std = the no-guess density from bete_density.rds; guess = the LBA race to gb (dlba.norm / plba.norm, DK-paper extract)"),
    generator = "tests/testthat/golden/mtlba/make_bete_guess_density.R",
    generated = format(Sys.time()), R = R.version.string)),
  "bete_guess_density.rds")
cat("wrote bete_guess_density.rds:", nrow(cells), "cells\n")
print(summary(cells$dens))
