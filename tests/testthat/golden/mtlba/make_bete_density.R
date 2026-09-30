# Stage M0 fixture (b): BETE's n1PDF_MTR_1 on a grid, balance of evidence (BE) and threshold
# count (TC), K = 1..4 ratings per choice, rho = 0 and rho != 0. Both-positive drift convention.
#   Run from this directory:  Rscript make_bete_density.R     (needs mvtnorm; <= 8 cores)
# Reference code: ref_extract/bete_model_MTR_extract.R (verbatim, see extract_refs.R).
#
# One row per (parameter set, rho, K, rule, BETE response R in 1..2K, dt). For each response
# BETE's get_dL gives the kernel: stopping accumulator `sa`, its stopping threshold index
# `tt`, the other accumulator `nsa` and its band (`lt`, `ut`) -- indices into the threshold
# matrix c(DT0 = 0, DT1, .., b). The call is exactly the one in BETE's log_likelihood_mt:
#   n1PDF_MTR_1(rt, pars[c(sa, nsa), ], dl = th[nsa, lt], du = th[nsa, ut], b = th[sa, tt])
# BETE's response scale R = 1..2K runs from high-confidence accumulator 1 to high-confidence
# accumulator 2; translating it to EMC2's (R, RR) is Stage M1's job (the response lookup).
library(mvtnorm); library(parallel)
source("mtlba_ref_lib.R")
source("ref_extract/bete_model_MTR_extract.R")
cores <- 8

grid <- expand.grid(set = names(mtlba_sets), rho = mtlba_rho, K = 1:4, rule = c("BE", "TC"),
                    stringsAsFactors = FALSE)
cells <- do.call(rbind, mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; p <- mtlba_sets[[g$set]]; th <- mtlba_thresholds(p, g$K)
  dL <- get_dL(2 * g$K, g$rule)
  pars <- cbind(A = p$A, t0 = p$t0, v = p$v, sv = p$sv, r = g$rho)
  do.call(rbind, lapply(seq_len(nrow(dL)), function(j) {
    d <- dL[j, ]; sa <- d[["sa"]]; nsa <- d[["nsa"]]
    bs <- th[sa, d[["tt"]]]; dl <- th[nsa, d[["lt"]]]; du <- th[nsa, d[["ut"]]]
    rt <- p$t0 + mtlba_dt
    data.frame(g, R = d[["RR"]], sa = sa, tt = d[["tt"]], nsa = nsa, lt = d[["lt"]],
               ut = d[["ut"]], bs = bs, dl = dl, du = du, dt = mtlba_dt, rt = rt,
               dens = vapply(rt, function(x)
                 n1PDF_MTR_1(x, pars[c(sa, nsa), ], dl = dl, du = du, b = bs), 0),
               row.names = NULL)
  }))
}, mc.cores = cores))
stopifnot(!anyNA(cells$dens), nrow(cells) == nrow(grid) / 4 * sum(2 * 1:4) * length(mtlba_dt))

# lookup tables, even (2K) and odd (2K-1, don't know = the two middle responses merged)
dL <- list()
for (rule in c("BE", "TC")) for (NR in 2:8) dL[[paste0(rule, NR)]] <- get_dL(NR, rule)
thresholds <- lapply(mtlba_sets, function(p)
  setNames(lapply(1:4, function(K) mtlba_thresholds(p, K)), paste0("K", 1:4)))

saveRDS(list(
  cells = cells, sets = mtlba_sets, thresholds = thresholds, dL = dL,
  provenance = list(
    what = "n1PDF_MTR_1(rt, pars, dl, du, b): BETE closed form from truncated bivariate normal moments (recintab0 / mvtnorm::pmvnorm)",
    source = "EMC2 branch origin/BETE (9eab6fd0), R/model_MTR.R, verbatim extract in ref_extract/",
    notes = c(
      "pars columns A, t0, v, sv, r with r already on the (-1, 1) correlation scale",
      "thresholds: d_k = A + B * (k/K)^pow, b = A + B (sets in mtlba_ref_lib.R); all above A",
      "don't know (odd number of responses): BETE's likelihood AVERAGES the two middle terms; the plan's likelihood is their SUM. Only the separate terms are stored here (rows R = K and R = K + 1 of the even table).",
      "absolute accuracy of BETE vs the independent references is ~1e-14 (see ref_density.rds); far-tail values are not relatively accurate"),
    generator = "tests/testthat/golden/mtlba/make_bete_density.R",
    mvtnorm = as.character(packageVersion("mvtnorm")),
    generated = format(Sys.time()), R = R.version.string)),
  "bete_density.rds")
cat("wrote bete_density.rds:", nrow(cells), "cells\n")
print(summary(cells$dens))
