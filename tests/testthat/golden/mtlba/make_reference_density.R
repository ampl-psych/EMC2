# Stage M0 fixture (c): independent high-accuracy references for the kernel, and the check
# that BETE's closed form (fixture b) agrees with them.
#   Run from this directory, after make_bete_density.R:
#   Rscript make_reference_density.R      (base R only; <= 8 cores; mvtnorm only for the
#                                          normaliser comparison)
# For every cell of bete_density.rds:
#   dens_1d = integrate() on the one-dimensional form (kernel_1d, closed-form band)
#   dens_2d = brute-force two-dimensional integral of the joint drift density (kernel_2d,
#             no conditional-normal algebra, non-negative integrand)
# both with integrate(rel.tol = 1e-13) on pieces split at the integrand's kinks and peaks.
library(parallel)
source("mtlba_ref_lib.R")
cores <- 8
b <- readRDS("bete_density.rds")
cells <- b$cells

# normaliser P(both drifts > 0) per (set, rho): 1-D, 2-D, and mvtnorm
pp <- unique(cells[, c("set", "rho")])
pp$ppos_1d <- pp$ppos_2d <- pp$ppos_pmvnorm <- NA_real_
for (i in seq_len(nrow(pp))) {
  p <- mtlba_sets[[pp$set[i]]]; rho <- pp$rho[i]
  pp$ppos_1d[i] <- ppos_1d(p$v, p$sv, rho)
  pp$ppos_2d[i] <- ppos_2d(p$v, p$sv, rho)
  S <- matrix(c(p$sv[1]^2, rep(rho * prod(p$sv), 2), p$sv[2]^2), 2)
  pp$ppos_pmvnorm[i] <- mvtnorm::pmvnorm(lower = c(0, 0), upper = c(Inf, Inf), mean = p$v,
                                         sigma = S)
}

# the distinct kernels (BE and TC, and different K, share many)
key <- c("set", "rho", "sa", "bs", "dl", "du", "dt")
u <- unique(cells[, key])
res <- mclapply(seq_len(nrow(u)), function(i) {
  x <- u[i, ]; p <- mtlba_sets[[x$set]]; s <- x$sa; o <- 3 - s
  j <- which(pp$set == x$set & pp$rho == x$rho)
  c(kernel_1d(x$dt, p$A[s], x$bs, p$v[s], p$sv[s], p$A[o], x$dl, x$du, p$v[o], p$sv[o],
              x$rho, ppos = pp$ppos_1d[j]),
    kernel_2d(x$dt, p$A[s], x$bs, p$v[s], p$sv[s], p$A[o], x$dl, x$du, p$v[o], p$sv[o],
              x$rho, ppos = pp$ppos_2d[j]))
}, mc.cores = cores)
stopifnot(all(lengths(res) == 2))
u$dens_1d <- sapply(res, `[`, 1); u$dens_2d <- sapply(res, `[`, 2)
m <- merge(cbind(cells, .ord = seq_len(nrow(cells))), u, by = key, sort = FALSE)
m <- m[order(m$.ord), ]; m$.ord <- NULL
stopifnot(nrow(m) == nrow(cells), all(m$dens == cells$dens))
names(m)[names(m) == "dens"] <- "dens_bete"
m <- m[, c(setdiff(names(cells), "dens"), "dens_bete", "dens_1d", "dens_2d")]
rownames(m) <- NULL

agree <- function(a, b, sub) c(max_abs = max(abs(a - b)[sub]),
                               max_rel = max((abs(a - b) / pmax(abs(b), 1e-300))[sub]))
body <- m$dens_2d > 1e-5           # where a relative comparison is meaningful in double
agreement <- list(
  n_cells = nrow(m), n_distinct_kernels = nrow(u), n_body = sum(body),
  all_cells = rbind(`1d_vs_2d` = agree(m$dens_1d, m$dens_2d, TRUE),
                    bete_vs_2d = agree(m$dens_bete, m$dens_2d, TRUE),
                    bete_vs_1d = agree(m$dens_bete, m$dens_1d, TRUE))[, "max_abs", drop = FALSE],
  body_cells_dens_gt_1e5 = rbind(`1d_vs_2d` = agree(m$dens_1d, m$dens_2d, body),
                    bete_vs_2d = agree(m$dens_bete, m$dens_2d, body),
                    bete_vs_1d = agree(m$dens_bete, m$dens_1d, body)),
  by_rho_zero = rbind(rho0 = agree(m$dens_bete, m$dens_2d, body & m$rho == 0),
                      rho_non0 = agree(m$dens_bete, m$dens_2d, body & m$rho != 0)),
  ppos = c(max_abs_1d_vs_2d = max(abs(pp$ppos_1d - pp$ppos_2d)),
           max_abs_pmvnorm_vs_2d = max(abs(pp$ppos_pmvnorm - pp$ppos_2d))),
  min_dens = c(bete = min(m$dens_bete), `1d` = min(m$dens_1d), `2d` = min(m$dens_2d)))
print(agreement)
# the stage's requirement: (b) and (c) agree, and the two references agree to 1e-10
stopifnot(agreement$all_cells[, "max_abs"] < 1e-10,
          agreement$body_cells_dens_gt_1e5[, "max_rel"] < 1e-9,
          agreement$ppos < 1e-10, m$dens_2d >= 0)

saveRDS(list(
  cells = m, ppos = pp, agreement = agreement, sets = mtlba_sets,
  provenance = list(
    what = "kernel density by integrate(): dens_1d (1-D form, closed-form band) and dens_2d (brute-force 2-D); dens_bete copied from bete_density.rds",
    use = c("dens_2d is the arbiter: non-negative integrand, no cancellation, valid in the far tail",
            "dens_1d and dens_bete are accurate in ABSOLUTE terms (~1e-14) but not relatively in the far tail (band = difference of two near-equal numbers); dens_1d can be slightly negative there"),
    generator = "tests/testthat/golden/mtlba/make_reference_density.R + mtlba_ref_lib.R",
    generated = format(Sys.time()), R = R.version.string)),
  "ref_density.rds")
cat("wrote ref_density.rds\n")
