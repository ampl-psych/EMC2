# Stage M0 fixture (a): the don't-know paper's n1PDF.DK on a grid. EXACT PARITY TARGET for the
# uncorrelated MTLBA with K = 2 and dont_know = TRUE.
#   Run from this directory:  Rscript make_dk_density.R      (needs no packages)
# Reference code: ref_extract/dk_LBA-MTR_dists_extract.R (verbatim, see extract_refs.R).
# Only cells with the don't-know threshold d >= A (plan decision: all thresholds above A).
# Drift convention: independent positive-truncated normals (posdrift = TRUE), which is the
# both-positive convention at rho = 0.
source("mtlba_ref_lib.R")
source("ref_extract/dk_LBA-MTR_dists_extract.R")
new <- environment(); old <- new.env()
sys.source("ref_extract/dk_MTR-LBA_dists_extract.R", envir = old)  # older d >= A only version
environment(old$n1PDF.DK) <- environment(old$n1PDFfixedt0.DK) <- old
for (f in c("dlba.norm", "plba.norm", "dlba.norm.core", "plba.norm.core"))
  assign(f, get(f), envir = old)

# accumulator order (1, 2); internal DMC parameters A, b, d, mean_v, sd_v, t0
dk_sets <- list(
  paper_check = list(A = c(.5, .5), b = c(2, 2),     d = c(1.5, 1.5), v = c(2, 1),
                     sv = c(1, 1),    t0 = .2),   # the check example in the DMC file
  d_equals_A  = list(A = c(1, 1),   b = c(2, 2),     d = c(1, 1),     v = c(1, 0),
                     sv = c(1, 1),    t0 = .1),   # boundary case d = A
  asym        = list(A = c(.4, .8), b = c(1.6, 1.5), d = c(.9, 1.2),  v = c(3, 1),
                     sv = c(.7, 1.2), t0 = .15),
  qD          = list(A = c(.6, .6), b = c(1.5, 2),   d = plogis(c(.5, 1.5)) * c(1.5, 2),
                     v = c(1.2, 2.2), sv = c(1, .8), t0 = .25),  # lba_DK_qD.R: d = plogis(qD)*b
  neg_mean    = list(A = c(1, .3),  b = c(1.5, 1.8), d = c(1.2, .9),  v = c(.8, -.3),
                     sv = c(1, .6),   t0 = .3)
)
dt_grid <- c(.05, .1, .2, .35, .5, .75, 1, 1.5, 2.5, 4)

dens <- function(fun, p, rt, resp) {
  o <- if (resp == "2") 2:1 else 1:2   # n1 order: the winning accumulator first
  fun(rt, r2 = if (resp == "DK") "DK" else "", A = p$A[o], b = p$b[o], d = p$d[o],
      mean_v = p$v[o], sd_v = p$sv[o], t0 = p$t0)
}
cells <- do.call(rbind, lapply(names(dk_sets), function(nm) {
  p <- dk_sets[[nm]]; stopifnot(p$d >= p$A, p$d < p$b)
  do.call(rbind, lapply(c("1", "2", "DK"), function(resp) {
    rt <- p$t0 + dt_grid
    x <- data.frame(set = nm, resp = resp, dt = dt_grid, rt = rt,
      dens = dens(n1PDF.DK, p, rt, resp),
      dens_old_version = dens(old$n1PDF.DK, p, rt, resp),
      # the DK density must not depend on which accumulator is listed first
      dens_swapped = if (resp == "DK") {
        n1PDF.DK(rt, r2 = "DK", A = p$A[2:1], b = p$b[2:1], d = p$d[2:1],
                 mean_v = p$v[2:1], sd_v = p$sv[2:1], t0 = p$t0) } else NA_real_)
    # independent check: the plan's kernel at rho = 0 (mtlba_ref_lib.R)
    k <- function(s, dl, du) vapply(dt_grid, function(t) { o <- 3 - s
      kernel_1d(t, p$A[s], p$b[s], p$v[s], p$sv[s], p$A[o], dl, du, p$v[o], p$sv[o], 0) }, 0)
    x$kernel_ref <- switch(resp, "1" = k(1, 0, p$d[2]), "2" = k(2, 0, p$d[1]),
                           "DK" = k(1, p$d[2], p$b[2]) + k(2, p$d[1], p$b[1]))
    x
  }))
}))
# EMC2 coding (plan F3): R = winning accumulator, RR = 2 definitive / 1 don't know;
# don't-know trials carry the placeholder R = first response level.
cells$R  <- ifelse(cells$resp == "2", 2L, 1L)
cells$RR <- ifelse(cells$resp == "DK", 1L, 2L)

# Far-tail cells (dens <= 1e-8): the DMC closed form loses relative accuracy there (differences
# of normal CDFs; one DK cell is slightly negative), so they are parity targets in ABSOLUTE
# terms only. Flagged, kept.
cells$tail <- cells$dens <= 1e-8

# total probability over the three responses, per set
mass <- sapply(dk_sets, function(p) {
  f <- function(resp) integrate(function(dt) dens(n1PDF.DK, p, dt + p$t0, resp), 0, Inf,
                                rel.tol = 1e-10)$value
  c(p1 = f("1"), p2 = f("2"), pDK = f("DK"))
})
mass <- rbind(mass, total = colSums(mass))

checks <- list(
  max_abs_new_vs_old_version = max(abs(cells$dens - cells$dens_old_version)),
  max_abs_DK_order_swap = max(abs(cells$dens - cells$dens_swapped), na.rm = TRUE),
  max_abs_vs_kernel_ref = max(abs(cells$dens - cells$kernel_ref)),
  max_rel_vs_kernel_ref_not_tail = with(cells[cells$dens > 1e-8, ],
    max(abs(dens - kernel_ref) / dens)),
  n_tail_cells = sum(cells$dens <= 1e-8),
  mass = mass)
print(checks)
stopifnot(checks$max_abs_new_vs_old_version < 1e-14, checks$max_abs_DK_order_swap < 1e-14,
          checks$max_abs_vs_kernel_ref < 1e-12, abs(mass["total", ] - 1) < 1e-8)

saveRDS(list(
  cells = cells[, c("set", "resp", "R", "RR", "dt", "rt", "dens", "kernel_ref", "tail")],
  sets = dk_sets, checks = checks,
  provenance = list(
    what = "n1PDF.DK(t, r2, A, b, d, mean_v, sd_v, t0, posdrift = TRUE) on a grid",
    source = "DMC_DK_OSF/rod/dmc/models/LBA-MTR/dists.R (https://osf.io/download/akzxw/), verbatim extract in ref_extract/",
    columns = c(resp = "'1'/'2' = that accumulator wins with a definitive response, 'DK' = don't know (either wins, loser above d)",
                R = "EMC2 coding: winning accumulator; placeholder 1 on DK rows",
                RR = "EMC2 coding: 2 = definitive, 1 = don't know",
                dens = "the parity target (DK value = SUM of the two terms)",
                kernel_ref = "mtlba_ref_lib.R kernel_1d at rho = 0, independent check",
                tail = "dens <= 1e-8: compare in absolute terms only (DMC closed form not relatively accurate there)"),
    generator = "tests/testthat/golden/mtlba/make_dk_density.R",
    generated = format(Sys.time()), R = R.version.string)),
  "dk_density.rds")
cat("wrote dk_density.rds:", nrow(cells), "cells\n")
