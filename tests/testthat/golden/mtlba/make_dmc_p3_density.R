# Generator of golden/mtlba/dmc_p3_density.rds (Stage M7c, 2026-10-07): the density of the
# manuscript's DMC model MTRp3 (drifts bivariate normal conditioned on not both negative, the
# lowest confidence threshold d1 free below A, per-accumulator confidence guessing) from the
# manuscript's own DMC code (rating-work/mtlba/stageM4/dmc_copy, model MTRp3_censor.R, dists.R
# functions n1PDFMTRp3 / n1PDFMTRp3.CENSOR), with the 5-point Simpson rule over the winner's
# start point replaced by 4 panels x 16 Gauss-Legendre nodes (stageM4/dists_quad.R, the
# accurate rule of the M4 accuracy map) so the values are accurate to ~1e-12. Censored
# masses are integrate()d at rel.tol 1e-10 over the sum of the six response densities.
# Parameter sets: the posterior means of DMC's fits of r94 subject 1 (model 5, no guessing;
# model 1, guessing) rounded to 4 digits, and an uncorrelated set with d1 below A on both
# accumulators. Rows: accumulators NEW, OLD; the response is the first accumulator of the
# stated order ("n1" = NEW responds, "n2" = OLD responds).
M4 <- "/Users/andrew/Documents/Projects/EMCpackage/rating-work/mtlba/stageM4"
source(file.path(M4, "dmc_setup.R")); load_dmc_model("MTRp3_censor.R"); source(file.path(M4, "dists_quad.R"))
gl <- function(n) { i <- seq_len(n - 1); b <- i / sqrt(4 * i^2 - 1); J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b; e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2) }
g <- gl(16); M4QUAD <<- list(x = as.vector(outer((g$x + 1) / 2, 0:3, "+")) / 4, w = rep(g$w, 4) / 8)
sets <- list(
  m5 = list(A = c(.3037, .4577), b = c(.8069, .8711), fb = c(.3835, .5233), pf = c(0, 0), p2 = c(.5, .5), p3 = c(.5, .5),
            d1 = c(.1255, .0484), d2 = c(.4451, .6708), t0 = .3251, mean_v = c(1.3327, -.3391), sd_v = c(1, 1.2),
            corr_v = -.8144),
  m1 = list(A = c(.2490, .4678), b = c(.6909, .8370), fb = c(.3197, .5270), pf = c(.0013, .0432), p2 = c(.0143, .0143),
            p3 = c(.9805, .9805), d1 = c(.0247, .0436), d2 = c(.3293, .6189), t0 = .3572, mean_v = c(1.2618, -.4137),
            sd_v = c(1, 1), corr_v = -.8403),
  u  = list(A = c(.5, .4), b = c(1.4, 1.2), fb = c(.9, .7), pf = c(.3, .1), p2 = c(.4, .4), p3 = c(.5, .5),
            d1 = c(.3, .2), d2 = c(.8, .7), t0 = .2, mean_v = c(1.2, .3), sd_v = c(.8, 1.1), corr_v = 0))
rts <- c(.42, .5, .6, .75, .9, 1.2, 1.8, 3)
cells <- NULL
for (sn in names(sets)) for (ord in c("n1", "n2")) for (r2 in 1:3) {
  s <- sets[[sn]]
  o <- if (ord == "n1") 1:2 else 2:1
  dens <- n1PDFMTRp3(rts, rep(r2, length(rts)), s$A[o], s$b[o], s$fb[o], s$pf[o], s$p2[o], s$p3[o], s$d1[o], s$d2[o],
                     s$mean_v[o], s$sd_v[o], s$corr_v, s$t0)
  cells <- rbind(cells, data.frame(set = sn, order = ord, r2 = r2, rt = rts, dens = dens))
}
# censored masses: P(t0 < rt < CF) and P(rt > CS)
dmc_sum <- function(x, s) {
  out <- matrix(ncol = 6, nrow = length(x)); t <- rep(x, times = 3); r2 <- rep(1:3, each = length(x))
  out[, 1:3] <- n1PDFMTRp3(t, r2, s$A, s$b, s$fb, s$pf, s$p2, s$p3, s$d1, s$d2, s$mean_v, s$sd_v, s$corr_v, s$t0)
  o <- 2:1
  out[, 4:6] <- n1PDFMTRp3(t, r2, s$A[o], s$b[o], s$fb[o], s$pf[o], s$p2[o], s$p3[o], s$d1[o], s$d2[o], s$mean_v[o], s$sd_v[o], s$corr_v, s$t0)
  apply(out, 1, sum)
}
cens <- NULL
for (sn in names(sets)) {
  s <- sets[[sn]]
  for (CF in c(.45, .55)) cens <- rbind(cens, data.frame(set = sn, kind = "fast", bound = CF,
    mass = integrate(dmc_sum, s$t0, CF, s = s, rel.tol = 1e-10, subdivisions = 2000L)$value))
  for (CS in c(1.5, 2.5)) cens <- rbind(cens, data.frame(set = sn, kind = "slow", bound = CS,
    mass = integrate(dmc_sum, CS, Inf, s = s, rel.tol = 1e-10, subdivisions = 2000L)$value))
}
saveRDS(list(sets = sets, cells = cells, cens = cens,
             note = "DMC MTRp3 (manuscript code, 4 x 16 Gauss-Legendre start-point rule); see make_dmc_p3_density.R"),
        "/Users/andrew/Documents/Projects/EMCpackage/EMCrelease/EMC2-rating/tests/testthat/golden/mtlba/dmc_p3_density.rds")
print(summary(cells$dens)); print(cens)
