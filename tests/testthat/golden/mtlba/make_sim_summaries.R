# Stage M0 fixture (d): simulator summaries, 1e6 trials each, from the reference simulators
#   rlba.DK                 (don't-know paper, ref_extract/dk_LBA-MTR_dists_extract.R)
#   rLBA_BE / rLBA_TC       (EMC2 branch BETE, ref_extract/bete_model_MTR_extract.R)
# Run from this directory:  Rscript make_sim_summaries.R    (needs msm; <= 8 cores; ~15 min)
#
# Substitutions for packages that are not installed here (nothing else is changed):
#  * rtnorm  (rlba.DK): msm::rtnorm, the function DMC used.
#  * rtmvnorm (rLBA_BE / rLBA_TC call tmvtnorm::rtmvnorm(n = 1, mean, sigma, lower = c(0, 0))):
#    the rejection sampler below -- draw from the untruncated bivariate normal until both
#    components are positive, which is also tmvtnorm's default algorithm ("rejection").
# Each summary also stores the EXPECTED response probabilities from the independent 1-D
# reference density (kernel_1d integrated over time) and the z scores of the simulated
# proportions against them.
library(parallel)
source("mtlba_ref_lib.R")
source("ref_extract/dk_LBA-MTR_dists_extract.R")
source("ref_extract/bete_model_MTR_extract.R")
rtnorm <- msm::rtnorm
rtmvnorm <- function(n = 1, mean, sigma, lower = rep(0, 2)) {
  stopifnot(n == 1, length(mean) == 2)
  s1 <- sqrt(sigma[1, 1]); s2 <- sqrt(sigma[2, 2]); r <- sigma[1, 2] / (s1 * s2)
  c2 <- sqrt(1 - r^2)
  repeat {
    z <- rnorm(2)
    x <- c(mean[1] + s1 * z[1], mean[2] + s2 * (r * z[1] + c2 * z[2]))
    if (x[1] > lower[1] && x[2] > lower[2]) return(matrix(x, 1))
  }
}
# One correction to the reference, for rLBA_BE only. When accumulator 1 LOSES, BETE bins its
# evidence with accumulator 2's threshold b as the top break (`pmats[2,!looser2,"b"]`), so
# with b1 > b2 any evidence between b2 and b1 gets rating NA. The top break should be the
# loser's own b. rLBA_BE_fixed is the verbatim source with that index changed from 2 to 1 (two
# places, done by text substitution below); it is identical to the original whenever
# b1 <= b2 (checked below on the `sym` set, same seed) and is what the BE summaries use.
local({
  src <- readLines("ref_extract/bete_model_MTR_extract.R")
  a <- grep("^# ---- rLBA_BE ", src); b <- grep("^# ---- rLBA_TC ", src)
  blk <- src[(a + 1):(b - 1)]
  hit <- grep('pmats[2,!looser2,"b"]', blk, fixed = TRUE)
  stopifnot(length(hit) == 2)
  blk <- sub('pmats[2,!looser2,"b"]', 'pmats[1,!looser2,"b"]', blk, fixed = TRUE)
  blk <- sub("^rLBA_BE <- ", "rLBA_BE_fixed <- ", blk)
  eval(parse(text = blk), envir = globalenv())
})
cores <- 8; n <- 1e6; chunk <- 1e5; probs <- c(.1, .3, .5, .7, .9)

summarise <- function(R, rt, levels) {
  R <- factor(R, levels = levels)
  list(n = length(rt), n_NA = sum(is.na(R)),
       p = as.vector(table(R)) / length(rt),
       mean_rt = as.vector(tapply(rt, R, mean)),
       q = t(sapply(levels, function(l) quantile(rt[R == l], probs, names = FALSE))),
       levels = levels, probs = probs)
}
# P(response) = integral over decision time of the kernel
p_kernel <- function(p, th, rho, sa, tt, lt, ut, ppos) {
  o <- 3 - sa
  f <- function(ts) vapply(ts, function(t)
    kernel_1d(t, p$A[sa], th[sa, tt], p$v[sa], p$sv[sa], p$A[o], th[o, lt], th[o, ut],
              p$v[o], p$sv[o], rho, ppos = ppos), 0)
  pw_int(f, c(0, .05, .1, .2, .4, .8, 1.6, 3.2, 6.4, 12.8, 50), rel.tol = 1e-9, abs.tol = 1e-14) +
    integrate(f, 50, Inf, rel.tol = 1e-7, abs.tol = 1e-13, stop.on.error = FALSE)$value
}
zscore <- function(p_sim, p_exp, n) (p_sim - p_exp) / sqrt(p_exp * (1 - p_exp) / n)

# ---- rlba.DK ---------------------------------------------------------------------------
dk_sets <- readRDS("dk_density.rds")$sets
dk <- mclapply(seq_along(dk_sets), function(i) {
  p <- dk_sets[[i]]; set.seed(7100 + i)
  sim <- do.call(rbind, lapply(seq_len(n / chunk), function(j)
    rlba.DK(chunk, A = p$A, b = p$b, d = p$d, t0 = p$t0, mean_v = p$v, sd_v = p$sv)))
  resp <- ifelse(sim$response2 == 2, "DK", as.character(sim$response))
  out <- summarise(resp, sim$rt, c("1", "2", "DK"))
  th <- cbind(DT0 = 0, DT1 = p$d, b = p$b); pp <- ppos_1d(p$v, p$sv, 0)
  out$p_expected <- c(p_kernel(p, th, 0, 1, 3, 1, 2, pp), p_kernel(p, th, 0, 2, 3, 1, 2, pp),
    p_kernel(p, th, 0, 1, 3, 2, 3, pp) + p_kernel(p, th, 0, 2, 3, 2, 3, pp))
  out$z <- zscore(out$p, out$p_expected, out$n)
  out$seed <- 7100 + i; out$pars <- p
  out
}, mc.cores = cores)
names(dk) <- names(dk_sets)
cat("rlba.DK done\n")

# ---- rLBA_BE / rLBA_TC ------------------------------------------------------------------
grid <- expand.grid(set = names(mtlba_sets), rho = c(0, -.7, .8), K = 2:4,
                    rule = c("BE", "TC"), stringsAsFactors = FALSE)
grid$seed <- 8200 + seq_len(nrow(grid))
run_bete <- function(i) {
  g <- grid[i, ]; p <- mtlba_sets[[g$set]]; K <- g$K; th <- mtlba_thresholds(p, K)
  # parameter matrix as BETE's Ttransform leaves it: two rows per trial (accumulator 1, 2)
  one <- cbind(v = p$v, sv = p$sv, b = th[, "b"], A = p$A, t0 = p$t0, r = g$rho,
               th[, paste0("DT", 1:(K - 1)), drop = FALSE], gp = 0)
  pars <- one[rep(1:2, chunk), ]
  rfun <- if (g$rule == "BE") rLBA_BE_fixed else rLBA_TC
  set.seed(g$seed)
  sim <- do.call(rbind, lapply(seq_len(n / chunk), function(j) rfun(lR = 1:(2 * K), pars)))
  Rn <- as.numeric(as.character(sim$R))
  out <- summarise(Rn, sim$rt, 1:(2 * K))
  dL <- get_dL(2 * K, g$rule); pp <- ppos_1d(p$v, p$sv, g$rho)
  out$p_expected <- vapply(seq_len(nrow(dL)), function(j) p_kernel(p, th, g$rho,
    dL[j, "sa"], dL[j, "tt"], dL[j, "lt"], dL[j, "ut"], pp), 0)
  out$z <- zscore(out$p, out$p_expected, out$n)
  # don't-know version (odd number of responses, 2K - 1): BETE's simulators only relabel,
  # levels(out$R) <- c(1:nr, nr:(2*nr-1)), i.e. responses K and K+1 become one response K.
  Rdk <- ifelse(Rn > K, Rn - 1, Rn)
  out$dont_know <- summarise(Rdk, sim$rt, 1:(2 * K - 1))
  out$dont_know$p_expected <- c(out$p_expected[seq_len(K - 1)], sum(out$p_expected[K + 0:1]),
                                out$p_expected[-seq_len(K + 1)])
  out$dont_know$z <- zscore(out$dont_know$p, out$dont_know$p_expected, out$n)
  out$config <- g; out$pars <- one
  out
}
bete <- mclapply(seq_len(nrow(grid)), run_bete, mc.cores = cores, mc.preschedule = FALSE)
names(bete) <- with(grid, paste(rule, K, set, rho, sep = "_"))
stopifnot(!sapply(bete, inherits, "try-error"), !sapply(dk, inherits, "try-error"))

# the corrected simulator equals the verbatim one when b1 <= b2; the verbatim one gives NA
# ratings when b1 > b2
fix_check <- local({
  f <- function(fun, set, K = 3, rho = .4) {
    p <- mtlba_sets[[set]]; th <- mtlba_thresholds(p, K)
    one <- cbind(v = p$v, sv = p$sv, b = th[, "b"], A = p$A, t0 = p$t0, r = rho,
                 th[, paste0("DT", 1:(K - 1)), drop = FALSE], gp = 0)
    set.seed(99); fun(lR = 1:(2 * K), one[rep(1:2, 2e4), ])
  }
  list(identical_when_b_equal = identical(f(rLBA_BE, "sym"), f(rLBA_BE_fixed, "sym")),
       identical_when_b1_lt_b2 = identical(f(rLBA_BE, "lowv"), f(rLBA_BE_fixed, "lowv")),
       verbatim_NA_when_b1_gt_b2 = sum(is.na(f(rLBA_BE, "asym")$R)),
       fixed_NA_when_b1_gt_b2 = sum(is.na(f(rLBA_BE_fixed, "asym")$R)))
})
print(fix_check)
stopifnot(fix_check$identical_when_b_equal, fix_check$identical_when_b1_lt_b2,
          fix_check$fixed_NA_when_b1_gt_b2 == 0)

zs <- c(unlist(lapply(dk, `[[`, "z")), unlist(lapply(bete, `[[`, "z")))
ps <- c(unlist(lapply(dk, `[[`, "p_expected")), unlist(lapply(bete, `[[`, "p_expected")))
po <- c(unlist(lapply(dk, `[[`, "p")), unlist(lapply(bete, `[[`, "p")))
big <- ps * n >= 10
checks <- list(
  n_NA_responses = c(dk = sum(sapply(dk, `[[`, "n_NA")), bete = sum(sapply(bete, `[[`, "n_NA"))),
  expected_mass_range = range(c(sapply(dk, function(x) sum(x$p_expected)),
                                sapply(bete, function(x) sum(x$p_expected)))),
  # z scores only where the expected count is at least 10 (normal approximation)
  n_z = sum(big), max_abs_z = max(abs(zs[big])), n_abs_z_gt_3 = sum(abs(zs[big]) > 3),
  z_quantiles = quantile(zs[big], c(.005, .025, .5, .975, .995)),
  z_sd = sd(zs[big]),
  small_cells = cbind(expected_n = ps[!big] * n, observed_n = po[!big] * n),
  max_abs_p_diff = max(abs(po - ps)),
  fix_check = fix_check)
print(checks)

saveRDS(list(
  dk = dk, bete = bete, grid = grid, checks = checks,
  provenance = list(
    what = "1e6-trial summaries per configuration: response proportions p, mean_rt, RT quantiles q (rows = responses, columns = probs); rt includes t0",
    responses = c(dk = "'1'/'2' = definitive response by that accumulator, 'DK' = don't know",
                  bete = "BETE response scale 1..2K (see bete_density.rds$dL); $dont_know = the 2K-1 response version"),
    sources = "ref_extract/ (verbatim); rtnorm = msm::rtnorm; rtmvnorm = rejection sampler in the generator; rLBA_BE with the loser's own b as the top rating break (see generator)",
    expected = "p_expected = kernel_1d (mtlba_ref_lib.R) integrated over decision time; z = (p - p_expected)/binomial SE",
    generator = "tests/testthat/golden/mtlba/make_sim_summaries.R",
    msm = as.character(packageVersion("msm")),
    generated = format(Sys.time()), R = R.version.string, RNGkind = RNGkind())),
  "sim_summaries.rds")
cat("wrote sim_summaries.rds\n")
