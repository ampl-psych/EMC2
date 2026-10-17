# PROVES: the same comparison over three seed pairs, all subjects, with on-on and off-off null pairs and a pooled 3-vs-3 comparison (run 08 first, from the same directory).
# PRINTED 8b257e3d: pooled 3v3 |z| max mu 0.85, sigma2 1.66, alpha 1.51 (none above 2); single pairs and the off-off null pair both reach |z| 2.6-3.9 on a few alphas; sd ratio 0.92-1.18; median ESS on/off alpha 8486/6262, sigma2 13004/10911.
suppressPackageStartupMessages(library(EMC2))
source_lines <- readLines("08_posterior_onoff.R")
eval(parse(text = source_lines[4:16]))   # data, design, sc(), one()
stats <- function(e, selection) {
  g <- get_pars(e, selection = selection, stage = "sample", return_mcmc = TRUE, merge_chains = FALSE, flatten = TRUE, by_subject = selection == "alpha")
  m <- do.call(cbind, lapply(g, function(x) as.matrix(do.call(rbind, x))))
  ess <- unlist(ess_summary(e, selection = selection, stage = "sample", stat = NULL, flatten = TRUE, by_subject = selection == "alpha"))
  data.frame(mean = colMeans(m), sd = apply(m, 2, sd), ess = as.numeric(ess)[seq_len(ncol(m))])
}
cmp <- function(a, b, what, label) {
  z <- (a$mean - b$mean) / sqrt(a$sd^2 / a$ess + b$sd^2 / b$ess); r <- a$sd / b$sd
  cat(sprintf("%-12s %-7s %3d params | |z| max %.2f, n>2: %d, mean z %+.2f | sd ratio %.2f-%.2f | median ESS %.0f / %.0f\n",
              label, what, length(z), max(abs(z)), sum(abs(z) > 2), mean(z), min(r), max(r), median(a$ess), median(b$ess)))
}
saved <- readRDS("posterior_onoff_emc.rds")
on11 <- saved$on; off11 <- saved$off
on12 <- one(TRUE, 12); off12 <- one(FALSE, 12); off13 <- one(FALSE, 13); on13 <- one(TRUE, 13)
for (sel in c("mu", "sigma2", "alpha")) {
  s <- lapply(list(on11 = on11, off11 = off11, on12 = on12, off12 = off12, on13 = on13, off13 = off13), stats, selection = sel)
  cmp(s$on11, s$off11, sel, "on11-off11"); cmp(s$on12, s$off12, sel, "on12-off12"); cmp(s$on13, s$off13, sel, "on13-off13")
  cmp(s$off11, s$off13, sel, "off11-off13 (null)"); cmp(s$on11, s$on13, sel, "on11-on13 (null)")
  # pooled: mean over the three ON runs vs mean over the three OFF runs
  pool <- function(l) data.frame(mean = rowMeans(sapply(l, `[[`, "mean")), sd = rowMeans(sapply(l, `[[`, "sd")), ess = rowSums(sapply(l, `[[`, "ess")))
  cmp(pool(s[c("on11","on12","on13")]), pool(s[c("off11","off12","off13")]), sel, "pooled 3v3")
}
