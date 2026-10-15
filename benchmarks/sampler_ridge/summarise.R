folder <- "benchmarks/sampler_ridge/results"
results <- setNames(lapply(c("current", "curvature", "directions"), function(name)
  readRDS(file.path(folder, paste0(name, ".rds")))),
  c("Current", "Curvature only", "Automatic directions"))

comparison <- do.call(rbind, lapply(names(results), function(name) {
  fits <- results[[name]]
  ridge <- do.call(rbind, lapply(fits, `[[`, "summary"))
  pars <- do.call(rbind, lapply(fits, `[[`, "diagnostics"))
  data.frame(method = name, seconds = median(vapply(fits, `[[`, 0, "seconds")),
    max_rhat = max(pars$rhat), ridge_ess = median(ridge$ess),
    ridge_ess_second = median(ridge$ess_second),
    converged = sum(vapply(fits, function(x) max(x$diagnostics$rhat) < 1.05, TRUE)))
}))
print(comparison, row.names = FALSE)
writeLines(capture.output(print(comparison, row.names = FALSE)),
           file.path(folder, "summary.txt"))

png(file.path(folder, "traces.png"), width = 1500, height = 600, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), las = 1)
col <- adjustcolor(c("#3366AA", "#DD7744", "#339966", "#AA55AA"), .65)
ylim <- range(results[[1]][[1]]$ridge[, , 1], results[[3]][[1]]$ridge[, , 1])
for (name in names(results)[c(1, 3)]) {
  matplot(results[[name]][[1]]$ridge[, , 1], type = "l", lty = 1, col = col,
          ylim = ylim, xlab = "Draw", ylab = "Population mean", main = name)
}
dev.off()

weak <- lapply(c("current", "directions"), function(name)
  readRDS(file.path(folder, paste0("weak-", name, ".rds"))))
print(do.call(rbind, lapply(weak, function(x) data.frame(
  method = x$summary$sampler[1], max_rhat = max(x$diagnostics$rhat),
  ridge_ess = min(x$summary$ess), seconds = x$seconds))), row.names = FALSE)
