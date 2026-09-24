# Posterior-shift screens (60 data sets x 400 trials per S cell, 2000 draws) of the four NLE cells.
#   Rscript nle_phase6_screen.R <reg_s402 card.rds> [cores] [out.rds]
suppressPackageStartupMessages(library(EMC2))
args <- commandArgs(TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "nle_phase6_cells.R"))
cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
cat(sprintf("EMC2 %s @ %s; cores %d; BLAS %s\n", as.character(utils::packageVersion("EMC2")),
            find.package("EMC2"), cores, utils::sessionInfo()$BLAS))
cells <- nle_phase6_cells(args[1])
out <- list()
for (nm in names(cells)) {
  cat("\n====", nm, "====\n")
  out[[nm]] <- nn_posterior_shift(cells[[nm]], n_datasets = 60, n_trials = 400, n_draws = 2000, cores = cores)
  print(out[[nm]])
  if (length(args) >= 3) saveRDS(out, args[3])
}
