# Cluster bundle (no run) for one of the NLE cells:
#   Rscript nle_phase6_bundle.R <cell name> <out dir> [reg_s402 card.rds] [trials] [replicates]
# e.g. reg_s402, whose control is the cap256w_c4 analytic DDM control (identical design and
# prior, checked here against cap256w_c4's cell), so it is bundled with run_control = FALSE.
suppressPackageStartupMessages(library(EMC2))
args <- commandArgs(TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "nle_phase6_cells.R"))
nm <- args[1]; dir <- args[2]
card <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL
trials <- if (length(args) >= 4) as.integer(args[4]) else 400L
reps <- if (length(args) >= 5) as.integer(args[5]) else 500L
cells <- nle_phase6_cells(card)
cell <- cells[[nm]]
if (nm == "reg_s402") {
  base <- cells$cap256w_c4
  stopifnot(identical(names(cell$mean), names(base$mean)), identical(cell$mean, base$mean),
            identical(cell$sd, base$sd), identical(cell$control_design$Ffactors, base$control_design$Ffactors),
            identical(sampled_pars(cell$control_design), sampled_pars(base$control_design)))
  cat("reg_s402 cell: design and prior identical to cap256w_c4's (its control is shared)\n")
}
nn_sbc_cell(cell, trials = trials, replicates = reps, run_control = FALSE, archive_dir = dir, run = FALSE,
            cores_per_chain = 8, cores_for_chains = 1,
            info = list(branch = "dev-nle", commit = system("git -C ~/Documents/Projects/EMCpackage/EMCrelease/EMC2-nle rev-parse --short HEAD", intern = TRUE),
                        why = paste("Phase 6 SERVER BATCH:", nm, "network cell, no local control (shares the cap256w_c4 DDM control)"),
                        where = "run on tux26/Gadi", topic = "NLE migration"))
cat("bundle written to", dir, "\n")
