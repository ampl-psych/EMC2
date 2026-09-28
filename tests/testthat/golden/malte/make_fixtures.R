# Builds the fixtures of tests/testthat/test-nn-malte.R from the NLE project
# (~/Downloads/julian/nle, stage M1): the three random-weight cards converted
# from Malte Lueken's flows by flow/scripts/malte_to_card.py, stored as .rds,
# and for every 10th row of their golden.csv (float64 values of his
# evaluate_pdf_sf) the values of the NLE R port (port/R/flow_race.R).
#   Rscript tests/testthat/golden/malte/make_fixtures.R ~/Downloads/julian/nle
nle <- path.expand(commandArgs(trailingOnly = TRUE)[1])
out <- "tests/testthat/golden/malte"
source(file.path(nle, "port/R/flow_race.R"))
for (nm in c("rdm_final", "rdm_plain", "crdm_final")) {
  card <- file.path(nle, "malte/fixtures", nm, "card.json")
  saveRDS(jsonlite::fromJSON(card, simplifyDataFrame = FALSE), file.path(out, paste0(nm, ".rds")), compress = "xz")
  fl <- flow_load(card)
  g <- utils::read.csv(file.path(nle, "malte/fixtures", nm, "golden.csv"))
  g <- g[seq(1, nrow(g), by = 10), ]
  ctx <- unlist(fl$context_names)
  g$port_log_pdf <- g$port_cdf <- NA_real_
  for (i in seq_len(nrow(g))) {
    ev <- flow_eval(fl, unlist(g[i, ctx]), g$dt[i])
    g$port_log_pdf[i] <- ev$log_pdf; g$port_cdf[i] <- ev$cdf
  }
  cat(sprintf("%s: %d rows, R port vs golden.csv max |log pdf diff| %.1e\n", nm, nrow(g),
              max(abs(g$port_log_pdf - g$log_pdf))))
  saveRDS(g[c(ctx, "dt", "log_pdf", "log_sf_exact", "port_log_pdf", "port_cdf")],
          file.path(out, paste0(nm, "_golden.rds")), compress = "xz")
}
