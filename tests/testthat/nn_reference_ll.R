# The compiled neural likelihood of this EMC2 on fixed data and parameter draws,
# in a form two installations can compare (tests/testthat/test-nn-reference.R).
#
# Shared by the test process and a subprocess running another EMC2 library:
#   Rscript -e 'source("nn_reference_ll.R")' <inputs.rds> <out.rds>
# <inputs.rds> is a list of cases, each with the data (`dat`), the proposal matrix
# (`P`) and the model's recipe (`model`, `formula`, `constants`, `factors`,
# `Rlevels`, and `card` for a card registered by path). The design and the dadm
# are built here, in the library under test, so what is compared is the
# likelihood code from the design onwards: sums, trial-wise values, and the
# multithreaded backend. Nothing random happens in this file.

nn_reference_cases <- function(card = NULL) {
  ddm <- function(model, formula, constants = c(s = 0), card = NULL)
    list(model = model, formula = formula, constants = constants, card = card,
         factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"))
  race <- function(levels)
    list(model = "RDMnn", formula = list(v ~ lR, B ~ 1, t0 ~ 1, A ~ 1), constants = c(s = log(1)),
         factors = list(subjects = 1, S = levels), Rlevels = levels)
  cases <- list(
    ddm_default = ddm("DDMnn", list(v ~ S, a ~ 1, t0 ~ 1, Z ~ 1, SZ ~ 1, sv ~ 1)),
    rdm_2acc   = race(c("a", "b")),
    rdm_3acc   = race(c("a", "b", "c")))
  if (!is.null(card))
    cases$card <- ddm("card", list(v ~ S, a ~ 1, t0 ~ 1, Z ~ 1, SZ ~ 1, sv ~ 1),
                      constants = NULL, card = card)
  cases
}

nn_reference_model <- function(case) {
  switch(case$model,
         DDMnn = EMC2::DDMnn,
         RDMnn = EMC2::RDMnn,
         card = EMC2::register_nn_model(case$card, twin = EMC2::DDM),
         stop("unknown model ", case$model))
}

nn_reference_design <- function(case)
  suppressMessages(EMC2::design(factors = case$factors, Rlevels = case$Rlevels,
                                model = nn_reference_model(case), formula = case$formula,
                                constants = case$constants, report_p_vector = FALSE))

# The true parameter vector each case's data are simulated from (test side).
nn_reference_p <- function(des, case) {
  p <- EMC2::sampled_pars(des, doMap = FALSE)
  if (case$model == "RDMnn") {
    p[] <- log(1)
    p[grepl("^v_lR", names(p))] <- log(1.8)
    p[c("B", "t0", "A")] <- log(c(1, .25, .3))
  } else {
    vals <- c(v = 1, v_Sb = -.5, a = log(1.2), t0 = log(.25), Z = 0, SZ = qnorm(.3), sv = log(.5))
    p[] <- vals[names(p)]
  }
  p
}

# Likelihoods of one case in the loaded EMC2: sums, trial-wise, both backends.
nn_reference_ll <- function(case) {
  des <- nn_reference_design(case)
  emc <- suppressMessages(EMC2::make_emc(case$dat, des, type = "single", n_chains = 1,
                                         compress = FALSE, verbose = FALSE, rt_resolution = .001))
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model
  manager <- get("calc_ll_manager", asNamespace("EMC2"))
  old <- options(emc.ll_backend = "multiprocess")
  on.exit(options(old))
  out <- list(sum = as.numeric(manager(case$P, dadm, model)),
              trialwise = manager(case$P, dadm, model, return_trialwise = TRUE))
  options(emc.ll_backend = "multithreaded", emc.n_threads = 2)
  out$sum_mt <- as.numeric(manager(case$P, dadm, model))
  out$trialwise_mt <- manager(case$P, dadm, model, return_trialwise = TRUE)
  out$n_rows <- nrow(dadm)
  out
}

nn_reference_run <- function(inputs) {
  res <- lapply(inputs, nn_reference_ll)
  attr(res, "library") <- find.package("EMC2")
  attr(res, "version") <- as.character(utils::packageVersion("EMC2"))
  res
}

# As a subprocess (the test sets NN_REFERENCE_MAIN=1): <inputs.rds> <out.rds>
if (identical(Sys.getenv("NN_REFERENCE_MAIN"), "1")) {
  a <- commandArgs(TRUE)
  if (length(a) != 2L) stop("usage: Rscript -e 'source(\"nn_reference_ll.R\")' <inputs.rds> <out.rds>")
  saveRDS(nn_reference_run(readRDS(a[1])), a[2])
}
