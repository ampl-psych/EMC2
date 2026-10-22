# The compiled neural likelihood of this EMC2 against another installed EMC2
# (the reference, e.g. dev-nle-ondev) on the same data and parameter draws:
# every shipped card, and a card given by path. Skipped unless
#   EMC2_NLE_REF_LIB = the library holding the reference EMC2 (put first on R_LIBS
#                      of a subprocess; its dependencies come from this session's
#                      .libPaths())
#   EMC2_NLE_CARD    = a card (<tag>_weights.json or .rds) for the fifth case;
#                      without it the shipped cards alone are compared.
# Sums, trial-wise values and the multithreaded backend must agree to 1e-12
# (absolute, log-likelihood units). The data and the proposals are drawn once
# here; each library builds the design and dadm itself (nn_reference_ll.R).

ref_lib <- Sys.getenv("EMC2_NLE_REF_LIB")
card <- Sys.getenv("EMC2_NLE_CARD")
have_ref <- nzchar(ref_lib) && dir.exists(file.path(ref_lib, "EMC2"))
have_card <- nzchar(card) && file.exists(card)
source(test_path("nn_reference_ll.R"), local = TRUE)

reference_inputs <- function(cases, n_trials = 120, n = 40, n_out = 3, seed = 2026) {
  set.seed(seed)
  for (nm in names(cases)) {
    des <- nn_reference_design(cases[[nm]])
    p <- nn_reference_p(des, cases[[nm]])
    cases[[nm]]$dat <- make_data(p, des, n_trials = n_trials)
    P <- matrix(p, n, length(p), byrow = TRUE, dimnames = list(NULL, names(p)))
    P <- P + matrix(rnorm(length(P), sd = .1), n)
    P[seq_len(n_out), 1] <- P[seq_len(n_out), 1] + 20        # out of bounds
    cases[[nm]]$P <- P
  }
  cases
}

reference_subprocess <- function(inputs) {
  inp <- tempfile(fileext = ".rds"); out <- tempfile(fileext = ".rds")
  saveRDS(inputs, inp)
  libs <- paste(c(ref_lib, .libPaths()), collapse = .Platform$path.sep)
  log <- system2(file.path(R.home("bin"), "Rscript"),
                 c("-e", shQuote(sprintf("source(%s)", shQuote(test_path("nn_reference_ll.R")))),
                   shQuote(inp), shQuote(out)),
                 env = c(paste0("R_LIBS=", libs), "R_PROFILE_USER=/dev/null", "NN_REFERENCE_MAIN=1"),
                 stdout = TRUE, stderr = TRUE)
  if (!file.exists(out)) stop("the reference subprocess failed:\n", paste(log, collapse = "\n"))
  readRDS(out)
}

# 1e-12 on each trial's log-likelihood (values of order 1 to 10); on a sum over
# trials 1e-12 of its magnitude (a few hundred: an absolute 1e-12 is about ten
# ulp there, and the two engines add the trials in a different order).
expect_same_likelihood <- function(here, ref, label, tol = 1e-12) {
  expect_identical(dim(here$trialwise), dim(ref$trialwise), label = paste(label, "trialwise dim"))
  expect_true(all(is.finite(here$sum)), label = paste(label, "finite sums"))
  for (f in c("sum", "sum_mt", "trialwise", "trialwise_mt")) {
    scale <- if (startsWith(f, "sum")) pmax(1, abs(ref[[f]])) else 1
    d <- max(abs(here[[f]] - ref[[f]]) / scale)
    expect_lte(d, tol, label = sprintf("%s %s: max |diff| %.2e (%s)", label, f, d,
                                       if (startsWith(f, "sum")) "relative" else "absolute"))
  }
}

test_that("the compiled likelihood of every shipped card equals the reference library's at 1e-12", {
  skip_if_not(have_ref, "EMC2_NLE_REF_LIB not set")
  inputs <- reference_inputs(nn_reference_cases())
  here <- nn_reference_run(inputs)
  ref <- reference_subprocess(inputs)
  expect_true(startsWith(attr(ref, "library"), normalizePath(ref_lib)),
              label = paste("the subprocess loaded", attr(ref, "library")))
  expect_false(identical(attr(ref, "library"), attr(here, "library")))
  for (nm in names(inputs)) expect_same_likelihood(here[[nm]], ref[[nm]], nm)
})

test_that("the compiled likelihood of the card EMC2_NLE_CARD equals the reference library's at 1e-12", {
  skip_if_not(have_ref, "EMC2_NLE_REF_LIB not set")
  skip_if_not(have_card, "EMC2_NLE_CARD not set")
  inputs <- reference_inputs(nn_reference_cases(card)["card"])
  here <- nn_reference_run(inputs)
  ref <- reference_subprocess(inputs)
  expect_same_likelihood(here$card, ref$card, basename(card))
})
