## SBC robustness: a replicate whose prior draw cannot be simulated (make_data
## returns FALSE for out-of-bound parameters) is reported as failed, and a
## try-error replicate does not break assembly.
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

test_that("run_SBC_subject reports an unsimulable prior draw as failed", {
  des <- design(model = SSEXG, factors = list(subjects = 1, S = c("left", "right")),
                Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
                formula = list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1,
                               gf ~ 1, tf ~ 1), report_p_vector = FALSE)
  pn <- names(sampled_pars(des))
  pmean <- c(mu = log(.45), mu_lMTRUE = .3, sigma = log(.05), tau = log(.08), muS = log(.2),
             sigmaS = log(.03), tauS = log(.04), gf = qnorm(.03), tf = qnorm(.1))[pn]
  pri <- prior(des, type = "single", pmean = pmean, psd = setNames(rep(.5, length(pn)), pn))
  # tf far below its lower bound (0.001 on the probability scale) -> make_data() FALSE
  bad <- pmean; bad["tf"] <- qnorm(1e-6)
  prior_alpha <- rbind(bad)
  td <- tempfile("sbc_"); dir.create(td)
  dots <- list(max_tries = 1, compress = FALSE, rt_resolution = 1e-12, cores_per_chain = 1,
               functions = list(SSD = make_ssd(staircase = FALSE, values = c(.2, .3))))
  r <- withCallingHandlers(
    EMC2:::run_SBC_subject(1, des, prior_alpha, trials = 20, pri, dots, td),
    warning = function(w) {
      expect_match(conditionMessage(w), "out of model bounds|Data simulation failed")
      invokeRestart("muffleWarning")
    })
  expect_true(isTRUE(r$failed))
  expect_true(file.exists(file.path(td, "p_vector_rep1.Rdata")))
  unlink(td, recursive = TRUE)
})

test_that("SBC assembly skips failed and try-error replicates", {
  good <- list(rank = c(a = .3, b = .7), med = c(a = 0, b = 1), bias = c(a = .1, b = -.1),
               coverage = c(a = TRUE, b = FALSE))
  reps <- list("1" = good, "2" = list(rank = NULL, failed = TRUE),
               "3" = structure("boom", class = "try-error"), "4" = good)
  SBC <- EMC2:::.sbc_assemble_single(reps)
  expect_equal(attr(SBC, "recovered_reps"), c(1L, 4L))
  expect_equal(nrow(SBC$rank$alpha), 2)
})
