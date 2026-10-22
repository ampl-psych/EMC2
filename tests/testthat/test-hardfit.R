# hardfit(): the record it writes, the early stop and the continuation;
# diagnose(): the outcome classes on constructed draws. The one-block primitive
# (run_emc() with an iteration-only criterion) is tested in test-sample-kernel.R.

test_that("hardfit() runs the stages, writes the record and continues", {
  withr::local_seed(123)
  withr::local_options(emc.adapt_converge = FALSE)
  emc <- tiny_lnr()
  expect_message(h <- suppressWarnings(hardfit(emc, iter = 40, step_size = 20, max_tries = 3, burn = list(max_tries = 2),
                                               min_ess = NULL, max_subj_spread = NULL, cull = FALSE, abort_gd = NULL,
                                               cores_for_chains = 1, particle_factor = 10, verbose = TRUE)),
                 "Running sample stage")
  r <- h[[1]]$stop_record
  expect_true(all(c("criteria", "met", "tries", "rule", "stats", "culls", "discards", "aborted", "restarts", "trace", "previous", "time") %in% names(r)))
  expect_gte(unname(chain_n(h)[1, "sample"]), 40)
  expect_equal(nrow(r$trace), r$tries); expect_equal(r$trace$draws[nrow(r$trace)], unname(chain_n(h)[1, "sample"]))
  expect_equal(r$criteria$max_gd, 1.1); expect_null(r$criteria$min_ess); expect_null(r$criteria$cull); expect_null(r$criteria$abort_gd)
  expect_null(r$previous); expect_equal(r$restarts, 0L); expect_equal(r$culls, 0L)
  expect_true(is.finite(r$stats$ll_gap)); expect_true(is.list(r$flat)); expect_s3_class(r$prior, "data.frame")
  expect_identical(get_last_stage(h), "sample")
  expect_s3_class(h, "emc")
  d <- diagnose(h)
  expect_true(d$class %in% c("met", "ess_held", "ridge_only", "region_mode"))
  # a continuation adds `iter` draws, keeps the kernel and nests the earlier record
  n1 <- unname(chain_n(h)[1, "sample"])
  h2 <- suppressWarnings(hardfit(h, iter = 20, step_size = 20, max_tries = 2, min_ess = NULL, max_subj_spread = NULL,
                                 cull = FALSE, abort_gd = NULL, cores_for_chains = 1, particle_factor = 10, verbose = FALSE))
  expect_gte(unname(chain_n(h2)[1, "sample"]), n1 + 20)
  expect_identical(h2[[1]]$sample_kernel, h[[1]]$sample_kernel)
  expect_identical(h2[[1]]$stop_record$previous$trace, r$trace)
  # the restart after a cull: one draw kept, reburn draws under the same kernel, one draw kept
  rs <- EMC2:::cull_restart(h2, 20, function(e, n) run_emc(e, "sample", stop_criteria = list(iter = n), step_size = n, max_tries = 1,
                                                        cores_for_chains = 1, verbose = FALSE, particle_factor = 10, trim = FALSE))
  expect_equal(unname(chain_n(rs)[, "sample"]), rep(1, 2)); expect_identical(rs[[1]]$sample_kernel, h2[[1]]$sample_kernel)
  expect_false(identical(rs[[1]]$samples$theta_mu[, rs[[1]]$samples$idx], h2[[1]]$samples$theta_mu[, h2[[1]]$samples$idx]))
  # validation
  expect_error(hardfit(emc, abort_gd = 1), "above 1")
  expect_error(hardfit(emc, min_ess = -1), "min_ess")
  expect_error(hardfit(emc, gd_quantile = 0), "gd_quantile")
  expect_error(hardfit(emc, cull = list(bogus = 1)), "unknown")
  expect_error(hardfit(emc, stage = "warm"), "stage must be")
})

test_that("the early stop ends the sample stage at its check and is recorded", {
  withr::local_seed(123)
  withr::local_options(emc.adapt_converge = FALSE)
  emc <- tiny_lnr()
  w <- capture_warnings(h <- hardfit(emc, iter = 1000, step_size = 20, max_tries = 4, burn = list(max_tries = 2), max_gd = 1.0001,
                                     min_ess = NULL, max_subj_spread = NULL, cull = FALSE, abort_gd = 1.0001, abort_try = 2,
                                     cores_for_chains = 1, particle_factor = 10, verbose = FALSE))
  expect_true(any(grepl("stopped early at check 2", w))); expect_false(any(grepl("see diagnose", w)))
  r <- h[[1]]$stop_record
  expect_false(r$met); expect_equal(r$aborted$try, 2); expect_gt(r$aborted$rule, 1.0001); expect_equal(r$tries, 2)
  expect_equal(unname(chain_n(h)[1, "sample"]), 40)
  expect_equal(diagnose(h)$class, "aborted")
  # without the early stop the same settings run to max_tries, unmet, every draw kept
  w <- capture_warnings(m <- hardfit(emc, iter = 1000, step_size = 20, max_tries = 3, burn = list(max_tries = 2), max_gd = 1.0001,
                                     min_ess = NULL, max_subj_spread = NULL, cull = FALSE, abort_gd = NULL,
                                     cores_for_chains = 1, particle_factor = 10, verbose = FALSE))
  expect_true(any(grepl("see diagnose", w)))
  expect_false(m[[1]]$stop_record$met); expect_null(m[[1]]$stop_record$aborted)
  expect_equal(unname(chain_n(m)[1, "sample"]), 60); expect_length(m[[1]]$stop_record$discards, 0)
})

test_that("diagnose() reads the outcome classes from the draws and the record", {
  n <- unname(chain_n(samples_LNR)[1, "sample"]); it <- which(samples_LNR[[1]]$samples$stage == "sample")
  base <- function(...) diagnose(samples_LNR, max_gd = 10, ...)
  expect_equal(base(min_ess = NULL)$class, "met")
  expect_equal(base(min_ess = 1e6)$class, "ess_held")
  expect_false(base(min_ess = 1e6)$falling)
  # a chain shifted throughout on every group mean: the rule fails along identified directions
  shift <- samples_LNR; shift[[1]]$samples$theta_mu[, it] <- shift[[1]]$samples$theta_mu[, it] + 20
  ds <- diagnose(shift, min_ess = NULL)
  expect_false(ds$rule$met); expect_gt(ds$stats$rhat_ident, 1.1); expect_equal(ds$class, "region_mode")
  # chains apart along a ridge (a + b pinned, a - b free): the identified direction agrees
  set.seed(5); ridge <- samples_LNR; p <- nrow(ridge[[1]]$samples$theta_mu)
  for(i in seq_along(ridge)){
    u <- rnorm(n, i * 3, 1)
    m <- ridge[[i]]$samples$theta_mu
    m[1, it] <- u + rnorm(n, 0, .05); m[2, it] <- -u + rnorm(n, 0, .05)
    if(p > 2) for(j in 3:p) m[j, it] <- rnorm(n)
    ridge[[i]]$samples$theta_mu <- m
  }
  dr <- diagnose(ridge, min_ess = NULL)
  expect_false(dr$rule$met); expect_lt(dr$stats$rhat_ident, 1.1); expect_gt(dr$stats$rhat_b, 1.5)
  expect_equal(dr$class, "ridge_only")
  # the record decides aborted and met, and the trace decides falling
  ab <- shift; ab[[1]]$stop_record <- list(aborted = list(try = 20, rule = 2, abort_gd = 1.5), met = FALSE)
  expect_equal(diagnose(ab)$class, "aborted")
  fell <- shift; fell[[1]]$stop_record <- list(met = FALSE, trace = data.frame(try = 1:4, rule = c(2, 1.8, 1.6, 1.4)))
  df <- diagnose(fell, min_ess = NULL)
  expect_equal(df$class, "region_mode"); expect_true(df$falling); expect_match(df$next_step, "continuation")
  expect_output(print(df), "still falling")
})

test_that("prior tension flags a group mean far from its prior mean, and is a report, not a criterion", {
  d <- diagnose(samples_LNR, max_gd = 10, min_ess = NULL)
  expect_s3_class(d$prior, "data.frame")
  expect_equal(d$prior$parameter, rownames(samples_LNR[[1]]$samples$theta_mu))
  expect_equal(d$prior$prior_sd, rep(1, 4)); expect_false(any(d$prior$flagged))
  expect_false(grepl("Prior tension", paste(capture.output(print(d)), collapse = "\n")))
  # the prior mean of m moved 10 SDs away: m flagged with the sign of the shift, the class unchanged
  far <- samples_LNR; far[[1]]$prior$theta_mu_mean["m"] <- far[[1]]$prior$theta_mu_mean["m"] + 10
  dfar <- diagnose(far, max_gd = 10, min_ess = NULL)
  expect_equal(dfar$class, d$class)
  expect_equal(dfar$prior$parameter[dfar$prior$flagged], "m"); expect_lt(dfar$prior$z[dfar$prior$parameter == "m"], -9)
  expect_match(dfar$next_step, "posterior mean of m is more than 3 prior SDs")
  expect_output(print(dfar), "Prior tension .*: m -")
  expect_false(any(diagnose(far, max_gd = 10, min_ess = NULL, prior_se = 20)$prior$flagged))
  # a vector prior variance, and no prior mean at all
  vec <- far; vec[[1]]$prior$theta_mu_var <- diag(far[[1]]$prior$theta_mu_var)
  expect_equal(diagnose(vec, max_gd = 10, min_ess = NULL)$prior, dfar$prior)
  none <- samples_LNR; none[[1]]$prior$theta_mu_mean <- NULL
  expect_null(diagnose(none, max_gd = 10, min_ess = NULL)$prior)
  # no group level (type single): one row per parameter and subject, read against the parameter's prior
  single <- far
  for(i in seq_along(single)){ single[[i]]$samples$theta_mu <- NULL; single[[i]]$samples$theta_var <- NULL }
  single[[1]]$type <- "single"
  ps <- diagnose(single, max_gd = 10, min_ess = NULL)$prior
  expect_equal(nrow(ps), 4 * dim(samples_LNR[[1]]$samples$alpha)[2])
  expect_true(all(grepl("^m:", ps$parameter[ps$flagged])))
})
