test_that("SSEXG C++ likelihood matches R and responds to go and stop parameters", {
  des <- ss_test_design("exgaussian")
  expect_stop_cpp_equal_r(des, make_stop_cpp_p(des))
})

test_that("SSEXG and SSRDEX serial and two-thread likelihoods agree", {
  skip_if(testthat::is_parallel(),
          "Run this file serially to avoid nested test-worker/OpenMP parallelism.")
  withr::local_options(emc.ll_backend = "multiprocess", emc.n_threads = 2)
  for (go in c("exgaussian", "racing_diffusion")) {
    des <- ss_test_design(go)
    p <- make_stop_cpp_p(des)
    emc <- make_emc(ss_test_data(), des, type = "single", n_chains = 1,
                    compress = FALSE, verbose = FALSE)
    proposals <- rbind(p, p - .01, p + .01)
    options(emc.ll_backend = "multiprocess")
    serial <- calc_ll_manager(proposals, emc[[1]]$data[[1]], emc[[1]]$model,
                               return_trialwise = TRUE)
    options(emc.ll_backend = "multithreaded")
    threaded <- withCallingHandlers(
      calc_ll_manager(proposals, emc[[1]]$data[[1]], emc[[1]]$model,
                       return_trialwise = TRUE),
      warning = function(w) {
        if (grepl("OpenMP not available", conditionMessage(w), fixed = TRUE))
          skip("OpenMP is not compiled in; run this check on an OpenMP-enabled build.")
      }
    )
    expect_true(all(is.finite(serial)) && all(is.finite(threaded)))
    expect_identical(threaded, serial)
  }
})
