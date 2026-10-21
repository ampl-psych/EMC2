RNGkind("L'Ecuyer-CMRG")
set.seed(123)

test_that("information criteria and model weights agree with fixed likelihoods", {
  y <- c(-1.2, -.7, -.3, .4, .6)
  dat <- data.frame(subjects = factor(c("a", "a", "b", "b", "b")),
                    R = factor("response"), rt = .2 + exp(y))
  des <- design(data = dat, model = LNR, formula = list(m ~ 1),
                 constants = c(s = log(.7), t0 = log(.2)), report_p_vector = FALSE)
  draws <- rbind(a = c(-1.5, -1.1, -.8, -.6), b = c(-.5, .1, .2, .4))
  models <- lapply(c(A = 0, B = .35), function(shift) {
    alpha <- draws + shift
    # One accumulator: lognormal density, including its Jacobian.
    ll <- vapply(seq_along(y), function(i)
      dnorm(y[i], alpha[as.character(dat$subjects[i]), ], .7, log = TRUE) - y[i],
      numeric(ncol(alpha)))
    subject_ll <- rowsum(t(ll), dat$subjects)
    mean_ll <- rowsum(dnorm(y, rowMeans(alpha)[dat$subjects], .7, log = TRUE) - y,
                       dat$subjects)[, 1]
    emc <- make_emc(dat, des, type = "single", n_chains = 2,
                    compress = FALSE, rt_resolution = 1e-12)
    for (chain in 1:2) {
      idx <- (chain - 1) * 2 + 1:2
      emc[[chain]]$samples <- list(
        alpha = array(alpha[, idx], c(1, 2, 2), dimnames = list("m", rownames(alpha), NULL)),
        subj_ll = subject_ll[, idx], stage = rep("sample", 2), idx = 2)
    }
    list(emc = emc, ll = ll, deviance = -2 * subject_ll, mean_deviance = -2 * mean_ll)
  })

  for (best in c(FALSE, TRUE)) {
    expected <- t(vapply(models, function(model) {
      meanD <- mean(colSums(model$deviance))
      Dmean <- sum(model$mean_deviance)
      minima <- apply(model$deviance, 1, min)
      if (best) minima <- pmin(minima, model$mean_deviance)
      minD <- sum(minima)
      pD <- meanD - if (best) minD else Dmean
      waic <- -2 * sum(log(colMeans(exp(model$ll))) - apply(model$ll, 2, var))
      c(WAIC = waic, DIC = meanD + pD, BPIC = meanD + 2*pD,
        EffectiveN = pD, meanD = meanD, Dmean = Dmean, minD = minD)
    }, numeric(7)))
    actual <- compare(lapply(models, `[[`, "emc"), BayesFactor = FALSE,
                      use_best_fit = best, print_summary = FALSE)
    expect_equal(as.matrix(actual[, colnames(expected)]), expected, tolerance = 1e-8)
    for (metric in c("WAIC", "DIC", "BPIC")) {
      weights <- exp(-expected[, metric] / 2)
      expect_equal(actual[[paste0("w", metric)]], unname(weights / sum(weights)), tolerance = 1e-8)
    }
  }
})

test_that("savage-dickey", {
  expect_snapshot(
    round(hypothesis(samples_LNR, parameter = "m", do_plot = F, H0 = -1), 2))
  expect_snapshot(
    round(hypothesis(samples_LNR, fun = function(d) d["m"] - d["m_lMd"],
                  H0 = -0.5, do_plot = F), 2))
})
