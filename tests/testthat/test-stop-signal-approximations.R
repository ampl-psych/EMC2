approx_ss_design <- function(go = "exgaussian", approximation = NULL,
                             varying = NULL) {
  formula <- if (go == "exgaussian") {
    list(mu ~ lM, sigma ~ 1, tau ~ 1)
  } else {
    list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1)
  }
  formula <- c(formula, list(muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))
  if (!is.null(varying)) {
    # Identical race parameters across stimuli isolate threshold memoization.
    formula[[1]] <- if (go == "exgaussian") mu ~ 1 else v ~ 1
    formula <- c(formula, list(stats::as.formula(paste(varying, "~ S"))))
  }
  constants <- if (is.null(varying)) numeric(0) else
    stats::setNames(c(0, if (varying == "exg_normal_ratio") 3 else .1),
                    c(varying, paste0(varying, "_Sright")))
  suppressMessages(design(
    model = function() withCallingHandlers(
      stop_signal(go = go, stop = "exgaussian", exgaussian_approximation = approximation),
      warning = function(w) {
        if (startsWith(conditionMessage(w), "Using the ex-Gaussian normal approximation"))
          invokeRestart("muffleWarning")
      }),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) d$S == d$lR, formula = formula, constants = constants,
    report_p_vector = FALSE))
}

approx_ss_pars <- function(des) {
  p <- sampled_pars(des, doMap = FALSE)
  values <- c(mu = log(.65), mu_lMTRUE = log(.5/.65), sigma = log(.05),
              tau = log(.08), v = log(1.8), v_lMTRUE = log(.35), B = log(.8),
              A = log(.1), t0 = log(.2), s = 0,
              muS = log(.2), sigmaS = log(.03), tauS = log(.06),
              gf = qnorm(.05), tf = qnorm(.1),
              exg_normal_ratio = 0, exg_normal_ratio_Sright = 3,
              exg_shifted_exp_sigma = 0, exg_shifted_exp_sigma_Sright = .1)
  stopifnot(all(names(p) %in% names(values)))
  p[] <- values[names(p)]
  p
}

approx_ss_data <- function() {
  data.frame(subjects = factor(rep(1, 6)), trials = 1:6,
             S = factor(rep("left", 6), levels = c("left", "right")),
             R = factor(c("left", "right", "left", NA, "left", NA),
                        levels = c("left", "right")),
             rt = c(.52, .72, .62, NA, .60, NA),
             SSD = c(Inf, Inf, .35, .2, .3, .35))
}

approx_ss_ll <- function(dat, des, p, reference = FALSE, trialwise = FALSE) {
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                  compress = FALSE, verbose = FALSE))
  m <- emc[[1]]$model
  if (reference) {
    info <- m()
    info$c_name <- NULL
    m <- function() info
  }
  as.numeric(calc_ll_manager(matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
                             emc[[1]]$data[[1]], m, return_trialwise = trialwise))
}

test_that("SSEXG and SSRDEX honor optional approximations in fitting and censoring", {
  settings <- list(normal = list(normal_ratio = 3),
                   exponential = list(shifted_exponential_sigma = .1))
  for (go in c("exgaussian", "racing_diffusion")) {
    dat <- approx_ss_data()
    baseline <- approx_ss_design(go)
    baseline_ll <- approx_ss_ll(dat, baseline, approx_ss_pars(baseline))
    explicit_zero <- approx_ss_design(go, list(normal_ratio = 0,
                                              shifted_exponential_sigma = 0))
    expect_identical(approx_ss_ll(dat, explicit_zero, approx_ss_pars(explicit_zero)),
                     baseline_ll)
    for (setting in settings) {
      des <- approx_ss_design(go, setting)
      p <- approx_ss_pars(des)
      enabled <- approx_ss_ll(dat, des, p)
      expect_gt(abs(enabled - baseline_ll), 1e-4)
      for (bounds in list(NULL, list(UC = .68), list(LC = .56, UC = .68))) {
        data <- if (is.null(bounds)) dat else
          do.call(make_missing, c(list(dat, rt_resolution = NULL), bounds))
        expect_lt(abs(approx_ss_ll(data, des, p) -
                        approx_ss_ll(data, des, p, reference = TRUE)), 1e-6)
      }
    }
  }
})

test_that("SS integral memoization distinguishes row-varying approximation settings", {
  for (go in c("exgaussian", "racing_diffusion")) {
    for (threshold in c("exg_normal_ratio", "exg_shifted_exp_sigma")) {
      des <- approx_ss_design(go, varying = threshold)
      p <- approx_ss_pars(des)
      dat <- approx_ss_data()[rep(4, 4), ]
      dat$trials <- 1:4
      dat$S[] <- rep(c("left", "right"), 2)
      mixed <- approx_ss_ll(dat, des, p, trialwise = TRUE)
      separate <- vapply(seq_len(nrow(dat)), function(i)
        approx_ss_ll(dat[i, , drop = FALSE], des, p), 0.0)
      expect_equal(mixed, separate, tolerance = 1e-12)
      expect_gt(abs(mixed[1] - mixed[2]), 1e-4)
    }
  }
})
