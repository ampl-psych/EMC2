make_minimal_emc <- function(trend, covariate1 = c(-.5, 0, .5, 1, 2), ...) {
  this_design <- design(
    factors = list(subjects = 1, S = 1), Rlevels = 1,
    covariates = c("covariate1", names(list(...))),
    matchfun = function(x) x$S == x$R, trend = trend,
    formula = list(m ~ 1, s ~ 1, t0 ~ 1), model = LNR
  )
  dat <- data.frame(subjects = factor(1), S = factor(1), R = factor(1),
                    rt = 1, covariate1 = covariate1, ...)
  make_emc(dat, this_design, type = "single", n_chains = 1, compress = FALSE)
}

# lin_incr ----------------------------------------------------------------
trend_lin_incr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'lin_incr')))
kernel_pars <- NULL
emc <- make_minimal_emc(trend_lin_incr)
expected_output <- c(-.5, 0, .5, 1, 2)
test_that("lin_incr_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})

# lin_decr ----------------------------------------------------------------
trend_lin_decr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'lin_decr')))
kernel_pars <- NULL
emc <- make_minimal_emc(trend_lin_decr)
expected_output <- -c(-.5, 0, .5, 1, 2)
test_that("lin_decr_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})

# exp_decr ----------------------------------------------------------------
trend_exp_decr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'exp_decr')))
kernel_pars <- c('m.d_ed'=1)
covariate1 <- c(-.5, 0, .5, 1, 2)
emc <- make_minimal_emc(trend_exp_decr, covariate1 = covariate1)
expected_output <- exp(-exp(kernel_pars) * covariate1)
test_that("exp_decr_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})

# exp_incr ----------------------------------------------------------------
trend_exp_incr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'exp_incr')))
kernel_pars <- c('m.d_ei'=1)
covariate1 <- c(-.5, 0, .5, 1, 2)
emc <- make_minimal_emc(trend_exp_incr, covariate1 = covariate1)
expected_output <- 1-exp(-exp(kernel_pars) * covariate1)
test_that("exp_incr_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})

# slin_incr / slin_decr ---------------------------------------------------
# k = +/- min(1, k_sat * c); non-finite covariates (NA, Inf: e.g. SSD on go trials) give 0
trend_slin_incr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'slin_incr')))
kernel_pars <- c('m.k_sat'=log(2))
covariate1 <- c(0, .1, .3, 1, Inf, NA)
emc <- make_minimal_emc(trend_slin_incr, covariate1 = covariate1)
expected_output <- ifelse(is.finite(covariate1), pmin(1, exp(kernel_pars) * covariate1), 0)
test_that("slin_incr_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
})
trend_slin_decr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'slin_decr')))
emc <- make_minimal_emc(trend_slin_decr, covariate1 = covariate1)
test_that("slin_decr_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(-expected_output))
})

# pow_decr ----------------------------------------------------------------
trend_pow_decr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'pow_decr')))
kernel_pars <- c('m.d_pd'=1)
covariate1 <- c(-.5, 0, .5, 1, 2)
emc <- make_minimal_emc(trend_pow_decr, covariate1 = covariate1)
expected_output = (1 + covariate1)^(-exp(kernel_pars))
test_that("pow_decr_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})

# pow_incr ----------------------------------------------------------------
trend_pow_incr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'pow_incr')))
kernel_pars <- c('m.d_pi'=1)
covariate1 <- c(-.5, 0, .5, 1, 2)
emc <- make_minimal_emc(trend_pow_incr, covariate1 = covariate1)
expected_output <- 1-(1 + covariate1)^(-exp(kernel_pars))
test_that("pow_incr_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})

# poly2: Quadratic polynomial: k = d1 * c + d2 * c^2 -----------
trend_poly2 <- make_trend(make_base('m', 'add', make_kernel('covariate1', 'poly2')))
kernel_pars <- c('m.d1'=1, 'm.d2'=1)
covariate1 <- c(-.5, 0, .5, 1, 2)
emc <- make_minimal_emc(trend_poly2, covariate1 = covariate1)
expected_output <- kernel_pars[[1]]*covariate1+kernel_pars[[2]]*covariate1^2
test_that("poly2_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})

# poly3: Cubic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3 ----------
trend_poly3 <- make_trend(make_base('m', 'add', make_kernel('covariate1', 'poly3')))
kernel_pars <- c('m.d1'=1, 'm.d2'=1, 'm.d3'=1)
covariate1 <- c(-.5, 0, .5, 1, 2)
emc <- make_minimal_emc(trend_poly3, covariate1 = covariate1)
expected_output <- kernel_pars[[1]]*covariate1+kernel_pars[[2]]*covariate1^2+kernel_pars[[3]]*covariate1^3
test_that("poly3_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})

# poly4: Quartic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3 + d4 * c^4 ---------
trend_poly4 <- make_trend(make_base('m', 'add', make_kernel('covariate1', 'poly4')))
kernel_pars <- c('m.d1'=1, 'm.d2'=1, 'm.d3'=1, 'm.d4'=.1)
covariate1 <- c(-.5, 0, .5, 1, 2)
emc <- make_minimal_emc(trend_poly4, covariate1 = covariate1)
expected_output <- kernel_pars[[1]]*covariate1+kernel_pars[[2]]*covariate1^2+kernel_pars[[3]]*covariate1^3+kernel_pars[[4]]*covariate1^4
test_that("poly4_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})


# LEARNING RULES ----------------------------------------------------------
# include NA check

# delta --------
trend_delta <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'delta')))
kernel_pars <- c('m.q0'=0.5, 'm.alpha'=qnorm(0.2))
covariate1 <- c(NA, 1, NA, 1, NA)
emc <- make_minimal_emc(trend_delta, covariate1 = covariate1)
expected_output <- matrix(c(0.5, 0.5, 0.6, 0.6, 0.68)) # manually computed for this specific covariate + parameter vector
test_that("delta_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})


# delta 2 LR --------
trend_delta2lr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'delta2lr')))
kernel_pars <- c('m.q0'=0.5, 'm.alphaPos'=qnorm(0.5), 'm.alphaNeg' = qnorm(0.25))
covariate1 <- c(NA, 1, NA, 0, NA)
emc <- make_minimal_emc(trend_delta2lr, covariate1 = covariate1)
expected_output <- matrix(c(0.5, 0.5, 0.75, 0.75, 0.5625)) # manually computed for this specific covariate + parameter vector
test_that("delta2lr_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})


## test resetting
trend_delta <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'delta', kernel_args = list(q_reset_column='do_reset'))))
kernel_pars <- c('m.q0'=0, 'm.alpha'=qnorm(0.2))
covariate1 <- c(1, 1, 1, 1, 1)
emc <- make_minimal_emc(trend_delta, covariate1 = covariate1, do_reset=c(F,F,T,F,F))
expected_output <- matrix(c(0, 0.2, 0, 0.2, 0.36)) # manually computed for this specific covariate + parameter vector
test_that("delta reset_Rcpp", {
  expect_equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output),
               tolerance = 1e-8)
})



# Rescorla-Wagner uses the summed expectation of active cues.
trend_rw <- make_trend(make_base("m", "add",
  make_kernel(paste0("covariate", 1:4), "rescorlawagner",
    kernel_args = list(q_reset_column = "do_reset"))))
covariates <- matrix(c(NA, NA, 1, 1,
                       1, 1, NA, NA,
                       NA, 0, 0, NA,
                       1, NA, NA, 1,
                       1, 1, 1, 1,
                       0, 0, 0, 0), ncol = 4, byrow = TRUE)
emc <- make_minimal_emc(
  trend_rw, covariate1 = covariates[, 1], covariate2 = covariates[, 2],
  covariate3 = covariates[, 3], covariate4 = covariates[, 4],
  do_reset = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)
)
kernel_pars <- c(m.q0 = 0, m.alpha = qnorm(.2))
expected_output <- matrix(c(0, 0, 0, 0,
                            0, 0, .2, .2,
                            .2, .2, .2, .2,
                            .2, .12, .12, .2,
                            0, 0, 0, 0,
                            .2, .2, .2, .2), ncol = 4, byrow = TRUE)
test_that("rescorlawagner_Rcpp", {
  expect_equal(unname(apply_kernel(kernel_pars, emc)), expected_output,
               tolerance = 1e-8)
})

# ----------------------------------------------------------------------------~
# BETA-BINOMIAL, DBM, TPM LEARNING RULES                                  -----
# ----------------------------------------------------------------------------~

# NB using helper functions defined in tests/testthat/helper-kernels.R that are
# independent of the EMC2 trends framework

# shared covariate vector (true mean roughly 0.3; includes NA)
covariate1 <- c(
  0, 1, 0, 0, 0, 0, 1, NA, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, NA, 0, 1, 0, 0
)

make_minimal_emc_wrap <- function(kernel_name, output_code = 1L) {
  make_minimal_emc(
    trend = make_trend(
      make_base(
        target_parameter = "m",
        type = "lin",
        kernel = make_kernel(
          cov_names = "covariate1",
          type = kernel_name
        ),
        kernel_output = output_code
      )
    ),
    covariate1 = covariate1
  )
}

output_codes <- 1L
output_types <- c("mean")

# Beta-binomial (basic)
emc <- setNames(lapply(
  X = output_codes,
  FUN = function(x) {
    make_minimal_emc_wrap("beta_binomial", output_code = x)
  }
), output_types)
snapshot_list <- vector(mode = "list", length = 3L)
names(snapshot_list) <- c("uniform", "uniform_reset", "informative")

# uniform prior: Beta shape parameters = 1
kernel_pars <- c("m.a0" = log(1), "m.b0" = log(1))
true_out <- sapply(
  X = output_types,
  FUN = function(output_type) {
    beta_binomial(
      x = covariate1, a0 = exp(kernel_pars[1]), b0 = exp(kernel_pars[2]),
      output_type = output_type
    )
  }
)
emc_out <- sapply(
  X = emc,
  FUN = function(x) {as.numeric(apply_kernel(kernel_pars, x))}
)
snapshot_list[["uniform"]] <- list(
  manual = true_out,
  emc = emc_out
)

# test belief resetting
reset_column <- rep(FALSE, times = length(covariate1))
reset_idx <- floor(length(covariate1)/2)
reset_column[reset_idx] <- TRUE

emc_reset <- make_minimal_emc(
  trend = make_trend(
    make_base(
      target_parameter = "m",
      type = "lin",
      kernel = make_kernel(
        cov_names = "covariate1",
        type = "beta_binomial",
        kernel_args = list(belief_reset_column = "do_reset")
      ),
      kernel_output = 1L
    )
  ),
  covariate1 = covariate1,
  do_reset = reset_column
)

true_out <- beta_binomial(
  x = covariate1, a0 = exp(kernel_pars[1]), b0 = exp(kernel_pars[2]),
  belief_reset = reset_column
)
emc_out <- as.numeric(apply_kernel(kernel_pars, emc_reset))
snapshot_list[["uniform_reset"]] <- list(
  manual = true_out,
  emc = emc_out
)

# informative prior: shape1 = mean * scale; shape2 = (1-mean) * scale
kernel_pars <- c(
  "m.a0" = log(mean(covariate1, na.rm = TRUE) * 10),
  "m.b0" = log((1 - mean(covariate1, na.rm = TRUE)) * 10)
)
true_out <- sapply(
  X = output_types,
  FUN = function(output_type) {
    beta_binomial(
      x = covariate1, a0 = exp(kernel_pars[1]), b0 = exp(kernel_pars[2]),
      output_type = output_type
    )
  }
)
emc_out <- sapply(
  X = emc,
  FUN = function(x) {as.numeric(apply_kernel(kernel_pars, x))}
)
snapshot_list[["informative"]] <- list(
  manual = true_out,
  emc = emc_out
)

test_that("beta_binomial_Rcpp", {
  # Compare the compiled kernel with the independent R calculation.
  for (case in snapshot_list) expect_equal(case$emc, case$manual)
  reset <- snapshot_list$uniform_reset$emc
  expect_identical(reset[1], reset[reset_idx])
  expect_snapshot(snapshot_list)
})


# Beta-binomial (exponential decay)
emc <- lapply(
  X = output_codes,
  FUN = function(x) {
    make_minimal_emc_wrap("beta_binomial_decay", output_code = x)
  }
)
names(emc) <- output_types

# uniform prior with decay = 4: Updates are weighted by exp(-1/4); implies half-life of 4*log(2)=2.8 trials
kernel_pars <- c("m.a0" = log(1), "m.b0" = log(1), "m.decay" = log(4))
true_out <- sapply(
  X = output_types,
  FUN = function(output_type) {
    beta_binomial(
      x = covariate1, a0 = exp(kernel_pars[1]), b0 = exp(kernel_pars[2]),
      decay = exp(kernel_pars[3]),
      output_type = output_type
    )
  }
)
emc_out <- sapply(
  X = emc,
  FUN = function(x) {as.numeric(apply_kernel(kernel_pars, x))}
)
test_that("beta_binomial_decay_Rcpp", {
  expect_equal(emc_out, true_out)
  expect_snapshot(list(manual = true_out, emc = emc_out))
})


# Beta-binomial (sliding window)
emc <- lapply(
  X = output_codes,
  FUN = function(x) {
    make_minimal_emc_wrap("beta_binomial_window", output_code = x)
  }
)
names(emc) <- output_types
snapshot_list <- vector(mode = "list", length = 2L)
names(snapshot_list) <- c("uniform", "uniform_reset")

# uniform prior with memory window of 6 trials
kernel_pars <- c("m.a0" = log(1), "m.b0" = log(1), "m.window" = log(6))
true_out <- sapply(
  X = output_types,
  FUN = function(output_type) {
    beta_binomial(
      x = covariate1, a0 = exp(kernel_pars[1]), b0 = exp(kernel_pars[2]),
      window = exp(kernel_pars[3]),
      output_type = output_type
    )
  }
)
emc_out <- sapply(
  X = emc,
  FUN = function(x) {as.numeric(apply_kernel(kernel_pars, x))}
)
snapshot_list[["uniform"]] <- list(
  manual = true_out,
  emc = emc_out
)

# test belief resetting
emc_reset <- make_minimal_emc(
  trend = make_trend(
    make_base(
      target_parameter = "m",
      type = "lin",
      kernel = make_kernel(
        cov_names = "covariate1",
        type = "beta_binomial_window",
        kernel_args = list(belief_reset_column = "do_reset")
      ),
      kernel_output = 1L
    )
  ),
  covariate1 = covariate1,
  do_reset = reset_column
)

true_out <- beta_binomial(
  x = covariate1, a0 = exp(kernel_pars[1]), b0 = exp(kernel_pars[2]),
  window = exp(kernel_pars[3]),
  belief_reset = reset_column
)
emc_out <- as.numeric(apply_kernel(kernel_pars, emc_reset))
snapshot_list[["uniform_reset"]] <- list(
  manual = true_out,
  emc = emc_out
)

test_that("beta_binomial_window_Rcpp", {
  for (case in snapshot_list) expect_equal(case$emc, case$manual)
  reset <- snapshot_list$uniform_reset$emc
  expect_identical(reset[1], reset[reset_idx])
  expect_snapshot(snapshot_list)
})


# Dynamic Belief Model
emc <- lapply(
  X = output_codes,
  FUN = function(x) {
    make_minimal_emc_wrap("dbm", output_code = x)
  }
)
names(emc) <- output_types
snapshot_list <- vector(mode = "list", length = 4L)
names(snapshot_list) <- c("cp_0", "cp_1", "cp_real", "cp_real_reset")

# change point probability of zero, output should be equivalent to Beta binomial
kernel_pars <- c(
  "m.cp" =  qnorm(1e-12),
  "m.mu0" = qnorm(mean(covariate1, na.rm = TRUE)),
  "m.s0" = log(10)
)
true_out <- sapply(
  X = output_types,
  FUN = function(output_type) {
    dbm(
      x = covariate1, cp = pnorm(kernel_pars[1]), mu0 = pnorm(kernel_pars[2]),
      s0 = exp(kernel_pars[3])
    )
  }
)
emc_out <- sapply(
  X = emc,
  FUN = function(x) {as.numeric(apply_kernel(kernel_pars, x))}
)

emc_out_bb <- sapply(
  X = output_codes,
  FUN = function(x) {
    as.numeric(
      apply_kernel(
        kernel_pars = c(
          "m.a0" = log(mean(covariate1, na.rm = TRUE) * 10),
          "m.b0" = log((1 - mean(covariate1, na.rm = TRUE)) * 10)
        ),
        emc = make_minimal_emc_wrap("beta_binomial", output_code = x)
      )
    )
  }
)
colnames(emc_out_bb) <- output_types

snapshot_list[["cp_0"]] <- list(
  manual = true_out,
  emc = emc_out,
  emc_bb = emc_out_bb
)

# change point probability of one, output should be constant across trials
kernel_pars <- c(
  "m.cp" =  qnorm(1 - 1e-12),
  "m.mu0" = qnorm(mean(covariate1, na.rm = TRUE)),
  "m.s0" = log(10)
)
true_out <- sapply(
  X = output_types,
  FUN = function(output_type) {
    dbm(
      x = covariate1, cp = pnorm(kernel_pars[1]), mu0 = pnorm(kernel_pars[2]),
      s0 = exp(kernel_pars[3])
    )
  }
)
emc_out <- sapply(
  X = emc,
  FUN = function(x) {as.numeric(apply_kernel(kernel_pars, x))}
)

fixed_beta_out <- matrix(nrow = length(covariate1), ncol = length(output_types))
colnames(fixed_beta_out) <- output_types
fixed_beta_out[ , "mean"] <- rep(
  x = beta_mean(
    a = unname(pnorm(kernel_pars[2])) * unname(exp(kernel_pars[3])),
    b = (1 - unname(pnorm(kernel_pars[2]))) * unname(exp(kernel_pars[3]))
  ),
  length = length(covariate1)
)

snapshot_list[["cp_1"]] <- list(
  manual = true_out,
  emc = emc_out,
  fixed_beta = fixed_beta_out
)


# change point probability of 0.2
kernel_pars <- c(
  "m.cp" =  qnorm(0.2),
  "m.mu0" = qnorm(mean(covariate1, na.rm = TRUE)),
  "m.s0" = log(10)
)
true_out <- sapply(
  X = output_types,
  FUN = function(output_type) {
    dbm(
      x = covariate1, cp = pnorm(kernel_pars[1]), mu0 = pnorm(kernel_pars[2]),
      s0 = exp(kernel_pars[3])
    )
  }
)
emc_out <- sapply(
  X = emc,
  FUN = function(x) {as.numeric(apply_kernel(kernel_pars, x))}
)

snapshot_list[["cp_real"]] <- list(
  manual = true_out,
  emc = emc_out
)

# test belief resetting
emc_reset <- make_minimal_emc(
  trend = make_trend(
    make_base(
      target_parameter = "m",
      type = "lin",
      kernel = make_kernel(
        cov_names = "covariate1",
        type = "dbm",
        kernel_args = list(belief_reset_column = "do_reset")
      ),
      kernel_output = 1L
    )
  ),
  covariate1 = covariate1,
  do_reset = reset_column
)

true_out <- dbm(
  x = covariate1, cp = pnorm(kernel_pars[1]), mu0 = pnorm(kernel_pars[2]),
  s0 = exp(kernel_pars[3]),
  belief_reset = reset_column
)
emc_out <- as.numeric(apply_kernel(kernel_pars, emc_reset))
snapshot_list[["cp_real_reset"]] <- list(
  manual = true_out,
  emc = emc_out
)

test_that("dbm_Rcpp", {
  for (case in snapshot_list) expect_equal(case$emc, case$manual)
  # At the two limits, DBM reduces to cumulative learning or a fixed prior.
  expect_equal(snapshot_list$cp_0$emc, snapshot_list$cp_0$emc_bb, tolerance = 1e-6)
  expect_equal(snapshot_list$cp_1$emc, snapshot_list$cp_1$fixed_beta, tolerance = 1e-5)
  reset <- snapshot_list$cp_real_reset$emc
  expect_identical(reset[1], reset[reset_idx])
  expect_snapshot(snapshot_list)
})


# minimal Transition Probability Model (TPM)

# NB using Matlab toolbox (copyright 2016 Florent Meyniel & Maxime Maheu;
# https://github.com/florentmeyniel/MinimalTransitionProbsModel) as ground truth.
# This toolbox doesn't accomodate missing obs, hence we change the covariate1
# vector as follows:
covariate1 <- c(
  1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 1
) - 1

# Matlab script:
# clear; close('all');
# try cd('MarkovReview'); catch, end;
# addpath('Tools');
# addpath('IdealObserversCode');
# s = [1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 1];
#
# % Set parameters
# in.s            = s;                % sequence
# in.learned      = 'transition';     % estimate transition
# in.jump         = 1;                % estimate with jumps
# in.mode         = 'HMM';            % use the HMM (not sampling) algorithm
# in.opt.pJ       = 0.2;              % a priori probability that a jump occur at each outcome
# n               = 100;              % resolution of the univariate probability grid
# in.opt.pgrid    = linspace(0,1,n);  % estimation probability grid
# in.opt.Alpha0   = ones(n)/(n^2);    % uniform prior on transition probabilities
# in.verbose      = 1;                % to check that no default values are used.
#
# % Compute the observer
# out = IdealObserver(in);
# 1 - out.p1_mean'

true_out <- c(
  0.5000, 0.5000, 0.6088, 0.4374, 0.3561, 0.3108, 0.4541, 0.4895, 0.3863,
  0.3364, 0.3051, 0.4382, 0.4741, 0.3767, 0.3971, 0.5030, 0.4020, 0.3810,
  0.5145, 0.4139, 0.3721, 0.5201, 0.3392, 0.5822, 0.4502, 0.3565, 0.5490,
  0.3248, 0.5993, 0.4633
)

emc <- make_minimal_emc(
  trend = make_trend(
    make_base(
      target_parameter = "m",
      type = "lin",
      kernel = make_kernel(
        cov_names = "covariate1",
        type = "tpm"
      ),
      kernel_output = 1L
    )
  ),
  covariate1 = covariate1
)

emc <- make_minimal_emc_wrap("tpm", output_code = 1L)
kernel_pars <- c("m.cp" =  qnorm(0.2), "m.a0" = log(1), "m.b0" = log(1))
emc_out <- as.numeric(apply_kernel(kernel_pars, emc))

# NB after thorough debugging, it turns out that the output from the Matlab toolbox
# is offset by one index relative to our output. Specifically, in Matlab,
# out.p1_mean(t) is the model's forecast for trial t+1, computed from data through
# trial t. This is in contrast to our implementation of Beta-binomial / DBM-like
# kernels: our output for trial t is the forecast for trial t, computed before
# seeing trial t (i.e., a true one-step-ahead prediction).
# After taking this into account, our output is correct:

test_that("tpm_Rcpp", {
  expect_equal(emc_out[-1], true_out[-length(true_out)], tolerance = 1e-3)
  expect_snapshot(list(Matlab = true_out, emc = emc_out))
})
