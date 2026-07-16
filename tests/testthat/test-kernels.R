# Convenience functions ---------------------------------------------------------
make_tiny_design <- function(trend, ...) {
  cov_names <- trend$kernels[[1]]$cov_names
  for(nm in names(list(...))) cov_names <- c(cov_names, nm)

  return(design(factors = list(subjects = 1, S = 1), Rlevels = 1,
                covariates = cov_names, matchfun = function(x) x$S==x$R, trend = trend,
                formula = list(m ~ 1, s ~ 1, t0 ~ 1),
                model = LNR, report_p_vector = FALSE))
}

make_minimal_emc <- function(trend, n_trials=5, covariate1=1:5, ...) {
  ## This *ONLY* creates an EMC object with the covariates. Kernel pars are ignored
  this_design <- make_tiny_design(trend, ...)
  p_vector <- sampled_pars(this_design, doMap = FALSE)  # no kernel pars

  opts <- list(...)

  covariates <- data.frame(covariate1=covariate1)
  for(nm in names(opts)) {
    covariates[,nm] <- opts[[nm]]
  }

  dat <- make_data(p_vector, this_design, n_trials = n_trials, covariates = covariates)


  emc <- make_emc(dat, this_design, type='single')
  return(emc)
}


# lin_incr ----------------------------------------------------------------
trend_lin_incr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'lin_incr')))
#trend_lin_incr <- make_trend(par_names = "m",
#                             cov_names = 'covariate1',
#                             kernels = 'lin_incr')
kernel_pars <- NULL
emc <- make_minimal_emc(trend_lin_incr)
expected_output <- 1:5
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("lin_incr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})

# lin_decr ----------------------------------------------------------------
trend_lin_decr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'lin_decr')))
kernel_pars <- NULL
emc <- make_minimal_emc(trend_lin_decr)
expected_output <- -(1:5)
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("lin_decr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})

# exp_decr ----------------------------------------------------------------
trend_exp_decr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'exp_decr')))
kernel_pars <- c('m.d_ed'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_exp_decr, covariate1 = covariate1)
expected_output <- exp(-exp(kernel_pars) * covariate1)
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("exp_decr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})

# exp_incr ----------------------------------------------------------------
trend_exp_incr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'exp_incr')))
kernel_pars <- c('m.d_ei'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_exp_incr, covariate1 = covariate1)
expected_output <- 1-exp(-exp(kernel_pars) * covariate1)
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("exp_incr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})

# pow_decr ----------------------------------------------------------------
trend_pow_decr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'pow_decr')))
kernel_pars <- c('m.d_pd'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_pow_decr, covariate1 = covariate1)
expected_output = (1 + covariate1)^(-exp(kernel_pars))
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("pow_decr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})

# pow_incr ----------------------------------------------------------------
trend_pow_incr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'pow_incr')))
kernel_pars <- c('m.d_pi'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_pow_incr, covariate1 = covariate1)
expected_output <- 1-(1 + covariate1)^(-exp(kernel_pars))
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("pow_incr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})

# poly2: Quadratic polynomial: k = d1 * c + d2 * c^2 -----------
trend_poly2 <- make_trend(make_base('m', 'add', make_kernel('covariate1', 'poly2')))
kernel_pars <- c('m.d1'=1, 'm.d2'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_poly2, covariate1 = covariate1)
expected_output <- kernel_pars[[1]]*covariate1+kernel_pars[[2]]*covariate1^2
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("poly2_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})

# poly3: Cubic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3 ----------
trend_poly3 <- make_trend(make_base('m', 'add', make_kernel('covariate1', 'poly3')))
kernel_pars <- c('m.d1'=1, 'm.d2'=1, 'm.d3'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_poly3, covariate1 = covariate1)
expected_output <- kernel_pars[[1]]*covariate1+kernel_pars[[2]]*covariate1^2+kernel_pars[[3]]*covariate1^3
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("poly3_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})

# poly4: Quartic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3 + d4 * c^4 ---------
trend_poly4 <- make_trend(make_base('m', 'add', make_kernel('covariate1', 'poly4')))
kernel_pars <- c('m.d1'=1, 'm.d2'=1, 'm.d3'=1, 'm.d4'=.1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_poly4, covariate1 = covariate1)
expected_output <- kernel_pars[[1]]*covariate1+kernel_pars[[2]]*covariate1^2+kernel_pars[[3]]*covariate1^3+kernel_pars[[4]]*covariate1^4
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("poly4_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})


# LEARNING RULES ----------------------------------------------------------
# include NA check

# delta --------
trend_delta <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'delta')))
kernel_pars <- c('m.q0'=0.5, 'm.alpha'=qnorm(0.2))
covariate1 <- c(NA, 1, NA, 1, NA)
emc <- make_minimal_emc(trend_delta, covariate1 = covariate1)
expected_output <- matrix(c(0.5, 0.5, 0.6, 0.6, 0.68)) # manually computed for this specific covariate + parameter vector
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("delta_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})


# delta 2 LR --------
trend_delta2lr <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'delta2lr')))
kernel_pars <- c('m.q0'=0.5, 'm.alphaPos'=qnorm(0.5), 'm.alphaNeg' = qnorm(0.25))
covariate1 <- c(NA, 1, NA, 0, NA)
emc <- make_minimal_emc(trend_delta2lr, covariate1 = covariate1)
expected_output <- matrix(c(0.5, 0.5, 0.75, 0.75, 0.5625)) # manually computed for this specific covariate + parameter vector
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("delta2lr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})


# delta 2 kernel ----------------------------------------------------------
# trend_delta2kernel <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'delta2kernel', base='lin')
# kernel_pars <- c('m.q0'=0.8, 'm.alphaFast'=qnorm(0.50), 'm.propSlow' = qnorm(0.10), 'm.dSwitch'=qnorm(0.1))
# covariate1 <- c(1, 1, 1, 1, 0, 0, 0, 0)
# ## use delta rules to construct target output.
# # this test is contingent on delta rules being implemented correctly,
# # but if that's not the case, an error is thrown above
# delta_fast <- EMC2:::run_delta(q0=rep(kernel_pars[1], length(covariate1)),
#                                alpha=rep(kernel_pars[2], length(covariate1)),
#                                covariate = covariate1)
# delta_slow <- EMC2:::run_delta(q0=rep(kernel_pars[1], length(covariate1)),
#                                alpha=rep(kernel_pars[2]*kernel_pars[3], length(covariate1)),
#                                covariate = covariate1)
# do_switch <- abs(delta_fast - delta_slow) > kernel_pars[4]
# expected_output <- delta_slow
# expected_output[do_switch] <- delta_fast[do_switch]
# emc <- make_minimal_emc(trend_delta2kernel, n_trials=8, covariate1 = covariate1)
# all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
# test_that("delta2kernel_Rcpp", {
#   expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
# })


## test resetting
trend_delta <- make_trend(make_base('m', 'lin', make_kernel('covariate1', 'delta', kernel_args = list(q_reset_column='do_reset'))))
kernel_pars <- c('m.q0'=0, 'm.alpha'=qnorm(0.2))
covariate1 <- c(1, 1, 1, 1, 1)
emc <- make_minimal_emc(trend_delta, covariate1 = covariate1, do_reset=c(F,F,T,F,F))
expected_output <- matrix(c(0, 0.2, 0, 0.2, 0.36)) # manually computed for this specific covariate + parameter vector
all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))
test_that("delta_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc)))
})



## Rescorla-wagner
trend_rw <-make_trend(make_base('m', 'add', make_kernel(c('covariate1','covariate2','covariate3','covariate4'),
                                                        'rescorlawagner', kernel_args = list(q_reset_column='do_reset'))))

# make_trend(par_names = "m",
#                        cov_names = list(c('covariate1','covariate2','covariate3','covariate4')),
#                        kernels = 'rescorlawagner', base='add',
#                        kernel_args=list(q_reset_column='do_reset'))

covariate_matrix <- matrix(c(NA,NA,  1,  1,
                             1,  1, NA, NA,
                             NA, 0,  0, NA,
                             1,  NA, NA, 1,
                             0,  0,  0,  0), nrow=5, byrow=TRUE)
covariates <- data.frame(covariate_matrix)
colnames(covariates) <- paste0('covariate', 1:4)
covariates$do_reset <- rep(FALSE, 5)
emc <- make_minimal_emc(trend_delta,
                        covariate1 = covariates$covariate1,
                        covariate2 = covariates$covariate2,
                        covariate3 = covariates$covariate3,
                        covariate4 = covariates$covariate4,
                        do_reset=c(F,F,F,T,F))
kernel_pars <- c('m.q0'=0, 'm.alpha'=0.2)
apply_kernel(kernel_pars, emc)
## looks ok...
# all.equal(matrix(apply_kernel(kernel_pars, emc)), matrix(expected_output))


# BETA-BINOMIAL, DBM, TPM LEARNING RULES --------------------------------------

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
    n_trials = length(covariate1),
    covariate1 = covariate1
  )
}

# TODO: add non-default output codes
# output_codes <- as.integer(1:4)
# output_types <- c("mean", "mode", "surprise", "log-precision")
output_codes <- 1L
output_types <- c("mean")

# Beta-binomial (basic)
emc <- lapply(
  X = output_codes,
  FUN = function(x) {
    make_minimal_emc_wrap("beta_binomial", output_code = x)
  }
)
names(emc) <- output_types
snapshot_list <- vector(mode = "list", length = 2L)
names(snapshot_list) <- c("uniform", "informative")

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
all.equal(emc_out, true_out)
snapshot_list[["uniform"]] <- list(
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
all.equal(emc_out, true_out)
snapshot_list[["informative"]] <- list(
  manual = true_out,
  emc = emc_out
)

# conclude with test snapshot
test_that("beta_binomial_Rcpp", {expect_snapshot(snapshot_list)})


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
all.equal(emc_out, true_out)
# conclude with test snapshot
test_that("beta_binomial_decay_Rcpp", {
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
all.equal(emc_out, true_out)
# conclude with test snapshot
test_that("beta_binomial_window_Rcpp", {
  expect_snapshot(list(manual = true_out, emc = emc_out))
})


# Dynamic Belief Model
emc <- lapply(
  X = output_codes,
  FUN = function(x) {
    make_minimal_emc_wrap("dbm", output_code = x)
  }
)
names(emc) <- output_types
snapshot_list <- vector(mode = "list", length = 3L)
names(snapshot_list) <- c("cp_0", "cp_1", "cp_real")

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
all.equal(emc_out, true_out)

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

all.equal(emc_out, emc_out_bb, tolerance = 1e-6)

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
all.equal(emc_out, true_out)

fixed_beta_out <- matrix(nrow = length(covariate1), ncol = length(output_types))
colnames(fixed_beta_out) <- output_types
fixed_beta_out[ , "mean"] <- rep(
  x = beta_mean(
    a = unname(pnorm(kernel_pars[2])) * unname(exp(kernel_pars[3])),
    b = (1 - unname(pnorm(kernel_pars[2]))) * unname(exp(kernel_pars[3]))
  ),
  length = length(covariate1)
)
# fixed_beta_out[ , "mode"] <- rep(
#   x = beta_mode(
#     a = unname(pnorm(kernel_pars[2])) * unname(exp(kernel_pars[3])),
#     b = (1 - unname(pnorm(kernel_pars[2]))) * unname(exp(kernel_pars[3]))
#   ),
#   length = length(covariate1)
# )
# fixed_beta_out[ , "surprise"] <- shannon_surprise(fixed_beta_out[ , "mean"], covariate1)
# fixed_beta_out[ , "log-precision"] <- rep(
#   x = beta_logprecision(
#     a = unname(pnorm(kernel_pars[2])) * unname(exp(kernel_pars[3])),
#     b = (1 - unname(pnorm(kernel_pars[2]))) * unname(exp(kernel_pars[3]))
#   ),
#   length = length(covariate1)
# )

all.equal(emc_out, fixed_beta_out, tolerance = 1e-5)

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
all.equal(emc_out, true_out)

snapshot_list[["cp_real"]] <- list(
  manual = true_out,
  emc = emc_out
)
# conclude with test snapshot
test_that("dbm_Rcpp", {expect_snapshot(snapshot_list)})


# # Custom kernel -- only C -------------------------------------------------
# This cannot be part of that-test :-( but we can still use it for manual tests...
# Write a custom kernel to a separate file
# tf <- tempfile(fileext = ".cpp")
# writeLines(c(
#   "// [[Rcpp::depends(EMC2)]]",
#   "#include <Rcpp.h>",
#   "#include \"EMC2/userfun.hpp\"",
#   "",
#   "// Example: two params (a, b) and two inputs (covariate1, t0)",
#   "Rcpp::NumericVector custom_kernel(Rcpp::NumericMatrix kernel_pars, Rcpp::NumericMatrix input) {",
#   "  int n = input.nrow();",
#   "  Rcpp::NumericVector out(n, 0.0);",
#   "  for (int i = 0; i < n; ++i) {",
#   "    double a = (kernel_pars.ncol() > 0) ? kernel_pars(i, 0) : 0.0;",
#   "    double b = (kernel_pars.ncol() > 1) ? kernel_pars(i, 1) : 0.0;",
#   "    double in1 = input(i, 0);  // covariate1",
#   "    double in2 = input(i, 1);  // t0",
#   "    if ((i % 2) == 0) out[i] = (Rcpp::NumericVector::is_na(in1) ? 0.0 : in1) + a;",
#   "    else              out[i] = (Rcpp::NumericVector::is_na(in2) ? 0.0 : in2) * b;",
#   "  }",
#   "  return out;",
#   "}",
#   "",
#   "// Export pointer maker for registration",
#   "// [[Rcpp::export]]",
#   "SEXP EMC2_make_custom_kernel_ptr();",
#   "EMC2_MAKE_PTR(custom_kernel)"
# ), tf)
#
# # Register with parameter names, transforms, and a default base
# ck <- register_kernel(
#   kernel_parameters = c("a", "b"),
#   file = tf,
#   transforms = c(a = "identity", b = "pnorm")
# #  base = "add"
# )
# kernel <- make_kernel('custom', 'covariate1', par_input='t0', custom_kernel = ck)
# trend_custom <- make_trend(make_base('add', target_parameter='m', kernel=kernel))
#

##
# kernel_pars <- c('m.a'=0.1, 'm.b'=1)
# covariate1 <- c(NA, 1, 2, 0, 20, 2, 1)
# t0 <- matrix(rep(0.2, length(covariate1)))
# colnames(t0) <- 't0'
# emc <- make_minimal_emc(trend_custom, covariate1 = covariate1, n_trials=7)
# expected_output <- matrix(NA, nrow=length(covariate1))
# a <- rep(kernel_pars['m.a'], length(covariate1))
# b <- rep(kernel_pars['m.b'], length(covariate1))
# for(i in 1:length(covariate1)) {
#   if((i%%2) == 1) {  # NB: Rcpp does 0-based indexing, so i%%2 == 0 corresponds to (i+1)%%2 here
#     expected_output[i] <- ifelse(is.na(covariate1[i]), 0, covariate1[i])+a[i]
#   } else {
#     expected_output[i] <- ifelse(is.na(t0[i]), 0, t0[i])*pnorm(b)[i] # note that we need to pnorm b here!
#   }
# }
# all.equal(matrix(apply_kernel(kernel_pars, emc, input_pars=t0)), matrix(expected_output))



