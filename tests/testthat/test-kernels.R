# Convenience functions ---------------------------------------------------------
make_tiny_design <- function(trend) {
  cov_names <- trend[[1]]$covariate
  return(design(factors = list(subjects = 1, S = 1), Rlevels = 1,
                covariates = cov_names, matchfun = function(x) x$S==x$R, trend = trend,
                formula = list(m ~ 1, s ~ 1, t0 ~ 1),
                model = LNR))
}

make_minimal_emc <- function(trend, n_trials=5, covariate1=1:5) {
  ## This *ONLY* creates an EMC object with the covariates. Kernel pars are ignored
  this_design <- make_tiny_design(trend)
  p_vector <- sampled_pars(this_design, doMap = FALSE)  # no kernel pars
  # n_kp <- length(kernel_pars)
  # if(n_kp>0) {
  #   p_vector[names(kernel_pars)] <- kernel_pars
  # }

  dat <- make_data(p_vector, this_design, n_trials = n_trials, covariates = data.frame(covariate1 = covariate1))
  emc <- make_emc(dat, this_design, type='single')
  return(emc)
}


# lin_incr ----------------------------------------------------------------
trend_lin_incr <- make_trend(par_names = "m",
                             cov_names = 'covariate1',
                             kernels = 'lin_incr')
kernel_pars <- NULL
emc <- make_minimal_emc(trend_lin_incr)
expected_output <- 1:5
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("lin_incr_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("lin_incr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})

# lin_decr ----------------------------------------------------------------
trend_lin_decr <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'lin_decr')
kernel_pars <- NULL
emc <- make_minimal_emc(trend_lin_decr)
expected_output <- -(1:5)
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("lin_decr_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("lin_decr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})

# exp_decr ----------------------------------------------------------------
trend_exp_decr <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'exp_decr')
kernel_pars <- c('m.d_ed'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_exp_decr, covariate1 = covariate1)
expected_output <- exp(-kernel_pars * covariate1)
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("exp_decr_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("exp_decr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})

# exp_incr ----------------------------------------------------------------
trend_exp_incr <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'exp_incr')
kernel_pars <- c('m.d_ei'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_exp_incr, covariate1 = covariate1)
expected_output <- 1-exp(-kernel_pars * covariate1)
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("exp_incr_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("exp_incr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})

# pow_decr ----------------------------------------------------------------
trend_pow_decr <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'pow_decr')
kernel_pars <- c('m.d_pd'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_pow_decr, covariate1 = covariate1)
expected_output = (1 + covariate1)^(-kernel_pars)
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("pow_decr_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("pow_decr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})

# pow_incr ----------------------------------------------------------------
trend_pow_incr <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'pow_incr')
kernel_pars <- c('m.d_pi'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_pow_incr, covariate1 = covariate1)
expected_output <- 1-(1 + covariate1)^(-kernel_pars)
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("pow_incr_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("pow_incr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})

# poly2: Quadratic polynomial: k = d1 * c + d2 * c^2 -----------
trend_poly2 <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'poly2', base='add')
kernel_pars <- c('m.d1'=1, 'm.d2'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_poly2, covariate1 = covariate1)
expected_output <- kernel_pars[[1]]*covariate1+kernel_pars[[2]]*covariate1^2
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("poly2_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("poly2_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})

# poly3: Cubic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3 ----------
trend_poly3 <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'poly3', base='add')
kernel_pars <- c('m.d1'=1, 'm.d2'=1, 'm.d3'=1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_poly3, covariate1 = covariate1)
expected_output <- kernel_pars[[1]]*covariate1+kernel_pars[[2]]*covariate1^2+kernel_pars[[3]]*covariate1^3
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("poly3_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("poly3_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})

# poly4: Quartic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3 + d4 * c^4 ---------
trend_poly4 <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'poly4', base='add')
kernel_pars <- c('m.d1'=1, 'm.d2'=1, 'm.d3'=1, 'm.d4'=.1)
covariate1 <- 1:5
emc <- make_minimal_emc(trend_poly4, covariate1 = covariate1)
expected_output <- kernel_pars[[1]]*covariate1+kernel_pars[[2]]*covariate1^2+kernel_pars[[3]]*covariate1^3+kernel_pars[[4]]*covariate1^4
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("poly4_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("poly4_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})


# LEARNING RULES ----------------------------------------------------------
# include NA check

# delta --------
trend_delta <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'delta', base='add')
kernel_pars <- c('m.q0'=0.5, 'm.alpha'=0.2)
covariate1 <- c(NA, 1, NA, 1, NA)
emc <- make_minimal_emc(trend_delta, covariate1 = covariate1)
expected_output <- matrix(c(0.5, 0.5, 0.6, 0.6, 0.68)) # manually computed for this specific covariate + parameter vector
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("delta_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("delta_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})


# delta 2 LR --------
trend_delta2lr <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'delta2lr', base='add')
kernel_pars <- c('m.q0'=0.5, 'm.alphaPos'=0.5, 'm.alphaNeg' = 0.25)
covariate1 <- c(NA, 1, NA, 0, NA)
emc <- make_minimal_emc(trend_delta2lr, covariate1 = covariate1)
expected_output <- matrix(c(0.5, 0.5, 0.75, 0.75, 0.5625)) # manually computed for this specific covariate + parameter vector
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output)) # WRONG - but why?!
test_that("delta2lr_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("delta2lr_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})


# delta 2 kernel ----------------------------------------------------------
trend_delta2kernel <- make_trend(par_names = "m", cov_names = 'covariate1', kernels = 'delta2kernel', base='add')
kernel_pars <- c('m.q0'=0.8, 'm.alphaFast'=0.50, 'm.propSlow' = 0.10, 'm.dSwitch'=0.1)
covariate1 <- c(1, 1, 1, 1, 0, 0, 0, 0)
## use delta rules to construct target output.
# this test is contingent on delta rules being implemented correctly,
# but if that's not the case, an error is thrown above
delta_fast <- EMC2:::run_delta(q0=rep(kernel_pars[1], length(covariate1)),
                               alpha=rep(kernel_pars[2], length(covariate1)),
                               covariate = covariate1)
delta_slow <- EMC2:::run_delta(q0=rep(kernel_pars[1], length(covariate1)),
                               alpha=rep(kernel_pars[2]*kernel_pars[3], length(covariate1)),
                               covariate = covariate1)
do_switch <- abs(delta_fast - delta_slow) > kernel_pars[4]
expected_output <- delta_slow
expected_output[do_switch] <- delta_fast[do_switch]
emc <- make_minimal_emc(trend_delta2kernel, n_trials=8, covariate1 = covariate1)
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="R")), matrix(expected_output))
all.equal(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")), matrix(expected_output))
test_that("delta2kernel_R", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="R")))
})
test_that("delta2kernel_Rcpp", {
  expect_snapshot(matrix(apply_kernel(kernel_pars, emc, mode="Rcpp")))
})



# Custom kernel -- only C -------------------------------------------------
# This cannot be part of that-test :-( but we can still use it for manual tests...
# # Write a custom kernel to a separate file
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
# ct <- register_trend(
#   trend_parameters = c("a", "b"),
#   file = tf,
#   transforms = c(a = "identity", b = "pnorm"),
#   base = "add"
# )
#
# # Use in a trend (note par_input to add t0 as an input column)
# trend_custom <- make_trend(
#   par_names  = "m",
#   cov_names  = "covariate1",
#   kernels    = "custom",
#   par_input  = list("t0"),
#   phase      = "pretransform",
#   bases      = NULL,         # uses ct$base (here: add)
#   custom_trend = ct
# )
#
# ##
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
#     expected_output[i] <- ifelse(is.na(t0[i]), 0, t0[i])*b[i]
#   }
# }
# all.equal(matrix(apply_kernel(kernel_pars, emc, input_pars=t0, mode="R")), matrix(expected_output))
# all.equal(matrix(apply_kernel(kernel_pars, emc, input_pars=t0, mode="Rcpp")), matrix(expected_output))
# test_that("customkernel_R", {
#   expect_snapshot(matrix(apply_kernel(kernel_pars, emc, input_pars=t0, mode="R")))
# })
# test_that("customkernel_Rcpp", {
#   expect_snapshot(matrix(apply_kernel(kernel_pars, emc, input_pars=t0, mode="Rcpp")))
# })
#




##  "    double a = (kernel_pars.ncol() > 0) ? kernel_pars(i, 0) : 0.0;",
#"    double b = (kernel_pars.ncol() > 1) ? kernel_pars(i, 1) : 0.0;",
#"    double in1 = input(i, 0);  // covariate1",
#"    double in2 = input(i, 1);  // t0",
#"    if ((i % 2) == 0) out[i] = (Rcpp::NumericVector::is_na(in1) ? 0.0 : in1) + a;",
#"    else              out[i] = (Rcpp::NumericVector::is_na(in2) ? 0.0 : in2) * b;",




# rm(list=ls())
# library(EMC2)
#
# apply_kernel <- function(kernel_pars, emc, subject=1, input_pars=NULL, trend_n=1, mode='Rcpp') {
#   ##
#   dadm <- emc[[1]]$data[[1]]
#   trend_list <- emc[[1]]$model()$trend
#   if(length(trend_list) > 1) {
#     warning(paste0('Multiple trends found - applying trend number ', trend_n))
#   }
#   trend <- trend_list[[trend_n]]
#   trend_par <- names(trend_list)[[trend_n]]
#
#   # extract kernel pars
#   if(trend$kernel %in% c("lin_incr", "lin_decr")) {
#     trend_pars <- matrix(0, nrow=nrow(dadm))
#     colnames(trend_pars) <- 'PLACEHOLDER'
#   } else {
#     trend_pars <- matrix(rep(kernel_pars, each=nrow(dadm)), ncol=length(kernel_pars), byrow=FALSE)
#     colnames(trend_pars) <- names(kernel_pars)
#   }
#
#   # Add base par -- first check if part of input_pars
#   if(trend_par %in% colnames(input_pars)) {
#     trend_pars <- cbind(input_pars[trend_par], trend_pars)
#     colnames(trend_pars)[1] <- trend_par
#   } else {
#     base_par <- trend_help(base=trend$base, do_return=TRUE)$default_pars
#     if(length(base_par) > 0) {
#       trend_pars <- cbind(0, trend_pars)
#       colnames(trend_pars)[1] <- trend$trend_pnames[1]
#     }
#   }
#
#
#   # all pars
#   pars_full <- trend_pars
#   if(!is.null(input_pars)) {
#     pars_full <- cbind(pars_full[,!colnames(pars_full)%in%colnames(input_pars)], input_pars)
#   }
#
#   # Define output parameter - 0 except when it's part of pars_full
#   if(trend_par %in% colnames(pars_full)) {
#     param <- pars_full[,trend_par]
#   } else {
#     param <- rep(0, nrow(dadm))
#   }
#
#
#   out_c <- out_R <- NULL
#   if(mode %in% c('Rcpp', 'compare')) {
#     out_c <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars, pars_full = pars_full, return_kernel = TRUE)
#   } else if(mode %in% c('R', 'compare')) {
#     out_R <- EMC2:::run_trend(dadm = dadm, trend=trend, param=param, trend_pars=trend_pars, pars_full = pars_full, return_kernel = TRUE)
#   }
#
#   if(mode == 'Rcpp') return(out_c)
#   if(mode == 'R') return(out_R)
#   if(mode == 'compare') all.equal(out_c, out_R)
# }
#
