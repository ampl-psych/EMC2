# Per-factor factor-disturbance variance priors (a_d / b_d).
# The sampler (gibbs_step_SEM), the prior sampler (get_prior_SEM) and the
# bridge density (bridge_group_and_prior_and_jac_SEM) must all index a_d/b_d
# per factor group and stay bit-identical to the scalar behaviour when a scalar
# (or a constant vector) is supplied.

ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
matchfun <- function(d) d$S == d$lR

make_sem_emc <- function(n_sub = 3, n_chains = 1) {
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[seq_len(n_sub)], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, matchfun = matchfun,
                formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                contrasts = list(m = list(lM = ADmat)), report_p_vector = FALSE)
  # F1 & F2 form a correlated group of size 2 (IW branch); F3 is a singleton
  # (inverse-gamma branch) -> both delta_inv branches are exercised.
  ss <- make_sem_structure(data = dat, design = des,
                           lambda_specs = list(F1 = c("m", "m_lMd"), F2 = c("s"), F3 = c("t0")),
                           factor_groups = list(c("F1", "F2")))
  make_emc(dat, des, type = "SEM", sem_settings = ss, n_chains = n_chains,
           compress = FALSE)
}

test_that("normalise_factor_var_prior recycles, matches names, and validates", {
  fg <- c(1, 1, 2); fn <- c("F1", "F2", "F3")
  # scalar recycles to n_factors, unnamed
  expect_identical(normalise_factor_var_prior(0.3, "b_d", 3, fn, fg), rep(0.3, 3))
  # length-n_factors vector passes through in column order
  expect_identical(normalise_factor_var_prior(c(1, 2, 3), "b_d", 3, fn, fg), c(1, 2, 3))
  # named vector is reordered to column order
  expect_identical(
    normalise_factor_var_prior(c(F3 = 0.01, F1 = 1, F2 = 0.5), "b_d", 3, fn, fg),
    c(1, 0.5, 0.01))
  # wrong length, negative, and NA all error informatively
  expect_error(normalise_factor_var_prior(c(1, 2), "b_d", 3, fn, fg), "length-3")
  expect_error(normalise_factor_var_prior(c(1, -2, 3), "b_d", 3, fn, fg), "positive")
  expect_error(normalise_factor_var_prior(c(1, NA, 3), "b_d", 3, fn, fg), "NA")
  # a_d must be homogeneous within the F1/F2 group
  expect_error(
    normalise_factor_var_prior(c(5, 8, 5), "a_d", 3, fn, fg, homogeneous_in_group = TRUE),
    "constant within each correlated")
  # b_d may differ within a group (heterogeneous IW scale diagonal is valid)
  expect_silent(normalise_factor_var_prior(c(1, 0.5, 0.01), "b_d", 3, fn, fg,
                                           homogeneous_in_group = FALSE))
})

test_that("make_emc normalises scalar a_d/b_d to per-factor vectors", {
  emc <- make_sem_emc()
  expect_length(emc[[1]]$prior$a_d, emc[[1]]$n_factors)
  expect_length(emc[[1]]$prior$b_d, emc[[1]]$n_factors)
  expect_true(all(emc[[1]]$prior$b_d == 0.3))
})

test_that("gibbs_step_SEM rejects malformed per-factor a_d/b_d", {
  # These fire during prior normalisation, before any sampling state is used,
  # so a freshly-made (unsampled) emc is enough to test them.
  emc <- make_sem_emc()
  s <- emc[[1]]
  n_pars <- sum(!s$nuisance)
  alpha <- matrix(0, n_pars, s$n_subjects,
                  dimnames = list(s$par_names[!s$nuisance], NULL))
  s_wrong <- s; s_wrong$prior$b_d <- c(1, 2)
  expect_error(gibbs_step_SEM(s_wrong, alpha), "length-3")
  s_neg <- s; s_neg$prior$b_d <- c(1, -2, 3)
  expect_error(gibbs_step_SEM(s_neg, alpha), "positive")
  s_het <- s; s_het$prior$a_d <- c(5, 8, 5)   # F1,F2 grouped
  expect_error(gibbs_step_SEM(s_het, alpha), "constant within each correlated")
})

test_that("bridge density is bit-identical for scalar vs constant-vector prior", {
  # Hand-built minimal info so this needs no fitting. Column layout of
  # proposals_group: theta_mu(3), eps_log(3), delta group1 chol(3), group2(1).
  n_pars <- 3; n_factors <- 3; n_subjects <- 4
  par_names <- c("p1", "p2", "p3"); factor_names <- c("F1", "F2", "F3")
  Lambda_mat <- matrix(0, n_pars, n_factors, dimnames = list(par_names, factor_names))
  diag(Lambda_mat) <- 1
  B_mat <- matrix(0, n_factors, n_factors, dimnames = list(factor_names, factor_names))
  K_mat <- matrix(0, n_pars, 0, dimnames = list(par_names, NULL))
  G_mat <- matrix(0, n_factors, 0, dimnames = list(factor_names, NULL))
  sem_settings <- list(Lambda_mat = Lambda_mat, B_mat = B_mat, K_mat = K_mat,
                       G_mat = G_mat, factor_groups = c(1, 1, 2),
                       covariates = data.frame(matrix(ncol = 0, nrow = n_subjects)))
  mk_info <- function(a_d, b_d) list(
    prior = list(theta_mu_mean = rep(0, n_pars), theta_mu_var = rep(1, n_pars),
                 lambda_var = rep(.7, n_factors), K_var = numeric(0),
                 G_var = rep(.5, n_factors), B_var = rep(.5, n_factors),
                 a_d = a_d, b_d = b_d, a_e = rep(5, n_pars), b_e = rep(.3, n_pars)),
    n_pars = n_pars, n_factors = n_factors, n_subjects = n_subjects,
    sem_settings = sem_settings)

  set.seed(99); n_iter <- 5
  proposals_group <- cbind(matrix(rnorm(n_iter * n_pars), n_iter),   # theta_mu
                           matrix(rnorm(n_iter * n_pars), n_iter),   # eps_log
                           matrix(rnorm(n_iter * 3), n_iter),        # delta group1
                           matrix(rnorm(n_iter * 1), n_iter))        # delta group2
  proposals_list <- list(matrix(rnorm(n_iter * n_pars * n_subjects), n_iter))

  d_scalar <- bridge_group_and_prior_and_jac_SEM(proposals_group, proposals_list, mk_info(5, 0.3))
  d_vector <- bridge_group_and_prior_and_jac_SEM(proposals_group, proposals_list,
                                                 mk_info(rep(5, 3), rep(0.3, 3)))
  expect_identical(d_scalar, d_vector)
  # heterogeneous a_d within the correlated group errors on the bridge side too
  expect_error(
    bridge_group_and_prior_and_jac_SEM(proposals_group, proposals_list, mk_info(c(5, 8, 5), 0.3)),
    "constant within each correlated")
})

test_that("gibbs_step_SEM runs with a non-constant per-factor lambda_var (byrow expansion)", {
  skip_on_cran()
  emc <- make_sem_emc(n_sub = 3, n_chains = 2)
  emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 2),
                 verbose = FALSE, cores_for_chains = 1)
  s <- emc[[1]]
  n_pars <- sum(!s$nuisance)
  set.seed(1)
  alpha <- matrix(rnorm(n_pars * s$n_subjects), n_pars, s$n_subjects,
                  dimnames = list(s$par_names[!s$nuisance], NULL))
  # Distinct per-factor loading-prior variances must be applied column-wise
  # (byrow), not scrambled by column-major recycling; the step should run and
  # return finite loadings.
  s$prior$lambda_var <- c(2, 0.5, 0.1)
  set.seed(7); d <- gibbs_step_SEM(s, alpha)
  expect_true(all(is.finite(d$lambda)))
})

test_that("gibbs_step_SEM draws are bit-identical for scalar vs constant-vector b_d", {
  skip_on_cran()
  emc <- make_sem_emc(n_sub = 3, n_chains = 2)
  emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 2),
                 verbose = FALSE, cores_for_chains = 1)
  s <- emc[[1]]
  n_pars <- sum(!s$nuisance)
  set.seed(1)
  alpha <- matrix(rnorm(n_pars * s$n_subjects), n_pars, s$n_subjects,
                  dimnames = list(s$par_names[!s$nuisance], NULL))

  s1 <- s; s1$prior$a_d <- 5;              s1$prior$b_d <- 0.3
  s2 <- s; s2$prior$a_d <- rep(5, s$n_factors); s2$prior$b_d <- rep(0.3, s$n_factors)
  set.seed(42); d1 <- gibbs_step_SEM(s1, alpha)
  set.seed(42); d2 <- gibbs_step_SEM(s2, alpha)
  expect_identical(d1, d2)

  # A genuinely per-factor b_d must change the delta_inv draw -> the vector path
  # actually reaches the sampler (stand-in for the full recovery test).
  s3 <- s; s3$prior$b_d <- c(1, 0.5, 0.01)
  set.seed(42); d3 <- gibbs_step_SEM(s3, alpha)
  expect_false(isTRUE(all.equal(d3$delta_inv, d1$delta_inv)))
  expect_true(all(is.finite(d3$delta_inv)))
})

# ---------------------------------------------------------------------------
# Non-zero loading prior means (lambda_mean) + the 3c bridge-side lambda_var
# unwind-order fix. lambda_mean is a conjugate Gaussian prior mean on the free
# loadings; lambda_mean = 0 must be bit-identical to the historical behaviour.
# ---------------------------------------------------------------------------

make_sem_emc_free <- function(n_sub = 4, n_chains = 1) {
  # F1 loads on m (fixed 1) + m_lMd (free) + s (free); F2 on t0 (fixed 1).
  # -> 2 free loadings across 2 factors, so per-factor and per-entry priors bite.
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[seq_len(n_sub)], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, matchfun = matchfun,
                formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                contrasts = list(m = list(lM = ADmat)), report_p_vector = FALSE)
  ss <- make_sem_structure(data = dat, design = des,
                           lambda_specs = list(F1 = c("m", "m_lMd", "s"), F2 = c("t0")))
  list(dat = dat, des = des, ss = ss,
       emc = make_emc(dat, des, type = "SEM", sem_settings = ss,
                      n_chains = n_chains, compress = FALSE))
}

# Minimal hand-built bridge info: 3 params, 2 factors, 3 free loadings, no
# covariates. Column layout of proposals_group: theta_mu(3), v_lambda(3),
# eps_log(3), delta group1(1), delta group2(1).
bridge_lambda_fixture <- function() {
  n_pars <- 3; n_factors <- 2; n_subjects <- 4
  pn <- c("p1", "p2", "p3"); fn <- c("F1", "F2")
  Lam <- matrix(0, n_pars, n_factors, dimnames = list(pn, fn))
  Lam[1, 1] <- 1; Lam[2, 1] <- Inf; Lam[3, 1] <- Inf; Lam[3, 2] <- Inf
  ss <- list(Lambda_mat = Lam,
             B_mat = matrix(0, n_factors, n_factors, dimnames = list(fn, fn)),
             K_mat = matrix(0, n_pars, 0, dimnames = list(pn, NULL)),
             G_mat = matrix(0, n_factors, 0, dimnames = list(fn, NULL)),
             factor_groups = c(1, 2),
             covariates = data.frame(matrix(ncol = 0, nrow = n_subjects)))
  mk_info <- function(lambda_mean = NULL, lambda_var = rep(.7, n_factors)) {
    pr <- list(theta_mu_mean = rep(0, n_pars), theta_mu_var = rep(1, n_pars),
               lambda_var = lambda_var, K_var = numeric(0), G_var = rep(.5, n_factors),
               B_var = rep(.5, n_factors), a_d = 5, b_d = .3,
               a_e = rep(5, n_pars), b_e = rep(.3, n_pars))
    if (!is.null(lambda_mean)) pr$lambda_mean <- lambda_mean
    list(prior = pr, n_pars = n_pars, n_factors = n_factors, n_subjects = n_subjects,
         sem_settings = ss)
  }
  set.seed(5); n_iter <- 6
  pg <- cbind(matrix(rnorm(n_iter * 3), n_iter), matrix(rnorm(n_iter * 3), n_iter),
              matrix(rnorm(n_iter * 3), n_iter), matrix(rnorm(n_iter * 1), n_iter),
              matrix(rnorm(n_iter * 1), n_iter))
  pl <- list(matrix(rnorm(n_iter * n_pars * n_subjects), n_iter))
  list(mk_info = mk_info, pg = pg, pl = pl, n_pars = n_pars, n_factors = n_factors)
}

test_that("normalise_lambda_matrix_prior handles scalar / per-factor / per-entry", {
  pn <- c("p1", "p2", "p3"); fn <- c("F1", "F2")
  # scalar -> whole matrix
  expect_equal(normalise_lambda_matrix_prior(1, "lambda_mean", 3, 2, pn, fn),
               matrix(1, 3, 2))
  # per-factor vector -> byrow (column f constant)
  expect_equal(normalise_lambda_matrix_prior(c(0.5, 2), "lambda_var", 3, 2, pn, fn, positive = TRUE),
               matrix(c(0.5, 0.5, 0.5, 2, 2, 2), 3, 2))
  # named per-factor vector reordered to column order
  expect_equal(normalise_lambda_matrix_prior(c(F2 = 2, F1 = 0.5), "lambda_var", 3, 2, pn, fn, positive = TRUE),
               matrix(c(0.5, 0.5, 0.5, 2, 2, 2), 3, 2))
  # per-entry matrix passes through; named dimnames reordered
  m <- matrix(c(0, 1, -1, 0, 0, 1), 3, 2, dimnames = list(pn, fn))
  expect_equal(normalise_lambda_matrix_prior(m, "lambda_mean", 3, 2, pn, fn), unname(m))
  expect_equal(normalise_lambda_matrix_prior(m[c(3, 1, 2), c(2, 1)], "lambda_mean", 3, 2, pn, fn),
               unname(m))
  # validation
  expect_error(normalise_lambda_matrix_prior(matrix(0, 3, 3), "lambda_mean", 3, 2, pn, fn),
               "must be 3 x 2")
  expect_error(normalise_lambda_matrix_prior(c(1, 2, 3), "lambda_mean", 3, 2, pn, fn),
               "scalar, a length-2")
  expect_error(normalise_lambda_matrix_prior(c(1, -2), "lambda_var", 3, 2, pn, fn, positive = TRUE),
               "positive")
  # lambda_mean allows negative values (signs are meaningful)
  expect_silent(normalise_lambda_matrix_prior(c(-1, 2), "lambda_mean", 3, 2, pn, fn))
})

test_that("get_prior_SEM defaults lambda_mean to 0", {
  emc <- make_sem_emc_free()$emc
  expect_identical(emc[[1]]$prior$lambda_mean, 0)
})

test_that("bridge density is bit-identical for lambda_mean = 0, absent, and matrix", {
  fx <- bridge_lambda_fixture()
  b_absent <- bridge_group_and_prior_and_jac_SEM(fx$pg, fx$pl, fx$mk_info(lambda_mean = NULL))
  b_zero   <- bridge_group_and_prior_and_jac_SEM(fx$pg, fx$pl, fx$mk_info(lambda_mean = 0))
  b_mat0   <- bridge_group_and_prior_and_jac_SEM(fx$pg, fx$pl,
                fx$mk_info(lambda_mean = matrix(0, fx$n_pars, fx$n_factors)))
  expect_identical(b_absent, b_zero)
  expect_identical(b_absent, b_mat0)
})

test_that("bridge lambda_var: per-factor vector == equivalent per-entry matrix (3c fix)", {
  fx <- bridge_lambda_fixture()
  b_vec <- bridge_group_and_prior_and_jac_SEM(fx$pg, fx$pl,
             fx$mk_info(lambda_mean = 0, lambda_var = c(0.5, 2)))
  b_mat <- bridge_group_and_prior_and_jac_SEM(fx$pg, fx$pl,
             fx$mk_info(lambda_mean = 0,
                        lambda_var = matrix(c(0.5, 0.5, 0.5, 2, 2, 2), fx$n_pars, fx$n_factors)))
  expect_identical(b_vec, b_mat)
  # a per-entry lambda_mean must change the density
  b_mean <- bridge_group_and_prior_and_jac_SEM(fx$pg, fx$pl,
              fx$mk_info(lambda_mean = matrix(c(0, 1, -1, 0, 0, 1), fx$n_pars, fx$n_factors)))
  expect_false(isTRUE(all.equal(b_mean, b_vec)))
})

test_that("gibbs_step_SEM rejects malformed lambda_mean (before sampling state is used)", {
  s <- make_sem_emc_free()$emc[[1]]
  n_pars <- sum(!s$nuisance)
  alpha <- matrix(0, n_pars, s$n_subjects, dimnames = list(s$par_names[!s$nuisance], NULL))
  s_dim <- s; s_dim$prior$lambda_mean <- matrix(0, n_pars, s$n_factors + 1)
  expect_error(gibbs_step_SEM(s_dim, alpha), paste0("must be ", n_pars, " x ", s$n_factors))
  s_len <- s; s_len$prior$lambda_mean <- c(1, 2, 3)   # not scalar, not n_factors
  expect_error(gibbs_step_SEM(s_len, alpha), "per-factor vector")
})

test_that("make_emc warns when lambda_mean is set but there are no free loadings", {
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:3], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, matchfun = matchfun,
                formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                contrasts = list(m = list(lM = ADmat)), report_p_vector = FALSE)
  ss <- make_sem_structure(data = dat, design = des,
                           lambda_specs = list(F1 = c("m"), F2 = c("t0")))  # all fixed
  expect_equal(sum(ss$Lambda_mat == Inf), 0)
  expect_warning(
    make_emc(dat, des, type = "SEM", sem_settings = ss,
             prior = prior(des, type = "SEM", sem_settings = ss, lambda_mean = 1),
             n_chains = 2, compress = FALSE),
    "no free .* loadings")
})

test_that("gibbs_step_SEM: lambda_mean = 0 bit-identical; vector==matrix; mean shifts loadings", {
  skip_on_cran()
  fr <- make_sem_emc_free(n_sub = 4, n_chains = 2)
  emc <- run_emc(fr$emc, "preburn", stop_criteria = list(iter = 2),
                 verbose = FALSE, cores_for_chains = 1)
  s <- emc[[1]]
  n_pars <- sum(!s$nuisance); n_factors <- s$n_factors
  set.seed(1)
  alpha <- matrix(rnorm(n_pars * s$n_subjects), n_pars, s$n_subjects,
                  dimnames = list(s$par_names[!s$nuisance], NULL))

  # Backward compatibility: 0, absent, and a constant-0 matrix all coincide.
  s0 <- s; s0$prior$lambda_mean <- 0
  sA <- s; sA$prior$lambda_mean <- NULL
  sM <- s; sM$prior$lambda_mean <- matrix(0, n_pars, n_factors)
  set.seed(42); d0 <- gibbs_step_SEM(s0, alpha)
  set.seed(42); dA <- gibbs_step_SEM(sA, alpha)
  set.seed(42); dM <- gibbs_step_SEM(sM, alpha)
  expect_identical(d0, dA)
  expect_identical(d0, dM)

  # per-factor lambda_var vector == equivalent per-entry matrix (Gibbs side)
  sv  <- s; sv$prior$lambda_var  <- c(0.5, 2)
  sMv <- s; sMv$prior$lambda_var <- matrix(c(0.5, 0.5, 0.5, 0.5, 2, 2, 2, 2), n_pars, n_factors)
  set.seed(7); dv  <- gibbs_step_SEM(sv, alpha)
  set.seed(7); dMv <- gibbs_step_SEM(sMv, alpha)
  expect_identical(dv, dMv)

  # A per-entry lambda_mean reaches the sampler (changes the loading draw).
  sMe <- s; sMe$prior$lambda_mean <- matrix(c(0, 1, -1, 0, 0, 0, 0, 0), n_pars, n_factors)
  set.seed(42); dMe <- gibbs_step_SEM(sMe, alpha)
  expect_false(isTRUE(all.equal(dMe$lambda, d0$lambda)))

  # Regulariser/conjugacy direction: shrinking toward 1 lifts the free-loading
  # posterior mean relative to the mean-0 prior (moment check over many draws).
  fi <- which(s$sem_settings$Lambda_mat == Inf, arr.ind = TRUE)[1, ]
  draw_free <- function(sx) sapply(1:300, function(z) {
    set.seed(z); gibbs_step_SEM(sx, alpha)$lambda[fi[1], fi[2]] })
  s_m0 <- s; s_m0$prior$lambda_mean <- 0; s_m0$prior$lambda_var <- 0.25
  s_m1 <- s; s_m1$prior$lambda_mean <- 1; s_m1$prior$lambda_var <- 0.25
  expect_gt(mean(draw_free(s_m1)), mean(draw_free(s_m0)))
})
