RNGkind("L'Ecuyer-CMRG")
set.seed(123)

ADmat_sem <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
matchfun_sem <- function(d) d$S == d$lR

dat_sem <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:3], ]
dat_sem$subjects <- droplevels(dat_sem$subjects)
subject_scores <- setNames(c(-1, 0, 1), levels(dat_sem$subjects))
dat_sem$impulsiveness <- subject_scores[dat_sem$subjects]

design_sem <- design(
  data = dat_sem,
  model = LNR,
  matchfun = matchfun_sem,
  formula = list(m ~ lM, s ~ 1, t0 ~ 1),
  contrasts = list(m = list(lM = ADmat_sem))
)

test_that("make_sem_structure allows covariates as factor indicators", {
  sem_settings <- make_sem_structure(
    data = dat_sem,
    design = design_sem,
    covariate_cols = "impulsiveness",
    lambda_specs = list(Speed = c("m", "m_lMd", "impulsiveness"))
  )

  expect_equal(sem_settings$Lambda_mat["m", "Speed"], 1)
  expect_equal(sem_settings$Lambda_mat["m_lMd", "Speed"], Inf)
  expect_equal(sem_settings$Lambda_cov_mat["impulsiveness", "Speed"], Inf)
  expect_equal(colnames(sem_settings$lambda_covariates), "impulsiveness")
  expect_equal(mean(sem_settings$lambda_covariates$impulsiveness), 0)
})

test_that("make_sem_structure rejects overlap between lambda indicators and K/G regressors", {
  expect_error(
    make_sem_structure(
      data = dat_sem,
      design = design_sem,
      covariate_cols = "impulsiveness",
      lambda_specs = list(Speed = c("m", "impulsiveness")),
      k_specs = list(t0 = "impulsiveness")
    ),
    "cannot also appear in 'k_specs'"
  )
})

test_that("SEM group step conditions alpha means on observed covariate indicators", {
  subjects <- c("S1", "S2", "S3")
  pars <- c("a_Speed1", "a_Speed2")
  sem_settings <- list(
    Lambda_mat = matrix(c(1, Inf), nrow = 2, ncol = 1,
                        dimnames = list(pars, "Speed")),
    Lambda_cov_mat = matrix(Inf, nrow = 1, ncol = 1,
                            dimnames = list("impulsiveness", "Speed")),
    B_mat = matrix(0, nrow = 1, ncol = 1, dimnames = list("Speed", "Speed")),
    K_mat = matrix(0, nrow = 2, ncol = 0, dimnames = list(pars, NULL)),
    G_mat = matrix(0, nrow = 1, ncol = 0, dimnames = list("Speed", NULL)),
    factor_groups = 1L,
    covariates = data.frame(row.names = subjects),
    lambda_covariates = data.frame(
      impulsiveness = c(-1, 0, 1),
      row.names = subjects
    )
  )

  sampler <- list(
    nuisance = c(FALSE, FALSE),
    par_names = pars,
    n_subjects = length(subjects)
  )
  sampler <- add_info_SEM(sampler, sem_settings = sem_settings)
  sampler$samples <- sample_store_SEM(
    data = data.frame(subjects = factor(subjects, levels = subjects)),
    par_names = pars,
    iters = 1,
    is_nuisance = sampler$nuisance,
    sem_settings = sampler$sem_settings
  )

  start <- get_startpoints_SEM(sampler, start_mu = c(0, 0), start_var = diag(2))
  sampler$samples$theta_mu[, 1] <- start$tmu
  sampler$samples$theta_var[, , 1] <- start$tvar
  sampler$samples$lambda[, , 1] <- start$lambda
  sampler$samples$lambda_cov[, , 1] <- start$lambda_cov
  sampler$samples$B[, , 1] <- start$B
  sampler$samples$K[, , 1] <- start$K
  sampler$samples$G[, , 1] <- start$G
  sampler$samples$epsilon_inv[, 1] <- start$epsilon_inv
  sampler$samples$epsilon_cov_inv[, 1] <- start$epsilon_cov_inv
  sampler$samples$delta_inv[, , 1] <- start$delta_inv
  sampler$samples$eta[, , 1] <- start$eta
  sampler$samples$idx <- 1

  alpha <- matrix(
    c(-0.5, 0.1,
      0.0, -0.1,
      0.6, 0.2),
    nrow = 2
  )

  set.seed(1)
  step <- gibbs_step_SEM(sampler, alpha)

  expect_equal(dim(step$subj_mu), c(2, 3))
  expect_equal(dim(step$proposal_var), c(2, 2))
  expect_false(isTRUE(all.equal(step$subj_mu[, 1], step$subj_mu[, 3])))
})
