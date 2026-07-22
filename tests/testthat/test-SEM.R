ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
matchfun=function(d)d$S==d$lR

dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:3],]
dat$subjects <- droplevels(dat$subjects)

make_LNR_design <- function(){
  design(data = dat, model = LNR, matchfun = matchfun,
         formula = list(m~lM, s~1, t0~1),
         contrasts = list(m=list(lM=ADmat)), report_p_vector = FALSE)
}

# Lambda_mat and K_mat are applied to the parameter vector by position, so a row
# order that disagrees with sampled_pars() silently estimates loadings on the
# wrong parameters unless add_info_SEM realigns them.
test_that("SEM row constraints are aligned by name, not position", {
  design_LNR <- make_LNR_design()
  sem_settings <- make_sem_structure(data = dat, design = design_LNR,
                                     lambda_specs = list(F1 = c("m","m_lMd","s","t0")))

  permuted <- sem_settings
  permuted$Lambda_mat <- sem_settings$Lambda_mat[c("t0","s","m_lMd","m"), , drop = FALSE]

  SEM <- make_emc(dat, design_LNR, type = "SEM", sem_settings = permuted,
                  n_chains = 1, compress = FALSE)

  expect_identical(rownames(SEM[[1]]$sem_settings$Lambda_mat), SEM[[1]]$par_names)
  expect_identical(SEM[[1]]$sem_settings$Lambda_mat, sem_settings$Lambda_mat)
})

test_that("SEM row constraints follow the fitted joint design order", {
  # make_sem_structure gets the joint design list in a different order than the
  # one being fitted, which is how the row permutation arises in practice.
  designs_fitted <- list(taskA = make_LNR_design(), taskB = make_LNR_design())
  designs_sem <- designs_fitted[c("taskB", "taskA")]

  sem_settings <- make_sem_structure(data = dat, design = designs_sem,
                                     lambda_specs = list(F1 = c("taskA|m", "taskB|m")))
  expect_identical(rownames(sem_settings$Lambda_mat), names(sampled_pars(designs_sem)))

  SEM <- make_emc(list(taskA = dat, taskB = dat), designs_fitted, type = "SEM",
                  sem_settings = sem_settings, n_chains = 1, compress = FALSE)

  Lambda_mat <- SEM[[1]]$sem_settings$Lambda_mat
  expect_identical(rownames(Lambda_mat), SEM[[1]]$par_names)
  # The loadings named in lambda_specs must remain the free/fixed ones after alignment
  expect_equal(Lambda_mat["taskA|m", "F1"], 1)
  expect_equal(Lambda_mat["taskB|m", "F1"], Inf)
  expect_true(all(Lambda_mat[setdiff(rownames(Lambda_mat), c("taskA|m","taskB|m")), "F1"] == 0))
})

# A Lambda_mat with no Inf entries is a valid spec: the loadings are all fixed
# (e.g. at 1) and only the factor variances/covariances are estimated.
test_that("constrain_lambda keeps fixed non-zero loadings", {
  Lambda_mat <- matrix(c(1, 1, 0, 0), ncol = 1,
                       dimnames = list(c("m","m_lMd","s","t0"), "F1"))
  lambda <- array(rnorm(4 * 1 * 3), dim = c(4, 1, 3),
                  dimnames = list(rownames(Lambda_mat), "F1", NULL))

  out <- constrain_lambda(lambda, Lambda_mat)

  expect_true(all(out[c("m","m_lMd"), 1, ] == 1))
  expect_true(all(out[c("s","t0"), 1, ] == 0))
})

test_that("constrain_lambda leaves estimated entries alone and zeroes fixed zeroes", {
  Lambda_mat <- matrix(c(1, Inf, 0, Inf), ncol = 1,
                       dimnames = list(c("m","m_lMd","s","t0"), "F1"))
  lambda <- array(rnorm(4 * 1 * 3), dim = c(4, 1, 3),
                  dimnames = list(rownames(Lambda_mat), "F1", NULL))

  out <- constrain_lambda(lambda, Lambda_mat)

  expect_true(all(out["m", 1, ] == 1))
  expect_true(all(out["s", 1, ] == 0))
  expect_identical(out[c("m_lMd","t0"), 1, ], lambda[c("m_lMd","t0"), 1, ])
})

test_that("constrain_lambda is unchanged for Inf/0 constraint matrices", {
  # variant_standard passes a par_groups mask with no fixed non-zero entries
  constraintMat <- matrix(0, 4, 4)
  constraintMat[1:2, 1:2] <- Inf
  vars <- array(rnorm(4 * 4 * 2), dim = c(4, 4, 2))

  out <- constrain_lambda(vars, constraintMat)

  expect_identical(out[1:2, 1:2, ], vars[1:2, 1:2, ])
  expect_true(all(out[3:4, , ] == 0))
  expect_true(all(out[, 3:4, ] == 0))
})

test_that("SEM constraints with unknown parameter names are rejected", {
  design_LNR <- make_LNR_design()
  sem_settings <- make_sem_structure(data = dat, design = design_LNR,
                                     lambda_specs = list(F1 = c("m","m_lMd","s","t0")))
  rownames(sem_settings$Lambda_mat) <- c("m","m_lMd","s","NOT_A_PAR")

  expect_error(
    make_emc(dat, design_LNR, type = "SEM", sem_settings = sem_settings,
             n_chains = 1, compress = FALSE),
    "Lambda_mat rownames do not match the sampled parameters"
  )
})

test_that("SEM constraints with the wrong number of rows are rejected", {
  design_LNR <- make_LNR_design()
  sem_settings <- make_sem_structure(data = dat, design = design_LNR,
                                     lambda_specs = list(F1 = c("m","m_lMd","s","t0")))
  sem_settings$Lambda_mat <- sem_settings$Lambda_mat[1:3, , drop = FALSE]

  expect_error(
    make_emc(dat, design_LNR, type = "SEM", sem_settings = sem_settings,
             n_chains = 1, compress = FALSE),
    "Lambda_mat has 3 rows but the model has 4 sampled parameters"
  )
})
