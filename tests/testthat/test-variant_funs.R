for (variant in c("single", "standard", "diagonal", "blocked", "diagonal-gamma",
                  "SEM", "infnt_factor")) {
  test_that(paste(variant, "initializes and stores valid samples"), {
    withr::local_seed(123)
    dat <- data.frame(subjects = factor(rep(1:4, each = 3)),
                      C = factor(rep(1:3, 4)), R = factor(1),
                      rt = .2 + exp(c(-.6, .1, .8, -.4, .3, .5, -.2, .5, .9, -.3, .2, .7)))
    des <- design(data = dat, model = LNR, formula = list(m ~ 0 + C),
                  constants = c(s = 0, t0 = log(.2)), report_p_vector = FALSE)
    args <- switch(variant,
      blocked = list(par_groups = c(1, 1, 2)),
      SEM = list(sem_settings = make_sem_structure(data = dat, design = des,
                   lambda_specs = list(F1 = c("m_C1", "m_C2")))),
      infnt_factor = list(n_factors = 2),
      list())
    emc <- do.call(make_emc, c(list(dat, des, type = variant, n_chains = 2), args))
    emc <- init_chains(emc, particles = 5, cores_for_chains = 1)
    emc <- run_emc(emc, "preburn", stop_criteria = list(iter = 5),
                   particle_factor = 10, cores_for_chains = 1)
    for (chain in emc) {
      samples <- chain$samples
      expect_equal(dim(samples$alpha), c(3, 4, 6))
      expect_equal(dim(samples$subj_ll), c(4, 6))
      expect_equal(samples$idx, 6)
      expect_true(all(is.finite(samples$alpha)) && all(is.finite(samples$subj_ll)))
      if (variant == "single") next
      expect_equal(dim(samples$theta_mu), c(3, 6))
      expect_equal(dim(samples$theta_var), c(3, 3, 6))
      expect_true(all(is.finite(samples$theta_mu)))
      expect_true(all(vapply(2:6, function(i)
        min(eigen(samples$theta_var[, , i], symmetric = TRUE, only.values = TRUE)$values) > 0,
        logical(1))))
      if (variant %in% c("diagonal", "blocked", "diagonal-gamma")) {
        groups <- if (variant == "blocked") c(1, 1, 2) else 1:3
        covariance <- samples$theta_var[, , 2:6]
        expect_true(all(covariance[rep(outer(groups, groups, "!="), 5)] == 0))
      }
      if (variant == "SEM") {
        expect_true(all(samples$lambda[1, 1, ] == 1))
        expect_true(all(samples$lambda[3, 1, ] == 0))
      }
    }
  })
}
