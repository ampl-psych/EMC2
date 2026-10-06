# Conjugate normal model on sufficient statistics, whose posterior is known: parameter j of subject s
# has likelihood N(ybar_sj | alpha_sj, 1 / T_j), i.e. precision T_j and linear term T_j * ybar_sj.
conj_ll <- function(pars, dadm, ...) -0.5 * sum(dadm$T * (dadm$ybar - pars)^2)

conj_emc_from <- function(ybar, T, pars = colnames(ybar)) {
  dat <- do.call(rbind, lapply(seq_len(nrow(ybar)), function(s)
    data.frame(subjects = s, par = pars, T = T, ybar = ybar[s, ])))
  dat$subjects <- factor(dat$subjects)
  des <- design(model = conj_ll, custom_p_vector = pars, report_p_vector = FALSE)
  make_emc(dat, des, type = "standard", n_chains = 2, compress = FALSE)
}

slow_tests <- function() Sys.getenv("EMC2_SLOW_TESTS") == "true"
