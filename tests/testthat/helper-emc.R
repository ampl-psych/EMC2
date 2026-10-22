# A four-subject, two-trial LNR with one free parameter: the smallest object
# the sampler runs through every stage.
tiny_lnr <- function(n_chains = 2){
  dat <- data.frame(subjects = factor(rep(1:4, each = 2)), R = factor(1),
                    rt = .2 + exp(c(-.7, -.2, .1, .4, .3, .8, -.4, .6)))
  des <- design(data = dat, model = LNR, formula = list(m ~ 1),
                constants = c(s = 0, t0 = log(.2)), report_p_vector = FALSE)
  make_emc(dat, des, type = "standard", n_chains = n_chains)
}

# Every stage in a few draws, 20 of them in sample
fit_tiny_lnr <- function(emc){
  criteria <- list(preburn = list(iter = 10), burn = list(iter = 10),
                   adapt = list(iter = 20, min_unique = 5), sample = list(iter = 20))
  fit(emc, cores_for_chains = 1, stop_criteria = criteria, verbose = FALSE,
      particle_factor = 10, step_size = 20, trim = FALSE)
}
