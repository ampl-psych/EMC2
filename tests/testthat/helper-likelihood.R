# Fixed coefficient draws on the sampled (unmapped) scale.
likelihood_parameters <- function(centre, n = 16, sd = .25) {
  withr::local_seed(123, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion")
  draws <- matrix(rnorm(n * length(centre), sd = sd), ncol = length(centre))
  colnames(draws) <- names(centre)
  rbind(unname(centre), sweep(draws, 2, centre, "+"))
}

model_likelihood <- function(emc, parameters) {
  vapply(emc[[1]]$data, function(dadm) {
    EMC2:::calc_ll_manager(parameters, dadm, emc[[1]]$model)
  }, numeric(nrow(parameters)))
}

circular_parameters <- list(
  CDM = c(v = log(1.1), theta = qnorm(.55), a = 0, t0 = log(.2), s = 0, sv = log(.1)),
  SDM = c(v = log(1.1), theta1 = qnorm(.35), theta2 = qnorm(.6), a = 0, t0 = log(.2), s = 0, sv = log(.1)),
  HSDM = c(v = 0, theta1 = qnorm(.32), theta2 = qnorm(.57), theta3 = qnorm(.66), a = 0, t0 = log(.2), s = 0, sv = log(.1)),
  PSDM = c(v = 0, theta = qnorm(.38), a = 0, t0 = log(.2), s = 0, sv = log(.1)),
  PHSDM = c(v = 0, theta1 = qnorm(.33), theta2 = qnorm(.58), a = 0, t0 = log(.2), s = 0, sv = log(.1))
)

circular_design <- function(name) {
  design(model = get(name, asNamespace("EMC2")), factors = list(subjects = 1),
    Rlevels = "resp", report_p_vector = FALSE,
    formula = lapply(names(circular_parameters[[name]]), function(p) reformulate("1", p)))
}
