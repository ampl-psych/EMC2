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
