# Simple normal likelihood, the first half of the parameters are the means,
# the second half the variances. We exp them to ensure positive variances
custom_ll <- function(pars, dadm, ...){
  idx <- dadm$S == "left"
  ll_left <- dnorm(dadm$rt[idx], pars[1], exp(pars[3]))
  ll_right <- dnorm(dadm$rt[!idx], pars[2], exp(pars[4]))
  return(sum(log(c(ll_left, ll_right))))
}


library(EMC2)
pars <- c("muL", "muR", "sdL", "sdR")


design <- design(model = custom_ll, custom_p_vector = pars)

# using nuisance = 3:4
custom_emc <- make_emc(get_data(samples_LNR), design, type = "standard", nuisance_non_hyper = 3:4)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)
test_that("custom_ll", {
  expect_snapshot(init_chains(custom_emc, cores_for_chains = 1, particles = 10))
})
