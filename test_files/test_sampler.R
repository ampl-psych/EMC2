rm(list=ls())
devtools::load_all()
print(load("test_files/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)
head(dat)
dat$age <- rnorm(nrow(dat))
# NB: This data has been truncated at 0.25s and 1.5s

# Average rate = intercept, and rate d = difference (match-mismatch) contrast
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
ADmat

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat

undebug(make_dm)
design_B <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Fcovariates = list(age = dat$age),
  Flist=list(v~age*E,sv~1,B~age,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)

p_vector <- sampled_p_vector(design_B)
p_vector[1:length(p_vector)] <- abs(rnorm(length(p_vector)))

mapped_par(p_vector, design_B)

# prior <- get_prior_single(design = design_B)
# plot_prior(prior)

prior <- list(
  theta_mu_mean = 1:5,
  theta_mu_var = diag(c(5:1))
) # This way we're using default priors for the nuisance parameters

dat2 <- dat[dat$subjects %in% unique(dat$subjects)[1:4],]
dat2$subjects <- droplevels(dat2$subjects)

# Nuisance non hyper = non hierarchically estimated parameters
samplers_test <- make_samplers(dat2, design_B, type = "single")
samplers_test <- run_emc(samplers_test, cores_per_chain = 4, verbose = T)

microbenchmark::microbenchmark(
  EMC2:::run_bridge_sampling(samplers_test, cores_per_prop = 4, cores_for_props = 2),
  EMC2:::run_bridge_sampling(samplers_test, cores_per_prop = 2, cores_for_props = 4), times = 5
)



devtools::load_all()
attr(samplers[[1]], "variant_funs") <- EMC2:::get_variant_funs("diagonal")
debug(run_bridge_sampling)
EMC2:::run_bridge_sampling(samplers, subfilter = 300)
debug(compare_IC)
test <- compare_IC(list(s1 = samplers_good, s2 = samplers), BayesFactor = T, cores_per_prop = 3)



debug(gd_pmwg)
undebug(as_Mcmc)
test <- gd_pmwg(samplers, selection = "random", filter = "adapt")

test2 <- gd_pmwg(samplers_test, selection = "alpha", filter = "preburn")




debug(IS2)
samplers <- run_IS2(samplers, IS_samples = 50, n_cores = 14)

undebug(plot_density)
test <- plot_density(samplers, selection = "mu")

debug(IS2)
samplers <- run_IS2(samplers, filter = "burn")


debug(test_adapted)
samplers <- run_adapt(samplers, cores_for_chains = 3, cores_per_chain = 2, verbose = T, min_unique = 50)
samplers <- run_sample(samplers, cores_for_chains = 3, cores_per_chain = 2, verbose = T)
# samplers <- run_emc(samplers, cores_per_chain = 5, cores_for_chains = 1, verbose = T)
?run_sample

# nuisance = hierarchically estimated parameters, but no covariances or other relationships estimated
# grouped pars = pars estimated the same across participants.
samplers <- make_samplers(dat, design_B, nuisance = c(6,7), grouped_pars = 5, type = "infnt_factor")

# we could also specify a prior for these:
prior <- list(
  theta_mu_mean = 1:4,
  theta_mu_var = c(4:1), # Type infinite factor requires a vector of variances not a matrix. I'll make better docs about this.
  prior_nuis = list(
    theta_mu_mean = 1:2,
    theta_mu_var = rep(1,2)
  ),
  prior_grouped = list(
    theta_mu_mean = .5,
    prior_var = .3
  )
)
samplers <- make_samplers(dat, design_B, nuisance_non_hyper = c("B", "t0"), grouped_pars = "v", type = "infnt_factor")

samplers <- auto_burn(samplers, cores_per_chain = 4, cores_for_chains = 3, verbose = T)
debug(test_adapted)
samplers <- run_adapt(samplers, cores_per_chain = 4, verbose = T, min_unique = 100)

plot_chains(samplers, filter = "burn")

plot_acfs(samplers)

test <- plot_density(samplers)
test <- plot_density(samplers, selection = "mu")
check_run(samplers, interactive = F)

pp <- post_predict(samplers, n_cores = 12)

plot_fit(dat, pp)

plot_defective_density(dat)


