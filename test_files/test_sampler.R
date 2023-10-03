rm(list=ls())
print(load("test_files/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)
head(dat)
# NB: This data has been truncated at 0.25s and 1.5s

# Average rate = intercept, and rate d = difference (match-mismatch) contrast
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
ADmat

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat

# Here we fit a series of models
# We'll first build several plausible models and then do model selection
# Only B affected by E
design_B <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Flist=list(v~lM,sv~1,B~E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)

# prior <- get_prior_single(design = design_B)
# plot_prior(prior)

prior <- list(
  theta_mu_mean = 1:7,
  theta_mu_var = diag(c(7:1))
) # This way we're using default priors for the nuisance parameters

dat2 <- dat[dat$subjects %in% unique(dat$subjects)[1:4],]
dat2$subjects <- droplevels(dat2$subjects)

# Nuisance non hyper = non hierarchically estimated parameters
devtools::load_all()
samplers <- make_samplers(dat2, design_B, type = "blocked", prior = prior, par_groups = 1:7)
debug(EMC2:::init)
samplers <- run_samplers(samplers, stage= "preburn", iter = 20, step_size = 10, max_gd = 1.1, verbose = T, cores_for_chains = 1, cores_per_chain = 1)

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

devtools::load_all()
debug(run_stage)
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


