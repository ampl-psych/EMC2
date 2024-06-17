devtools::load_all()
# First define a design:
design_DDMaE <- make_design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))

# prior <- make_prior_new(design_DDMaE)

# Then create a p_vector:
p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
          t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))

# debug(make_prior_new)
# debug(ask_user_prior)
# prior <- make_prior_new(design_DDMaE)

prior <- make_prior_new(design_DDMaE, pmean = p_vector[-1], psd = c("sv" = .5),
                        fill_default = TRUE)


devtools::load_all()
debug(plot_prior_new)
plot_prior_new(prior, design_DDMaE, selection = "covariance", N = 1e3, flatten = T)


qdesign_DDMaEv <- make_design(data = forstmann,model=DDM,
                            formula =list(v~0+S*E,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                            constants=c(s=log(1)))
# devtools::load_all()
# debug(make_prior_new)
prior_new <- make_prior_new(design_DDMaEv, pmean = p_vector, update = prior,
                        psd = c("v_Eneutral" = -100),
                        fill_default = TRUE)

plot_prior_new(prior_new, design_DDMaEv, selection = "variance", N = 5e4,
               mapped = F)


devtools::load_all()
debug(plot_prior_new)
plot_prior_new(prior, design_DDMaE, selection = "alpha", N = 1e3)

plot_prior_new(prior, design_DDMaE, selection = "covariance")
plot_prior_new(prior, design_DDMaE, selection = "mu")
plot_prior_new(prior, design_DDMaE, selection = "variance")
plot_prior(prior, design_DDMaE, selection = "correlation")

# This will map the parameters of the p_vector back to the design
# undebug(mapped_par)
# debug(make_data)
mapped_par(p_vector,design_DDMaE)


debug(as_mcmc_new)
plot_chains_new(samplers_LNR, selection = "LL")

debug(plot_pars_new)
plot_pars_new(samplers_LNR, plot_prior = F, selection = "LL")


devtools::load_all()
posterior_summary_new(samplers_LNR, selection = "alpha")



test <- es_summary_new(samplers_LNR, selection = "alpha", stat = "mean", stat_only = F)



debug(get_summary_stat)
gd_summary_new(samplers_LNR, selection = "alpha", omit_mpsrf = T)

debug(get_summary_stat)

devtools::load_all()

debug(summary.emc)
summary.emc(samplers_LNR)


debug(get_summary_stat)
get_summary_stat(samplers_LNR, fun = c(effectiveSize, get_posterior_quantiles), map = TRUE, stat_name = c("ESS", "1", "3"),
                 probs = c(.2, .5, .8))


posterior_summary_new(samplers_LNR, selection = "alpha", stat_name = 1:3)

debug(make_nice_summary)

# To do:
# Check add_constants
# Fix post_predict with hyper
# Email Michael, change schedule
# check run
devtools::load_all()
debug(savage_dickey)
savage_dickey(samplers_LNR, parameter = "m", selection = "alpha", subject = 1)


debug(plot_pars_new)
plot_pars_new(samplers_LNR, selection = "alpha", show_chains = TRUE, by_subject = FALSE, plot_prior = TRUE,
              use_prior_lim = TRUE, use_par = c("m", "s"))

debug(as_mcmc_new)

devtools::load_all()
debug(as_mcmc_new)
debug(as_mcmc_new)
plot_pars_new(samplers_LNR, selection = "alpha", show_chains = TRUE, by_subject = T,
              plot_prior = T, use_prior_lim = T, subject = 1:2, map = TRUE,
              use_par = c('t0', 'm_TRUE'))


plot_pars_new(samplers_LNR, selection = "covariance", show_chains = TRUE,
              plot_prior = T, use_prior_lim = T,
              use_par = c('t0'))



plot_chains(samplers_LNR, selection = "mu", subfilter = 1:40)



devtools::load_all()
debug(as_mcmc_new)

MCMC <- as_mcmc_new(samplers_LNR, selection = "covariance", use_par = c("s.m"), flatten = TRUE)

debug(filter_emc)
plot_chains_new(samplers_LNR, selection = "mu",flatten = TRUE, map = F,
                subfilter = seq(1, 40, by = 2))

devtools::load_all()
debug(get_objects_standard)
plot_pars_new(samplers_LNR, selection = "alpha", show_chains = FALSE)
plot_chains(samplers_LNR, selection = "mu")

debug(plot_chains)
plot_chains(samplers_LNR, selection = "alpha", layout = NULL)


