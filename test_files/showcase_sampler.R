rm(list = ls())

library(EMC2)

# Renames and Overview -----------------------------------------------------------------
### First arguments
# filter -> stage
# subfilter -> filter
# samplers -> emc
# variance -> sigma2
# full variance covariance matrix -> Sigma

# Function renames
# run_emc -> fit
# check_run -> check
# parameters_data_frame -> parameters
# plot_recovery -> recovery
# post_predict -> predict
# savage_dickey -> hypothesis
# p_test -> credible
# make_design -> design
# make_prior -> prior
# make_samplers -> make_emc


# New/massively changed functions
# get_data
# gd_summary
# ess_summary
# posterior_summary
# profile_plot
# subset
# get_pars #BACKBONE OF ALL (the version used in ll version is now called
# get_pars_matrix)

# Got rid of:
# attributes model_list and data_list and a whole lot of backend functions/code

############# IMPORTANT ########################
#### Almost every argument of get_pars works for every new function, and almost
#### all can be arbitrarily combined ###

# OK let's get going:
# It pays off to look thoroughly at ?get_pars, defaults are mentioned in the
# function description
?get_pars

# let's go over each option separately focusing on those relevant for plotting.
# The type argument is only for other functions, and so is true_pars
# the get_pars return_mcmc argument is more for other functions and not
# applicable for plot, when TRUE it returns an mcmc object

# We will use a small built in emc object that only contains 50 ''sample''
# stage samples, so can't set differently here
plot(samples_LNR, selection = "mu", stage = "adapt")

# Samples can be sub-setted in three ways
plot(samples_LNR, selection = "mu", filter = 20) # will exclude the first 20
plot(samples_LNR, selection = "mu", filter = c(1,2,3,20:30)) # only takes these
plot(samples_LNR, selection = "mu", length.out = 20)

# They can also be thinned
plot(samples_LNR, selection = "mu", thin = 2)

# Use map = TRUE to map
plot(samples_LNR, selection = "mu", map = TRUE)

# The by subject argument determines how the plots/stats are presented for
# selection = 'alpha'
plot(samples_LNR, selection = "alpha", by_subject = FALSE)
plot(samples_LNR, selection = "alpha", by_subject = TRUE)

# merge_chains concatenates chains
plot(samples_LNR, selection = "alpha", merge_chains = T)

# subject picks out certain subjects
plot(samples_LNR, selection = "alpha", subject = 1:2, by_subject = T)
plot(samples_LNR, selection = "alpha", subject = 1:2, by_subject = F)
plot(samples_LNR, selection = "alpha", subject = "as1t", by_subject = F)

# flatten determines how 3-d samples are presented (excluding selection = 'alpha')
plot(samples_LNR, selection = "correlation", flatten = FALSE)
plot(samples_LNR, selection = "correlation", flatten = TRUE)

# remove_dup removes duplicate samples, notice how the panels of correlations
# gets smaller and smaller
plot(samples_LNR, selection = "correlation", remove_dup = TRUE)

# remove constants isn't relevant for plotting, but can be turned off
plot(samples_LNR, selection = "correlation", remove_constants = FALSE)

# use_par can be used to select certain parameters
plot(samples_LNR, selection = "mu", use_par = c("m", "s"))
# or sub-groups of parameters, her it will only plot 2 panels, the correlations
# with m and s
plot(samples_LNR, selection = "correlation", use_par = c("m", "s"))
# also works with flatten
plot(samples_LNR, selection = "correlation", use_par = c("s.m"), flatten = T)
# and with map
plot(samples_LNR, selection = "mu", map = T, use_par = "m_TRUE")

# ????
plot(samples_LNR, selection = "mu", map = T, use_par = "m_TRUE", add_recalculated = TRUE)

# chain we can use to only take some chains
plot(samples_LNR, selection = "mu", chain = 1:2)
plot(samples_LNR, selection = "mu", chain = 3)

# Plotting functions now also take graphical parameters of base R with a few
# exceptions (e.g., plot_fit)

# For example putting amny options together
plot(samples_LNR, selection = "alpha", chain = 1:2, map = TRUE,
     use_par = c("m_TRUE", "s"),subject = 1:2, filter = 20:50, length.out = 20,
     by_subject = TRUE, col = c("darkgreen", "purple"),
     lwd = 2, lty = 1, main = "wow this is flexible")

# Something we are still working on is having the plot arguments passable as a
# list, one for each plot. So for now you can't set separate ylim for each plot,
# for example.

# Ok having discussed most of the new argument functionality let's go over the new functions

# statistics --------------------------------------------------------------
# Only the functions plot, check, and summary take length(selection) > 1
check(samples_LNR) # plots the worst parameter of each selection
summ <- summary(samples_LNR) # By default prints the first subject
# all subjects are stored in here though
summ

# Default is mu
posterior_summary(samples_LNR)
posterior_summary(samples_LNR, selection = "correlation")
gd_summary(samples_LNR)
gd_summary(samples_LNR, selection = "correlation", stat = "mean")
# can be any function that takes a vector, even if they don't make sense ;)
gd_summary(samples_LNR, selection = "correlation", stat = "sum")
gd_summary(samples_LNR, stat_only = T, selection = "Sigma")

# This function works the same
ess_summary(samples_LNR)


# plotting ----------------------------------------------------------------
# plot_pars now only returns the contraction
# to get the quantiles use `posterior_summary`
# to get the recovery use `recovery`
contract <- plot_pars(samples_LNR)
contract
plot_pars(samples_LNR, plot_prior = F)

# can also set random plotting arguments
plot_pars(samples_LNR, col = "orange", xlim = c(-2, 2), lty = 3, lwd = 2,
          prior_plot_args = list(col = "purple", lty = 1))

# With selection is alpha, the prior is now a rmvnorm draw using the group-level
# mean and group-level variance
plot_pars(samples_LNR, selection = "alpha")

# Can now also plot all subjects, which in this case is just 3
plot_pars(samples_LNR, selection = "alpha", all_subjects = T, use_par = c("s", "t0"))
# ... now takes arguments for plotting, get_pars and density
plot_pars(samples_LNR, selection = "alpha", all_subjects = T, adjust = 5)


# make_data, and design ---------------------------------------------------------------

# besides a p_vector and a p_mat make_data now also takes an emc object
# of which the medians and design will be used.
sim_dat <- make_data(samples_LNR)

pmedian <- posterior_summary(samples_LNR, probs = .5)[[1]]
# We can use these to do revamped profile plotting
design <- attr(samples_LNR, "design_list")[[1]]
profile_plot(sim_dat, design, pmedian)
# It's quite flexible
profile_plot(sim_dat, design, pmedian, use_par = c("m", "t0"),
             p_min = c(t0 = -1.5), p_max = c(t0 = -1.75))


# prior -------------------------------------------------------------------
# let's use that design to make a prior, here enter at the prompt to by
# setting fill_default to be F
prior_input <- prior(design, fill_default = F)

# or use ask to only fill those
prior_input <- prior(design, ask = "Sigma")

# We can also plot our prior, again the get_pars() and plot arguments work here
plot_prior(prior_input, design)
plot_prior(prior_input, design, use_par = c("s", "t0"), col = "green")


# recovery ----------------------------------------------------------------
# first make an emc object
sim_emc <- make_emc(sim_dat, design, prior_list = prior_input)
sim_emc <- fit(sim_emc, cores_per_chain = 3)

# use the new subset to make it equally long as samples_LNR
sim_emc <- subset(sim_emc, length.out = 50)

# Now we can use the improved recovery function
recovery(sim_emc, pmedian)

# Can also supply old samplers used to make_data
recovery(sim_emc, samples_LNR)
recovery(sim_emc, samples_LNR, selection = "alpha")
recovery(sim_emc, samples_LNR, selection = "alpha", col = "darkgreen")
recovery(sim_emc, samples_LNR, selection = "correlation", col = "darkgreen",
         ci_plot_args = list(angle = 60, col = "brown", lwd = .5))

# this now outputs the miss
recovered <- recovery(sim_emc, samples_LNR, selection = "alpha")
recovered
recovered <- recovery(sim_emc, samples_LNR, selection = "alpha", by_subject = T)
recovered

recovered <- recovery(sim_emc, samples_LNR, selection = "mu")
recovered

# or for a vector instead of emc object
recovered <- recovery(sim_emc, pmedian, selection = "mu")
recovered


# more plots --------------------------------------------------------------

# A ton of new plots coming soon




# design_specification ----------------------------------------------------

# New feature I'm working on will come here


