rm(list = ls())
library(EMC2)
set.seed(1)

# Set up a contrast, general lexical sensitivity and drift bias
LexMat <- cbind(d = c(-1, 1))

# See parameters and transformations of the DDM
?DDM

# Here we simulate for 1 subject
design_lex <- design(factors = list(Lex = c("Non-Word", "Word"), subjects = 1),
                     covariates = "Freq",
                     Rlevels = c("Non-Word", "Word"),
                     formula = list(v ~ Lex + Freq, a ~ 1, t0 ~ 1, Z ~ 1, sv ~ 1, s ~ Lex),
                     contrast = list(Lex = LexMat),
                     constants = c(s = log(1)),
                     model = DDM)

mapped_pars(design_lex)

sampled_pars(design_lex)

p_vector <- sampled_pars(design_lex)
p_vector[] <- c(.1, 1.5, .2, log(1.1), log(.3), qnorm(.55), log(.3), .2)


mapped_pars(design_lex, p_vector)

# plot a representation of the DDM for this design
plot_design(design_lex, p_vector = p_vector, factors = list(v = "Lex"))

dat <- make_data(parameters = p_vector, design = design_lex, n_trials = 100)
# We get a little message that our covariate is randomly imputed, let's fix that!
# To add covariates just add it as an argument to the data
# Usually frequency is quite tailed

# frequency is not defined for non-words
word_frequency <- rgamma(sum(dat$Lex == "Word"), 5, .1)
hist(word_frequency)

# To get rid of skew in a covariate we can use the log of it
word_frequency <- log(word_frequency)

# Also standardize to retain intercept as average effect
word_frequency <- scale(word_frequency)
# Create vector that stores frequency for words and non-words
frequency <- numeric(nrow(dat))
frequency[dat$Lex == "Word"] <- word_frequency

dat <- make_data(parameters = p_vector, design = design_lex, n_trials = 100,
                 covariates = list(Freq = frequency))

plot_density(dat, factors = "Lex")

# Priors be mindful of transformations!
?DDM

# Default priors are N(0,1)
# Priors on effects left at 0, except for drift rate, with the assumption
# that on average people are correct. Other priors seem kind of reasonable, but
# we urge you to think about prior expectations yourself!
prior_lex <- prior(pmean = c(v = 0, v_Lexd = 2, v_freq = 0, a = log(1), t0 = log(.25),
                             Z = qnorm(.5), sv = log(.3), s_Lexd = 0),
                   design = design_lex,
                   type = "single")

emc <- make_emc(dat, design_lex, prior_list = prior_lex, type = "single")
emc <- fit(emc, fileName = "DDM.RData")

summary(emc)
# Decent recovery, need more trials for t0 and s_Lexd
plot_pars(emc, true_pars = p_vector, use_prior_lim = F)

pp <- predict(emc)
plot_cdf(dat, pp, factors = "Lex")
# This plot is a bit odd since frequency is only non-zero for words
# so quantile 2 contains all the non-words
plot_cdf(dat, pp, factors = "Freq")

