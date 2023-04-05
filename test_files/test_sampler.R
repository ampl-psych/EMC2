rm(list = ls())
library(devtools)
# Rcpp only works as fast with install like this unfortunately, if this throws errors, restarting Rstudio worked for me
devtools::install("~/Documents/UVA/2022/EMC2")
load_all()
# library(EMC2)
library(statmod)
library(corrplot)
library(factor.switching)
library(colorspace)

n_pars <- 6
n_subjects <- 20
n_trials <- 100
qs <- n_pars:1
qs[as.logical(qs %% 2)] <- -qs[as.logical(qs %% 2)]
covmat <- toeplitz(qs/n_pars)

is_zero <- sample(1:sum(lower.tri(covmat)), n_pars * n_pars /5)
# covmat[lower.tri(covmat)][is_zero] <- rnorm(n_pars * n_pars /5, sd = .2)
covmat[lower.tri(covmat)] <- covmat[lower.tri(covmat)]/1.5
covmat[upper.tri(covmat)] <- t(covmat)[upper.tri(covmat)]
corrplot(covmat)

data1 <- mvtnorm::rmvnorm(n_subjects, mean = seq(-5, 5, length.out = n_pars), sigma =  covmat)
new_df1 <- data.frame()
cov_person <- MCMCpack::riwish(1.5*n_pars, diag(n_pars))
for(i in 1:nrow(data1)){
  tmp <- data.frame(mvtnorm::rmvnorm(n_trials, data1[i,], cov_person))
  tmp$subjects <- i
  new_df1 <- rbind(new_df1, tmp)
}

data2 <- mvtnorm::rmvnorm(n_subjects, mean = seq(-5, 5, length.out = n_pars), sigma =  covmat)
new_df2 <- data.frame()
for(i in 1:nrow(data2)){
  tmp <- data.frame(mvtnorm::rmvnorm(n_trials, data2[i,], cov_person))
  tmp$subjects <- i
  new_df2 <- rbind(new_df2, tmp)
}

custom_ll <- function(dadm, pars){
  sum(dnorm(t(as.matrix(dadm[,1:length(pars)])), mean = pars, sd = 1, log = T))
}

design_normal <- make_design(model = custom_ll, custom_p_vector = paste0("X", 1:n_pars))

# first speed test:
samplers <- make_samplers(list(new_df1, new_df2), list(design_normal, design_normal))
samplers <- run_samplers(samplers, stage = "preburn", iter = 100, cores_per_chain =10, cores_for_chains =1, verbose = T, verboseProgress = T,
                         n_blocks = 2)
samplers <- run_samplers(samplers, stage = "burn", max_gd = 1.2, cores_per_chain =4, cores_for_chains =3, verbose = T, verboseProgress = T,
                         n_blocks = 2)
samplers <- run_samplers(samplers, stage = "adapt", max_gd = 1.2, cores_per_chain =4, cores_for_chains =3, verbose = T, verboseProgress = T,
                         n_blocks = 2)
samplers <- run_samplers(samplers, stage = "sample", max_gd = 1.2, cores_per_chain =4, cores_for_chains =3, verbose = T, verboseProgress = T,
                         n_blocks = 2)
samplers <- run_emc(samplers, cores_per_chain =1, cores_for_chains = 1, verbose = T)


samplers <- run_samplers(samplers, stage = "preburn", iter = 200, cores_per_chain =10, cores_for_chains =1, useC = F, verbose = T, verboseProgress = T)


samplers <- auto_burn(samplers, cores_per_chain =10, cores_for_chains =1, useC = F, verbose = T, verboseProgress = T)
adapted <- run_adapt(samplers, cores_per_chain =3, cores_for_chains = 3, useC = F, verbose = T, min_unique = 50)
debug(create_eff_proposals)
debug(new_particle)
sampled <- run_sample(adapted, cores_per_chain = 1, cores_for_chains = 3, useC = F, verbose = T, verboseProgress = T)

samplers_infnt <- make_samplers(new_df, design_normal, type = "infnt_factor")
samplers_infnt <- run_emc(samplers_infnt, cores_per_chain =4, cores_for_chains = 3, useC = F, verbose = T)

samplers_diag <- make_samplers(new_df, design_normal, type = "diagonal")
samplers_diag <- run_sample(samplers_diag, cores_per_chain =4, cores_for_chains = 3, useC = F, verbose = T)



# first let's run init separately

library(factor.switching)
library(EMC2)
library(corrplot)
samples_fa <- merge_samples(samplers_infnt)
samples_std <- merge_samples(samplers)
samples_diag <- merge_samples(samplers_diag)

alpha_diag <- samples_diag$samples$alpha[,,samples_diag$samples$stage == "sample"]
cor_diag <- array(0, dim = c(samples_diag$n_pars, samples_diag$n_pars, dim(alpha_diag)[3]))
for(i in 1:dim(alpha_diag)[3]){
  cor_diag[,,i] <- cor(t(alpha_diag[,,i]))
}
corrs_with_numbers(samples_std, do_corr = T)
corrs_with_numbers(samples_fa, do_corr = T)
corrs_with_numbers(corrs = cor_diag, do_corr = T)

corrplot(covmat[1:10, 1:10], addCoef.col ="black", number.cex = .75, tl.col = "black")

loadings <- samples_fa$samples$theta_lambda[,1:5,samples_fa$samples$stage == "sample"]

samples_fa_sml <- samples_fa
samples_fa_sml$samples$theta_lambda <- samples_fa_sml$samples$theta_lambda[,1:2,]
for(i in 1:dim(samples_fa_sml$samples$theta_var)[3]){
  samples_fa_sml$samples$theta_var[,,i] <- samples_fa_sml$samples$theta_lambda[,,i] %*% t(samples_fa_sml$samples$theta_lambda[,,i]) + diag(1/samples_fa$samples$theta_sig_err_inv[,i])
}

corrs_with_numbers(samples_fa_sml)

max_factors <- 7

# align loadings ----------------------------------------------------------
samples <- samples_fa
idx <- which(samples$samples$stage == "sample")
idx <- idx[-c(1:200)]

loadings_recov <- aperm(samples$samples$theta_lambda[,1:max_factors,idx], c(2,1,3))
lambda_mcmc <- t(matrix(loadings_recov, prod(dim(loadings_recov)[1:2]), dim(loadings_recov)[3]))
expanded <- expand.grid(1:max_factors, 1:samples$n_pars)
colnames(lambda_mcmc) <- paste0("LambdaV", expanded[,2], "_", expanded[,1])


library(factor.switching)
res <- rsp_exact(lambda_mcmc)

plot(res)
probs <- .95
factors_out <- c(which(credible.region(res$lambda_reordered_mcmc, probs = probs)[[1]][1,] > 0),
                 which(credible.region(res$lambda_reordered_mcmc, probs = probs)[[1]][2,] < 0))
factors_out <- unique(as.numeric(sub("^[^_]*_", "", names(factors_out))))
factors_out

lambda_fixed <- array(0, dim = c(nrow(loadings), max_factors, length(idx)))
for(i in 1:length(idx)){lambda_fixed[,,i] <- t(matrix(t(res$lambda_reordered_mcmc)[,i], nrow = max_factors, byrow = F))}

samples_fixed <- samples_fa
samples_fixed$samples$theta_lambda[,,samples_fixed$samples$stage == "sample"] <- lambda_fixed

debug(corrs_with_numbers)
corrs_with_numbers(loadings = lambda_fixed)
#
# IS_samples <- attr(samplers, "IS_samples")
# group_sample1 <- IS_samples$sub_and_group
#
# log_ratios <- -1 * group_sample1
# r_eff <- relative_eff(exp(-log_ratios))
# psis_result <- psis(log_ratios[1,,], r_eff = r_eff)
#
# group_sample <- group_sample1[,1:500,]
# group_sample1[group_sample1 < -5000] <- -5000
# group_sample1[group_sample1 > 5000] <- 5000

library(loo)

# Note that the useC argument is a bit deceptive. It uses the C dists in any case, but if TRUE will also use the C log_likelihood_race
microbenchmark::microbenchmark(
  run_IS2(samplers, filter = "burn", IS_samples = 10, max_particles = 1000, n_cores = 10),
  run_IS2(samplers, filter = "burn", IS_samples = 10, max_particles = 1000, n_cores = 10, useC = F), times = 20
) # Literally twice as fast!

# Then compare
samplers <- make_samplers(dat_single, design_RDM, type = "single")
samplersC <- auto_burn(samplers, verbose = T, min_es = 1000, useC = T)
samplersR <- auto_burn(samplers, verbose = T, min_es = 1000, useC = F)
# samplers <- run_samplers(samplers, stage = "preburn", iter = 50, verbose = T, cores_per_chain = 1, cores_for_chains = 1, useC = T)
# samplers <- run_samplers(samplers, stage = "burn", iter = 50, verbose = T, cores_per_chain = 1, cores_for_chains = 1, useC = T)


samplersC_merg <- merge_samples(samplersC)
samplersR_merg <- merge_samples(samplersR)



library(ggplot2)
library(cowplot)

N <- 1500
plots <- list()
for(i in 1:nrow(samplersR_merg$samples$alpha)){
  R <- samplersR_merg$samples$alpha[i,1,( samplersR_merg$samples$idx - N):  samplersR_merg$samples$idx]
  C <- samplersC_merg$samples$alpha[i,1,( samplersC_merg$samples$idx - N):  samplersC_merg$samples$idx]
  df <- data.frame(par_val = c(R, C), group = rep(c("R", "C"), each = N + 1))
  ggp <- ggplot(df, aes(x = par_val, colour = group, fill = group)) +
    geom_density(alpha = .1) + ggtitle(samplersR_merg$par_names[i]) + theme_bw()
  plots[[i]] <- ggp
}
plot_grid(plots[[1]], plots[[2]])
plot_grid(plots[[3]], plots[[4]])
plot_grid(plots[[5]], plots[[6]])
plot_grid(plots[[7]])



