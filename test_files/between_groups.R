rm(list = ls())

devtools::load_all()
load("~/Documents/UVA/2023/BetweenSubsApplications/Manning2022/Manning2022.RData")

data <- data[as.numeric(data$subjects) %% 2 != 0,]

source("test_files/sampling_lm.R")
source("test_files/helps_lm.R")
source("test_files/utils_lm.R")

# Remove NULL responses
data <- data[data$keypress != 0,]

data$keypress[data$keypress == 1] <- "left"
data$keypress[data$keypress == 2] <- "right"

colnames(data)[c(5,6,7)] <- c("R", "rt", "subjects")
data$R <- factor(data$R)
data$subjects <- factor(data$subjects)
data$cond <- factor(data$cond)
data$group <- factor(data$group)


aggregate(EEG_respslope ~ cond, data = data, FUN = mean)


dat_coh <- data[data$task == "coherence",]
dat_coh$slope_mean <- NA
dat_coh$slope_diff <- NA

for(sub in unique(dat_coh$subjects)){
  idx <- dat_coh$subjects == sub
  mean_sub <- mean(dat_coh$EEG_respslope[idx])
  diff_sub <- 2*(dat_coh$EEG_respslope[dat_coh$cond == "easy"][idx][1] - mean_sub)
  dat_coh$slope_mean[idx] <- mean_sub
  dat_coh$slope_diff[idx] <- diff_sub
}
form <- list(v~cond + (cond|subjects), a~1+ (1|subjects), t0 ~ 1 + (1|subjects),
             Z ~ 1 + (1|subjects), s ~1, DP ~ 1, SZ ~ 1, sv ~1, st0 ~ 1)

# NAs in slope mean are treated as 0!!
design <- make_design_lm(form, ddata = dat_coh, model = DDM,
                       constants = c(sv = log(0), s = log(1), st0 = log(0), DP = qnorm(.5),
                                     SZ = qnorm(0), st0 = log(0)))

samplers <- make_samplers_lm(dat_coh, design)
# Things to do:
# Check whether replacing NAs with 0 is alright for DMs
# Check the adding of between variables out of the random variables in get_conditionals (happens with nesting)!!!!
# Finding which are intercepts with 0 + in pmwgs$is_intercept

source("test_files/sampling_lm.R")

create_eff_proposals_between <- function(samples){
  idx_subtract <- min(300, samples$idx/2)
  idx <- round(samples$idx - idx_subtract):samples$idx
  all_samples <- rbind(samples$fixed[,idx], samples$random[,idx])
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- stats::cov(t(all_samples))
  n_pars <- nrow(samples$fixed)
  if(any(eigen(var_tilde, only.values = T)$values < 1e-08)){
    return(list(eff_mu = mu_tilde[1:n_pars], eff_var = as.matrix(nearPD(var_tilde[1:n_pars, 1:n_pars])$mat)))
  }
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$random[,samples$idx]))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

debug(create_eff_proposals_between)
samplers <- run_samplers(samplers, stage = "burn", verbose = T, cores_per_chain = 10, cores_for_chains = 1, fileName = "test_groups.RData", mean_gd = 1.1)

source("test_files/sampling_lm.R")
source("test_files/helps_lm.R")
source("test_files/utils_lm.R")

samplers <- run_adapt(samplers, verbose = T, cores_per_chain = 4, cores_for_chains = 3, particle_factor = 30)
samplers <- run_sample(samplers, verbose = T, cores_per_chain = 4, cores_for_chains = 3, iter = 700)

plot_chains(samplers, selection = "random", filter = "sample", subfilter = 50)


devtools::load_all()
source("test_files/sampling_lm.R")
source("test_files/helps_lm.R")
source("test_files/utils_lm.R")

samplers <- run_sample(samplers, cores_per_chain = 4, cores_for_chains = 3, verbose = T)

samplers_vEEG_main <- make_samplers(dat_coh, design, type = "lm",
                                    formula = list(v ~ slope_mean))

samplers_vEEG_inter <- make_samplers(dat_coh, design, type = "lm",
                                     formula = list(v ~ slope_mean*group))

samplers_vdiffEEG_main <- make_samplers(dat_coh, design, type = "lm",
                                        formula = list(v_diff ~ slope_diff))

samplers_vdiffEEG_inter <- make_samplers(dat_coh, design, type = "lm",
                                         formula = list(v_diff ~ slope_diff*group))



samplers_a <- make_samplers(dat_coh, design, type = "lm",
                            formula = list(a ~ group))
samplers_v <- make_samplers(dat_coh, design, type = "lm",
                            formula = list(v ~ group))
samplers_v_diff <- make_samplers(dat_coh, design, type = "lm",
                                 formula = list(v_cond1 ~ group))
samplers_t0 <- make_samplers(dat_coh, design, type = "lm",
                             formula = list(t0 ~ group))
samplers_Z <- make_samplers(dat_coh, design, type = "lm",
                            formula = list(Z ~ group))

debug(make_samplers)
devtools::load_all("~/Documents/UVA/2023/EMC2_Andrew/")
samplers_null <- make_samplers(dat_coh, design, type = "lm")

samplers_null <- auto_burn(samplers_null, cores_for_chains = 1, verbose = T, cores_per_chain = 15)
samplers_null <- run_adapt(samplers_null, cores_for_chains = 1, cores_per_chain = 15, min_unique = 25,
                           step_size = 25, verbose = T)

samplers_null <- run_IS2(samplers_null, filter = "burn", IS_samples = 5)


samplers_a <- run_emc(samplers_a, cores_per_chain = 8, verbose = T, fileName = "a.RData")
samplers_v <- run_emc(samplers_v, cores_per_chain = 8, verbose = T, fileName = "v.RData")
samplers_v_diff <- run_emc(samplers_v, cores_per_chain = 8, verbose = T, fileName = "v_diff.RData")
samplers_t0 <- run_emc(samplers_t0, cores_per_chain = 8, verbose = T, fileName = "t0.RData")
samplers_Z <- run_emc(samplers_Z, cores_per_chain = 8, verbose = T, fileName = "Z.RData")

samplers_a <- run_IS2(samplers_a, n_cores = 25)
samplers_v <- run_IS2(samplers_v, n_cores = 25)
samplers_v_diff <- run_IS2(samplers_v_diff, n_cores = 25)
samplers_t0 <- run_IS2(samplers_t0, n_cores = 25)
samplers_Z <- run_IS2(samplers_Z, n_cores = 25)




