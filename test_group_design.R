rm(list = ls())
# devtools::install()
# library(EMC2)
devtools::load_all("~/Documents/2025/EMC2/")
load("~/Documents/2025/EMC2_tutorial/Rfiles/data_ELP.RData")

data <- data[data$subjects %in% unique(data$subjects)[seq(1, 200, by = 10)],]
data$uni <- factor(data$uni)
data$subjects <- factor(data$subjects)
table(data$uni)


data$log_freq[is.na(data$log_freq)] <- 0
data$R <- factor(data$R)
data$S <- factor(data$S)
Smat <- cbind(d = c(-1, 1))

data <- data[,!colnames(data) %in% c("item", "length")]

# Full analysis -----------------------------------------------------------
full_design <- design(list(v ~ S*log_freq, a ~ 1, t0 ~ 1, Z ~ 1, sv ~ 1),
                      data = data, model = DDM,
                      contrasts = list(S = Smat),
                      constants = c(v_log_freq = 0))

sampled_pars(full_design)

data$age <- data$age - mean(data$age)
group_des <- group_design(formula = list(v_Sd ~ age, `v_Sd:log_freq` ~ age,
                                         t0 ~ age + (1|uni), a ~ age + (1|uni)),
                          data = data, subject_design = full_design)
pri <- prior(design = full_design, group_design = group_des)
plot(pri, selection = "s")
#
emc <- make_emc(data, design = full_design, group_design = group_des)
#
#
#
debug(get_startpoints_standard)
emc <- init_chains(emc, cores_for_chains = 1, cores_per_chain = 3)

emc <- fit(emc, cores_for_chains = 1, cores_per_chain = 5)
