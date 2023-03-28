rm(list = ls())
library(EMC2)
library(formula.tools)

load("test_files/Juanita.RData")

sub_dat <- sub_dat[sub_dat$C != "miss",]
dat <- sub_dat[,c("s", "Target", "SLength", "SType", "R", "RT", "group")]
names(dat) <- c("subjects", "target", "length", "type", "R", "rt", "group")
dat$target[dat$target == "LeftResponse"] <- "left"
dat$target[dat$target == "RightResponse"] <- "right"
dat$R[dat$R %in% c(211, 201)] <- "left"
dat$R[dat$R %in% c(212, 202)] <- "right"

dat$type[!dat$type %in% c("EHD", "LD")] <- "STD"
dat <- dat[dat$type != "EHD",]

dat$rt <- dat$rt/1000

dat$R <- factor(dat$R)
dat$subjects <- factor(dat$subjects)
dat$target <- factor(dat$target)
dat$group <- factor(dat$group)
dat$type <- factor(dat$type, levels = c("STD", "LD"))


design_lbaB <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),target=levels(dat$target), type = levels(dat$type)),
  Rlevels=levels(dat$R),matchfun=function(d)d$target==d$lR,
  Flist=list(v~lM*type,sv~1,B~1,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)

samplers_inter_group <- make_samplers(dat, design_lbaB, type = "lm",
  formula = list(v_lMTRUE ~ group, v_typeLD ~ group, v_lMTRUE:typeLD ~ group))
samplers_inter_group <- run_emc(samplers_inter_group, cores_per_chain = 5, verbose = T)
save(samplers_inter_group, file = "inter_group.RData")

samplers_main_group <- make_samplers(dat, design_lbaB, type = "lm",
                                     formula = list(v_lMTRUE ~ group, v_typeLD ~ group))
samplers_main_group <- run_emc(samplers_main_group, cores_per_chain = 5)
save(samplers_main_group, file = "main_group.RData")

samplers_mainB_group <- make_samplers(dat, design_lbaB, type = "lm", n_chains = 3,
                                      formula = list(v_lMTRUE ~ group, v_typeLD ~ group, B ~ group))
samplers_mainB_group <- run_emc(samplers_mainB_group, cores_per_chain = 5)
save(samplers_mainB_group, file = "main_Bgroup.RData")

samplers_interB_group <- make_samplers(dat, design_lbaB, type = "lm", n_chains = 3,
                                       formula = list(v_lMTRUE ~ group, v_typeLD ~ group, v_lMTRUE:typeLD ~ group, B ~ group))
samplers_interB_group <- run_emc(samplers_interB_group, cores_per_chain = 5)
save(samplers_interB_group, file = "inter_Bgroup.RData")

samplers <- samplers_interB_group
theta <- get_theta(samplers[[1]], mapped = T)

rowMeans(summ_old)
rowMeans(summ_young)
rowMeans(summ_SZ)

rowMeans(theta$B)
rowMeans(theta$v_lMTRUE)
rowMeans(theta$`v_lMTRUE:typeLD`)


samplers_main_group <- run_IS2(samplers_main_group, n_cores = 15, max_particles = 2500, IS_samples = 500)
save(samplers_main_group, file = "main_group.RData")

samplers_interB_group <- run_IS2(samplers_interB_group, n_cores = 15, max_particles = 2500, IS_samples = 500)
save(samplers_interB_group, file = "inter_Bgroup.RData")

samplers_mainB_group <- run_IS2(samplers_mainB_group, n_cores = 15, max_particles = 2500, IS_samples = 500)
save(samplers_mainB_group, file = "main_Bgroup.RData")

samplers_1mainB_group <- make_samplers(dat, design_lbaB, type = "lm", n_chains = 3,
                                       formula = list(v_lMTRUE ~ group, B ~ group))
samplers_1mainB_group <- run_emc(samplers_1mainB_group, cores_per_chain = 5, useC = F, verbose = T)
save(samplers_1mainB_group, file = "main_Bgroup.RData")

samplers_mainB_group <- run_IS2(samplers_mainB_group, n_cores = 15, max_particles = 2500, IS_samples = 500)
save(samplers_mainB_group, file = "main_Bgroup.RData")

samplers_1mainB_group <- run_IS2(samplers_1mainB_group, IS_samples = 500, max_particles = 2500, n_cores = 15)
save(samplers_1mainB_group, file = "1main_Bgroup.RData")


pp <- post_predict(samplers)
plot_fit(dat, pp, signalFactor = "target", factors = c("group", "target"), layout = c(2,3))

