rm(list=ls())
install.packages("EMC2")

# test_in_package <- F
#
# if(test_in_package){
#   library(devtools)
#   load_all()
# } else{
#   library(Rcpp)
#   all_files <- list.files("R")
#   all_files <- all_files[!all_files %in% "EMC2-package.R"]
#   for(file in all_files) source(paste0("R/", file))
#   sourceCpp("~/Documents/UVA/2022/test_rcpp/test_Niek.cpp")
# }


load("~/Documents/UVA/2022/EMC_test/PNAS.RData")

dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Average rate = intercept, and rate d = difference (match-mismatch) contrast
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
ADmat

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat

design_RDM <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,s~1,B~E,A~1,t0~1),
  constants=c(s=log(1)),
  model=rdmB)


# Test single subject
dat_single <- dat[which(dat$subjects %in% (unique(dat$subjects)[1])),]
dat_single <- droplevels(dat_single)

samplers <- make_samplers(dat_single, design_B, type = "single", n_chains = 3)
samplers <- auto_burn(samplers, verbose = T, min_es = 1000)
# samplers <- run_samplers(samplers, stage = "preburn", iter = 50, verbose = T, cores_per_chain = 1, cores_for_chains = 1)

samplersC_merg <- merge_samples(samplers)
samplersR_merg <- merge_samples(samplers)
# microbenchmark::microbenchmark(
#   samplers <- run_samplers(samplers, stage = "preburn", iter = 50, verbose = T, cores_per_chain = 1, cores_for_chains = 1),
#   times = 15
# )


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
