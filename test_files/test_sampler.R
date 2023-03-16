rm(list=ls())
library(devtools)
load_all()

load("~/Documents/UVA/2022/EMC_test/PNAS.RData")

dat <- data[,c("s","E","S","R","RT")]
dat$response <- dat$R
dat$response[dat$response == "right"] <- "upper"
dat$response[dat$response == "left"] <- 'lower'
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Average rate = intercept, and rate d = difference (match-mismatch) contrast
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
ADmat

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat

ll_func <- function(pars, dadm){
  out <- rtdists::ddiffusion(dadm, a = exp(pars["a"]), t0 = exp(pars["t0"]),
                      v = pars["v"])
  bad <- (out<1e-10)|(!is.finite(out))
  out[bad] <- 1e-10
  return(sum(log(pmax(out, 1e-10))))
}

design_B <- make_design(custom_p_vector = c("a", "t0", "v"), model=ll_func)


# Test single subject
dat_single <- dat[which(dat$subjects %in% (unique(dat$subjects)[1:4])),]
dat_single <- droplevels(dat_single)

samplers <- make_samplers(dat_single, design_B, type = "standard", n_chains = 3)
samplers <- auto_burn(samplers, verbose = T, cores_per_chain = 4, cores_for_chains = 3)

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

