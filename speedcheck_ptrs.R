library(EMC2)
rm(list=ls())
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR
trend <- make_trend(par_names = "m", cov_names = list(c("covariate1")),
                    kernels = "delta", ffill_na=TRUE)
design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c('covariate1'),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)
##mapped_pars(design_base)
p_vector <- sampled_pars(design_base, doMap = FALSE)

p_vector[1:7] <- c(-1, 1.5, log(1), log(.2), 1, .5, qnorm(.2))

covariate1 <- c(NA, 1, NA, NA, 1, NA, NA, NA, NA, NA)#, 1, NA, 1, 1, NA, rep(NA, 10))

#debug(EMC2:::run_kernel)
dat <- make_data(p_vector, design_base, n_trials = 5, covariates = data.frame(covariate1 = covariate1), return_trialwise_parameters=TRUE)
# cbind(attr(dat, 'trialwise_parameters'), rep(covariate1, each=2))


## Test code for comparing the output of run_trend between R and Rcpp
# signature:
# NumericVector run_trend_rcpp(DataFrame data, List trend, NumericVector param, NumericMatrix trend_pars, NumericMatrix pars_full) {
emc <- make_emc(dat, design_base, type='single')
dadm <- emc[[1]]$data[[1]]
trend <- emc[[1]]$model()$trend$m
param <- rep(0,  nrow(dadm))
trend_pars <- matrix(rep(c(p_vector[5:7]), each=nrow(dadm)), ncol=3, byrow=FALSE)
trend_pars[,3] <- pnorm(trend_pars[,3])
pars_full <- matrix(rep(p_vector, each=nrow(dadm)), ncol=length(p_vector), byrow=FALSE)
pars_full[,ncol(pars_full)] <- pnorm(pars_full[,ncol(pars_full)])
cv1_updated <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
                                     pars_full = pars_full, return_kernel = TRUE, use_ptrs = FALSE)
cv1_updated

cv1_updated_ptr <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
                                         pars_full = pars_full, return_kernel = TRUE, use_ptrs = TRUE)

cv1_updated_noptr <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
                                           pars_full = pars_full, return_kernel = TRUE, use_ptrs = FALSE)


n_trials <- 10000
covariate1 <- rnorm(n_trials)
covariate2 <- rnorm(n_trials)
covariate3 <- rnorm(n_trials)
covariate2[sample(1:n_trials, 50)] <- NA

trend <- make_trend(par_names = "m", cov_names = list(c("covariate1","covariate2","covariate3")),
                    kernels = "delta", ffill_na=TRUE)
design_base <- design(factors = list(subjects = 1, S = 1:2),
                      Rlevels = 1:2,
                      covariates = c("covariate1","covariate2","covariate3"),
                      matchfun = matchfun,
                      trend = trend,
                      formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                      contrasts = list(lM = ADmat),
                      model = LNR)
dat <- make_data(p_vector, design_base, n_trials = n_trials,
                 covariates = data.frame(covariate1 = covariate1,
                                         covariate2 = covariate2,
                                         covariate3 = covariate3
                 ), return_trialwise_parameters=TRUE)
emc <- make_emc(dat, design_base, type='single')


emc <- make_emc(dat, design_base, type='single')
dadm <- emc[[1]]$data[[1]]
trend <- emc[[1]]$model()$trend$m
param <- rep(0,  nrow(dadm))
trend_pars <- matrix(rep(c(p_vector[5:7]), each=nrow(dadm)), ncol=3, byrow=FALSE)
trend_pars[,3] <- pnorm(trend_pars[,3])
pars_full <- matrix(rep(p_vector, each=nrow(dadm)), ncol=length(p_vector), byrow=FALSE)
pars_full[,ncol(pars_full)] <- pnorm(pars_full[,ncol(pars_full)])

cv1_updated_ptr <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
                                         pars_full = pars_full, return_kernel = TRUE, use_ptrs = TRUE)

cv1_updated_noptr <- EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
                                           pars_full = pars_full, return_kernel = TRUE, use_ptrs = FALSE)

microbenchmark::microbenchmark(noptrs=EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
                                                          pars_full = pars_full, return_kernel = TRUE, use_ptrs = FALSE),
                               ptrs=EMC2:::run_trend_rcpp(data = dadm, trend=trend, param=param, trend_pars=trend_pars,
                                                     pars_full = pars_full, return_kernel = TRUE, use_ptrs = TRUE)
                               )
