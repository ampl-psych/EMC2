dGNG <- design(Rlevels = c("left","right"),
               factors=list(subjects=1,S=c("left","right")),
               functions=list(
                 TIMEOUT=function(d)rep(2.5,nrow(d)),
                 Rnogo=function(d)factor(rep("left",nrow(d)),levels=c("left","right")), # no go response level
                 Rgo=function(d)factor(rep("right",nrow(d)),levels=c("left","right"))), # go response level
               formula=list(v~S,a~1, Z~1, t0~1),
               model=DDMGNG)

p_vector <- sampled_pars(dGNG,doMap=F)
p_vector[1:5] <- c(0,1,log(1),qnorm(.5),log(.4))

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

sdat <- make_data(p_vector,dGNG,n_trials=10)
emc <- make_emc(sdat,dGNG,type="single", compress = FALSE, n_chains = 1)

test_that("DDMGNG", {
  expect_snapshot(init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples)
  expect_snapshot(make_data(p_vector,dGNG,n_trials=10))
})


dprobit <- design(Rlevels = c("left","right"),
                  factors=list(subjects=1,S=c("left","right")),
                  formula=list(mean ~ 0+S, sd ~ 1,threshold ~ 1),
                  matchfun=function(d)d$S==d$lR,
                  constants=c(sd=log(1),threshold=0),
                  model=SDT)

p_vector <- sampled_pars(dprobit,doMap=F)
p_vector[1:2] <- c(-1,1)
sdat <- make_data(p_vector,dprobit,n_trials=10)
emc <- make_emc(sdat,dprobit,type="single", compress = F, n_chains = 1)

test_that("probit", {
  expect_snapshot(init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples)
  expect_snapshot(make_data(p_vector,dprobit,n_trials=10))
})
