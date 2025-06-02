RNGkind("L'Ecuyer-CMRG")
set.seed(123)

lIfun <- function(d) factor(rep(2,nrow(d)),levels=1:2)

# All staircase (indicated by NA)
mySSD_function <- function(d)SSD_function(d,SSD=NA,p=.25)

designSSexG <- design(model=SSexG,
                      factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
                      matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
                      functions=list(lI=lIfun, SSD = mySSD_function),
                      formula=list(mu~lM,sigma~1,tau~1,muS~1,sigmaS~1,tauS~1, gf~1,tf~1)
)

p_vector <- sampled_pars(designSSexG,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(.6), log(.8), log(0.05), log(0.2),
                                  log(0.2),log(0.03), log(0.05),
                                  qnorm(.1),qnorm(.1))


dat <- make_data(p_vector, designSSexG, n_trials = 10,staircase=TRUE)
emc <- make_emc(dat, designSSexG, type = "single")


test_that("exG", {
  expect_snapshot(
    init_chains(emc, particles = 10, cores_for_chains = 1)[[1]]$samples
  )
})
