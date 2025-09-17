RNGkind("L'Ecuyer-CMRG")
set.seed(123)

lIfun <- function(d) factor(rep(2,nrow(d)),levels=1:2)

# All staircase (indicated by NA)
mySSD_function <- function(d)SSD_function(d,SSD=NA,p=.25)

designSSexG <- design(model=SSEXG,
                      factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
                      matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
                      functions=list(lI=lIfun),
                      covariates = "SSD",
                      formula=list(m~lM,sigma~1,k~1,mS~1,sigmaS~1,kS~1, gf~1,tf~1)
)

p_vector <- sampled_pars(designSSexG,doMap = FALSE)
p_values <- c(
  m = log(0.65),
  sigma = log(0.8),
  k = qnorm(0.07),
  mS = log(0.23),
  sigmaS = log(0.2),
  kS = qnorm(0.13),
  gf = qnorm(.1),
  tf = qnorm(.1)
)
p_vector[names(p_values)] <- p_values


dat <- make_data(p_vector, designSSexG, n_trials = 10,staircase=TRUE, functions = list(SSD = mySSD_function))
emc <- make_emc(dat, designSSexG, type = "single")


test_that("exG", {
  expect_snapshot(
    init_chains(emc, particles = 10, cores_for_chains = 1)[[1]]$samples
  )
})
