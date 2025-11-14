RNGkind("L'Ecuyer-CMRG")
set.seed(123)

lIfun <- function(d) factor(rep(2,nrow(d)),levels=1:2)

# All staircase (indicated by NA)
mySSD_function <- make_ssd()

designSSexG <- design(model=SSEXG,
                      factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
                      matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
                      functions=list(lI=lIfun),
                      covariates = "SSD",
                      formula=list(mu~lM,sigma~1,tau~1,muS~1,sigmaS~1,tauS~1, gf~1,tf~1)
)

p_vector <- sampled_pars(designSSexG,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(.6), log(.8), log(0.05), log(0.2),
                                  log(0.2),log(0.03), log(0.05),
                                  qnorm(.1),qnorm(.1))

dat <- make_data(p_vector, designSSexG, n_trials = 10,
                 functions = list(SSD = mySSD_function))
emc <- make_emc(dat, designSSexG, type = "single")


test_that("exG", {
  expect_snapshot(
    init_chains(emc, particles = 10, cores_for_chains = 1)[[1]]$samples
  )
})

test_that("staircase resets per subject", {
  design_multi <- design(model = SSEXG,
                         factors = list(subjects = 1:3, S = c("left", "right")),
                         Rlevels = c("left", "right"),
                         matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
                         functions = list(lI = lIfun),
                         formula = list(mu ~ lM, sigma ~ 1, tau ~ 1,
                                         muS ~ 1, sigmaS ~ 1, tauS ~ 1,
                                         gf ~ 1, tf ~ 1))

  p_multi <- sampled_pars(design_multi, doMap = FALSE)
  p_multi[1:length(p_multi)] <- c(log(.6), log(.8), log(0.05), log(0.2),
                                  log(0.2), log(0.03), log(0.05),
                                  qnorm(.1), qnorm(.1))

  dat_multi <- make_data(p_multi, design_multi, n_trials = 400,
                         functions = list(SSD = make_ssd()))

  finite_ssd <- is.finite(dat_multi$SSD)
  if (!all(tapply(finite_ssd, dat_multi$subjects, any))) {
    skip("No stop trials generated for at least one subject")
  }
  first_ssd <- tapply(dat_multi$SSD[finite_ssd], dat_multi$subjects[finite_ssd], `[`, 1)
  expect_equal(as.numeric(first_ssd), rep(0.25, length(first_ssd)))
})
