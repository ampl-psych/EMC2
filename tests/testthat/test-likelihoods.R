# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR

dat <- expand.grid(subjects = factor(1:2), S = factor(c("left", "right")),
                   E = factor(c("speed", "neutral", "accuracy"),
                              levels = c("speed", "neutral", "accuracy")),
                   R = factor(c("left", "right")))
dat$rt <- rep(c(.35, .45, .6, .8), length.out = nrow(dat))

# First race models
design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                            formula=list(m~lM + E,s~1,t0~1),
                            contrasts=list(m=list(lM=ADmat)))

p_LNR <- c(m=-1,m_lMd=2, m_Eneutral = .1, m_Eaccuracy = .1, s = log(.5), t0=log(.2))

design_LBA <- design(data = dat,model=LBA,matchfun=matchfun,
                            formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
                            contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))

p_LBA <- c(v=2,v_lMd=1, sv_lMTRUE = log(.5), B = log(.4), B_Eneutral = .1,
           B_Eaccuracy = .2, B_lRright = 0, A = log(.2), t0=log(.2))


design_RDM <- design(data = dat,model=RDM,matchfun=matchfun,
                            formula=list(v~lM,s~lM,B~E+lR,A~1,t0~1),
                            contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))

p_RDM <- c(v=log(2),v_lMd=1, s_lMTRUE = log(.8), B = log(.4), B_Eneutral = .1,
           B_Eaccuracy = .2, B_lRright = 0, A = log(.2), t0=log(.2))

# Also DDM
design_DDM <- design(data = dat,model=DDM,
                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1, st0~1),
                            constants=c(s=log(1)))

p_DDM <- c(v_Sleft = -2, v_Sright = 2, a = log(2), a_Eneutral = .3, a_Eaccuracy = .3, t0 = log(.2),
           Z = qnorm(.5), sv = log(.5), SZ = qnorm(.3), st0 = log(.05))

LNR_s <- make_emc(dat, design_LNR, type = "single", compress = FALSE, n_chains = 1)
LBA_s <- make_emc(dat, design_LBA, type = "single", compress = FALSE, n_chains = 1)
RDM_s <- make_emc(dat, design_RDM, type = "single", compress = FALSE, n_chains = 1)
DDM_s <- make_emc(dat, design_DDM, type = "single", compress = FALSE, n_chains = 1)

test_that("LNR likelihood", {
  parameters <- likelihood_parameters(p_LNR)
  parameters <- rbind(parameters, invalid = p_LNR)
  parameters[nrow(parameters), "t0"] <- log(10)
  ll <- model_likelihood(LNR_s, parameters)
  expected <- apply(parameters, 1, function(p) {
    effect <- c(speed = 0, neutral = p[["m_Eneutral"]], accuracy = p[["m_Eaccuracy"]])
    m <- p[["m"]] + effect[as.character(dat$E)]
    difference <- ifelse(dat$S == dat$R, .5, -.5) * p[["m_lMd"]]
    rt <- dat$rt - exp(p[["t0"]])
    trial_ll <- dlnorm(rt, m + difference, exp(p[["s"]]), log = TRUE) +
      plnorm(rt, m - difference, exp(p[["s"]]), lower.tail = FALSE, log.p = TRUE)
    tapply(pmax(log(1e-10), trial_ll), dat$subjects, sum)
  })
  expect_equal(unname(ll), unname(t(expected)), tolerance = 1e-8)
  # Repeated observations exercise the compressed likelihood's expansion.
  compressed <- make_emc(rbind(dat, dat), design_LNR,
                         type = "single", compress = TRUE, n_chains = 1)
  expect_equal(model_likelihood(compressed, parameters), 2 * ll)
})

test_that("LBA likelihood", {
  parameters <- likelihood_parameters(p_LBA)
  parameters <- rbind(parameters, invalid = p_LBA)
  parameters[nrow(parameters), "t0"] <- log(10)
  expect_snapshot(round(model_likelihood(LBA_s, parameters), 6))
})

test_that("RDM likelihood", {
  parameters <- likelihood_parameters(p_RDM)
  # A = 0 has a separate density/CDF implementation.
  parameters[1:8, "A"] <- -Inf
  parameters <- rbind(parameters, invalid = p_RDM)
  parameters[nrow(parameters), "t0"] <- log(10)
  expect_snapshot(round(model_likelihood(RDM_s, parameters), 6))
})

test_that("DDM likelihood", {
  parameters <- likelihood_parameters(p_DDM)
  parameters <- rbind(parameters, invalid = p_DDM)
  parameters[nrow(parameters), "t0"] <- log(10)
  expect_snapshot(round(model_likelihood(DDM_s, parameters), 6))
})

design_WDM <- design(data = dat, model = DDM,
                      formula = list(v ~ 0 + S, a ~ E, t0 ~ 1, Z ~ 1),
                      constants = c(s = 0, sv = -Inf, SZ = -Inf, st0 = -Inf))
p_WDM <- p_DDM[!names(p_DDM) %in% c("sv", "SZ", "st0")]
WDM_s <- make_emc(dat, design_WDM, type = "single", compress = FALSE, n_chains = 1)

test_that("WDM likelihood without trial variability", {
  parameters <- likelihood_parameters(p_WDM)
  expect_snapshot(round(model_likelihood(WDM_s, parameters), 6))
})
