RNGkind("L'Ecuyer-CMRG")
set.seed(123)

matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR) |
  (d$lR=="pm" & as.numeric(d$S)>2)
designLBA <- design(
  factors=list(subjects=1,S=c("left","right","leftpm","rightpm"),RACE=2:3),
  Rlevels=c("left","right","pm"),
  matchfun=matchfun,
  model=LBA,constants=c(v_RACE3=0),
  formula=list(v~RACE*lM,B~1,t0~1,A~1),
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(2), log(4), log(1), log(2), log(0.2),log(.5))

# Make square data so can remove pm in RACE = 2
template <- make_data(p_vector,designLBA,n_trials=10)
template <- template[!(template$RACE==2 & (template$S %in% c("leftpm","rightpm"))),]
dat <- make_data(p_vector,designLBA,data=template)

emc <- make_emc(dat, designLBA, type = "single", compress = F)

test_that("RACE", {
  expect_snapshot(init_chains(emc, particles = 10, cores_for_chains = 1))
})
