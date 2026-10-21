matchfun <- function(d) as.numeric(d$S)==as.numeric(d$lR) |
  (d$lR=="pm" & as.numeric(d$S)>2)
designLBA <- design(
  factors=list(subjects=1,S=c("left","right","leftpm","rightpm"),RACE=2:3),
  Rlevels=c("left","right","pm"),
  matchfun=matchfun,
  model=LBA,constants=c(v_RACE3=0),
  functions=list(missingness=function(d) ifelse(d$RACE==2 & d$lR=="pm", 0, NA)),
  formula=list(v~RACE*lM,B~1,t0~1,A~1),
)

p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(2), log(4), log(1), log(2), log(0.2),log(.5))

dat <- data.frame(
  subjects = factor(1),
  S = factor(c("left", "right", "left", "right", "leftpm", "rightpm"),
             levels = c("left", "right", "leftpm", "rightpm")),
  RACE = factor(c(2, 2, 3, 3, 3, 3)),
  R = factor(c("left", "right", "right", "left", "pm", "right"),
             levels = c("left", "right", "pm")),
  rt = c(.6, .8, .7, 1.1, .9, .65)
)
emc <- make_emc(dat, designLBA, type = "single", compress = FALSE, n_chains = 1)

test_that("RACE likelihood", {
  parameters <- likelihood_parameters(p_vector)
  expect_snapshot(round(model_likelihood(emc, parameters), 6))
})
