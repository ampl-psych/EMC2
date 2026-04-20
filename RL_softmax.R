## Full-feedback RL + softmax example.
## For chosen-only feedback, fitting is possible with NA-coded unchosen outcomes,
## but simulation still needs a manual outer loop because rewards depend on the sampled choice.
devtools::load_all(".", quiet = TRUE)

set.seed(123)

reward_map <- function(dadm, cov_names) {
  out <- sapply(cov_names, function(col) {
    option <- sub("^reward_", "", col)
    as.numeric(dadm$lR == option)
  })
  if (!is.matrix(out)) out <- matrix(out, ncol = length(cov_names))
  colnames(out) <- cov_names
  out
}

n_trials <- 200
reward_prob <- c(left = 0.75, right = 0.25)

reward_schedule <- data.frame(
  reward_left = rbinom(n_trials, size = 1, prob = reward_prob["left"]),
  reward_right = rbinom(n_trials, size = 1, prob = reward_prob["right"])
)

trend_rl <- make_trend(
  par_names = "utility",
  cov_names = list(c("reward_left", "reward_right")),
  kernels = "delta",
  bases = "lin",
  maps = list(choice = reward_map),
  at = "lR"
)

design_rl_softmax <- design(
  Rlevels = c("left", "right"),
  factors = list(subjects = 1),
  covariates = c("reward_left", "reward_right"),
  formula = list(utility ~ 1),
  constants = c(utility = 0),
  trend = trend_rl,
  model = multinomial_logit,
)

pars <- sampled_pars(design_rl_softmax)
pars[c("utility.w_choice", "utility.q0", "utility.alpha")] <- c(
  4,
  0.5,
  qnorm(0.1)
)

dat <- make_data(
  pars,
  design_rl_softmax,
  n_trials = n_trials,
  covariates = reward_schedule,
  return_trialwise_parameters = TRUE
)

emc <- make_emc(
  dat,
  design_rl_softmax,
  type = "single",
  n_chains = 2
)

emc <- fit(emc)

print(head(dat, 10))
print(tail(attr(dat, "trialwise_parameters"), 10))

pars
credint(emc)
