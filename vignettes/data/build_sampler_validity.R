# Builds vignettes/data/sampler-validity.RData for vignette("sampler-validity"):
# one 25-parameter single-subject LNR data set (about 500 trials) fitted with
# the legacy sampler (options(emc.sampler = "legacy")) and with the current
# one, plus the Laplace SDs at the posterior mode. About 2 minutes on 3 cores.
# Run from the package root with the branch installed:
#   Rscript vignettes/data/build_sampler_validity.R [library path]
args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(path.expand(args[1]), .libPaths()))
library(EMC2)
suppressPackageStartupMessages(library(mvtnorm))

matchfun <- function(d) (d$NWS == "new" & d$lR == "NEW") | (d$NWS != "new" & d$lR == "OLD")
des <- design(factors = list(subjects = 1, NWS = c("new", "weak", "strong"), WF = c("low", "high")),
              Rlevels = c("NEW", "OLD"), model = LNR, matchfun = matchfun,
              formula = list(m ~ 0 + NWS:WF:lM, s ~ 0 + NWS:WF:lM, t0 ~ 1),
              report_p_vector = FALSE)
pn <- names(sampled_pars(des))
# subject means of a real participant (Ratcliff & Starns, 2013 data; DMC fit)
truth <- c("m_NWSnew:WFlow:lMFALSE" = -0.0322, "m_NWSweak:WFlow:lMFALSE" = 0.1606,
  "m_NWSstrong:WFlow:lMFALSE" = 0.6663, "m_NWSnew:WFhigh:lMFALSE" = -0.3178,
  "m_NWSweak:WFhigh:lMFALSE" = -0.0628, "m_NWSstrong:WFhigh:lMFALSE" = -0.0539,
  "m_NWSnew:WFlow:lMTRUE" = -0.9061, "m_NWSweak:WFlow:lMTRUE" = -1.0254,
  "m_NWSstrong:WFlow:lMTRUE" = -1.2573, "m_NWSnew:WFhigh:lMTRUE" = -0.5986,
  "m_NWSweak:WFhigh:lMTRUE" = -0.8527, "m_NWSstrong:WFhigh:lMTRUE" = -0.9149,
  "s_NWSnew:WFlow:lMFALSE" = -0.1493, "s_NWSweak:WFlow:lMFALSE" = -0.0386,
  "s_NWSstrong:WFlow:lMFALSE" = 0.1036, "s_NWSnew:WFhigh:lMFALSE" = -0.2355,
  "s_NWSweak:WFhigh:lMFALSE" = -0.1724, "s_NWSstrong:WFhigh:lMFALSE" = -0.2418,
  "s_NWSnew:WFlow:lMTRUE" = -0.3291, "s_NWSweak:WFlow:lMTRUE" = -0.0401,
  "s_NWSstrong:WFlow:lMTRUE" = -0.1724, "s_NWSnew:WFhigh:lMTRUE" = -0.2883,
  "s_NWSweak:WFhigh:lMTRUE" = -0.132, "s_NWSstrong:WFhigh:lMTRUE" = -0.1426,
  t0 = -0.9032)[pn]
pmean <- setNames(numeric(length(pn)), pn); psd <- setNames(rep(1, length(pn)), pn)
pmean[grep("lMTRUE", pn)] <- -0.7
pmean[grep("^s_", pn)] <- log(.8); psd[grep("^s_", pn)] <- .5
pmean["t0"] <- log(.3); psd["t0"] <- .5
pri <- prior(des, type = "single", pmean = pmean, psd = psd)

set.seed(2026)
dat <- make_data(truth, des, n_trials = 84)      # 6 cells x 84 = 504 trials

fit_one <- function(sampler) {
  options(emc.sampler = sampler)
  on.exit(options(emc.sampler = NULL))
  emc <- make_emc(dat, des, type = "single", prior_list = pri, rt_resolution = 1e-9, n_chains = 3)
  fit(emc, cores_for_chains = 3, verbose = FALSE, trim = FALSE)
}
fits <- list(legacy = fit_one("legacy"), current = fit_one(NULL))

# Laplace SDs: numerical Hessian of log-likelihood + log prior at the mode
dadm <- fits$current[[1]]$data[[1]]; mod <- fits$current[[1]]$model
ll <- function(p) as.vector(EMC2:::calc_ll_manager(matrix(p, 1, dimnames = list(NULL, pn)), dadm, mod))
obj <- function(p) -(ll(p) + dmvnorm(p, pmean, diag(psd^2), log = TRUE))
X0 <- t(get_pars(fits$current, selection = "alpha", stage = "sample", return_mcmc = FALSE, merge_chains = TRUE)[, 1, ])
mode <- optim(colMeans(X0), obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12))$par
k <- length(mode); h <- 1e-3; H <- matrix(0, k, k)
for (a in 1:k) for (b in a:k) {
  ea <- replace(numeric(k), a, h); eb <- replace(numeric(k), b, h)
  H[a, b] <- H[b, a] <- (obj(mode + ea + eb) - obj(mode + ea - eb) - obj(mode - ea + eb) + obj(mode - ea - eb)) / (4 * h^2)
}
laplace <- list(mode = mode, sd = sqrt(diag(solve(H))), ll_mode = ll(mode))

summarise_fit <- function(emc) {
  X <- t(get_pars(emc, selection = "alpha", stage = "sample", return_mcmc = FALSE, merge_chains = TRUE)[, 1, ])
  D <- -2 * unlist(lapply(emc, function(ch) ch$samples$subj_ll[1, ch$samples$stage == "sample"]))
  A <- emc[[1]]$samples$alpha[, 1, emc[[1]]$samples$stage == "sample"]
  moved <- colSums(abs(A[, -1] - A[, -ncol(A)])) > 0
  pm <- attr(emc[[1]]$samples, "pm_settings")[[1]][[1]]
  list(pars = data.frame(par = pn, truth = truth, post_mean = colMeans(X), post_sd = apply(X, 2, sd),
                         laplace_sd = laplace$sd, z = (colMeans(X) - truth) / apply(X, 2, sd)),
       excess = mean(D) + 2 * laplace$ll_mode, k = k, move_rate = mean(moved),
       epsilon = pm$epsilon, mix = pm$mix, n_particles = pm$n_particles,
       iters = chain_n(emc)[1, ],
       max_rhat = max(gd_summary(emc, selection = "alpha", stage = "sample", print_summary = FALSE)))
}
summaries <- lapply(fits, summarise_fit)
save(fits, summaries, laplace, truth, pmean, psd, file = "vignettes/data/sampler-validity.RData")
for (nm in names(summaries)) with(summaries[[nm]], cat(nm, ": excess/k", round(excess / k, 2),
  " median SD/Laplace", round(median(pars$post_sd / pars$laplace_sd), 2), " move", round(move_rate, 2),
  " eps", paste(round(epsilon, 2), collapse = "/"), " Rhat", round(max_rhat, 3), "\n"))
