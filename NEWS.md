# EMC2 (development)

## Sampler

-   `adapt` and `sample` use an exact particle step: each iteration is either an ensemble move around the current value or importance resampling with proposals that do not depend on it. The previous step put proposals centred on the current value into one importance-weighted batch, which is not a valid MCMC kernel: with 25-35 parameters posterior SDs came out 1.25-1.7 times too wide, and chains with much data per subject got stuck. `preburn` and `burn` keep it as a fast search. `options(emc.sampler = "legacy")` restores the previous sampler for comparison.
-   The sample stage runs one fixed kernel, built from the end of `adapt` and tuned in a 100-iteration tail of it (step size and number of particles); nothing adapts once draws are kept. For models with a group level `adapt` waits until its last 250 draws have converged (largest Rhat of the subject parameters below 1.2, or no longer improving, at most 1000 iterations); `options(emc.adapt_converge = FALSE)` switches this off.
-   For models with a group level the subject step's proposals follow the current group mean and covariance, through each subject's likelihood precision, so they no longer freeze when a group variance collapses. Where the group level is stable the precision comes from the subject's recent draws. `lik_precision()` brackets its finite-difference step, and gives a parameter whose likelihood is floored at a model bound (e.g. the DDM's `sv`, `SZ`) no precision instead of a spurious curvature.
-   Hierarchical (`type = "standard"`) fits run an interweaving sweep after the group step (Yu & Meng, 2011): a scale move per parameter (group SD and every subject's deviation together) and, without a group-level design, a location move (group mean and every subject together). The moves are proposed on each subject's quadratic likelihood surrogate and accepted against the exact likelihoods, so the chain targets the exact posterior. This removes the slow mixing of weakly identified parameters in the hierarchical funnel. Step sizes adapt in `adapt` and are frozen for `sample`; moves that add little are dropped. `options(emc.scale_move = FALSE)` turns the sweep off; its state is in `attr(samples, "scale_move")`.
-   One definition of Rhat: `gd_summary()`, `summary()`, plots and `fit()`'s `burn` and `sample` stop rules use the split-Rhat (each chain halved, Gelman-Rubin over the halves, no correction or transform) instead of the point estimate of `coda::gelman.diag(transform = TRUE)`. Values are lower than coda's, most for the highest ones (a subject whose posterior has a shoulder that one chain visits reads about 1.05 instead of 1.15-1.2), so burn and sample stop sooner. The stop rule is otherwise unchanged, including dropping the first third of the stage's draws when that lowers the largest Rhat, and a check is one pass over the stored draws (half a second with 518 subjects). New option `stop_criteria$gd_quantile` (e.g. `.99`, off by default) judges the subject level by that quantile of its Rhats instead of the largest.
-   `fit()` re-draws a numerically singular group-covariance draw up to 3 times before the chain errors (`on_singular`; `list(max_retries = 0)` errors at once).
-   On Windows the `emc.*` options (including `emc.sampler`) reach the processes that run the chains.
-   SBC: `get_gamma()` / `get_lims()` take the binomial quantile from `cumsum(dbinom())` instead of `qbinom()`, which is numerically wrong on some R builds.
-   New `vignette("sampler-validity")`: why the previous sampler was replaced, what the current one does, and the evidence that it samples the right posterior (exact references, Stan's NUTS on nine hard cases, and forstmann's five standard models).
-   Tests: the sampler's long statistical checks run only with `EMC2_SLOW_TESTS=true`; the default run keeps fast unit tests of the same code.

# EMC2 3.4.1

## New features (cens_trunc2-SS-dEXG3mu)

-   New built-in trend kernels `slin_incr` (`k = min(1, k_sat * c)`) and `slin_decr` (`k = -min(1, k_sat * c)`); non-finite covariates give 0. With a `lin` base and `design(transform = list(func = c(<target>.w = "exp")))` they give a rise (or fall) that saturates, e.g. the dEXG3 stop-signal model (`muS = muS(0) + d * min(1, k * SSD)`).
-   `make_ssd()` generators can be given to `design(functions = list(SSD = ...))`: they are pre-trial design functions that run the staircase trial by trial in `make_data()` (`conditional_on_data = FALSE`), so parameters that depend on SSD see the realised SSD. `make_ssd()` gains a `UC` argument for the late-response rule. `make_data(functions = )` keeps the vectorised staircase for models without a trend on SSD.
-   `get_data()` keeps `make_ssd()` columns (observed SSDs are data), so `predict(conditional_on_data = TRUE)` predicts at the observed delays.

## Bug fixes

-   Small makevars corrections for new CRAN checks

## New features

-   Trend parameters are now also returned with map = TRUE

# EMC2 3.4.0

## Bug fixes

-   Fixed some maths that was wrong about the ECDF plot for SBC

## New features

-   Finished general trends implementation (stay tuned; tutorial still on the way)

# EMC2 3.3.0

## Bug fixes

-   Addressed a major issue in `map = TRUE`, so that now mapping of population level parameters is now done correctly. See issue #119

## New features

-   More general support for `group_design` (tutorial on the way), also using `map = TRUE`

-   Start of a more general trends implementation (stay tuned; tutorial still on the way)

-   Some more SEM/FA functionality added i.e. `rotate_loadings`

# EMC2 3.2.1

## Bug fix

-   Added a warning that in hierarchical models, using `map = TRUE` in functions like `map = TRUE` or `map = TRUE` does not return the population-level marginal mean and variances on the original scale for group-level parameters. See issue #119

# EMC2 3.2.0

## New features

-   group_design specification (tutorial coming up)

-   trends specification (tutorial also coming up)

-   broader continuous covariates support with `map = TRUE`

-   fMRI joint modelling (tutorial on <https://osf.io/preprints/psyarxiv/rhfk3_v1>)

-   Made changes to how lR was used in compressed likelihood for race models. Run `update2version` to reuse older samples

## Bug fixes

-   Made legend including or excluding more flexible for `plot_cdf`, `plot_density` and `plot_stat`

# EMC2 3.1.1

## New features

-   added thin to fit/run_emc which can either be set to TRUE to automatically thin based on ESS, or on a numeric to only keep 1/x samples

-   added probit/SDT model for bimanual choices

## Bug fixes

-   Rare bug in sampling removed

-   Small bug fixes in plot_data to make it more flexible

-   cleared up argumentation of run_emc/fit

# EMC2 3.1.0

## New features

-   model_averaging function, which allows you to compare evidence for an effect across a set of models

## Bug fixes

-   Small hotfix in which creating proposals and the start of burn would sometimes fail for large number of subjects

-   Patched up old error in which model bounds weren't considered in data generation

-   Fixed error in which compare_subject would return IC for whole dataset for every subject.

# EMC2 3.0.0

## New features

-   IMPORTANT: to keep your old samples compatible with current EMC2, run update2version(<name of old samples>)

-   IMPORTANT: Design and prior are now also S3 methods with their own S3 classes, see EMC2 paper

-   plot_fit is deprecated and has branched of into plot_density, plot_cdf and plot_stat

-   sampled_p_vector is deprecated and is now named sampled_pars

-   Added a design_plot function, which makes a plot of the proposed accumulation process

-   Sampling is completely reworked. Adaptive tuning of the number of particles and more stable convergence

## Bug Fixes

-   Fixed rare case where conditional MVN would break

-   Fixed bug in predict on joint models

-   Suppressed unwanted print statements in DDM estimation/prediction

-   Added more checks to a wide array of functions to ensure proper input format

# EMC2 2.1.0

## New features

-   Added a website with vignettes, changelog and a reference

-   Added `run_sbc()` function to perform simulation-based calibration for a design

-   Added `prior_help()` to get more information on the prior for a certain `type`

-   Changed DDM implementation, which is faster and more accurate

-   Bridge sampling now also works for `type = "blocked"`

## Bug Fixes

-   Made bridge sampling for inverse-gamma and inverse-wishart more robust

-   Made `prior()` function work more generally
