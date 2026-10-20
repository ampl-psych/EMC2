# EMC2 (development)

## Sampler

-   `make_emc()` defaults to 6 chains for hierarchical types and 3 for `type = "single"` (`n_chains = NULL`; a number given is used as it is). Previously 3. The default is a recommendation, not a validity requirement: a fit that met the stop rule with 3 chains stands.
-   Sample-stage stopping reads every criterion (the iteration target, `min_es`, the opt-in criteria below and the flat test) on one set of draws. The first third of the stage is dropped only when doing so makes the stop pass, and each drop is recorded; otherwise every draw is kept and the stage grows. `max_tries` now ends an unmet sample stage on its own. The progress line prints the stop rule's value (what `max_gd` is compared with), alongside the plain maximum when that differs. `stop_criteria$min_es` is no longer partially matched by `min_ess`.
-   New opt-in sample-stage criteria in `stop_criteria`, each off unless given: `min_ess` (smallest effective size over the group means and log group SDs, summed over chains), `max_ll_gap` (best minus worst chain in median summed subject log-likelihood) and `max_subj_spread` (largest range over chains of a subject's median log-likelihood). `abort_gd` with `abort_try` (default 20) ends a sample stage early when, at that check, the stop is unmet and the rule's Rhat exceeds `abort_gd`. The fit records how the sample stage ended (criteria on, met, tries, the rule's value, last statistics, discards, culls, `aborted`) in `emc[[1]]$stop_record`.
-   New opt-in `cull` (burn and sample): a chain held in a different region is respawned from a donor chain, with the donor's kernel settings. A chain is flagged when it sits more than `chain_gap` (200) below the best chain on the summed subject log-likelihood, at z above `chain_z` (5) within-chain spreads, for `persistence` (2) checks, and either more than `lp_guard_below` (1000) below or also lower on its own log posterior. Burn compares from 250 iterations, re-seeds up to twice, and re-seeds any chain still apart at burn exit. Sample compares from 500 iterations; a cull restarts the stage after `reburn` (100) iterations, at most once. Culls are recorded in `emc[[1]]$culls`.
-   New opt-in `max_flat_loc` (sample): the group means, subject parameters, summed and per-subject log-likelihoods must not drift between the stage's first and last thirds by more than that many interquartile ranges (0.3 is a good value); `flat_mcse` optionally passes a drift within that many Monte Carlo errors of zero (default 0, off). `chain_diagnostics()` reports the stop statistics, the cull comparison and the flat test for any fitted object.
-   `thin_to_ess()` thins the sample stage of a fitted object to keep about a target effective sample size (default 3000) without discarding information from slower-mixing parameters, recording what it did in `emc[[1]]$thinning`.
-   The joint location move's curvature-guided directions are built once per stage (and once more when the sample kernel is tuned) instead of at every step, so their step sizes adapt continuously across the steps of burn and adapt. The sample kernel's directions are rebuilt from adapt's tail, where the step restarts at 1 and is frozen at its average over the tail's last 80 iterations. The acceptance ratio is unchanged, and the sample-stage kernel stays fixed. A prior covariance of the group means that cannot be Cholesky-factorised falls back to an eigenvalue square root, and unusable curvature falls back to prior-scaled axes, instead of stopping with an Armadillo error.
-   `adapt` and `sample` now use exact ensemble and importance-resampling kernels, fixing biased posterior spreads and stuck chains. `preburn` and `burn` retain the search kernel. `options(emc.sampler = "legacy")` restores the previous sampler for comparison.
-   Independent particle draws now follow their mixture probabilities without forcing a draw from every component. Parameter blocks update sequentially, and shared likelihoods are counted once.
-   Forked sampling now advances L'Ecuyer-CMRG streams between batches and uses separate participant substreams.
-   Sample proposals are built from `adapt`, tuned for 100 more adapt iterations, then fixed. Hierarchical adaptation checks subject-level Rhat over the last 250 draws (target 1.2; stops after three checks without improvement or 1000 iterations). Disable this check with `options(emc.adapt_converge = FALSE)`.
-   Hierarchical subject proposals follow the current group mean and covariance using estimated likelihood precision. `lik_precision()` brackets finite-difference steps and avoids spurious curvature at model bounds.
-   Standard hierarchical models add exact scale moves and, without a group-level design, location moves to improve funnel mixing. Tuning ends before sampling. Disable with `options(emc.scale_move = FALSE)`; state is stored in `attr(samples, "scale_move")`.
-   Summaries, plots and stopping rules now share split-Rhat without coda's correction or transform. Stop checks are computed in one pass. `stop_criteria$gd_quantile` optionally uses a subject-level quantile; other parameters still use their largest Rhat.
-   Singular group-covariance draws are retried up to 3 times; `on_singular = list(max_retries = 0)` disables retries.
-   Windows workers now receive the `emc.*` options.
-   SBC plot bands use directly computed binomial quantiles and root-finding for simultaneous coverage.
-   The long posterior-preservation test runs locally and is skipped on CRAN and CI.

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
