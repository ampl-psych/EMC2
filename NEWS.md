# EMC2 (development)

## Sampler (dev-sampler)

-   The particle step used in every stage put proposal components centred on the chain's current value into one importance-weighted batch, which is not a valid MCMC kernel: in 25-35 dimensions the draws followed the proposal rather than the posterior (SD 1.25-1.4x the Laplace SD in converged fits, stuck chains with 2000+ trials per subject; hierarchical fits affected too). `adapt` and `sample` now use an exact kernel (per iteration either an ensemble move around the current value or an importance-resampling step with the independent components), the local step size is adapted to a 3% "better-than-current" target with a windowed acceptance count, independence proposals are used at their estimated scale, and the `diag(.5)` fallback for an unusable chain covariance is gone. `preburn`/`burn` keep the old step as a search. `options(emc.sampler = "legacy")` restores the previous sampler for comparison. Calibration on a 48-cell grid, the MTLNR SBC and the Ratcliff-Starns Table 1 fits: `vignette("sampler-validity")`.
-   The sample stage now runs one fixed kernel. Its proposals are built from the `adapt` draws and its step size, mixing weights and number of particles are tuned in a 100-iteration tail of `adapt`; nothing adapts once draws are kept. Before, the proposals were re-estimated at every step from the chain's last 250 draws and the step size kept adapting, which is not a valid sampler: on a conjugate hierarchical model with weakly identified parameters, chains started at the truth sank into the funnel (group SD median 0.004 against 0.16 exact). For models with a group level the local proposal and one independence proposal now follow the current group mean and covariance (subject likelihood precision, by finite differences, plus the current group precision), so the subject step no longer freezes when a group variance collapses. On that model long runs now reproduce the exact posterior from a healthy and from a collapsed start. `options(emc.sampler = "legacy")` still restores the old sampler. `vignette("sampler-validity")`, sections 4 and 8, is corrected accordingly.
-   Hierarchical (`type = "standard"`) fits with few trials per weakly identified parameter per subject used to stop at `fit()`'s maximum number of tries with the group variances of those parameters collapsed (Rhat 1.2-2.4): the centred alternation of the subject step and the group step cannot leave the funnel's neck (autocorrelation time of a log group SD about 500 iterations even with an exact subject draw). After the group step, `adapt` and `sample` now run an interweaving sweep (Yu & Meng, 2011): for each parameter a scale move that rescales the group SD, every subject's deviation from the group mean and the prior's auxiliary variable together, and a location move that shifts the group mean and every subject together. The sweep runs on the quadratic approximation of each subject's likelihood that the subject step already has and is then accepted against the exact likelihoods (one evaluation per subject per iteration), so it is exact whatever the approximation. Step sizes are tuned in the tail of `adapt` and frozen for `sample`. On a conjugate model with the exact posterior known, the autocorrelation time of the weak parameters' log group SD falls from about 500 to about 10 and of their group means to 4-5, and `fit()` converges; `options(emc.scale_move = FALSE)` turns the sweep off. The location move only runs without a group-level design. See `vignette("sampler-validity")` section 7 for how it works and section 8 for two open issues found in regression testing (a step-size runaway on a weakly-identified DDM parameter, and a real-data hierarchical fit where the sweep worsens subject-level mixing with no funnel present) that are under active investigation and not yet fixed in this release; the sweep's own ~20% cost estimate also does not generalise past models shaped like the one it was tuned on, and can run 1.5-11x per sample iteration elsewhere.
-   The interweaving sweep's step sizes now adapt on the realised acceptance of a move (accepted on the quadratic surrogate *and* its block then passed the exact likelihood check), not on the surrogate's acceptance alone, which tuned the step to the surrogate whatever its quality. On forstmann's DDM the surrogate for `sv` and `SZ` was wrong by three orders of magnitude (below), its acceptance was a coin toss on the proposal's sign at any step, the step ran to 6-11 log-SD units, and the exact check rejected nearly every such block, leaving those group SDs to the slow centred random walk (Rhat 5.3 at `fit()`'s maximum tries). The exact check's rejections now bound the step by what the exact likelihoods accept. Where the surrogate is exact (a Gaussian likelihood) every block passes and nothing changes: the conjugate model's draws are bit-identical. `attr(samples, "scale_move")` gains `acc_real` / `acc_loc_real` (realised acceptances; `acc_in` / `acc_loc` keep counting the surrogate's). The conjugate model, the `k25` funnel reproducer and the grid cells nearest their margins are unchanged in convergence and mixing. **forstmann itself is not fixed by this** (`sv` Rhat 2.4 at maximum tries): the step now shrinks instead of running away, because the surrogate for `sv`/`SZ` is wrong by three orders of magnitude -- `lik_precision()`'s finite-difference step search oscillates between a step on the flat stretch and one that straddles the DDM's likelihood floor at its parameter bounds (`sv` < .01, `SZ` < .01), and reads the floor as curvature. That also narrows the subject step's local proposal for those parameters to ~0.01 and is the pre-existing cause of the `sv` convergence failure on this data set. A fix to that search is prototyped and pending (`rating-work/sampler/hier/stageH5/REPORT.md`).
-   `options(emc.ll_backend = "multithreaded")` no longer hangs a fit on Linux. With gcc's OpenMP runtime (libgomp), a process that has run a threaded likelihood and then forks leaves its children waiting for ever, so a threaded likelihood in the R session itself (`compare()`, WAIC, a profile plot, a stage run with `cores_for_chains = 1` and `cores_per_chain = 1`) followed by a fit with parallel chains in the same session hung without an error. The R session now releases its OpenMP thread pool after each threaded likelihood (about 0.1 ms per call at 4 threads, 0.6 ms at 16; forked chains and subject workers keep theirs and lose nothing). Where the pool cannot be released, processes forked afterwards use the serial likelihood, with a warning. Nothing changes on macOS, on Windows or with Intel's OpenMP runtime, which were not affected. `emc2_build_info()` reports the OpenMP runtime in use.

## New features (dev-rating)

-   New model `MTLNR()`, the correlated multiple-threshold log-normal race of Reynolds, Kvam, Osth & Heathcote (2020), for choice and response-time data with an additional ordered confidence (or other) rating. `MTLNR(n_ratings = K)` sets the number of rating categories; ratings are read from the design's `RR` column (1 = lowest, `K` = highest) and reported on the natural threshold scale with `add_recalculated = TRUE`. See `vignette("rating-models")`.
-   New rating-data helpers, shared across future rating models: `rating_summary()` and `plot_ratings()` (response proportions and RT quantiles by folded response), `zroc()` and `plot_zroc()` (z-transformed ROC).

## New features (dev-nle)

-   `register_nn_model()` accepts two more kinds of neural likelihood, `"regression_joint"` (direct-regression MLPs) and `"mlp_joint"` (likelihood approximation networks such as HSSM's LANs). Both run in the compiled likelihood, like the flows, with the weights loaded once. `inst/scripts/onnx_to_card.py` converts an ONNX MLP (tanh hidden layers) into a model card; EMC2 does not depend on onnxruntime.

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
