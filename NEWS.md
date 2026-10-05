# EMC2 (development)

## Sampler (dev-sampler)

-   The particle step used in every stage put proposal components centred on the chain's current value into one importance-weighted batch, which is not a valid MCMC kernel: in 25-35 dimensions the draws followed the proposal rather than the posterior (SD 1.25-1.4x the Laplace SD in converged fits, stuck chains with 2000+ trials per subject; hierarchical fits affected too). `adapt` and `sample` now use an exact kernel (per iteration either an ensemble move around the current value or an importance-resampling step with the independent components), the local step size is adapted to a 3% "better-than-current" target with a windowed acceptance count, independence proposals are used at their estimated scale, and the `diag(.5)` fallback for an unusable chain covariance is gone. `preburn`/`burn` keep the old step as a search. `options(emc.sampler = "legacy")` restores the previous sampler for comparison. Calibration on a 48-cell grid, the MTLNR SBC and the Ratcliff-Starns Table 1 fits: `vignette("sampler-validity")`.
-   The sample stage now runs one fixed kernel. Its proposals are built from the `adapt` draws and its step size, mixing weights and number of particles are tuned in a 100-iteration tail of `adapt`; nothing adapts once draws are kept. Before, the proposals were re-estimated at every step from the chain's last 250 draws and the step size kept adapting, which is not a valid sampler: on a conjugate hierarchical model with weakly identified parameters, chains started at the truth sank into the funnel (group SD median 0.004 against 0.16 exact). For models with a group level the local proposal and one independence proposal now follow the current group mean and covariance (subject likelihood precision, by finite differences, plus the current group precision), so the subject step no longer freezes when a group variance collapses. On that model long runs now reproduce the exact posterior from a healthy and from a collapsed start. `options(emc.sampler = "legacy")` still restores the old sampler. `vignette("sampler-validity")`, sections 4 and 8, is corrected accordingly.
-   Hierarchical (`type = "standard"`) fits with few trials per weakly identified parameter per subject used to stop at `fit()`'s maximum number of tries with the group variances of those parameters collapsed (Rhat 1.2-2.4): the centred alternation of the subject step and the group step cannot leave the funnel's neck (autocorrelation time of a log group SD about 500 iterations even with an exact subject draw). After the group step, `adapt` and `sample` now run an interweaving sweep (Yu & Meng, 2011): for each parameter a scale move that rescales the group SD, every subject's deviation from the group mean and the prior's auxiliary variable together, and a location move that shifts the group mean and every subject together. The sweep runs on the quadratic approximation of each subject's likelihood that the subject step already has and is then accepted against the exact likelihoods (one evaluation per subject per iteration), so it is exact whatever the approximation. Step sizes are tuned in the tail of `adapt` and frozen for `sample`. On a conjugate model with the exact posterior known, the autocorrelation time of the weak parameters' log group SD falls from about 500 to about 10 and of their group means to 4-5, and `fit()` converges; `options(emc.scale_move = FALSE)` turns the sweep off. The location move only runs without a group-level design. See `vignette("sampler-validity")` section 7 for how it works and section 8 for two open issues found in regression testing (a step-size runaway on a weakly-identified DDM parameter, and a real-data hierarchical fit where the sweep worsens subject-level mixing with no funnel present) that are under active investigation and not yet fixed in this release; the sweep's own ~20% cost estimate also does not generalise past models shaped like the one it was tuned on, and can run 1.5-11x per sample iteration elsewhere.
-   The interweaving sweep's step sizes now adapt on the realised acceptance of a move (accepted on the quadratic surrogate *and* its block then passed the exact likelihood check), not on the surrogate's acceptance alone, which tuned the step to the surrogate whatever its quality. On forstmann's DDM the surrogate for `sv` and `SZ` was wrong by three orders of magnitude (below), its acceptance was a coin toss on the proposal's sign at any step, the step ran to 6-11 log-SD units, and the exact check rejected nearly every such block, leaving those group SDs to the slow centred random walk (Rhat 5.3 at `fit()`'s maximum tries). The exact check's rejections now bound the step by what the exact likelihoods accept. Where the surrogate is exact (a Gaussian likelihood) every block passes and nothing changes: the conjugate model's draws are bit-identical. `attr(samples, "scale_move")` gains `acc_real` / `acc_loc_real` (realised acceptances; `acc_in` / `acc_loc` keep counting the surrogate's). The conjugate model and the `k25` funnel cells are unchanged in convergence and mixing; a move shares the fate of its block, so the target of the realised acceptance is .3 times a running block acceptance (`r_block` in the same attribute: exponentially weighted, starting at 1, frozen with the steps, and never taken below .35, the block acceptance under which the sweep is split into more blocks). With a fixed target of .3 through blocks that pass about half the time every step in a block shrank to a quarter of what the surrogate allows, which cost group-SD mixing on the grid's hardest cell (`m35_standard_500`, 35 parameters, 8 subjects; five replicates per build: 9-12 group SDs with an autocorrelation time over 40 against 2-6 before, largest Rhat of a log group SD 1.20-1.33 against 1.09-1.21); with the scaled target that cell is back where it was (3-7 parameters, 1.10-1.20, median autocorrelation time 13 against 15) and its subject-level Rhat is lower than under either earlier rule (1.08-1.16 against 1.13-1.65), the `k25` cells are unchanged, and forstmann's `sv` and `SZ` group SDs mix about twice as fast with no runaway (steps of 0.3-2.0). On its own this did not fix forstmann (`sv` Rhat 2.4 at maximum tries: the step shrank to nothing instead of running away), which exposed the real fault, next.
-   The subject likelihood precision (`lik_precision()`, the quadratic surrogate that the subject step's local proposal and the interweaving sweep both run on) could be wrong by three orders of magnitude for a parameter whose likelihood is flat on one side and floored on the other: its finite-difference step search oscillated for ever between a step on the flat stretch (grow x10) and one across the floor (shrink x0.1), and used whichever side the last round landed on. The DDM returns `min_ll` per trial outside its bounds (`sv`, `SZ` < .01), forstmann's subjects sit a log unit or two above that floor, and the floor was read as a curvature of 10^3-10^4 with a gradient to match: the local proposal for a subject's `SZ` had SD ~0.01 against a posterior 0.2-0.5 wide (the cause of the `sv` convergence failure on this data set that predates the sweep), and every sweep move on `sv`/`SZ` failed its exact check. The search now brackets the step and bisects in log scale, and a parameter whose step never settles gets a flat surrogate (no precision, no gradient). On forstmann `SZ` now converges (Rhat 1.03 at 5000 draws) and `sv` goes from Rhat 5.3 to 1.2-1.7 at `fit()`'s maximum tries and 1.18 at 5000 draws, with every chain visiting the whole of its posterior (every other parameter is below 1.01, in half the time); what is left on `sv` is the slow random walk of its group mean across a two-log-unit plateau where the likelihood is indifferent, which is the model and prior rather than the sampler. Where the likelihood is quadratic (the conjugate model) nothing changes, and on the LNR cells of the grid the surrogate is identical entry for entry.
-   Hierarchical fits with many subjects: the fixed sample-stage kernel was built from `adapt` draws taken while the chains still disagreed (adapt stopped on `min_unique` after about 100 iterations), and the few subjects whose kernel was poor set the fit's Rhat. On the Eisenberg Simon data (RDM, 518 subjects) every one of 20 fits stopped above 1.1 at `fit()`'s maximum tries, with a minimum subject-level ESS of 14 per 1000 draws. Four changes. (1) `adapt` also waits until, over the last 250 iterations (the draws the kernel is built from), the largest Rhat of any subject's parameter is below 1.2, checked from 250 iterations on and given up after three checks without improvement or at 1000; `options(emc.adapt_converge = FALSE)` restores the old rule. (2) For `type = "standard"` the particle step's conditional proposal takes the subject likelihood's precision from the window's draws (their inverse covariance minus the group precision, with the finite-difference estimate kept where the draws show only the prior) instead of the curvature at one point, which was two to three times too narrow for skewed subject posteriors. It does this only where the group precision is stable over the window (largest coefficient of variation below .5); in a funnel it is not, the subtraction is wrong, and the finite-difference precision is kept. (3) The interweaving sweep's steps start at each move's conditional scale, and a gate drops, for the rest of the fit, moves that add little to the group-level Gibbs step once the sweep is split into one block per parameter (with several hundred subjects the surrogate's error adds up and moves of well-identified parameters only pass at tiny steps). (4) Each block's exact check (one likelihood per subject) is timed in-process against forked and the faster way is kept: forking cost more than the likelihoods themselves. On Eisenberg with `fit()` defaults: subject-level Rhat 1.10 [1.08, 1.16] over five runs (3 of 5 below 1.1), ESS 5% quantile / minimum 211 / 48, 63 minutes (the sweep used to take 150); agreement with a Stan (NUTS) reference on every group-level quantity, which the earlier sampler that passed this data set (`002a5c6f`, re-estimating its proposals during sampling) did not have (it was off on 10 of 16 group means and SDs, and too narrow). On the replicated ladder (two `k25` funnel cells, `m35_standard_500`, forstmann) subject-level mixing is as good or better (conditional-proposal acceptance about doubled, fewer particles) and runs are 2-4 times faster, except forstmann (DDM), 1.3 times slower because its `sv` plateau keeps adapt running. `vignette("sampler-validity")`, sections 7-8.
-   One definition of Rhat, and a stop rule that reads it in one pass. `fit()`'s `burn` and `sample` stop rules, `gd_summary()` and `summary()` now use the split-Rhat (each chain halved; the potential scale reduction factor of Gelman & Rubin over the half-chains) instead of the point estimate of `coda::gelman.diag(transform = TRUE)`. coda's estimate multiplies by a degrees-of-freedom correction that is large whenever the chains' variances differ, so a subject whose posterior has a shoulder that one chain happens to visit read 1.15-1.2 where the split-Rhat reads 1.05, and with several hundred subjects the largest of thousands of such values decided the fit: on the Eisenberg Simon data (518 subjects, 4144 subject-level parameters) the sample stage's criterion was met by 3 of 20 runs within `fit()`'s 20 checks, although every run's posterior agreed with a Stan (NUTS) reference. Reported Rhat values are therefore somewhat lower than before, mostly for the parameters with the highest values. The rule itself is unchanged: checked every `step_size` draws, the largest Rhat of the selection below 1.1 with at least 1000 draws kept, and on a failed check the first third of the stage's draws is dropped if the largest Rhat is lower without it (chains still arriving at the posterior; replays of the rule on runs with their burn-in put back in front showed that this discard is needed -- chains that drift together pass Rhat < 1.1 -- and that the same rule on the split-Rhat removes such a transient as well as before). A check is now one pass over the stored draws (0.5 s on Eisenberg against 58-108 s, which was more than the sampling between two checks). Replayed on 20 saved Eisenberg runs the new rule stops on its criterion in 18, after a median 1380 draws instead of 2000; in five new `fit()` runs it did so in 3 (as the old rule did in its five), in 37 minutes against 63, with every group-level quantity within 2.4 Monte-Carlo standard errors of Stan. The two `k25` funnel cells and `m35_standard_500` run 1.6-1.9 times faster with the same or better Rhat and the same effective sizes per draw (`m35`: 5 of 5 runs below 1.1 against 3 of 5); forstmann's DDM is 1.25 times faster and otherwise within its run-to-run spread (1 of 5 runs below 1.1 under either rule: its `sv` plateau needs more draws than `fit()`'s limit); the conjugate model's exact posterior is recovered as before. New option `stop_criteria$gd_quantile` (e.g. `.99`, not the default): the subject-level parameters are judged by that quantile of their Rhats instead of the largest, every other selected parameter still by the largest, and the discard still by the largest; with it the criterion does not tighten as the number of subjects grows. See `?fit`, `?gd_summary` and `vignette("sampler-validity")`, section 7.
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
