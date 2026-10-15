# EMC2 (development)

## Sampler

-   The sampler runs in C++. Chains and subjects are OpenMP threads (`cores_for_chains` chains at a time, `cores_per_chain` subjects each), so `OMP_NUM_THREADS` does not cap them, and with Intel's runtime `KMP_AFFINITY=disabled` is needed or every thread shares one CPU. Models whose likelihood must stay on the R thread (`CDM`, `PSDM`, `PHSDM`, custom trend kernels, builds with `PNORM_MODE=0`) run one chain and one subject at a time, with one warning that the core settings are ignored. Likelihoods implemented in R can no longer be sampled.
-   `adapt` and `sample` now use exact ensemble and importance-resampling kernels, fixing biased posterior spreads and stuck chains. `preburn` and `burn` retain the search kernel.
-   Independent particle draws now follow their mixture probabilities without forcing a draw from every component. Parameter blocks update sequentially, and shared likelihoods are counted once.
-   Each chain draws from its own generator, seeded from R's. `set.seed()` before `fit()` or `run_emc()` determines the draws, a saved object continues its own stream (kept in `attr(samples, "rng")`), and the core settings do not change the draws, so the same seed on the same object reproduces a run.
-   Sample proposals are built from `adapt`, tuned for 100 more adapt iterations, then fixed. Hierarchical adaptation checks subject-level Rhat over the last 250 draws (target 1.2; stops after three checks without improvement or 1000 iterations). Disable this check with `options(emc.adapt_converge = FALSE)`.
-   A stage can be run again on any saved or finished object, including more `adapt` after its tail has tuned the sample kernel, so a fit killed at any point is resumed with `fit()` on the saved file.
-   Hierarchical subject proposals follow the current group mean and covariance using estimated likelihood precision. `lik_precision()` brackets finite-difference steps and avoids spurious curvature at model bounds.
-   Every hierarchical type (`standard`, `diagonal`, `infnt_factor`, `SEM`) adds exact scale moves and, without a group-level design, location moves to improve funnel mixing. Tuning ends before sampling. Disable with `options(emc.scale_move = FALSE)`; state is stored in `attr(samples, "scale_move")`.
-   Group-level corrections: `infnt_factor` refreshes the shrinkage product after every `delta` draw and draws the factor scores before the group mean (a valid partially collapsed order); its starting `delta` uses the rate parameterisation. `SEM` regressions account for fixed non-zero loadings and paths (a marker loading fixed at 1), the subject level is conditional on the factor scores and each subject's covariates (`theta_var` holds the population covariance), cyclic `B` gets its Jacobian correction, the start is centred on the prior mean, and a factor group of `d` factors needs `a_d > d - 1` (otherwise its prior is improper, and the fit stops before sampling).
-   Summaries, plots and stopping rules now share split-Rhat without coda's correction or transform. Stop checks are computed in one pass. `stop_criteria$gd_quantile` optionally uses a subject-level quantile; other parameters still use their largest Rhat.
-   Singular group-level draws, including the `SEM` and `infnt_factor` decompositions, are retried up to 3 times; `on_singular = list(max_retries = 0)` disables retries. `on_singular$on_exhausted` is accepted from earlier scripts with a warning that it is ignored.
-   A parameter block whose component has no likelihood, a samples array too small for the draws, or `is_blocked`/`par_group` shorter than the parameters stop with an error instead of sampling the prior or corrupting memory.
-   Fork safety on Linux with GNU libgomp: the pool of OpenMP threads is released when a threaded call returns, so `fit()` in forked processes (`run_sbc()`, `mclapply()`) after a fit in the session no longer hangs. Intel's and LLVM's runtimes need nothing and are left alone. `emc2_build_info()` reports the runtime.
-   Build: `$(FLIBS)` is no longer linked (the package compiles no Fortran; on macOS without gfortran the link failed), every object depends on every header (an in-place rebuild after a header change no longer mixes stale and new objects), and the sources compile with Intel icpc 2021.
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
