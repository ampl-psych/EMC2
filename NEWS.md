# EMC2 3.4.1

## New features (dSSfix)

-   `make_kernel()` gains `centre`: the kernel output is centred over the rows it applies to (those with a finite covariate), so the target parameter becomes its value at the *average* covariate instead of at covariate zero. This removes the collinearity between a target parameter and its base weight that otherwise arises whenever the covariate varies over a narrow range away from zero -- stop-signal delays being the motivating case. Rows with a non-finite covariate carry no trend and stay at zero, and centring is per subject.
-   New kernels `sat_incr` / `sat_decr`: the same saturating shape as `slin_incr` / `slin_decr` but parameterised by the saturation point, `k = +/- min(1, c / s_sat)`. `s_sat` is in the units of the covariate ("where the plateau starts"), which is the quantity data can speak to and a far easier thing to put a prior on than a rate. Bound it to the observed covariate range with a `pnorm` transform (`design(transform = list(func = c(<target>.s_sat = "pnorm"), lower = ..., upper = ...))`).

## Bug fixes

-   Every non-sequential kernel (`lin_*`, `exp_*`, `pow_*`, `poly2/3/4`) now treats a non-finite covariate as "no covariate on this row" and contributes no trend there, as the `slin_*` kernels already did. Previously `design()` refused them outright if the covariate contained `NA`, which made them unusable on stop-signal delays (`SSD` is `Inf` on go trials, and `NA` while data are simulated trial by trial).

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
