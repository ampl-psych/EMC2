# EMC2 (development)

## Neural likelihoods

-   The scaled DDM card (DDMnle: a flow over the unit-scale decision time given `kappa`, `eta`, `Z`, `SZ`, with `a`, `s`, `t0` and `st0` applied exactly) is a card kind of its own, `"flow_scaled"`. `register_nn_model(card, twin = DDM)` reads it directly: its parameterisation, bounds, zero exceptions and simulator are the DDM's, its R path is the NLE project's reference (`ddm_sc_density`, `SZ` as the proportion), and its compiled path the evaluator that used to need a hand-built registration from a scaffold card. The two agree to 1e-11 over prior draws in every variability regime.
-   Cards carry a format version (`card_version`; absent = 0). EMC2 reads versions 0 and 1 and refuses any other by name. Version 1 adds the structured `scaling`, `support`/`mass_tolerance`/`beyond_support` (read by `nn_total_mass()`), `control_limits` for race cards (the analytic control's `s` range and `v b / s^2` overflow limit; `nn_cell()` and `nn_in_box()` warn when a prior passes them) and a `validation` record. `vignette("neural-likelihood-cards")` is the specification.
-   `nn_golden_check()` evaluates a card on its golden vectors (`<tag>_golden.csv` beside it; the shipped artefacts' files are now installed with them) and reports the largest difference per output, `PASS` or `FAIL` in one line. `inst/scripts/nle_verify.R` runs it for every shipped artefact and every card given, plus the compiled-vs-R and oracle checks of a scaled card, and exits non-zero on any failure.
-   Card fields are read by exact name. `$` on a list matches a field by prefix, so a card field `preset` was read as `pre` and the card refused; `register_nn_model()` and `nn_golden_check()` now read every card, member and golden-file field with `[[name, exact = TRUE]]`, and a test registers a card carrying extra fields named like prefixes of fields EMC2 reads.
-   `nn_golden_check()` takes its default tolerance from the card's `validation$golden$tolerance` (1e-9 without one; the converter records 1e-4 for a LAN converted from a float32 ONNX file). `inst/scripts/nle_verify.R` checks `mlp_joint` cards too: it registers one with the card's optional `p_types` field or `--p-types=v=0,a=1,...` and checks it to the card's own tolerance.
-   `inst/scripts/onnx_to_card.py` is now a pointer to `nleflow convert` (the converter moved into the nleflow Python package of the NLE repository); its test in `test-nn-mlp.R` is skipped with that message and the converter is tested in nleflow (`tests/test_convert.py`).
-   `nn_sz_convert()` converts between the two start-point-variability conventions (EMC2's proportion; the width `2 SZ min(Z, 1 - Z)` that `DDM()` passes to its likelihood); the conventions are a documented contract with a test at each boundary.
-   The training-region rule of the validation kit is advisory: `nn_cell()`, `nn_in_box()` and `run_sbc()` warn when a prior leaves the network's training region within `k` sd, and no longer refuse (`nn_in_box(error = TRUE)` for the old behaviour).
-   `make_emc()` warns when a neural-likelihood design samples `t0` and the response times are floored to a resolution coarser than 1 ms (the default `rt_resolution = 1/60` biases t0 by about -1.9 posterior sd at N = 100, even for the analytic likelihood; the flooring applies whether or not the data are compressed).
-   On R < 4.5 (no `tools::sha256sum`) the shipped artefacts are pinned by md5 lines in `inst/extdata/flownn/MANIFEST`; before, every shipped artefact failed to load there with a misleading "does not match the sha256" error.
-   `?DDMnn` records the state of the `ddm_cap256w_c4` t0 shift: no correction is applied, an unresolved shift of about 0.002 (log scale) at N = 800 remains, and the scaled card replaces the artefact for the variability models.
-   `register_nn_model()` reads a card's `fixed` field (the values of twin parameters held fixed in training) as the default `exceptions`, so such a card registers, and `nle_verify.R` checks it, from the card alone; an exception that contradicts the card is refused.

## New features (dev-nle)

-   `register_nn_model()` accepts two more kinds of neural likelihood, `"regression_joint"` (direct-regression MLPs) and `"mlp_joint"` (likelihood approximation networks such as HSSM's LANs). Both run in the compiled likelihood, like the flows, with the weights loaded once. `inst/scripts/onnx_to_card.py` converts an ONNX MLP (tanh hidden layers) into a model card; EMC2 does not depend on onnxruntime.

## Sampler

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
