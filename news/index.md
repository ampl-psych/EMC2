# Changelog

## EMC2 3.3.0

### Bug fixes

- Addressed a major issue in `map = TRUE`, so that now mapping of
  population level parameters is now done correctly. See issue
  [\#119](https://github.com/ampl-psych/EMC2/issues/119)

### New features

- More general support for `group_design` (tutorial on the way), also
  using `map = TRUE`

- Start of a more general trends implementation (stay tuned; tutorial
  still on the way)

- Some more SEM/FA functionality added i.e. `rotate_loadings`

## EMC2 3.2.1

CRAN release: 2025-09-22

### Bug fix

- Added a warning that in hierarchical models, using `map = TRUE` in
  functions like `map = TRUE` or `map = TRUE` does not return the
  population-level marginal mean and variances on the original scale for
  group-level parameters. See issue
  [\#119](https://github.com/ampl-psych/EMC2/issues/119)

## EMC2 3.2.0

CRAN release: 2025-06-27

### New features

- group_design specification (tutorial coming up)

- trends specification (tutorial also coming up)

- broader continuous covariates support with `map = TRUE`

- fMRI joint modelling (tutorial on
  <https://osf.io/preprints/psyarxiv/rhfk3_v1>)

- Made changes to how lR was used in compressed likelihood for race
  models. Run `update2version` to reuse older samples

### Bug fixes

- Made legend including or excluding more flexible for `plot_cdf`,
  `plot_density` and `plot_stat`

## EMC2 3.1.1

CRAN release: 2025-04-07

### New features

- added thin to fit/run_emc which can either be set to TRUE to
  automatically thin based on ESS, or on a numeric to only keep 1/x
  samples

- added probit/SDT model for bimanual choices

### Bug fixes

- Rare bug in sampling removed

- Small bug fixes in plot_data to make it more flexible

- cleared up argumentation of run_emc/fit

## EMC2 3.1.0

CRAN release: 2025-03-10

### New features

- model_averaging function, which allows you to compare evidence for an
  effect across a set of models

### Bug fixes

- Small hotfix in which creating proposals and the start of burn would
  sometimes fail for large number of subjects

- Patched up old error in which model bounds weren’t considered in data
  generation

- Fixed error in which compare_subject would return IC for whole dataset
  for every subject.

## EMC2 3.0.0

CRAN release: 2025-03-07

### New features

- IMPORTANT: to keep your old samples compatible with current EMC2, run
  update2version()

- IMPORTANT: Design and prior are now also S3 methods with their own S3
  classes, see EMC2 paper

- plot_fit is deprecated and has branched of into plot_density, plot_cdf
  and plot_stat

- sampled_p_vector is deprecated and is now named sampled_pars

- Added a design_plot function, which makes a plot of the proposed
  accumulation process

- Sampling is completely reworked. Adaptive tuning of the number of
  particles and more stable convergence

### Bug Fixes

- Fixed rare case where conditional MVN would break

- Fixed bug in predict on joint models

- Suppressed unwanted print statements in DDM estimation/prediction

- Added more checks to a wide array of functions to ensure proper input
  format

## EMC2 2.1.0

CRAN release: 2024-10-14

### New features

- Added a website with vignettes, changelog and a reference

- Added
  [`run_sbc()`](https://ampl-psych.github.io/EMC2/reference/run_sbc.md)
  function to perform simulation-based calibration for a design

- Added
  [`prior_help()`](https://ampl-psych.github.io/EMC2/reference/prior_help.md)
  to get more information on the prior for a certain `type`

- Changed DDM implementation, which is faster and more accurate

- Bridge sampling now also works for `type = "blocked"`

### Bug Fixes

- Made bridge sampling for inverse-gamma and inverse-wishart more robust

- Made [`prior()`](https://ampl-psych.github.io/EMC2/reference/prior.md)
  function work more generally
