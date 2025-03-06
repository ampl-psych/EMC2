# EMC2 3.0.0

## New features

* IMPORTANT: to keep your old samples compatible with current EMC2, run update2version(<name of old samples>)

* IMPORTANT: Design and prior are now also S3 methods with their own S3 classes, see EMC2 paper

* plot_fit is deprecated and has branched of into plot_density, plot_cdf and plot_stat

* sampled_p_vector is deprecated and is now named sampled_pars

* Added a design_plot function, which makes a plot of the proposed accumulation process

* Sampling is completely reworked. Adaptive tuning of the number of particles 
and more stable convergence

## Bug Fixes

* Fixed rare case where conditional MVN would break

* Fixed bug in predict on joint models

* Suppressed unwanted print statements in DDM estimation/prediction

* Added more checks to a wide array of functions to ensure proper input format

# EMC2 2.1.0

## New features 

* Added a website with vignettes, changelog and a reference

* Added `run_sbc()` function to perform simulation-based calibration for a design

* Added `prior_help()` to get more information on the prior for a certain `type`

* Changed DDM implementation, which is faster and more accurate

* Bridge sampling now also works for `type = "blocked"`

## Bug Fixes

* Made bridge sampling for inverse-gamma and inverse-wishart more robust

* Made `prior()` function work more generally
