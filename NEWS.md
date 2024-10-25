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
