# Bayes Factors

returns the Bayes Factor for two models

## Usage

``` r
get_BayesFactor(MLL1, MLL2)
```

## Arguments

- MLL1:

  Numeric. Marginal likelihood of model 1. Obtained with
  [`run_bridge_sampling()`](https://ampl-psych.github.io/EMC2/reference/run_bridge_sampling.md)

- MLL2:

  Numeric. Marginal likelihood of model 2. Obtained with
  [`run_bridge_sampling()`](https://ampl-psych.github.io/EMC2/reference/run_bridge_sampling.md)

## Value

The BayesFactor for model 1 over model 2

## Examples

``` r
# \donttest{
# Normally one would compare two different models
# Here we use two times the same model:
M1 <- M0 <- run_bridge_sampling(samples_LNR, both_splits = FALSE, cores_for_props = 1)
get_BayesFactor(M1, M0)
#> [1] 1
# }
```
