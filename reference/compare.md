# Information Criteria and Marginal Likelihoods

Returns the BPIC/DIC or marginal deviance (-2\*marginal likelihood) for
a list of samples objects.

## Usage

``` r
compare(
  sList,
  stage = "sample",
  filter = NULL,
  use_best_fit = TRUE,
  BayesFactor = TRUE,
  cores_for_props = 4,
  cores_per_prop = 1,
  print_summary = TRUE,
  digits = 0,
  digits_p = 3,
  ...
)
```

## Arguments

- sList:

  List of samples objects

- stage:

  A string. Specifies which stage the samples are to be taken from
  `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`

- filter:

  An integer or vector. If it's an integer, iterations up until the
  value set by `filter` will be excluded. If a vector is supplied, only
  the iterations in the vector will be considered.

- use_best_fit:

  Boolean, defaults to `TRUE`, uses the minimal or mean likelihood
  (whichever is better) in the calculation, otherwise always uses the
  mean likelihood.

- BayesFactor:

  Boolean, defaults to `TRUE`. Include marginal likelihoods as estimated
  using WARP-III bridge sampling. Usually takes a minute per model added
  to calculate

- cores_for_props:

  Integer, how many cores to use for the Bayes factor calculation, here
  4 is the default for the 4 different proposal densities to evaluate,
  only 1, 2 and 4 are sensible.

- cores_per_prop:

  Integer, how many cores to use for the Bayes factor calculation if you
  have more than 4 cores available. Cores used will be cores_for_props
  \* cores_per_prop. Best to prioritize cores_for_props being 4 or 2

- print_summary:

  Boolean (default `TRUE`), print table of results

- digits:

  Integer, significant digits in printed table for information criteria

- digits_p:

  Integer, significant digits in printed table for model weights

- ...:

  Additional, optional arguments

## Value

Matrix of effective number of parameters, mean deviance, deviance of
mean, DIC, BPIC, Marginal Deviance (if `BayesFactor=TRUE`) and
associated weights.

## Examples

``` r
# \donttest{
compare(list(samples_LNR), cores_for_props = 1)
#>     MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
#> 1 -571   1 -610    1 -571     1         39  -648  -656 -687
# Typically we would define a list of two (or more) different models:
# # Here the full model is an emc object with the hypothesized effect
# # The null model is an emc object without the hypothesized effect
# design_full <- design(data = forstmann,model=DDM,
#                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#                            constants=c(s=log(1)))
# # Now without a ~ E
# design_null <- design(data = forstmann,model=DDM,
#                            formula =list(v~0+S,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
#                            constants=c(s=log(1)))
#
# full_model <- make_emc(forstmann, design_full)
# full_model <- fit(full_model)
#
# null_model <- make_emc(forstmann, design_null)
# null_model <- fit(null_model)
# sList <- list(full_model, null_model)
# # By default emc uses 4 cores to parallelize marginal likelihood estimation across proposals
# # So cores_per_prop = 3 results in 12 cores used.
# compare(sList, cores_per_prop = 3)
# }
```
