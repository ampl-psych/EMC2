# Information Criteria For Each Participant

Returns the BPIC/DIC based model weights for each participant in a list
of samples objects

## Usage

``` r
compare_subject(
  sList,
  stage = "sample",
  filter = 0,
  use_best_fit = TRUE,
  print_summary = TRUE,
  digits = 3,
  return_summary = FALSE,
  n_cores = 1,
  subject = NULL
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

  Boolean, defaults to `TRUE`, use minimal likelihood or mean likelihood
  (whichever is better) in the calculation, otherwise always uses the
  mean likelihood.

- print_summary:

  Boolean (defaults to `TRUE`) print tables of model weight results

- digits:

  Integer, significant digits in printed table

- return_summary:

  Return tables of model weight results

- n_cores:

  Number of cores for parallel processing

- subject:

  Used to select subset of subjects (integer or character vector)

## Value

List of matrices for each subject of effective number of parameters,
mean deviance, deviance of mean, DIC, BPIC and associated weights.

## Examples

``` r
# For a broader illustration see `compare`.
# Here we just take two times the same model, but normally one would compare
# different models
compare_subject(list(m0 = samples_LNR, m1 = samples_LNR))
#> ...$DIC
#>      wDIC_m0 wDIC_m1
#> as1t     0.5     0.5
#> bd6t     0.5     0.5
#> bl1t     0.5     0.5
#> 
#> $BPIC
#>      wBPIC_m0 wBPIC_m1
#> as1t      0.5      0.5
#> bd6t      0.5      0.5
#> bl1t      0.5      0.5
#> 
#> 
#> Winners
#>      m0
#> DIC   3
#> BPIC  3
```
