# Effective Sample Size

Returns the effective sample size (ESS) of the selected parameter type.
Full range of possible samples manipulations described in `get_pars`.

## Usage

``` r
# S3 method for class 'emc'
ess_summary(
  emc,
  selection = "mu",
  stat = "min",
  stat_only = FALSE,
  digits = 1,
  ...
)

ess_summary(emc, ...)
```

## Arguments

- emc:

  An emc object

- selection:

  A Character vector. Indicates which parameter types to check (e.g.,
  `alpha`, `mu`, `sigma2`, `correlation`).

- stat:

  A string. Should correspond to a function that can be applied to a
  vector, which will be performed on the vector/rows or columns of the
  matrix of the parameters

- stat_only:

  Boolean. If `TRUE` will only return the result of the applied stat
  function, otherwise returns both the stat result and the result of the
  function on all parameters.

- digits:

  Integer. How many digits to round the output to

- ...:

  Optional additional arguments that can be passed to `get_pars`

## Value

A matrix or vector of ESS values for the selected parameter type.

## Examples

``` r
ess_summary(samples_LNR, selection = "alpha")
#>        as1t  bd6t  bl1t   min
#> m     387.5 139.8 181.1 139.8
#> m_lMd 150.0 150.0 192.9 150.0
#> s     180.0 599.6 252.9 180.0
#> t0    402.3 158.4 256.2 158.4
#> min   150.0 139.8 181.1 139.8
```
