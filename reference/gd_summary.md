# Gelman-Rubin Statistic

Returns the Gelman-Rubin diagnostics (otherwise known as the R-hat) of
the selected parameter type; i.e. the ratio of between to within MCMC
chain variance.

## Usage

``` r
# S3 method for class 'emc'
gd_summary(
  emc,
  selection = "mu",
  omit_mpsrf = TRUE,
  stat = "max",
  stat_only = FALSE,
  digits = 3,
  ...
)

gd_summary(emc, ...)
```

## Arguments

- emc:

  An emc object

- selection:

  A Character vector. Indicates which parameter types to check (e.g.,
  `alpha`, `mu`, `sigma2`, `correlation`).

- omit_mpsrf:

  Boolean. If `TRUE` also returns the multivariate point scale reduction
  factor (see
  [`?coda::gelman.diag`](https://rdrr.io/pkg/coda/man/gelman.diag.html)).

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

A matrix or vector of R-hat values for the selected parameter type.

## Details

See: Gelman, A and Rubin, DB (1992) Inference from iterative simulation
using multiple sequences, *Statistical Science*, 7, 457-511.

Full range of possible samples manipulations described in `get_pars`.

## Examples

``` r
gd_summary(samples_LNR, selection = "correlation", stat = "mean", flatten = TRUE)
#>             m_lMd.m   s.m  t0.m s.m_lMd t0.m_lMd  t0.s  mean
#> correlation   1.014 1.012 0.997   1.024    1.025 1.013 1.014
```
