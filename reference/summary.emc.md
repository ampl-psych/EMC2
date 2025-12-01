# Summary Statistics for emc Objects

Computes quantiles, `Rhat` and `ESS` for selected model parameters.

## Usage

``` r
# S3 method for class 'emc'
summary(
  object,
  selection = c("mu", "sigma2", "alpha"),
  probs = c(0.025, 0.5, 0.975),
  digits = 3,
  ...
)
```

## Arguments

- object:

  An object of class `emc`

- selection:

  A character string indicating the parameter type Defaults to `mu`,
  `sigma2`, and `alpha`. See below for more information.

- probs:

  The quantiles to be computed. Defaults to the the 2.5%, 50% and 97.5%
  quantiles.

- digits:

  An integer specifying rounding of output.

- ...:

  Optional arguments that can be passed to `get_pars`

## Value

A list of summary output.

## Details

Note that if `selection = alpha` and `by_subject = TRUE` (default) is
used, summary statistics are computed at the individual level. to the
console but summary statistics for all subjects are returned by the
function.

If an emc object that has not been run with `fit` yet is supplied,
summary of the design will be returned.
