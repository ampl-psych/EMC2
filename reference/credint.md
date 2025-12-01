# Posterior Quantiles

Returns the quantiles of the selected parameter type. Full range of
possible samples manipulations described in `get_pars`.

## Usage

``` r
# S3 method for class 'emc.prior'
credint(
  x,
  selection = "mu",
  probs = c(0.025, 0.5, 0.975),
  digits = 3,
  N = 1000,
  covariates = NULL,
  ...
)

# S3 method for class 'emc'
credint(x, selection = "mu", probs = c(0.025, 0.5, 0.975), digits = 3, ...)

credint(x, ...)
```

## Arguments

- x:

  An emc or emc.prior object

- selection:

  A Character vector. Indicates which parameter types to check (e.g.,
  `alpha`, `mu`, `sigma2`, `correlation`).

- probs:

  A vector. Indicates which quantiles to return from the posterior.

- digits:

  Integer. How many digits to round the output to

- N:

  An integer. Number of samples to use for the quantile calculation
  (only for prior.emc objects)

- covariates:

  A list of covariates to use for the quantile calculation (only for
  prior.emc objects)

- ...:

  Optional additional arguments that can be passed to `get_pars`

## Value

A list of posterior quantiles for each parameter group in the selected
parameter type.

## Examples

``` r
credint(samples_LNR)
#> $mu
#>         2.5%    50%  97.5%
#> m     -1.222 -0.968 -0.775
#> m_lMd -0.844 -0.514 -0.121
#> s     -0.939 -0.533 -0.028
#> t0    -2.010 -1.627 -1.137
#> 
```
