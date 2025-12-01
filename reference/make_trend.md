# Create a trend specification for model parameters

Create a trend specification for model parameters

## Usage

``` r
make_trend(
  par_names,
  cov_names = NULL,
  kernels,
  bases = NULL,
  shared = NULL,
  trend_pnames = NULL,
  phase = "premap",
  par_input = NULL,
  at = NULL,
  custom_trend = NULL
)
```

## Arguments

- par_names:

  Character vector specifying which parameters to apply trend to

- cov_names:

  Character vector specifying which covariates to use for each trend

- kernels:

  Character vector specifying which kernel function to use for each
  trend

- bases:

  Optional character vector specifying which base function to use for
  each trend

- shared:

  Named list with entries the parameter names to be shared and the names
  the new names of the shared parameter.

- trend_pnames:

  Optional character vector specifying custom parameter names

- phase:

  Character vector (length 1 or `length(par_names)`) specifying the
  phase for each trend entry; one of "premap", "pretransform", or
  "posttransform". Defaults to "premap".

- par_input:

  Optional character vector(s) of parameter names to use as additional
  inputs for the trend

- at:

  If NULL (default), trend is applied everywhere. If a factor name
  (e.g., "lR"), trend is applied only to entries corresponding to the
  first level of that factor.

- custom_trend:

  A trend registered with `register_trend`

## Value

A list containing the trend specifications for each parameter

## Examples

``` r
# Put trend on B and v parameters
trend <- make_trend(
  par_names = c("B", "v"),
  cov_names = "strial",
  kernels = c("exp_incr", "poly3"),
  phase = "premap",
  shared = list(shrd = list("B.B0", "v.d1"))
)
get_trend_pnames(trend)
#> [1] "B.w"    "B.d_ei" "v.d2"   "v.d3"   "shrd"  
```
