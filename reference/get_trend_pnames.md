# Get parameter types from trend object

Get parameter types from trend object

## Usage

``` r
get_trend_pnames(trend)
```

## Arguments

- trend:

  A trend object created by make_trend()

## Value

A character vector of parameter names used in the trend

## Examples

``` r
trend <- make_trend(par_names = "v", cov_names = "trial", kernels = "exp_incr")
get_trend_pnames(trend)
#> [1] "v.w"    "v.d_ei"
```
