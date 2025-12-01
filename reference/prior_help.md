# Prior Specification Information

Prints information associated with the prior for certain 'type'

## Usage

``` r
prior_help(type)
```

## Arguments

- type:

  A character string indicating which 'type' of model to run (e.g.
  'standard' or 'single')

## Value

Invisible return with a list of all the information that is also printed

## Examples

``` r
prior_help('diagonal')
#> type: mu 
#> Group-level mean 
#> 
#>   Hyperparameters: 
#>   - theta_mu_mean : mean of the group-level mean prior 
#>   - theta_mu_var : variance of the group-level mean prior 
#> 
#> type: Sigma 
#> Group-level covariance matrix 
#> 
#>   Hyperparameters: 
#>   - v : degrees of freedom on the group-level (co-)variance prior, 2 leads to uniform correlations. Single value 
#>   - A : scale on the group-level variance prior, larger values lead to larger variances 
#> 
```
