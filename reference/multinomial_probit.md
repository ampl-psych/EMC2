# Multinomial Probit Response Model

Unordered choice among response alternatives using latent Gaussian
utility noise. Assumes independent Gaussian shocks across alternatives,
which implies correlated pairwise differences in the implied comparison
space.

## Usage

``` r
multinomial_probit()
```

## Value

A model list with all the necessary functions to sample

## Details

Model parameters are: utility (unbounded).

## Examples

``` r
dmnp <- design(
  Rlevels = c("left", "right"),
  factors = list(subjects = 1, S = c("left", "right")),
  formula = list(utility ~ 0 + S),
  matchfun = function(d) d$S == d$lR,
  model = multinomial_probit
)
#> 
#>  Sampled Parameters: 
#> [1] "utility_Sleft"  "utility_Sright"
#> 
#>  Design Matrices: 
#> $utility
#>      S utility_Sleft utility_Sright
#>   left             1              0
#>  right             0              1
#> 
```
