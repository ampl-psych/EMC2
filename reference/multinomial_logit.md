# Multinomial Logit Response Model

Unordered choice among response alternatives using a softmax over latent
utilities.

## Usage

``` r
multinomial_logit()
```

## Value

A model list with all the necessary functions to sample

## Details

Model parameters are: utility (unbounded).

## Examples

``` r
dmnl <- design(
  Rlevels = c("left", "right"),
  factors = list(subjects = 1, S = c("left", "right")),
  formula = list(utility ~ 0 + S),
  matchfun = function(d) d$S == d$lR,
  model = multinomial_logit
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
