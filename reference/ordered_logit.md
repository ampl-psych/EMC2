# Ordered Logit Response Model

Ordered responses based on a latent logistic evidence variable. Binary
response models are the 2-category special case.

## Usage

``` r
ordered_logit()
```

## Value

A model list with all the necessary functions to sample

## Details

The `scale` parameter is the latent logistic scale, and thresholds
follow the same ordered `K - 1` parameterization as
[`ordered_probit()`](https://ampl-psych.github.io/EMC2/reference/ordered_probit.md),
with derived ordered thresholds returned in `cut_expanded`.

## Examples

``` r
dord <- design(
  Rlevels = c("left", "right"),
  factors = list(subjects = 1, S = c("left", "right")),
  formula = list(location ~ 0 + S, scale ~ 1, cut ~ 1),
  matchfun = function(d) d$S == d$lR,
  constants = c(scale = log(1)),
  model = ordered_logit
)
#> 
#>  Sampled Parameters: 
#> [1] "location_Sleft"  "location_Sright" "cut"            
#> 
#>  Design Matrices: 
#> $location
#>      S location_Sleft location_Sright
#>   left              1               0
#>  right              0               1
#> 
#> $scale
#>  scale
#>      1
#> 
#> $cut
#>  cut
#>    1
#> 
```
