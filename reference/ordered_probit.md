# Ordered Probit Response Model

Ordered responses based on a latent Gaussian evidence variable. Binary
response models are the 2-category special case.

## Usage

``` r
ordered_probit()
```

## Value

A model list with all the necessary functions to sample

## Details

Model parameters are: location (unbounded), scale (log scale), and cut
(first threshold free, remaining thresholds modeled as positive
increments).

The `cut ~ 1` specification yields the standard flexible `K - 1`
threshold parameterization from the ordinal regression literature. The
final response category has an implicit upper cut at `Inf`. Internally,
`Ttransform` leaves sampled `cut` values unchanged and returns the
ordered thresholds in a derived `cut_expanded` column used by the
likelihood and random-number generator.

## Examples

``` r
dord <- design(
  Rlevels = c("left", "right"),
  factors = list(subjects = 1, S = c("left", "right")),
  formula = list(location ~ 0 + S, scale ~ 1, cut ~ 1),
  matchfun = function(d) d$S == d$lR,
  constants = c(scale = log(1)),
  model = ordered_probit
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
