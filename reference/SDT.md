# Gaussian Signal Detection Theory Model for Binary Responses

Discrete binary choice based on continuous Gaussian latent, with no rt
(rt must be set to NA in data).

## Usage

``` r
SDT()
```

## Value

A model list with all the necessary functions to sample

## Details

Model parameters are: mean (unbounded) sd (log scale) and threshold
(unbounded).

For identifiability in one condition two parameters must be fixed
(conventionally mean=0 and sd = 1). When used with data that records
only accuracy (so reponse bias cannot be evaluated) a single threshold
must be assumed and fixed (e.g., threshold = 0).

At present this model is not fully implemented in C, but as its
likelihood requires only pnorm evaluation it is quite fast.

## Examples

``` r
dprobit <- design(Rlevels = c("left","right"),
           factors=list(subjects=1,S=c("left","right")),
           formula=list(mean ~ 0+S, sd ~ 1,threshold ~ 1),
           matchfun=function(d)d$S==d$lR,
           constants=c(sd=log(1),threshold=0),
           model=SDT)
#> 
#>  Sampled Parameters: 
#> [1] "mean_Sleft"  "mean_Sright"
#> 
#>  Design Matrices: 
#> $mean
#>      S mean_Sleft mean_Sright
#>   left          1           0
#>  right          0           1
#> 
#> $sd
#>  sd
#>   1
#> 
#> $threshold
#>  threshold
#>          1
#> 

p_vector <- sampled_pars(dprobit)
```
