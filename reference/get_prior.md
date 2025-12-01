# Get Prior

Extracts prior from an emc object

## Usage

``` r
# S3 method for class 'emc'
get_prior(emc)

get_prior(emc)
```

## Arguments

- emc:

  an emc object

## Value

A prior with class emc.prior

## Examples

``` r
get_prior(samples_LNR)
#> Mean and variance of the prior on the transformed parameters: 
#> m ~ 𝑁(0, 1)
#> m_lMd ~ 𝑁(0, 1)
#> s ~ 𝑁(0, 1)
#> t0 ~ 𝑁(0, 1)
#> 
#> For detailed info use summary(<prior>)
```
