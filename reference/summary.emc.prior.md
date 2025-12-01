# Summary method for emc.prior objects

Prints a summary of the prior specification, including descriptions of
the prior types and their associated hyperparameters.

## Usage

``` r
# S3 method for class 'emc.prior'
summary(object, ...)
```

## Arguments

- object:

  An object of class 'emc.prior' containing prior specifications

- ...:

  Additional arguments passed to other methods (not currently used)

## Value

Invisibly returns NULL. Called for its side effect of printing the
summary.

## See also

[`prior`](https://ampl-psych.github.io/EMC2/reference/prior.md) for
creating prior objects

## Examples

``` r
# Take a prior object
prior <- get_prior(samples_LNR)
summary(prior)
#> mu  -  Group-level mean 
#> 
#> mean of the group-level mean prior : 
#>     m m_lMd     s    t0 
#>     0     0     0     0 
#> variance of the group-level mean prior : 
#>       m m_lMd s t0
#> m     1     0 0  0
#> m_lMd 0     1 0  0
#> s     0     0 1  0
#> t0    0     0 0  1
#> 
#> Sigma  -  Group-level covariance matrix 
#> 
#> degrees of freedom on the group-level (co-)variance prior, 2 leads to uniform correlations. Single value : 
#> [1] 2
#> scale on the group-level variance prior, larger values lead to larger variances : 
#>     m m_lMd     s    t0 
#>   0.3   0.3   0.3   0.3 
#> 
```
