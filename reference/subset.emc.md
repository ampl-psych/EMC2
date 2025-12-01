# Shorten an emc Object

Shorten an emc Object

## Usage

``` r
# S3 method for class 'emc'
subset(
  x,
  stage = "sample",
  filter = NULL,
  thin = 1,
  keep_stages = FALSE,
  length.out = NULL,
  ...
)
```

## Arguments

- x:

  an emc object

- stage:

  A character string. Indicates from which sampling stage(s) to take the
  samples from (i.e. `preburn`, `burn`, `adapt`, `sample`)

- filter:

  Integer or numeric vector. If an integer is supplied, iterations up
  until that integer are removed. If a vector is supplied, the
  iterations within the range are kept.

- thin:

  An integer. By how much to thin the chains

- keep_stages:

  Boolean. If `TRUE`, will not remove samples from unselected stages.

- length.out:

  Integer. Alternatively to thinning, you can also select a desired
  length of the MCMC chains, which will be thinned appropriately.

- ...:

  additional optional arguments

## Value

A shortened emc object

## Examples

``` r
subset(samples_LNR, length.out = 10)
#> Iterations: 
#>      preburn burn adapt sample
#> [1,]       0    0     0     10
#> [2,]       0    0     0     10
#> [3,]       0    0     0     10
#> 
#> Subjects: 
#> [1] "as1t" "bd6t" "bl1t"
#> 
#> Parameters: 
#> [1] "m"     "m_lMd" "s"     "t0"   
```
