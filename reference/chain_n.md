# MCMC Chain Iterations

Returns a matrix with the number of samples per chain for each stage
that is present in the emc object (i.e., `preburn`, `burn`, `adapt`,
`sample`). The number of rows of the matrix reflects the number of
chains and the number of columns the number of sampling stages.

## Usage

``` r
chain_n(emc)
```

## Arguments

- emc:

  A list, the output of
  [`fit()`](https://ampl-psych.github.io/EMC2/reference/fit.md).

## Value

A matrix

## Examples

``` r
chain_n(samples_LNR)
#>      preburn burn adapt sample
#> [1,]       0    0     0     50
#> [2,]       0    0     0     50
#> [3,]       0    0     0     50
```
