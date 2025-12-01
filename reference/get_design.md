# Get Design

Extracts design from an emc object

## Usage

``` r
# S3 method for class 'emc.prior'
get_design(x)

# S3 method for class 'emc'
get_design(x)

get_design(x)
```

## Arguments

- x:

  an `emc` or `emc.prior` object

## Value

A design with class emc.design

## Examples

``` r
get_design(samples_LNR)
#> m ~ lM 
#> s ~ 1 
#> t0 ~ 1 
```
