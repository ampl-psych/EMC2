# Rotate loadings based on posterior median

This function rotates factor loadings using a rotation function based on
the posterior median.

## Usage

``` r
rotate_loadings(emc, rot_fun)
```

## Arguments

- emc:

  An 'emc' object containing factor analysis results

- rot_fun:

  a rotation function for factor loadings, see also `GPArotation`W

## Value

An 'emc' object with rotated factor loadings
