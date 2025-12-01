# Cut Factors Based on Credible Loadings

This function removes factors that do not have more than one credible
loading based on the specified confidence interval.

## Usage

``` r
cut_factors(emc, CI = 95)
```

## Arguments

- emc:

  An 'emc' object containing factor analysis results

- CI:

  Numeric. Confidence interval percentage (default is 95)

## Value

An 'emc' object with factors that don't meet the credibility criterion
removed
