# Plot the ECDF Difference in SBC Ranks

Plots the difference in observed cumulative rank statistics and the
expected cumulative distribution of a uniform distribution. The blue
shaded areas indicate the 95% credible interval.

## Usage

``` r
plot_sbc_ecdf(ranks, layout = NA, add_stats = TRUE, main = "")
```

## Arguments

- ranks:

  A list of named dataframes of the rank statistic

- layout:

  Optional. A numeric vector specifying the layout using
  `par(mfrow = layout)`

- add_stats:

  Boolean. Should coverage, bias and precision be included in the
  figure.

- main:

  Optional. A character specifying plot title.

## Value

No returns
