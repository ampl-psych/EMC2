# Plot the Histogram of the Observed Rank Statistics of SBC

Note that this plot is dependent on the number of bins, and a more
general visualization is to use `plot_sbc_ecdf`

## Usage

``` r
plot_sbc_hist(ranks, bins = 10, layout = NA, add_stats = TRUE)
```

## Arguments

- ranks:

  A list of named dataframes of the rank statistic

- bins:

  An integer specifying the number of bins to use when plotting the
  histogram

- layout:

  Optional. A numeric vector specifying the layout using
  `par(mfrow = layout)`

- add_stats:

  Boolean. Should coverage, bias and precision be included in the
  figure.

## Value

No returns
