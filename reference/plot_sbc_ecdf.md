# Plot the ECDF Difference in SBC Ranks

Plots F_hat(z) - z where F_hat is the empirical CDF of the (normalized)
ranks and z is the Uniform(0,1) CDF. The shaded band is the simultaneous
(1 - gamma) envelope for F_hat(z) - z.

## Usage

``` r
plot_sbc_ecdf(ranks, layout = NA, add_stats = TRUE, main = NULL, K = 500)
```

## Arguments

- ranks:

  A list of named dataframes of the rank statistic (raw or normalized)

- layout:

  Optional. A numeric vector specifying the layout using
  `par(mfrow = layout)`

- add_stats:

  Boolean. Should coverage, bias and precision be included in the
  figure.

- main:

  Optional. A character specifying plot title.

- K:

  Optional. Effective sample size of the MCMC that produced the ranks.

## Value

No returns
