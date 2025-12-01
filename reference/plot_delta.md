# Plot Difference of Cumulative Distribution Functions

Plots panels of differences in cumulative distribution functions (CDFs)
between conditions specified by the delta factor in the data.
Optionally, posterior and/or prior predictive delta functions can be
overlaid.

## Usage

``` r
plot_delta(
  input,
  post_predict = NULL,
  prior_predict = NULL,
  subject = NULL,
  quants = c(0.025, 0.975),
  functions = NULL,
  factors = NULL,
  delta_factor = "R",
  n_cores = 1,
  n_post = 50,
  layout = NA,
  to_plot = c("data", "posterior", "prior")[1:2],
  use_lim = c("data", "posterior", "prior")[1:2],
  legendpos = c("topleft"),
  posterior_args = list(),
  prior_args = list(),
  add_percentiles = c(1:9) * 10,
  rev_delta = FALSE,
  ...
)
```

## Arguments

- input:

  Either an `emc` object or a data frame, or a *list* of such objects.

- post_predict:

  Optional posterior predictive data (matching columns) or *list*
  thereof.

- prior_predict:

  Optional prior predictive data (matching columns) or *list* thereof.

- subject:

  Subset the data to a single subject (by index or name).

- quants:

  Numeric vector of credible interval bounds (e.g. `c(0.025, 0.975)`).

- functions:

  A function (or list of functions) that create new columns in the
  datasets or predictives

- factors:

  Character vector of factor names to aggregate over; defaults to
  plotting full data set ungrouped by factors if `NULL`.

- delta_factor:

  The name of the factor to delta

- n_cores:

  Number of CPU cores to use if generating predictives from an `emc`
  object.

- n_post:

  Number of posterior draws to simulate if needed for predictives.

- layout:

  Numeric vector used in `par(mfrow=...)`; use `NA` for auto-layout.

- to_plot:

  Character vector: any of `"data"`, `"posterior"`, `"prior"`.

- use_lim:

  Character vector controlling which source(s) define `xlim`.

- legendpos:

  Character vector controlling the positions of the legends

- posterior_args:

  Optional list of graphical parameters for posterior lines/ribbons.

- prior_args:

  Optional list of graphical parameters for prior lines/ribbons.

- add_percentiles:

  Vector of integers giving percentiles to plot as points, NULL stops
  plotting.

- rev_delta:

  If FALSE (the default) the first level of the defective factor is
  subtracted from the second, if TRUE this is reversed.

- ...:

  Other graphical parameters for the real data lines.

## Value

Returns `NULL` invisibly.

## Examples

``` r
# Plot delta function for data only, not that the delta_factor must have two
# levels.
# fortsmann_speed_accuracy <- forstmann[forstmann$E!="neutral",]
# fortsmann_speed_accuracy$E <- droplevels(fortsmann_speed_accuracy$E)
# plot_delta(fortsmann_speed_accuracy, to_plot = "data")
#
# Plot with posterior predictions
# plot_delta(samples_LNR, to_plot = c("data","posterior"), n_post=10)
#
# Or a list of multiple emc objects ...
```
