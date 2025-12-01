# Plot conditional accuracy functions

Plots panels of conditional accuracy functions (CAFs, one for each level
of caf_factor on the same panel). Accuracy is calculated with smoothing
box car filter on percentile ranges, 0..X, 1..(X+1), ... , (100-X+1)..
Inf, where 1 \< X \<= 50. Optionally, posterior and/or prior predictive
CAFs can be overlaid.

## Usage

``` r
plot_caf(
  input,
  post_predict = NULL,
  prior_predict = NULL,
  subject = NULL,
  quants = c(0.025, 0.975),
  functions = NULL,
  factors = NULL,
  caf_factor = NULL,
  n_cores = 1,
  n_post = 50,
  layout = NA,
  to_plot = c("data", "posterior", "prior")[1:2],
  use_lim = c("data", "posterior", "prior")[1:2],
  legendpos = c("bottomleft", "bottomright"),
  posterior_args = list(),
  prior_args = list(),
  accuracy_function = function(d) d$S == d$R,
  smooth_window = 5,
  which_plot = 1:2,
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

- caf_factor:

  The name of within-panel factor

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

- accuracy_function:

  Accuracy score, default: function(d) d\$S==d\$R,

- smooth_window, :

  range of RT over which calculate accuracy, default 5

- which_plot:

  which of levels of caf_factor to plot, default is both i.e,.
  which_plot = 1:2

- ...:

  Other graphical parameters for the real data lines.

## Value

Returns `NULL` invisibly.

## Examples

``` r
# Plot conditional accuracy function for data only,
# NB: the caf_factor must have two levels levels.
# forstmann_speed_accuracy <- forstmann[forstmann$E!="neutral",]
# forstmann_speed_accuracy$E <- droplevels(forstmann_speed_accuracy$E)
# plot_caf(forstmann_speed_accuracy, caf_factor="E",factors="S", smooth_window=10)
#
# Or a list of multiple emc objects ...
```
