# Plot Inhibition Functions

Plots panels of the inhibition functions (probability of responding;
Pr(R)) for each level of specified factors of stop-signal data as a
function of user-defined SSD bins/categories. Optionally, posterior
and/or prior predictive inhibition functions can be overlaid.

## Usage

``` r
plot_ss_if(
  input,
  post_predict = NULL,
  prior_predict = NULL,
  probs = seq(0, 1, length.out = 5),
  factors = NULL,
  within_plot = NULL,
  use_global_quantiles = FALSE,
  subject = NULL,
  quants = c(0.025, 0.975),
  functions = NULL,
  n_cores = 1,
  n_post = 50,
  layout = NA,
  to_plot = c("data", "posterior", "prior")[1:2],
  use_lim = c("data", "posterior", "prior")[1:2],
  legendpos = c("topleft", "bottomright"),
  posterior_args = list(),
  prior_args = list(),
  ...
)
```

## Arguments

- input:

  Either an emc object or a stop-signal data frame, or a *list* of such
  objects. SSD column in data required.

- post_predict:

  Optional posterior predictive data (matching columns) or *list*
  thereof.

- prior_predict:

  Optional prior predictive data (matching columns) or list thereof.

- probs:

  Numeric vector of probabilities with values on the unit interval that
  defines SSD bins/categories.

- factors:

  Character vector of factor names to aggregate over; defaults to
  plotting full data set ungrouped by factors if NULL.

- within_plot:

  Character indicating factor for which inhibition functions are plotted
  in the same panel

- use_global_quantiles:

  If set to `TRUE`, SSDs are pooled over participants before calculating
  percentiles, so the same absolute SSD range is used to get Pr(R) for
  each participant, and then these probabilities are averaged over
  participants.

- subject:

  Subset the data to a single subject (by index or name).

- quants:

  Numeric vector of credible interval bounds (e.g. c(0.025, 0.975)).

- functions:

  A function (or list of functions) that create new columns in the
  datasets or predictives

- n_cores:

  Number of CPU cores to use if generating predictives from an emc
  object.

- n_post:

  Number of posterior draws to simulate if needed for predictives.

- layout:

  Numeric vector used in par(mfrow=...); use NA for auto-layout.

- to_plot:

  Character vector: any of "data", "posterior", "prior".

- use_lim:

  Character vector controlling which source(s) define axis limits

- legendpos:

  Character vector controlling the positions of the legends

- posterior_args:

  Optional list of graphical parameters for posterior lines/ribbons.

- prior_args:

  Optional list of graphical parameters for prior lines/ribbons.

- ...:

  Other graphical parameters for the real data lines.

## Value

Returns NULL invisibly

## Details

Per default, the SSD-categories are defined in terms of the percentiles
of the SSD distribution for each participant, and then averaged over
participants (see `use_global_quantiles`).

If credible regions are not plotted, the data is plotted with error bars
(plus/minus the standard error per SSD bin/category)
