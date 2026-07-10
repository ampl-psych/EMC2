# Plots trends over time

Plots the trend for selected parameters of a model. Can be used either
with a p_vector, or trial-wise parameters or covariates obtained from
predict(). Optionally overlays a ground-truth trajectory for recovery
evaluation.

## Usage

``` r
plot_trend(
  input_data,
  emc,
  par_name,
  true_trajectories = NULL,
  subject = NULL,
  filter = NULL,
  on_x_axis = "trials",
  pp_shaded = TRUE,
  group_average = FALSE,
  layout = NULL,
  ...
)
```

## Arguments

- input_data:

  A p_vector or posterior predictives compatible with the provided emc
  object. If posterior predictives, should be the output of
  `predict(..., return_trialwise_parameters = TRUE)`.

- emc:

  An emc object.

- par_name:

  Parameter name (or covariate name) to plot.

- true_trajectories:

  Optional. A data frame with columns `subjects`, the column named by
  `on_x_axis`, and at least `par_name`. If supplied, the true trajectory
  is overlaid as a dashed black line, and `ylim` expands to include it.

- subject:

  Subject to plot. If NULL, plots all subjects separately. Ignored when
  `group_average = TRUE`.

- filter:

  Optional function that takes a data frame and returns a logical vector
  indicating which rows to include. For race models, should select one
  accumulator, e.g. `function(d) d$lR == 'left'`.

- on_x_axis:

  Column name to plot on the x-axis. Default `'trials'`.

- pp_shaded:

  Logical. If TRUE, plots 95% credible interval as a shaded area.
  Otherwise plots separate lines for each posterior predictive
  iteration. Only applicable when `input_data` is posterior predictives.

- group_average:

  Logical. If TRUE, plots a single panel averaged across subjects rather
  than one panel per subject.

- layout:

  Optional integer vector of length 2 to override the automatic mfrow
  layout. Ignored when `group_average = TRUE`.

- ...:

  Optional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns the credible interval data frame when `input_data` is
posterior predictives, otherwise NULL.
