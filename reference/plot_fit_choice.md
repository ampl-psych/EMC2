# Plot Choice Model Fit

Plot observed and predicted fit for choice-only response models.

## Usage

``` r
plot_fit_choice(
  input,
  post_predict = NULL,
  prior_predict = NULL,
  subject = NULL,
  quants = c(0.025, 0.975),
  functions = NULL,
  factors = NULL,
  signalFactor = "S",
  n_cores = 1,
  n_post = 50,
  layout = NA,
  style = c("prob", "cumulative", "roc"),
  to_plot = c("data", "posterior", "prior")[1:2],
  legendpos = "topright",
  posterior_args = list(),
  prior_args = list(),
  zROC = FALSE,
  qfun = qnorm,
  lim = NULL,
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

- signalFactor:

  The factor defining signal and noise classes for ROC plots.

- n_cores:

  Number of CPU cores to use if generating predictives from an `emc`
  object.

- n_post:

  Number of posterior draws to simulate if needed for predictives.

- layout:

  Numeric vector used in `par(mfrow=...)`; use `NA` for auto-layout.

- style:

  A string indicating which choice fit plot to draw: `"prob"`,
  `"cumulative"`, or `"roc"`.

- to_plot:

  Character vector: any of `"data"`, `"posterior"`, `"prior"`.

- legendpos:

  Character vector controlling the positions of the legends

- posterior_args:

  Optional list of graphical parameters for posterior lines/ribbons.

- prior_args:

  Optional list of graphical parameters for prior lines/ribbons.

- zROC:

  Boolean; if `TRUE`, plot a z-transformed ROC.

- qfun:

  Quantile function used when `zROC = TRUE`.

- lim:

  Optional common limits for ROC or zROC plots.

- ...:

  Other graphical parameters for the real data lines.

## Value

Invisibly returns a list with the plotted summaries for each source.

## Details

The default `style = "prob"` compares observed response probabilities to
posterior and/or prior predictive intervals. For ordered response
models, `style = "cumulative"` plots cumulative response probabilities.
For SDT-like two-signal designs, `style = "roc"` plots observed and
predictive ROC or zROC curves.

## Examples

``` r
# dmnl <- design(
#   Rlevels = c("left", "right", "up"),
#   factors = list(subjects = 1, S = c("left", "right", "up")),
#   formula = list(utility ~ lM),
#   contrasts = list(utility = list(lM = matrix(c(-1/2, 1/2), ncol = 1))),
#   matchfun = function(d) d$S == d$lR,
#   model = multinomial_logit
# )
# dat <- make_data(c(utility = 0, utility_lM1 = 2), dmnl, n_trials = 40)
# plot_fit_choice(dat, style = "prob", factors = "S")
```
