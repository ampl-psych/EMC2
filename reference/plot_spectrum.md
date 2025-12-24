# Plot Empirical and Posterior Predictive Power Spectra

Computes and plots the empirical power spectrum, with optional overlay
of spectra from posterior predictive simulations. All customization of
the spectral estimator is done through `spectrum.args`, while plot
appearance is controlled via `...`.

## Usage

``` r
plot_spectrum(
  dat,
  pp = NULL,
  plot.log = TRUE,
  spectrum.args = list(),
  trial_duration = NULL,
  ...
)
```

## Arguments

- dat:

  A data frame containing empirical reaction time data, with at least
  columns `subjects` and `rt`.

- pp:

  Optional posterior predictive data in the same format as `dat`,
  including `subjects`, `postn`, and `trials`.

- plot.log:

  Logical. Whether to log-transform frequencies and power before
  plotting. This does not affect the call to
  [`spectrum`](https://rdrr.io/r/stats/spectrum.html): use
  `spectrum.args$list(log = TRUE)` to request log spectral density from
  the estimator itself.

- spectrum.args:

  A named list of arguments forwarded directly to
  [`spectrum`](https://rdrr.io/r/stats/spectrum.html) inside
  [`get_power_spectra`](https://ampl-psych.github.io/EMC2/reference/get_power_spectra.md).

- trial_duration:

  Optional duration of a trial in seconds. If supplied, the x-axis is
  labeled in human-readable time units. Otherwise the x-axis is in (log)
  frequencies of 1/trial.

- ...:

  Additional graphical parameters passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns a list containing the empirical spectrum and, if
posterior predictive data is supplied, the posterior predictive spectra
and their mean.
