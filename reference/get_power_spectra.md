# Compute Power Spectra With Optional Subject-Level Aggregation

Computes power spectral density estimates using
[`spectrum`](https://rdrr.io/r/stats/spectrum.html), optionally
aggregated across subjects or posterior predictive samples. All
arguments intended for the underlying spectral estimator should be
supplied through `spectrum.args`.

## Usage

``` r
get_power_spectra(data, by.postn = FALSE, spectrum.args = list())
```

## Arguments

- data:

  A data frame with reaction time data. Must contain `subjects` and
  `rt`, and for posterior predictive data optionally `postn` and
  `trials`.

- by.postn:

  Logical. If `TRUE`, compute a separate spectrum for each posterior
  predictive draw and each posterior sample index.

- spectrum.args:

  A named list of arguments passed directly to
  [`spectrum`](https://rdrr.io/r/stats/spectrum.html). These override
  the defaults internally used in this function. Useful for customizing
  smoothing spans, detrending, tapering, and so on. Defaults:
  `list(spans=c(3, 5), detrend=FALSE, demean=TRUE, log=FALSE, taper=0)`.
  By default, we run `spectrum` without `log`, and log-transform while
  plotting

## Value

Either a data frame with columns `freq` and `power`, or (if
`by.postn = TRUE`) a list with frequency vector and a matrix of spectra
across posterior samples.

## Details

The function organizes the data by subject (and optionally posterior
sample index), computes spectra individually, interpolates spectra to a
common frequency grid if needed, and averages them appropriately.
