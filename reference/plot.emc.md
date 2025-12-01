# Plot Function for emc Objects

Makes trace plots for model parameters.

## Usage

``` r
# S3 method for class 'emc'
plot(
  x,
  stage = "sample",
  selection = c("mu", "sigma2", "alpha"),
  layout = NA,
  ...
)
```

## Arguments

- x:

  An object of class `emc`

- stage:

  A character string indicating the sampling stage to be summarized. Can
  be `preburn`, `burn`, `adapt`, or `sample`.

- selection:

  A character vector indicating the parameter group(s). Defaults to
  `mu`, `sigma2`, and `alpha`.

- layout:

  A vector indicating which layout to use as in par(mfrow = layout). If
  NA, will automatically generate an appropriate layout.

- ...:

  Optional arguments that can be passed to `get_pars` or `plot.default`
  (see [`par()`](https://rdrr.io/r/graphics/par.html))

## Value

A trace/acf plot of the selected MCMC chains

## Details

If an emc object that has not been run with `fit` yet is supplied prior
plots will be returned.

## Examples

``` r
plot(samples_LNR)






# Or trace autocorrelation for the second subject:
plot(samples_LNR, subject = 2, selection = "alpha")


# Can also plot the trace of for example the group-level correlation:
plot(samples_LNR, selection = "correlation", col = c("green", "purple", "orange"), lwd = 2)



```
