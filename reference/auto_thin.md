# Automatically Thin an emc Object

Uses the effective sample size of `selection` to determine how much to
optimally thin an emc object

## Usage

``` r
# S3 method for class 'emc'
auto_thin(emc, stage = "sample", selection = c("alpha", "mu"), ...)

auto_thin(emc, stage = "sample", selection = c("alpha", "mu"), ...)
```

## Arguments

- emc:

  an emc object.

- stage:

  A character string. Indicates from which sampling stage(s) to take the
  samples from (i.e. `preburn`, `burn`, `adapt`, `sample`)

- selection:

  Which parameter types (i.e. 'alpha' or 'mu' to consider when
  determining the effective sample size)

- ...:

  additional optional arguments
