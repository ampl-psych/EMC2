# Summary method for emc.design objects

Prints a summary of the design object, including sampled parameters and
design matrices. For continuous covariates just prints one row, instead
of all covariates.

## Usage

``` r
# S3 method for class 'emc.design'
summary(object, ...)
```

## Arguments

- object:

  An object of class `emc.design` containing the design to summarize

- ...:

  Additional arguments (not used)

## Value

Invisibly returns the design matrices
