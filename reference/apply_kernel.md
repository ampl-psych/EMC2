# Apply a kernel implied in an emc object

This function extracts the appropriate trend, and applies its implied
kernel using either the pure R implementation or the Rcpp version. When
`mode = "compare"`, the function checks whether the two implementations
produce identical output.

## Usage

``` r
apply_kernel(
  kernel_pars,
  emc,
  subject = 1,
  input_pars = NULL,
  trend_n = 1,
  mode = "Rcpp"
)
```

## Arguments

- kernel_pars:

  A named vector of kernel parameters on the natural scale. Use `NULL`
  for kernels that do not require parameters.

- emc:

  An `emc`.

- subject:

  Subject index for which to apply the kernel. Defaults to `1`.

- input_pars:

  Optional parameter matrix containing externally supplied parameter
  values (e.g., trend parameters). Only needed for custom kernels.

- trend_n:

  Integer specifying which trend to apply when multiple trends exist in
  the model. Defaults to `1`. A warning is issued if the model contains
  more than one trend.

- mode:

  Character string specifying which implementation to use:

  `"R"`

  :   Use the pure R implementation.

  `"Rcpp"`

  :   Use the Rcpp implementation (default).

## Value

Returns a kernel matrix produced by the corresponding implementation.

## Details

Applies the trend-specific kernel associated with an `emc` model to the
a subject's data and returns the resulting kernel matrix.
