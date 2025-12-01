# Register a custom C++ trend kernel

Compiles and registers a user-provided C++ function that maps per-trial
kernel parameters and inputs to a numeric vector. The C++ function must
have signature: NumericVector f(NumericMatrix trend_pars, NumericMatrix
input) and provide an exported pointer creator using EMC2_MAKE_PTR.

## Usage

``` r
register_trend(trend_parameters, file, transforms = NULL, base = "add")
```

## Arguments

- trend_parameters:

  Character vector of kernel parameter names (in order).

- file:

  Path to the C++ file implementing the custom kernel. The file should
  include EMC2/userfun.hpp and define a pointer creator (via
  EMC2_MAKE_PTR) that is exported to R.

- transforms:

  Optional named character vector or list mapping each custom kernel
  parameter name to a transform name (e.g., "identity", "exp", "pnorm").
  Length must match `trend_parameters`. If unnamed but the correct
  length, the order is assumed to match `trend_parameters`.

- base:

  Default base to use when creating trends with this custom kernel if no
  `bases` argument is supplied to `make_trend`. One of
  c("lin","exp_lin","centered","add","identity"). Default "add".

## Value

An object to pass to `make_trend(custom_trend=...)`, carrying the
pointer, parameter names, default base, and optional transform mapping.
