# Posterior Credible Interval Tests

Modeled after `t.test`, returns the credible interval of the parameter
or test and what proportion of the posterior distribution (or the
difference in posterior distributions in case of a two sample test)
overlaps with mu. For a one sample test provide `x` and for two sample
also provide `y`. Note that for comparisons within one model, we
recommend using
[`hypothesis()`](https://ampl-psych.github.io/EMC2/reference/hypothesis.md)
if the priors were well chosen.

## Usage

``` r
# S3 method for class 'emc'
credible(
  x,
  x_name = NULL,
  x_fun = NULL,
  x_fun_name = "fun",
  selection = "mu",
  y = NULL,
  y_name = NULL,
  y_fun = NULL,
  y_fun_name = "fun",
  x_subject = NULL,
  y_subject = NULL,
  mu = 0,
  alternative = c("less", "greater")[1],
  probs = c(0.025, 0.5, 0.975),
  digits = 2,
  p_digits = 3,
  print_table = TRUE,
  ...
)

credible(x, ...)
```

## Arguments

- x:

  An emc object

- x_name:

  A character string. Name of the parameter to be tested for `x`

- x_fun:

  Function applied to the MCMC chains to create variable to be tested.

- x_fun_name:

  Name to give to quantity calculated by `x_fun`

- selection:

  A character string designating parameter type (e.g. `alpha` or
  `covariance`)

- y:

  A second emc object

- y_name:

  A character string. Name of the parameter to be tested for `y`

- y_fun:

  Function applied to the MCMC chains to create variable to be tested.

- y_fun_name:

  Name to give to quantity calculated by `y_fun`

- x_subject:

  Integer or name selecting a subject

- y_subject:

  Integer or name selecting a subject

- mu:

  Numeric. `NULL` value for single sample test if `y` is not supplied
  (default 0)

- alternative:

  `less` or `greater` determining direction of test probability

- probs:

  Vector defining quantiles to return.

- digits:

  Integer, significant digits for estimates in printed results

- p_digits:

  Integer, significant digits for probability in printed results

- print_table:

  Boolean (defaults to `TRUE`) for printing results table

- ...:

  Additional optional arguments that can be passed to `get_pars`

## Value

Invisible results table with no rounding.

## Examples

``` r
{
# Run a credible interval test (Bayesian ''t-test'')
credible(samples_LNR, x_name = "m")
# We can also compare between two sets of emc objects

# # Now without a ~ E
# design_null <- design(data = forstmann,model=DDM,
#                            formula =list(v~0+S,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
#                            constants=c(s=log(1)))
#
# null_model <- make_emc(forstmann, design_null)
# null_model <- fit(null_model)
# credible(x = null_model, x_name = "a", y = full_model, y_name = "a")
#
# # Or provide custom functions:
# credible(x = full_model, x_fun = function(d) d["a_Eaccuracy"] - d["a_Eneutral"])
}
#>           m mu
#> 2.5%  -1.22 NA
#> 50%   -0.97  0
#> 97.5% -0.77 NA
#> attr(,"less")
#> [1] 1
```
