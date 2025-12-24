# Specify a Design and Model

This function combines information regarding the data, type of model,
and the model specification.

## Usage

``` r
design(
  formula = NULL,
  factors = NULL,
  Rlevels = NULL,
  model,
  data = NULL,
  contrasts = NULL,
  matchfun = NULL,
  constants = NULL,
  covariates = NULL,
  functions = NULL,
  report_p_vector = TRUE,
  custom_p_vector = NULL,
  trend = NULL,
  transform = NULL,
  bound = NULL,
  ...
)
```

## Arguments

- formula:

  A list. Contains the design formulae in the format
  `list(y ~ x, a ~ z)`.

- factors:

  A named list containing all the factor variables that span the design
  cells and that should be taken into account by the model. The name
  `subjects` must be used to indicate the participant factor variable,
  also in the data.

  Example:
  `list(subjects=levels(dat$subjects), condition=levels(dat$condition))`

- Rlevels:

  A character vector. Contains the response factor levels. Example:
  `c("right", "left")`

- model:

  A function, specifies the model type. Choose from the drift diffusion
  model ([`DDM()`](https://ampl-psych.github.io/EMC2/reference/DDM.md),
  `DDMt0natural()`), the log-normal race model
  ([`LNR()`](https://ampl-psych.github.io/EMC2/reference/LNR.md)), the
  linear ballistic model
  ([`LBA()`](https://ampl-psych.github.io/EMC2/reference/LBA.md)), the
  racing diffusion model
  ([`RDM()`](https://ampl-psych.github.io/EMC2/reference/RDM.md),
  `RDMt0natural()`), or define your own model functions.

- data:

  A data frame. `data` can be used to automatically detect `factors`,
  `Rlevels` and `covariates` in a dataset. The variable `R` needs to be
  a factor variable indicating the response variable. Any numeric column
  except `trials` and `rt` are treated as covariates, and all remaining
  factor variables are internally used in `factors`.

- contrasts:

  Optional. A named list specifying a design matrix. Example for
  supplying a customized design matrix:
  `list(lM = matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))))`

- matchfun:

  A function. Only needed for race models. Specifies whether a response
  was correct or not. Example: `function(d)d$S==d$lR` where lR refers to
  the latent response factor.

- constants:

  A named vector that sets constants. Any parameter in `sampled_pars`
  can be set constant.

- covariates:

  Names of numeric covariates.

- functions:

  List of functions to create new factors based on those in the factors
  argument. These new factors can then be used in `formula`.

- report_p_vector:

  Boolean. If TRUE (default), it returns the vector of parameters to be
  estimated.

- custom_p_vector:

  A character vector. If specified, a custom likelihood function can be
  supplied.

- trend:

  A trend list, as made by
  [`make_trend`](https://ampl-psych.github.io/EMC2/reference/make_trend.md)

- transform:

  A list with custom transformations to be applied to the parameters of
  the model, if the conventional transformations aren't desired. See
  [`DDM()`](https://ampl-psych.github.io/EMC2/reference/DDM.md) for an
  example of such transformations

- bound:

  A list with custom bounds to be applied to the parameters of the
  model, if the conventional bound aren't desired. see
  [`DDM()`](https://ampl-psych.github.io/EMC2/reference/DDM.md) for an
  example of such bounds. Bounds are used to set limits to the
  likelihood landscape that cannot reasonable be achieved with
  `transform`

- ...:

  Additional, optional arguments

## Value

A design list.

## Examples

``` r
# load example dataset
dat <- forstmann

# create a function that takes the latent response (lR) factor (d) and returns a logical
# defining the correct response for each stimulus. Here the match is simply
# such that the S factor equals the latent response factor
matchfun <- function(d)d$S==d$lR

# When working with lM and lR, it can be useful to design  an
# "average and difference" contrast matrix. For binary responses, it has a
# simple canonical form
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))

# Create a design for a linear ballistic accumulator model (LBA) that allows
# thresholds to be a function of E and lR. The final result is a 9 parameter model.
design_LBABE <- design(data = dat,model=LBA,matchfun=matchfun,
                            formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
                            contrasts=list(v=list(lM=ADmat)),
                            constants=c(sv=log(1)))
#> 
#>  Sampled Parameters: 
#> [1] "v"           "v_lMdiff"    "sv_lMTRUE"   "B"           "B_Eneutral" 
#> [6] "B_Eaccuracy" "B_lRright"   "A"           "t0"         
#> 
#>  Design Matrices: 
#> $v
#>     lM v v_lMdiff
#>   TRUE 1      0.5
#>  FALSE 1     -0.5
#> 
#> $sv
#>     lM sv sv_lMTRUE
#>   TRUE  1         1
#>  FALSE  1         0
#> 
#> $B
#>         E    lR B B_Eneutral B_Eaccuracy B_lRright
#>     speed  left 1          0           0         0
#>     speed right 1          0           0         1
#>   neutral  left 1          1           0         0
#>   neutral right 1          1           0         1
#>  accuracy  left 1          0           1         0
#>  accuracy right 1          0           1         1
#> 
#> $A
#>  A
#>  1
#> 
#> $t0
#>  t0
#>   1
#> 
```
