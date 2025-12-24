# Create a trend specification for model parameters

Create a trend specification for model parameters

## Usage

``` r
make_trend(
  par_names,
  cov_names = NULL,
  kernels,
  bases = NULL,
  shared = NULL,
  trend_pnames = NULL,
  phase = "premap",
  par_input = NULL,
  at = "lR",
  maps = NULL,
  custom_trend = NULL,
  ffill_na = NULL
)
```

## Arguments

- par_names:

  Character vector specifying which parameters to apply trend to

- cov_names:

  Character vector specifying which covariates to use for each trend

- kernels:

  Character vector specifying which kernel function to use for each
  trend

- bases:

  Optional character vector specifying which base function to use for
  each trend

- shared:

  Named list with entries the parameter names to be shared and the names
  the new names of the shared parameter.

- trend_pnames:

  Optional character vector specifying custom parameter names

- phase:

  Character vector (length 1 or `length(par_names)`) specifying the
  phase for each trend entry; one of "premap", "pretransform", or
  "posttransform". Defaults to "premap".

- par_input:

  Optional character vector(s) of parameter names to use as additional
  inputs for the trend

- at:

  If NULL, trend is applied to every row in the `dadm`. If a factor name
  (e.g., "lR"), trend is applied only to entries corresponding to the
  first level of that factor, and fed forward to the other levels of
  that factor. Defaults to "lR". For DDMs, `at` should be set to NULL.

- maps:

  List of functions that create matrices with which to multiply the
  covariates before applying the base. See details.

- custom_trend:

  A trend registered with `register_trend`

- ffill_na:

  Determines how missing covariate values are handled. If `TRUE`,
  missing values are forward-filled using the last known non-`NA` value
  after applying the kernel. If `FALSE`, trials with missing covariates
  contribute `0` instead. The default (NULL) is interpreted as `TRUE`
  for delta-rule models and `FALSE` otherwise.

## Value

A list containing the trend specifications for each parameter

## Details

The `maps` argument accepts one or more functions that translate
trial-level covariates into accumulator-specific predictors.

Example of a minimal map function:

    advantage_map <- function(dadm, cov_names) {
      lS      <- paste0('cov', ifelse(dadm$lR == 'left',  dadm$cov_left,  dadm$cov_right))
      lSother <- paste0('cov', ifelse(dadm$lR == 'right', dadm$cov_left,  dadm$cov_right))
      plus  <- sapply(cov_names, function(x) ifelse(lS      == x,  1, 0))
      minus <- sapply(cov_names, function(x) ifelse(lSother == x, -1, 0))
      plus + minus
    }

Multiple maps may be supplied, in which case the model will create a
separate base parameter for each map. See more examples below.

## Examples

``` r
# Put trend on B and v parameters
trend <- make_trend(
  par_names = c("B", "v"),
  cov_names = "strial",
  kernels = c("exp_incr", "poly3"),
  phase = "premap",
  shared = list(shrd = list("B.B0", "v.d1"))
)
get_trend_pnames(trend)
#> [1] "B.w"    "B.d_ei" "v.d2"   "v.d3"   "shrd"  


# Using covariate maps

# Covariate maps allow you to specify how trial-by-trial covariates influence
# model parameters for each accumulator. The example below uses a simple data
# frame with two trials. `cov_left` and `cov_right` specify which covariates
# correspond to the left and right accumulators on each trial. `S` indicates the
# correct response, and `cov1`–`cov4` contain the actual covariate values.

data <- data.frame(
  subjects = rep(1, 2),
  S        = c('left', 'right'),
  cov_left = c('1', '4'),
  cov_right= c('3', '2'),
  rt       = c(1.2, 0.8),
  R        = factor(c('left', 'right')),
  cov1     = c(1, NA),
  cov2     = c(NA, 1),
  cov3     = c(NA, 1),
  cov4     = c(1, 1)
)

# A covariate map function must take `dadm` and `cov_names` as inputs and return
# a matrix of size (nrow(dadm), length(cov_names)), coding how each covariate
# contributes to each accumulator.

advantage_map <- function(dadm, cov_names) {

  # Which stimulus does the accumulator correspond to on each trial?
  lS <- paste0('cov', ifelse(dadm$lR == 'left', dadm$cov_left, dadm$cov_right))

  # Which stimulus does the *other* accumulator correspond to?
  lSother <- paste0('cov', ifelse(dadm$lR == 'right', dadm$cov_left, dadm$cov_right))

  # Build indicator matrices
  map_plus1 <- sapply(cov_names, function(col) ifelse(lS     == col,  1, 0))
  map_minus1<- sapply(cov_names, function(col) ifelse(lSother == col, -1, 0))

  map_plus1 + map_minus1
}

# A covariate map function can be supplied to make_trend(), which creates the mapping
# specification for the model for each participant. Here, a single map ('differences') is provided.

trend <- make_trend(
  par_names = 'v',
  kernels   = 'delta',
  bases     = 'lin',
  cov_names = list(c('cov1', 'cov2', 'cov3', 'cov4')),
  maps      = list('differences' = advantage_map),
  at        = 'lR'
)

design_RDM <- design(
  model  = RDM,
  data   = data,
  formula= list(B ~ 1, v ~ 1, t0 ~ 1),
  trend  = trend
)
#> Intercept formula added for trend_pars: v.w_differences, v.q0, v.alpha
#> Parameter(s) A, s not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "B"               "v"               "t0"              "v.w_differences"
#> [5] "v.q0"            "v.alpha"        
#> 
#>  Design Matrices: 
#> $B
#>  B
#>  1
#> 
#> $v
#>  v
#>  1
#> 
#> $t0
#>  t0
#>   1
#> 
#> $v.w_differences
#>  v.w_differences
#>                1
#> 
#> $v.q0
#>  v.q0
#>     1
#> 
#> $v.alpha
#>  v.alpha
#>        1
#> 
#> $A
#>  A
#>  1
#> 
#> $s
#>  s
#>  1
#> 

emc <- make_emc(data, design_RDM, type = 'single')
#> Because the model contains a delta rule, data will not be compressed.
#> Processing data set 1

# The resulting covariate maps for each subject are attached to the `dadm`:
attr(emc[[1]]$data[[1]], 'covariate_maps')
#> $differences
#>      cov1 cov2 cov3 cov4
#> [1,]    1    0   -1    0
#> [2,]   -1    0    1    0
#> [3,]    0   -1    0    1
#> [4,]    0    1    0   -1
#> 
# And to confirm that this mapping is correct, compare with the corresponding `dadm`
emc[[1]]$data[[1]]
#>   subjects     S cov_left cov_right  rt     R cov1 cov2 cov3 cov4 trials    lR
#> 1        1  left        1         3 1.2  left    1   NA   NA    1      1  left
#> 2        1  left        1         3 1.2  left    1   NA   NA    1      1 right
#> 3        1 right        4         2 0.8 right   NA    1    1    1      2  left
#> 4        1 right        4         2 0.8 right   NA    1    1    1      2 right
#>   winner
#> 1   TRUE
#> 2  FALSE
#> 3  FALSE
#> 4   TRUE

# You can also provide multiple covariate maps. Each additional map introduces
# a separate base parameter. For example, the following `sum_map` is suitable
# for RL-ARD–type models:

sum_map <- function(dadm, cov_names) {
  # Which stimulus does the accumulator correspond to on each trial?
  lS <- paste0('cov', ifelse(dadm$lR == 'left', dadm$cov_left, dadm$cov_right))
  # Which stimulus does the *other* accumulator correspond to?
  lSother <- paste0('cov', ifelse(dadm$lR == 'right', dadm$cov_left, dadm$cov_right))

  # Indicator matrices (note: both are added rather than subtracted)
  map_this  <- sapply(cov_names, function(col) ifelse(lS      == col, 1, 0))
  map_other <- sapply(cov_names, function(col) ifelse(lSother == col, 1, 0))

  map_this + map_other
}

trend <- make_trend(
  par_names = 'v',
  kernels   = 'delta',
  bases     = 'lin',
  cov_names = list(c('cov1', 'cov2', 'cov3', 'cov4')),
  maps      = list('differences' = advantage_map,
                   'sums'        = sum_map),
  at        = 'lR'
)

design_RDM <- design(
  model  = RDM,
  data   = data,
  formula= list(B ~ 1, v ~ 1, t0 ~ 1),
  trend  = trend
)
#> Intercept formula added for trend_pars: v.w_differences, v.w_sums, v.q0, v.alpha
#> Parameter(s) A, s not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "B"               "v"               "t0"              "v.w_differences"
#> [5] "v.w_sums"        "v.q0"            "v.alpha"        
#> 
#>  Design Matrices: 
#> $B
#>  B
#>  1
#> 
#> $v
#>  v
#>  1
#> 
#> $t0
#>  t0
#>   1
#> 
#> $v.w_differences
#>  v.w_differences
#>                1
#> 
#> $v.w_sums
#>  v.w_sums
#>         1
#> 
#> $v.q0
#>  v.q0
#>     1
#> 
#> $v.alpha
#>  v.alpha
#>        1
#> 
#> $A
#>  A
#>  1
#> 
#> $s
#>  s
#>  1
#> 

emc <- make_emc(data, design_RDM, type = 'single')
#> Because the model contains a delta rule, data will not be compressed.
#> Processing data set 1

# Now the dadm contains two covariate maps, and the model includes two
# corresponding base parameters (e.g., v.w1 and v.w2):
attr(emc[[1]]$data[[1]], 'covariate_maps')
#> $differences
#>      cov1 cov2 cov3 cov4
#> [1,]    1    0   -1    0
#> [2,]   -1    0    1    0
#> [3,]    0   -1    0    1
#> [4,]    0    1    0   -1
#> 
#> $sums
#>      cov1 cov2 cov3 cov4
#> [1,]    1    0    1    0
#> [2,]    1    0    1    0
#> [3,]    0    1    0    1
#> [4,]    0    1    0    1
#> 

```
