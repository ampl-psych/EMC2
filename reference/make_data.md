# Simulate Data

Simulates data based on a model design and a parameter vector
(`p_vector`) by one of two methods:

1.  Creating a fully crossed and balanced design specified by the
    design, with number of trials per cell specified by the `n_trials`
    argument

2.  Using the design of a data frame supplied, which allows creation of
    unbalanced and other irregular designs, and replacing previous data
    with simulated data

## Usage

``` r
make_data(
  parameters,
  design = NULL,
  n_trials = NULL,
  data = NULL,
  expand = 1,
  staircase = NULL,
  functions = NULL,
  ...
)
```

## Arguments

- parameters:

  parameter vector used to simulate data. Can also be a matrix with one
  row per subject (with corresponding row names) or an emc object with
  sampled parameters (in which case posterior medians of `alpha` are
  used to simulate data)

- design:

  Design list created by
  [`design()`](https://ampl-psych.github.io/EMC2/reference/design.md)

- n_trials:

  Integer. If `data` is not supplied, number of trials to create per
  design cell

- data:

  Data frame. If supplied, the factors are taken from the data.
  Determines the number of trials per level of the design factors and
  can thus allow for unbalanced designs

- expand:

  Integer. Replicates the `data` (if supplied) expand times to increase
  number of trials per cell.

- staircase:

  Default NULL, used with stop-signal paradigm simulation to specify a
  staircase algorithm. If non-null and a list then passed through as is,
  if not it is assigned the default list structure:
  list(p=.25,SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf)

- functions:

  List of functions you want to apply to the data generation.

- ...:

  Additional optional arguments

## Value

A data frame with simulated data

## Details

To create data for multiple subjects see `?make_random_effects()`.

## Examples

``` r
# First create a design
design_DDMaE <- design(factors = list(S = c("left", "right"),
                                           E = c("SPD", "ACC"),
                                           subjects = 1:30),
                            Rlevels = c("left", "right"), model = DDM,
                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                            constants=c(s=log(1)))
#> Parameter(s) st0 not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "v_Sleft"  "v_Sright" "a"        "a_EACC"   "t0"       "Z"        "sv"      
#> [8] "SZ"      
#> 
#>  Design Matrices: 
#> $v
#>      S v_Sleft v_Sright
#>   left       1        0
#>  right       0        1
#> 
#> $a
#>    E a a_EACC
#>  SPD 1      0
#>  ACC 1      1
#> 
#> $t0
#>  t0
#>   1
#> 
#> $s
#>  s
#>  1
#> 
#> $Z
#>  Z
#>  1
#> 
#> $sv
#>  sv
#>   1
#> 
#> $SZ
#>  SZ
#>   1
#> 
#> $st0
#>  st0
#>    1
#> 
# Then create a p_vector:
parameters <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_EACC=log(2), t0=log(.2),
              Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))

# Now we can simulate data
data <- make_data(parameters, design_DDMaE, n_trials = 30)

# We can also simulate data based on a specific dataset
design_DDMaE <- design(data = forstmann,model=DDM,
                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                            constants=c(s=log(1)))
#> Parameter(s) st0 not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "v_Sleft"     "v_Sright"    "a"           "a_Eneutral"  "a_Eaccuracy"
#> [6] "t0"          "Z"           "sv"          "SZ"         
#> 
#>  Design Matrices: 
#> $v
#>      S v_Sleft v_Sright
#>   left       1        0
#>  right       0        1
#> 
#> $a
#>         E a a_Eneutral a_Eaccuracy
#>     speed 1          0           0
#>   neutral 1          1           0
#>  accuracy 1          0           1
#> 
#> $t0
#>  t0
#>   1
#> 
#> $s
#>  s
#>  1
#> 
#> $Z
#>  Z
#>  1
#> 
#> $sv
#>  sv
#>   1
#> 
#> $SZ
#>  SZ
#>   1
#> 
#> $st0
#>  st0
#>    1
#> 
parameters <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
              t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))

data <- make_data(parameters, design_DDMaE, data = forstmann)
```
