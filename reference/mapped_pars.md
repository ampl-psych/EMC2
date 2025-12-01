# Parameter Mapping Back to the Design Factors

Maps parameters of the cognitive model back to the experimental design.
If p_vector is left unspecified will print a textual description of the
mapping. Otherwise the p_vector can be created using
[`sampled_pars()`](https://ampl-psych.github.io/EMC2/reference/sampled_pars.md).
The returned matrix shows whether/how parameters differ across the
experimental factors.

## Usage

``` r
mapped_pars(
  x,
  p_vector = NULL,
  model = NULL,
  digits = 3,
  remove_subjects = TRUE,
  covariates = NULL,
  ...
)

# S3 method for class 'emc.design'
mapped_pars(
  x,
  p_vector = NULL,
  model = NULL,
  digits = 3,
  remove_subjects = TRUE,
  covariates = NULL,
  ...
)

# S3 method for class 'emc.prior'
mapped_pars(
  x,
  p_vector = NULL,
  model = NULL,
  digits = 3,
  remove_subjects = TRUE,
  covariates = NULL,
  ...
)

# S3 method for class 'emc'
mapped_pars(
  x,
  p_vector = NULL,
  model = NULL,
  digits = 3,
  remove_subjects = TRUE,
  covariates = NULL,
  ...
)
```

## Arguments

- x:

  an `emc`, `emc.prior` or `emc.design` object

- p_vector:

  Optional. Specify parameter vector to get numeric mappings. Must be in
  the form of `sampled_pars(design)`

- model:

  Optional model type (if not already specified in `design`)

- digits:

  Integer. Will round the output parameter values to this many decimals

- remove_subjects:

  Boolean. Whether to include subjects as a factor in the design

- covariates:

  Covariates specified in the design can be included here.

- ...:

  optional arguments

## Value

Matrix with a column for each factor in the design and for each model
parameter type (`p_type`).

## Examples

``` r
# First define a design:
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
mapped_pars(design_DDMaE)
#> $v 
#>   S      
#>  left   : v_Sleft
#>  right  : v_Sright
#> 
#> $a 
#>   E         
#>  speed     : exp(a)
#>  neutral   : exp(a + a_Eneutral)
#>  accuracy  : exp(a + a_Eaccuracy)
#> 
# Then create a p_vector:
p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
          t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
# This will map the parameters of the p_vector back to the design
mapped_pars(design_DDMaE, p_vector)
#>          E     S  v   a  sv  t0 st0 s   Z  SZ    z   sz
#> 1    speed  left -2 1.0 0.5 0.2   0 1 0.5 0.5 0.50 0.50
#> 2  neutral  left -2 1.5 0.5 0.2   0 1 0.5 0.5 0.75 0.75
#> 3 accuracy  left -2 2.0 0.5 0.2   0 1 0.5 0.5 1.00 1.00
#> 4    speed right  2 1.0 0.5 0.2   0 1 0.5 0.5 0.50 0.50
#> 5  neutral right  2 1.5 0.5 0.2   0 1 0.5 0.5 0.75 0.75
#> 6 accuracy right  2 2.0 0.5 0.2   0 1 0.5 0.5 1.00 1.00
```
