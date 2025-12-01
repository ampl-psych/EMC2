# Specify Priors for the Chosen Model

These values are entered manually by default but can be recycled from
another prior (given in the `update` argument).

## Usage

``` r
prior(
  design,
  type = NULL,
  group_design = NULL,
  update = NULL,
  do_ask = NULL,
  fill_default = TRUE,
  ...
)
```

## Arguments

- design:

  Design list for which a prior is constructed, typically the output of
  [`design()`](https://ampl-psych.github.io/EMC2/reference/design.md)

- type:

  Character. What type of group-level model you plan on using i.e.
  `diagonal`

- group_design:

  An `emc.group_design` object created with
  [`group_design()`](https://ampl-psych.github.io/EMC2/reference/group_design.md)

- update:

  Prior list from which to copy values

- do_ask:

  Character. For which parameter types or hyperparameters to ask for
  prior specification, i.e. `Sigma`, `mu` or `loadings` for factor
  models, but `theta_mu_mean` or `A` also works.

- fill_default:

  Boolean, If `TRUE` will fill all non-specified parameters, and
  parameters outside of `do_ask`, to default values

- ...:

  Either values to prefill, i.e. `theta_mu_mean = c(1:6)`, or additional
  arguments such as `n_factors = 2`

## Value

A prior list object

## Details

Where a value is not supplied, the user is prompted to enter numeric
values (or functions that evaluate to numbers).

To get the prior help use `prior_help(type)`. With `type` e.g.
'diagonal'.

## Examples

``` r
# First define a design for the model
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
# Then set up a prior using prior
p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
                     t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
psd <- c(v_Sleft=1,v_Sright=1,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
                     t0=.4,Z=1,sv=.4,SZ=1)
# Here we left the variance prior at default
prior_DDMaE <- prior(design_DDMaE,mu_mean=p_vector,mu_sd=psd)
# Also add a group-level variance prior:
pscale <- c(v_Sleft=.6,v_Sright=.6,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
                             t0=.2,Z=.5,sv=.4,SZ=.3)
df <- .4
prior_DDMaE <- prior(design_DDMaE,mu_mean=p_vector,mu_sd=psd, A = pscale, df = df)
# If we specify a new design
design_DDMat0E <- design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~E, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))
#> Parameter(s) st0 not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#>  [1] "v_Sleft"      "v_Sright"     "a"            "a_Eneutral"   "a_Eaccuracy" 
#>  [6] "t0"           "t0_Eneutral"  "t0_Eaccuracy" "Z"            "sv"          
#> [11] "SZ"          
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
#>         E t0 t0_Eneutral t0_Eaccuracy
#>     speed  1           0            0
#>   neutral  1           1            0
#>  accuracy  1           0            1
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
# We can easily update the prior
prior_DDMat0E <- prior(design_DDMat0E, update = prior_DDMaE)
```
