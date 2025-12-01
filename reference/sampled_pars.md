# Get Model Parameters from a Design

Makes a vector with zeroes, with names and length corresponding to the
model parameters of the design.

## Usage

``` r
sampled_pars(
  x,
  group_design = NULL,
  doMap = FALSE,
  add_da = FALSE,
  all_cells_dm = FALSE,
  data = NULL
)

# S3 method for class 'emc.design'
sampled_pars(
  x,
  group_design = NULL,
  doMap = FALSE,
  add_da = FALSE,
  all_cells_dm = FALSE,
  data = NULL
)

# S3 method for class 'emc.group_design'
sampled_pars(
  x,
  group_design = NULL,
  doMap = FALSE,
  add_da = FALSE,
  all_cells_dm = FALSE,
  data = NULL
)

# S3 method for class 'emc.prior'
sampled_pars(
  x,
  group_design = NULL,
  doMap = FALSE,
  add_da = FALSE,
  all_cells_dm = FALSE,
  data = NULL
)

# S3 method for class 'emc'
sampled_pars(
  x,
  group_design = NULL,
  doMap = FALSE,
  add_da = FALSE,
  all_cells_dm = FALSE,
  data = NULL
)
```

## Arguments

- x:

  an `emc.design` object made with
  [`design()`](https://ampl-psych.github.io/EMC2/reference/design.md) or
  an `emc` object.

- group_design:

  an `emc.group_design` object made with
  [`group_design()`](https://ampl-psych.github.io/EMC2/reference/group_design.md)

- doMap:

  logical. If `TRUE` will also include an attribute `map` with the
  design matrices that perform the mapping back to the design

- add_da:

  Boolean. Whether to include the relevant data columns in the map
  attribute

- all_cells_dm:

  Boolean. Whether to include all levels of a factor in the mapping
  attribute, even when one is dropped in the design

- data:

  A data frame to be included for accurate covariate mapping in
  summary.design

## Value

Named vector.

## Examples

``` r
# First define a design
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
# Then for this design get which cognitive model parameters are sampled:
sampled_pars(design_DDMaE)
#>     v_Sleft    v_Sright           a  a_Eneutral a_Eaccuracy          t0 
#>           0           0           0           0           0           0 
#>           Z          sv          SZ 
#>           0           0           0 
```
