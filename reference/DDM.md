# The Diffusion Decision Model

Model file to estimate the Diffusion Decision Model (DDM) in EMC2.

## Usage

``` r
DDM()
```

## Value

A model list with all the necessary functions for EMC2 to sample

## Details

Model files are almost exclusively used in
[`design()`](https://ampl-psych.github.io/EMC2/reference/design.md).

Default values are used for all parameters that are not explicitly
listed in the `formula` argument of
[`design()`](https://ampl-psych.github.io/EMC2/reference/design.md).They
can also be accessed with `DDM()$p_types`.

|               |               |                   |             |                                                 |                                                      |
|---------------|---------------|-------------------|-------------|-------------------------------------------------|------------------------------------------------------|
| **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping**                                     | **Interpretation**                                   |
| *v*           | \-            | \[-Inf, Inf\]     | 1           |                                                 | Mean evidence-accumulation rate (drift rate)         |
| *a*           | log           | \[0, Inf\]        | log(1)      |                                                 | Boundary separation                                  |
| *t0*          | log           | \[0, Inf\]        | log(0)      |                                                 | Non-decision time                                    |
| *s*           | log           | \[0, Inf\]        | log(1)      |                                                 | Within-trial standard deviation of drift rate        |
| *Z*           | probit        | \[0, 1\]          | qnorm(0.5)  | *z* = *Z* x *a*                                 | Relative start point (bias)                          |
| *SZ*          | probit        | \[0, 1\]          | qnorm(0)    | *sz* = 2 x *SZ* x min(*a* x *Z*, *a* x (1-*Z*)) | Relative between-trial variation in start point      |
| *sv*          | log           | \[0, Inf\]        | log(0)      |                                                 | Between-trial standard deviation of drift rate       |
| *st0*         | log           | \[0, Inf\]        | log(0)      |                                                 | Between-trial variation (range) in non-decision time |

`a`, `t0`, `sv`, `st0`, `s` are sampled on the log scale because these
parameters are strictly positive, `Z`, `SZ` and `DP` are sampled on the
probit scale because they should be strictly between 0 and 1.

`Z` is estimated as the ratio of bias to one boundary where 0.5 means no
bias. `DP` comprises the difference in non-decision time for each
response option.

Conventionally, `s` is fixed to 1 to satisfy scaling constraints.

See Ratcliff, R., & McKoon, G. (2008). The diffusion decision model:
theory and data for two-choice decision tasks. *Neural computation,
20*(4), 873-922. doi:10.1162/neco.2008.12-06-420.

## Examples

``` r
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
# For all parameters that are not defined in the formula, default values are assumed
# (see Table above).
```
