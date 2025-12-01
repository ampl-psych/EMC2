# Contrast Enforcing Decreasing Estimates

Each level will be estimated as a reduction from the previous level

## Usage

``` r
contr.decreasing(n)
```

## Arguments

- n:

  an integer. The number of items for which to create the contrast.

## Value

a contrast matrix.

## Examples

``` r
{
design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.decreasing),
formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
constants=c(s=log(1)))
}
#> Parameter(s) st0 not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "v"        "v_Sright" "a"        "a_E2"     "a_E3"     "t0"       "Z"       
#> [8] "sv"       "SZ"      
#> 
#>  Design Matrices: 
#> $v
#>      S v v_Sright
#>   left 1        0
#>  right 1        1
#> 
#> $a
#>         E a a_E2 a_E3
#>     speed 1    1    1
#>   neutral 1    1    0
#>  accuracy 1    0    0
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
```
