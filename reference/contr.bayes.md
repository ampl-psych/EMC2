# Contrast Enforcing Equal Prior Variance on each Level

Typical contrasts impose different levels of marginal prior variance for
the different levels. This contrast can be used to ensure that each
level has equal marginal priors (Rouder, Morey, Speckman, & Province;
2012).

## Usage

``` r
contr.bayes(n)
```

## Arguments

- n:

  An integer. The number of items for which to create the contrast

## Value

A contrast matrix.

## Examples

``` r
{
design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.bayes),
formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
constants=c(s=log(1)))
}
#> Parameter(s) st0 not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "v"        "v_Sright" "a"        "a_E1"     "a_E2"     "t0"       "Z"       
#> [8] "sv"       "SZ"      
#> 
#>  Design Matrices: 
#> $v
#>      S v v_Sright
#>   left 1        0
#>  right 1        1
#> 
#> $a
#>         E a       a_E1       a_E2
#>     speed 1  0.0000000  0.8164966
#>   neutral 1 -0.7071068 -0.4082483
#>  accuracy 1  0.7071068 -0.4082483
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
