# The Log-Normal Race Model

Model file to estimate the Log-Normal Race Model (LNR) in EMC2.

## Usage

``` r
LNR()
```

## Value

A model list with all the necessary functions for EMC2 to sample

## Details

Model files are almost exclusively used in
[`design()`](https://ampl-psych.github.io/EMC2/reference/design.md).

Default values are used for all parameters that are not explicitly
listed in the `formula` argument of
[`design()`](https://ampl-psych.github.io/EMC2/reference/design.md).They
can also be accessed with `LNR()$p_types`.

|               |               |                   |             |             |                    |
|---------------|---------------|-------------------|-------------|-------------|--------------------|
| **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
| *m*           | \-            | \[-Inf, Inf\]     | 1           |             | Scale parameter    |
| *s*           | log           | \[0, Inf\]        | log(1)      |             | Shape parameter    |
| *t0*          | log           | \[0, Inf\]        | log(0)      |             | Non-decision time  |

Because the LNR is a race model, it has one accumulator per response
option. EMC2 automatically constructs a factor representing the
accumulators `lR` (i.e., the latent response) with level names taken
from the `R` column in the data.

In [`design()`](https://ampl-psych.github.io/EMC2/reference/design.md),
`matchfun` can be used to automatically create a latent match (`lM`)
factor with levels `FALSE` (i.e., the stimulus does not match the
accumulator) and `TRUE` (i.e., the stimulus does match the accumulator).
This is added internally and can also be used in the model formula,
typically for parameters related to the rate of accumulation (see the
example below).

Rouder, J. N., Province, J. M., Morey, R. D., Gomez, P., & Heathcote, A.
(2015). The lognormal race: A cognitive-process model of choice and
latency with desirable psychometric properties. *Psychometrika, 80*,
491-513. https://doi.org/10.1007/s11336-013-9396-3

## Examples

``` r
# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR
# We now construct our design, with v ~ lM and the contrast for lM the ADmat.
design_LNRmE <- design(data = forstmann,model=LNR,matchfun=matchfun,
                       formula=list(m~lM + E,s~1,t0~1),
                       contrasts=list(m=list(lM=ADmat)))
#> 
#>  Sampled Parameters: 
#> [1] "m"           "m_lMd"       "m_Eneutral"  "m_Eaccuracy" "s"          
#> [6] "t0"         
#> 
#>  Design Matrices: 
#> $m
#>     lM        E m m_lMd m_Eneutral m_Eaccuracy
#>   TRUE    speed 1   0.5          0           0
#>  FALSE    speed 1  -0.5          0           0
#>   TRUE  neutral 1   0.5          1           0
#>  FALSE  neutral 1  -0.5          1           0
#>   TRUE accuracy 1   0.5          0           1
#>  FALSE accuracy 1  -0.5          0           1
#> 
#> $s
#>  s
#>  1
#> 
#> $t0
#>  t0
#>   1
#> 
# For all parameters that are not defined in the formula, default values are assumed
# (see Table above).
```
