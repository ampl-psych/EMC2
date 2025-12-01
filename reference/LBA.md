# The Linear Ballistic Accumulator model

Model file to estimate the Linear Ballistic Accumulator (LBA) in EMC2.

## Usage

``` r
LBA()
```

## Value

A model list with all the necessary functions for EMC2 to sample

## Details

Model files are almost exclusively used in
[`design()`](https://ampl-psych.github.io/EMC2/reference/design.md).

Default values are used for all parameters that are not explicitly
listed in the `formula` argument of
[`design()`](https://ampl-psych.github.io/EMC2/reference/design.md).They
can also be accessed with `LBA()$p_types`.

|               |               |                   |             |               |                                                       |
|---------------|---------------|-------------------|-------------|---------------|-------------------------------------------------------|
| **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping**   | **Interpretation**                                    |
| *v*           | \-            | \[-Inf, Inf\]     | 1           |               | Mean evidence-accumulation rate                       |
| *A*           | log           | \[0, Inf\]        | log(0)      |               | Between-trial variation (range) in start point        |
| *B*           | log           | \[0, Inf\]        | log(1)      | *b* = *B*+*A* | Distance from *A* to *b* (response threshold)         |
| *t0*          | log           | \[0, Inf\]        | log(0)      |               | Non-decision time                                     |
| *sv*          | log           | \[0, Inf\]        | log(1)      |               | Between-trial variation in evidence-accumulation rate |

All parameters are estimated on the log scale, except for the drift rate
which is estimated on the real line.

Conventionally, `sv` is fixed to 1 to satisfy scaling constraints.

The *b* = *B* + *A* parameterization ensures that the response threshold
is always higher than the between trial variation in start point of the
drift rate.

Because the LBA is a race model, it has one accumulator per response
option. EMC2 automatically constructs a factor representing the
accumulators `lR` (i.e., the latent response) with level names taken
from the `R` column in the data.

The `lR` factor is mainly used to allow for response bias, analogous to
`Z` in the DDM. For example, in the LBA, response thresholds are
determined by the *B* parameters, so `B~lR` allows for different
thresholds for the accumulator corresponding to left and right stimuli
(e.g., a bias to respond left occurs if the left threshold is less than
the right threshold). For race models, the
[`design()`](https://ampl-psych.github.io/EMC2/reference/design.md)
argument `matchfun` can be provided, a function that takes the `lR`
factor (defined in the augmented data (d) in the following function) and
returns a logical defining the correct response. In the example below,
the match is simply such that the `S` factor equals the latent response
factor: `matchfun=function(d)d$S==d$lR`. Then `matchfun` is used to
automatically create a latent match (`lM`) factor with levels `FALSE`
(i.e., the stimulus does not match the accumulator) and `TRUE` (i.e.,
the stimulus does match the accumulator). This is added internally and
can also be used in model formula, typically for parameters related to
the rate of accumulation.

Brown, S. D., & Heathcote, A. (2008). The simplest complete model of
choice response time: Linear ballistic accumulation. *Cognitive
Psychology, 57*(3), 153-178.
https://doi.org/10.1016/j.cogpsych.2007.12.002

## Examples

``` r
# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR
# We now construct our design, with v ~ lM and the contrast for lM the ADmat.
design_LBABE <- design(data = forstmann,model=LBA,matchfun=matchfun,
                       formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
                       contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
#> 
#>  Sampled Parameters: 
#> [1] "v"           "v_lMd"       "sv_lMTRUE"   "B"           "B_Eneutral" 
#> [6] "B_Eaccuracy" "B_lRright"   "A"           "t0"         
#> 
#>  Design Matrices: 
#> $v
#>     lM v v_lMd
#>   TRUE 1   0.5
#>  FALSE 1  -0.5
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
# For all parameters that are not defined in the formula, default values are assumed
# (see Table above).
```
