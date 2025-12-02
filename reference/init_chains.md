# Initialize Chains

Adds a set of start points to each chain. These start points are sampled
from a user-defined multivariate normal across subjects.

## Usage

``` r
init_chains(
  emc,
  start_mu = NULL,
  start_var = NULL,
  particles = 1000,
  cores_per_chain = 1,
  cores_for_chains = length(emc),
  ...
)
```

## Arguments

- emc:

  An emc object made by
  [`make_emc()`](https://ampl-psych.github.io/EMC2/reference/make_emc.md)

- start_mu:

  A vector. Mean of multivariate normal used in proposal distribution

- start_var:

  A matrix. Variance covariance matrix of multivariate normal used in
  proposal distribution. Smaller values will lead to less deviation
  around the mean.

- particles:

  An integer. Number of starting values

- cores_per_chain:

  An integer. How many cores to use per chain. Parallelizes across
  participant calculations.

- cores_for_chains:

  An integer. How many cores to use to parallelize across chains.
  Default is the number of chains.

- ...:

  optional additional arguments

## Value

An emc object

## Examples

``` r
# \donttest{
# Make a design and an emc object
design_DDMaE <- design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1),
                           constants=c(s=log(1)))
#> Parameter(s) SZ, Z, st0, sv not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "v_Sleft"     "v_Sright"    "a"           "a_Eneutral"  "a_Eaccuracy"
#> [6] "t0"         
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
#> $SZ
#>  SZ
#>   1
#> 
#> $Z
#>  Z
#>  1
#> 
#> $st0
#>  st0
#>    1
#> 
#> $sv
#>  sv
#>   1
#> 

DDMaE <- make_emc(forstmann, design_DDMaE, compress = FALSE)
#> Processing data set 1
# set up our mean starting points (same used across subjects).
mu <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
       t0=log(.2))
# Small variances to simulate start points from a tight range
var <- diag(0.05, length(mu))
# Initialize chains, 4 cores per chain, and parallelizing across our 3 chains as well
# so 4*3 cores used.
DDMaE <- init_chains(DDMaE, start_mu = mu, start_var = var,
                     cores_per_chain = 1, cores_for_chains = 1, particles = 3)
# Afterwards we can just use fit
# DDMaE <- fit(DDMaE, cores_per_chain = 4)
# }
```
