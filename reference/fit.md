# Model Estimation in EMC2

General purpose function to estimate models specified in EMC2.

## Usage

``` r
# S3 method for class 'emc'
fit(
  emc,
  stage = NULL,
  iter = 1000,
  stop_criteria = NULL,
  search_width = 1,
  step_size = 100,
  verbose = TRUE,
  fileName = NULL,
  particle_factor = 50,
  cores_per_chain = 1,
  cores_for_chains = length(emc),
  max_tries = 20,
  thin = FALSE,
  ...
)

fit(emc, ...)
```

## Arguments

- emc:

  An emc object created with `make_emc`, or a path to where the emc
  object is stored.

- stage:

  A string. Indicates which stage to start the run from, either
  `preburn`, `burn`, `adapt` or `sample`. If unspecified, it will run
  the subsequent stage (if there is one).

- iter:

  An integer. Indicates how many iterations to run in the sampling
  stage.

- stop_criteria:

  A list. Defines the stopping criteria and for which types of
  parameters these should hold. See the details and examples section.

- search_width:

  A double. Tunes target acceptance probability of the MCMC process.
  This fine-tunes the width of the search space to obtain the desired
  acceptance probability. 1 is the default width, increases lead to
  broader search.

- step_size:

  An integer. After each step, the stopping requirements as specified by
  `stop_criteria` are checked and proposal distributions are updated.
  Defaults to 100.

- verbose:

  Logical. Whether to print messages between each step with the current
  status regarding the `stop_criteria`.

- fileName:

  A string. If specified, will auto-save emc object at this location on
  every iteration.

- particle_factor:

  An integer. `particle_factor` multiplied by the square root of the
  number of sampled parameters determines the number of particles used.

- cores_per_chain:

  An integer. How many cores to use per chain. Parallelizes across
  participant calculations. Only available on Linux or Mac OS. For
  Windows, only parallelization across chains (`cores_for_chains`) is
  available.

- cores_for_chains:

  An integer. How many cores to use across chains. Defaults to the
  number of chains. The total number of cores used is equal to
  `cores_per_chain` \* `cores_for_chains`.

- max_tries:

  An integer. How many times should it try to meet the finish conditions
  as specified by `stop_criteria`? Defaults to 20. `max_tries` is
  ignored if the required number of iterations has not been reached yet.

- thin:

  A boolean. If `TRUE` will automatically thin the MCMC samples, closely
  matched to the ESS. Can also be set to a double, in which case 1/thin
  of the chain will be removed (does not have to be an integer).

- ...:

  Additional optional arguments

## Value

An emc object

## Details

`stop_criteria` is either a list of lists with names of the stages, or a
single list in which case its assumed to be for the sample `stage` (see
examples). The potential stop criteria to be set are:

`selection` (character vector): For which parameters the `stop_criteria`
should hold

`mean_gd` (numeric): The mean Gelman-Rubin diagnostic across all
parameters in the selection

`max_gd` (numeric): The max Gelman-Rubin diagnostic across all
parameters in the selection

`min_unique` (integer): The minimum number of unique samples in the MCMC
chains across all parameters in the selection

`min_es` (integer): The minimum number of effective samples across all
parameters in the selection

`omit_mpsrf` (Boolean): Whether to include the multivariate point-scale
reduction factor in the Gelman-Rubin diagnostic. Default is `FALSE`.

`iter` (integer): The number of MCMC samples to collect.

The estimation is performed using particle-metropolis within-Gibbs
sampling. For sampling details see:

Gunawan, D., Hawkins, G. E., Tran, M.-N., Kohn, R., & Brown, S. (2020).
New estimation approaches for the hierarchical linear ballistic
accumulator model. *Journal of Mathematical Psychology* ,96, 102368.
doi.org/10.1016/j.jmp.2020.102368

Stevenson, N., Donzallaz, M. C., Innes, R. J., Forstmann, B., Matzke,
D., & Heathcote, A. (2024). EMC2: An R Package for cognitive models of
choice. doi.org/10.31234/osf.io/2e4dq

## Examples

``` r
# \donttest{
# Define a design first
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR

# Drop most subjects for speed
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)

design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                     formula=list(m~lM,s~1,t0~1),
                     contrasts=list(m=list(lM=ADmat)))
#> 
#>  Sampled Parameters: 
#> [1] "m"     "m_lMd" "s"     "t0"   
#> 
#>  Design Matrices: 
#> $m
#>     lM m m_lMd
#>   TRUE 1   0.5
#>  FALSE 1  -0.5
#> 
#> $s
#>  s
#>  1
#> 
#> $t0
#>  t0
#>   1
#> 
# Before fit can be called, we first need to make an emc object
LNR_s <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2)
#> Processing data set 1
#> Likelihood speedup factor: 16 (104 unique trials)
# Run fit, here illustrating how to use stop_criteria (also for speed purposes)
LNR_s <- fit(LNR_s, cores_for_chains = 1, stop_criteria = list(
  preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
  sample = list(iter = 25, max_gd = 2)), verbose = FALSE, particle_factor = 30, step_size = 25)
# }
```
