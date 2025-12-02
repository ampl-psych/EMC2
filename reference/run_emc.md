# Fine-Tuned Model Estimation

Although typically users will rely on `fit`, this function can be used
for more fine-tuned specification of estimation needs. The function will
throw an error if a stage is skipped, the stages have to be run in order
("preburn", "burn", "adapt", "sample"). More details can be found in the
`fit` help files
([`?fit`](https://ampl-psych.github.io/EMC2/reference/fit.md)).

## Usage

``` r
run_emc(
  emc,
  stage,
  stop_criteria,
  search_width = 1,
  step_size = 100,
  verbose = FALSE,
  verboseProgress = FALSE,
  fileName = NULL,
  particle_factor = 50,
  cores_per_chain = 1,
  cores_for_chains = length(emc),
  max_tries = 20,
  n_blocks = 1,
  thin = FALSE,
  trim = TRUE,
  r_cores = 1
)
```

## Arguments

- emc:

  An emc object

- stage:

  A string. Indicates which stage is to be run, either `preburn`,
  `burn`, `adapt` or `sample`

- stop_criteria:

  A list. Defines the stopping criteria and for which types of
  parameters these should hold. See
  [`?fit`](https://ampl-psych.github.io/EMC2/reference/fit.md).

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
  status regarding the stop_criteria.

- verboseProgress:

  Logical. Whether to print a progress bar within each step or not. Will
  print one progress bar for each chain and only if cores_for_chains =
  1.

- fileName:

  A string. If specified will autosave emc at this location on every
  iteration.

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
  number of chains. the total number of cores used is equal to
  `cores_per_chain` \* `cores_for_chains`.

- max_tries:

  An integer. How many times should it try to meet the finish conditions
  as specified by stop_criteria? Defaults to 20. max_tries is ignored if
  the required number of iterations has not been reached yet.

- n_blocks:

  An integer. Number of blocks. Will block the parameter chains such
  that they are updated in blocks. This can be helpful in extremely
  tough models with a large number of parameters.

- thin:

  A boolean. If `TRUE` will automatically thin the MCMC samples, closely
  matched to the ESS. Can also be set to a double, in which case 1/thin
  of the chain will be removed (does not have to be an integer).

- trim:

  A boolean. If `TRUE` will automatically remove redundant samples (i.e.
  from preburn, burn, adapt).

- r_cores:

  An integer for number of cores to use in R-based likelihood
  calculations, default 1.

## Value

An emc object

## Examples

``` r
# \donttest{
# First define a design
design_in <- design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1),
                           constants=c(s=log(1)))
#> Parameter(s) SZ, st0, sv not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "v_Sleft"     "v_Sright"    "a"           "a_Eneutral"  "a_Eaccuracy"
#> [6] "t0"          "Z"          
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
#> $SZ
#>  SZ
#>   1
#> 
#> $st0
#>  st0
#>    1
#> 
#> $sv
#>  sv
#>   1
#> 
# Then make the emc, we've omitted a prior here for brevity so default priors will be used.
emc <- make_emc(forstmann, design_in, compress = FALSE)
#> Processing data set 1

# Now for example we can specify that we only want to run the "preburn" phase
# for MCMC 10 iterations
# emc <- run_emc(emc, stage = "preburn", stop_criteria = list(iter = 10), cores_for_chains = 1)
# }
```
