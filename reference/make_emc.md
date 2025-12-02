# Make an emc Object

Creates an emc object by combining the data, prior, and model
specification into a `emc` object that is needed in
[`fit()`](https://ampl-psych.github.io/EMC2/reference/fit.md).

## Usage

``` r
make_emc(
  data,
  design,
  model = NULL,
  type = "standard",
  n_chains = 3,
  compress = TRUE,
  rt_resolution = 1/60,
  prior_list = NULL,
  group_design = NULL,
  par_groups = NULL,
  ...
)
```

## Arguments

- data:

  A data frame, or a list of data frames. Needs to have the variable
  `subjects` as participant identifier.

- design:

  A list with a pre-specified design, the output of
  [`design()`](https://ampl-psych.github.io/EMC2/reference/design.md).

- model:

  A model list. If none is supplied, the model specified in
  [`design()`](https://ampl-psych.github.io/EMC2/reference/design.md) is
  used.

- type:

  A string indicating whether to run a `standard` group-level,
  `blocked`, `diagonal`, `factor`, or `single` (i.e., non-hierarchical)
  model.

- n_chains:

  An integer. Specifies the number of mcmc chains to be run (has to be
  more than 1 to compute `rhat`).

- compress:

  A Boolean, if `TRUE` (i.e., the default), the data is compressed to
  speed up likelihood calculations.

- rt_resolution:

  A double. Used for compression, response times will be binned based on
  this resolution.

- prior_list:

  A named list containing the prior. Default prior created if `NULL`.
  For the default priors, see `?get_prior_{type}`.

- group_design:

  A design for group-level mappings, made using
  [`group_design()`](https://ampl-psych.github.io/EMC2/reference/group_design.md).

- par_groups:

  A vector. Indicates which parameters are allowed to correlate. Could
  either be a list of character vectors of covariance blocks. Or a
  numeric vector, e.g., `c(1,1,1,2,2)` means the covariances of the
  first three and of the last two parameters are estimated as two
  separate blocks.

- ...:

  Additional, optional arguments.

## Value

An uninitialized emc object

## Examples

``` r
dat <- forstmann

# function that takes the lR factor (named diff in the following function) and
# returns a logical defining the correct response for each stimulus. In this
# case the match is simply such that the S factor equals the latent response factor.
matchfun <- function(d)d$S==d$lR

# design an "average and difference" contrast matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))

# specify design
design_LBABE <- design(data = dat,model=LBA,matchfun=matchfun,
formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
#> 
#>  Sampled Parameters: 
#> [1] "v"           "v_lMdiff"    "sv_lMTRUE"   "B"           "B_Eneutral" 
#> [6] "B_Eaccuracy" "B_lRright"   "A"           "t0"         
#> 
#>  Design Matrices: 
#> $v
#>     lM v v_lMdiff
#>   TRUE 1      0.5
#>  FALSE 1     -0.5
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

# specify priors
pmean <- c(v=1,v_lMdiff=1,sv_lMTRUE=log(.5), B=log(.5),B_Eneutral=log(1.5),
           B_Eaccuracy=log(2),B_lRright=0, A=log(0.25),t0=log(.2))
psd <- c(v=1,v_lMdiff=0.5,sv_lMTRUE=.5,
         B=0.3,B_Eneutral=0.3,B_Eaccuracy=0.3,B_lRright=0.3,A=0.4,t0=.5)
prior_LBABE <- prior(design_LBABE, type = 'standard',pmean=pmean,psd=psd)

# create emc object
LBABE <- make_emc(dat,design_LBABE,type="standard",  prior=prior_LBABE,
                  compress = FALSE)
#> Processing data set 1
```
