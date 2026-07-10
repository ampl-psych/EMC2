# Add information criteria to emc object

Adds DIC, BPIC, and optionally MLL values as attributes to an emc
object. Can be useful to offload computational burden.

## Usage

``` r
add_ICs_MLL(
  emc,
  stage = "sample",
  filter = NULL,
  use_best_fit = TRUE,
  BayesFactor = TRUE,
  cores_for_props = 4,
  cores_per_prop = 1,
  ...
)
```

## Arguments

- emc:

  List of samples objects

- stage:

  A string. Specifies which stage the samples are to be taken from
  `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`

- filter:

  An integer or vector. If it's an integer, iterations up until the
  value set by `filter` will be excluded. If a vector is supplied, only
  the iterations in the vector will be considered.

- use_best_fit:

  Boolean, defaults to `TRUE`, uses the minimal or mean likelihood
  (whichever is better) in the calculation, otherwise always uses the
  mean likelihood.

- BayesFactor:

  Boolean, defaults to `TRUE`. Include marginal likelihoods as estimated
  using WARP-III bridge sampling. Usually takes a minute per model added
  to calculate

- cores_for_props:

  Integer, how many cores to use for the Bayes factor calculation, here
  4 is the default for the 4 different proposal densities to evaluate,
  only 1, 2 and 4 are sensible.

- cores_per_prop:

  Integer, how many cores to use for the Bayes factor calculation if you
  have more than 4 cores available. Cores used will be cores_for_props
  \* cores_per_prop. Best to prioritize cores_for_props being 4 or 2

- ...:

  Additional, optional arguments

## Value

An `emc` object with new attributes 'ICs' and 'MLL'

## Examples

``` r
# \donttest{
samples_with_ICs <- add_ICs_MLL(samples_LNR, cores_for_props = 1)
attr(samples_with_ICs, 'MLL')
#> [1] 285.3428
attr(samples_with_ICs, 'ICs')
#>        DIC       BPIC EffectiveN      meanD      Dmean       minD 
#> -620.95182 -606.17010   14.78171 -635.73353 -647.74821 -650.51524 

# Pre-computed MLLs and ICs are extracted when using compare():
compare(sList=list(samples_with_ICs, samples_LNR), cores_for_props=1)
#>     MD   wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
#> 1 -571 0.433 -621  0.5 -606   0.5         15  -636  -648 -651
#> 2 -571 0.567 -621  0.5 -606   0.5         15  -636  -648 -651
# Returns the same MD (barring noise), BPIC, DIC for both emc objects, as expected -
# but the first is extracted, the second computed in compare().
# }
```
