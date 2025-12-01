# Estimating Marginal Likelihoods Using WARP-III Bridge Sampling

Uses bridge sampling that matches a proposal distribution to the first
three moments of the posterior distribution to get an accurate estimate
of the marginal likelihood. The marginal likelihood can be used for
computing Bayes factors and posterior model probabilities.

## Usage

``` r
run_bridge_sampling(
  emc,
  stage = "sample",
  filter = NULL,
  repetitions = 1,
  cores_for_props = 4,
  cores_per_prop = 1,
  both_splits = TRUE,
  ...
)
```

## Arguments

- emc:

  An emc object with a set of converged samples

- stage:

  A character indicating which stage to use, defaults to `sample`

- filter:

  An integer or vector. If integer, it will exclude up until that
  integer. If vector it will include everything in that range.

- repetitions:

  An integer. How many times to repeat the bridge sampling scheme. Can
  help get an estimate of stability of the estimate.

- cores_for_props:

  Integer. Warp-III evaluates the posterior over 4 different proposal
  densities. If you have the CPU, 4 cores will do this in parallel, 2 is
  also already helpful.

- cores_per_prop:

  Integer. Per density we can also parallelize across subjects. Eventual
  cores will be `cores_for_props` \* `cores_per_prop`. For efficiency
  users should prioritize cores_for_props being 4.

- both_splits:

  Boolean. Bridge sampling uses a proposal density and a target density.
  We can estimate the stability of our samples and therefore MLL
  estimate, by running 2 bridge sampling iterations The first one uses
  the first half of the samples as the proposal and the second half as
  the target, the second run uses the opposite. If this is is set to
  `FALSE`, it will only run bridge sampling once and it will instead do
  an odd-even iterations split to get a more reasonable estimate for
  just one run.

- ...:

  Additional, optional more in-depth hyperparameters

## Value

A vector of length repetitions which contains the marginal log
likelihood estimates per repetition

## Details

If not enough posterior samples were collected using
[`fit()`](https://ampl-psych.github.io/EMC2/reference/fit.md), bridge
sampling can be unstable. It is recommended to run
`run_bridge_sampling()` several times with the `repetitions` argument
and to examine how stable the results are.

It can be difficult to converge bridge sampling for exceptionally large
models, because of a large number of subjects (\> 100) and/or cognitive
model parameters.

For a practical introduction:

Gronau, Q. F., Heathcote, A., & Matzke, D. (2020). Computing Bayes
factors for evidence-accumulation models using Warp-III bridge sampling.
*Behavior research methods*, 52(2), 918-937.
doi.org/10.3758/s13428-019-01290-6

For mathematical background:

Meng, X.-L., & Wong, W. H. (1996). Simulating ratios of normalizing
constants via a simple identity: A theoretical exploration. *Statistica
Sinica*, 6, 831-860.
http://www3.stat.sinica.edu.tw/statistica/j6n4/j6n43/j6n43.htm

Meng, X.-L., & Schilling, S. (2002). Warp bridge sampling. *Journal of
Computational and Graphical Statistics*, 11(3), 552-586.
doi.org/10.1198/106186002457

## Examples

``` r
# \donttest{
# After `fit` has converged on a specific model
# We can take those samples and calculate the marginal log-likelihood for them
MLL <- run_bridge_sampling(samples_LNR, cores_for_props = 1, both_splits = FALSE)
# This will run on 2*4 cores (since 4 is the default for ``cores_for_props``)
# }
```
