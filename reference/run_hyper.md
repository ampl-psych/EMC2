# Run a Group-level Model.

Separate function for running only the group-level model. This can be
useful in a two-step analysis. Works similar in functionality to
make_emc, except also does the fitting and returns an emc object that
works with most posterior checking tests (but not the data
generation/posterior predictives).

## Usage

``` r
run_hyper(
  type = "standard",
  data,
  prior = NULL,
  iter = 1000,
  n_chains = 3,
  ...
)
```

## Arguments

- type:

  A string indicating whether to run a `standard` group-level,
  `blocked`, `diagonal`, `factor`, or `single` (i.e., non-hierarchical)
  model.

- data:

  A data frame, or a list of data frames. Needs to have the variable
  `subjects` as participant identifier.

- prior:

  an emc.prior object.

- iter:

  Number of MCMC samples to collect.

- n_chains:

  An integer. Specifies the number of mcmc chains to be run (has to be
  more than 1 to compute `rhat`).

- ...:

  Additional, optional arguments.

## Value

an emc object with only group-level samples
