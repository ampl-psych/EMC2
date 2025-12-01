# Simulation-Based Calibration

Runs SBC for an EMC2 model and associated design. Returns normalized
rank (between 0 and 1) and prior samples. For hierarchical models the
group-level mean and the (implied) group-level (co-)variance are
returned. For non-hierarchical models only the subject-level parameters
rank is returned.

## Usage

``` r
run_sbc(
  design_in,
  prior_in,
  replicates = 250,
  trials = 100,
  n_subjects = 30,
  plot_data = FALSE,
  verbose = TRUE,
  fileName = NULL,
  ...
)
```

## Arguments

- design_in:

  An emc design list. The design of the model to be used in SBC

- prior_in:

  An emc prior list. The prior for the design to be used in SBC

- replicates:

  Integer. The number of samples to draw from the prior

- trials:

  Integer. The number of trials of the simulated data (per subject)

- n_subjects:

  Integer. Only used for hierarchical models. The number of subjects to
  be used in data generation of each replicate

- plot_data:

  Boolean. Whether to plot the data simulated (aggregated across
  subjects)

- verbose:

  Verbose. Whether to print progress related messages

- fileName:

  Character. Highly recommended, saves temporary results to the fileName

- ...:

  A list of optional additional arguments that can be passed to `fit`
  and `make_emc`

## Value

The ranks and prior samples. For hierarchical models also the
prior-generated subject-level parameters.
