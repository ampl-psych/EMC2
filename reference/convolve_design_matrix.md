# Convolve Events with HRF to Construct Design Matrices

This function convolves events with the HRF to construct design matrices
for fMRI analysis.

## Usage

``` r
convolve_design_matrix(
  timeseries,
  events,
  factors = NULL,
  contrasts = NULL,
  covariates = NULL,
  add_constant = TRUE,
  hrf_model = "glover",
  cell_coding = NULL,
  scale = TRUE,
  high_pass = TRUE,
  high_pass_model = "cosine",
  cut_off = 1e-12
)
```

## Arguments

- timeseries:

  A data frame containing fMRI time series data with columns 'subjects',
  'run', 'time', and at least one ROI column

- events:

  A data frame containing event information with required columns
  `subjects`, `run`, `onset`, `duration`, `event_type`, and `modulation`

- factors:

  A named list mapping factor names to event types

- contrasts:

  A named list of contrast matrices for each factor

- covariates:

  A character vector of event types to include as covariates

- add_constant:

  A boolean specifying whether a 1 should be included to the design
  matrix post convolution

- hrf_model:

  A character string specifying the HRF model to use ('glover', 'spm',
  'glover + derivative', or 'spm + derivative')

- cell_coding:

  A character vector of factor names to use cell coding for

- scale:

  A boolean indicating whether to scale the design matrix.

- high_pass:

  Logical indicating whether to apply high-pass filtering.
  Alternatively, specifying 'add' adds the regressors to the design
  matrix

- high_pass_model:

  Character indicating which type of high-pass filtering to apply
  ('cosine', 'poly')

- cut_off:

  A numeric value specifying the cutoff for the high-pass filter

## Value

A list containing the design matrices

## Examples

``` r
# Generate a simple example timeseries
ts <- data.frame(
  subjects = rep(1, 100),
  run = rep(1, 100),
  time = seq(0, 99),
  ROI1 = rnorm(100)
)

# Generate example events
events <- data.frame(
  subjects = rep(1, 4),
  run = rep(1, 4),
  onset = c(10, 30, 50, 70),
  duration = rep(0.5, 4),
  event_type = c("hard", "easy", "hard", "easy"),
  modulation = c(1, 1, 1, 1)
)

# Build design matrices
design_matrices <-  convolve_design_matrix(
  timeseries = ts,
  events = events,
  factors = list(difficulty = c("hard", "easy")),
  contrasts = list(difficulty = matrix(c(-1, 1)))
)
#>   event_type subjects run onset duration modulation   regressor
#> 1       hard        1   1    10      0.5          1 difficulty1
#> 2       easy        1   1    30      0.5         -1 difficulty1
#> Filtering out high_pass noise, make sure you also use high_pass_filter(<timeseries>)
```
