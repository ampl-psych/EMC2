# Plot Design

Makes design illustration by plotting simulated data based on the design

## Usage

``` r
# S3 method for class 'emc.design'
plot_design(
  x,
  data = NULL,
  factors = NULL,
  plot_factor = NULL,
  n_data_sim = 10,
  p_vector = NULL,
  functions = NULL,
  ...
)

# S3 method for class 'emc.prior'
plot_design(
  x,
  data = NULL,
  factors = NULL,
  plot_factor = NULL,
  n_data_sim = 10,
  p_vector = NULL,
  functions = NULL,
  ...
)

plot_design(
  x,
  data = NULL,
  factors = NULL,
  plot_factor = NULL,
  n_data_sim = 10,
  p_vector = NULL,
  functions = NULL,
  ...
)

# S3 method for class 'emc'
plot_design(
  x,
  data = NULL,
  factors = NULL,
  plot_factor = NULL,
  n_data_sim = 10,
  p_vector = NULL,
  functions = NULL,
  ...
)
```

## Arguments

- x:

  An `emc` or `emc.prior` object containing the design to plot

- data:

  Optional data to overlay on the design plot

- factors:

  Factors to use for varying parameters

- plot_factor:

  Optional. Make separate plots for each level of this factor

- n_data_sim:

  If data is provided, number of simulated datasets to generate for the
  plot. Default is 10.

- p_vector:

  Only needed when x is an `emc.design` object, which parameters to use
  for data generation.

- functions:

  A named list of functions that create additional columns in the data.

- ...:

  Additional arguments to pass to `make_design_plot`

## Value

No return value. Just plots the design
