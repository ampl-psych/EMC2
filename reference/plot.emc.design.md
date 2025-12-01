# Plot method for emc.design objects

Makes design illustration by plotting simulated data based on the design

## Usage

``` r
# S3 method for class 'emc.design'
plot(
  x,
  p_vector,
  data = NULL,
  factors = NULL,
  plot_factor = NULL,
  n_data_sim = 10,
  functions = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `emc.design` containing the design to plot

- p_vector:

  A named vector of parameter values to use for data generation

- data:

  Optional data frame to overlay on the design plot. If NULL, data will
  be simulated.

- factors:

  Character vector. Factors to use for varying parameters in the plot

- plot_factor:

  Optional character. Make separate plots for each level of this factor

- n_data_sim:

  Integer. If data is NULL, number of simulated datasets to generate for
  the plot. Default is 10.

- functions:

  Optional named list of functions that create additional columns in the
  data

- ...:

  Additional arguments passed to `make_design_plot`

## Value

No return value, called for side effect of plotting
