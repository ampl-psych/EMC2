# Factor diagram plot \#Makes a factor diagram plot. Heavily based on the fa.diagram function of the `psych` package.

Factor diagram plot \#Makes a factor diagram plot. Heavily based on the
fa.diagram function of the `psych` package.

## Usage

``` r
factor_diagram(
  emc = NULL,
  stage = "sample",
  loadings = NULL,
  standardize = TRUE,
  simple = FALSE,
  only_cred = TRUE,
  cut = 0,
  nice_names = NULL,
  factor_names = NULL,
  sort = TRUE,
  adj = 1,
  main = NULL,
  cex = NULL
)
```

## Arguments

- emc:

  An emc object

- stage:

  Character. The stage from which to take the samples

- loadings:

  An array of loadings. Can be alternatively supplied if emc is not
  supplied

- standardize:

  Boolean. Whether to standardize the loadings

- simple:

  Boolean. Whether the factor diagram should be simplified for visual
  clarity.

- only_cred:

  Boolean. Whether to only plot the credible loadings

- cut:

  Numeric. Mean loadings beneath this number will be excluded.

- nice_names:

  Character vector. Alternative names to give the parameters

- factor_names:

  Character vector. Names to give the different factors

- sort:

  Boolean. Whether to sort the paramaters before plotting for visual
  clarity.

- adj:

  Integer. Adjust to adjust loading values positions in the diagram if
  illegible.

- main:

  Character vector. Title of the plot

- cex:

  Integer. Font size
