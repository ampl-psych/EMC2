# Make SEM Diagram

Make SEM Diagram

## Usage

``` r
make_SEM_diagram(
  emc,
  plot_values = TRUE,
  cred_only = FALSE,
  par_names = NULL,
  cut = NULL,
  ...
)
```

## Arguments

- emc:

  an emc object

- plot_values:

  whether to plot the values or just the nodes/edges

- cred_only:

  whether to only plot credible values

- par_names:

  optional, if specified will overwrite the parameter names with
  user-defined names

- cut:

  optional. A numeric value factor loadings smaller than this value will
  be omitted.

- ...:

  optional additional arguments passed to render_graph from DiagrammeR

## Value

Invisibly returns a DiagrammeR graph object
