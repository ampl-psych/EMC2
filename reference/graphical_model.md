# Graphical Model

Draws a probabilistic graphical model (PGM) using circles, squares,
arrows, and plates for the standard EMC hierarchical model.

## Usage

``` r
graphical_model(emc)
```

## Arguments

- emc:

  an emc object of type `standard`, `blocked`, or `diagonal`

## Value

Invisibly returns a DiagrammeR `grViz` graph object
