# Merge Samples

Merges samples from all chains as one unlisted object.

## Usage

``` r
merge_chains(emc)
```

## Arguments

- emc:

  An emc object, commonly the output of
  [`fit()`](https://ampl-psych.github.io/EMC2/reference/fit.md)

## Value

An unlisted emc object with all chains merged

## Details

Note that all sampling stages are included in the merged output,
including iterations from the `preburn`, `burn`, and `adapt` stages.
`merge_chains(emc)$samples$stage` shows the corresponding sampling
stages.
