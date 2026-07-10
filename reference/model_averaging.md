# Model Averaging

Computes model weights and a Bayes factor by comparing two groups of
models based on their Information Criterion (IC) values. The function
works with either numeric vectors or data frames containing multiple IC
measures (e.g., MD, BPIC, DIC).

## Usage

``` r
model_averaging(IC_for, IC_against)
```

## Arguments

- IC_for:

  A numeric vector or the output of `compare`

- IC_against:

  A numeric vector or the output of `compare`

## Value

A `data.frame` with the following columns:

- `wFor`:

  The aggregated weight of the models in favor.

- `wAgainst`:

  The aggregated weight of the models against.

- `Factor`:

  The Bayes factor (ratio of `wFor` to `wAgainst`).

If `IC_for` is a data frame, a matrix with rows corresponding to each IC
measure is returned.

## Details

When provided with numeric vectors, it computes the weights for the two
groups by first converting the IC values into relative weights and then
normalizing them. When provided with a data frame, it assumes that the
data frame is the output of a call to `compare` and applies averaging to
each IC metric

## Examples

``` r
# First set up some example models (normally these would be alternative models)
samples_LNR2 <- subset(samples_LNR, length.out = 45)
samples_LNR3 <- subset(samples_LNR, length.out = 40)
samples_LNR4 <- subset(samples_LNR, length.out = 35)

# Run compare on them, BayesFactor = F is set for speed.
ICs <- compare(list(S1 = samples_LNR, S2 = samples_LNR2,
                    S3 = samples_LNR3, S4 = samples_LNR4), BayesFactor = FALSE)
#>     DIC  wDIC BPIC wBPIC EffectiveN meanD Dmean minD
#> S1 -621 0.193 -606 0.157         15  -636  -648 -651
#> S2 -621 0.208 -606 0.176         15  -636  -648 -651
#> S3 -621 0.213 -606 0.182         15  -636  -648 -651
#> S4 -622 0.387 -608 0.485         14  -636  -648 -650

# Model averaging can either be done with a vector of ICs:
model_averaging(ICs$BPIC[1:2], ICs$BPIC[2:4])
#>        wFor  wAgainst    Factor
#> 1 0.2834696 0.7165304 0.3956141

# Or the output of compare:
model_averaging(ICs[1:2,], ICs[3:4,])
#>           wFor  wAgainst    Factor
#> BPIC 0.3333656 0.6666344 0.5000725
#> DIC  0.4008657 0.5991343 0.6690749
```
