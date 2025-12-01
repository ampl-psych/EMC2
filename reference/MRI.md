# GLM model for fMRI data

Creates a model specification for fMRI data using a normal distribution.
This model assumes that the observed BOLD signal follows a normal
distribution with a mean determined by the design matrix and betas, and
a standard deviation parameter for noise.

## Usage

``` r
MRI()
```

## Value

A list containing model specification

## Details

The model uses a normal distribution to model fMRI BOLD signals. Beta
parameters represent the effect sizes for different conditions, and the
sd parameter represents the standard deviation of the noise.

The log-likelihood function centers the predicted values by subtracting
the mean, which helps with model identifiability.

## Examples

``` r
# Create a normal MRI model specification
model_spec <- MRI()

# Access model parameters
model_spec$p_types
#> beta   sd 
#>    0    0 
```
