# Create an AR(1) GLM model for fMRI data

This function creates a model specification for MRI data with an AR(1)
error structure. The model includes beta parameters for the design
matrix, a rho parameter for the autocorrelation, and a standard
deviation parameter for the noise.

## Usage

``` r
MRI_AR1()
```

## Value

A list containing the model specifications

## Details

The AR(1) model accounts for temporal autocorrelation in the data, where
each timepoint is correlated with the previous timepoint according to
the rho parameter.

## Examples

``` r
# Create an AR(1) GLM model for fMRI data
model_spec <- MRI_AR1()

# Access model parameters
model_spec$p_types
#>      beta       rho        sd 
#> 0.0000000 0.5003989 0.0000000 
```
