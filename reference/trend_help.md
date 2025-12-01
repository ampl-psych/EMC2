# Get help information for trend kernels and bases

Get help information for trend kernels and bases

## Usage

``` r
trend_help(kernel = NULL, base = NULL, ...)
```

## Arguments

- kernel:

  Character string specifying the kernel type to get information about

- base:

  Character string specifying the base type to get information about

- ...:

  Additional arguments

## Value

Formatted trend information

## Examples

``` r
# Get information about exponential increasing kernel
trend_help(kernel = "exp_incr")
#> Description: 
#> Increasing exponential kernel: k = 1 - exp(-d_ei * c) 
#>  
#> Default transformations (in order): 
#> list(d_ei = "exp")
#>  
#> Available bases, first is the default: 
#> lin, exp_lin, centered 
#>  

# Get information about linear base
trend_help(base = "lin")
#> Description: 
#> Linear base: parameter + w * k 
#>  
#> Default transformations: 
#> list(w = "identity")
#>  

# Return available kernel and base types
trend_help()
#> Available kernels:
#>   custom: Custom C++ kernel: provided via register_trend().
#>   lin_decr: Decreasing linear kernel: k = -c
#>   lin_incr: Increasing linear kernel: k = c
#>   exp_decr: Decreasing exponential kernel: k = exp(-d_ed * c)
#>   exp_incr: Increasing exponential kernel: k = 1 - exp(-d_ei * c)
#>   pow_decr: Decreasing power kernel: k = (1 + c)^(-d_pd)
#>   pow_incr: Increasing power kernel: k = 1 - (1 + c)^(-d_pi)
#>   poly2: Quadratic polynomial: k = d1 * c + d2 * c^2
#>   poly3: Cubic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3
#>   poly4: Quartic polynomial: k = d1 * c + d2 * c^2 + d3 * c^3 + d4 * c^4
#>   delta: Standard delta rule kernel: k = q[i].
#>          Updates q[i] = q[i-1] + alpha * (c[i-1] - q[i-1]).
#>          Parameters: q0 (initial value), alpha (learning rate).
#>   delta2kernel: Dual kernel delta rule: k = q[i].
#>           Combines fast and slow learning rates
#>           and switches between them based on dSwitch.
#>           Parameters: q0 (initial value), alphaFast (fast learning rate),
#>           propSlow (alphaSlow = propSlow * alphaFast), dSwitch (switch threshold).
#>   delta2lr: Dual learning rate delta rule: k = q[i].
#>           Like the standard delta rule, but with separate
#>           learning rates for positive and negative prediction errors.
#>           Parameters: q0 (initial value), alphaPos (learning rate for positive PEs),
#>           alphaNeg (learning rate for negative PEs).
#> 
#> Available base types:
#>   lin: Linear base: parameter + w * k
#>   exp_lin: Exponential linear base: exp(parameter) + exp(w) * k
#>   centered: Centered mapping: parameter + w*(k - 0.5)
#>   add: Additive base: parameter + k
#>   identity: Identity base: k
#> 
#> Phase options:
#>   premap: Trend is applied before parameter mapping. This means the trend parameters
#>           are mapped first, then used to transform cognitive model parameters before 
#>           their mapping.
#>   pretransform: Trend is applied after parameter mapping but before transformations.
#>                 Cognitive model parameters are mapped first, then trend is applied, 
#>                 followed by transformations.
#>   posttransform: Trend is applied after both mapping and transformations.
#>                  Cognitive model parameters are mapped and transformed first, 
#>                  then trend is applied.
```
