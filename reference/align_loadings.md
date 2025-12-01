# Reorder MCMC Samples of Factor Loadings

This function reorders MCMC samples of factor loadings to address the
label switching problem in Bayesian factor analysis. It implements a
parallelized version of the code and algorithm proposed by Papastamoulis
and Ntzoufras (2022)

## Usage

``` r
align_loadings(
  emc = NULL,
  lambda = NULL,
  n_cores = 1,
  verbose = TRUE,
  rotate_fun = NULL
)
```

## Arguments

- emc:

  an 'emc' object of type `infnt_factor`.

- lambda:

  Needs to be supplied if emc is not supplied. Array of factor loadings
  with dimensions p (variables) x q (factors) x n (MCMC iterations)

- n_cores:

  Number of cores for parallel processing

- verbose:

  Logical; whether to print progress information

- rotate_fun:

  A function that returns an orthogonally rotated factor loadings
  matrix. If NULL uses `varimax`

## Value

A list containing:

- lambda_reordered:

  Array of reordered loadings

- lambda_reordered_mcmc:

  Array of reordered loadings as MCMC object

- lambda_hat:

  Matrix of mean loadings after reordering

- v_vectors:

  Matrix of permutation vectors

- c_vectors:

  Matrix of sign-switching vectors

## References

Papastamoulis, P., & Ntzoufras, I. (2022). On the identifiability of
Bayesian factor analytic models. *Statistical Computing*, 32(2), 1-29.
doi: 10.1007/s11222-022-10084-4

## Examples

``` r
# This function works natively with emc objects, but also factor arrays:
# Simulate a small example with 5 variables, 2 factors, and 10 MCMC iterations
set.seed(123)
p <- 5  # Number of variables
q <- 2  # Number of factors
n <- 10 # Number of MCMC iterations

# Create random factor loadings with label switching
lambda <- array(0, dim = c(p, q, n))
for (i in 1:n) {
  # Generate base loadings
  base_loadings <- matrix(rnorm(p*q, 0, 0.5), p, q)
  base_loadings[1:3, 1] <- abs(base_loadings[1:3, 1]) + 0.5  # Strong loadings on factor 1
  base_loadings[4:5, 2] <- abs(base_loadings[4:5, 2]) + 0.5  # Strong loadings on factor 2

  # Randomly switch labels and signs
  if (runif(1) > 0.5) {
    # Switch factor order
    base_loadings <- base_loadings[, c(2, 1)]
  }
  if (runif(1) > 0.5) {
    # Switch sign of factor 1
    base_loadings[, 1] <- -base_loadings[, 1]
  }
  if (runif(1) > 0.5) {
    # Switch sign of factor 2
    base_loadings[, 2] <- -base_loadings[, 2]
  }

  lambda[,,i] <- base_loadings
}

# Align the loadings
result <- align_loadings(lambda = lambda, verbose = TRUE, n_cores = 1)
#> * iteration: 1 
#>    -----  objective function = 35.634 ----- 
#> 
#> * iteration: 2 
#>    -----  objective function = 21.451 ----- 
#> 
#> * iteration: 3 
#>    -----  objective function = 21.451 ----- 
#> 

# Examine the aligned loadings
print(result)
#> , , 1
#> 
#>            [,1]       [,2]
#> [1,]  0.2689155 1.12774886
#> [2,]  0.4272338 0.49891516
#> [3,]  1.4257358 0.06418314
#> [4,] -0.3765649 0.75551960
#> [5,] -0.2925759 0.66412558
#> 
#> , , 2
#> 
#>            [,1]       [,2]
#> [1,]  1.8080348 -0.1169574
#> [2,]  0.7331361 -0.2514198
#> [3,]  0.4786492 -0.5549789
#> [4,] -0.2503449  0.9816466
#> [5,]  0.7959448  0.3938997
#> 
#> , , 3
#> 
#>           [,1]       [,2]
#> [1,] 0.8153272 -0.6374533
#> [2,] 0.8607901  0.5587857
#> [3,] 1.3562341  0.1033301
#> [4,] 0.4701434  0.6113424
#> [5,] 0.1534571  0.9381935
#> 
#> , , 4
#> 
#>            [,1]        [,2]
#> [1,]  0.4581829 -0.70278658
#> [2,]  0.6929112  0.09575379
#> [3,]  2.1511153  0.01839815
#> [4,] -0.2764209  0.64124171
#> [5,]  0.3401396  1.07819467
#> 
#> , , 5
#> 
#>             [,1]        [,2]
#> [1,]  0.70117199 -0.02414626
#> [2,]  0.73295327 -0.03175625
#> [3,]  0.89952737  0.67170494
#> [4,] -0.03305273  0.61341156
#> [5,]  0.14435905  1.25632763
#> 
#> , , 6
#> 
#>            [,1]      [,2]
#> [1,]  0.1735129 1.4821360
#> [2,]  1.4728791 0.2815132
#> [3,]  0.3442323 0.5730832
#> [4,] -0.8623041 0.7475942
#> [5,] -0.8021993 0.4166953
#> 
#> , , 7
#> 
#>            [,1]        [,2]
#> [1,]  1.5633438 0.007026856
#> [2,]  0.6114172 0.667016952
#> [3,]  1.6443336 0.232574070
#> [4,]  0.2408506 1.194891313
#> [5,] -0.4781115 0.496004302
#> 
#> , , 8
#> 
#>             [,1]        [,2]
#> [1,]  0.32566050 -0.75455692
#> [2,]  0.73340200 -0.08634494
#> [3,]  0.53530423 -0.17334955
#> [4,] -0.02470959  0.66019441
#> [5,]  0.07100331  0.67582055
#> 
#> , , 9
#> 
#>             [,1]      [,2]
#> [1,] -0.02332887 0.9850349
#> [2,]  0.69697339 0.4366192
#> [3,]  1.23211360 0.3720901
#> [4,] -0.78206467 0.4580670
#> [5,]  0.42756135 1.1866916
#> 
#> , , 10
#> 
#>           [,1]        [,2]
#> [1,] 0.4762522 -0.62745314
#> [2,] 1.1839210 -0.18135240
#> [3,] 1.2118446  0.04758179
#> [4,] 0.6811685  0.37261935
#> [5,] 0.6309906  0.26789701
#> 
```
