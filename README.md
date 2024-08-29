# EMC2: Extended Models of Choice 2:

The `R` package `EMC2` provides tools to perform Bayesian hierarchical analyses of the following cognitive models: Diffusion Decision Model (DDM), Linear Ballistic Accumulator Model (LBA), Racing Diffusion Model (RDM), and Lognormal Racing Model (LNR). Specifically, the package provides functionality for specifying individual model designs, estimating the models, examining convergence as well as model fit through posterior prediction methods. It also includes various plotting functions and relative model comparison methods such as Bayes factors. In addition, users can specify their own likelihood function and perform non-hierarchical estimation. The package uses particle metropolis Markov chain Monte Carlo sampling. For hierarchical models, it uses efficient Gibbs sampling at the population level and supports a variety of covariance structures, extending the work of Gunawan and colleagues (2020). 
## Installation

To install the R package, and its dependencies you can use: 

`remotes::install_github("ampl-psych/EMC2",dependencies=TRUE)`

## Workflow Overview

Pictured below are the four phases of an `EMC2`cognitive model analysis with associated functions (in `courier` font).  

&nbsp; 

<img width="800" alt="workflow-emc2" src="https://github.com/user-attachments/assets/ed96de6c-17b9-47f5-9153-cd19cff47d34">

&nbsp;

For details, please see: 

Stevenson, N., Donzallaz, M. C., Innes, R. J., Forstmann, B., Matzke, D., & Heathcote, A. (2024, January 30). EMC2: An R Package for cognitive models of choice. https://doi.org/10.31234/osf.io/2e4dq

## Bug Reports, Contributing, and Feature Requests

If you come across any bugs, or have ideas for extensions of `EMC2`, you can add them as an issue [here](https://github.com/ampl-psych/EMC2/issues). If you would like to contribute to the package's code, please submit a pull request.

## References

Stevenson, N., Donzallaz, M. C., Innes, R. J., Forstmann, B., Matzke, D., & Heathcote, A. (2024, January 30). EMC2: An R Package for cognitive models of choice. https://doi.org/10.31234/osf.io/2e4dq

Gunawan, D., Hawkins, G. E., Tran, M. N., Kohn, R., & Brown, S. D. (2020). New estimation approaches for the hierarchical Linear Ballistic Accumulator model. *Journal of Mathematical Psychology, 96,* 102368. https://doi.org/10.1016/j.jmp.2020.102368


<!-- badges: start -->

[![R-CMD-check](https://github.com/ampl-psych/EMC2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ampl-psych/EMC2/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->
