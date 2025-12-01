# LNR Model of Forstmann Data (First 3 Subjects)

An emc object with a limited number of samples and subjects of the
Forstmann dataset. The object is a nested list of lenght three, each
list containing the MCMC samples of the respective chain. The MCMC
samples are stored in the samples element.

## Usage

``` r
samples_LNR
```

## Format

An emc object. An emc object is a list with a specific structure and
elements, as outlined below.

- data:

  A list of dataframes, one for each subject included

- par_names:

  A character vector containing the model parameter names

- n_pars:

  The number of parameters in the model

- n_subjects:

  The number of unique subject ID's in the data

- model:

  A list containing the model functions

- nuisance:

  A logical vector indicating which parameters are nuisance parameters

- subjects:

  A vector containing the unique subject ID's

- type:

  The type of model e.g., "standard" or "diagonal"

- prior:

  A list that holds the prior for `theta_mu` (the model parameters).
  Contains the mean (`theta_mu_mean`), covariance matrix
  (`theta_mu_var`), degrees of freedom (`v`), and scale (`A`) and
  inverse covariance matrix (`theta_mu_invar`)

- samples:

  A list with defined structure containing the samples, see the Samples
  Element section for more detail

- sampler_nuis:

  A sampler list for nuisance parameters (in this case there are none),
  similarly structured to the overall samples list of one of the MCMC
  chains.

## Source

<https://www.pnas.org/doi/10.1073/pnas.0805903105>

## Samples Element

The samples element of a emc object contains the different types of
samples estimated by EMC2. These include the three main types of samples
`theta_mu`, `theta_var` and `alpha` as well as a number of other items
which are detailed here.

- theta_mu:

  samples used for estimating the model parameters (group level), an
  array of size (n_pars x n_samples)

- theta_var:

  samples used for estimating the parameter covariance matrix, an array
  of size (n_pars x n_pars x n_samples)

- alpha:

  samples used for estimating the subject random effects, an array of
  size (n_pars x n_subjects x n_samples)

- stage:

  A vector containing what PMwG stage each sample was drawn in

- subj_ll:

  The winning particles log-likelihood for each subject and sample

- a_half:

  Mixing weights used during the Gibbs step when creating a new sample
  for the covariance matrix

- last_theta_var_inv:

  The inverse of the last samples covariance matrix

- idx:

  The index of the last sample drawn
