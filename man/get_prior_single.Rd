% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_single.R
\name{get_prior_single}
\alias{get_prior_single}
\title{Prior specification or prior sampling for single subject estimation.}
\usage{
get_prior_single(
  prior = NULL,
  n_pars = NULL,
  sample = TRUE,
  N = 1e+05,
  type = "alpha",
  design = NULL,
  map = FALSE
)
}
\arguments{
\item{prior}{A named list containing the prior mean (theta_mu_mean) and
variance (theta_mu_var). Default prior created if NULL}

\item{n_pars}{Argument used by the sampler, best left NULL. In user case inferred from the design}

\item{sample}{Whether to sample from the prior. Default is TRUE}

\item{N}{How many samples to draw from the prior, default 1e5}

\item{type}{FIX ME}

\item{design}{The design obtained from \code{make_design}, required when map = TRUE}

\item{map}{Boolean, default TRUE reverses malformation used by model to make
sampled parameters unbounded}
}
\value{
A list with single entry named "alpha" of samples from the prior (if sample = TRUE) or else a prior object
}
\description{
With this type of estimation, we assume that one, or multiple, subjects are
estimated without any hierarchical constraint. We need to specify a prior
with a multivariate normal from, by providing specifying prior$theta_mu_mean
a vector with an entry for each parameter, and a prior covariance matrix
prior$theta_mu_var, with default list(theta_mu_mean = rep(0, n_pars),
theta_mu_var = diag(rep(0, n_pars))).
}
