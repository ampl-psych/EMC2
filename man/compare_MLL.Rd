% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{compare_MLL}
\alias{compare_MLL}
\title{Calculate a table of model probabilities based for a list of samples objects
based on samples of marginal log-likelihood (MLL) added to these objects by
run_IS2. Probabilities estimated by a bootstrap ath picks a vector of MLLs,
one for each model in the list randomly with replacement nboot times,
calculates model probabilities and averages}
\usage{
compare_MLL(mll, nboot = 1e+05, digits = 2, print_summary = TRUE)
}
\arguments{
\item{mll}{List of samples objects with IS_samples attribute added by by run_IS2}

\item{nboot}{Integer number of bootstrap samples, the default (1e5) usually
gives stable results at 2 decimal places.}

\item{digits}{Integer, significant digits in printed table}

\item{print_summary}{Boolean (default TRUE) print table of results}
}
\value{
Vector of model probabilities with names from samples list.
}
\description{
Calculate a table of model probabilities based for a list of samples objects
based on samples of marginal log-likelihood (MLL) added to these objects by
run_IS2. Probabilities estimated by a bootstrap ath picks a vector of MLLs,
one for each model in the list randomly with replacement nboot times,
calculates model probabilities and averages
}
