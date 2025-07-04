% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factor_analysis.R
\name{factor_diagram}
\alias{factor_diagram}
\title{Factor diagram plot
#Makes a factor diagram plot. Heavily based on the fa.diagram function of the \code{psych} package.}
\usage{
factor_diagram(
  emc = NULL,
  stage = "sample",
  loadings = NULL,
  standardize = TRUE,
  simple = FALSE,
  only_cred = TRUE,
  cut = 0,
  nice_names = NULL,
  factor_names = NULL,
  sort = TRUE,
  adj = 1,
  main = NULL,
  cex = NULL
)
}
\arguments{
\item{emc}{An emc object}

\item{stage}{Character. The stage from which to take the samples}

\item{loadings}{An array of loadings. Can be alternatively supplied if emc is not supplied}

\item{standardize}{Boolean. Whether to standardize the loadings}

\item{simple}{Boolean. Whether the factor diagram should be simplified for visual clarity.}

\item{only_cred}{Boolean. Whether to only plot the credible loadings}

\item{cut}{Numeric. Mean loadings beneath this number will be excluded.}

\item{nice_names}{Character vector. Alternative names to give the parameters}

\item{factor_names}{Character vector. Names to give the different factors}

\item{sort}{Boolean. Whether to sort the paramaters before plotting for visual clarity.}

\item{adj}{Integer. Adjust to adjust loading values positions in the diagram if illegible.}

\item{main}{Character vector. Title of the plot}

\item{cex}{Integer. Font size}
}
\description{
Factor diagram plot
#Makes a factor diagram plot. Heavily based on the fa.diagram function of the \code{psych} package.
}
