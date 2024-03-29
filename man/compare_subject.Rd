% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{compare_subject}
\alias{compare_subject}
\title{IC-based model weights for each participant in a list of samples objects}
\usage{
compare_subject(
  sList,
  filter = "sample",
  subfilter = 0,
  use_best_fit = TRUE,
  print_summary = TRUE,
  digits = 3
)
}
\arguments{
\item{sList}{List of samples objects}

\item{filter}{A string. Specifies which stage you want to plot.}

\item{subfilter}{An integer or vector. If integer it will exclude up until}

\item{use_best_fit}{Boolean, default TRUE use best of minD and Dmean in
calculation otherwise always use Dmean (see compare)}

\item{print_summary}{Boolean (default TRUE) print table of results}

\item{digits}{Integer, significant digits in printed table}
}
\value{
List of tables for each subject of effective number of parameters,
mean deviance, deviance of mean, DIC, BPIC (and optionally MD), and associated weights.
}
\description{
IC-based model weights for each participant in a list of samples objects
}
