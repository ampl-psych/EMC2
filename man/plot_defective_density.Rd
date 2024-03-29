% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_defective_density}
\alias{plot_defective_density}
\title{Plot defective densities for each subject and cell.}
\usage{
plot_defective_density(
  data,
  subject = NULL,
  factors = NULL,
  layout = NULL,
  mfcol = FALSE,
  xlim = NULL,
  bw = "nrd0",
  adjust = 1,
  correct_fun = NULL,
  rt = "top",
  accuracy = "topright"
)
}
\arguments{
\item{data}{data frame with at least subjects (subjects factor) R (response factor)
and rt (response time) columns,and optionally other factor columns with any name except
subjects, R, rt or trials.}

\item{subject}{Integer or string selecting a subject (default NULL = all).}

\item{factors}{character vector of factor names in design (default NULL = all).}

\item{layout}{2-vector specifying par(mfrow) or par(mfcol) (default NULL use current).}

\item{mfcol}{Boolean, default TRUE use mfcol else mfrow.}

\item{xlim}{x-axis limit for all cells (default NULL = scale per cell).}

\item{bw}{number or string bandwidth for density (default "nrd0").}

\item{adjust}{density function bandwidth adjust parameter.}

\item{correct_fun}{function scoring accuracy using columns in data.}

\item{rt}{legend function position string for mean RT (default "top)}

\item{accuracy}{legend function position string for accuracy (default "topright")}
}
\value{
Invisibly if correct_fun specified a subject accuracy vector
}
\description{
Each panel contains a set of densities (i.e., densities for each possible
possible response) that are defective (i.e., have areas potentially less
than 1, where for all responses the area sums to 1).
}
