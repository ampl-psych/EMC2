% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/design.R
\name{plot.emc.design}
\alias{plot.emc.design}
\title{Plot method for emc.design objects}
\usage{
\method{plot}{emc.design}(
  x,
  p_vector,
  data = NULL,
  factors = NULL,
  plot_factor = NULL,
  n_data_sim = 10,
  functions = NULL,
  ...
)
}
\arguments{
\item{x}{An object of class \code{emc.design} containing the design to plot}

\item{p_vector}{A named vector of parameter values to use for data generation}

\item{data}{Optional data frame to overlay on the design plot. If NULL, data will be simulated.}

\item{factors}{Character vector. Factors to use for varying parameters in the plot}

\item{plot_factor}{Optional character. Make separate plots for each level of this factor}

\item{n_data_sim}{Integer. If data is NULL, number of simulated datasets to generate for the plot. Default is 10.}

\item{functions}{Optional named list of functions that create additional columns in the data}

\item{...}{Additional arguments passed to \code{make_design_plot}}
}
\value{
No return value, called for side effect of plotting
}
\description{
Makes design illustration by plotting simulated data based on the design
}
