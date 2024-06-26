% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{pairs_posterior}
\alias{pairs_posterior}
\title{Plot within-chain correlations}
\usage{
pairs_posterior(
  samplers,
  filter = "sample",
  thin = 1,
  subfilter = 0,
  mapped = FALSE,
  selection = c("alpha", "mu", "variance", "covariance", "correlation")[1],
  scale.subjects = TRUE,
  use_par = NULL,
  do_plot = TRUE,
  maxp = 500
)
}
\arguments{
\item{samplers}{An EMC2 samplers object}

\item{filter}{A character. Specifies which sampling stage should be plotted.
Defaults to the sampling stage \code{sample}.}

\item{thin}{An integer. Only iterations that are a multiple of \code{thin} are kept.}

\item{subfilter}{Integer or numeric vector. If an integer is supplied, iterations
up until that integer within the sampling stage \code{filter} are kept. If a vector is supplied, the iterations
within the range are kept.}

\item{mapped}{Boolean. If \code{TRUE}, plots the parameters as mapped back to the factor levels of the design,
if \code{FALSE}, the sampled parameters are plotted}

\item{selection}{A Character string. Indicates which parameter type to plot (\code{alpha}, \code{mu}, \code{variance}, \code{covariance}, \code{correlation}).
\code{LL} will plot the log-likelihood chains.}

\item{scale.subjects}{Boolean. To standardize each participant with \code{selection = "alpha"},
by subtracting the mean and divding by the standard deviation. This ensures the plot has every participant on the same scale.}

\item{use_par}{Character vector of names of parameters to plot (default \code{NULL} plots everything)}

\item{do_plot}{Boolean. Whether to plot the pairs plot, if \code{FALSE}, only the correlations
are returned.}

\item{maxp}{Integer for maximum number of iterations used (defaults to 500).
If number of samples in stage or selection exceeds maxp, a random subset will be taken of size maxp.}
}
\value{
Invisibly returns a matrix with the correlations between the parameters.
}
\description{
Plots within-chain parameter correlations (upper triangle) and corresponding
scatterplots (lower triangle) to visualize parameter sloppiness.
}
\details{
If \code{selection = alpha} the parameter chains are concatenated across participants,
(after standardizing if \code{scale.subjects = TRUE}) and then correlated.
}
\examples{
\dontrun{
# Plot the sloppiness for the individual-level subjects
pairs_posterior(samplers, selection = "alpha")

# We can also choose group-level parameters and subsets of the parameter space
pairs_posterior(samplers, use = c("v", "B", "t0"), selection = "variance")
}
}
