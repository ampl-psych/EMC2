% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bridge_sampling.R
\name{run_bridge_sampling}
\alias{run_bridge_sampling}
\title{Estimating Marginal likelihoods using WARP-III bridge sampling}
\usage{
run_bridge_sampling(
  samplers,
  filter = "sample",
  subfilter = 0,
  repetitions = 1,
  cores_for_props = 4,
  cores_per_prop = 1,
  both_splits = T,
  maxiter = 5000,
  r0 = 1e-05,
  tol1 = 1e-10,
  tol2 = 1e-06,
  hyper_only = F
)
}
\arguments{
\item{samplers}{A list with a set of converged samplers}

\item{filter}{A character indicating which stage to use, default is sample}

\item{subfilter}{An integer or vector indicating whether to exclude up to (integer case), or which to exclude}

\item{repetitions}{An integer. How many times to repeat the bridge sampling scheme. Can help get an estimate of stability of the estimate.}

\item{cores_for_props}{Integer. Warp-III evaluates the posterior over 4 different proposal densities. If you have the power, 4 cores will do this in parallel, 2 also helps, 3 not really.}

\item{cores_per_prop}{Integer. Per density we can also parallelize across subjects. Eventual cores will be cores_for_props * cores_per_prop. Prioritize cores_for_props being 4.}

\item{both_splits}{Boolean. Bridge sampling uses a proposal density and a target density. We can estimate the stability of our samples and therefore MLL estimate, by running 2 bridge sampling iteratations.
The first one uses the first half of the samples as the proposal and the second half as the target, the second run uses the opposite. If this is is FALSE, it will only run bridge samplign once and
it will instead do a odd-even iterations split to get a more reasonable estimate for just one run.}

\item{maxiter}{Integer, how many iterations to run in a post-processing scheme, best left at default}

\item{r0}{Numeric. Hyperparameter in a post-processing scheme, best left at default}

\item{tol1}{Numeric. Hyperparameter in a post-processing scheme, best left at default}

\item{tol2}{Numeric. Hyperparameter in a post-processing scheme, best left at default}

\item{hyper_only}{Boolean. Experimental, to bridge sample only considering the group-level.}
}
\value{
A vector of length repetitions which contains the Marginal Log Likelihood estimates per repetition
}
\description{
Estimating Marginal likelihoods using WARP-III bridge sampling
}
