% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{run_emc}
\alias{run_emc}
\title{Model estimation in EMC2}
\usage{
run_emc(
  samplers,
  stage = NULL,
  iter = 1000,
  stop_criteria = NULL,
  p_accept = 0.8,
  step_size = 100,
  verbose = TRUE,
  verboseProgress = FALSE,
  fileName = NULL,
  particles = NULL,
  particle_factor = 50,
  cores_per_chain = 1,
  cores_for_chains = length(samplers),
  max_tries = 20,
  n_blocks = 1
)
}
\arguments{
\item{samplers}{An EMC2 samplers object created with \code{make_samplers},
or a path to where the samplers object is stored.}

\item{stage}{A string. Indicates which stage to start the run from, either \code{preburn}, \code{burn}, \code{adapt} or \code{sample}.
If unspecified, it will run the subsequent stage (if there is one).}

\item{iter}{An integer. Indicates how many iterations to run in the sampling stage.}

\item{stop_criteria}{A list. Defines the stopping criteria and for which types
of parameters these should hold. See the details and examples section.}

\item{p_accept}{A double. The target acceptance probability of the MCMC process.
This fine-tunes the width of the search space to obtain the desired acceptance probability. Defaults to .8}

\item{step_size}{An integer. After each step, the stopping requirements as specified
by \code{stop_criteria} are checked and proposal distributions are updated. Defaults to 100.}

\item{verbose}{Logical. Whether to print messages between each step with the current status regarding the \code{stop_criteria}.}

\item{verboseProgress}{Logical. Whether to print a progress bar within each step or not.
Will print one progress bar for each chain and only if \code{cores_for_chains = 1}.}

\item{fileName}{A string. If specified, will auto-save samplers object at this location on every iteration.}

\item{particles}{An integer. How many particles to use, default is \code{NULL} and
\code{particle_factor} is used instead. If specified, \code{particle_factor} is overwritten.}

\item{particle_factor}{An integer. \code{particle_factor} multiplied by the square
root of the number of sampled parameters determines the number of particles used.}

\item{cores_per_chain}{An integer. How many cores to use per chain. Parallelizes across
participant calculations. Only available on Linux or Mac OS. For Windows, only
parallelization across chains (\code{cores_for_chains}) is available.}

\item{cores_for_chains}{An integer. How many cores to use across chains.
Defaults to the number of chains. The total number of cores used is equal to \code{cores_per_chain} * \code{cores_for_chains}.}

\item{max_tries}{An integer. How many times should it try to meet the finish
conditions as specified by \code{stop_criteria}? Defaults to 20. \code{max_tries} is
ignored if the required number of iterations has not been reached yet.}

\item{n_blocks}{An integer. Number of blocks. Will block the parameter chains such that they are
updated in blocks. This can be helpful in extremely tough models with a large number of parameters.}
}
\value{
A list of samplers
}
\description{
General purpose function to estimate models specified in EMC2.
}
\details{
\code{stop_criteria} is either a list of lists with names of the stages,
or a single list in which case its assumed to be for the sample \code{stage} (see examples).
The potential stop criteria to be set are:

\code{selection} (character vector): For which parameters the \code{stop_criteria} should hold

\code{mean_gd} (numeric): The mean Gelman-Rubin diagnostic across all parameters in the selection

\code{max_gd} (numeric): The max Gelman-Rubin diagnostic across all parameters in the selection

\code{min_unique} (integer): The minimum number of unique samples in the MCMC chains across all parameters in the selection

\code{min_es} (integer): The minimum number of effective samples across all parameters in the selection

\code{omit_mpsrf} (Boolean): Whether to include the multivariate point-scale reduction factor in the Gelman-Rubin diagnostic. Default is \code{FALSE}.

\code{iter} (integer): The number of MCMC samples to collect.

The estimation is performed using particle-metropolis within-Gibbs sampling.
For sampling details see:

Gunawan, D., Hawkins, G. E., Tran, M.-N., Kohn, R., & Brown, S. (2020).
New estimation approaches for the hierarchical linear ballistic accumulator model.
\emph{Journal of Mathematical Psychology} ,96, 102368. doi.org/10.1016/j.jmp.2020.102368

Stevenson, N., Donzallaz, M. C., Innes, R. J., Forstmann, B., Matzke, D., & Heathcote, A. (2024).
EMC2: An R Package for cognitive models of choice. doi.org/10.31234/osf.io/2e4dq
}
\examples{
\dontrun{
# First define a design
design_DDMaE <- make_design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))
# Then make the samplers, we've omitted a prior here for brevity so default priors will be used.
samplers <- make_samplers(forstmann, design)

# With the samplers object we can start sampling by simply calling run_emc
samplers <- run_emc(samplers, fileName = "intermediate_save_location.RData")

# For particularly hard models it pays off to increase the ``particle_factor``
# and, although to a lesser extent, lower ``p_accept``.
samplers <- run_emc(samplers, particle_factor = 100, p_accept = .6)

# Example of how to use the stop_criteria:
samplers <- run_emc(samplers, stop_criteria = list(mean_gd = 1.1, max_gd = 1.5,
            selection = c('alpha', 'variance'), omit_mpsrf = TRUE, min_es = 1000).
# In this case the stop_criteria are set for the sample stage, which will be
# run until the mean_gd < 1.1, the max_gd < 1.5 (omitting the multivariate psrf)
# and the effective sample size > 1000,
# for both the individual-subject parameters ("alpha")
# and the group-level variance parameters.

# For the unspecified stages in the ``stop_criteria`` the default values
# are assumed which are found in Stevenson et al. 2024 <doi.org/10.31234/osf.io/2e4dq>

# Alternatively, you can also specify the stop_criteria for specific stages by creating a
# list of lists:
samplers <- run_emc(samplers, stop_criteria = list("burn" = list(mean_gd = 1.1, max_gd = 1.5,
            selection = c('alpha')), "adapt" = list(min_unique = 300)))
}
}
