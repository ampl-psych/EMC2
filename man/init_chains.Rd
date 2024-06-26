% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{init_chains}
\alias{init_chains}
\title{Initialize chains}
\usage{
init_chains(
  samplers,
  start_mu = NULL,
  start_var = NULL,
  particles = 1000,
  cores_per_chain = 1,
  cores_for_chains = length(samplers)
)
}
\arguments{
\item{samplers}{List of chains made by \code{make_samplers()}}

\item{start_mu}{A vector. Mean of multivariate normal used in proposal distribution}

\item{start_var}{A matrix. Variance covariance matrix of multivariate normal used in proposal distribution.
Smaller values will lead to less deviation around the mean.}

\item{particles}{An integer. Number of starting values}

\item{cores_per_chain}{An integer. How many cores to use per chain. Parallelizes across participant calculations.}

\item{cores_for_chains}{An integer. How many cores to use to parallelize across chains. Default is the number of chains.}
}
\value{
A samplers object
}
\description{
Adds a set of start points to each chain. These start points are sampled from a user-defined multivariate
normal across subjects.
}
\examples{
\dontrun{
# Make a design and a samplers object
design_DDMaE <- make_design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))

DDMaE <- make_samplers(forstmann, design_DDMaE)
# set up our mean starting points (same used across subjects).
mu <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
       t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
# Small variances to simulate start points from a tight range
var <- diag(0.05, length(mu))
# Initialize chains, 4 cores per chain, and parallelizing across our 3 chains as well
# so 4*3 cores used.
DDMaE <- init_chains(DDMaE, start_mu = p_vector, start_var = var, cores_per_chain = 4)
# Afterwards we can just use run_emc
DDMaE <- run_emc(DDMaE, cores_per_chain = 4)
}
}
