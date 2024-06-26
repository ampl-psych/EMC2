% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_standard.R
\name{get_prior_standard}
\alias{get_prior_standard}
\title{Prior specification or prior sampling for standard estimation.}
\usage{
get_prior_standard(
  prior = NULL,
  n_pars = NULL,
  sample = TRUE,
  N = 1e+05,
  type = "mu",
  design = NULL,
  map = FALSE
)
}
\arguments{
\item{prior}{A named list that can contain the prior mean (\code{theta_mu_mean}) and
variance (\code{theta_mu_var}) on the group-level mean, or the scale (\code{A}), or degrees of freedom (\code{v})
for the group-level variance-covariance matrix. For \code{NULL} entries, a default prior gets created.}

\item{n_pars}{Often inferred from the design, but if \code{design = NULL}, \code{n_pars}
will be used to determine the size of prior.}

\item{sample}{Boolean, defaults to \code{TRUE}, sample from the prior or simply return the prior specifications?}

\item{N}{How many samples to draw from the prior, the default is 1e5}

\item{type}{Character. If \code{sample = TRUE}, what prior to sample from. Options:
\code{"mu"}, \code{"variance"}, \code{"covariance"} \code{"full_var"}, \code{"alpha"}.}

\item{design}{The design obtained from \code{make_design()}, required when \code{map = TRUE}}

\item{map}{Boolean, defaults to \code{TRUE}. If \code{sample = TRUE}, the implied prior is sampled.
This includes back-transformations for naturally bounded parameters such as
the non-decision time and an inverse mapping from the design matrix back to the
cells of the design. If \code{FALSE}, the transformed, unmapped, parameters are used.
Note that \code{map} does not affect the prior used in the sampling process.}
}
\value{
A list with a single entry of type of samples from the prior (if sample = TRUE) or else a prior object
}
\description{
To get the default prior for a created design: get_prior_standard(design = design, sample = FALSE)
}
\details{
For details see Huang, A., & Wand, M. P. (2013). Simple marginally noninformative
prior distributions for covariance matrices. \emph{Bayesian Analysis}, 8, 439-452. https://doi.org/10.1214/13-BA815.

Note that if \code{sample = FALSE}, prior$theta_mu_invar (the inverse of the prior covariance matrix on the group-level mean) is returned,
which is only used for computational efficiency
}
\examples{
\dontrun{
# First define a design for the model
design_DDMaE <- make_design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))
# Now get the default prior
prior <- get_prior_standard(design = design_DDMaE, sample = FALSE)
# We can change values in the default prior or use make_prior
# Then we can get samples from this prior e.g.
samples <- get_prior_standard(prior = prior, design = design_DDMaE,
  sample = TRUE, type = "mu")
}
}
