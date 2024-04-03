% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_blocked.R
\name{get_prior_blocked}
\alias{get_prior_blocked}
\title{Prior specification or prior sampling for blocked estimation.}
\usage{
get_prior_blocked(
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
\item{prior}{A named list that can contain the prior mean (theta_mu_mean) and
variance (theta_mu_var) on the group-level mean, or the scale (A), or degrees of freedom (df)
for the group-level variance-covariance matrix.
Default prior created for NULL entries.}

\item{n_pars}{Often inferred from the design, but if design = NULL will be used to determine size of prior.}

\item{sample}{Whether to sample from the prior or to simply return the prior. Default is TRUE,}

\item{N}{How many samples to draw from the prior, default 1e5}

\item{type}{character. If sample = TRUE, what prior to sample from. Options: "mu", "variance", "covariance" "full_var", "alpha".}

\item{design}{The design obtained from \code{make_design}, required when map = TRUE}

\item{map}{Boolean, default TRUE reverses malformation used by model to make
sampled parameters unbounded}
}
\value{
A list with a single entry of type of samples from the prior (if sample = TRUE) or else a prior object
}
\description{
Works analogous to \code{get_prior_standard}. Blocks of the covariance matrix to estimate are only considered in sampling.
For details see Huang and Wand, 2013.
To get the default prior for a created design: get_prior_diag(design = design, sample = FALSE)
}
\examples{
\dontrun{
# First define a design for the model
design_DDMaE <- make_design(data = forstmann,model=DDM,
                           formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                           constants=c(s=log(1)))
# Now get the default prior
prior <- get_prior_blocked(design = design_DDMaE, sample = FALSE)
# We can change values in the default prior or use make_prior
# Then we can get samples from this prior e.g.
samples <- get_prior_blocked(prior = prior, design = design_DDMaE,
  sample = TRUE, type = "mu")
}
}