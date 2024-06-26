% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{make_samplers}
\alias{make_samplers}
\title{\code{make_samplers()}}
\usage{
make_samplers(
  data,
  design,
  model = NULL,
  type = "standard",
  n_chains = 3,
  compress = TRUE,
  rt_resolution = 0.02,
  prior_list = NULL,
  grouped_pars = NULL,
  par_groups = NULL,
  n_factors = NULL,
  constraintMat = NULL,
  ...
)
}
\arguments{
\item{data}{A data frame, or a list of data frames. Needs to have the variable \code{subjects} as participant identifier.}

\item{design}{A list with a pre-specified design, the output of \code{make_design()}.}

\item{model}{A model list. If none is supplied, the model specified in \code{make_design()} is used.}

\item{type}{A string indicating whether to run a \code{standard} group-level, \code{blocked}, \code{diagonal}, \code{factor}, or \code{single} (i.e., non-hierarchical) model.}

\item{n_chains}{An integer. Specifies the number of mcmc chains to be run (has to be more than 1 to compute \code{rhat}).}

\item{compress}{A Boolean, if \code{TRUE} (i.e., the default), the data is compressed to speed up likelihood calculations.}

\item{rt_resolution}{A double. Used for compression, response times will be binned based on this resolution.}

\item{prior_list}{A named list containing the prior. Default prior created if \code{NULL}. For the default priors, see \verb{?get_prior_\{type\}}.}

\item{grouped_pars}{An integer vector. Parameters on this location of the vector of parameters are treated as constant across subjects}

\item{par_groups}{A vector. Only to be specified with type \code{blocked}, e.g., \code{c(1,1,1,2,2)} means the covariances
of the first three and of the last two parameters are estimated as two separate blocks.}

\item{n_factors}{An integer. Only to be specified with type \code{factor}.}

\item{constraintMat}{A matrix of rows equal to the number of estimated parameters, and columns equal to the number of factors. Only to be specified with type factor.
If \code{NULL}, default settings as specified in Innes et al. 2022 will be used.}

\item{...}{Additional, optional arguments.}
}
\value{
A list
}
\description{
initializes the sampling process by combining the data, prior,
and model specification into a \code{samplers} object that is needed in \code{run_emc()}.
}
\examples{
dat <- forstmann

# function that takes the lR factor (named diff in the following function) and
# returns a logical defining the correct response for each stimulus. In this
# case the match is simply such that the S factor equals the latent response factor.
matchfun <- function(d)d$S==d$lR

# design an "average and difference" contrast matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))

# specify design
design_LBABE <- make_design(data = dat,model=LBA,matchfun=matchfun,
formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))

# specify priors
pmean <- c(v=1,v_lMdiff=1,sv_lMTRUE=log(.5), B=log(.5),B_Eneutral=log(1.5),
           B_Eaccuracy=log(2),B_lRright=0, A=log(0.25),t0=log(.2))
psd <- c(v=1,v_lMdiff=0.5,sv_lMTRUE=.5,
         B=0.3,B_Eneutral=0.3,B_Eaccuracy=0.3,B_lRright=0.3,A=0.4,t0=.5)
prior_LBABE <- make_prior(design_LBABE, type = 'standard',pmean=pmean,psd=psd)

# create samplers object
LBABE <- make_samplers(dat,design_LBABE,type="standard",  prior=prior_LBABE)

}
