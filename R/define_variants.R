# Creates and initializes a sample storage object based on the specified type
# Parameters:
#   data: Input data
#   par_names: Names of parameters
#   type: Type of sampler ("standard", "single", "blocked", etc.)
#   iters: Number of iterations
#   stage: Sampling stage ("init" or other)
#   is_nuisance: Boolean vector indicating which parameters are nuisance parameters
sample_store <- function(data, par_names, type, iters = 1, stage = "init", is_nuisance = rep(F, length(par_names)), ...) {
  switch(type,
    "standard" = sample_store_standard(data, par_names, iters, stage, is_nuisance = is_nuisance, ...),
    "single" = sample_store_base(data, par_names, iters, stage, is_nuisance = is_nuisance, ...),
    "blocked" = sample_store_standard(data, par_names, iters, stage, is_nuisance = is_nuisance, ...),
    "diagonal" = sample_store_standard(data, par_names, iters, stage, is_nuisance = is_nuisance, ...),
    "factor" = sample_store_factor(data, par_names, iters, stage, is_nuisance = is_nuisance, ...),
    "infnt_factor" = sample_store_infnt_factor(data, par_names, iters, stage, is_nuisance = is_nuisance, ...),
    "SEM" = sample_store_SEM(data, par_names, iters, stage, is_nuisance = is_nuisance, ...),
    "diagonal-gamma" = sample_store_diag_gamma(data, par_names, iters, stage, is_nuisance = is_nuisance, ...),
    stop("Invalid type specified")
  )
}

# Adds prior and model information to the sampler object
# Parameters:
#   sampler: The sampler object to update
#   prior: Prior distribution information
#   type: Type of sampler
add_info <- function(sampler, prior, type, ...) {
  switch(type,
    "standard" = add_info_standard(sampler, prior, ...),
    "single" = add_info_single(sampler, prior, ...),
    "blocked" = add_info_blocked(sampler, prior, ...),
    "diagonal" = add_info_diag(sampler, prior, ...),
    "factor" = add_info_factor(sampler, prior, ...),
    "infnt_factor" = add_info_infnt_factor(sampler, prior, ...),
    "SEM" = add_info_SEM(sampler, prior, ...),
    "diagonal-gamma" = add_info_diag_gamma(sampler, prior, ...),
    stop("Invalid type specified")
  )
}

# Generates starting points for the MCMC chain
# Parameters:
#   sampler: The sampler object
#   start_mu: Initial mean values
#   start_var: Initial variance values
#   type: Type of sampler
get_startpoints <- function(sampler, start_mu, start_var, type, ...) {
  switch(type,
    "standard" = get_startpoints_standard(sampler, start_mu, start_var, ...),
    "single" = get_startpoints_single(sampler, start_mu, start_var, ...),
    "blocked" = get_startpoints_blocked(sampler, start_mu, start_var, ...),
    "diagonal" = get_startpoints_diag(sampler, start_mu, start_var, ...),
    "factor" = get_startpoints_factor(sampler, start_mu, start_var, ...),
    "infnt_factor" = get_startpoints_infnt_factor(sampler, start_mu, start_var, ...),
    "SEM" = get_startpoints_SEM(sampler, start_mu, start_var, ...),
    "diagonal-gamma" = get_startpoints_diag_gamma(sampler, start_mu, start_var, ...),
    stop("Invalid type specified")
  )
}

# Retrieves group-level parameters based on sampler type
# Parameters:
#   parameters: Current parameter values
#   s: subject
#   type: Type of sampler
get_group_level <- function(parameters, s, type, ...) {
  switch(type,
    "standard" = get_group_level_standard(parameters, s, ...),
    "single" = get_group_level_single(parameters, s, ...),
    "blocked" = get_group_level_standard(parameters, s, ...),
    "diagonal" = get_group_level_standard(parameters, s, ...),
    "factor" = get_group_level_standard(parameters, s, ...),
    "infnt_factor" = get_group_level_standard(parameters, s, ...),
    "SEM" = get_group_level_SEM(parameters, s, ...),
    "diagonal-gamma" = get_group_level_standard(parameters, s, ...),
    stop("Invalid type specified")
  )
}

# Updates the samples storage with new MCMC samples
# Parameters:
#   samples: Current samples storage
#   group_level: Group-level parameters
#   proposals: Proposed parameter values
#   j: Current iteration
#   n_pars: Number of parameters
#   type: Type of sampler
fill_samples <- function(samples, group_level, proposals, j = 1, n_pars, type, ...) {
  switch(type,
    "standard" = fill_samples_standard(samples, group_level, proposals, j, n_pars, ...),
    "single" = fill_samples_RE(samples, proposals, j, n_pars, ...),
    "blocked" = fill_samples_standard(samples, group_level, proposals, j, n_pars, ...),
    "diagonal" = fill_samples_standard(samples, group_level, proposals, j, n_pars, ...),
    "factor" = fill_samples_factor(samples, group_level, proposals, j, n_pars, ...),
    "infnt_factor" = fill_samples_infnt_factor(samples, group_level, proposals, j, n_pars, ...),
    "SEM" = fill_samples_SEM(samples, group_level, proposals, j, n_pars, ...),
    "diagonal-gamma" = fill_samples_diag_gamma(samples, group_level, proposals, j, n_pars, ...),
    stop("Invalid type specified")
  )
}

# Performs a Gibbs sampling step
# Parameters:
#   sampler: The sampler object
#   alpha: Acceptance probability
#   type: Type of sampler
gibbs_step <- function(sampler, alpha, type, ...) {
  switch(type,
    "standard" = gibbs_step_standard(sampler, alpha, ...),
    "single" = gibbs_step_single(sampler, alpha, ...),
    "blocked" = gibbs_step_blocked(sampler, alpha, ...),
    "diagonal" = gibbs_step_diag(sampler, alpha, ...),
    "factor" = gibbs_step_factor(sampler, alpha, ...),
    "infnt_factor" = gibbs_step_infnt_factor(sampler, alpha, ...),
    "SEM" = gibbs_step_SEM(sampler, alpha, ...),
    "diagonal-gamma" = gibbs_step_diag_gamma(sampler, alpha, ...),
    stop("Invalid type specified")
  )
}

# Calculates conditional distributions for parameters
# Parameters:
#   s: Sampler object
#   samples: Current samples
#   n_pars: Number of parameters
#   iteration: Current iteration
#   idx: Index for current parameter
#   type: Type of sampler
get_conditionals <- function(s, samples, n_pars, iteration = NULL, idx, type, ...) {
  switch(type,
    "standard" = get_conditionals_standard(s, samples, n_pars, iteration, idx, ...),
    "single" = get_conditionals_single(s, samples, n_pars, iteration, idx, ...),
    "blocked" = get_conditionals_blocked(s, samples, n_pars, iteration, idx, ...),
    "diagonal" = get_conditionals_diag(s, samples, n_pars, iteration, idx, ...),
    "factor" = get_conditionals_factor(s, samples, n_pars, iteration, idx, ...),
    "infnt_factor" = get_conditionals_infnt_factor(s, samples, n_pars, iteration, idx, ...),
    "SEM" = get_conditionals_SEM(s, samples, n_pars, iteration, idx, ...),
    "diagonal-gamma" = get_conditionals_diag(s, samples, n_pars, iteration, idx, ...),
    stop("Invalid type specified")
  )
}

# Adds additional information for bridge sampling
# Parameters:
#   info: Information object
#   samples: MCMC samples
#   type: Type of sampler
bridge_add_info <- function(info, samples, type, ...) {
  switch(type,
    "standard" = bridge_add_info_standard(info, samples, ...),
    "single" = bridge_add_info_single(info, samples, ...),
    "blocked" = bridge_add_info_blocked(info, samples, ...),
    "diagonal" = bridge_add_info_diag(info, samples, ...),
    "factor" = bridge_add_info_factor(info, samples, ...),
    "SEM" = bridge_add_info_SEM(info, samples, ...),
    stop("Invalid type specified")
  )
}

# Computes group-level parameters, prior densities, and Jacobian determinants for bridge sampling
# Parameters:
#   proposals_group: Group-level proposals
#   proposals_list: List of all proposals
#   info: Additional information
#   type: Type of sampler
bridge_group_and_prior_and_jac <- function(proposals_group, proposals_list, info, type, ...) {
  switch(type,
    "standard" = bridge_group_and_prior_and_jac_standard(proposals_group, proposals_list, info, ...),
    "single" = bridge_group_and_prior_and_jac_single(proposals_group, proposals_list, info, ...),
    "blocked" = bridge_group_and_prior_and_jac_blocked(proposals_group, proposals_list, info, ...),
    "diagonal" = bridge_group_and_prior_and_jac_diag(proposals_group, proposals_list, info, ...),
    "factor" = bridge_group_and_prior_and_jac_factor(proposals_group, proposals_list, info, ...),
    "SEM" = bridge_group_and_prior_and_jac_SEM(proposals_group, proposals_list, info, ...),
    stop("Invalid type specified")
  )
}

# Adds group-level samples to the complete samples collection
# Parameters:
#   all_samples: Complete collection of samples
#   samples: New samples to add
#   idx: Index for current group
#   type: Type of sampler
bridge_add_group <- function(all_samples, samples, idx, type, ...) {
  switch(type,
    "standard" = bridge_add_group_standard(all_samples, samples, idx, ...),
    "single" = bridge_add_group_single(all_samples, samples, idx, ...),
    "blocked" = bridge_add_group_blocked(all_samples, samples, idx, ...),
    "diagonal" = bridge_add_group_diag(all_samples, samples, idx, ...),
    "factor" = bridge_add_group_factor(all_samples, samples, idx, ...),
    "SEM" = bridge_add_group_SEM(all_samples, samples, idx, ...),
    stop("Invalid type specified")
  )
}

# Calculates group-level information criteria
# Parameters:
#   sampler: The sampler object
#   type: Type of sampler
group_IC <- function(sampler, type, filter = NULL, stage = "sample", ...) {
  switch(type,
    "standard" = group__IC_standard(sampler, filter = filter, stage = stage, ...),
    "single" = group__IC_single(sampler, filter = filter, stage = stage, ...),
    "blocked" = group__IC_standard(sampler, filter = filter, stage = stage, ...),
    "diagonal" = group__IC_standard(sampler, filter = filter, stage = stage, ...),
    "factor" = group__IC_standard(sampler, filter = filter, stage = stage, ...),
    "infnt_factor" = group__IC_standard(sampler, filter = filter, stage = stage, ...),
    "SEM" = group__IC_SEM(sampler, filter = filter, stage = stage, ...),
    "diagonal-gamma" = group__IC_standard(sampler, filter = filter, stage = stage, ...),
    stop("Invalid type specified")
  )
}

# Filters and processes MCMC samples
# Parameters:
#   samples: MCMC samples to filter
#   type: Type of sampler
filtered_samples <- function(samples, type, ...) {
  switch(type,
    "standard" = filtered_samples_standard(samples, ...),
    "single" = filtered_samples_single(samples, ...),
    "blocked" = filtered_samples_standard(samples, ...),
    "diagonal" = filtered_samples_standard(samples, ...),
    "factor" = filtered_samples_factor(samples, ...),
    "infnt_factor" = filtered_samples_infnt_factor(samples, ...),
    "SEM" = filtered_samples_SEM(samples, ...),
    "diagonal-gamma" = filtered_samples_standard(samples, ...),
    stop("Invalid type specified")
  )
}

