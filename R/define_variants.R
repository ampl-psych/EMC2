sample_store <- function(data, par_names, type, iters = 1, stage = "init", ...) {
  switch(type,
    "standard" = sample_store_standard(data, par_names, iters, stage, ...),
    "single" = sample_store_base(data, par_names, iters, stage, ...),
    "blocked" = sample_store_standard(data, par_names, iters, stage, ...),
    "diagonal" = sample_store_standard(data, par_names, iters, stage, ...),
    "infnt_factor" = sample_store_infnt_factor(data, par_names, iters, stage, ...),
    "SEM" = sample_store_SEM(data, par_names, iters, stage, ...),
    "diagonal-gamma" = sample_store_diag_gamma(data, par_names, iters, stage, ...),
    stop("Invalid type specified")
  )
}

add_info <- function(sampler, prior, type, ...) {
  switch(type,
    "standard" = add_info_standard(sampler, prior, ...),
    "single" = add_info_single(sampler, prior, ...),
    "infnt_factor" = add_info_infnt_factor(sampler, prior, ...),
    "SEM" = add_info_SEM(sampler, prior, ...),
    "diagonal-gamma" = add_info_diag_gamma(sampler, prior, ...),
    stop("Invalid type specified")
  )
}

bridge_add_info <- function(info, samples, type, ...) {
  switch(type,
    "standard" = bridge_add_info_standard(info, samples, ...),
    "single" = bridge_add_info_single(info, samples, ...),
    "SEM" = bridge_add_info_SEM(info, samples, ...),
    stop("Invalid type specified")
  )
}

bridge_group_and_prior_and_jac <- function(proposals_group, proposals_list, info, type, ...) {
  switch(type,
    "standard" = bridge_group_and_prior_and_jac_standard(proposals_group, proposals_list, info, ...),
    "single" = bridge_group_and_prior_and_jac_single(proposals_group, proposals_list, info, ...),
    "SEM" = bridge_group_and_prior_and_jac_SEM(proposals_group, proposals_list, info, ...),
    stop("Invalid type specified")
  )
}

bridge_add_group <- function(all_samples, samples, idx, type, ...) {
  switch(type,
    "standard" = bridge_add_group_standard(all_samples, samples, idx, ...),
    "single" = bridge_add_group_single(all_samples, samples, idx, ...),
    "SEM" = bridge_add_group_SEM(all_samples, samples, idx, ...),
    stop("Invalid type specified")
  )
}

group_IC <- function(sampler, type, filter = NULL, stage = "sample", ...) {
  switch(type,
    "standard" = group__IC_standard(sampler, filter = filter, stage = stage, ...),
    "single" = group__IC_single(sampler, filter = filter, stage = stage, ...),
    "infnt_factor" = group__IC_standard(sampler, filter = filter, stage = stage, ...),
    "SEM" = group__IC_SEM(sampler, filter = filter, stage = stage, ...),
    "diagonal-gamma" = group__IC_standard(sampler, filter = filter, stage = stage, ...),
    stop("Invalid type specified")
  )
}

filtered_samples <- function(samples, type, ...) {
  switch(type,
    "standard" = filtered_samples_standard(samples, ...),
    "single" = filtered_samples_single(samples, ...),
    "infnt_factor" = filtered_samples_infnt_factor(samples, ...),
    "SEM" = filtered_samples_SEM(samples, ...),
    "diagonal-gamma" = filtered_samples_standard(samples, ...),
    stop("Invalid type specified")
  )
}
