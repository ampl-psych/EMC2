.oo_model_list <- function(model) {
  if (is.function(model)) {
    return(model())
  }
  model
}

.oo_expanded_designs <- function(dadm, row_idx = NULL) {
  designs <- attr(dadm, "designs")
  if (is.null(designs)) {
    stop("dadm must have a 'designs' attribute")
  }

  out <- lapply(designs, function(design_mat) {
    expand_idx <- attr(design_mat, "expand")
    if (is.null(expand_idx)) {
      expand_idx <- seq_len(nrow(design_mat))
    }

    expanded <- design_mat[expand_idx, , drop = FALSE]
    if (!is.null(row_idx)) {
      expanded <- expanded[row_idx, , drop = FALSE]
    }
    expanded
  })

  names(out) <- names(designs)
  out
}

.oo_identity_transform <- function(param_names) {
  list(
    func = stats::setNames(rep("identity", length(param_names)), param_names),
    lower = stats::setNames(rep(-Inf, length(param_names)), param_names),
    upper = stats::setNames(rep(Inf, length(param_names)), param_names)
  )
}

.oo_reorder_public_pars <- function(pars, model_list) {
  base_names <- names(model_list$p_types)
  ord <- c(intersect(base_names, colnames(pars)), setdiff(colnames(pars), base_names))
  if (identical(ord, colnames(pars))) {
    return(pars)
  }

  ok <- attr(pars, "ok")
  pars <- pars[, ord, drop = FALSE]
  attr(pars, "ok") <- ok
  pars
}

.oo_particle_matrix <- function(p, dadm, keep_all_columns = FALSE) {
  sampled_p_names <- attr(dadm, "sampled_p_names")
  p_names <- attr(dadm, "p_names")

  if (is.null(sampled_p_names)) {
    stop("dadm must have a 'sampled_p_names' attribute")
  }

  if (is.data.frame(p)) {
    p <- as.matrix(p)
  }

  if (is.null(dim(p))) {
    if (is.null(names(p))) {
      target_names <- if (keep_all_columns && !is.null(p_names)) p_names else sampled_p_names
      if (length(p) != length(target_names)) {
        stop("Unnamed p vector must have length ", length(target_names))
      }
      names(p) <- target_names
    }

    target_names <- if (keep_all_columns && !is.null(p_names)) p_names else sampled_p_names
    if (!all(target_names %in% names(p))) {
      stop("p names must include: ", paste(target_names, collapse = ", "))
    }

    out <- matrix(p[target_names], nrow = 1)
    colnames(out) <- target_names
    return(out)
  }

  if (is.null(colnames(p))) {
    if (ncol(p) != length(sampled_p_names)) {
      stop("p matrix must have columns for: ", paste(sampled_p_names, collapse = ", "))
    }
    colnames(p) <- sampled_p_names
  }

  target_names <- if (keep_all_columns && !is.null(p_names)) p_names else sampled_p_names

  if (all(target_names %in% colnames(p))) {
    p <- p[, target_names, drop = FALSE]
  } else if (!keep_all_columns && !is.null(p_names) && all(p_names %in% colnames(p)) && all(sampled_p_names %in% p_names)) {
    p <- p[, sampled_p_names, drop = FALSE]
  } else {
    stop("p matrix columns must include: ", paste(target_names, collapse = ", "))
  }

  p
}

get_pars_oo <- function(p, dadm, model,
                        pretransformed = FALSE,
                        constants_included = FALSE,
                        return_kernel_matrix = FALSE,
                        return_all_pars = FALSE,
                        kernel_output_codes = 1L) {
  model_list <- .oo_model_list(model)
  particle_matrix <- .oo_particle_matrix(p, dadm, keep_all_columns = constants_included)
  constants <- attr(dadm, "constants")
  if (constants_included || is.null(constants)) {
    constants <- NA
  }
  pretransforms <- if (pretransformed) {
    .oo_identity_transform(colnames(particle_matrix))
  } else {
    model_list$pre_transform
  }

  call_one <- function(cur_particles, cur_dadm, row_idx = NULL) {
    get_pars_c_wrapper_oo(
      particle_matrix = cur_particles,
      data = cur_dadm,
      constants = constants,
      designs = .oo_expanded_designs(dadm, row_idx),
      bounds = model_list$bound,
      transforms = model_list$transform,
      pretransforms = pretransforms,
      trend = model_list$trend,
      return_kernel_matrix = return_kernel_matrix,
      return_all_pars = return_all_pars,
      kernel_output_codes = kernel_output_codes
    )
  }

  if (!("subjects" %in% names(dadm)) || nrow(particle_matrix) == 1L) {
    return(call_one(particle_matrix, dadm))
  }

  if (is.null(rownames(particle_matrix))) {
    stop("Multi-row p matrix must have rownames for every subject in dadm")
  }

  subj_levels <- levels(dadm$subjects)
  used_subjects <- subj_levels[subj_levels %in% as.character(dadm$subjects)]

  if (!all(used_subjects %in% rownames(particle_matrix))) {
    stop("p matrix must have rows named for every subject in dadm")
  }

  particle_matrix <- particle_matrix[used_subjects, , drop = FALSE]

  pieces <- vector("list", length(used_subjects))
  row_ids <- vector("list", length(used_subjects))

  for (i in seq_along(used_subjects)) {
    row_idx <- dadm$subjects == used_subjects[i]
    row_ids[[i]] <- which(row_idx)
    pieces[[i]] <- call_one(
      cur_particles = particle_matrix[i, , drop = FALSE],
      cur_dadm = dadm[row_idx, , drop = FALSE],
      row_idx = row_idx
    )
  }

  out <- matrix(NA_real_, nrow = nrow(dadm), ncol = ncol(pieces[[1]]))
  colnames(out) <- colnames(pieces[[1]])
  for (i in seq_along(pieces)) {
    out[row_ids[[i]], ] <- pieces[[i]]
  }
  out
}

get_pars_matrix_oo <- function(p_vector, dadm, model) {
  model_list <- .oo_model_list(model)
  pars <- get_pars_oo(p_vector, dadm, model_list)
  pars <- model_list$Ttransform(pars, dadm)
  pars <- add_bound(pars, model_list$bound, dadm$lR)
  .oo_reorder_public_pars(pars, model_list)
}
