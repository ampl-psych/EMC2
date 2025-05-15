#' Create Group-Level Design Matrices
#'
#' Creates design matrices for group-level parameters based on subject-level design and formulas.
#' This function is used for hierarchical modeling to specify how subject-level parameters
#' vary across groups or conditions.
#'
#' @param formula A list of formulas specifying the relationship between subject-level parameters
#'        and group-level predictors. Each formula should have a subject-level parameter on the left-hand side
#'        and group-level predictors on the right-hand side.
#' @param data The same data as used in the subject-level design. Must include a 'subjects' column.
#' @param subject_design An emc.design object containing the subject-level design.
#' @param contrasts Optional list of contrast matrices to be used for categorical predictors.
#'
#' @return A list of design matrices, one for each parameter specified in the formula. The intercept is
#' automatically included as the group-level mean and is omitted from the design matrices.
#'
#' @details Here it is important to consider the interpretation of the group-level mean. This allows
#' one to add covariates/group-level factors to the model. However, mu, the group-level mean, is still
#' included for all parameters. Mu represents the intercept in the design matrix, this intercept is always
#' added to the group-level model. Therefore, to keep the interpretation of mu as the group-level mean,
#' it is important to ensure that the design matrix has a mean of zero. If not, this function will throw a
#' warning. For some unbalanced designs, this is unavoidable and the warning can be ignored.
#'
#' @examples
#' # Create subject-level design 
#' subj_design <- design(data = forstmann, model = DDM,
#'                       formula = list(v ~ S, a ~ E, t0 ~ 1),
#'                       contrasts = list(S = contr.helmert))
#' # Add some age covariate and roughly demeans
#' # Demeaning is important to ensure that the interpretation of the group-level intercept
#' # is the mean of the group (i.e., 'mu' still represents the group-level mean)
#' forstmann$age <- as.numeric(forstmann$subjects) -mean(as.numeric(forstmann$subjects))
#' # Create fake group column
#' forstmann$group <- ifelse(forstmann$subjects %in%
#'               unique(forstmann$subjects)[seq(1, 19, 2)], "A", "B")
#'
#' # Create group-level design matrices
#' group_des <- group_design(
#'   formula = list(v_S1 ~ age + group, a ~ age),
#'   data = forstmann,
#'   subject_design = subj_design,
#'   contrasts = list(group = contr.bayes)
#' )
#' # Then you can make the emc object with
#' emc <- make_emc(forstmann, subj_design, compress = FALSE, group_design = group_des)
#' @export
group_design <- function(formula, data, subject_design, contrasts = NULL){
  par_names <- names(sampled_pars(subject_design))

  # Extract dependent variables (left hand side) from formula
  lhs_terms <- unlist(lapply(formula, function(x) as.character(terms(x)[[2]])))

  # Check if all dependent variables are in par_names
  if (!all(lhs_terms %in% par_names)) {
    invalid_terms <- lhs_terms[!lhs_terms %in% par_names]
    stop(paste0("Parameter(s) ", paste0(invalid_terms, collapse=", "),
                " in formula not found in subject design parameters"))
  }

  # Extract all variables from formula, both left and right sides
  rhs_terms <- unique(unlist(lapply(formula, function(x) {
    # Get all variables from right hand side of formula
    all.vars(terms(x)[[3]])
  })))

  # Check if all variables are in data
  if (length(rhs_terms) > 0 && !all(rhs_terms %in% names(data))) {
    missing_vars <- rhs_terms[!rhs_terms %in% names(data)]
    stop(paste0("Variable(s) ", paste0(missing_vars, collapse=", "),
                " in formula not found in data"))
  }

  # Check if any factor has multiple levels per subject
  for (var in rhs_terms) {
    if (is.factor(data[[var]])) {
      # Get count of unique levels per subject
      level_counts <- tapply(data[[var]], data$subjects, function(x) length(unique(x)))

      # Check if any subject has more than one level
      if (any(level_counts > 1)) {
        problematic_subjects <- names(level_counts[level_counts > 1])
        stop(paste0("Factor '", var, "' has multiple levels per subject. ",
                    "First problematic subject: ", problematic_subjects[1], " with ",
                    level_counts[problematic_subjects[1]], " levels. ",
                    "Group-level design requires exactly one level per subject for each factor."))
      }
    }
  }

  # Create an empty list to store design matrices for each formula
  design_matrices <- list()

  # Process each formula separately
  for (i in seq_along(formula)) {
    current_formula <- formula[[i]]

    # Get current formula's left-hand side (parameter name)
    param_name <- as.character(terms(current_formula)[[2]])

    # Get all variables used in the current formula (right-hand side)
    current_vars <- all.vars(terms(current_formula)[[3]])

    # Subset data to include only relevant variables
    subset_data <- data[, c("subjects", current_vars), drop = FALSE]

    # Aggregate data by subject
    agg_data <- setNames(
      data.frame(unique(subset_data$subjects)),
      "subjects"
    )

    # For each variable in the formula, add it to the aggregated data
    for (var in current_vars) {
      # Get unique values per subject
      var_values <- aggregate(
        subset_data[[var]],
        by = list(subjects = subset_data$subjects),
        FUN = function(x) x[1]
      )

      # Add to aggregated data
      agg_data[[var]] <- var_values$x

      # Preserve factor levels if applicable
      if (is.factor(subset_data[[var]])) {
        agg_data[[var]] <- factor(agg_data[[var]], levels = levels(subset_data[[var]]))
      }
    }

    # Apply model.matrix with the current formula
    # Instead of constructing a formula string, use the original formula with modified LHS
    rhs_formula <- terms(current_formula)[[3]]
    dm <- build_design(current_formula, agg_data, contrasts.arg = contrasts)
    # # Check if the design matrix has an intercept
    # has_intercept <- "(Intercept)" %in% colnames(dm)
    #
    # if (!has_intercept) {
    #   stop("Intercept-less design matrix not supported yet")
    # }
    #
    # # Drop the intercept column if present
    # if (has_intercept) {
    #   dm <- dm[, colnames(dm) != "(Intercept)", drop = FALSE]
    # }
    is_int <- grepl("(Intercept)", colnames(dm), fixed = TRUE)
    colnames(dm)[is_int] <- param_name
    colnames(dm)[!is_int] <- paste0(param_name, "_", colnames(dm)[!is_int])
    # Store in the list
    design_matrices[[param_name]] <- dm
    # # Check if overall design matrix mean is zero
    # if (ncol(dm) > 0) {
    #   design_mean <- mean(as.matrix(dm))
    #   if (abs(design_mean) > 1e-1) {  # Small threshold for numerical precision
    #     warning(paste0("Design matrix for parameter '", param_name, "' does not have mean zero. ",
    #                    "For factors, consider using zero-sum contrast matrices (e.g., contr.bayes, contr.helmert). ",
    #                    "For covariates, consider centering them. ",
    #                    "This ensures the intercept can be interpreted as the group-level mean."))
    #   }
    # }
  }
  class(design_matrices) <- "emc.group_design"
  return(design_matrices)
}

#' @rdname sampled_pars
#' @export
sampled_pars.emc.group_design <- function(x,model=NULL,doMap=FALSE, add_da = FALSE, all_cells_dm = FALSE){
  if(is.null(x)) return(logical(0))
  par_names <- unlist(lapply(x, colnames))
  par_names <- setNames(rep(0, length(par_names)), par_names)
  return(par_names)
}


add_group_par_names <- function(pars, group_des){
  if(is.null(group_des)) return(pars)
  out <- c()
  for(par in pars){
    if(par %in% names(group_des)){
      out <- c(out, names(sampled_pars.emc.group_design(group_des[par])))
    } else{
      out <- c(out, par)
    }
  }
  return(out)
}

n_additional_group_pars <- function(group_des){
  if(is.null(group_des)) return(0)
  return(length(sampled_pars(group_des)) - length(group_des))
}

add_group_design <- function(par_names, group_designs = NULL, n_subjects) {
  #
  # par_names:   character vector of parameter names (length p)
  # group_designs: named list of existing design matrices, e.g. $"param1" => (n x m1), etc.
  #                    some par_names might not appear here.
  # n_subjects:    integer, number of subjects (rows).
  #
  # Returns a new list, same length as par_names, each entry an (n_subjects x m_k) matrix
  # If param was not in group_designs, we create a single column of 1's of dimension (n_subjects x 1).

  out_list <- vector("list", length(par_names))
  names(out_list) <- par_names

  for (k in seq_along(par_names)) {
    pname <- par_names[k]
    if (!is.null(group_designs[[pname]])) {
      out_list[[k]] <- group_designs[[pname]]
    } else {
      # no existing design => just a column of 1's
      out_list[[k]] <- matrix(1, nrow = n_subjects, ncol = 1)
    }
  }
  return(out_list)
}


build_design <- function(formula, data, contrasts.arg = NULL) {

  ## ---------- helpers --------------------------------------------------
  is_bar <- function(x) is.call(x) && (identical(x[[1]], quote(`|`)) ||
                                         identical(x[[1]], quote(`||`)))
  nobars <- function(term) {                                  # drop bar terms
    if (is_bar(term)) return(NULL)
    if (!is.call(term)) return(term)
    as.call(c(term[[1]], lapply(as.list(term)[-1], nobars)))
  }
  bars_to_plus <- function(term) {                            # (a|g) → a + g
    if (is_bar(term))
      return(as.call(list(quote(`+`), bars_to_plus(term[[2]]), bars_to_plus(term[[3]]))))
    if (!is.call(term)) return(term)
    as.call(c(term[[1]], lapply(as.list(term)[-1], bars_to_plus)))
  }
  find_bars <- function(term) if (is_bar(term)) list(term) else
    if (is.call(term)) unlist(lapply(as.list(term)[-1], find_bars), FALSE)

  ## ---------- 0. pre‑filter contrasts ----------------------------------
  if (!is.null(contrasts.arg)) {
    keep <- intersect(names(contrasts.arg), names(data))
    contrasts.arg <- if (length(keep)) contrasts.arg[keep] else NULL
  }

  ## ---------- 1. safe model frame --------------------------------------
  rhs_safe <- bars_to_plus(formula[[3]])
  mf <- model.frame(as.formula(call("~", rhs_safe)), data, na.action = NULL)

  ## ---------- 2. fixed‑effect matrix -----------------------------------
  rhs_fixed <- nobars(formula[[3]])
  tt_fixed  <- terms(as.formula(call("~", rhs_fixed)))
  fml_fixed <- reformulate(attr(tt_fixed, "term.labels"),
                           intercept = attr(tt_fixed, "intercept"))
  X <- suppressWarnings(model.matrix(fml_fixed, mf, contrasts.arg = contrasts.arg))

  ## ---------- 3. random‑effect matrix ----------------------------------
  Z_parts <- list()
  for (bar in find_bars(formula[[3]])) {

    ## -- 3a. grouping incidence ----------------------------------------
    gname <- deparse(bar[[3]])
    g     <- as.factor(mf[[gname]])
    J     <- model.matrix(~0 + g)
    colnames(J) <- paste0(gname, levels(g))              # idA, idB, …

    ## -- 3b. slope design ---------------------------------------------
    lhs <- bar[[2]]
    tt  <- terms(as.formula(paste("~", deparse(lhs))))
    slope_labels <- attr(tt, "term.labels")
    has_int <- attr(tt, "intercept") == 1

    mm <- if (length(slope_labels)) {
      fml_slope <- reformulate(slope_labels, intercept = FALSE)
      need <- if (is.null(contrasts.arg)) character(0)
      else intersect(names(contrasts.arg), all.vars(fml_slope))
      model.matrix(fml_slope, mf,
                   contrasts.arg = if (length(need)) contrasts.arg[need] else NULL)
    } else matrix(, nrow(mf), 0)

    ## -- 3c. assemble block -------------------------------------------
    Z <- if (has_int) J else NULL
    if (ncol(mm)) {
      Zs <- do.call(cbind,
                    lapply(seq_len(ncol(mm)),
                           \(k) sweep(J, 1, mm[, k], `*`)))
      colnames(Zs) <- outer(paste0(gname, levels(g)),
                            colnames(mm),
                            paste, sep = ":")            # idA:x, …
      Z <- if (is.null(Z)) Zs else cbind(Z, Zs)
    }
    Z_parts[[length(Z_parts) + 1]] <- Z
  }
  Zbig <- if (length(Z_parts)) do.call(cbind, Z_parts) else matrix(0, nrow(mf), 0)

  ## ---------- 4. combined dense matrix ---------------------------------
  cbind(X, Zbig)
}

