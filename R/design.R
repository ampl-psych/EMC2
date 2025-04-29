#' Specify a Design and Model
#'
#' This function combines information regarding the data, type of model, and
#' the model specification.
#'
#' @param formula A list. Contains the design formulae in the
#' format `list(y ~ x, a ~ z)`.
#' @param factors A named list containing all the factor variables that span
#' the design cells and that should be taken into account by the model.
#' The name `subjects` must be used to indicate the participant factor variable,
#' also in the data.
#'
#' Example: `list(subjects=levels(dat$subjects), condition=levels(dat$condition))`
#'
#' @param Rlevels A character vector. Contains the response factor levels.
#' Example: `c("right", "left")`
#' @param model A function, specifies the model type.
#' Choose from the drift diffusion model (`DDM()`, `DDMt0natural()`),
#' the log-normal race model (`LNR()`), the linear ballistic model (`LBA()`),
#' the racing diffusion model (`RDM()`, `RDMt0natural()`), or define your own
#' model functions.
#' @param data A data frame. `data` can be used to automatically detect
#'  `factors`, `Rlevels` and `covariates` in a dataset. The variable `R` needs
#'  to be a factor variable indicating the response variable. Any numeric column
#'  except `trials` and `rt` are treated as covariates, and all remaining factor
#'  variables are internally used in `factors`.
#' @param contrasts Optional. A named list specifying a design matrix.
#' Example for supplying a customized design matrix:
#' `list(lM = matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))))`
#' @param matchfun A function. Only needed for race models. Specifies whether a
#' response was correct or not. Example: `function(d)d$S==d$lR` where lR refers
#' to the latent response factor.
#' @param constants A named vector that sets constants. Any parameter in
#' `sampled_pars` can be set constant.
#' @param covariates Names of numeric covariates.
#' @param functions List of functions to create new factors based on those in
#' the factors argument. These new factors can then be used in `formula`.
#' @param report_p_vector Boolean. If TRUE (default), it returns the vector of
#' parameters to be estimated.
#' @param custom_p_vector A character vector. If specified, a custom likelihood
#' function can be supplied.
#' @param transform A list with custom transformations to be applied to the parameters of the model,
#' if the conventional transformations aren't desired.
#' See `DDM()` for an example of such transformations
#' @param bound A list with custom bounds to be applied to the parameters of the model,
#' if the conventional bound aren't desired.
#' see `DDM()` for an example of such bounds. Bounds are used to set limits to
#' the likelihood landscape that cannot reasonable be achieved with `transform`
#' @param ... Additional, optional arguments
#'
#' @return A design list.
#' @examples
#'
#' # load example dataset
#' dat <- forstmann
#'
#' # create a function that takes the latent response (lR) factor (d) and returns a logical
#' # defining the correct response for each stimulus. Here the match is simply
#' # such that the S factor equals the latent response factor
#' matchfun <- function(d)d$S==d$lR
#'
#' # When working with lM and lR, it can be useful to design  an
#' # "average and difference" contrast matrix. For binary responses, it has a
#' # simple canonical form
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))
#'
#' # Create a design for a linear ballistic accumulator model (LBA) that allows
#' # thresholds to be a function of E and lR. The final result is a 9 parameter model.
#' design_LBABE <- design(data = dat,model=LBA,matchfun=matchfun,
#'                             formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
#'                             contrasts=list(v=list(lM=ADmat)),
#'                             constants=c(sv=log(1)))
#' @export
#'
#'
design <- function(formula = NULL,factors = NULL,Rlevels = NULL,model,data=NULL,
                       contrasts=NULL,matchfun=NULL,constants=NULL,covariates=NULL,
                       functions=NULL,report_p_vector=TRUE, custom_p_vector = NULL,
                   transform = NULL, bound = NULL, ...){

  optionals <- list(...)
  if(!is.null(optionals$trend)){
    trend <- optionals$trend
  } else {
    trend <- NULL
  }
  if(!is.null(optionals$pre_transform)){
    pre_transform <- optionals$pre_transform
  } else {
    pre_transform <- NULL
  }

  if(any(names(factors) %in% c("trial", "R", "rt", "lR", "lM"))){
    stop("Please do not use any of the following names within factors argument: trial, R, rt, lR, lM")
  }

  if(any(grepl("_", names(factors)))){
    stop("_ in variable names detected. Please refrain from using any underscores.")
  }

  if(!is.null(custom_p_vector)){

    model_list <- function(){list(log_likelihood = model)}
    if(!is.null(list(...)$rfun)){
      model_list$rfun <- list(...)$rfun
    }
    design <- list(Flist = formula, model = model_list, Ffactors = factors)

    attr(design, "sampled_p_names") <-custom_p_vector
    attr(design, "custom_ll") <- TRUE
    class(design) <- "emc.design"
    return(design)
  }
  if (!is.null(data)) {
    if(!"subjects" %in% colnames(data)) stop("make sure subjects identifier is present in data")
    data$subjects <- factor(data$subjects)
    facs <- lapply(data,levels)
    nfacs <- facs[unlist(lapply(facs,is.null))]
    facs <- facs[!unlist(lapply(facs,is.null))]
    Rlevels <- facs[["R"]]
    factors <- facs[names(facs)!="R"]
    nfacs <- nfacs[!(names(nfacs) %in% c("trials","rt"))]
    if (length(nfacs)>0) covariates <- nfacs
  }
  # Check if all parameters in the model are specified in the formula
  nams <- unlist(lapply(formula,function(x) as.character(stats::terms(x)[[2]])))
  if (!all(sort(names(model()$p_types)) %in% sort(nams)) & is.null(custom_p_vector)){
    p_types <- model()$p_types
    not_specified <- sort(names(p_types))[!sort(names(p_types)) %in% sort(nams)]
    message(paste0("Parameter(s) ", paste0(not_specified, collapse = ", "), " not specified in formula and assumed constant."))
    additional_constants <- p_types[not_specified]
    names(additional_constants) <- not_specified
    constants <- c(constants, additional_constants[!names(additional_constants) %in% names(constants)])
    for(add_constant in not_specified) formula[[length(formula)+ 1]] <- as.formula(paste0(add_constant, "~ 1"))
  }

  design <- list(Flist=formula,Ffactors=factors,Rlevels=Rlevels,
                 Clist=contrasts,matchfun=matchfun,constants=constants,
                 Fcovariates=covariates,Ffunctions=functions,model=model)
  class(design) <- "emc.design"
  p_vector <- sampled_pars(design,model)
  lhs_terms <- unlist(lapply(formula, function(x) as.character(stats::terms(x)[[2]])))

  # Check if any terms are not in model parameters
  if (!is.null(formula) && !all(lhs_terms %in% names(model()$p_types))) {
    invalid_terms <- lhs_terms[!lhs_terms %in% names(model()$p_types)]
    stop(paste0("Parameter(s) ", paste0(invalid_terms, collapse=", "),
                " in formula not found in model p_types"))
  }
  model_list <- model()
  model_list$transform <- fill_transform(transform,model)
  model_list$bound <- fill_bound(bound,model)
  model_list$pre_transform <- fill_transform(pre_transform, model = model, p_vector = p_vector, is_pre = TRUE)
  model <- function(){return(model_list)}
  design$model <- model
  attr(design,"p_vector") <- p_vector
  if (report_p_vector) {
    summary(design)
  }
  return(design)
}

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
#' @details Here it is important to consider the interpreation of the group-level mean. This allows
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
  lhs_terms <- unlist(lapply(formula, function(x) as.character(stats::terms(x)[[2]])))

  # Check if all dependent variables are in par_names
  if (!all(lhs_terms %in% par_names)) {
    invalid_terms <- lhs_terms[!lhs_terms %in% par_names]
    stop(paste0("Parameter(s) ", paste0(invalid_terms, collapse=", "),
                  " in formula not found in subject design parameters"))
  }

  # Extract all variables from formula, both left and right sides
  rhs_terms <- unique(unlist(lapply(formula, function(x) {
    # Get all variables from right hand side of formula
    all.vars(stats::terms(x)[[3]])
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
    param_name <- as.character(stats::terms(current_formula)[[2]])

    # Get all variables used in the current formula (right-hand side)
    current_vars <- all.vars(stats::terms(current_formula)[[3]])

    # Check if this is an intercept-only formula
    if (length(current_vars) == 0 && deparse(stats::terms(current_formula)[[3]]) == "1") {
      warning(paste0("Formula '", param_name, " ~ 1' detected. Note that the intercept is automatically ",
                     "included as the group-level mean and does not need to be specified."))
    }

    # Subset data to include only relevant variables
    subset_data <- data[, c("subjects", current_vars), drop = FALSE]

    # Aggregate data by subject
    agg_data <- stats::setNames(
      data.frame(unique(subset_data$subjects)),
      "subjects"
    )

    # For each variable in the formula, add it to the aggregated data
    for (var in current_vars) {
      # Get unique values per subject
      var_values <- stats::aggregate(
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
    rhs_formula <- stats::terms(current_formula)[[3]]
    if (!is.null(contrasts)) {
      # Filter contrasts to only include variables in the current formula
      formula_vars <- all.vars(rhs_formula)
      filtered_contrasts <- contrasts[names(contrasts) %in% formula_vars]

      if (length(filtered_contrasts) > 0) {
        dm <- stats::model.matrix(stats::reformulate(termlabels = deparse(rhs_formula)),
                                data = agg_data,
                                contrasts.arg = filtered_contrasts)
      } else {
        dm <- stats::model.matrix(stats::reformulate(termlabels = deparse(rhs_formula)),
                                data = agg_data)
      }
    } else {
      dm <- stats::model.matrix(stats::reformulate(termlabels = deparse(rhs_formula)),
                              data = agg_data)
    }

    # Check if the design matrix has an intercept
    has_intercept <- "(Intercept)" %in% colnames(dm)

    if (!has_intercept) {
      stop("Intercept-less design matrix not supported yet")
    }

    # Drop the intercept column if present
    if (has_intercept) {
      dm <- dm[, colnames(dm) != "(Intercept)", drop = FALSE]
    }

    # Store in the list
    design_matrices[[param_name]] <- dm

    # Check if overall design matrix mean is zero
    if (ncol(dm) > 0) {
      design_mean <- mean(as.matrix(dm))
      if (abs(design_mean) > 1e-1) {  # Small threshold for numerical precision
        warning(paste0("Design matrix for parameter '", param_name, "' does not have mean zero. ",
                      "For factors, consider using zero-sum contrast matrices (e.g., contr.bayes, contr.helmert). ",
                      "For covariates, consider centering them. ",
                      "This ensures the intercept can be interpreted as the group-level mean."))
      }
    }
  }

  return(design_matrices)
}

#' Contrast Enforcing Equal Prior Variance on each Level
#'
#' Typical contrasts impose different levels of marginal prior variance for the different levels.
#' This contrast can be used to ensure that each level has equal marginal priors (Rouder, Morey, Speckman, & Province; 2012).
#'
#' @param n An integer. The number of items for which to create the contrast
#'
#' @return A contrast matrix.
#' @export
#' @examples{
#' design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.bayes),
#' formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' constants=c(s=log(1)))
#' }
contr.bayes <- function(n) {
  if (length(n) <= 1L) {
    if (is.numeric(n) && length(n) == 1L && n > 1L)
      levels <- seq_len(n)
    else stop("not enough degrees of freedom to define contrasts")
  }
  else levels <- n
  levels <- as.character(levels)
  n <- length(levels)
  cont <- diag(n)
  a <- n
  I_a <- diag(a)
  J_a <- matrix(1, nrow = a, ncol = a)
  Sigma_a <- I_a - J_a/a
  cont <- eigen(Sigma_a)$vectors[,seq_len(a-1), drop = FALSE]
  return(cont)
}


#' Contrast Enforcing Increasing Estimates
#'
#' Each level will be estimated additively from the previous level
#'
#' @param n an integer. The number of items for which to create the contrast.
#'
#' @return a contrast matrix.
#' @export
#' @examples{
#' design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.increasing),
#' formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' constants=c(s=log(1)))
#' }
contr.increasing <- function(n)
{
  if (length(n) <= 1L) {
    if (is.numeric(n) && length(n) == 1L && n > 1L)
      levels <- seq_len(n)
    else stop("not enough degrees of freedom to define contrasts")
  }
  else levels <- n
  levels <- as.character(levels)
  n <- length(levels)
  contr <- matrix(0,nrow=n,ncol=n-1,dimnames=list(NULL,2:n))
  contr[lower.tri(contr)] <- 1
  contr
}

#' Contrast Enforcing Decreasing Estimates
#'
#' Each level will be estimated as a reduction from the previous level
#'
#' @param n an integer. The number of items for which to create the contrast.
#'
#' @return a contrast matrix.
#' @export
#' @examples{
#' design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.decreasing),
#' formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' constants=c(s=log(1)))
#' }
contr.decreasing <- function(n) {
  out <- contr.increasing(n)
  out[dim(out)[1]:1,]
}

#' Anova Style Contrast Matrix
#'
#' Similar to `contr.helmert`, but then scaled to estimate differences between conditions. Use in `design()`.
#'
#' @param n An integer. The number of items for which to create the contrast
#'
#' @return A contrast matrix.
#' @export
#' @examples{
#' design_DDMaE <- design(data = forstmann,model=DDM, contrasts = list(E = contr.anova),
#' formula =list(v~S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' constants=c(s=log(1)))
#' }

contr.anova <- function(n) {
  if (length(n) <= 1L) {
    if (is.numeric(n) && length(n) == 1L && n > 1L)
       levels <- seq_len(n) else
         stop("not enough degrees of freedom to define contrasts")
  } else levels <- n
  levels <- as.character(levels)
  n <- length(levels)
  contr <- stats::contr.helmert(n)
  contr/rep(2*apply(abs(contr),2,max),each=dim(contr)[1])
}

add_accumulators <- function(data,matchfun=NULL,simulate=FALSE,type="RACE", Fcovariates=NULL) {
  if (!is.factor(data$R)) stop("data must have a factor R")
  factors <- names(data)[!names(data) %in% c("R","rt","trials",Fcovariates)]
if (type=="DDM") {
    datar <- cbind(data,lR=factor(rep(levels(data$R)[1],dim(data)[1]),
      levels=levels(data$R)),lM=factor(rep(TRUE,dim(data)[1])))
  }
  if (type %in% c("RACE","SDT")) {
    nacc <- length(levels(data$R))
    datar <- cbind(do.call(rbind,lapply(1:nacc,function(x){data})),
                   lR=factor(rep(levels(data$R),each=dim(data)[1]),levels=levels(data$R)))
    datar <- datar[order(rep(1:dim(data)[1],nacc),datar$lR),]
    if (!is.null(matchfun)) {
      lM <- matchfun(datar)
      if (!is.factor(lM)){
        datar$lM <- factor(lM)
      } else{
        datar$lM <- factor(lM,levels=levels(lM))
      }
    }
    # Advantage NAFC
    nam <- unlist(lapply(strsplit(dimnames(datar)[[2]],"lS"),function(x)x[[1]]))
    islS <- nam ==""
    if (any(islS)) {
      if (sum(islS) != length(levels(data$R)))
        stop("The number of lS columns in the data must equal the length of Rlevels")
      lR <- unlist(lapply(strsplit(dimnames(datar)[[2]],"lS"),function(x){
        if (length(x)==2) x[[2]] else NULL}))
      if (!all(lR %in% levels(data$R)))
        stop("x  in lSx must be in Rlevels")
      if (any(names(datar=="lSmagnitude")))
        stop("Do not use lSmagnitude as a factor name")
      lSmagnitude <- as.character(datar$lR)
      for (i in levels(datar$lR)) {
        isin <- datar$lR==i
        lSmagnitude[isin] <- as.character(datar[isin,paste0("lS",i)])
      }
      factors <- factors[!(factors %in% dimnames(datar)[[2]][islS])]
      # datar <- datar[,!islS]
      datar$lSmagnitude <- as.numeric(lSmagnitude)
    }
  }
  if (type %in% c("BE","TC")) {
    datar <- cbind(do.call(rbind,lapply(1:2,function(x){data})),
      lR=factor(rep(1:2,each=dim(data)[1]),levels=1:2))
    if (!is.null(matchfun)) {
      lM <- matchfun(datar)
      if (any(is.na(lM)) || !(is.logical(lM)))
        stop("matchfun not scoring properly")
      datar$lM <- factor(lM)
    }
  }
  row.names(datar) <- NULL
  if (simulate) datar$rt <- NA else {
    R <- datar$R
    R[is.na(R)] <- levels(datar$lR)[1]

    if (type %in% c("BE","TC")) datar$winner <- NA else
      datar$winner <- datar$lR==R
    # datar$winner[is.na(datar$winner)] <- FALSE
  }
  # # sort cells together
  # if ("trials" %in% names(data)){
  #   if(length(factors) > 1){
  #     datar[order(apply(datar[,c(factors)],1,paste,collapse="_"), as.numeric(datar$trials),as.numeric(datar$lR)),]
  #   } else{
  #     datar[order(datar[,c(factors)], as.numeric(datar$trials),as.numeric(datar$lR)),]
  #   }
  # }
  # else{
  #   if(length(factors) > 1){
  #     datar[order(apply(datar[,c(factors)],1,paste,collapse="_"), as.numeric(datar$lR)),]
  #   } else{
  #     datar[order(datar[,c(factors)], as.numeric(datar$lR)),]
  #   }
  # }
  datar
}

design_model_custom_ll <- function(data, design, model){
  if (!is.factor(data$subjects)) {
    data$subjects <- factor(data$subjects)
    warning("subjects column was converted to a factor")
  }
  dadm <- data
  model_input <- model
  attr(dadm, "model") <- function(){
    return(list(log_likelihood = model_input))
  }
  attr(dadm,"sampled_p_names") <- attr(design, "sampled_p_names")
  attr(dadm, "custom_ll") <- TRUE
  return(dadm)

}


compress_dadm <- function(da,designs,Fcov,Ffun)
    # out keeps only unique rows in terms of all parameters design matrices
    # R, lR and rt (at given resolution) from full data set
  {
    nacc <- length(unique(da$lR))
    # contract output
    cells <- paste(
      apply(do.call(cbind,lapply(designs,function(x){
        apply(x[attr(x,"expand"),,drop=FALSE],1,paste,collapse="_")})
      ),1,paste,collapse="+"),da$subjects,da$R,da$lR,da$rt,sep="+")
    # Make sure that if row is included for a trial so are other rows
    if (!is.null(Fcov)) {
      if (is.null(names(Fcov))) nFcov <- Fcov else nFcov <- names(Fcov)
      cells <- paste(cells,apply(da[,nFcov,drop=FALSE],1,paste,collapse="+"),sep="+")
    }
    if (!is.null(Ffun))
      cells <- paste(cells,apply(da[,Ffun,drop=FALSE],1,paste,collapse="+"),sep="+")

    if (nacc>1) cells <- paste0(rep(apply(matrix(cells,nrow=nacc),2,paste0,collapse="_"),
                                    each=nacc),rep(1:nacc,times=length(cells)/nacc),sep="_")

    contract <- !duplicated(cells)
    out <- da[contract,,drop=FALSE]
    attr(out,"contract") <- contract
    attr(out,"expand") <- as.numeric(factor(cells,levels=unique(cells)))
    lR1 <- da$lR==levels(da$lR)[[1]]
    attr(out,"expand_winner") <- as.numeric(factor(cells[lR1],levels=unique(cells[lR1])))
    attr(out,"s_expand") <- da$subjects
    attr(out,"designs") <- lapply(designs,function(x){
      attr(x,"expand") <- attr(x,"expand")[contract]; x})

    # indices to use to contract further ignoring rt then expand back
    cells_nort <- paste(
      apply(do.call(cbind,lapply(designs,function(x){
        apply(x[attr(x,"expand"),,drop=FALSE],1,paste,collapse="_")})
      ),1,paste,collapse="+"),da$subjects,da$R,da$lR,sep="+")[contract]
    attr(out,"unique_nort") <- !duplicated(cells_nort)
    attr(out,"expand_nort") <- as.numeric(factor(cells_nort,levels=unique(cells_nort)))

    # cells_nort <- paste(
    #   apply(do.call(cbind,lapply(designs,function(x){
    #     apply(x[attr(x,"expand"),,drop=FALSE],1,paste,collapse="_")})
    #   ),1,paste,collapse="+"),da$subjects,da$R,da$lR,sep="+")[contract]
    # attr(out,"unique_nort") <- !duplicated(cells_nort)
    # # Only first level WHY????
    # cells <- cells[da$lR==levels(da$lR)[1]]
    # cells_nort <- cells_nort[out$lR==levels(out$lR)[1]]
    # attr(out,"expand_nort") <- as.numeric(factor(cells_nort,
    #    levels=unique(cells_nort)))[as.numeric(factor(cells,levels=unique(cells)))]

    # indices to use to contract ignoring rt and response (R), then expand back
    cells_nortR <- paste(apply(do.call(cbind,lapply(designs,function(x){
      apply(x[attr(x,"expand"),,drop=FALSE],1,paste,collapse="_")})),1,paste,collapse="+"),
      da$subjects,da$lR,sep="+")[contract]
    attr(out,"unique_nortR") <- !duplicated(cells_nortR)
    attr(out,"expand_nortR") <- as.numeric(factor(cells_nortR,levels=unique(cells_nortR)))

    # # indices to use to contract ignoring rt and response (R), then expand back
    # cells_nortR <- paste(apply(do.call(cbind,lapply(designs,function(x){
    #   apply(x[attr(x,"expand"),,drop=FALSE],1,paste,collapse="_")})),1,paste,collapse="+"),
    #   da$subjects,da$lR,sep="+")[contract]
    # attr(out,"unique_nortR") <- !duplicated(cells_nortR)
    # # Only first level WHY????
    # cells_nortR <- cells_nortR[out$lR==levels(out$lR)[1]]
    # attr(out,"expand_nortR") <- as.numeric(factor(cells_nortR,
    #    levels=unique(cells_nortR)))[as.numeric(factor(cells,levels=unique(cells)))]

    # Lower censor
    if (!any(is.na(out$rt))) { # Not a choice only model
      winner <- out$lR==levels(out$lR)[[1]]
      ok <- out$rt[winner]==-Inf
      if (any(ok)) {
        ok[ok] <- 1:sum(ok)
        attr(out,"expand_lc") <- ok[attr(out,"expand_winner")] + 1
      }
      # Upper censor
      ok <- out$rt[winner]==Inf
      if (any(ok)) {
        ok[ok] <- 1:sum(ok)
        attr(out,"expand_uc") <- ok[attr(out,"expand_winner")] + 1
      }
    }
    out
}

check_rt <- function(b,d,upper=TRUE)
  # Check bounds respected if present
{
  if (!all(sort(levels(d$subjects))==sort(names(b))))
    stop("Bound vector must have same names as subjects")
  d <- d[!is.na(d$rt),]
  d <- d[is.finite(d$rt),]
  bound <- d$subjects
  levels(bound) <- unlist(b)[levels(bound)]
  if (upper)
    ok <- all(d$rt < as.numeric(as.character(bound))) else
      ok <- all(d$rt > as.numeric(as.character(bound)))
  if (!all(ok)) stop("Bound not respected in data")
}

rt_check_function <- function(data){
  # Truncation
  if (!is.null(attr(data,"UT"))) {
    if (length(attr(data,"UT"))==1 && is.null(names(attr(data,"UT"))))
      attr(data,"UT") <- stats::setNames(rep(attr(data,"UT"),length(levels(data$subjects))),
                                         levels(data$subjects))
    check_rt(attr(data,"UT"),data)
  }
  if (!is.null(attr(data,"LT"))) {
    if (length(attr(data,"LT"))==1 && is.null(names(attr(data,"LT"))))
      attr(data,"LT") <- stats::setNames(rep(attr(data,"LT"),length(levels(data$subjects))),
                                         levels(data$subjects))
    if (any(attr(data,"LT")<0)) stop("Lower truncation cannot be negative")
    check_rt(attr(data,"LT"),data,upper=FALSE)
  }
  if (!is.null(attr(data,"UT")) & !is.null(attr(data,"LT"))) {
    DT <- attr(data,"UT") - attr(data,"LT")
    if (!is.null(DT) && any(DT<0)) stop("UT must be greater than LT")
  }

  # Censoring
  if (!is.null(attr(data,"UC"))) {
    if (length(attr(data,"UC"))==1 && is.null(names(attr(data,"UC"))))
      attr(data,"UC") <- stats::setNames(rep(attr(data,"UC"),length(levels(data$subjects))),
                                         levels(data$subjects))
    check_rt(attr(data,"UC"),data)
    if (!is.null(attr(data,"UT")) && attr(data,"UT") < attr(data,"UC"))
      stop("Upper censor must be less than upper truncation")
  }
  if (!is.null(attr(data,"LC"))) {
    if (length(attr(data,"LC"))==1 && is.null(names(attr(data,"LC"))))
      attr(data,"LC") <- stats::setNames(rep(attr(data,"LC"),length(levels(data$subjects))),
                                         levels(data$subjects))
    if (any(attr(data,"LC")<0)) stop("Lower censor cannot be negative")
    check_rt(attr(data,"LC"),data,upper=FALSE)
    if (!is.null(attr(data,"LT")) && attr(data,"LT") > attr(data,"LC"))
      stop("Lower censor must be greater than lower truncation")
  }
  if (any(data$rt[!is.na(data$rt)]==-Inf) & is.null(attr(data,"LC")))
    stop("Data must have an LC attribute if any rt = -Inf")
  if (any(data$rt[!is.na(data$rt)]==Inf) & is.null(attr(data,"UC")))
    stop("Data must have an UC attribute if any rt = Inf")
  if (!is.null(attr(data,"UC"))) check_rt(attr(data,"UC"),data)
  if (!is.null(attr(data,"LC"))) check_rt(attr(data,"LC"),data,upper=FALSE)
  if (!is.null(attr(data,"UC")) & !is.null(attr(data,"LC"))) {
    DC <- attr(data,"UC") - attr(data,"LC")
    if (!is.null(DC) && any(DC<0)) stop("UC must be greater than LC")
  }
}


design_model <- function(data,design,model=NULL,
                         add_acc=TRUE,rt_resolution=0.02,verbose=TRUE,
                         compress=TRUE,rt_check=TRUE, add_da = FALSE, all_cells_dm = FALSE)
{
  if (is.null(model)) {
    if (is.null(design$model))
      stop("Model must be supplied if it has not been added to design")
    model <- design$model
  }
  if (model()$type=="SDT") rt_check <- FALSE
  if(grepl("MRI", model()$type)){
    dadm <- data
    attr(dadm, "design_matrix") <- attr(design, "design_matrix")
    # attr(design, "design_matrix") <- NULL
    p_names <- names(model()$p_types)
    attr(dadm,"p_names") <- p_names
    sampled_p_names <- p_names[!(p_names %in% names(design$constants))]
    attr(dadm,"sampled_p_names") <- sampled_p_names
    return(dadm)
  }
  if (any(names(model()$p_types) %in% names(data)))
    stop("Data cannot have columns with the same names as model parameters")
  if (!is.factor(data$subjects)) {
    data$subjects <- factor(data$subjects)
    warning("subjects column was converted to a factor")
  }

  if (!any(names(data)=="trials")) data$trials <- 1:dim(data)[1]
  if(rt_check){rt_check_function(data)}
  if (!add_acc) da <- data else
    da <- add_accumulators(data,design$matchfun,type=model()$type,Fcovariates=design$Fcovariates)
  order_idx <- order(da$subjects)
  da <- da[order_idx,] # fixes different sort in add_accumulators depending on subject type

  if (!is.null(design$Ffunctions)) for (i in names(design$Ffunctions)) {
    newF <- stats::setNames(data.frame(design$Ffunctions[[i]](da)),i)
    da <- cbind.data.frame(da,newF)
  }

  if (is.null(model()$p_types) | is.null(model()$Ttransform))
    stop("p_types and Ttransform must be supplied")
  if (!all(unlist(lapply(design$Flist,class))=="formula"))
    stop("Flist must contain formulas")
  nams <- unlist(lapply(design$Flist,function(x)as.character(stats::terms(x)[[2]])))
  names(design$Flist) <- nams
  if (is.null(design$Clist)) design$Clist=list(stats::contr.treatment)
  if (!is.list(design$Clist)) stop("Clist must be a list")
  pnames <- names(model()$p_types)
  if (!is.list(design$Clist[[1]])[1]){
    design$Clist <- stats::setNames(lapply(1:length(pnames),
                                           function(x)design$Clist),pnames)
  } else {
   missing_p_types <- pnames[!(pnames %in% names(design$Clist))]
   if (length(missing_p_types)>0) {
     nok <- length(design$Clist)
      for (i in 1:length(missing_p_types)) {
        design$Clist[[missing_p_types[i]]] <- list(stats::contr.treatment)
        names(design$Clist)[nok+i] <- missing_p_types[i]
      }
    }
  }
  for (i in pnames) attr(design$Flist[[i]],"Clist") <- design$Clist[[i]]

  out <- lapply(design$Flist,make_dm,da=da,Fcovariates=design$Fcovariates, add_da = add_da, all_cells_dm = all_cells_dm)
  if (!is.null(rt_resolution) & !is.null(da$rt)) da$rt <- round(da$rt/rt_resolution)*rt_resolution
  if (compress){
    dadm <- compress_dadm(da,designs=out, Fcov=design$Fcovariates,Ffun=names(design$Ffunctions))
    # Change expansion names
    # attr(dadm,"expand_all") <- attr(dadm,"expand")
    attr(dadm,"expand") <- attr(dadm,"expand_winner")
    attr(dadm,"expand_winner") <- NULL
  }  else {
    dadm <- da
    attr(dadm,"designs") <- out
    attr(dadm,"s_expand") <- da$subjects
    # attr(dadm,"expand_all") <- 1:nrow(dadm)
    if(is.null(dadm$lR)){
      attr(dadm,"expand") <- 1:nrow(dadm)
    } else{
      attr(dadm,"expand") <- 1:(nrow(dadm)/length(unique(dadm$lR)))
    }
  }
  p_names <-  unlist(lapply(out,function(x){dimnames(x)[[2]]}),use.names=FALSE)
  bad_constants <- names(design$constants)[!(names(design$constants) %in% p_names)]
  if (length(bad_constants) > 0)
    stop("Constant(s) ",paste(bad_constants,collapse=" ")," not in design")

  # Pick out constants
  sampled_p_names <- p_names[!(p_names %in% names(design$constants))]
  attr(dadm,"p_names") <- p_names
  attr(dadm,"sampled_p_names") <- sampled_p_names
  if (model()$type=="DDM") nunique <- dim(dadm)[1] else
    nunique <- dim(dadm)[1]/length(levels(dadm$lR))
  if (verbose & compress) message("Likelihood speedup factor: ",
  round(dim(da)[1]/dim(dadm)[1],1)," (",nunique," unique trials)")

  attr(dadm,"model") <- model
  attr(dadm,"constants") <- design$constants

  if (add_acc) {
    attr(dadm, "ok_dadm_winner") <- is.finite(dadm$rt) & dadm$winner
    attr(dadm, "ok_dadm_looser") <- is.finite(dadm$rt) & !dadm$winner
    attr(dadm, "ok_da_winner") <- attr(dadm, "ok_dadm_winner")[attr(dadm,"expand")]
    attr(dadm, "ok_da_looser") <- attr(dadm, "ok_dadm_looser")[attr(dadm,"expand")]
  }
  attr(dadm,"ok_trials") <- is.finite(data$rt)
  attr(dadm,"s_data") <- data$subjects
  if (model()$type %in% c("BE","TC"))
    attr(dadm,"dL") <- get_dL(length(levels(dadm$R)),model()$type)
  dadm
}


make_dm <- function(form,da,Clist=NULL,Fcovariates=NULL, add_da = FALSE, all_cells_dm = FALSE)
  # Makes a design matrix based on formula form from augmented data frame da
{

  compress_dm <- function(dm, da = NULL, all_cells_dm = FALSE)
    # out keeps only unique rows, out[attr(out,"expand"),] gets back original.
  {
    cells <- apply(dm,1,paste,collapse="_")
    ass <- attr(dm,"assign")
    contr <- attr(dm,"contrasts")
    if(!is.null(da)){
      dups <- duplicated(paste0(cells, apply(da, 1, paste0, collapse = "_")))
    } else{
      dups <- duplicated(cells)
    }
    out <- dm[!dups,,drop=FALSE]
    if(!is.null(da) & !all_cells_dm){
      if(nrow(da) != 0){
        out <- cbind(da[!dups,colnames(da) != "subjects",drop=FALSE], out)
      }
    }
    attr(out,"expand") <- as.numeric(factor(cells,levels=unique(cells)))
    attr(out,"assign") <- ass
    attr(out,"contrasts") <- contr
    out
  }

  if (is.null(Clist)) Clist <- attr(form,"Clist")
  pnam <- stats::terms(form)[[2]]
  da[[pnam]] <- 1
  for (i in names(Clist)) if (i %in% names(da)) {
    if (!is.factor(da[[i]]))
      stop(i," must be a factor (design factors has a parameter name?)")
    levs <- levels(da[[i]])
    nl <- length(levs)
    if (class(Clist[[i]])[1]=="function")
      stats::contrasts(da[[i]]) <- do.call(Clist[[i]],list(n=levs)) else {
        if (!is.matrix(Clist[[i]]) || dim(Clist[[i]])[1]!=nl) {
          if (all(levs %in% row.names(Clist[[i]]))) # design with missing cells
            Clist[[i]] <- Clist[[i]][levs,] else
            stop("Clist for ",i," not a ",nl," row matrix")
        } else dimnames(Clist[[i]])[[1]] <- levs
        stats::contrasts(da[[i]],how.many=dim(Clist[[i]])[2]) <- Clist[[i]]
      }
  }
  out <- stats::model.matrix(form,da)
  if (dim(out)[2]==1) dimnames(out)[[2]] <- as.character(pnam) else {
    if (attr(stats::terms(form),"intercept")!=0) {
      cnams <- paste(pnam,dimnames(out)[[2]][-1],sep="_")
      dimnames(out)[[2]] <- c(pnam,cnams)
    } else dimnames(out)[[2]] <- paste(pnam,dimnames(out)[[2]],sep="_")
  }
  if(add_da){
    da <- da[,all.vars(form)[-1], drop = F]
    out <- compress_dm(out, da, all_cells_dm)
  } else{
    out <- compress_dm(out)
  }
  return(out)
}

# data generation

# Used in make_data and make_emc
add_trials <- function(dat)
  # Add trials column, 1:n for each subject
{
  n <- table(dat$subjects)
  if (!any(names(dat)=="trials")) dat <- cbind.data.frame(dat,trials=NA)
  for (i in names(n)) dat$trials[dat$subjects==i] <- 1:n[i]
  dat
}

dm_list <- function(dadm)
  # Makes data model into subjects list for use by likelihood
  # Assumes each subject has the same design.
{

  sub_design <- function(designs,isin)
    lapply(designs,function(x) {
      attr(x,"expand") <- attr(x,"expand")[isin]
      x
    })


  model <- attr(dadm,"model")
  p_names <- attr(dadm,"p_names")
  sampled_p_names <- attr(dadm,"sampled_p_names")
  designs <- attr(dadm,"designs")
  expand <- attr(dadm,"expand")
  s_expand <- attr(dadm,"s_expand")
  unique_nort <- attr(dadm,"unique_nort")
  expand_nort <- attr(dadm,"expand_nort")
  unique_nortR <- attr(dadm,"unique_nortR")
  expand_nortR <- attr(dadm,"expand_nortR")
  # ok_trials <- attr(dadm,"ok_trials")
  # ok_dadm_winner <- attr(dadm,"ok_dadm_winner")
  # ok_dadm_looser <- attr(dadm,"ok_dadm_looser")
  # ok_da_winner <- attr(dadm,"ok_da_winner")
  # ok_da_looser <- attr(dadm,"ok_da_looser")
  # expand_uc <- attr(dadm,"expand_uc")
  # expand_lc <- attr(dadm,"expand_lc")
  dms_mri <- attr(dadm, "design_matrix")

  # winner on expanded dadm
  expand_winner <- attr(dadm,"expand")
  # subjects for first level of lR in expanded dadm
  slR1=dadm$subjects[expand][dadm$lR[expand]==levels(dadm$lR)[[1]]]

  dl <- stats::setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                        levels(dadm$subjects))
  for (i in levels(dadm$subjects)) {
    isin <- dadm$subjects==i         # dadm
    dl[[i]] <- dadm[isin,]
    dl[[i]]$subjects <- factor(as.character(dl[[i]]$subjects))
    if(is.null(attr(dadm, "custom_ll"))){

      isin1 <- s_expand==i             # da
      isin2 <- attr(dadm,"s_data")==i  # data
      if(length(isin2) > 0){
        attr(dl[[i]],"expand") <- expand_winner[isin2]-min(expand_winner[isin2]) + 1
      }
      attr(dl[[i]],"model") <- NULL
      attr(dl[[i]],"p_names") <- p_names
      attr(dl[[i]],"sampled_p_names") <- sampled_p_names
      attr(dl[[i]],"designs") <- sub_design(designs,isin)
      # if(!is.null(expand)) attr(dl[[i]],"expand_all") <- expand[isin1]-min(expand[isin1]) + 1
      attr(dl[[i]],"contract") <- NULL
      attr(dl[[i]],"expand_winner") <- NULL
      attr(dl[[i]],"ok_dadm_winner") <- NULL
      attr(dl[[i]],"ok_dadm_looser") <- NULL
      attr(dl[[i]],"ok_da_winner") <- NULL
      attr(dl[[i]],"ok_da_looser") <- NULL
      attr(dl[[i]],"ok_trials") <- NULL
      attr(dl[[i]],"s_data") <- NULL
      attr(dl[[i]],"s_expand") <- NULL
      attr(dl[[i]],"prior") <- NULL
      if(!is.null(dms_mri)){
        attr(dl[[i]], "designs") <- make_mri_sampling_design(dms_mri[[i]], sampled_p_names)
        attr(dl[[i]], "design_matrix") <- NULL
      }

      attr(dl[[i]], "unique_nort") <- NULL
      attr(dl[[i]], "unique_nortR") <- NULL
      attr(dl[[i]], "expand_nort") <- NULL
      attr(dl[[i]], "expand_nortR") <- NULL

      if (!is.null(attr(dadm,"LT"))){
        attr(dl[[i]],"LT") <- attr(dadm,"LT")[names(attr(dadm,"LT"))==i]
      }
      if (!is.null(attr(dadm,"UT"))){
        attr(dl[[i]],"UT") <- attr(dadm,"UT")[names(attr(dadm,"UT"))==i]
      }
      if (!is.null(attr(dadm,"LC"))){
        attr(dl[[i]],"LC") <- attr(dadm,"LC")[names(attr(dadm,"LC"))==i]
      }
      if (!is.null(attr(dadm,"UC"))){
        attr(dl[[i]],"UC") <- attr(dadm,"UC")[names(attr(dadm,"UC"))==i]
      }

      # TE and BE
      if (!is.null(attr(dadm,"dL")))
        attr(dl[[i]],"dL") <- attr(dadm,"dL")

    }
  }


  return(dl)
}

#' Update EMC Objects to the Current Version
#'
#' This function updates EMC objects created with older versions of the package to be compatible with the current version.
#'
#' @param emc An EMC object to update
#' @return An updated EMC object compatible with the current version
#' @examples
#' # Update the model to current version
#' updated_model <- update2version(samples_LNR)
#'
#' @export
update2version <- function(emc){
  # For older versions, ensure that the class is emc:
  class(emc) <- "emc"
  get_new_model <- function(old_model, pars){
    if(old_model()$c_name == "LBA"){
      model <- LBA
    } else if(old_model()$c_name == "DDM"){
      model <- DDM
    } else if(old_model()$c_name == "RDM"){
      model <- RDM
    } else if(old_model()$c_name == "LNR"){
      model <- LNR
    } else{
      stop("current model not supported for updating, sorry!!")
    }
    model_list <- model()
    model_list$transform <- fill_transform(transform = NULL,model)
    model_list$pre_transform <- fill_transform(transform = NULL, model = model, p_vector = pars, is_pre = TRUE)
    model_list$bound <- fill_bound(bound = NULL,model)
    model <- function(){return(model_list)}
    return(model)
  }

  design_list <- get_design(emc)
  if(is.null(emc[[1]]$type)){
    type <- attr(emc[[1]], "variant_funs")$type
    emc <- lapply(emc, FUN = function(x){
      x$type <- type
      return(x)
    })
  } else{
    type <- emc[[1]]$type
  }
  # Model used to be stored in data
  first_data <- emc[[1]]$data[[1]]
  if(is.null(emc[[1]]$model)){
    if(is.data.frame(first_data)){
      old_model <- attr(first_data, "model")
      new_model <- get_new_model(old_model, sampled_pars(design_list[[1]]))
      design_list[[1]]$model <- new_model
      emc[[1]]$model <- new_model
    } else{
      old_model <- lapply(first_data, function(x) attr(x, "model"))
      new_model <- mapply(get_new_model, old_model, lapply(design_list, sampled_pars))
      design_list <- mapply(function(x, y){
        x$model <- y
        return(list(x))
      }, design_list, new_model)
    }
    emc[[1]]$model <- new_model
  }
  prior_new <- emc[[1]]$prior
  attr(prior_new, "type") <- type
  prior_new <- prior(design_list, type, update = prior_new)
  class(prior_new) <- "emc.prior"
  emc <- lapply(emc, function(x){
    x$prior <- prior_new
    return(x)
  })
  class(emc) <- "emc"
  attr(emc, "design_list") <- NULL
  return(emc)
}


# Some s3 classes for design objects ---------------------------------------------------------
#' Parameter Mapping Back to the Design Factors
#'
#' Maps parameters of the cognitive model back to the experimental design. If p_vector
#' is left unspecified will print a textual description of the mapping.
#' Otherwise the p_vector can be created using ``sampled_pars()``.
#' The returned matrix shows whether/how parameters
#' differ across the experimental factors.
#'
#' @param x an `emc`, `emc.prior` or `emc.design` object
#' @param p_vector Optional. Specify parameter vector to get numeric mappings.
#' Must be in the form of ``sampled_pars(design)``
#' @param model Optional model type (if not already specified in ``design``)
#' @param digits Integer. Will round the output parameter values to this many decimals
#' @param ... optional arguments
#' @param remove_subjects Boolean. Whether to include subjects as a factor in the design
#' @param covariates Covariates specified in the design can be included here.
#' @return Matrix with a column for each factor in the design and for   each model parameter type (``p_type``).
#' @examples
#' # First define a design:
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' mapped_pars(design_DDMaE)
#' # Then create a p_vector:
#' p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'           t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' # This will map the parameters of the p_vector back to the design
#' mapped_pars(design_DDMaE, p_vector)
#'
#' @export
mapped_pars <- function(x, p_vector = NULL, model=NULL,
                        digits=3,remove_subjects=TRUE,
                        covariates=NULL,...)
  # Show augmented data and corresponding mapped parameter
{
  UseMethod("mapped_pars")
}

#' @rdname mapped_pars
#' @export
mapped_pars.emc.design <- function(x, p_vector = NULL, model=NULL,
                                   digits=3,remove_subjects=TRUE,
                                   covariates=NULL,...){
  if(is.null(x)) return(NULL)
  if(is.null(x$Ffactors)){
    x <- x[[1]]
  }
  if(!is.null(attr(x, "custom_ll"))){
    stop("Mapped_pars not available for this design type")
  }
  design <- x
  if(is.null(p_vector)){
    return(verbal_dm(design))
  }
  remove_RACE <- TRUE
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if (is.null(covariates))
    Fcovariates <- design$Fcovariates else
      Fcovariates <- covariates
  if (is.null(model)) if (is.null(design$model))
    stop("Must specify model as not in design") else model <- design$model
  if (remove_subjects) design$Ffactors$subjects <- design$Ffactors$subjects[1]
  if (!is.matrix(p_vector)) p_vector <- make_pmat(p_vector,design)
  dadm <- design_model(make_data(p_vector,design,n_trials=1,Fcovariates=Fcovariates),
                       design,model,rt_check=FALSE,compress=FALSE, verbose = FALSE)
  ok <- !(names(dadm) %in% c("subjects","trials","R","rt","winner"))
  out <- cbind(dadm[,ok],round(get_pars_matrix(p_vector,dadm, design$model()),digits))
  if (model()$type=="SDT")  out <- out[dadm$lR!=levels(dadm$lR)[length(levels(dadm$lR))],]
  if (model()$type=="DDM")  out <- out[,!(names(out) %in% c("lR","lM"))]
  if (any(names(out)=="RACE") && remove_RACE)
    out <- out[as.numeric(out$lR) <= as.numeric(as.character(out$RACE)),,drop=FALSE]
  return(out)
}



#' Get Model Parameters from a Design
#'
#' Makes a vector with zeroes, with names and length corresponding to the
#' model parameters of the design.
#'
#' @param x an `emc.design` object made with `design()` or an `emc` object.
#' @param model a model list. Defaults to the model specified in the design list.
#' @param doMap logical. If `TRUE` will also include an attribute `map`
#' with the design matrices that perform the mapping back to the design
#' @param add_da Boolean. Whether to include the relevant data columns in the map attribute
#' @param all_cells_dm Boolean. Whether to include all levels of a factor in the mapping attribute,
#' even when one is dropped in the design
#'
#'
#' @return Named vector.
#' @examples
#' # First define a design
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then for this design get which cognitive model parameters are sampled:
#' sampled_pars(design_DDMaE)
#'
#' @export
sampled_pars <- function(x,model=NULL,doMap=TRUE, add_da = FALSE, all_cells_dm = FALSE)
{
  UseMethod("sampled_pars")
}

#' @rdname sampled_pars
#' @export
sampled_pars.emc.design <- function(x,model=NULL,doMap=TRUE, add_da = FALSE, all_cells_dm = FALSE){
  design <- x
  if(is.null(design)) return(NULL)
  if("Flist" %in% names(design)){
    design <- list(design)
  }
  out <- c()
  map_list <- list()
  if(is.null(names(design))){
    names(design) <- as.character(1:length(design))
  }
  for(j in 1:length(design)){
    cur_name <- names(design)[j]
    cur_design <- design[[j]]
    if(!is.null(attr(cur_design, "custom_ll"))){
      pars <- numeric(length(attr(cur_design,"sampled_p_names")))
      if(length(design) != 1){
        map_list[[cur_name]] <- NA
        names(pars) <- paste(cur_name,  attr(cur_design,"sampled_p_names"), sep = "|")
      } else{
        names(pars) <- attr(cur_design,"sampled_p_names")
      }
      out <- c(out, pars)
      next
    }
    model <- cur_design$model
    if (is.null(model)) stop("Must supply model as not in design")

    if(grepl("MRI", model()$type)){
      pars <- model()$p_types
      if(length(design) != 1){
        names(pars) <- paste(cur_name,  names(pars), sep = "|")
        map_list[[cur_name]] <- NA
      }
      pars[1:length(pars)] <- 0
      out <- c(out, pars)
      next
    }

    Ffactors=c(cur_design$Ffactors,list(R=cur_design$Rlevels))
    data <- as.data.frame.table(array(dim=unlist(lapply(Ffactors,length)),
                                      dimnames=Ffactors))[,-(length(Ffactors)+1)]
    for (i in names(cur_design$Ffactors))
      data[[i]] <- factor(data[[i]],levels=cur_design$Ffactors[[i]])

    # if (!is.null(design$Ffunctions))
    #   data <- cbind.data.frame(data,data.frame(lapply(design$Ffunctions,function(f){f(data)})))

    if (!is.null(cur_design$Fcovariates)) {
      covs <- matrix(1,nrow=dim(data)[1],ncol=length(cur_design$Fcovariates),
                     dimnames=list(NULL,names(cur_design$Fcovariates)))
      data <- cbind.data.frame(data,covs)
    }
    dadm <- design_model(
      add_accumulators(data,matchfun=cur_design$matchfun,type=model()$type,Fcovariates=cur_design$Fcovariates),
      cur_design,model,add_acc=FALSE,verbose=FALSE,rt_check=FALSE,compress=FALSE, add_da = add_da,
      all_cells_dm = all_cells_dm)
    sampled_p_names <- attr(dadm,"sampled_p_names")
    if(length(design) != 1){
      map_list[[cur_name]] <- lapply(attributes(dadm)$designs,function(x){x[,,drop=FALSE]})
      sampled_p_names <- paste(cur_name, sampled_p_names, sep = "|")
    }
    out <- c(out, stats::setNames(numeric(length(sampled_p_names)),sampled_p_names))
    if(length(design) == 1){
      if (doMap) attr(out,"map") <- lapply(attributes(dadm)$designs,function(x){x[,,drop=FALSE]})
    }
  }
  if(length(design) != 1) attr(out, "map") <- map_list
  if(!add_da & any(duplicated(names(out)))) stop("duplicate parameter names found! Usually this happens when joint designs share indicator names")
  return(out)
}

#' Summary method for emc.design objects
#'
#' Prints a summary of the design object, including sampled parameters and design matrices.
#' For continuous covariates just prints one row, instead of all covariates.
#'
#' @param object An object of class `emc.design` containing the design to summarize
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the design matrices
#' @export
summary.emc.design <- function(object, ...){
  p_vector <- sampled_pars(object)
  cat("\n Sampled Parameters: \n")
  print(names(p_vector))
  cat("\n Design Matrices: \n")
  map_out <- sampled_pars(object,object$model, add_da = TRUE)
  print(attr(map_out, "map"), row.names = FALSE)
  return(invisible(map_out))
}

#' @export
print.emc.design <- function(x, ...){
  if("Ffactors" %in% names(x)){
    x <- list(x)
  }
  for(i in 1:length(x)){
    for(j in 1:length(x[[i]]$Flist)){
      cat(deparse(x[[i]]$Flist[j][[1]]), "\n")
    }
  }
}

#' @rdname plot_design
#' @export
plot_design.emc.design <- function(x, data = NULL, factors = NULL, plot_factor = NULL, n_data_sim = 10,
                            p_vector = NULL, functions = NULL, ...){
  if(is.null(p_vector)) stop("p_vector must be supplied if only the design is given")
  plot(x, p_vector, data = data, factors = factors, plot_factor = plot_factor, n_data_sim = n_data_sim,
       functions = functions, ...)
}


#' Plot method for emc.design objects
#'
#' Makes design illustration by plotting simulated data based on the design
#'
#' @param x An object of class `emc.design` containing the design to plot
#' @param p_vector A named vector of parameter values to use for data generation
#' @param data Optional data frame to overlay on the design plot. If NULL, data will be simulated.
#' @param factors Character vector. Factors to use for varying parameters in the plot
#' @param plot_factor Optional character. Make separate plots for each level of this factor
#' @param n_data_sim Integer. If data is NULL, number of simulated datasets to generate for the plot. Default is 10.
#' @param functions Optional named list of functions that create additional columns in the data
#' @param ... Additional arguments passed to `make_design_plot`
#' @return No return value, called for side effect of plotting
#' @export
plot.emc.design <- function(x, p_vector, data = NULL, factors = NULL, plot_factor = NULL, n_data_sim = 10,
                            functions = NULL, ...){
  if(!"Ffactors" %in% names(x)){
    if(length(x) != 1){
      stop("Current design type not supported for plotting")
    } else{
      x <- x[[1]]
    }
  }
  x$Ffunctions <- c(x$Ffunctions, functions)
  # Get a mapped parameter for each cell of the design
  pars <- mapped_pars(x, p_vector)
  if(is.null(data)){
    data <- vector("list", n_data_sim)
    # If no data is supplied generate some data sets
    for(i in 1:n_data_sim){
      data[[i]] <- make_data(p_vector, design = x, n_trials = 50)
    } # and bind them back together
    data <- do.call(rbind, data)
  }
  data <- data[!is.na(data$rt) & !is.infinite(data$rt),]
  data <- design_model(data, x, compress = FALSE, rt_resolution = 1e-15)

  if(is.null(x$model()$c_name)) stop("Current design type not supported for plotting")
  # if(x$model()$c_name == "LNR") stop("LNR designs not supported for plotting")
  type <- ifelse(x$model()$c_name == "DDM", "DDM", ifelse(x$model()$c_name == "LNR", "LNR", "race"))
  within_noise <- ifelse(x$model()$c_name == "LBA", FALSE, TRUE)
  # Split only relevant for DDM
  dots <- add_defaults(list(...), split = "R", within_noise = within_noise, plot_legend = TRUE)
  if(type != "DDM"){
    dots$split = NULL
    data <- data[data$winner,]
  }
  do.call(make_design_plot, c(list(data = data, pars = pars, factors = factors, main = dots$main,
                   plot_factor = plot_factor,
                   type = type), fix_dots(dots, make_design_plot)))
}


#' @exportS3Method
sampled_pars.default <- function(x,model=NULL,doMap=TRUE, add_da = FALSE, all_cells_dm = FALSE){
  if(is.null(x)) return(NULL)
  if(!is.null(attr(x, "custom_ll"))){
    pars <- numeric(length(attr(x,"sampled_p_names")))
    names(pars) <- attr(x,"sampled_p_names")
    return(pars)
  }
  if(!is.null(x$Ffactors)){
    x <- list(x)
    class(x) <- "emc.design"
  }
  out <- sampled_pars.emc.design(x, model = model, doMap = doMap, add_da = add_da, all_cells_dm = all_cells_dm)
  return(out)
}
