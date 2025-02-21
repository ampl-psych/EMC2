do_transform <- function(pars,transform)
  # pars is the parameter matrix, transform is a transform list
{
  isexp <- transform$func[colnames(pars)] == "exp"
  isprobit <- transform$func[colnames(pars)] == "pnorm"
  # Here instead of using isexp directly I use the column names in case v is replicated in pars (i.e. in map_mcmc)
  pars[,isexp] <- exp(sweep(pars[,isexp, drop = F], 2, transform$lower[colnames(pars)[isexp]], "-"))
  pars[,isprobit] <- pnorm(sweep(sweep(pars[,isprobit, drop = F], 2, transform$lower[colnames(pars)[isprobit]], "-")
                                 , 2, transform$upper[colnames(pars)[isprobit]]-transform$lower[colnames(pars)[isprobit]], "/"))
  pars
}

do_pre_transform <- function(p_vector,transform)
  # pars is the parameter matrix, transform is a transform list
{
  isexp <- transform$func[names(p_vector)] == "exp"
  isprobit <- transform$func[names(p_vector)] == "pnorm"
  # Here instead of using isexp directly I use the names in case v is replicated in p_vector (i.e. in map_mcmc)
  p_vector[isexp] <- exp(p_vector[isexp] - transform$lower[names(p_vector)[isexp]])
  p_vector[isprobit] <- pnorm((p_vector[isprobit] - transform$lower[names(p_vector)[isprobit]])/
                                (transform$upper[names(p_vector)[isprobit]]-transform$lower[names(p_vector)[isprobit]]))
  p_vector
}


#### Functions to look at parameters ----

#### Functions to look at parameters ----

map_p <- function(p,dadm)
  # Map p to dadm and returns matrix of mapped parameters
  # p is either a vector or a matrix (ncol = number of subjects) of p_vectors
  # dadm is a design matrix with attributes containing model information
{

  # Check if p is a matrix and validate column names match parameter names
  if ( is.matrix(p) ) {
    if (!all(sort(dimnames(p)[[2]])==sort(attr(dadm,"p_names"))))
      stop("p col.names must be: ",paste(attr(dadm,"p_names"),collapse=", "))
    if (!all(levels(dadm$subjects) %in% dimnames(p)[[1]]))
      stop("p must have rows named for every subject in dadm")
    p <- p[dadm$subjects,]
  } else if (!all(sort(names(p))==sort(attr(dadm,"p_names")))) # If p is vector, check names
    stop("p names must be: ",paste(attr(dadm,"p_names"),collapse=", "))

  # Get parameter names from model and create output matrix
  do_p <- names(attr(dadm,"model")()$p_types)
  pars <- matrix(nrow=nrow(dadm),ncol=length(do_p),dimnames=list(NULL,do_p))

  # If there are any trends do these first, they might be used later in mapping
  # Otherwise we're not applying the trend premap, but we are doing it pre-transform
  # So these trend parameters are post-map, pre-transform and have to be included in the pars output
  premap_idx <- rep(F, length(do_p))
  if(!is.null(attr(dadm,"model")()$trend) &&
     (attr(attr(dadm,"model")()$trend, "premap") || attr(attr(dadm,"model")()$trend, "pretransform"))){
    trend_names <- get_trend_pnames(attr(dadm,"model")()$trend)
    pretrend_idx <- do_p %in% trend_names
    if((attr(attr(dadm,"model")()$trend, "premap"))){
      # These can be removed from the pars matrix at the end
      # Since they are already used before the mapping
      premap_idx <- pretrend_idx
    }
    # Reorder parameters to make design matrix for trends first
    do_p <- c(do_p[pretrend_idx], do_p[!pretrend_idx])
  } else{
    pretrend_idx <- rep(F, length(do_p))
  }
  k <- 1
  # Loop through each parameter
  for (i in do_p) {
    cur_design <- attr(dadm,"designs")[[i]]
    # Handle vector vs matrix input differently
    if ( !is.matrix(p) ) {
      pm <- t(as.matrix(p[colnames(cur_design)]))
      pm <- pm[rep(1,nrow(pars)),,drop=FALSE]
    } else pm <- p[,colnames(cur_design),drop=FALSE]

    # Apply pre-mapped trends if they exist
    if (!is.null(attr(dadm,"model")()$trend) && attr(attr(dadm,"model")()$trend, "premap")) {
      trend <- attr(dadm,"model")()$trend
      isin <- names(trend) %in% colnames(pm)
      if (any(isin)){ # At this point the trend has already been mapped and transformed
        for (j in names(trend)[isin]) {
          cur_trend <- trend[[j]]
          # We can select the trend pars from the already update pars matrix
          trend_pars <- pars[,cur_trend$trend_pnames]
          pm[,j] <- run_trend(dadm, cur_trend, pm[,j], trend_pars)
        }
      }
    }

    # Apply design matrix and sum parameter effects
    tmp <- pm*cur_design[attr(cur_design,"expand"),,drop=FALSE]
    tmp[is.nan(tmp)] <- 0 # Handle 0 weight x Inf parameter cases
    tmp <- apply(tmp,1,sum)
    # If this is a premap trend parameter, transform it here already
    # We'll need it transformed later in this loop (for trending other parameters)
    if(k <= sum(pretrend_idx)){
      tmp <- as.matrix(tmp)
      colnames(tmp) <- i
      tmp <- do_transform(tmp, attr(dadm,"model")()$transform)
    }
    k <- k + 1
    pars[,i] <- tmp
  }
  # Return only non-trend parameters
  return(pars[,!premap_idx])
}


get_pars_matrix <- function(p_vector,dadm) {
  # Order:
  # 1 pretransform
  # 2 add constants
  # 3 map
  # # - if premap trend:
  # #   First make trend pars matrix
  # #   Apply trends premap and remove trend pars from pars matrix
  # # - map
  # 4 if pretransform trend:
  # # - apply trends and remove trend pars from pars matrix
  # 5 transform
  # 6 if posttransform trend:
  # # - apply trends and remove trend pars from pars matrix
  # 7 trial-wise transform
  # 8 bound

  # Niek should constants be included in pre_transform? I think not?
  p_vector <- do_pre_transform(p_vector, attr(dadm, "model")()$pre_transform)
  # If there's any premap trends, they're done in map_p
  pars <- map_p(add_constants(p_vector,attr(dadm,"constants")),dadm)
  if(!is.null(attr(dadm, "model")()$trend) && attr(attr(dadm, "model")()$trend, "pretransform")){
    # This runs the trend and afterwards removes the trend parameters
    pars <- prep_trend(dadm, attr(dadm, "model")()$trend, pars)
  }
  pars <- do_transform(pars, attr(dadm,"model")()$transform)
  if(!is.null(attr(dadm, "model")()$trend) && attr(attr(dadm, "model")()$trend, "posttransform")){
    # This runs the trend and afterwards removes the trend parameters
    pars <- prep_trend(dadm, attr(dadm, "model")()$trend, pars)
  }
  pars <- attr(dadm,"model")()$Ttransform(pars, dadm)
  pars <- add_bound(pars, attr(dadm,"model")()$bound)
  return(pars)
}


make_pmat <- function(p_vector,design)
  # puts vector form of p_vector into matrix form
{
  ss <- design$Ffactors$subjects
  matrix(rep(p_vector,each=length(ss)),nrow=length(ss),
         dimnames=list(ss,names(p_vector)))
}

add_constants <- function(p,constants)
{
  if (is.null(constants)) return(p)
  if (is.matrix(p)) {
    nams <- c(dimnames(p)[[2]],names(constants))
    p <- cbind(p,matrix(rep(constants,each=dim(p)[1]),nrow=dim(p)[1]))
    dimnames(p)[[2]] <- nams
    return(p)
  } else{
    return(c(p,constants))
  }

}

add_recalculated_pars <- function(pmat, model, cnams){
  modifiers <- unlist(lapply(strsplit(cnams,"_"),function(x){paste0(x[-1], collapse = "_")}))
  par_names <- colnames(pmat)
  unq_pars <- unique(par_names)
  par_table <- table(par_names)
  counts <- lapply(par_table, function(x) return(1:x))
  combn <- do.call(expand.grid, counts)
  colnames(combn) <- names(par_table)

  out <- list()
  modfs <- list()
  for(r in 1:nrow(combn)){
    pmat_in <- matrix(NA, nrow = nrow(pmat), ncol = length(unq_pars))
    colnames(pmat_in) <- unq_pars
    cur_modifiers <- setNames(numeric(length(unq_pars)), unq_pars)
    for(par in unq_pars){
      pmat_in[,par] <- pmat[,which(colnames(pmat) == par)[combn[r,par]]]
      cur_modifiers[par] <- modifiers[which(colnames(pmat) == par)[combn[r,par]]]
    }
    added <- model()$Ttransform(pmat_in)
    added <- added[, !(colnames(added) %in% colnames(pmat_in)), drop = F]
    attr(added, "ok") <- NULL
    out[[r]] <- added
    modfs[[r]] <- cur_modifiers
  }
  m_out <- matrix(0, nrow = nrow(pmat), ncol = 0)
  if(ncol(out[[1]]) == 0) return(NULL)
  for(i in 1:ncol(out[[1]])){
    cur_par <- lapply(out, function(x) x[,i, drop = F])
    not_dups <- !duplicated(cur_par)
    cur_combn <- combn[not_dups,]
    pars_vary <- colnames(cur_combn)[colMeans(cur_combn) != 1]
    to_add <- do.call(cbind, cur_par[not_dups])
    cur_modfs <- modfs[not_dups]
    fnams <- sapply(cur_modfs, function(x) paste0(unique(unlist(strsplit(x[pars_vary], split = "_"))), collapse = "_"))
    if(fnams != "") colnames(to_add) <- paste0(colnames(to_add)[1], "_", fnams)
    m_out <- cbind(m_out, to_add)
  }
  return(m_out)
}

map_mcmc <- function(mcmc,design,include_constants = TRUE, add_recalculated = FALSE, covariates = NULL)
  # Maps vector or matrix (usually mcmc object) of sampled parameters to native
  # model parameterization.
{
  if(is.null(design$model)){
    design <- design[[1]]
  }
  model <- design$model

  doMap <- function(mapi,pmat, covariates = NULL){
    if(!is.null(covariates)){
      if(nrow(covariates) != nrow(pmat)) covariates <- covariates[sample(1:nrow(covariates), nrow(pmat), replace = T),, drop = F]
      pmat <- cbind(pmat, covariates)
    }
    t(mapi %*% t(pmat[,dimnames(mapi)[[2]],drop=FALSE]))
  }

  get_p_types <- function(nams, reverse = FALSE){
    if(reverse){
      out <- unlist(lapply(strsplit(nams,"_"),function(x){x[[-1]]}))
    } else{
      out <- unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))
    }
    return(out)
  }


  map <- attr(sampled_pars(design, add_da = TRUE, all_cells_dm = TRUE),"map")

  constants <- design$constants
  if (!is.matrix(mcmc) & !is.array(mcmc)) mcmc <- t(as.matrix(mcmc))
  if(length(dim(mcmc)) == 2){
    is_matrix <- TRUE
    mcmc_array <- array(mcmc, dim = c(nrow(mcmc), 1, ncol(mcmc)))
    rownames(mcmc_array) <- rownames(mcmc)
  } else{
    mcmc_array <- mcmc
    is_matrix <- FALSE
  }
  mp <- mapped_pars(design,mcmc_array[,1,1],remove_RACE=FALSE, covariates = covariates)

  for(k in 1:ncol(mcmc_array)){
    mcmc <- t(mcmc_array[,k,])
    mcmc <- t(apply(mcmc, 1, do_pre_transform, design$model()$pre_transform))
    pmat <- add_constants(mcmc,constants)
    plist <- lapply(map,doMap,pmat=pmat, covariates)
    # Give mapped variables names and flag constant
    for (i in 1:length(plist)) {
      vars <- row.names(attr(terms(design$Flist[[i]]),"factors"))
      if(any(names(design$Fcovariates) %in% vars)){
        cov_vars <- vars[(vars %in% names(design$Fcovariates))]
        vars <- vars[!(vars %in% names(design$Fcovariates))]
        has_cov <- TRUE
      } else{
        has_cov <- FALSE
      }
      uniq <- !duplicated(apply(mp[,vars, drop = F],1,paste,collapse="_"))
      if (is.null(vars)) {
        colnames(plist[[i]]) <- names(plist)[i]
      }else {
        if(length(vars) == 1) colnames(plist[[i]]) <- vars
        else if(length(vars) == 2){
          colnames(plist[[i]]) <- paste(vars[1], apply(matrix(do.call(cbind, Map(paste0, vars[-1], mp[uniq, vars[-1]])))
                                                       ,1, paste0, collapse = "_"), sep = "_")
        } else{
          colnames(plist[[i]]) <- paste(vars[1], apply(do.call(cbind, Map(paste0, vars[-1], mp[uniq, vars[-1]]))
                                                       ,1, paste0, collapse = "_"), sep = "_")
        }
      }
      if(has_cov) colnames(plist[[i]]) <- paste(colnames(plist[[i]]), cov_vars, sep = "_")
      # if (is.matrix(plist[[i]])) isConstant <- c(isConstant, apply(plist[[i]],2,function(x){all(x[1]==x[-1])}))
    }
    pmat <- do.call(cbind,plist)
    cnams <- colnames(pmat)
    colnames(pmat) <- get_p_types(cnams)
    pmat[,] <- do_transform(pmat, model()$transform)[,1:length(cnams)]

    if(add_recalculated) extras <- add_recalculated_pars(pmat, model, cnams)
    colnames(pmat) <- cnams
    if(add_recalculated) pmat <- cbind(pmat, extras)
    is_constant <- apply(pmat,2,function(x){all(duplicated(x)[-1])})
    if(!include_constants){
      pmat <- pmat[,!is_constant, drop = F]
    }
    if(k == 1){
      out <- array(0, dim = c(ncol(pmat), ncol(mcmc_array), dim(mcmc_array)[3]))
      rownames(out) <- colnames(pmat)
      colnames(out) <- colnames(mcmc_array)
    }
    out[,k,] <- t(pmat)
  }
  if(is_matrix){
    out <- out[,1,]
  }
  attr(out,"isConstant") <- is_constant
  return(out)
}



fill_transform <- function(transform, model, p_vector,
                           supported=c("identity","exp","pnorm"),
                           has_lower=c("exp","pnorm"),has_upper=c("pnorm"),
                           is_pre = FALSE){
  if(!is.null(transform)){
    if (!all(transform$func %in% supported)){
      stop("Only ", paste(supported, collapse = ", "), " transforms supported")
    }
    if(!is_pre){
      if (!all(names(transform$func) %in% names(model()$p_types)))stop("transform on parameter not in the model p_types")
      if (!all(names(transform$lower) %in% names(model()$p_types)))stop("transform on parameter not in the model p_types")
      if (!all(names(transform$upper) %in% names(model()$p_types)))stop("transform on parameter not in the model p_types")
    } else{
      if (!all(names(transform$func) %in% names(p_vector))) stop("pre_transform on parameter not in the sampled_pars")
      if (!all(names(transform$lower) %in% names(p_vector))) stop("transform on parameter not in the model p_types")
      if (!all(names(transform$upper) %in% names(p_vector))) stop("transform on parameter not in the model p_types")
    }
  }
  model_list <- model()
  p_names <- names(model_list$p_types)
  if(is_pre){
    p_names <- names(p_vector)
  }
  filled_func <- filled_lower <- filled_upper <- setNames(rep(NA, length(p_names)), p_names)
  # First update from the model
  if(is_pre){
    filled_func[names(model_list$pre_transform$func)] <- model_list$pre_transform$func
    filled_lower[names(model_list$pre_transform$lower)] <- model_list$pre_transform$lower
    filled_upper[names(model_list$pre_transform$upper)] <- model_list$pre_transform$upper
  } else{
    filled_func[names(model_list$transform$func)] <- model_list$transform$func
    filled_lower[names(model_list$transform$lower)] <- model_list$transform$lower
    filled_upper[names(model_list$transform$upper)] <- model_list$transform$upper
  }
  # Now updates from the user transform (they override the model)
  filled_func[names(transform$func)] <- transform$func
  filled_lower[names(transform$lower)] <- transform$lower
  filled_upper[names(transform$upper)] <- transform$upper

  # Now fill in remainers
  filled_func[is.na(filled_func)] <- "identity"
  filled_lower[is.na(filled_lower) & (filled_func %in% has_lower)] <- 0
  filled_upper[is.na(filled_upper) & (filled_func %in% has_upper)] <- 1
  # Remainers must be identity and have no bounds
  filled_lower[is.na(filled_lower)] <- -Inf
  filled_upper[is.na(filled_upper)] <- Inf
  return(list(func=filled_func,lower=filled_lower,upper=filled_upper))
}

fill_bound <- function(bound, model) {
  filled_bound <- model()$bound
  if (!is.null(bound)) {
    if (names(bound)[1] != "minmax")
      stop("first entry of bound must be named minmax")
    if (!all(colnames(bound$minmax) %in% names(model()$p_types)))
      stop("minmax column names must correspond to parameter types")
    if (!is.null(bound$exception) &&
        (!all(names(bound$exception) %in% names(model()$p_types))))
      stop("exception names must correspond to parameter types")
    filled_bound$minmax[,colnames(bound$minmax)] <- bound$minmax
    if (!is.null(bound$exception)) {
      filled_bound$exception <- c(bound$exception,filled_bound$exception)
      filled_bound$exception <- filled_bound$exception[!duplicated(names(filled_bound$exception))]
    }
  }
  filled_bound
}

# This form used in random number generation
do_bound <- function(pars,bound) {
  tpars <- t(pars[,colnames(bound$minmax),drop=FALSE])
  ok <- tpars > bound$minmax[1,] & tpars < bound$minmax[2,]
  if (!is.null(bound$exception)) ok[names(bound$exception),] <-
    ok[names(bound$exception),] |
    (tpars[names(bound$exception),] == bound$exception)
  apply(ok,2,all)
}

# This form used in get_pars
add_bound <- function(pars,bound) {
  attr(pars,"ok") <- do_bound(pars,bound)
  pars
}

generate_design_equations <- function(design_matrix,
                                      factor_cols = NULL,
                                      numeric_cols = NULL,
                                      transform) {
  # 1. If user hasn't specified which columns are factors or numeric, guess:
  if (is.null(factor_cols)) {
    # We'll assume anything that is a factor or character is a "factor column"
    factor_cols <- names(design_matrix)[
      sapply(design_matrix, function(x) is.factor(x) || is.character(x))
    ]
  }
  if (is.null(numeric_cols)) {
    # We'll assume everything that is numeric is for the design (contrast) columns
    numeric_cols <- names(design_matrix)[
      sapply(design_matrix, is.numeric)
    ]
  }

  # 2. Build the "equation" strings from numeric cols
  make_equation_string <- function(row_i) {
    eq_terms <- c()
    for (colname in numeric_cols) {
      val <- as.numeric(row_i[[colname]])
      # Skip zeros
      if (abs(val) < 1e-15) next

      # If val is exactly +1 or -1, skip numeric part:
      if (abs(val - 1) < 1e-15) {
        # val == +1
        term_str <- paste0("+ ", colname)
      } else if (abs(val + 1) < 1e-15) {
        # val == -1
        term_str <- paste0("- ", colname)
      } else {
        # Some other numeric coefficient => e.g. "+ 0.5 * col"
        sign_str <- ifelse(val >= 0, "+ ", "- ")
        term_str <- paste0(sign_str, format(abs(val), digits = 3), " * ", colname)
      }
      eq_terms <- c(eq_terms, term_str)
    }

    if (length(eq_terms) == 0) {
      # If everything was zero
      return("0")
    }

    # Combine eq_terms
    eq_string <- paste(eq_terms, collapse = " ")
    # Remove the leading '+ ' if present
    sub("^\\+ ", "", eq_string)
  }

  # 3. Compute the equation string for each row
  eq_strings <- apply(design_matrix, 1, make_equation_string)

  # 4. Determine the widths needed for each factor column
  #    so everything lines up in columns nicely.
  col_widths <- sapply(factor_cols, function(fc) {
    # The width must accommodate either the column name or the largest factor level
    max(nchar(fc), max(nchar(as.character(design_matrix[[fc]])), na.rm = TRUE))
  })

  # We'll label the final column as "Equation"
  eq_header <- ""
  # We don't necessarily need to align the entire equation column itself
  # (since it may vary in length), but let's keep a short label for the header.

  # 5. Print the header row
  #    E.g. if factor_cols = c("S", "E"), we do:
  #    " S     E    : Equation"
  factor_header_str <- paste(mapply(function(fc, w) {
    sprintf("%-*s", w, fc)   # left-justify the column name
  }, factor_cols, col_widths),
  collapse = "  ")

  cat("  ", factor_header_str, "  ", eq_header, "\n", sep="")

  # Optionally, add a simple "rule" line or blank line:
  # (Uncomment if you like to separate header from rows)
  # cat(" ", paste(rep("-", nchar(factor_header_str) + nchar(eq_header) + 5),
  #               collapse=""), "\n", sep="")

  # 6. Print each row: factor columns in their spaces, then ": <equation>"
  for (i in seq_len(nrow(design_matrix))) {
    row_factor_str <- paste(mapply(function(fc, w) {
      # Each factor value is left-justified in the same width as the header
      sprintf("%-*s", w, as.character(design_matrix[[fc]][i]))
    }, factor_cols, col_widths),
    collapse = "  ")
    if(transform == "identity"){
      cat(" ", row_factor_str, "  : ", eq_strings[i], "\n", sep="")
    } else{
      cat(" ", row_factor_str, "  : ", transform, "(", eq_strings[i], ")", "\n", sep="")
    }

  }
}


verbal_dm <- function(design){
  map <- attr(sampled_pars(design,design$model, add_da = TRUE), "map")
  map_no_da <- attr(sampled_pars(design,design$model), "map")
  transforms <- design$model()$transform$func
  for(i in 1:length(map)){
    m <- map[[i]]
    if(ncol(m) == 1) next
    cat(paste0("$", names(map)[i]), "\n")
    mnd <- map_no_da[[i]]
    par_idx <- colnames(m) %in% colnames(mnd)
    generate_design_equations(m, colnames(m)[!par_idx], colnames(m)[par_idx],
                              transform = transforms[names(map)[i]])
    cat("\n")
  }
}




