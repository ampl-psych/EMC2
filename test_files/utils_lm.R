
# Adaptations from BayesFactor and other packages ------------------------------------

rowMultiply <- function(x, y){
  # This function takes two design matrices and multiplies them for interactions
  # As described in Rouder 2012
  sparse = is(x, "sparseMatrix") | is(y, "sparseMatrix")
  if (nrow(x) != nrow(y))
    stop("Unequal row numbers in row.multiply:", nrow(x),
         ", ", nrow(y))
  K <- Matrix::t(Matrix::KhatriRao(Matrix::t(x), Matrix::t(y))) # This is way quicker than the original BayesFactor way
  colnames(K) = as.vector(Matrix::t(outer(colnames(x), colnames(y),
                                  function(x, y) {
                                    paste(x, y, sep = "_")
                                  })))
  return(K)
}

contr.bayes <- function(n, contrasts = TRUE) {
  # This makes a balanced contrasts that assigns equal weight to a fixed effect factor
  if (length(n) <= 1L) {
    if (is.numeric(n) && length(n) == 1L && n > 1L)
      TRUE
    else stop("not enough degrees of freedom to define contrasts")
  } else n <- length(n)
  cont <- diag(n)
  if (contrasts) {
    a <- n
    I_a <- diag(a)
    J_a <- matrix(1, nrow = a, ncol = a)
    Sigma_a <- I_a - J_a/a
    cont <- eigen(Sigma_a)$vectors[,seq_len(a-1), drop = FALSE]
  }
  cont
}

model.Matrix <- function(object, data, contrasts.arg = NULL, drop.unused.levels = FALSE,
                         xlev = NULL){
  # Creates a sparse design matrix
  options(na.action='na.pass')
  m <- Matrix::sparse.model.matrix(object, data = data, contrasts.arg = contrasts.arg,
                                   drop.unused.levels = drop.unused.levels, xlev = xlev)
  new("dsparseModelMatrix", m, assign = attr(m, "assign"),
      contrasts = if (is.null(ctr <- attr(m, "contrasts")))
        list()
      else ctr)
  m[is.na(m)] <- 0
  return(m)
}

oneDM <- function (trm, data, dataTypes)
{
  # First figures out what type of effect we're dealing with
  # Then depending on the effect fixed (continuous or not) or random makes a design matrix
  # In case of interactions, first makes the DM for them separately then combines them using
  # Rowmultiply function.
  effects <- strsplit(trm, ":")[[1]]
  if (length(effects) == 1) {
    if (dataTypes[effects] == "fixed") {
      fmla = paste("~", effects)
      if(is.factor(data[,effects])){
        contrast.args <- list()
        contrast.args[[effects]] <- contr.bayes
        X = model.Matrix(formula(fmla), data = data, contrasts.arg = contrast.args)
      } else{
        X = model.Matrix(formula(fmla), data = data)
      }

      # colnames(X) = paste(effects, "_redu_", 1:ncol(X),
      #                     sep = "")
      X <- X[,-1, drop = F] # Remove intercept
    }
    else if (dataTypes[effects] == "random") {
      if(!is.factor(data[effects])) "Stop only factors can be random variables"
      fmla = paste("~", effects, "-1")
      X = model.Matrix(formula(fmla), data = data,)
    }
    return(X)
  }
  else {
    Xs = lapply(effects, function(trm, data, dataTypes) {
      oneDM(trm, data = data, dataTypes = dataTypes)
    }, data = data, dataTypes = dataTypes)
    X = Reduce(rowMultiply, x = Xs)
    X <- X[,Matrix::colSums(abs(X)) != 0, drop = F]
    return(X)
  }
}


# Custom code -------------------------------------------------------------

types_and_baseForm <- function(formula){
  # Based on a formula, finds out what main and interaction effects there are
  # And if they are fixed or random
  all_terms <- labels(terms(formula))
  rnd_idx <- grepl("\\|", all_terms)
  all_factors <- all.vars(delete.response(terms(formula)))
  terms_rnd <- all_terms[rnd_idx]
  types <- rep("fixed", length(all_factors))
  terms_out <- all_terms[!rnd_idx]
  types[!all_factors %in% all_terms] <- "random"
  if(any(rnd_idx)){
    # Random slopes can be written as interactions with the thing they slope on
    # Thus the bar can be interpreted as an interaction, thus we replace the bar with that
    for(term in terms_rnd){
      dep <- formula[[2]]
      b4_bar <- gsub("[|].*", "", term)
      b4_bar <- labels(terms(as.formula(paste0(dep, "~", b4_bar))))
      after_bar <- gsub(".*[|]", "", term)
      after_bar <- labels(terms(as.formula(paste0(dep, "~", after_bar))))
      if(length(b4_bar) > 0){
        for(b4 in b4_bar) terms_out <- c(terms_out, paste0(after_bar, ":", b4))
      }
      terms_out <- c(terms_out, after_bar)
    }
  }
  names(types) <- all_factors
  if(any(duplicated(terms_out))){
    warning("duplicate terms found, this usually means you specified
            redundant terms e.g. A + B + A*B")
  }
  return(list(types = types, terms = terms_out))
}

sepFixRandom <- function(Xs, types){
  # Separates the DMs into fixed and random effects
  fixed <- as.matrix(rep(1, nrow(Xs[[1]])))
  random <- NULL
  random_names <- names(types[types == "random"])
  for(X in Xs){
    is_random <- F
    for(name in random_names){
      # If any of the effects (even if it's an interaction) are random
      # The whole thing is treated as random for EMC2
      # Only for ease of use, mathematically they are mixed interactions stil
      if(any(grepl(name, colnames(X)))){
        is_random <- T
      }
    }
    if(!is_random){
      fixed <- cbind(fixed, X)
    } else{
      random <- cbind(random, X)
    }
  }
  colnames(fixed)[1] <- ""
  return(list(fixed = fixed, random = random))
}

get_parNames <- function(DM_fixed, DM_random, parameter, group_random_by = "subject"){
  # First find what are the fixed parameters
  # Then find out what are the parameters that are allowed to vary by subject
  # Then do some parsing to make a nice overview
  fixed_pars <- colnames(DM_fixed)
  to_add <- fixed_pars != ""
  fixed_pars[to_add] <- paste0(parameter, "_", fixed_pars[to_add])
  fixed_pars[!to_add] <- parameter
  random_pars <- colnames(DM_random)
  if(!is.null(random_pars)){
    seps <- strsplit(random_pars, "_")
    out <- numeric(length(seps))
    for(i in 1:length(seps)){
      bad_idx <- grepl(group_random_by, seps[[i]])
      if(all(bad_idx)){
        out[i] <- parameter
      } else{
        out[i] <- paste0(parameter, "_", paste0(seps[[i]][!bad_idx], collapse = "_"))
      }
    }
    uniques_in_group <- unique(out)
  } else{
    uniques_in_group <- NULL
  }

  return(list(random = paste0(parameter, "_", random_pars), fixed = fixed_pars,
              uniques_in_group = uniques_in_group))
}

compress_dadm_lm <- function(da,fixed_DM,random_DM, Fcov = NULL)
  # out keeps only unique rows in terms of all parameters design matrices
  # R, lR and rt (at given resolution) from full data set
{
  DM <- mapply(cbind, fixed_DM, random_DM)
  nacc <- length(unique(da$lR))
  # contract output
  cells <- paste(apply(do.call(cbind,lapply(DM,function(x){
    apply(x,1,paste,collapse="_")})
  ),1,paste,collapse="+"),da$subjects,da$R,da$lR,da$rt,sep="+")
  # Make sure that if row is included for a trial so are other rows
  if (!is.null(Fcov))
    cells <- paste(cells,apply(da[,names(Fcov),drop=FALSE],1,paste,collapse="+"),sep="+")
  if (nacc>1) cells <- paste0(rep(apply(matrix(cells,nrow=nacc),2,paste0,collapse="_"),
                                  each=nacc),rep(1:nacc,times=length(cells)/nacc),sep="_")
  contract <- !duplicated(cells)
  out <- da[contract,,drop=FALSE]
  fixed_DM <- lapply(fixed_DM, FUN = function(x) x[contract,,drop = F])
  random_DM <- lapply(random_DM, FUN = function(x) x[contract,,drop = F])
  attr(out,"contract") <- contract
  attr(out,"expand") <- as.numeric(factor(cells,levels=unique(cells)))
  lR1 <- da$lR==levels(da$lR)[[1]]
  attr(out,"expand_winner") <- as.numeric(factor(cells[lR1],levels=unique(cells[lR1])))
  attr(out,"s_expand") <- da$subjects

  # indices to use to contract further ignoring rt then expand back
  cells_nort <- paste(
    apply(do.call(cbind,lapply(DM,function(x){
      apply(x,1,paste,collapse="_")})
    ),1,paste,collapse="+"),da$subjects,da$R,da$lR,sep="+")[contract]
  attr(out,"unique_nort") <- !duplicated(cells_nort)
  # Only first level WHY????
  cells <- cells[da$lR==levels(da$lR)[1]]
  cells_nort <- cells_nort[out$lR==levels(out$lR)[1]]
  attr(out,"expand_nort") <- as.numeric(factor(cells_nort,
                                               levels=unique(cells_nort)))[as.numeric(factor(cells,levels=unique(cells)))]

  # indices to use to contract ignoring rt and response (R), then expand back
  cells_nortR <- paste(apply(do.call(cbind,lapply(DM,function(x){
    apply(x,1,paste,collapse="_")})),1,paste,collapse="+"),
    da$subjects,da$lR,sep="+")[contract]
  attr(out,"unique_nortR") <- !duplicated(cells_nortR)
  # Only first level WHY????
  cells_nortR <- cells_nortR[out$lR==levels(out$lR)[1]]
  attr(out,"expand_nortR") <- as.numeric(factor(cells_nortR,
                                                levels=unique(cells_nortR)))[as.numeric(factor(cells,levels=unique(cells)))]
  attr(out, "DM_fixed") <- fixed_DM
  attr(out, "DM_random") <- random_DM
  out
}

getGMap <- function(DM, use, data, constants){
  names_fixed <- colnames(DM$fixed)
  names_random <- colnames(DM$random)
  names_fixed <- names_fixed[!names_fixed %in% names(constants)]
  names_random <- names_random[!names_random %in% names(constants)]

  names_fixed <- gsub('[[:digit:]]+', '', names_fixed)
  # names_random <- gsub('[[:digit:]]+', '', names_random)

  g_fixed <- as.numeric(as.factor(names_fixed))

  names(g_fixed) <- names_fixed
  g_random <- get_GRandom(names_random, data, use)
  return(list(fixed = g_fixed, random = g_random))
}

get_GRandom <- function(names, data, use){
  # First order the random effects by order of appearance
  random_types <- names(use$types[use$types == "random"])
  n_random <- length(random_types)
  levs <- numeric(n_random)
  for(k in 1:n_random){
    levs[k] <- length(levels(data[,random_types[k]]))
  }
  names(levs) <- random_types
  levs <- names(levs[order(levs, decreasing = T)])
  seps <- strsplit(names, "_")
  out <- length(seps)
  for(i in 1:length(seps)){
    sep <- seps[[i]]
    idx <- unlist(lapply(levs, FUN = function(x, sep) return(any(grepl(x, sep))), sep))
    to_remove <- levs[idx][1]
    out[i] <- paste0(sep[!grepl(to_remove, sep)], collapse = "_")
  }
  g_random <- as.numeric(as.factor(out))
  names(g_random) <- out
  return(g_random)
}

oneG <- function(G){
  g_out <- numeric()
  for(g in G){
    if(length(g_out) == 0){
      g_out <- g
    } else{
      g_out <- c(g_out, g + max(g_out))
    }
  }
  names(g_out) <- unlist(lapply(G, names))
  return(g_out)
}

get_sub_idx <- function(s, par_names){
  seps <- strsplit(par_names, "_")
  idx <- unlist(lapply(seps, FUN = function(x) any(x == paste0("subjects", s))))
  return(idx)
}


