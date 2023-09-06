rm(list = ls())
library(lme4)
n_subjects <- 5
n_cond <- 3
n_trials <- 5
n_groups <- 2

subject <- rep(1:n_subjects, each = n_trials*n_cond)
cond <- rep(rep(1:n_cond, n_subjects), each = n_trials)
rt <- rnorm(n_subjects * n_trials * n_cond, mean = cond*subject)

df <- data.frame(subjects = subject, cond = cond, rt = rt)
for(i in 2:n_groups) df <- rbind(df, df)
df$group <- rep(1:n_groups, each = n_trials*n_cond*n_subjects)
df$subjects <- as.factor(df$subjects)
df$cond <- as.factor(df$cond)


types_and_baseForm <- function(formula){
  all_terms <- labels(terms(form))
  rnd_idx <- grepl("\\|", all_terms)
  all_factors <- all.vars(delete.response(terms(form)))
  terms_rnd <- all_terms[rnd_idx]
  types <- rep("fixed", length(all_factors))
  terms_out <- all_terms[!rnd_idx]
  types[!all_factors %in% all_terms] <- "random"
  if(any(rnd_idx)){
    for(term in terms_rnd){
      dep <- formula[[2]]
      b4_bar <- gsub("[|].*", "", term)
      b4_bar <- labels(terms(as.formula(paste0(dep, "~", b4_bar))))
      after_bar <- gsub(".*[|]", "", term) # Will break with nesting
      after_bar <- labels(terms(as.formula(paste0(dep, "~", after_bar))))
      terms_out <- c(terms_out, paste0(after_bar, ":", b4_bar))
      terms_out <- c(terms_out, after_bar)
    }
  }
  names(types) <- all_factors
  return(list(types = types, terms = terms_out))
}

form <- rt ~ cond + (cond|group:subjects) + (cond|group)

rowMultiply <- function (x, y)
{
  sparse = is(x, "sparseMatrix") | is(y, "sparseMatrix")
  if (nrow(x) != nrow(y))
    stop("Unequal row numbers in row.multiply:", nrow(x),
         ", ", nrow(y))
  K = sapply(1:nrow(x), function(n, x, y) {
    kronecker(x[n, ], y[n, ])
  }, x = x, y = y)
  K <- t(Matrix(as.vector(K), ncol = nrow(x), sparse = sparse))
  colnames(K) = as.vector(t(outer(colnames(x), colnames(y),
                                  function(x, y) {
                                    paste(x, y, sep = ".&.")
                                  })))
  return(K)
}

use <- types_and_baseForm(form)
trms <- use$terms
data <- df
dataTypes <- use$types
sparse <- TRUE
Xs = lapply(trms, function(trm, data, dataTypes) {
  BayesFactor:::oneDesignMatrix(trm, data = data, dataTypes = dataTypes,
                  sparse = sparse)
}, data = data, dataTypes = dataTypes)

m1 <- lmer(form, data)

getME(m1, "Z")


