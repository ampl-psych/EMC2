rm(list = ls())
library(BayesFactor)
n_subjects <- 5
n_cond <- 3
n_trials <- 10
n_groups <- 2

subject <- rep(1:n_subjects, each = n_trials*n_cond)
cond <- rep(rep(1:n_cond, n_subjects), each = n_trials)
rt <- rnorm(n_subjects * n_trials * n_cond, mean = cond*subject)

df <- data.frame(subjects = subject, cond = cond, rt = rt)
for(i in 2:n_groups) df <- rbind(df, df)
df$group <- rep(1:n_groups, each = n_trials*n_cond*n_subjects)


df$subjects <- as.factor(df$subjects)
df$group <- as.factor(df$group)

bf = lmBF(rt ~ cond* subjects*group, data = df, whichRandom=c("subjects", "subjects:cond", "subjects:group", "subjects:cond:group"), progress=FALSE)
X = model.matrix(bf)
X


data <- df
dataTypes <- c(subjects = "random", cond = "fixed", group = "fixed")
formula <- as.formula(rt ~ cond *group* subjects)
factors = BayesFactor:::fmlaFactors(formula, data)[-1]
nFactors = length(factors)


debug(BayesFactor:::createGMap)
if (nFactors == 0) {
  X = matrix(1, nrow(data), 1)
  gMap = c(intercept = NA)
} else {
  X = as.matrix(BayesFactor:::fullDesignMatrix(formula, data, dataTypes))
  gMap = BayesFactor:::createGMap(formula, data, dataTypes)
  X = cbind(1, X)
  gMap = c(intercept = NA, gMap)
}
attr(X, "gMap") = gMap


trms <- attr(terms(fmla, data = data), "term.labels")
sparse = any(dataTypes == "random")
Xs = lapply(trms, function(trm, data, dataTypes) {
  oneDesignMatrix(trm, data = data, dataTypes = dataTypes,
                  sparse = sparse)
}, data = data, dataTypes = dataTypes)
do.call("cbind", Xs)


terms(rt ~ cond * subjects, data = df)
