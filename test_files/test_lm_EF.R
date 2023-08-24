rm(list = ls())
library(BayesFactor)
n_subjects <- 5
n_cond <- 3
n_trials <- 3


subject <- rep(1:n_subjects, each = n_trials*n_cond)
cond <- rep(rep(1:n_cond, n_subjects), each = n_trials)
rt <- rnorm(n_subjects * n_trials * n_cond, mean = cond*subject)

df <- data.frame(subjects = subject, cond = cond, rt = rt)
df$subjects <- as.factor(df$subjects)
df$cond <- as.factor(df$cond)

bf = lmBF(rt ~ cond* subjects, data = df, whichRandom=c("subjects", "subjects:cond"), progress=FALSE)
X = model.matrix(bf)

terms(rt ~ cond * subjects, data = df)
