# Create subject-level design
subj_design <- design(data = forstmann, model = DDM,
                      formula = list(v ~ S, a ~ E, t0 ~ 1),
                      contrasts = list(S = contr.helmert))
# Add some age covariate and roughly demean
# Demeaning is important to ensure that the interpretation of the group-level intercept
# is the mean of the group (i.e., 'mu' still represents the group-level mean)
forstmann$age <- as.numeric(forstmann$subjects) -mean(as.numeric(forstmann$subjects))
# Create fake group column
forstmann$group <- ifelse(forstmann$subjects %in%
                            unique(forstmann$subjects)[seq(1, 19, 2)], "A", "B")

# Create group-level design matrices
group_des <- group_design(
  formula = list(v_S1 ~ age + group, a ~ age),
  data = forstmann,
  subject_design = subj_design,
  contrasts = list(group = contr.bayes)
)

pri <- prior(subj_design, group_design = group_des)
emc <- make_emc(forstmann, subj_design, group_design = group_des, prior_list = pri)

test_that("groups", {
  expect_snapshot(print(group_des))
  expect_snapshot(summary(group_des))
  expect_snapshot(credint(pri, group_design = group_des))
  expect_snapshot(credint(pri, group_design = group_des, selection = "beta"))
  expect_snapshot(init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples)
})
