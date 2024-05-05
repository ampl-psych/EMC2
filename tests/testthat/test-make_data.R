# First create a design
design_LNR <- make_design(factors = list(S = c("left", "right"),
                                           subjects = 1:3),
                            Rlevels = c("left", "right"), model = LNR,
                            formula =list(m~0+S,s~1, t0~1),
                            constants=c(s=log(1)))
# Then create a p_vector:
p_vector <- c(m_Sleft=-2,m_Sright=2,t0=log(.2))

# Now we can simulate data
data <- make_data(p_vector, design_LNR, trials = 10)

test_that("make_data", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  expect_snapshot(
    make_data(p_vector, design_LNR, trials = 10)
  )
})
