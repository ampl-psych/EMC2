# PROVES: the location-move prior ratio coded in JointMoves::draw equals log N(tmu + delta d_j | mean, P^-1) - log N(tmu | mean, P^-1), and the group term is invariant under the joint shift. Needs only mvtnorm.
# PRINTED: worst relative error over 200 random cases = 3.14e-14; group term change under joint shift = -3.55e-15.
# Numerical check of the location-move log prior ratio implemented in JointMoves::draw:
#   lp = -delta * <P d_j, tmu - mean> - .5 delta^2 d_j' P d_j
# against log N(tmu + delta d_j | mean, P^-1) - log N(tmu | mean, P^-1), for random inputs.
set.seed(42)
worst <- 0
for (rep in 1:200) {
  p <- sample(2:8, 1)
  A <- matrix(rnorm(p * p), p); S0 <- crossprod(A) + diag(p) * .1   # prior covariance
  P <- solve(S0)
  mean <- rnorm(p); tmu <- rnorm(p, mean, 2)
  # directions as the code builds them: root * V * diag(1/sqrt(max(lambda,1))) with M = root' H root + I
  root <- t(chol(S0)); B <- matrix(rnorm(p * p), p); H <- crossprod(B)
  M <- t(root) %*% H %*% root; diag(M) <- diag(M) + 1
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  D <- root %*% e$vectors %*% diag(1 / sqrt(pmax(e$values, 1)), p)
  PD <- P %*% D; curv <- colSums(D * PD)
  j <- sample(p, 1); delta <- rnorm(1, 0, 1)
  lp_code <- -delta * sum(PD[, j] * (tmu - mean)) - .5 * delta^2 * curv[j]
  lp_true <- mvtnorm::dmvnorm(tmu + delta * D[, j], mean, S0, log = TRUE) - mvtnorm::dmvnorm(tmu, mean, S0, log = TRUE)
  worst <- max(worst, abs(lp_code - lp_true) / max(1, abs(lp_true)))
}
cat(sprintf("location-move prior ratio: worst relative error over 200 random cases = %.2e\n", worst))
# The group term N(alpha_s | subj_mu_s, Sigma) is translation invariant when alpha_s and subj_mu_s shift together:
p <- 4; A <- matrix(rnorm(16), 4); Sg <- crossprod(A); a <- rnorm(4); m <- rnorm(4); sh <- rnorm(4)
cat(sprintf("group term change under joint shift = %.2e (should be 0)\n",
  mvtnorm::dmvnorm(a + sh, m + sh, Sg, log = TRUE) - mvtnorm::dmvnorm(a, m, Sg, log = TRUE)))
