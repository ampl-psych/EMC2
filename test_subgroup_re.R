source("R/group_design.R")
source("R/variant_standard.R")
source("R/get_type_objects.R")
source("R/sampling.R") # For run_hyper

library(mvtnorm)
library(MASS)

# 1. Define Ground Truth
# 1. Define Ground Truth
n_subjects <- 1000
n_groups <- 20
groups <- factor(rep(paste0("G", 1:n_groups), each = n_subjects / n_groups))
# Let's say we have 2 parameters: v, a
beta_true <- c(v = 1.0, a = 0.8) # Intercepts
# Random effects u (one per group)
# Model: v ~ 1 + (1|group)
# a ~ 1
# u_v ~ N(0, s_v^2)
s_v_true <- 0.3
u_v_true <- rnorm(n_groups, 0, s_v_true)
names(u_v_true) <- paste0("G", 1:n_groups)
Sigma_true <- diag(c(0.1, 0.1)) # Subject variance

# 2. Simulate Subjects (Alphas)
# alpha_i = beta + Z_i * u + epsilon_i
alphas <- matrix(NA, nrow = 2, ncol = n_subjects)
rownames(alphas) <- c("v", "a")
colnames(alphas) <- 1:n_subjects

for (i in 1:n_subjects) {
    g <- groups[i]
    # Fixed
    mu <- beta_true
    # Random offset for v
    # u_v has 2 elements.
    g_idx <- as.integer(g)
    mu["v"] <- mu["v"] + u_v_true[g_idx]
    # a has no random effect

    alphas[, i] <- mu + rmvnorm(1, sigma = Sigma_true)
}

# 3. Setup Design
design <- list(
    v ~ 1 + (1 | group),
    a ~ 1
)

# 4. Run Hyper
# We need to construct a "sampler" object structure expected by run_hyper?
# run_hyper signature: function(data_list, design, model, ...)
# Actually run_hyper takes "emc" objects?
# Let's check run_hyper in R/sampling.R or R/model_SS.R?
# I'll check usage in test_group_design.R

# Setup data frame for design
data <- data.frame(subjects = 1:n_subjects, group = groups)

# Create a mock "sampler" object or just pass alphas to run_hyper?
# run_hyper usually takes "samples" object.
# Let's use internal function `fit_hyper` if `run_hyper` is high level.
# Or looking at `run_hyper` code:
# It iterates stages.
# It expects an `emc` object.

# Let's rely on standard make_emc if possible, but we don't have full model here.
# Just testing Gibbs.
# `run_hyper` usually takes an initialized object.

# Let's manually invoke the Gibbs step in a loop.
# Initialize
n_iter <- 500
n_pars <- 2
par_names <- c("v", "a")
# Mocking infrastructure
sampled_pars <- function(x, ...) UseMethod("sampled_pars")
# We need to export/register it? No just defined in environment.

# Mock subject design
subject_design <- list(par_names = c("v", "a"))
# We need sampled_pars to return a NAMED vector/list where names are par names?
# Look at group_design.R line 50: par_names <- names(sampled_pars(subject_design))
# So sampled_pars must return something with names.
sampled_pars.emc.design <- function(x, ...) {
    out <- rep(1, length(x$par_names))
    names(out) <- x$par_names
    return(out)
}
class(subject_design) <- "emc.design"

# Also need to make sure sampled_pars.emc.group_design is available (from sourced file)
# But group_design.R calls it explicitly as sampled_pars.emc.group_design

# Pure R implementation of calculate_subject_means to avoid C++ dependency
calculate_subject_means <- function(group_designs, tmu) {
    n_subjects <- nrow(group_designs[[1]])
    n_pars <- length(group_designs)
    out <- matrix(NA, nrow = n_pars, ncol = n_subjects)
    rownames(out) <- names(group_designs)

    idx <- 1
    for (i in 1:n_pars) {
        dm <- group_designs[[i]]
        n_eff <- ncol(dm)
        beta <- tmu[idx:(idx + n_eff - 1)]
        out[i, ] <- as.vector(dm %*% beta)
        idx <- idx + n_eff
    }
    return(out)
}

group_designs <- group_design(design, data, subject_design = subject_design)

# Helper for Inverse Wishart
riwish <- function(v, S) {
    # Draw W ~ Wishart(v, solve(S))
    W <- stats::rWishart(1, df = v, Sigma = solve(S))[, , 1]
    return(solve(W))
}

# Build a mock sampler object
sampler <- list(
    samples = list(
        theta_mu = array(0, dim = c(2, 1)), # Start 0
        theta_var = array(diag(1, 2), dim = c(2, 2, 1)),
        a_half = array(1, dim = c(2, 1)),
        last_theta_var_inv = diag(1, 2),
        idx = 1
    ),
    group_designs = group_designs,
    par_names = par_names,
    nuisance = rep(FALSE, 2),
    is_blocked = c(TRUE, TRUE), # blocked
    par_group = c(1, 1), # one block
    prior = list(
        theta_mu_mean = rep(0, 2),
        theta_mu_invar = diag(0.1, 2), # Weak
        v = 2,
        A = rep(0.5, 2)
    )
)

# Storage
res_beta <- matrix(NA, n_iter, 2)
res_s <- matrix(NA, n_iter, 1) # s for v
res_u <- matrix(NA, n_iter, n_groups) # u for v (20 groups)

# Loop
curr_alpha <- alphas # Fixed data
# Initialize u, s in sampler (last sample)
sampler$samples$theta_u <- array(0, dim = c(n_groups, 1))
sampler$samples$theta_s <- array(0.1, dim = c(1, 1))
sampler$samples$a_half_s <- array(1, dim = c(1, 1))

cat("Starting Sampling...\n")
for (i in 1:n_iter) {
    # Call Gibbs Step
    out <- gibbs_step_standard(sampler, curr_alpha)

    # Store
    res_beta[i, ] <- out$tmu
    res_s[i, ] <- out$s
    res_u[i, ] <- out$u

    # Update sampler for next step
    sampler$samples$theta_mu <- as.matrix(out$tmu)
    sampler$samples$theta_var <- array(out$tvar, dim = c(2, 2, 1))
    sampler$samples$last_theta_var_inv <- out$tvinv
    sampler$samples$a_half <- as.matrix(out$a_half)
    sampler$samples$theta_u <- as.matrix(out$u)
    sampler$samples$theta_s <- as.matrix(out$s)
    sampler$samples$a_half_s <- as.matrix(out$a_half_s)
    sampler$samples$idx <- 1
}

# Results
cat("True Beta:", beta_true, "\n")
cat("Est Beta:", colMeans(res_beta[200:n_iter, ]), "\n")
cat("True s (SD):", s_v_true, "\n")
cat("Est s (SD):", mean(sqrt(res_s[200:n_iter, ])), "\n")
cat("Est s^2 (Var):", mean(res_s[200:n_iter, ]), "\n")
cat("True u (first 3):", head(u_v_true, 3), "\n")
cat("Est u (first 3):", head(colMeans(res_u[200:n_iter, ]), 3), "\n")

if (abs(mean(sqrt(res_s[200:n_iter, ])) - s_v_true) < 0.1) {
    cat("SUCCESS: Recovery acceptable.\n")
} else {
    cat("FAILURE: Recovery poor.\n")
}
