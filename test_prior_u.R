source("R/group_design.R")
source("R/variant_standard.R")
source("R/get_type_objects.R")
source("R/sampling.R")
source("R/priors.R")
source("R/design.R")
source("R/new_map.R")
source("R/objects.R")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
if (!requireNamespace("mvtnorm", quietly = TRUE)) install.packages("mvtnorm")
library(MASS)
library(mvtnorm)

# 1. Setup Design with Random Effects
n_subjects <- 10
n_groups <- 2
groups <- factor(rep(paste0("G", 1:n_groups), each = n_subjects / n_groups))
data <- data.frame(subjects = 1:n_subjects, group = groups)

# Mock design formula
design_formula <- list(
    v ~ 1 + (1 | group),
    a ~ 1
)

# Realistic subject design mock
subject_design <- list(
    par_names = c("v", "a"),
    Ffactors = list(subjects = factor(1:n_subjects)),
    model = function() list(type = "standard")
)
sampled_pars.emc.design <- function(x, ...) {
    out <- rep(1, length(x$par_names))
    names(out) <- x$par_names
    return(out)
}
sampled_pars.list <- function(x, ...) {
    if (inherits(x[[1]], "emc.design")) {
        return(sampled_pars(x[[1]], ...))
    }
    return(sampled_pars(x[[1]], ...))
}
class(subject_design) <- "emc.design"

# Group Design
group_des <- group_design(design_formula, data, subject_design = subject_design)

cat("Structure of group_des:\n")
print(str(group_des))
cat("Names: ", names(group_des), "\n")

# 2. Test get_prior_standard for 'u'
# We call get_prior_standard directly to bypass plot/get_pars wrapping issues.
cat("Testing get_prior_standard(selection='u')...\n")
prior_obj <- prior(
    design = list(subject_design), type = "standard", group_design = group_des,
    mu_mean = c(0, 0), mu_sd = c(1, 1)
)

samples_u <- get_prior_standard(prior = prior_obj, selection = "u", group_design = group_des, N = 10, design = list(subject_design))

if (!is.null(samples_u$theta_u)) {
    cat("Success: theta_u is present.\n")
    print(dim(samples_u$theta_u))
} else {
    cat("FAIL: theta_u is missing.\n")
}

# 3. Test get_prior_standard for 's'
cat("Testing get_prior_standard(selection='s')...\n")
samples_s <- get_prior_standard(prior = prior_obj, selection = "s", group_design = group_des, N = 10, design = list(subject_design))

if (!is.null(samples_s$theta_s)) {
    cat("Success: theta_s is present.\n")
    print(dim(samples_s$theta_s))
} else {
    cat("FAIL: theta_s is missing.\n")
}

if (!is.null(samples_s$a_half_s)) {
    cat("Success: a_half_s is present.\n")
    print(dim(samples_s$a_half_s))
} else {
    cat("FAIL: a_half_s is missing.\n")
}
