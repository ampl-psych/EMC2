# Test prior plots
source("R/group_design.R")
source("R/variant_standard.R")
source("R/get_type_objects.R")
source("R/sampling.R")
source("R/priors.R")
source("R/design.R")
files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in files) source(f)

if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
library(MASS)

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

# Prior
pri <- prior(
    design = list(subject_design), type = "standard", group_design = group_des,
    mu_mean = c(0, 0), mu_sd = c(1, 1)
)

# Mock C++ function
draw_alpha_from_design <- function(group_design, mu, var) {
    p <- length(group_design)
    if (p == 0) p <- nrow(mu)
    n_subjects <- 10 # Hardcoded to match test setup
    N <- ncol(mu)

    alpha <- array(NA_real_, dim = c(p, n_subjects, N))
    for (i in 1:N) {
        # Just draw random normal
        alpha[, , i] <- rnorm(p * n_subjects)
    }
    return(alpha)
}

cat("Testing plot(pri, selection='alpha')...\n")
tryCatch(
    {
        pdf(NULL)
        plot(pri, selection = "alpha", do_plot = FALSE, N = 10)
        dev.off()
        cat("Success: plot(pri, selection='alpha') ran.\n")
    },
    error = function(e) {
        cat("Error in plot alpha: ", e$message, "\n")
        print(e)
    }
)

cat("Testing plot(pri, selection='u')...\n")
tryCatch(
    {
        pdf(NULL)
        plot(pri, selection = "u", do_plot = FALSE, N = 10)
        dev.off()
        cat("Success: plot(pri, selection='u') ran.\n")
    },
    error = function(e) {
        cat("Error in plot u: ", e$message, "\n")
        print(e)
    }
)

cat("Testing plot(pri, selection='s')...\n")
tryCatch(
    {
        pdf(NULL)
        plot(pri, selection = "s", do_plot = FALSE, N = 10)
        dev.off()
        cat("Success: plot(pri, selection='s') ran.\n")
    },
    error = function(e) {
        cat("Error in plot s: ", e$message, "\n")
        print(e)
    }
)
