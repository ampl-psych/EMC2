# Define Structural Equation Model (SEM) Matrices

This function helps create the specification matrices (Lambda, B, K, G)
for an SEM. It takes a design object, data, factor names, covariate
column names, and list-based specifications for the paths to be
estimated. The subject-level parameter names for Lambda_mat and K_mat
rows are derived from `sampled_pars(design)`. It validates that
covariates are consistent per subject (subject column in `data` must be
named "subjects") and includes an aggregated subject-level covariate
data frame named `covariates` in the output list. For identifiability,
the first parameter listed in `lambda_specs` for each factor is fixed to
1.

## Usage

``` r
make_sem_structure(
  data = NULL,
  design,
  covariate_cols = NULL,
  lambda_specs = NULL,
  b_specs = NULL,
  k_specs = NULL,
  g_specs = NULL,
  fixed_value = 0
)
```

## Arguments

- data:

  A data frame containing a column named "subjects" and any covariate
  columns specified in `covariate_cols`.

- design:

  An emc.design object, as created by the
  [`design()`](https://ampl-psych.github.io/EMC2/reference/design.md)
  function. The parameter names for the SEM are derived from
  `names(sampled_pars(design))`.

- covariate_cols:

  Character vector or NULL. Column names in `data` to be used as
  covariates for K_mat and G_mat. If NULL, no covariates are processed.

- lambda_specs:

  A list defining factor loadings. The list names should be factor names
  and each element should be a character vector of parameter names (from
  `names(sampled_pars(design))`) that load onto that factor. The first
  parameter listed for each factor will be fixed to 1 for
  identifiability. Example:
  `list(Factor1 = c("v_Sleft", "a_Eneutral"), Factor2 = "t0")` Here,
  `Lambda_mat["v_Sleft", "Factor1"]` would be 1.

- b_specs:

  A list defining regressions among factors. List names are outcome
  factors, elements are character vectors of predictor factors. Example:
  `list(Factor2 = "Factor1", Factor3 = c("Factor1", "Factor2"))`

- k_specs:

  A list defining covariate effects on subject-level parameters. List
  names are parameter names (from `names(sampled_pars(design))`),
  elements are character vectors of covariate names (must be present in
  `covariate_cols` and thus in the processed `covariates` data frame).
  Example: `list(v_Sleft = "cov1", a_Eneutral = c("cov1", "cov2"))`

- g_specs:

  A list defining covariate effects on factors. List names are factor
  names, elements are character vectors of covariate names. Example:
  `list(Factor1 = "cov1", Factor2 = c("cov1", "cov2"))`

- fixed_value:

  Numeric. The value used for fixed paths in the matrices that are not
  set to 1 for identifiability or `Inf` for estimation. Default is 0.

## Value

A list containing:

- `Lambda_mat`: The factor loading matrix.

- `B_mat`: The matrix of regressions among factors.

- `K_mat`: The matrix of covariate effects on subject-level parameters.

- `G_mat`: The matrix of covariate effects on factors.

- `par_names`: The subject-level parameter names derived from
  `sampled_pars(design)`.

- `factor_names`: The provided SEM factor names.

- `covariates`: A data frame with one row per unique subject and columns
  for each covariate, containing the unique subject-level values. Column
  names are the covariate names.

## Examples

``` r
# Create a design object (simplified from design.R example)
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))
matchfun_example <- function(d) d$S==d$lR # Example match function

example_design_obj <- design(
  data = forstmann,
  model= LBA,
  matchfun=matchfun_example,
  formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
  contrasts=list(v=list(lM=ADmat)),
  constants=c(sv=log(1)),
)
#> 
#>  Sampled Parameters: 
#> [1] "v"           "v_lMdiff"    "sv_lMTRUE"   "B"           "B_Eneutral" 
#> [6] "B_Eaccuracy" "B_lRright"   "A"           "t0"         
#> 
#>  Design Matrices: 
#> $v
#>     lM v v_lMdiff
#>   TRUE 1      0.5
#>  FALSE 1     -0.5
#> 
#> $sv
#>     lM sv sv_lMTRUE
#>   TRUE  1         1
#>  FALSE  1         0
#> 
#> $B
#>         E    lR B B_Eneutral B_Eaccuracy B_lRright
#>     speed  left 1          0           0         0
#>     speed right 1          0           0         1
#>   neutral  left 1          1           0         0
#>   neutral right 1          1           0         1
#>  accuracy  left 1          0           1         0
#>  accuracy right 1          0           1         1
#> 
#> $A
#>  A
#>  1
#> 
#> $t0
#>  t0
#>   1
#> 

# SEM Factor names

# Make a copy of forstmann for example modification
forstmann_mod <- forstmann
set.seed(123) # for reproducibility
subj_trait_values <- stats::setNames(rnorm(length(levels(forstmann_mod$subjects))),
                                    levels(forstmann_mod$subjects))
forstmann_mod$SubjTrait <- subj_trait_values[forstmann_mod$subjects]

my_cov_cols <- c("SubjTrait")

lambda_example_specs <- list(
  Speed = c("v", "v_lMdiff"), # "v" will be fixed to 1
  Caution = c("B", "B_Eneutral", "B_Eaccuracy", "B_lRright", "A") # "B" fixed to 1
)
b_example_specs <- list(Caution = "Speed")
k_example_specs <- list(t0 = "SubjTrait") # "SubjTrait" must be in my_cov_cols
g_example_specs <- list(Speed = "SubjTrait")

sem_settings_definition <- make_sem_structure(
  data = forstmann_mod,
  design = example_design_obj,
  covariate_cols = my_cov_cols,
  lambda_specs = lambda_example_specs,
  b_specs = b_example_specs,
  k_specs = k_example_specs,
  g_specs = g_example_specs
)

print(sem_settings_definition$Lambda_mat)
#>             Speed Caution
#> v               1       0
#> v_lMdiff      Inf       0
#> sv_lMTRUE       0       0
#> B               0       1
#> B_Eneutral      0     Inf
#> B_Eaccuracy     0     Inf
#> B_lRright       0     Inf
#> A               0     Inf
#> t0              0       0
print(sem_settings_definition$B_mat)
#>         Speed Caution
#> Speed       0       0
#> Caution   Inf       0
print(sem_settings_definition$K_mat)
#>             SubjTrait
#> v                   0
#> v_lMdiff            0
#> sv_lMTRUE           0
#> B                   0
#> B_Eneutral          0
#> B_Eaccuracy         0
#> B_lRright           0
#> A                   0
#> t0                Inf
print(sem_settings_definition$G_mat)
#>         SubjTrait
#> Speed         Inf
#> Caution         0
print(head(sem_settings_definition$covariates))
#>     SubjTrait
#> 1 -0.56047565
#> 2 -0.23017749
#> 3  1.55870831
#> 4  0.07050839
#> 5  0.12928774
#> 6  1.71506499
```
