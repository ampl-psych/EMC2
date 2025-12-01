# Create Group-Level Design Matrices

Creates design matrices for group-level parameters based on
subject-level design and formulas. This function is used for
hierarchical modeling to specify how subject-level parameters vary
across groups or conditions.

## Usage

``` r
group_design(formula, data, subject_design, contrasts = NULL)
```

## Arguments

- formula:

  A list of formulas specifying the relationship between subject-level
  parameters and group-level predictors. Each formula should have a
  subject-level parameter on the left-hand side and group-level
  predictors on the right-hand side.

- data:

  The same data as used in the subject-level design. Must include a
  'subjects' column.

- subject_design:

  An emc.design object containing the subject-level design.

- contrasts:

  Optional list of contrast matrices to be used for categorical
  predictors.

## Value

A list of design matrices, one for each parameter specified in the
formula. The intercept is automatically included as the group-level mean
and is omitted from the design matrices.

## Details

Here it is important to consider the interpretation of the group-level
mean. This allows one to add covariates/group-level factors to the
model. However, mu, the group-level mean, is still included for all
parameters. Mu represents the intercept in the design matrix, this
intercept is always added to the group-level model. Therefore, to keep
the interpretation of mu as the group-level mean, it is important to
ensure that the design matrix has a mean of zero. If not, this function
will throw a warning. For some unbalanced designs, this is unavoidable
and the warning can be ignored.

## Examples

``` r
# Create subject-level design 
subj_design <- design(data = forstmann, model = DDM,
                      formula = list(v ~ S, a ~ E, t0 ~ 1),
                      contrasts = list(S = contr.helmert))
#> Parameter(s) SZ, Z, s, st0, sv not specified in formula and assumed constant.
#> 
#>  Sampled Parameters: 
#> [1] "v"           "v_S1"        "a"           "a_Eneutral"  "a_Eaccuracy"
#> [6] "t0"         
#> 
#>  Design Matrices: 
#> $v
#>      S v v_S1
#>   left 1   -1
#>  right 1    1
#> 
#> $a
#>         E a a_Eneutral a_Eaccuracy
#>     speed 1          0           0
#>   neutral 1          1           0
#>  accuracy 1          0           1
#> 
#> $t0
#>  t0
#>   1
#> 
#> $SZ
#>  SZ
#>   1
#> 
#> $Z
#>  Z
#>  1
#> 
#> $s
#>  s
#>  1
#> 
#> $st0
#>  st0
#>    1
#> 
#> $sv
#>  sv
#>   1
#> 
# Add some age covariate and roughly demeans
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
# Then you can make the emc object with
emc <- make_emc(forstmann, subj_design, compress = FALSE, group_design = group_des)
#> Processing data set 1
```
