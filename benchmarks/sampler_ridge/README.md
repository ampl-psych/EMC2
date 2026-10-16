# Correlated ridge comparison

The target uses EMC2's compiled Gaussian regression likelihood (`MRI`) with
23 participants, 24 observations each and 12 hierarchical coefficients
(11 regression coefficients and log noise SD). The population covariance is
estimated by the ordinary `standard` sampler.

For each of three seeds, a random orthogonal rotation creates nine sharply
identified directions and two unidentified directions. Predictor strength is
100. All coefficients have the same population-mean prior, with variance 1.
Four chains start at different locations along the unidentified directions.
Both samplers receive exactly the same saved initial states. No parameter
groups, directions or proposal covariances are supplied to the sampler.

Writing the two unidentified directions as the columns of U, XU = 0.
Therefore the marginal posterior of U' mu is exactly N(0, I). This provides
a reference distribution without fitting a second model.

Each run uses 50 preburn, 100 burn, 600 adapt, 100 kernel-tuning and 2,000
retained iterations per chain. The fixed warmup budget replaces the adaptive
stopping rule for this comparison. Particle counts and all other tuning use
the package defaults. Four chains run on four cores. Timings include warmup
and proposal construction, excluding the shared initial-state construction.

Three installed packages are compared in separate R sessions

- `current` is commit `ea74eb2defd12a585793b2a118f4ed9231d9aab3`.
- `curvature` is that commit with only the four-corner mixed derivative in
  `src/sampler/likelihood.cpp`.
- `directions` includes both that derivative and the directional joint moves.

The new joint move replaces coordinate-wise translations with translations
along estimated curvature directions. For a direction d and a symmetric
normal increment z, it proposes mu' = mu + zd and alpha_s' = alpha_s + zd
for every participant. The participant deviations alpha_s - mu do not change.
The acceptance ratio contains the exact likelihood change and the population
mean prior change. Directions and step sizes are learned during adaptation,
cached, and frozen during sampling. The particle and group updates are unchanged.

The original mixed derivative used two diagonal corners and the axial
evaluations. With finite steps and nonlinear parameters, this can introduce
spurious curvature along a flat direction. The four-corner central difference
avoids that error in this target. It changes proposal construction only.

## Run

From the repository root, load an installed package and source `benchmark.R`.
Create the shared inputs once

```r
source("benchmarks/sampler_ridge/benchmark.R")
dir.create("benchmarks/sampler_ridge/results", showWarnings = FALSE)
problems <- lapply(41:43, function(seed)
  ridge_model(seed, p = 11, subjects = 23, weak = 2))
saveRDS(problems, "benchmarks/sampler_ridge/results/problems.rds")
```

For each installed version, start a fresh R session with its library first
in `.libPaths()`, then run

```r
source("benchmarks/sampler_ridge/benchmark.R")
problems <- readRDS("benchmarks/sampler_ridge/results/problems.rds")
label <- "current" # or curvature / directions
out <- lapply(problems, run_ridge, label = label, draws = 2000)
saveRDS(out, paste0("benchmarks/sampler_ridge/results/", label, ".rds"))
```

Each result includes the full fit, per-parameter Rhat and bulk ESS, and the
projections onto the known reference directions. Generated results are ignored
by Git.

An additional check uses `ridge_model(44, p = 11, subjects = 23, weak = 2,
weak_strength = .1)`. Save that input as `weak-problem.rds`, run both samplers
with `draws = 2000`, and save their results as `weak-current.rds` and
`weak-directions.rds`. These directions are weakly identified; their posterior
is no longer the standard normal reference above.

Run `Rscript benchmarks/sampler_ridge/summarise.R` to produce the table and traces.

## Results

Seconds are medians across the three fits. ESS is the median bulk ESS across
the six weak-direction projections. Maximum Rhat covers all population means,
variances and participant coefficients across all three fits.

| Sampler | Seconds | Maximum Rhat | Weak-direction ESS | Fits with Rhat < 1.05 |
|---|---:|---:|---:|---:|
| Current | 9.38 | 4.286 | 4.4 | 0/3 |
| Curvature correction only | 6.83 | 4.158 | 4.4 | 0/3 |
| Automatic directions | 7.01 | 1.020 | 1,839 | 3/3 |

For the exact standard normal reference, the new sampler's six posterior
means range from -0.019 to 0.002 and SDs from 0.970 to 0.989.
The additional weakly identified target gives maximum Rhat 3.762 with the
current sampler and 1.018 with automatic directions.

The existing sample-kernel test passes, including the added check that
directions stay fixed. Continuing a saved fit also preserves the directions,
step sizes and particle kernel exactly.

## DDM5 timing and implementation check

The translation now caches P d and d' P d for every direction d, where P is
the population-mean prior precision. Its prior ratio needs one dot product.
It updates the group state only on acceptance and reuses the shift buffer.
Scale moves write back only the changed coefficient. The directions remain
cached across sampling batches.

The optimized implementation produced identical saved draws to the first
directional implementation on the first ridge dataset, a standard model with
a correlated mean prior, SEM and an infinite-factor model. The existing
sample-kernel continuation test also passes.

`ddm5.R` compares commit `ea74eb2d` against the optimized directional sampler.
It uses the saved DDM5 dataset from `ddm5_delta90_4x2/setup.rds`, with 30
participants, 3,000 trials, diagonal population covariance, unchanged priors
and bounds, and `rt_resolution = 1e-12`. Both runs use the same saved initial
states, four chains with two cores each, 600 warmup iterations (including 100
kernel-tuning iterations), and 1,000 retained iterations. The warmup length is
fixed for the timing comparison. Both packages are installed builds on macOS.

| Metric | Before | After |
|---|---:|---:|
| Warmup seconds | 5.010 | 5.102 |
| Sampling seconds | 3.100 | 3.313 |
| Total seconds, including setup | 8.172 | 8.479 |
| Minimum bulk ESS/s, means | 21.91 | 21.65 |
| Minimum bulk ESS/s, variances | 20.39 | 19.95 |
| Minimum bulk ESS/s, participants | 26.76 | 9.72 |
| Maximum Rhat | 1.02250 | 1.05005 |

The new sampler took 3.8% longer in this repetition. Participant 23's `sv`
had the lowest ESS after the change. These fits were not extended.

Create the common input once with the baseline package loaded

```r
source("benchmarks/sampler_ridge/ddm5.R")
saveRDS(prepare_ddm5(), "benchmarks/sampler_ridge/results/ddm5_input.rds")
```

In separate R sessions using each installed version, set `label` to `before`
or `after` and run

```r
source("benchmarks/sampler_ridge/ddm5.R")
input <- readRDS("benchmarks/sampler_ridge/results/ddm5_input.rds")
result <- run_ddm5(input, label)
saveRDS(result, paste0("benchmarks/sampler_ridge/results/ddm5_", label, ".rds"))
```

## Bisset prior comparison

`bisset_priors.R` resumes the same three saved Bisset chains under the default
prior, with the population means of log `S1.sigma` and log `S1.tau` centred at
`log(.05)` instead of zero, and with those centres plus `v = 8` instead of 2.
The population-mean prior SD remains 1 and `A` remains .3. All other priors,
data, likelihood inputs and the full population covariance model are unchanged.

Each arm adds 100 adaptation iterations, 100 discarded kernel-tuning iterations
and 200 sample iterations. The short comparison uses 25 particles in adaptation
and 25–30 during kernel tuning and sampling, with three chains and eight cores
per chain. Diagnostics use only the new sample iterations. The input snapshot
already contains 2,000 burn and 800 adaptation draws per chain.

```r
source("benchmarks/sampler_ridge/bisset_priors.R")
compare_bisset_priors("benchmarks/sampler_ridge/results/bisset_priors/input.rds",
                      "benchmarks/sampler_ridge/results/bisset_priors")
```

| Prior | Mean S1.sigma variance | Mean S1.tau variance | Maximum Rhat, means | Maximum Rhat, variances | Maximum Rhat, participants |
|---|---:|---:|---:|---:|---:|
| Default | 97.46 | 176.97 | 2.111 | 1.815 | 2.392 |
| Centres at log(.05) | 37.06 | 24.36 | 1.918 | 2.120 | 2.447 |
| Centres at log(.05), v = 8 | 8.05 | 4.39 | 1.691 | 2.084 | 2.406 |

Variances are for the sampled log coefficients. Their means here are averages
over the 200 new sample iterations, not estimates from converged fits. The
average population correlation between these coefficients fell from .989
under the default prior to .926 after recentering and .822 with `v = 8`.
None of the short continuations met Rhat < 1.05. The changed priors reduced the
extreme variance scales but did not rescue convergence within this budget.
The original R26 fit was not modified.

## Joint covariance shears

The shear adds a multiple of coefficient k to coefficient j, for every
participant, and transforms the covariance consistently. A partner is drawn
uniformly from the same covariance block. There are no model-specific pairs.
It runs alongside the existing particle, group, scale and translation moves
in burn, adapt and sample. It applies to the standard covariance blocks.

Write K = Sigma^-1, P for the population-mean prior precision, and m for its
prior mean. Let w_s be the cached likelihood curvature for coefficient j in
participant s. Without group regressors, the pivot offset is

```
o = sum_s w_s (mu_k - alpha_ks) / (sum_s w_s + P_jj)
delta ~ N(0, h_j^2 / (2 v K_jj / a_k + o^2 P_jj))
alpha'_js = alpha_js + delta (alpha_ks - mu_k + o)
mu'_j = mu_j + delta o
Sigma' = M Sigma M',     M = I + delta e_j e_k'
```

The pivot is weighted toward participants with sharply measured coefficients. With
group regressors, o is zero, deviations use each participant's group mean,
and those means stay fixed. The population auxiliary variables a also stay
fixed. The multiplier h adapts during warmup and is frozen during sampling.

The transformation has determinant one. The participant Gaussian densities
and Jacobians cancel. The reverse proposal has the same variance because
K_jj, a_k, alpha_ks, mu_k and o do not change. The full log acceptance ratio is

```
sum_s [log L_s(alpha'_s) - log L_s(alpha_s)]
  + (v / a_k) (2 delta K_jk - delta^2 K_jj)
  - delta o [P (mu - m)]_j - .5 delta^2 o^2 P_jj
```

The final two terms are absent with group regressors. The likelihood is
evaluated exactly. The approximation only determines the pivot. Covariance
and precision updates use row and column operations; no factorization is
needed for the shear. Existing likelihood caches and proposal buffers are
reused, and curvature weights and their denominators are cached.

`tails.R` supplies a second cheap target. Two exponential regression effects
have flat negative log-scale tails. Nine ordinary regression effects and
the noise parameter bring the total to 12 coefficients. There are 23
participants and a fully estimated covariance. Four chains start at different
variance scales, with the two tail coefficients tightly coupled to an observed
coefficient. The sampler receives no information about these roles.

Both installed samplers use the same saved initial states, 50 preburn, 100
burn, 400 adapt, 100 kernel-tuning and 4,000 sample iterations. The latter
were collected in two consecutive batches of 2,000, with the sample kernel
fixed. Baseline is commit `ee05e309`. Rhat and minimum bulk ESS below cover
all population means, variances, correlations and participant coefficients.

| Seed | Baseline max Rhat | Shear max Rhat | Baseline min ESS | Shear min ESS |
|---|---:|---:|---:|---:|
| 61 | 1.211 | 1.017 | 14.9 | 298.5 |
| 62 | 1.055 | 1.020 | 56.8 | 295.1 |
| 63 | 1.272 | 1.033 | 11.4 | 104.6 |

Median total time was 24.03 seconds before and 26.45 seconds after. Median
minimum ESS/s rose from 0.62 to 11.16 when including correlations, and from
3.55 to 11.16 for means, variances and participant coefficients alone. All
three new fits meet Rhat < 1.05, including correlations; no baseline fit
passes that complete check. The first 2,000-draw results are retained separately.

Create inputs once with `tail_model(seed)` for seeds 61–63, then call
`run_tails(input, label)` in a separate R session for each installed version.
Continue each result with

```r
more <- run_tails(list(emc = result$fit, seed = seed), label, stages = "sample")
more$seconds <- more$seconds + result$seconds
```

`check_shear.R` compares the acceptance ratio against full normal and
inverse-Wishart densities, checks reverse-proposal invariance, and samples
a constant-likelihood target with known Huang–Wand prior marginals. The
largest density-ratio discrepancy was 9.24e-13. The prior quantiles agree
within .02 on the checked scales. The sample-kernel continuation test also
passes, including the frozen shear step sizes.

The Gaussian ridge reference also remains correct with shears. Its two
standard-normal projections have means -0.008 and 0.018, SDs 0.990 and 0.990,
and maximum Rhat 1.003. Covariance-block and group-regression checks pass.
The DDM5 check produces identical draws before and after, taking 8.38 and
8.53 seconds respectively. Its diagonal covariance receives no shear moves.

The separate R26 Bisset check uses the same saved states, default priors,
100 adaptation, 100 kernel-tuning and 200 sample iterations as the earlier
short prior comparison, with three chains and eight cores per chain. Its
results are separate from the original fit and the installed package.

| Maximum Rhat | Baseline | Shear |
|---|---:|---:|
| Population means | 2.111 | 1.581 |
| Population variances | 1.815 | 1.734 |
| Participant coefficients | 2.392 | 2.708 |

The Bisset continuation does not meet Rhat < 1.05. Its worst remaining
disagreement is participant 20's d2.mu and d2.mu_smatchd. This short check
does not establish a solution to the real-data convergence problem.

### Legacy Bisset sampler

`bisset_legacy.R` tests the `AccumulatR` branch at `7d31b163` on R26.
It uses `Bisset2014/analysis.R` through `make_emc()`, with the same data,
model, default priors and `rt_resolution = 1e-3`. Three chains start fresh
with seed 415 and L'Ecuyer-CMRG, using 3 × 23 cores and the default particle
factor 50. The fixed budget is 50 preburn, 2,000 burn, 800 adaptation and
1,000 sample iterations per chain. Draws are retained without automatic
convergence trimming. Diagnostics use rank-normalized split Rhat and bulk
ESS over the last 500 draws of each stage.

Install a separate copy of that branch after applying `bisset_legacy.patch`.
The patch only connects its likelihood bridge to the current AccumulatR C++
API; the sampler is unchanged. The likelihood agrees exactly with the
current bridge at 138 saved participant coefficient vectors, and the
population prior parameters agree exactly.

With the current EMC2 installation, create the likelihood reference using

```r
source("bisset_legacy.R")
bisset_likelihood_reference("prior_check/adapt-current.rds", "legacy_check")
```

Then use a separate R session with the legacy installation first in
`.libPaths()` and run

```r
source("bisset_legacy.R")
run_bisset_legacy(normalizePath("analysis.R"), "legacy_check")
```

The original R26 fit and normal package installation are unchanged.

This run was stopped at the user's request during burn, after 700 saved burn
iterations. It did not reach adaptation or sampling. Over its last 500 burn
draws, S1.sigma variance chain means were 48.9, 12.1 and 16.0 (Rhat 1.491);
S1.tau variance means were 25.5, 2.1 and 8.0 (Rhat 1.547).
