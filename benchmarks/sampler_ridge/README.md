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
