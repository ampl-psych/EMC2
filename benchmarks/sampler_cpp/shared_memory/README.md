# DDM5 shared-memory benchmark

R26, installed EMC2, four chains and 30 participants. Each run continues the same saved DDM5 sample-stage state for 1,000 draws per chain, in ten batches of 100. The proposals, particle counts, likelihood, priors and initial chain RNG states are identical across methods. No Stan or DDM6 fits were launched.

Times are medians of three runs, including worker/context setup and sampling. Common R startup, package loading in the parent and reading the input fit are excluded. Peak RAM is proportional set size (PSS), summed over the parent and every worker in a separate profiling run. Timings are measured without memory profiling.

| Method | 4 × 2 cores, seconds | Peak MiB | 4 × 15 cores, seconds | Peak MiB |
|---|---:|---:|---:|---:|
| Forked workers | 10.88 | 1347 | 24.00 | 1824 |
| Persistent workers, full objects | 12.02 | 1660 | 25.66 | 1661 |
| Persistent workers, samples only | 9.29 | 1345 | 23.81 | 1345 |
| Shared C++ threads | 7.06 | 273 | 3.61 | 275 |

At 4 × 15 cores, the shared engine takes 3.60 seconds for sampling itself, versus 21.58 seconds for persistent workers returning samples only. Including setup, the speedup is 6.6× and peak RAM falls by 80%. At 4 × 2 cores, sampling speed is similar; the shared engine saves worker startup and RAM.

The shared engine runs chains and participant likelihood tasks in one OpenMP pool. Fixed likelihood inputs are shared and cached for the fit; each chain owns its RNG, proposal buffers and likelihood scratch. R object creation, progress and interrupt handling stay on the main thread. Likelihoods that need the R API run on that thread.

The particle, group and joint-update equations are unchanged. The sampler now uses a native RNG per chain, seeded from R. RNG state is stored with the samples and the compiled cache is rebuilt when resuming a saved fit. If the C++ library differs, each chain starts a fresh RNG stream from R's current seed. Fixed seeds therefore produce different draws from the earlier R RNG.

All 32 benchmark outputs contain identical draws and final RNG states. Tests cover serial/threaded equality, frozen proposals, saved resumes, restarting the RNG for a different C++ library, all supported sampler types, single-chain initialization and worker errors. DDM5 resumes were checked on macOS and Linux, and after starting OpenMP in the parent. The normal, gamma and multinomial generators were checked against their analytic moments using one million draws.

## Reproduce

The scripts use `/home/nsteven/2026/StanCompare/emc2_shared` on R26, with the package installed in its `lib` directory. The original fit is read from `results/hier_matched/rep01/ddm5_emc2.rds`. Run `prepare.R`, then `python3 measure.py`, then `summarise.R` to benchmark the installed shared-memory sampler. Each run starts in a fresh R process. The tables above retain the measured comparisons.

Raw timings and the table are saved locally in `../results/shared_memory/measurements.csv` and `../results/shared_memory/summary.csv`. The original R26 runs and installed packages remain unchanged.

## macOS, 4 × 2 cores

Same DDM5 state, 1,000 draws per chain, three timing runs and a separate RAM run. Both methods use four chains and two cores within each chain. RAM is peak combined RSS across the parent and workers.

| Method | Sampling seconds | Including setup | Peak MiB |
|---|---:|---:|---:|
| Persistent workers, samples only | 2.67 | 3.97 | 1585 |
| Shared C++ threads | 2.57 | 2.58 | 358 |

Sampling speed is similar. Shared threads save worker startup and use 77% less combined RSS. All eight runs produced identical samples and final RNG states.

Results are in `../results/shared_memory_macos/results/`. The scripts accept `EMC2_BENCHMARK_ROOT`, `EMC2_BENCHMARK_LIB` and `EMC2_BENCHMARK_FIT` environment variables. This run used that results folder as the root, `/tmp/emc2-shared-lib` and `/tmp/emc2-shared-ddm5.rds`. On macOS, `measure.py` uses 4 × 2 cores.
