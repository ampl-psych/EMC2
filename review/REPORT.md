# Review of joint_move at 8b257e3d (commits 88668e92, 6bea6e85, 8b257e3d)

Branch `review/joint_move-8b257e3d`; scripts and logs for every claim are in `review/scripts/` and `review/logs/`. HANDOVER.md covers parts 2-3.

Reviewed on Linux (Ubuntu 24.04, R 4.3.3, gcc 13.3, GNU libgomp). EMC2 was built clean from 8b257e3d and, for side-by-side checks, from PR #202's head 3f24a807. Every behavioural statement below was run against those two builds unless marked "by reading".

Not available to the review: the S3 reference material (review_part2.md, from_conflicteams.md, REPORT_jm.md, aws/README.md), the CI job logs and artifacts (only each job's conclusion was readable), and the EMC2-SBC-archive repository. Finding numbers 1-14 therefore come from the task description and PR #202's commit messages, not from review_part2.md itself.

## (a) Verdict on the new move: correct, and exact as a Metropolis-Hastings step

**What changed.** Before 88668e92 the location move shifted one parameter j of every subject and the group mean by the same scalar delta, through `GroupSampler::affine(state, j, 0, delta)`. Now `JointMoves` (src/sampler/runner.cpp) builds, once per stage, a p x p matrix of directions D and proposes, for each column j, a shift `s = delta * D[, j]` with `delta ~ N(0, loc[j]^2)` applied to every subject's whole parameter vector and to the group mean: `alpha[, s] += s` for all s, `tmu += s`, `subj_mu[, s] += s`. The scale move is unchanged and still goes through `affine`. The moves now also run in `burn` (R/variant_standard.R `scale_move_stages`, R/fitting.R `add_proposals` builds `lik_prec` in burn too), not only in adapt and sample.

**Ratio the code implements.** For the translation the group term N(alpha_s | subj_mu_s, Sigma) is invariant (alpha_s and subj_mu_s move together) and a translation has unit Jacobian, so only the likelihoods and the prior on the group mean change. With prior N(mean, P^-1):

    log p(tmu + s) - log p(tmu) = -delta * <P D[,j], tmu - mean> - 0.5 * delta^2 * D[,j]' P D[,j]

The code's `prior_directions = P %*% D`, `curvature[j] = D[,j]' P D[,j]` and `lp = -delta * dot(prior_directions.col(j), tmu - mean) - .5 * delta^2 * curvature[j]` are exactly this; the likelihood term is `sum_s (ll_new[s] - ll[s])` from a full re-evaluation of every subject. The proposal is symmetric (delta and -delta have the same density; D and loc are fixed within a stage and frozen in sample), so the acceptance `log U < lp + sum(ll_new - ll)` is the exact MH acceptance. Checked numerically against `mvtnorm::dmvnorm` differences for 200 random (P, D, tmu, delta): worst relative error 3e-14; the group-term change under a joint shift is 0 to machine precision.

**Prior used.** P is `theta_mu_invar` when present, otherwise `inv_sympd(theta_mu_var)`; `mean` is `theta_mu_mean` (recycled). This matches what GroupSampler (standard/diagonal/blocked/diagonal-gamma), SEMSampler (`1 / theta_mu_var`) and InfiniteFactorSampler (`1 / theta_mu_var`) use for their own group-mean updates, so the move is consistent with the Gibbs steps for every hierarchical type. For SEM and infnt_factor the subject means are `tmu + ...` so the joint shift also leaves their group terms unchanged (by reading).

**Which parameters move.** All p sampled parameters of all n subjects plus the p group means, p times per sweep (one move per direction), after the p scale moves of the same sweep. State is updated sequentially so each move sees the previous one's accepted state. Types with a group-level design keep the location move off (`location_ = !present("group_designs")`), verified: with `group_design()` the attribute has adapted `step` but `loc` untouched and no `directions`. Verified on: standard, diagonal, diagonal-gamma, blocked (directions built, loc adapted in all four).

**Directions.** D = L V diag(1/sqrt(max(lambda, 1))) where L L' = theta_mu_var, and V, lambda are the eigenpairs of I + L' (sum_s H_s) L with H_s the finite-difference likelihood precision of subject s from `lik_prec`. In the prior-whitened coordinates this is the curvature of the joint-translation target, so a unit step along D[, j] changes the log target by about 0.5, which is why `loc` restarts at 1 whenever D is rebuilt. This only shapes the proposal; validity does not depend on it.

**Thread safety and RNG.** `draw()` touches no R object (the Rcpp `settings` list is only read in the constructor, on the R thread, and written in `finish()`, also on the R thread). D, P D, curvature, mean are per-chain and read-only after construction; `shift`, `proposals`, `values` are per-chain buffers; the likelihood evaluation uses the chain's own `Likelihoods` context, as the particle step does. Each chain's `Random` draws one `rnorm` and one `runif` per move on the chain's own thread, as before. Nothing is shared between chains. Cores 1x1, 1x4 and 2x2 give bit-identical draws (LNR and DDM).

**The finite-difference change** in 88668e92 (src/sampler/likelihood.cpp) replaces the two-point mixed partial with the standard four-point formula `-(f++ - f+- - f-+ + f--) / (4 h_j h_k)`; it is second-order accurate where the old one mixed in the axis points' errors, and costs 2p(p-1) instead of p(p-1) extra evaluations per subject per rebuild.

**Verified by sampling.** (i) The move fires: after `burn` the `scale_move` attribute carries a 4x4 `directions` and `loc` adapted away from 1; after the adapt tail `loc`, `step` and `directions` are frozen and identical before and after a sample stage (as `test-sample-kernel.R` now also asserts); draws differ with `options(emc.scale_move = FALSE)` from the same seed, so the moves act in sample. (ii) ON/OFF posterior comparison, 4-subject LNR (forstmann), 3 chains, preburn/burn/adapt as in fit() and 4000 sample iterations, three seed pairs with the moves on (default) and off (`options(emc.scale_move = FALSE)`), plus two same-setting pairs as the null reference. Pooling the three runs per setting, the mean differences on mu (4), sigma2 (4) and alpha (16) parameters have |z| at most 0.85 / 1.66 / 1.51 against the ESS-based Monte Carlo error (none above 2; the single-pair comparisons and the OFF-vs-OFF null pairs both reach |z| 2.6-3.9 on a few alphas, so that is the resolution of the test, not a bias), posterior sd ratios on/off 0.92-1.18, and the moves raise the ESS (pooled median alpha ESS 8486 vs 6262, sigma2 13004 vs 10911). No evidence against exactness at this resolution; the SBC in part 2 is the definitive test.

**Caveats found (none affects exactness):**

1. `arma::chol(covariance, "lower")` and `arma::eig_sym(...)` in the JointMoves constructor are the throwing/unchecked forms. A prior covariance that is not positive definite (a zero prior variance, or a user matrix with a tiny eigenvalue) stops the stage with Armadillo's "chol(): decomposition failed" from inside `native_run_chains`, not an EMC2 message, and a failed `eig_sym` would leave `directions` empty and index out of bounds in `draw()`. Guard both (fall back to D = L, i.e. no curvature, or to the identity) and name the prior in the error.
2. Because `create_lik_prec()` clears `directions` at every `add_proposals()` call, D and `loc` restart (loc = 1) at every step of burn and adapt, including the 100-iteration tail that fixes the sample kernel: the frozen `loc` is the mean of log loc over iterations 21-100 of a window that started at 1. It converged to 2.5-4.5 in the runs here (the directions are small, so delta * D is a small move), so the kernel is usable, but the frozen value is a running-adaptation average, not a converged one. Carrying `directions` across steps within a stage (rebuild only when the stage changes), or starting `loc` from its previous value when D is rebuilt, would give the tail a warm start.
3. Nothing in the new code logs the location move's acceptance; only `loc` is visible. For the SBC report I will read `attr(samples, "scale_move")` (`directions` present, `loc` moved from 1) as the evidence that the move fired.

## (b) Behavioural checks, side by side

| Check | 8b257e3d (joint_move) | 3f24a807 (PR #202) |
|---|---|---|
| set.seed() before the first run of a fresh object changes the draws; same seed reproduces | yes / yes | yes / yes |
| set.seed() before a run on an already initialised object (resume, or a second run) changes the chain's own stream | **no**: draws and the stored RNG state are identical under seeds 1 and 2 (`native_run_chains` alone, and 5 more preburn iterations through `run_emc`). Only stages that rebuild proposals in R differ by seed (burn), through the R-side proposal construction. | yes (Random::mix) |
| cores 1x1 vs 2x2 bit-identical (LNR, 2 subjects) | yes | yes |
| DDM with probit Z: 1x1 == 1x4 == 2x2, and subject threads used | yes; 8 subjects, 20 preburn iterations: 3.5 s / 2.8 s / 2.1 s elapsed, user/elapsed 1.0 -> 1.75 -> 1.9 for LNR | yes |
| Resume at every stage boundary, incl. more adapt after adapt's tail and adapt from a sample object (#202's test-resume.R, NOT_CRAN=true) | 6/6 pass. Niek's 6bea6e85 truncates mix/epsilon/prop_performance in `check_mix()` etc., which is enough. | 6/6 pass |
| Settings/kernel error rather than an Armadillo message when particle settings mismatch the kernel | not present: `ParticleSampler` no longer checks the lengths; with the R-side truncation the mismatch cannot normally arise, so the message is only lost for hand-edited objects | present |
| `on_singular$on_exhausted` from earlier scripts | **error**: "unknown on_singular field(s): on_exhausted" (test-resume.R's first test fails on this expectation, 5/6 pass) | accepted with a warning |
| Threaded fit, then threaded fits in forked children (mclapply) return | yes, 0 s wait, children return idx 11 as expected (#202's test-omp-fork.R second test: 2/2 pass). Niek's src/openmp.h releases libgomp's pool on return (RAII, dladdr check), `configure` probes `omp_pause_resource_all` (OMP_PAUSE 1 here). | yes |
| `omp_runtime()` / `omp_release_pool()` / build-info "OMP runtime" line | not present (test-omp-fork.R's first test errors: object 'omp_runtime' not found) | present |
| CDM warns once that it runs serially when cores are requested | **no warning at all**; it runs serially silently | warns, naming the reason |
| Singular-covariance retry covers SEM/infnt_factor decompositions | yes (6bea6e85: non-throwing chol/solve with "singular ..." messages) | yes |
| R CMD check --as-cran (local, --no-manual --ignore-vignettes, CRAN index stubbed because CRAN is unreachable here) | 0 ERROR; 2 WARNING (no `checkbashisms` script here; locale not available) and 5 NOTE (vignette index because vignettes were skipped; DiagrammeR not installed; installed size 36.5 MB with Ubuntu's -g flags; clock; Ubuntu's non-portable flag) all environmental. 217 tests pass, 87 skipped (NOT_CRAN unset). `roxygen2::roxygenise()` changes nothing substantive (one blank line in run_emc.Rd). devtools::check() not run separately; it adds roxygen to the same check. | CI green on 3 OSes (run 38724643897) |
| Test suite with NOT_CRAN=true on Linux | **58 failing expectations, all snapshot mismatches**, 343 pass. The flat init_chains()/make_data() snapshots (choice, circular_models, group, joint, model_functions, trend) were regenerated on macOS and fail on Linux because std::gamma_distribution differs between libc++ and libstdc++ (the reason #202 made them per-OS); the kept Linux variants (Linux/fit, Linux/fMRI, Linux/variant_funs) are stale against Niek's draws; S3.md differs only by this environment's ASCII locale. | per-OS snapshots from the runners, green |
| CI on joint_move | R-CMD-check failed on ubuntu, windows and macos for every commit (8ad80ee1, 88668e92, 8b257e3d); annotation says only "R CMD check found ERRORs". The Linux run above is the likely content. | green at 3f24a807 |

## (c) PR #202 against joint_move, mechanism by mechanism

`git diff origin/joint_move_pr origin/joint_move` (tests and snapshots aside) touches 36 files. Niek's three commits re-implement most of #202 differently; this is what each fix looks like on each branch and whether the difference loses anything.

| Topic (#202 commit) | #202 | joint_move (Niek) | Loses anything? | Recommendation |
|---|---|---|---|---|
| set.seed() drives resumed runs (091f53cb, `Random::mix`: engine reseeded from stored state + 8 words from R at every run start) | yes | not ported; stored mt19937_64 state alone drives every run after init, set.seed() has no effect on the chain stream | **yes**: the behaviour in the table above. Repeated fits of one object cannot be given different streams except by re-initialising; three production fits of one case need different objects or `init_chains` seeds. (The ConflictEAMs label/attempt seed rule must be checked against this; the scripts are in S3.) | carry over to joint_move (one method and one call in the Chain constructor) |
| RNG format = C++ library family, not version (b5612cc9) | yes | yes (8b257e3d, same code) | no | done |
| Fork safety with libgomp (091f53cb: `omp_release_pool()` exported, called by RAII `PoolRelease` in the four threaded entry points and in `calc_ll_multithreaded`; `omp_runtime()`; configure probe `EMC2_OMP_PAUSE`) | yes | 6bea6e85: src/openmp.h `OpenMPPoolRelease` RAII in the same five places; the libgomp test (dladdr) moved into the destructor; configure probe kept (with a macOS variant linking -lomp). `omp_runtime()`, `omp_release_pool()`, the build-info runtime line and test-omp-fork.R dropped | functionally no (fork test passes). Loses the diagnostic (`emc2_build_info()` no longer says which OpenMP runtime is linked, which the AWS logs would otherwise record) and the test | keep Niek's header; carry over `omp_runtime()` + the build-info line + test-omp-fork.R (its first test adjusted to whatever helper is exported) |
| Thread-safe pnorm pre-transform (091f53cb: `PNORM_STD` in `fill_particle`; `safe = false` only for PNORM_MODE=0) | yes | 8b257e3d goes further: PNORM_MODE 0 (R::pnorm) removed from configure and pnorm_utils.h, `R::pnorm` removed from transform_utils.cpp too, so every build is thread-safe | no (mode 0 was the "exact but serial" option; Hart is 1e-15) | done |
| Serial-likelihood warning naming the reason (`unsafe_reason`, `warned_serial`) | yes | dropped | **yes** (silent serial run for CDM/PSDM/PHSDM, custom trend kernels) | carry over (small) |
| Resume at any stage (3ca07e3b: R-side truncation in `check_mix/epsilon/prop_performance`; length checks in `ParticleSampler`; `on_exhausted` accepted with warning; test-resume.R) | yes | 6bea6e85 ports the R-side truncation only | resume works (6/6). Loses the C++ length check (message) and `on_exhausted` compatibility (error for old scripts) and the test | carry over `on_exhausted` (3 lines) and test-resume.R; the length check is optional |
| Particle selection steps back from a zero-weight particle (091f53cb) | `while (selected > 0 && !(weights[selected] > 0)) --selected;` | 6bea6e85: `threshold >= weights[selected]` and `while (weights[selected] == 0) --selected;` without the `selected > 0` guard | no in practice (the walk back always meets the max-weight particle, weight 1), but the unguarded loop reads before index 0 if it ever does not | add the guard when carrying the rest over |
| `Likelihoods::one` throws for a block whose component no part carries; `shared_likelihood_index()` (091f53cb) | both | `shared_likelihood_index()` ported (6bea6e85); the throw dropped | yes: a mis-specified block samples its prior silently again | carry over (4 lines) |
| Bounds checks (091f53cb: `SampleStore::put` array length, `is_blocked`/`par_group` length, history index range, context chain count) | yes | all dropped | yes (memory safety/diagnostics only; no failure seen) | carry over |
| SEM/infnt_factor singular decompositions retried (091f53cb) | yes | yes (6bea6e85, same code minus comments) | no | done |
| Build: FLIBS (f7a0cac4 dropped it; 139da1d1 restored the LAPACK/BLAS/FLIBS triple except on macOS via `@LINALG_LIBS@`) | triple except macOS | 8b257e3d: the triple everywhere plus Accelerate on macOS, and `tools/check-linalg.sh` runs `R CMD SHLIB` at configure time and stops with install instructions when BLAS/LAPACK/Fortran do not link | no for CRAN ("portable use of BLAS_LIBS/LAPACK_LIBS" passes). macOS users without gfortran now get a clear configure error instead of a working build; CRAN's macOS toolchain has gfortran, so this is the conformant choice | accept Niek's |
| Header dependencies (f7a0cac4: `$(wildcard)`, GNU make declared in 139da1d1) | GNU make | 8b257e3d: `tools/header-dependencies.sh` appends one POSIX dependency line per header to the generated Makevars (and Makevars.win, now generated from Makevars.win.in by configure.win); `SystemRequirements: GNU make` removed | no (same effect, portable make) | accept Niek's |
| icpc-safe `static_cast<SEXP>` constructors (f7a0cac4) | yes | yes (6bea6e85) | no | done |
| Snapshots per OS (b69c7e4c, e9962cd8, 3f24a807) | per-OS variants from the runners, green | flat snapshots regenerated on macOS; Linux variants left stale | **yes**: CI red on all three OSes, 58 failing expectations on Linux | regenerate per OS from the runners once the sampler is final (see must-fix) |
| NEWS (50b5610c, d703d7a4) | describes the C++ sampler as it is | 8b257e3d keeps the older text: mentions `options(emc.sampler = "legacy")` (no such option on joint_move), "L'Ecuyer-CMRG streams", scale moves "for standard hierarchical models" (they run for every hierarchical type and now in burn), no entry for resume-at-any-stage, fork safety, set.seed semantics, the Intel KMP_AFFINITY note, the build changes, or the new curvature-guided translations | yes | carry over #202's NEWS text and add an entry for 88668e92 |
| `.Rbuildignore` `^\.git$`, `\.DS_Store$` | yes | dropped; src/.DS_Store deleted instead | no | drop |
| `.gitattributes` LF for build_info.h.win | yes | extended to configure*, cleanup*, tools/*.sh, src/Makevars* | no (better) | done |

Draw-based group precision, the "stability gate" (8b257e3d, src/sampler/likelihood.cpp): `likelihood_precision()` now uses the draw-based curvature only when the coefficient of variation of every marginal group precision over the window is below .5, and indexes the history through the selected stages consistently. By reading, the gate is what the task calls "R's stability gate".

## (d) Findings 6, 8, 11, 14 (left for Niek)

Without review_part2.md the mapping from number to text is the task's: "proposal precision without R's stability gate; non-PD proposal covariance on raw epsilon_inv; sparse checkUserInterrupt", one of the four cut off.

- Proposal precision without R's stability gate: **touched** by 8b257e3d (the CV < .5 gate above).
- Non-PD proposal covariance on raw epsilon_inv: **touched** by 8b257e3d (src/sampler/proposals.cpp conditions on `log(epsilon)` for every type, not only infnt_factor).
- Sparse checkUserInterrupt: **not touched**; `sample_chains()` still checks once per block of 10 iterations, so a DDM fit at 96 cores answers Ctrl-C after up to 10 iterations.
- The fourth: cannot say without the document.

## (e) Not yet on joint_move, with a recommendation for each

Carry over to joint_move (small, independent of the new move):
1. `Random::mix` and its call in the Chain constructor (set.seed() for resumed runs).
2. `resolve_on_singular()` accepting `on_exhausted` with a warning.
3. The serial-likelihood warning with `unsafe_reason`.
4. `Likelihoods::one` throwing when no part carries the component.
5. The four bounds checks (SampleStore::put, is_blocked/par_group, history index, context chain count) and the `selected > 0` guard in the two selection loops.
6. `omp_runtime()` and the "OMP runtime" line in `emc2_build_info()`.
7. tests/testthat/test-resume.R and test-omp-fork.R.
8. NEWS text from 50b5610c and d703d7a4, plus an entry for the curvature-guided translations and the burn-stage moves.

Drop (superseded by Niek's equivalent): the `@LINALG_LIBS@` configure logic, `$(wildcard)` header rule and `SystemRequirements: GNU make`, `PoolRelease` in sampler.cpp and the dladdr code in particle_ll.cpp, `.Rbuildignore` lines, the `ParticleSampler` length check (optional).

Regenerate rather than carry: the per-OS snapshot variants; #202's files are for #202's draws, and Niek's sampler draws differently (burn moves, selection change), so they would fail too. They have to come from the three CI runners after the last sampler change.

PR #202: close as superseded once items 1-8 are on joint_move (as a new small PR or direct commits); its diff against joint_move is now dominated by snapshots and re-implemented build code, and it is unmergeable ("dirty") as it stands.

## (f) Must fix in the new commits before a merge to dev

1. CI is red on all three OSes at every joint_move commit; the Linux failures are the snapshots. Regenerate per OS from the runners (the #202 recipe) after the last change that alters draws, and restore `variant = Sys.info()[["sysname"]]` on the init_chains()/make_data() snapshots.
2. Guard the JointMoves constructor (chol and eig_sym) so a non-PD prior covariance gives an EMC2 error or a fallback, not an Armadillo message or an out-of-bounds read.
3. NEWS as above (it currently documents an option that does not exist).
4. set.seed() on resumed runs (carry-over 1), because the documented semantics in #202's NEWS, and any seed rule in production scripts, assume it.

Not blocking but worth doing with the above: the warm start of `loc` across steps (caveat 2), and the `on_exhausted` compatibility if any production script still passes it.

## (g) Devtools::check()

`devtools::check(document = TRUE, args = c("--no-manual", "--ignore-vignettes", "--as-cran"))` on the 8b257e3d tree (same Linux box, CRAN index stubbed, duration 6 min): **1 ERROR, 2 WARNINGs, 4 NOTEs**. The ERROR is the test suite (devtools sets NOT_CRAN=true, so the snapshot tests run): `[ FAIL 58 | WARN 0 | SKIP 0 | PASS 343 ]`, every failure a snapshot mismatch as in section (b). The WARNINGs (no `checkbashisms` script here; locale not available) and NOTEs (DiagrammeR not installed; installed size 36.5 MB; clock; Ubuntu's non-portable compiler flag) are environmental. `roxygen2::roxygenise()` (7.3.1 here against 7.3.3 used by the package) rewrites nothing substantive (a blank line in run_emc.Rd and the locale rendering of a name in EMC2-package.Rd).

## (h) Measurements for the plan

- One hierarchical SBC replicate through the package's own worker (fit() defaults, 20 subjects, 100 trials per cell, 3 chains x 1 core, 4-vCPU Linux box): LBA (vignette design and prior) 2.1 min; DDM with probit Z (vignette design and prior, diagonal-gamma) 52 min.

## (i) Follow-up: branch `fix/joint-move-directions` (one commit on 8b257e3d)

Implements caveats 1 and 2 of section (a): `create_lik_prec()` clears the stored directions only when the stage changes (`add_proposals()` passes `stage != get_last_stage(emc)`), so the directions are built at a stage boundary and at the sample-kernel tail and `loc` adapts continuously within a stage; the JointMoves constructor uses the non-throwing `chol()` with an eigenvalue square root as fallback, skips non-finite likelihood precisions and falls back to prior-scaled axes if the whitened curvature cannot be decomposed. A regression test is added to test-sample-kernel.R and NEWS has an entry. Checks on the Linux build of that branch (scripts 11, 12, 04, 08, 09): directions identical across three burn steps (8b257e3d: rebuilt each step), rebuilt at burn->adapt and at the tail, kept in sample; `loc` reached 5.2-5.8 before the tail after 400 continuous uses where 8b257e3d restarted from 1 every 25; a zero prior variance runs (8b257e3d: "chol(): decomposition failed"); PR #202's resume test 6/6 and fork test pass as before; the three-seed on/off posterior comparison pooled over seeds has |z| at most 1.37 (mu), 0.94 (sigma2), 1.54 (alpha), none above 2, sd ratios 0.94-1.20, ESS higher with the move on (alpha 8861 vs 6262). Single seed pairs still reach |z| 2.6-3.8 on a few alphas, as the off-vs-off null pair does; the SBC is the definitive accuracy test of the per-stage rebuild.
