# MTLBA golden fixtures (Stage M0, 2026-09-30)

Reference values for the multiple-threshold LBA rating model (`MTLBA`, plan
`rating-clnr-plan.md`, stages M1 onward). Built before any package model code exists, from code
that is independent of the package. Nothing here is run by testthat; the tests written in M1+
read the `.rds` files.

Model convention the fixtures are for (plan section 1): two LBA accumulators, drifts bivariate
normal with correlation rho **truncated to both positive**, start points U(0, A), **all
confidence thresholds at or above A**.

The common object is the **kernel**: the density that the *stopping* accumulator reaches its
stopping threshold `bs` at decision time `dt` while the *other* accumulator's evidence lies in
the band `(dl, du)` (`dl = 0` = the lowest band). Balance of evidence (BE) stops at the choice
threshold b; threshold count (TC) stops at any threshold. A response likelihood is one kernel
(or, for don't know, the sum of two).

## Files

| file | what | generator |
|---|---|---|
| `dk_density.rds` | (a) don't-know paper `n1PDF.DK`, 5 parameter sets x 3 responses x 10 times, rho = 0, d >= A only. **Exact parity target.** | `make_dk_density.R` |
| `bete_density.rds` | (b) BETE `n1PDF_MTR_1`, 4 sets x 5 rho (0, -.7, -.3, .4, .8) x K = 1..4 x BE/TC x all 2K responses x 8 times = 6400 cells, plus BETE's response lookup tables `dL` | `make_bete_density.R` |
| `ref_density.rds` | (c) the same 6400 cells by two independent integrations (`dens_1d`, `dens_2d`) next to `dens_bete`, the normaliser P(both > 0) three ways, and the agreement table | `make_reference_density.R` |
| `sim_summaries.rds` | (d) 1e6-trial summaries (response proportions, mean RT, RT quantiles .1/.3/.5/.7/.9) from `rlba.DK` (5 sets) and BETE `rLBA_BE` / `rLBA_TC` (4 sets x rho 0, -.7, .8 x K = 2..4 x BE/TC, each also in its 2K-1 don't-know form), with expected proportions from the reference density | `make_sim_summaries.R` |
| `bete_guess_density.rds` | (e) **Stage M7.** BETE's likelihood with confidence guessing (`log_likelihood_mt`: trial-level mixture `(1 - gp) * ll + gp * LBAll * pR`), the 4 sets at rho = 0 x K = 1..4 x BE/TC x 3 guess settings x even and odd (don't-know) response sets x 8 times = 6720 cells, with the two components. Parity target of `MTLBA(guess = "trial")` only; the per-accumulator trigger model has no reference code with the both-positive convention | `make_bete_guess_density.R` |
| `mtlba_ref_lib.R` | parameter sets, grids, and the reference integrals `kernel_1d`, `kernel_2d`, `ppos_1d`, `ppos_2d` | (shared) |
| `ref_extract/*.R` | the reference functions, **verbatim** (text copied by srcref, with source line numbers and the md5 of the source file) | `extract_refs.R` |

Each `.rds` is a list with a `provenance` element (source, column meanings, generator, date,
R version).

To rebuild (from this directory, in this order; base R + mvtnorm + msm; at most 8 cores):

    Rscript extract_refs.R          # only if the reference tree changes; needs rating-work/mtlba/ref
    Rscript make_dk_density.R       # seconds
    Rscript make_bete_density.R     # seconds
    Rscript make_reference_density.R  # ~20 s
    Rscript make_sim_summaries.R    # ~6 min on 8 cores

## Provenance of the reference code

- **Don't-know paper** ("A dynamic model of deciding not to choose", https://osf.io/zxp5u):
  DMC code from `DMC_DK_OSF.zip` (https://osf.io/download/akzxw/),
  `rod/dmc/models/LBA-MTR/dists.R` (`dlba.norm`, `plba.norm`, `rlba.DK`, `n1PDF.DK`; this is the
  directory of the paper's final model `lba_DK_qD.R`). The older `MTR-LBA/dists.R` version of
  `n1PDF.DK` (valid only for d >= A) is extracted too and gives the same values on the grid
  (max difference 2e-16). Independent positive-truncated drifts (`posdrift = TRUE`) = the
  both-positive convention at rho = 0.
- **BETE**: EMC2 branch `origin/BETE` (commit 9eab6fd0, Oct 2025), `R/model_MTR.R`
  (`get_dL`, `rLBA_BE`, `rLBA_TC`, `recintab0` -- itself copied from CRAN package MomTrunc --
  and `n1PDF_MTR_1`). Needs `mvtnorm::pmvnorm`.
- **Independent references** (`mtlba_ref_lib.R`, written for this stage from the plan's
  equations, no code shared with either source):
  `kernel_1d` = `integrate()` over the stopping accumulator's drift with the closed-form band
  probability of the other accumulator; `kernel_2d` = brute-force nested `integrate()` of the
  joint bivariate-normal drift density times the (piecewise linear) probability over the other
  accumulator's start point, with no conditional-normal algebra. Both use `rel.tol = 1e-13` on
  pieces split at the integrand's kinks and peaks.

## What was checked (numbers from the generator runs)

- (a) `n1PDF.DK`: both versions agree to 2e-16; the don't-know density does not depend on the
  accumulator order; the three responses integrate to 1 in every set; agreement with `kernel_1d`
  at rho = 0 is 6e-15 absolute.
- (b) vs (c): over all 6400 cells the largest absolute differences are 3.4e-15 (`dens_1d` vs
  `dens_2d`) and 1.5e-14 (`dens_bete` vs either). Where the density exceeds 1e-5 (3818 cells)
  the largest relative differences are 6.4e-11 and 4.8e-10. The normaliser P(both > 0) agrees
  to 2e-16 three ways. So the requirement "independent references at 1e-10, and BETE agrees
  with them" holds in absolute terms everywhere and in relative terms in the body.
- (d) simulated response proportions against the reference density integrated over time:
  437 response cells with an expected count of at least 10 give z scores with sd 1.00 and
  largest |z| 3.6 (2 of 437 beyond 3; that configuration re-simulated with four new seeds gave
  |z| < 1.8, so chance), largest absolute difference in a proportion .0012; the expected
  proportions sum to 1 in every configuration; no NA responses. Holds for `rlba.DK`
  (15 cells, |z| <= 1.7), BE and TC, at rho = 0, -.7 and .8. Ten cells have expected counts
  below 10 (nearly all in set `sharp`, down to 5e-13): their quantiles are NA or very noisy.

## Things a user of these fixtures must know

1. **Far tails are absolute-accuracy only.** `n1PDF.DK`, `n1PDF_MTR_1` and `kernel_1d` all
   form a band probability as a difference of two nearly equal numbers, so where the density is
   below about 1e-8 they are right to ~1e-15 absolutely but not relatively (values of 1e-40
   differ by several percent between routes; one `n1PDF.DK` cell is -1e-29 and `dens_1d` can be
   -2e-15). `dk_density.rds` flags these cells (`tail`). `dens_2d` has a non-negative integrand
   and is the arbiter there. Test parity on the density with an absolute tolerance, or on the
   log density only where `dens > 1e-8`.
2. **Don't know = SUM.** `n1PDF.DK` returns the sum of the two lowest-rating terms (the plan's
   definition). BETE's likelihood averages them (a constant log 2 per don't-know trial).
   `bete_density.rds` stores the separate terms only.
3. **Response scales.** `dk_density.rds` has the EMC2 coding (`R`, `RR`; don't-know rows carry
   the placeholder `R = 1`, `RR = 1`). The BETE files use BETE's single response `R = 1..2K`
   (high-confidence accumulator 1 ... high-confidence accumulator 2) together with the kernel
   indices (`sa`, `tt`, `nsa`, `lt`, `ut`, 1-based into the threshold matrix
   `c(DT0 = 0, DT1, .., b)`), so a test can address the kernel directly. Turning BETE's R into
   EMC2's (R, RR) is M1's response lookup. Note for TC: on low-confidence responses the
   stopping accumulator is the one that is *behind* (e.g. K = 2, R = 2: accumulator 2 crosses
   its first threshold while accumulator 1 is already between its first threshold and b), so
   the choice is the non-stopping accumulator there.
4. **One correction to BETE's simulator.** `rLBA_BE` bins a losing accumulator 1's evidence
   with accumulator 2's b as the top break, which gives NA ratings when b1 > b2 (set `asym`).
   The BE summaries use the verbatim function with that one index corrected (done by text
   substitution in `make_sim_summaries.R`, which also checks the corrected function is
   identical to the original when b1 <= b2). `rLBA_TC`, `rlba.DK` and both density functions
   are used unmodified. BETE's `TC3PLBA`/`TC4PLBA` model lists call `rLBA_BE` as their
   simulator (not used here).
5. **Substituted samplers.** tmvtnorm is not installed on the build machine, so
   `tmvtnorm::rtmvnorm(n = 1, lower = c(0, 0))` is replaced by a rejection sampler from the
   untruncated bivariate normal (tmvtnorm's own default algorithm); `rtnorm` is `msm::rtnorm`.
6. Thresholds in the BETE grids are `d_k = A + B (k/K)^pow`, `b = A + B` -- proportions of B
   above A, the plan's parameterisation F2 (`mtlba_thresholds()`; the matrices are stored in
   `bete_density.rds$thresholds`).

The converted real data (don't-know Experiments 1 and 2, r94 unfiltered with per-participant
1st/99th RT percentiles) are outside the package, in `rating-work/mtlba/stageM0/data/`.

## dmc_p3_density.rds (Stage M7c, 2026-10-07)

The manuscript's DMC model `MTRp3` (the MTLBA's default convention: drift pair conditioned on
not both negative, lowest threshold free below A, per-accumulator confidence guessing) from the
manuscript's own DMC code (`rating-work/mtlba/stageM4/dmc_copy`, model `MTRp3_censor.R`,
`n1PDFMTRp3` / the censoring integrals of `n1PDFMTRp3.CENSOR`), with its 5-point Simpson rule
over the winner's start point replaced by 4 panels x 16 Gauss-Legendre nodes
(`stageM4/dists_quad.R`), so the values are accurate to ~1e-12. Three parameter sets (the
posterior means of DMC's r94 subject-1 fits without and with guessing, rounded, and an
uncorrelated set with both lowest thresholds below A), both response orders, three ratings,
eight rts (144 densities), and 12 censored masses (`integrate()` at rel.tol 1e-10). Generator:
`make_dmc_p3_density.R`. The suite compares `dMTLBA(drifts = "any")` and the stopping-time terms
with it (relative 1e-8 in the body: the two quadratures).
