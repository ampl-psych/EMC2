# The NLE project's four SBC cells (evaluation/sbc_specs.R: one participant, factor S =
# left/right, `trials` per S cell, priors well inside the training boxes), built with
# nn_cell() from the registrations. Sourced by nle_phase6_screen.R / nle_phase6_bundle.R.
#
#   source(system.file("scripts", "nle_phase6_cells.R", package = "EMC2"))
#   cells <- nle_phase6_cells(reg_card = "~/R/nle-cards/ddm_reg_s402.rds")
#
# reg_s402 is not shipped (a comparison baseline, no normalisation guarantee): give the
# path of the card made by inst/scripts/regress_to_card.R, or NULL to skip it. Its cell
# has EXACTLY cap256w_c4's design and prior, so the cap256w_c4 analytic DDM control is
# also its control.

nle_phase6_cells <- function(reg_card = NULL) {
  # DDM prior of sbc_specs.R, by name (s held at 1); order as its formula
  ddm_formula <- list(v ~ 1, a ~ 1, t0 ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1, st0 ~ 1)
  ddm_mean <- c(v = 1, a = log(.8), t0 = log(.3), Z = qnorm(.5), sv = log(.1), SZ = qnorm(.2), st0 = log(.05))
  ddm_sd <- c(v = .15, a = .15, t0 = .1, Z = .05, sv = .15, SZ = .15, st0 = .15)
  ddm_cell <- function(model, formula = ddm_formula, mean = ddm_mean, sd = ddm_sd)
    nn_cell(model, mean = mean, sd = sd, formula = formula,
            factors = list(S = c("left", "right")), Rlevels = c("left", "right"))
  cells <- list(
    cap256w_c4 = ddm_cell(DDMnn),
    st0zero = ddm_cell(function() DDMnn("ddm_st0zero"), formula = ddm_formula[-7],
                       mean = ddm_mean[-7], sd = ddm_sd[-7]),
    rdm_small = nn_cell(RDMnn, factors = list(S = c("left", "right")), Rlevels = c("left", "right"),
                        matchfun = function(d) d$S == d$lR,
                        formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
                        mean = c(v = log(1.5), v_lMTRUE = .4, B = log(1), A = log(.3), t0 = log(.2)),
                        sd = .1))
  if (!is.null(reg_card)) {
    # no `model` field in the card: the analytic twin is EMC2's DDM (its dfun/pfun are not
    # used; the simulator and the parameter defaults are)
    cells$reg_s402 <- ddm_cell(register_nn_model(path.expand(reg_card), twin = DDM))
  }
  cells
}
