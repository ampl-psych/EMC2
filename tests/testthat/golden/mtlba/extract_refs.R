# Stage M0 (2026-09-30): extract the reference functions the MTLBA golden fixtures are built
# from, VERBATIM (source text via srcref, nothing retyped), into ref_extract/ so the generator
# scripts in this directory are self-contained.
#
# Run from this directory:  Rscript extract_refs.R [path to rating-work/mtlba/ref]
# The reference tree is read-only and is not part of the package.
args <- commandArgs(trailingOnly = TRUE)
ref <- if (length(args)) args[1] else
  "~/Documents/Projects/EMCpackage/rating-work/mtlba/ref"
ref <- normalizePath(ref)

extract <- function(file, funs, out, header) {
  ex <- parse(file, keep.source = TRUE)
  sr <- attr(ex, "srcref")
  nm <- vapply(ex, function(e)
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]))
      as.character(e[[2]]) else NA_character_, "")
  txt <- c(header, paste0("# md5 of source file: ", tools::md5sum(file)), "")
  for (f in funs) {
    i <- which(nm == f)
    if (!length(i)) stop("not found: ", f, " in ", file)
    i <- i[length(i)]  # the last definition is the one in force when the file is sourced
    txt <- c(txt, sprintf("# ---- %s (source lines %d-%d) ----", f, sr[[i]][1], sr[[i]][3]),
             as.character(sr[[i]]), "")
  }
  writeLines(txt, out)
  cat("wrote", out, "\n")
}

dir.create("ref_extract", showWarnings = FALSE)

# (1) Don't-know paper (https://osf.io/zxp5u), DMC code from DMC_DK_OSF.zip. The model the paper
#     reports (lba_DK_qD.R) lives in LBA-MTR/, whose n1PDF.DK also handles d < A.
extract(file.path(ref, "dk_osf/DMC_DK_OSF/rod/dmc/models/LBA-MTR/dists.R"),
  c("dlba.norm.core", "plba.norm.core", "dlba.norm", "plba.norm",
    "rlba.DK", "n1PDFfixedt0.DK", "n1PDF.DK"),
  "ref_extract/dk_LBA-MTR_dists_extract.R",
  c("# VERBATIM extract (extract_refs.R) from the don't-know paper's DMC code:",
    "#   DMC_DK_OSF/rod/dmc/models/LBA-MTR/dists.R  (https://osf.io/download/akzxw/)",
    "# Do not edit."))
# (2) The older, simpler version in MTR-LBA/ (valid only for d >= A); used as a cross-check.
extract(file.path(ref, "dk_osf/DMC_DK_OSF/rod/dmc/models/MTR-LBA/dists.R"),
  c("n1PDFfixedt0.DK", "n1PDF.DK"),
  "ref_extract/dk_MTR-LBA_dists_extract.R",
  c("# VERBATIM extract (extract_refs.R) from the don't-know paper's DMC code:",
    "#   DMC_DK_OSF/rod/dmc/models/MTR-LBA/dists.R  (https://osf.io/download/akzxw/)",
    "# Do not edit."))
# (3) EMC2 branch origin/BETE (commit 9eab6fd0, Oct 2025), R/model_MTR.R
extract(file.path(ref, "bete/model_MTR.R"),
  c("get_dL", "get_pgm", "guess_mtr", "rLBA_BE", "rLBA_TC", "recintab0", "n1PDF_MTR_1"),
  "ref_extract/bete_model_MTR_extract.R",
  c("# VERBATIM extract (extract_refs.R) from EMC2 branch origin/BETE (9eab6fd0), R/model_MTR.R.",
    "# recintab0 is itself a copy of MomTrunc:::recintab0 (CRAN package MomTrunc, GPL-2).",
    "# Do not edit."))
