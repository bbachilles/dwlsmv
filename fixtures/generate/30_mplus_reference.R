# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 3: Mplus reference outputs
#
# Executed 2026-09-18 with Mplus 8.x on the development machine.
#
# Why a second reference at all. The v1 target statistic is T3, the
# scale-and-shift, and the spec records that it is BOTH lavaan's
# "scaled.shifted" and the Mplus WLSMV default (Gate 5). That claim rested on
# reading Asparouhov & Muthen (2010) Table 1 -- in a spec that had already
# caught itself citing the wrong source for this very quantity. Running both
# tools on the same data turns a reading of a paper into a measurement.
#
# Identification is deliberately matched to tier 2. lavaan's default is the
# marker method (first loading fixed to 1), so this fixes y1's loading to 1 and
# frees the factor variance rather than the reverse. Mplus's default is the same
# marker method (stated by the user, 2026-10-06; not measured); fixing the factor
# variance instead (fac@1) is an equivalent reparameterization, but the
# estimates are then NOT comparable element by element, which is an easy and
# very confusing mistake to make. Do not "simplify" this back to fac@1.
#
# PRECISION -- the thing to understand before using these files.
#
# Mplus prints estimates to three decimals in the .out. That is coarser than
# anything else this project compares against, so SAVEDATA RESULTS is used as
# well: it writes an E-format block with eight significant digits. The two
# sources are cross-checked against each other here rather than trusted
# separately.
#
# The RESULTS block is POSITIONAL and its layout is an Mplus implementation
# detail, not a documented stable contract -- in this run it is thresholds
# first, then loadings, whereas the parsed .out table lists loadings first.
# So the raw block is committed verbatim, and only the test statistic is
# extracted from it, by locating the entry that matches the parsed chi-square.
# That match is asserted, so a layout change breaks this script loudly instead
# of silently mislabelling a number.
#
# Needs: Mplus on PATH or at the default Windows install, and MplusAutomation.
# Run from the repository root:  Rscript fixtures/generate/30_mplus_reference.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("MplusAutomation", quietly = TRUE))
  stop("MplusAutomation is required for tier-3 fixtures.")

mplus_exe <- Sys.which("Mplus")
if (!nzchar(mplus_exe)) {
  cand <- "C:/Program Files/Mplus/Mplus.exe"
  if (file.exists(cand)) mplus_exe <- cand
}
if (!nzchar(mplus_exe) || !file.exists(mplus_exe))
  stop("Mplus executable not found. Tier 3 is optional; ",
       "tiers 1 and 2 do not need it.")

OUT <- file.path("fixtures", "expected")
written <- character(0)

# --- run Mplus in a scratch directory --------------------------------------
# Mplus resolves FILE IS relative to the input file and has historically been
# awkward with long absolute paths, so everything runs in one flat temp dir
# with bare filenames.
work <- file.path(tempdir(), "dwlsmv_mplus")
dir.create(work, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(file.path("fixtures", "data", "case_D.csv"))
write.table(dat, file.path(work, "case_D.dat"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

inp <- c(
  "TITLE: dwlsmv tier-3 reference -- case D, single factor, WLSMV",
  "DATA: FILE IS case_D.dat;",
  "VARIABLE:",
  "  NAMES ARE y1 y2 y3 y4;",
  "  CATEGORICAL ARE y1 y2 y3 y4;",
  "ANALYSIS:",
  "  ESTIMATOR = WLSMV;",
  "  PARAMETERIZATION = DELTA;",
  "MODEL:",
  # marker method, to match lavaan's default -- see the header
  "  fac BY y1@1 y2 y3 y4;",
  "  fac;",
  "OUTPUT: SAMPSTAT;",
  "SAVEDATA: RESULTS ARE results.dat;"
)
writeLines(inp, file.path(work, "case_D.inp"))

owd <- setwd(work)
status <- system2(mplus_exe, "case_D.inp", stdout = TRUE, stderr = TRUE)
setwd(owd)

# Mplus lowercases the output filename.
outfile <- file.path(work, "case_d.out")
if (!file.exists(outfile))
  stop("Mplus produced no output. Tail of its console output:\n",
       paste(utils::tail(status, 10), collapse = "\n"))

m <- MplusAutomation::readModels(outfile, quiet = TRUE)

if (length(m$errors))
  stop("Mplus reported errors:\n", paste(unlist(m$errors), collapse = "\n"))

# --- labelled estimates, three decimals, from the parsed .out ---------------
pars <- m$parameters$unstandardized
write.csv(pars, file.path(OUT, "mplus_case_D_parameters.csv"),
          row.names = FALSE, quote = FALSE)
written <- c(written, file.path(OUT, "mplus_case_D_parameters.csv"))

# --- polychoric correlations ------------------------------------------------
cr <- m$sampstat$correlations.vardiag
write.csv(cbind(label = rownames(cr), as.data.frame(cr)),
          file.path(OUT, "mplus_case_D_correlations.csv"),
          row.names = FALSE, quote = FALSE)
written <- c(written, file.path(OUT, "mplus_case_D_correlations.csv"))

# --- the high-precision RESULTS block, verbatim ------------------------------
res <- scan(file.path(work, "results.dat"), quiet = TRUE)
write.csv(data.frame(position = seq_along(res), value = res),
          file.path(OUT, "mplus_case_D_results_raw.csv"),
          row.names = FALSE, quote = FALSE)
written <- c(written, file.path(OUT, "mplus_case_D_results_raw.csv"))

# --- the test statistic, located rather than assumed ------------------------
chi_printed <- m$summaries$ChiSqM_Value
df_printed  <- m$summaries$ChiSqM_DF
p_printed   <- m$summaries$ChiSqM_PValue

hit <- which(abs(res - chi_printed) < 5e-4)
if (length(hit) != 1L)
  stop("could not locate the chi-square in the RESULTS block unambiguously ",
       "(", length(hit), " candidates). The block layout has ",
       "probably changed; ",
       "inspect mplus_case_D_results_raw.csv rather than trusting a position.")

chi_precise <- res[hit]
# df and p sit immediately after it in every layout seen so far -- asserted,
# not assumed, against the values parsed from the printed output.
if (abs(res[hit + 1L] - df_printed) > 1e-9)
  stop("the entry after the chi-square is not the degrees of freedom.")
if (abs(res[hit + 2L] - p_printed) > 5e-4)
  stop("the entry two after the chi-square is not the p-value.")

test <- data.frame(
  label = c("chisq.wlsmv", "df", "pvalue"),
  value = c(chi_precise, res[hit + 1L], res[hit + 2L])
)
write.csv(test, file.path(OUT, "mplus_case_D_test.csv"),
          row.names = FALSE, quote = FALSE)
written <- c(written, file.path(OUT, "mplus_case_D_test.csv"))

write_manifest(file.path("fixtures", "MANIFEST_mplus.json"), written,
               seeds = "none (deterministic given case_D.csv)",
               extra = list(packages = c(
                 paste("Mplus", m$summaries$Mplus.version),
                 paste("MplusAutomation",
                       as.character(utils::packageVersion("MplusAutomation"))),
                 R.version.string)))

cat(sprintf("%d Mplus reference files written.\n", length(written)))
cat(sprintf("chi-square: printed %.3f, precise %.8f, df %d\n",
            chi_printed, chi_precise, df_printed))
cat("Mplus prints three decimals; the precise value comes from SAVEDATA\n")
cat("RESULTS. Compare against tier 2 with that in mind -- see\n")
cat("the development record for the measured lavaan-vs-Mplus gap here.\n")
