# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- the .dta transfer channel
#
# Spec Section 7.1 step 2 requires the R-to-Stata file channel to be verified
# lossless BEFORE anything downstream is trusted, rather than assumed. The
# replicate machinery writes .dta because both tools read it natively, and a
# silent precision loss there would corrupt every paired comparison in the
# validation design while looking like an estimator disagreement.
#
# What is actually at risk, stated precisely. The ordinal payload is
# integer-valued, so it is nearly impossible to damage. The exposure is future:
# Section 7.1 records per-replicate quantities alongside the data, and the
# moment a non-integer column crosses this channel, "it worked for the ordinal
# data" stops being evidence. So this checks both:
#
#   case_D.dta          the real payload, integer ordinal
#   roundtrip_probe.dta doubles chosen to break a lossy channel
#
# The probe carries the values a careless format loses: full 17-digit doubles,
# subnormals, values near the double limits, negative zero, and fractions with
# no finite binary expansion. If the channel survives those it will survive a
# simulated dataset.
#
# WHY THE STATA SIDE IS A SEPARATE TEST. An R-to-R round-trip proves only that
# haven agrees with itself. The question is whether STATA reads back what R
# wrote, so the assertion that matters lives in tests/unit/test_dta_roundtrip.do
# and compares the .dta against the committed .csv inside Stata. This script
# only produces the files and checks the R half.
#
# The .dta files are committed so the Stata suite stays R-free, a project rule.
#
# Run from the repository root:  Rscript fixtures/generate/40_dta_roundtrip.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("haven", quietly = TRUE))
  stop("haven is required to write .dta fixtures.")

OUT_DATA <- file.path("fixtures", "data")
written <- character(0)

stopifnot_msg <- function(ok, msg) if (!ok) stop(msg, call. = FALSE)

# --- 1. the real payload ----------------------------------------------------
csv <- read.csv(file.path(OUT_DATA, "case_D.csv"))
dta <- file.path(OUT_DATA, "case_D.dta")
haven::write_dta(csv, dta)
back <- as.data.frame(haven::read_dta(dta))

stopifnot_msg(identical(dim(back), dim(csv)),
              "case_D: dimensions changed across the .dta round-trip")
stopifnot_msg(identical(names(back), names(csv)),
              "case_D: variable names changed across the .dta round-trip")
for (v in names(csv)) {
  stopifnot_msg(all(as.numeric(back[[v]]) == csv[[v]]),
                sprintf("case_D: column %s changed in the round-trip", v))
}
written <- c(written, dta)

# --- 2. the precision probe -------------------------------------------------
# Deliberately nastier than anything the study will produce.
#
# THE CHANNEL'S RANGE IS NARROWER THAN IEEE'S, discovered by this probe on its
# first run. Stata's largest double is 8.9884656743115785e+307, which in hex is
# +1.fffffffffffffX+3fe -- exactly one binade below IEEE's +1.fffffffffffffX+3ff.
# Stata reserves that entire top binade for its missing values (., .a through
# .z), so it gives up half the exponent range at the top to get them.
#
# .Machine$double.xmax therefore cannot cross this channel at all. haven
# refuses it outright:
#
#     Failed to insert value [3, 2]: A provided numeric value was outside the
#     range of representable values in the specified file format.
#
# That is the good outcome and worth stating: the channel fails LOUDLY at its
# limit rather than silently truncating, so an out-of-range value cannot become
# a wrong number further downstream. Nothing this project computes -- thresholds,
# correlations, chi-squares -- comes near 1e307, so the limit is academic here;
# it is recorded because "the transfer is lossless" is false as stated, and the
# next person to carry something unusual through it should know where the edge
# is.
#
# The probe now sits exactly ON the boundary rather than past it: maxdouble and
# its negative, and smallestdouble, which is the smallest NORMAL double. Stata
# advertises no subnormal support, so subnormals are deliberately not probed and
# should not be assumed to survive.
stata_maxdouble   <- 8.9884656743115785e+307   # c(maxdouble)
stata_smallestdbl <- 2.2250738585072014e-308   # c(smallestdouble)

probe <- data.frame(
  v1 = c(1 / 3, pi, exp(1), sqrt(2), 0.1),
  v2 = c(.Machine$double.eps, stata_smallestdbl, stata_maxdouble,
         -0.0, 1e-300),
  v3 = c(-1 / 7, 1234567890.123456, -9.87654321e-12, 2^53, 2^53 + 2),
  v4 = c(0, -1, -stata_maxdouble, 1e-300, 0.30000000000000004)
)

probe_dta <- file.path(OUT_DATA, "roundtrip_probe.dta")
haven::write_dta(probe, probe_dta)
probe_back <- as.data.frame(haven::read_dta(probe_dta))

for (v in names(probe)) {
  a <- probe[[v]]
  b <- as.numeric(probe_back[[v]])
  # identical() on doubles is the right test here: this is a bitwise question,
  # not a numerical-agreement one. A tolerance would defeat the purpose.
  stopifnot_msg(identical(a, b),
                sprintf("probe: column %s not bitwise after round-trip", v))
}
written <- c(written, probe_dta)

# The same values as text, at 17 significant digits, which round-trips an IEEE
# double exactly. This is what the Stata test compares the .dta against: a
# channel that agrees with itself proves nothing, so the reference has to reach
# Stata by a different route than the .dta does.
probe_csv <- file.path(OUT_DATA, "roundtrip_probe.csv")
txt <- as.data.frame(lapply(probe, function(x) sprintf("%.17g", x)),
                     stringsAsFactors = FALSE)
write.csv(txt, probe_csv, row.names = FALSE, quote = FALSE)

# Confirm the text route is itself exact before trusting it as the reference.
txt_back <- read.csv(probe_csv, colClasses = "character")
for (v in names(probe)) {
  stopifnot_msg(identical(as.numeric(txt_back[[v]]), probe[[v]]),
                sprintf("probe: column %s not exact via text route", v))
}
written <- c(written, probe_csv)

write_manifest(file.path("fixtures", "MANIFEST_dta.json"), written,
               seeds = "none (deterministic given case_D.csv and fixed probe values)",
               extra = list(packages = c(
                 paste("haven", as.character(utils::packageVersion("haven"))),
                 R.version.string)))

cat(sprintf("%d transfer-channel files written.\n", length(written)))
cat("R half verified: bitwise identical across write_dta/read_dta, and the\n")
cat("%.17g text route is exact. The Stata half is tests/unit/test_dta_roundtrip.do\n")
cat("-- an R-to-R round-trip only proves haven agrees with itself.\n")
