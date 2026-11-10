# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan on the marker-trap fixture
#
# Executed 2026-10-06 (lavaan 0.7.2): uv T3 59.938368 (64); lavaan's own marker
# fit converged to the ref solution (max gap 8.9e-7), i.e. lavaan is NOT trapped.
#
# E-016, "Marker-metric start trap and recovery through unit
# variance". The data are tier 1's fixtures/data/marker_trap.csv
# (95_marker_trap.R, fixture mode): MathAb by q1-q8 (binary), MathAtt by
# att1-att5 (5 categories, att3 and att4 reverse-keyed), MathAtt -> MathAb.
#
# Three lavaan 0.7.2 fits, all WLSMV / delta / scaled.shifted:
#   uv      unit variance: every loading free, psi(MathAtt) = 1 and the
#           DISTURBANCE variance psi(MathAb) = 1 (MathAb is endogenous)
#   marker  lavaan's default marker metric (att1, q1 at 1), lavaan's own starts
#   ref     the uv solution transformed to the marker metric, c_f = the
#           marker's uv loading: lambda/c_f, beta * c_B/c_A, psi_ff * c_f^2.
#           Computed here, not fitted. It is the reference marker solution
#           whether or not lavaan's own marker fit reaches it.
#
#   fixtures/expected/lavaan_marker_trap_sampstat.csv    label, value (s)
#   fixtures/expected/lavaan_marker_trap_uv.csv          label, est, se
#   fixtures/expected/lavaan_marker_trap_marker.csv      label, est, se
#   fixtures/expected/lavaan_marker_trap_ref.csv         label, est
#   fixtures/expected/lavaan_marker_trap_test.csv        label, value: T3 of
#       uv and marker, convergence of each, and the marker-vs-ref gap
#
# SELF-CHECKS, which stop the generator: the uv fit converged and is
# admissible (positive residual variances); the marker fit, whatever it did,
# is RECORDED, not required to match (a trapped lavaan is informative, not a
# reason to withhold the reference).
#
# Pinned options as 90_: bad.marker.crit = 0. Own manifest,
# MANIFEST_lavaan_marker_trap.json. LF on every platform (write_csv_lf).
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(lavaan))
source(file.path("fixtures", "generate", "00_common.R"))

DATA <- file.path("fixtures", "data", "marker_trap.csv")
OUT  <- file.path("fixtures", "expected")
MAN  <- file.path("fixtures", "MANIFEST_lavaan_marker_trap.json")
written <- character(0)

dat <- read.csv(DATA)
dat[] <- lapply(dat, ordered)
ord <- names(dat)
q <- sprintf("q%d", 1:8); a <- sprintf("att%d", 1:5)

fit_wlsmv <- function(model)
  lavaan::sem(model, data = dat, ordered = ord, estimator = "WLSMV",
              parameterization = "delta", test = "scaled.shifted",
              bad.marker.crit = 0)

m_uv <- paste0(
  "MathAb =~ ", paste0("NA*", q, collapse = " + "), "\n",
  "MathAtt =~ ", paste0("NA*", a, collapse = " + "), "\n",
  "MathAb ~ MathAtt\n",
  "MathAtt ~~ 1*MathAtt\n",
  "MathAb ~~ 1*MathAb\n")
m_mk <- paste0(
  "MathAb =~ ", paste(q, collapse = " + "), "\n",
  "MathAtt =~ ", paste(a, collapse = " + "), "\n",
  "MathAb ~ MathAtt\n")

fit_uv <- fit_wlsmv(m_uv)
stopifnot(lavInspect(fit_uv, "converged"))
th <- lavInspect(fit_uv, "theta")
stopifnot(all(diag(th) > 0))

fit_mk <- suppressWarnings(fit_wlsmv(m_mk))
mk_conv <- lavInspect(fit_mk, "converged")

lab <- function(pe) paste0(pe$lhs, pe$op, pe$rhs)
pe_uv <- parameterEstimates(fit_uv); pe_uv$label <- lab(pe_uv)
pe_mk <- parameterEstimates(fit_mk); pe_mk$label <- lab(pe_mk)
keep <- function(pe) pe[pe$op %in% c("=~", "~", "~~", "|") &
                          !(pe$op == "~~" & pe$lhs != pe$rhs) &
                          !(pe$op == "~~" & pe$lhs %in% ord), ]
pe_uv <- keep(pe_uv); pe_mk <- keep(pe_mk)

# --- the reference marker solution, transformed from uv --------------------------
est <- setNames(pe_uv$est, pe_uv$label)
cA <- est["MathAtt=~att1"]; cB <- est["MathAb=~q1"]
ref <- pe_uv[, c("label", "est")]
for (i in seq_len(nrow(ref))) {
  l <- ref$label[i]
  if (startsWith(l, "MathAb=~"))  ref$est[i] <- ref$est[i] / cB
  if (startsWith(l, "MathAtt=~")) ref$est[i] <- ref$est[i] / cA
  if (l == "MathAb~MathAtt")       ref$est[i] <- ref$est[i] * cB / cA
  if (l == "MathAtt~~MathAtt")     ref$est[i] <- cA^2
  if (l == "MathAb~~MathAb")       ref$est[i] <- cB^2
}
mk_est <- setNames(pe_mk$est, pe_mk$label)
gap <- max(abs(mk_est[ref$label] - ref$est))

t_uv <- fitMeasures(fit_uv, c("chisq.scaled", "df.scaled", "pvalue.scaled"))
t_mk <- fitMeasures(fit_mk, c("chisq.scaled", "df.scaled", "pvalue.scaled"))
cat(sprintf("uv: T3 %.6f df %d; marker: converged %s, T3 %.6f; max |marker - ref| %.3e\n",
            t_uv[1], as.integer(t_uv[2]), mk_conv, t_mk[1], gap))

s <- lavInspect(fit_uv, "wls.obs")
f <- file.path(OUT, "lavaan_marker_trap_sampstat.csv")
write_csv_lf(data.frame(label = names(s), value = as.numeric(s)), f); written <- c(written, f)
f <- file.path(OUT, "lavaan_marker_trap_uv.csv")
write_csv_lf(pe_uv[, c("label", "est", "se")], f); written <- c(written, f)
f <- file.path(OUT, "lavaan_marker_trap_marker.csv")
write_csv_lf(pe_mk[, c("label", "est", "se")], f); written <- c(written, f)
f <- file.path(OUT, "lavaan_marker_trap_ref.csv")
write_csv_lf(ref, f); written <- c(written, f)
f <- file.path(OUT, "lavaan_marker_trap_test.csv")
write_csv_lf(data.frame(
  label = c("uv.chisq.scaled", "uv.df.scaled", "uv.pvalue.scaled", "uv.converged",
            "marker.chisq.scaled", "marker.df.scaled", "marker.pvalue.scaled",
            "marker.converged", "marker.max_abs_gap_to_ref"),
  value = c(as.numeric(t_uv), 1, as.numeric(t_mk), as.numeric(mk_conv), gap)), f)
written <- c(written, f)

write_manifest(MAN, c(DATA, written),
  seeds = "none: deterministic given fixtures/data/marker_trap.csv",
  extra = list(packages = c(
    paste("lavaan", as.character(packageVersion("lavaan"))),
    "WLSMV, delta, scaled.shifted, bad.marker.crit = 0 (pinned)",
    sprintf("uv T3 %.6f; marker converged %s, T3 %.6f, max gap to ref %.3e",
            t_uv[1], mk_conv, t_mk[1], gap),
    R.version.string)))
cat(sprintf("%d file(s) written; %s written.\n", length(written), MAN))
