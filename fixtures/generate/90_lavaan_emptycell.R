# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan on data with a one-empty-cell pair
#
# O6 (decided and built 2026-10-06; spec Module 2, *One-empty-cell repair*;
# derivation @sec-o6). dwlsmv repairs a binary pair whose 2x2 table has
# exactly one empty cell as lavaan 0.7.2 does by default (zero.add = c(0.5, 0),
# zero.keep.margins = TRUE) and builds Gamma by construction (a), lavaan's.
# So lavaan's DEFAULT fit is the reference here, end to end: s, Gamma,
# theta-hat, SEs and T3. No other committed fixture contains such a pair.
#
# THE DATA are generated here, not by 10_fixtures.R, because they exist only
# for this reference: one factor, four binary items, N = 400, seeded. y1 is
# rare and y2 is set to 1 wherever y1 is 1, so the y1~~y2 table has exactly
# one empty cell (y1 = 1, y2 = 0) -- the same construction as
# tests/unit/test_dwlsmv.do section 10, but R's generator, so not the same
# numbers. Asserted below: that table has exactly one empty cell, no other
# pair has any, and no category is empty.
#
#   fixtures/data/emptycell.csv                 y1..y4, 0/1
#   fixtures/expected/lavaan_emptycell_sampstat.csv  label, value: wls.obs (s)
#   fixtures/expected/lavaan_emptycell_gamma.csv     label + D columns: gamma
#   fixtures/expected/lavaan_emptycell_parameters.csv label, est, se
#   fixtures/expected/lavaan_emptycell_test.csv       chisq.scaled, df.scaled,
#                                                     pvalue.scaled
#
# SELF-CHECKS, which stop the generator rather than write a fixture that does
# not test what it claims:
#   - lavaan's rho(y1, y2) is the keep-margins c = 0.5 repaired root: the
#     saturated inversion Phi2(a, b; rho) = (n11 + 0.5 sigma_e)/N, to 1e-8 (the
#     spec 7.2 correlation tier; lavaan's numerical fit, measured 1.3e-10 on
#     T1 in @sec-o6 level O5);
#   - without repair (zero.add = c(0, 0)) the pair's estimate is different
#     (by more than 1e-2, chosen after measuring 0.18), so the repair is what
#     this fixture exercises.
#
# Pinned options, as 20_ and 70_: estimator, parameterization, test,
# zero.add, zero.keep.margins, and bad.marker.crit = 0 (70_'s reason: lavaan
# otherwise may move the marker off a poor first indicator).
#
# Own manifest (MANIFEST_lavaan_emptycell.json), as 70_ and 80_ have.
# LF on every platform (write_csv_lf), as 80_.
#
# Run from the repository root:  Rscript fixtures/generate/90_lavaan_emptycell.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT  <- file.path("fixtures", "expected")
DATA <- file.path("fixtures", "data", "emptycell.csv")
MAN  <- file.path("fixtures", "MANIFEST_lavaan_emptycell.json")
SEED <- 20261006
written <- character(0)

# --- the data -----------------------------------------------------------------
set.seed(SEED)
N  <- 400
f  <- rnorm(N)
z  <- sapply(1:4, function(j) 0.7 * f + sqrt(1 - 0.49) * rnorm(N))
y1 <- as.integer(z[, 1] > 1)
y2 <- as.integer(z[, 2] > -0.5)
y2[y1 == 1] <- 1L
y3 <- as.integer(z[, 3] > 0)
y4 <- as.integer(z[, 4] > 0.3)
dat <- data.frame(y1 = y1, y2 = y2, y3 = y3, y4 = y4)

stopifnot(all(sapply(dat, function(v) all(table(factor(v, 0:1)) > 0))))
tab <- function(a, b) table(factor(dat[[a]], 0:1), factor(dat[[b]], 0:1))
nz  <- outer(1:4, 1:4, Vectorize(function(i, j) if (i < j) sum(tab(i, j) == 0) else 0L))
stopifnot(nz[1, 2] == 1, sum(nz) == 1)
n12 <- unclass(tab(1, 2))
cat("y1~~y2 table:", n12, "(column-major)\n")

write_csv_lf(dat, DATA)
written <- c(written, DATA)

# --- the fit ------------------------------------------------------------------
dato <- dat
dato[] <- lapply(dato, ordered)
fit_ec <- function(zadd)
  lavaan::sem("f =~ y1 + y2 + y3 + y4", data = dato, ordered = names(dato),
              estimator = "WLSMV", parameterization = "delta",
              test = "scaled.shifted", zero.add = zadd,
              zero.keep.margins = TRUE, bad.marker.crit = 0)
fit <- fit_ec(c(0.5, 0))
stopifnot(isTRUE(lavInspect(fit, "options")$zero.keep.margins),
          lavInspect(fit, "converged"))

# --- self-check 1: lavaan repaired y1~~y2 by keep-margins, c = 0.5 --------------
Phi2 <- function(a, b, r)
  integrate(function(x) dnorm(x) * pnorm((b - r * x) / sqrt(1 - r^2)),
            -Inf, a, rel.tol = 1e-12)$value
sg  <- matrix(c(1, -1, -1, 1), 2)
np  <- n12 + 0.5 * sg[which(n12 == 0)] * sg
a   <- qnorm(sum(n12[1, ]) / N)
b   <- qnorm(sum(n12[, 1]) / N)
rp  <- uniroot(function(r) Phi2(a, b, r) - np[1, 1] / N, c(-0.999999, 0.999999),
               tol = 1e-14)$root
rl  <- lavInspect(fit, "sampstat")$cov["y2", "y1"]
cat(sprintf("rho(y1,y2): lavaan %.12f, keep-margins inversion %.12f, gap %.2e\n",
            rl, rp, abs(rl - rp)))
stopifnot(abs(rl - rp) < 1e-8)

# --- self-check 2: the repair moves the estimate ------------------------------
# First written as "without repair the pair is at the bound (|rho| > 0.99)".
# That failed: on this table lavaan's unrepaired fit stops at 0.9698, short of
# its 0.999 cap (unlike L9d), on a likelihood that has no interior maximum.
# What the fixture needs is only that the repair acts, so the check is that the
# two fits differ; the floor 1e-2 is CHOSEN AFTER MEASURING (gap 0.18).
r0 <- suppressWarnings(lavInspect(fit_ec(c(0, 0)), "sampstat")$cov["y2", "y1"])
cat(sprintf("rho(y1,y2) without repair: %.6f (repaired %.6f)\n", r0, rl))
stopifnot(abs(r0 - rl) > 1e-2)

# --- the reference objects ----------------------------------------------------
s <- lavInspect(fit, "wls.obs")
f <- file.path(OUT, "lavaan_emptycell_sampstat.csv")
write_csv_lf(data.frame(label = names(s), value = as.numeric(s)), f)
written <- c(written, f)

G <- lavInspect(fit, "gamma")
stopifnot(identical(rownames(G), names(s)))
f <- file.path(OUT, "lavaan_emptycell_gamma.csv")
write_csv_lf(cbind(label = rownames(G), as.data.frame(unclass(G))), f)
written <- c(written, f)

se_type <- lavInspect(fit, "options")$se
pe <- parameterEstimates(fit)
pe$label <- paste0(pe$lhs, pe$op, pe$rhs)
f <- file.path(OUT, "lavaan_emptycell_parameters.csv")
write_csv_lf(pe[, c("label", "est", "se")], f)
written <- c(written, f)

ft <- fitMeasures(fit, c("chisq.scaled", "df.scaled", "pvalue.scaled"))
f <- file.path(OUT, "lavaan_emptycell_test.csv")
write_csv_lf(data.frame(label = names(ft), value = as.numeric(ft)), f)
written <- c(written, f)

write_manifest(MAN, written,
  seeds = paste("set.seed", SEED, "(data); fit deterministic given the data"),
  extra = list(packages = c(
    paste("lavaan", as.character(packageVersion("lavaan"))),
    paste("se type:", se_type),
    "zero.add = c(0.5, 0), zero.keep.margins = TRUE, bad.marker.crit = 0 (pinned)",
    sprintf("self-checks: rho(y1,y2) vs keep-margins inversion %.2e; without repair %.6f (repaired %.6f)", abs(rl - rp), r0, rl),
    R.version.string)))

cat(sprintf("%d file(s) written; %s written.\n", length(written), MAN))
