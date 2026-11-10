# ---------------------------------------------------------------------------
# Executed 2026-10-24 from the repository root (R 4.5.2, lavaan 0.7.2). The
# first run stopped at an exact-symmetry check on gamma (rounding, 6.7e-17
# relative); that sanity check was made relative (< 1e-12) and the second run
# held every assertion. The order check: wls.obs equals its sampstat pieces, 0.
#
# dwlsmv fixture generation -- tier 2: lavaan's conditional Stage 1
# (conditional.x = TRUE), one binary x, single group. Item 4a of the
# conditional-covariate work order (E-006, 2026-10-24); the
# derivation is @sec-condx.
#
# Written 2026-10-24 by the main session.
#
# Data (tier-1 style, base R, seeded): fixtures/data/condx_A.csv, N = 1500,
#   y1, y2 (3 categories), y3 (binary), x in {0, 1}. Generated from
#   @sec-condx's population: lambda = (0.8, 0.7, 0.6), gamma = 0.5 (F's
#   disturbance variance 1), P(x = 1) = 0.4, thresholds at x = 0
#   (-0.5, 0.6), (-0.2, 0.9), (0.3).
#
# lavaan is called with do.fit = FALSE: only Stage 1 is wanted, and no
# fitted start or estimate enters it. Written:
#   lavaan_condx_s.csv      label, value       lavaan's wls.obs, in its order
#   lavaan_condx_gamma.csv  row, col, value    lavaan's gamma (N x ACOV of s)
# Labels are lavaan's own (y1|t1 ... y1~x ... y1~~y2 ...). Before writing,
# wls.obs is asserted equal to its sampstat pieces (res.th, res.slopes, the
# lower triangle of res.cov), so the order is proved, not assumed.
#
# Run from the repository root:
#   Rscript fixtures/generate/99h_lavaan_condx.R
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))
if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

DAT <- file.path("fixtures", "data", "condx_A.csv")
OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_lavaan_condx.json")
SEED <- 20261024L
TOL_ORDER <- 1e-15   # wls.obs against its sampstat pieces (set before running)

set.seed(SEED)
N <- 1500L
lambda <- c(0.8, 0.7, 0.6); gamma <- 0.5; q <- 0.4
thr <- list(c(-0.5, 0.6), c(-0.2, 0.9), 0.3)
x <- rbinom(N, 1L, q)
eta <- gamma * x + rnorm(N)
Y <- sapply(seq_along(lambda), function(j) {
  ystar <- lambda[j] * eta + rnorm(N, 0, sqrt(1 - lambda[j]^2))
  findInterval(ystar, thr[[j]])
})
dat <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = x)
# every (category, group) cell non-empty: no separation (@sec-condx (b))
for (j in 1:3) stopifnot(all(table(factor(dat[[j]], 0:(length(thr[[j]]))), dat$x) > 0))
write_csv_lf(dat, DAT)

fit <- lavaan::sem("f =~ y1 + y2 + y3\n f ~ x", data = dat,
                   ordered = c("y1", "y2", "y3"), estimator = "WLSMV",
                   parameterization = "delta", conditional.x = TRUE, std.lv = TRUE,
                   zero.add = c(0.5, 0.0), zero.keep.margins = TRUE, do.fit = FALSE)
stopifnot(isTRUE(lavInspect(fit, "options")$conditional.x))

s  <- lavInspect(fit, "wls.obs")
ss <- lavInspect(fit, "sampstat")
R  <- ss$res.cov
pieces <- c(ss$res.th, ss$res.slopes[, "x"], R[lower.tri(R)])
gap <- max(abs(unname(s) - unname(pieces)))
cat(sprintf("wls.obs against sampstat pieces: %.2e (tol %.0e)\n", gap, TOL_ORDER))
stopifnot(gap <= TOL_ORDER, length(s) == 11L)

G <- lavInspect(fit, "gamma")
stopifnot(identical(rownames(G), names(s)), identical(colnames(G), names(s)),
          all(is.finite(G)))
asym <- max(abs(G - t(G))) / max(abs(G))
cat(sprintf("gamma asymmetry (relative): %.2e
", asym))
stopifnot(asym < 1e-12)   # rounding only; a sanity check, not a comparison

f_s <- file.path(OUT, "lavaan_condx_s.csv")
f_g <- file.path(OUT, "lavaan_condx_gamma.csv")
write_csv_lf(data.frame(label = names(s), value = sprintf("%.17g", unname(s))), f_s)
ij <- expand.grid(row = seq_along(s), col = seq_along(s))
write_csv_lf(data.frame(row = names(s)[ij$row], col = names(s)[ij$col],
                        value = sprintf("%.17g", G[cbind(ij$row, ij$col)])), f_g)
write_manifest(MAN, c(DAT, f_s, f_g), seeds = as.character(SEED),
  extra = list(packages = c(paste("lavaan", as.character(packageVersion("lavaan"))),
    "WLSMV, delta, conditional.x = TRUE, std.lv, do.fit = FALSE, zero.add = c(0.5, 0), zero.keep.margins = TRUE",
    sprintf("wls.obs against sampstat pieces tolerance %.0e (set before running)", TOL_ORDER),
    R.version.string)))
cat(sprintf("%s, %s, %s and %s written\n", DAT, f_s, f_g, MAN))
