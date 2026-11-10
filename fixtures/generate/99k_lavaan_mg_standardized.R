# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan standardizedSolution() at G = 2
#
# Module 11 D7 at G > 1 (spec 4.8 O17; derivation @sec-standardize-mg). The
# reference for the frozen map (P, Q-B (b)) and its delta-method SEs. Designs
# on committed data, no simulation here:
#   E    K = 4,4,4,2,2  F -> y1-y5           fixtures/data/mg_E.csv
#   B    K = 3,3,3,3,3  F -> y1-y5           fixtures/data/mg_B.csv
#   D    K = 4,2,2,2,2  F -> y1-y5           fixtures/data/mg_D.csv
#   B2F  K = 3,3,3,3,3  F -> y1-y3, H -> y4 y5, cov(F,H) free, mg_B.csv
# each at the four levels. B2F is here because a factor covariance row
# standardizes by two scales, which no one-factor design can test.
#
# Fits: the explicit level syntax of 98_ (never group.equal), its pinned
# options. Keys "ptype:group:row:col" as 98_: ptype 1 lambda (item, factor),
# 2 psi (factor, factor), 4 tau (item, k), 5 nu (item, 1), 6 omega (item, 1),
# 7 alpha (factor, 1), plus 8 theta (item, 1) in the standardized file only.
#
# TRAP -- the SEs carry lavaan's forward-difference error (80_'s header), about
# 1e-8 relative (derivation: 2.84e-8 at most). est_std has no such error.
#
# Outputs (fixtures/expected/):
#   lavaan_mg_std.csv         design, level, key, est_std, se
#   lavaan_mg_std_theta.csv   design, level, key, est    (lavaan's coef())
#   lavaan_mg_std_vcov.csv    design, level, key_i, key_j, value (vcov())
#
# Assertions before anything is written: converged; the identity
# lambda*^2 + Theta* = 1 on one-factor designs, 1e-12 absolute, which holds for
# the map, not for lavaan's arithmetic in particular (ID-P, @sec-standardize-mg).
#
# Run from the repository root:
#   Rscript fixtures/generate/99k_lavaan_mg_standardized.R
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))
if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_lavaan_mg_std.json")
LEVELS <- c("configural", "thresholds", "loadings", "intercepts")
p <- 5L
y <- paste0("y", seq_len(p))

DESIGNS <- list(
  E   = list(file = "mg_E.csv", K = c(4L, 4L, 4L, 2L, 2L), fac = list(f = 1:5)),
  B   = list(file = "mg_B.csv", K = rep(3L, 5L),           fac = list(f = 1:5)),
  D   = list(file = "mg_D.csv", K = c(4L, 2L, 2L, 2L, 2L), fac = list(f = 1:5)),
  B2F = list(file = "mg_B.csv", K = rep(3L, 5L),           fac = list(f = 1:3, h = 4:5)))

lev_syntax <- function(ds, level) {
  K <- ds$K; poly <- K > 2L
  tied_tau <- level != "configural"
  tied_lam <- level %in% c("loadings", "intercepts")
  om <- switch(level, configural = rep("fix", p), thresholds = ifelse(poly, "free", "fix"), rep("free", p))
  nu <- if (level %in% c("configural", "intercepts")) rep("fix", p) else rep("free", p)
  out <- character(0)
  for (fn in names(ds$fac)) {
    ii <- ds$fac[[fn]]
    ld <- if (tied_lam) sprintf("c(l%d,l%d)*%s", ii, ii, y[ii]) else sprintf("c(NA,NA)*%s", y[ii])
    out <- c(out, paste(fn, "=~", paste(ld, collapse = " + ")))
    phi <- level == "intercepts" || (level == "loadings" && any(poly[ii]))
    out <- c(out, sprintf("%s ~~ c(1,%s)*%s", fn, if (phi) "NA" else "1", fn),
             sprintf("%s ~ c(0,%s)*1", fn, if (level == "intercepts") "NA" else "0"))
  }
  if (length(ds$fac) == 2L) out <- c(out, "f ~~ c(NA,NA)*h")
  if (tied_tau) for (j in seq_len(p)) {
    cc <- seq_len(K[j] - 1L)
    out <- c(out, sprintf("%s | %s", y[j], paste(sprintf("c(t%d_%d,t%d_%d)*t%d", j, cc, j, cc, cc), collapse = " + ")))
  }
  paste(c(out, sprintf("%s ~*~ c(1,%s)*%s", y, ifelse(om == "free", "NA", "1"), y),
          sprintf("%s ~ c(0,%s)*1", y, ifelse(nu == "free", "NA", "0"))), collapse = "\n")
}

one_key <- function(op, lhs, rhs, g, fn) {
  fi <- function(x) match(x, fn)
  switch(op,
    "=~"  = sprintf("1:%d:%d:%d", g, match(rhs, y), fi(lhs)),
    "~~"  = if (lhs %in% fn && rhs %in% fn) { ij <- sort(c(fi(lhs), fi(rhs)), decreasing = TRUE); sprintf("2:%d:%d:%d", g, ij[1], ij[2]) }
            else if (lhs %in% y && lhs == rhs) sprintf("8:%d:%d:1", g, match(lhs, y)) else NA_character_,
    "|"   = sprintf("4:%d:%d:%d", g, match(lhs, y), as.integer(sub("t", "", rhs))),
    "~1"  = if (lhs %in% y) sprintf("5:%d:%d:1", g, match(lhs, y)) else sprintf("7:%d:%d:1", g, fi(lhs)),
    "~*~" = sprintf("6:%d:%d:1", g, match(lhs, y)),
    NA_character_)
}

std_rows <- list(); th_rows <- list(); v_rows <- list()
for (dn in names(DESIGNS)) {
  ds <- DESIGNS[[dn]]; fn <- names(ds$fac)
  dat <- read.csv(file.path("fixtures", "data", ds$file))
  dat[y] <- lapply(dat[y], ordered)
  for (lv in LEVELS) {
    fit <- lavaan::cfa(lev_syntax(ds, lv), data = dat, group = "g", ordered = y, estimator = "WLSMV",
                       parameterization = "delta", std.lv = TRUE, test = "scaled.shifted",
                       zero.add = c(0.5, 0.0), zero.keep.margins = TRUE, bad.marker.crit = 0)
    stopifnot(lavInspect(fit, "converged"))
    pt <- parTable(fit); fr <- pt[pt$free > 0, ]; fr <- fr[order(fr$free), ]
    kk <- vapply(seq_len(nrow(fr)), function(k) one_key(fr$op[k], fr$lhs[k], fr$rhs[k], fr$group[k], fn), "")
    stopifnot(!anyNA(kk))
    est <- coef(fit); V <- unname(vcov(fit))
    stopifnot(length(est) == length(kk), nrow(V) == length(kk), all(is.finite(est)), all(is.finite(V)))
    th_rows[[length(th_rows) + 1L]] <- data.frame(design = dn, level = lv, key = kk, est = unname(est))
    v_rows[[length(v_rows) + 1L]] <- data.frame(design = dn, level = lv, key_i = rep(kk, length(kk)),
                                                key_j = rep(kk, each = length(kk)), value = as.vector(V))
    ss <- standardizedSolution(fit, type = "std.all")
    ss <- ss[ss$op %in% c("=~", "|", "~~", "~1", "~*~"), ]
    sk <- vapply(seq_len(nrow(ss)), function(k) one_key(ss$op[k], ss$lhs[k], ss$rhs[k], ss$group[k], fn), "")
    stopifnot(!anyNA(sk), !anyDuplicated(sk))
    se <- ss$se; se[!is.finite(se)] <- 0
    stopifnot(all(is.finite(ss$est.std)))
    if (length(fn) == 1L) for (g in 1:2) for (j in seq_len(p)) {
      l <- ss$est.std[sk == sprintf("1:%d:%d:1", g, j)]; t <- ss$est.std[sk == sprintf("8:%d:%d:1", g, j)]
      stopifnot(abs(l^2 + t - 1) < 1e-12)
    }
    std_rows[[length(std_rows) + 1L]] <- data.frame(design = dn, level = lv, key = sk, est_std = ss$est.std, se = se)
    cat(sprintf("%-4s %-10s npar %d, %d standardized rows\n", dn, lv, length(kk), length(sk)))
  }
}

written <- character(0)
f <- file.path(OUT, "lavaan_mg_std.csv");       write_csv_lf(do.call(rbind, std_rows), f); written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_std_theta.csv"); write_csv_lf(do.call(rbind, th_rows), f);  written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_std_vcov.csv");  write_csv_lf(do.call(rbind, v_rows), f);   written <- c(written, f)
write_manifest(MAN, written, seeds = "none (deterministic given mg_E, mg_B, mg_D)",
  extra = list(packages = c(paste("lavaan", as.character(packageVersion("lavaan"))),
    "WLSMV, delta, scaled.shifted, std.lv, zero.add = c(0.5, 0), zero.keep.margins = TRUE, bad.marker.crit = 0 (pinned)",
    "explicit level syntax, never group.equal; standardizedSolution(type = 'std.all')",
    R.version.string)))
cat(sprintf("%d file(s) written; %s written.\n", length(written), MAN))
