# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- shared reference implementation
#
# This file is the R-side reference implementation of the Stage 1 quantities
# derived in spec/derivations/dwlsmv-derivations.qmd. It deliberately has NO
# package dependencies: Phi2 is computed from its one-dimensional integral
# representation rather than via mvtnorm, so these fixtures can be regenerated
# on any R installation and reproduce bit-for-bit.
#
# Everything here is verified in the derivation file. Do not change a formula
# here without changing it there, and re-running its verification chunks.
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

# Bivariate normal CDF, by the integral representation
#   Phi2(a,b;rho) = int_{-inf}^{a} phi(x) Phi((b - rho x)/sqrt(1-rho^2)) dx
# Accurate to about 1e-12.
Phi2 <- function(a, b, rho) {
  if (a == -Inf || b == -Inf) return(0)
  if (a ==  Inf && b ==  Inf) return(1)
  if (a ==  Inf) return(pnorm(b))
  if (b ==  Inf) return(pnorm(a))
  s <- sqrt(1 - rho^2)
  integrate(function(x) dnorm(x) * pnorm((b - rho * x) / s),
            -Inf, a, rel.tol = 1e-12)$value
}

# Bivariate normal density (Plackett: dPhi2/drho = phi2)
phi2 <- function(a, b, rho)
  exp(-(a^2 - 2 * rho * a * b + b^2) / (2 * (1 - rho^2))) /
  (2 * pi * sqrt(1 - rho^2))

# Threshold influence function for one observation in category `y`, given the
# cumulative probabilities `q` of one variable. Two values per threshold,
# switching sign exactly once at c = y (the "staircase" structure).
psi_tau <- function(y, q) ((seq_along(q) > y) - q) / dnorm(qnorm(q))

# The threshold block of Gamma.
#   Gamma[(j,c),(k,d)] = (P(Z_jc, Z_kd) - q_jc q_kd) / (phi_jc phi_kd)
# with the joint being min(q_jc, q_jd) within a variable and Phi2 across.
# No entry ever needs more than a bivariate normal CDF.
gamma_tau <- function(tau, q, var_id, rho_mat) {
  P <- length(tau)
  outer(seq_len(P), seq_len(P), Vectorize(function(a, b) {
    joint <- if (var_id[a] == var_id[b]) min(q[a], q[b])
             else Phi2(tau[a], tau[b], rho_mat[var_id[a], var_id[b]])
    (joint - q[a] * q[b]) / (dnorm(tau[a]) * dnorm(tau[b]))
  }))
}

# Labels matching lavaan's convention, which Section 4.2 of the spec adopts:
# thresholds "y1|t1" by variable then ascending index, correlations "y1~~y2"
# in column-major lower-triangle order.
threshold_labels <- function(n_thresh_per_var)
  unlist(lapply(seq_along(n_thresh_per_var), function(j)
    sprintf("y%d|t%d", j, seq_len(n_thresh_per_var[j]))))

correlation_labels <- function(p) {
  ij <- which(lower.tri(matrix(0, p, p)), arr.ind = TRUE)
  ij <- ij[order(ij[, "col"], ij[, "row"]), , drop = FALSE]
  sprintf("y%d~~y%d", ij[, "col"], ij[, "row"])
}

# Draw ordinal indicators from a single-factor model with equal loadings.
# `thr` is a list of threshold vectors, one per indicator.
simulate_onefactor <- function(n, lambda, thr) {
  f <- rnorm(n)
  sapply(seq_along(thr), function(j)
    findInterval(lambda * f + sqrt(1 - lambda^2) * rnorm(n), thr[[j]]))
}

# Cumulative indicator matrix: one column per threshold, in label order.
cumulative_indicators <- function(Y, thr)
  do.call(cbind, lapply(seq_along(thr), function(j)
    outer(Y[, j], seq_along(thr[[j]]), "<"))) * 1

# --- writing ----------------------------------------------------------------
# write.csv() opens its file in TEXT mode, so on Windows it emits CRLF and on
# Linux it emits LF. The tier-1 fixtures were produced on Linux and are
# committed with LF; .gitattributes marks them -text so git never touches them
# again. Regenerating on Windows with plain write.csv() therefore reproduces
# every value exactly and still rewrites every byte, and fixtures/verify.R
# reports CHANGED on files with nothing wrong with them.
#
# TELL: `git diff --numstat` says every line changed, and
# `git diff --ignore-all-space` says nothing changed.
#
# Opening the connection in binary mode ("wb") makes the newline LF on every
# platform, which is what the header above promises. Verified: with this,
# 10_fixtures.R regenerated under R 4.6.0 on Windows reproduces the
# Linux-generated case A-D bytes exactly.
write_csv_lf <- function(x, path, ...) {
  con <- file(path, open = "wb")
  on.exit(close(con))
  write.csv(x, con, row.names = FALSE, quote = FALSE, ...)
  invisible(path)
}

# --- provenance -------------------------------------------------------------
# Fixtures are committed to the repo, so the manifest has to make it possible
# to tell whether a file changed because the generator changed, the seed
# changed, or a package version moved underneath it.

write_manifest <- function(path, files, seeds, extra = list()) {
  info <- list(
    generated_utc  = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version      = R.version.string,
    platform       = R.version$platform,
    seeds          = seeds,
    files          = as.list(tools::md5sum(files)),
    packages       = if (length(extra$packages)) extra$packages else "none (base R only)"
  )
  # hand-rolled JSON so the generator stays dependency-free
  esc <- function(x) gsub('"', '\\\\"', as.character(x))
  fmt <- function(x, ind = "  ") {
    if (is.list(x)) paste0("{\n", paste(sprintf('%s  "%s": %s', ind, names(x),
        vapply(x, function(v) fmt(v, paste0(ind, "  ")), character(1))),
        collapse = ",\n"), "\n", ind, "}")
    else if (length(x) > 1) paste0("[", paste(sprintf('"%s"', esc(x)), collapse = ", "), "]")
    else sprintf('"%s"', esc(x))
  }
  writeLines(fmt(info), path)
  invisible(info)
}
