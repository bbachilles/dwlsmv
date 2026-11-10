# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan's scaled-shifted chi-square
# DIFFERENCE test for nested pairs on Case E (DIFFTEST, target 0.3.0)
#
# Executed 2026-10-08 (R 4.6.0, lavaan 0.7.2): every assertion held.
#
# Reference: lavTestLRT(method = "satorra.2000", A.method = "delta",
# scaled.shifted = TRUE), i.e. the nesting Jacobian estimated as
# h-hat = ginv(Delta1(theta1-hat)) Delta0(theta0-hat), with Delta1, W, Gamma
# and P all taken from the H1 fit (spec Module 12 *Evaluation point* and
# *The statistic, frozen*; derivation @sec-difftest). Every per-pair quantity
# is REBUILT here from lavaan's own Delta, Gamma and WLS.V and asserted against
# lavTestLRT before it is written. The rebuild takes (Delta1, Delta0, Gamma,
# Td) and no N, the routine's inputs (E-003, 2026-10-08).
#
# Negative contrast, never a reference: A.method = "exact" (the parameter-
# selection H). lavaan's exact variant fails on xload2; that row carries
# lavaan's NA with a note beside the rebuilt selection-form value.
#
# Options identical to 20_/90_: WLSMV, delta, scaled.shifted,
# zero.add = c(0.5, 0), zero.keep.margins = TRUE, bad.marker.crit = 0
# (pinned, as 70_/90_/96_). Own manifest (MANIFEST_lavaan_difftest.json).
#
# Pairs (marker metric, case_E.csv):
#   beta0   H1 f2 ~ f1                     vs H0 f2 ~ 0*f1      d = 1 (H0 false)
#   xload2  H1 + y4 on f1 and y3 on f2     vs H0 correlated     d = 2 (H0 true)
#
# Outputs (fixtures/expected/):
#   lavaan_case_E_difftest_fits.csv      one row per fit:
#     pair, role, N, q, df, T_raw, T3, scaling_factor, shift, trace_UG, trace_UG2
#   lavaan_case_E_difftest.csv           one row per pair x variant
#     (reference | negative_contrast):
#     pair, variant, d, T0, T1, Td, tr_M1, tr_M1sq, a_d, b_d, stat, p,
#     lavaan_stat, lavaan_p, note
#   lavaan_case_E_difftest_matrices.csv  long: pair, role, matrix, row, col,
#     rowlab, collab, value. Delta per fit (columns in lavaan's free order),
#     hhat per pair (rows H1's free order, columns H0's), Gamma ONCE
#     (pair = all, role = none): it is model-free, asserted identical
#     across all four fits. W_D, V1 and P1 are not stored: derivable.
#
# Run from the repository root: Rscript fixtures/generate/97_lavaan_difftest.R
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))
if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan is required for tier-2 fixtures.")
if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS is required (ginv).")
library(lavaan)

OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_lavaan_difftest.json")
written <- character(0)

# Rebuilt statistic against lavTestLRT, relative. Chosen after measuring:
# 2026-10-08 gaps 5.8e-14 (beta0) and 2.7e-11 (xload2), lavaan 0.7.2.
TOL_REBUILD <- 1e-9

datE <- read.csv(file.path("fixtures", "data", "case_E.csv"))
datE[] <- lapply(datE, ordered)

fit_E <- function(model)
  lavaan::sem(model, data = datE, ordered = names(datE), estimator = "WLSMV",
              parameterization = "delta", test = "scaled.shifted",
              zero.add = c(0.5, 0.0), zero.keep.margins = TRUE,
              bad.marker.crit = 0)

BASE <- "f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5 + y6"
pairs <- list(
  beta0  = list(d = 1,
    h1 = paste(BASE, "f2 ~ f1", sep = "\n"),
    h0 = paste(BASE, "f2 ~ 0*f1", sep = "\n")),
  xload2 = list(d = 2,
    h1 = "f1 =~ y1 + y2 + y3 + y4\n f2 =~ y4 + y5 + y6 + y3\n f1 ~~ f2",
    h0 = paste(BASE, "f1 ~~ f2", sep = "\n")))

# data.frame() silently drops a NULL column: every lavaan field goes through this
need <- function(x, what) {
  u <- unlist(x)
  if (is.null(x) || length(x) == 0L || anyNA(u) || (is.numeric(u) && any(!is.finite(u))))
    stop("lavaan field missing or non-finite: ", what)
  x
}
# lavaan 0.7 renamed A.method and scaled.shifted; the old names travel through
# "..." with only a warning, so pin the new ones
stopifnot(all(c("method", "a_method", "scaled_shifted") %in% names(formals(lavaan::lavTestLRT))))

flab <- function(f) {
  pt <- parTable(f); pt <- pt[pt$free > 0, ]; pt <- pt[order(pt$free), ]
  paste0(pt$lhs, pt$op, pt$rhs)
}
trc <- function(A) sum(diag(A))

# M1 from H1's Delta, diag(Gamma)^-1, Gamma and a nesting matrix H; no N
m1_traces <- function(D1, Gm, H) {
  WD <- (1 / diag(Gm)) * D1
  P  <- crossprod(D1, WD); Pi <- solve(P)
  V  <- Pi %*% crossprod(WD, Gm %*% WD) %*% Pi
  M  <- (P - P %*% H %*% solve(t(H) %*% P %*% H) %*% t(H) %*% P) %*% V
  c(tr = trc(M), tr2 = sum(M * t(M)))
}
t3d <- function(tr, Td, d) {
  a <- sqrt(d / tr[["tr2"]]); b <- d - a * tr[["tr"]]
  c(a = a, b = b, stat = a * Td + b, p = pchisq(a * Td + b, d, lower.tail = FALSE))
}
lrt <- function(f0, f1, am) {
  o <- tryCatch(suppressWarnings(lavTestLRT(f0, f1, method = "satorra.2000",
         a_method = am, scaled_shifted = TRUE)), error = function(e) e)
  if (inherits(o, "error")) return(list(stat = NA_real_, p = NA_real_,
    note = gsub("[,\n\r]+", " ", conditionMessage(o))))
  stopifnot(all(c("Chisq diff", "Df diff", "Pr(>Chisq)") %in% names(o)))
  r <- which(!is.na(o[["Df diff"]])); stopifnot(length(r) == 1L)
  list(stat = need(o[["Chisq diff"]][r], "Chisq diff"),
       p = need(o[["Pr(>Chisq)"]][r], "Pr(>Chisq)"), df = o[["Df diff"]][r], note = "")
}
long <- function(M, pair, role, name, rl, cl)
  data.frame(pair = pair, role = role, matrix = name,
             row = rep(seq_len(nrow(M)), ncol(M)), col = rep(seq_len(ncol(M)), each = nrow(M)),
             rowlab = rep(rl, ncol(M)), collab = rep(cl, each = nrow(M)), value = as.vector(M))

fit_rows <- list(); diff_rows <- list(); mat_rows <- list(); G_ref <- NULL; slab <- NULL
for (tag in names(pairs)) {
  p  <- pairs[[tag]]
  f0 <- fit_E(p$h0); f1 <- fit_E(p$h1)
  stopifnot(lavInspect(f0, "converged"), lavInspect(f1, "converged"))
  l0 <- flab(f0); l1 <- flab(f1)
  stopifnot(all(l0 %in% l1), length(l1) - length(l0) == p$d,
            fitMeasures(f0, "df") - fitMeasures(f1, "df") == p$d)

  for (role in c("h0", "h1")) {
    f  <- if (role == "h0") f0 else f1
    ss <- f@test$scaled.shifted
    Gm <- need(lavTech(f, "gamma")[[1]], "gamma")
    if (is.null(G_ref)) {
      G_ref <- Gm
      slab  <- need(rownames(lavInspect(f, "delta")), "delta rownames")
    } else stopifnot(identical(Gm, G_ref))
    # premise of W_D = diag(Gamma)^-1: printed, not asserted (no tolerance
    # measured yet); the rebuild assertion below is what binds
    Wv <- need(lavTech(f, "WLS.V")[[1]], "WLS.V")
    cat(sprintf("%-7s %s  max rel |diag(WLS.V) - 1/diag(Gamma)| = %.2e; offdiag max %.2e\n",
        tag, role, max(abs(diag(Wv) * diag(Gm) - 1)), max(abs(Wv[row(Wv) != col(Wv)]))))
    D  <- need(lavTech(f, "delta")[[1]], "delta")
    stopifnot(ncol(D) == length(if (role == "h0") l0 else l1))
    mat_rows[[length(mat_rows) + 1]] <- long(D, tag, role, "Delta", slab, flab(f))
    fit_rows[[length(fit_rows) + 1]] <- data.frame(
      pair = tag, role = role, N = need(lavInspect(f, "nobs"), "nobs"), q = ncol(D),
      df = need(ss$df, "ss$df"), T_raw = need(f@test$standard$stat, "standard$stat"),
      T3 = need(ss$stat, "ss$stat"), scaling_factor = need(ss$scaling.factor, "scaling.factor"),
      shift = need(ss$shift.parameter, "shift.parameter"),
      trace_UG = need(ss$trace.UGamma, "trace.UGamma"),
      trace_UG2 = need(ss$trace.UGamma2, "trace.UGamma2"))
  }

  D1 <- lavTech(f1, "delta")[[1]]; D0 <- lavTech(f0, "delta")[[1]]
  T0 <- f0@test$standard$stat; T1 <- f1@test$standard$stat; Td <- T0 - T1
  hh <- MASS::ginv(D1) %*% D0
  Hs <- matrix(0, length(l1), length(l0)); Hs[cbind(match(l0, l1), seq_along(l0))] <- 1
  mat_rows[[length(mat_rows) + 1]] <- long(hh, tag, "pair", "hhat", l1, l0)

  for (v in c("reference", "negative_contrast")) {
    tr  <- m1_traces(D1, G_ref, if (v == "reference") hh else Hs)
    st  <- t3d(tr, Td, p$d)
    lv  <- lrt(f0, f1, if (v == "reference") "delta" else "exact")
    if (v == "reference") {
      stopifnot(!is.na(lv$stat), lv$df == p$d,
                abs(st[["stat"]] / lv$stat - 1) < TOL_REBUILD)
    } else if (!is.na(lv$stat)) {
      stopifnot(abs(st[["stat"]] / lv$stat - 1) < TOL_REBUILD)
    }
    note <- if (v == "reference") "lavTestLRT satorra.2000 A.method=delta" else
      if (is.na(lv$stat)) paste("selection H; rebuilt only: lavaan exact failed:", lv$note) else
        "selection H (A.method=exact); not the reference"
    diff_rows[[length(diff_rows) + 1]] <- data.frame(
      pair = tag, variant = v, d = p$d, T0 = T0, T1 = T1, Td = Td,
      tr_M1 = tr[["tr"]], tr_M1sq = tr[["tr2"]], a_d = st[["a"]], b_d = st[["b"]],
      stat = st[["stat"]], p = st[["p"]], lavaan_stat = lv$stat, lavaan_p = lv$p,
      note = note)
  }
}
mat_rows[[length(mat_rows) + 1]] <- long(G_ref, "all", "none", "Gamma", slab, slab)

fits  <- do.call(rbind, fit_rows)
diffs <- do.call(rbind, diff_rows)
mats  <- do.call(rbind, mat_rows)
stopifnot(nrow(fits) == 4L, nrow(diffs) == 4L, all(fits$T_raw > 0),
          all(is.finite(as.matrix(diffs[, c("Td", "tr_M1", "tr_M1sq", "a_d", "b_d", "stat", "p")]))),
          all(is.finite(mats$value)))

f <- file.path(OUT, "lavaan_case_E_difftest_fits.csv");     write_csv_lf(fits, f);  written <- c(written, f)
f <- file.path(OUT, "lavaan_case_E_difftest.csv");          write_csv_lf(diffs, f); written <- c(written, f)
f <- file.path(OUT, "lavaan_case_E_difftest_matrices.csv"); write_csv_lf(mats, f);  written <- c(written, f)

write_manifest(MAN, written,
  seeds = "none (deterministic given case_E.csv)",
  extra = list(packages = c(paste("lavaan", as.character(packageVersion("lavaan"))),
    paste("MASS", as.character(packageVersion("MASS"))),
    "WLSMV, delta, scaled.shifted, zero.add = c(0.5, 0), zero.keep.margins = TRUE, bad.marker.crit = 0 (pinned)",
    "reference: lavTestLRT satorra.2000, A.method = delta, scaled.shifted; rebuilt from Delta, Gamma, Td (no N)",
    "negative contrast: A.method = exact (selection H); NA where lavaan fails",
    sprintf("rebuild tolerance %.0e relative, chosen after measuring", TOL_REBUILD),
    R.version.string)))
cat(sprintf("%d difftest reference file(s) written; %s written.\n", length(written), MAN))
