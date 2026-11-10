# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan two-group fits of design M with
# the structural path FREE in group 1 and ABSENT in group 2 (g1only), at the four levels (target
# 0.3.0; spec 4.8 O2 as widened 2026-10-11 after mgb-absent-prereg passed)
#
# Executed 2026-10-11 (R 4.6.0, lavaan 0.7.2): every assertion held on the first run.
#
# Data: the committed fixtures/data/mg_beta.csv (99_'s data, read, never
# regenerated; 99_ remains its only writer). The generating raw beta is .4 in BOTH
# groups, so group 2 is deliberately misspecified, as in mgb-absent-prereg.
#
# Fits: 99_'s level syntax with "f2 ~ c(NA,0)*f1" (free in group 1, fixed 0 in group 2)
# (@sec-mg-beta mgb_syntax, beta "g1only": the syntax mgb-absent-prereg
# showed to agree exactly with our parameter table at every level). Options as 99_.
#
# Keys as 99_; one beta key, 3:1:2:1 (group 1 only).
#
# Outputs (fixtures/expected/), 99_'s layouts:
#   lavaan_mg_betag1_fits.csv, lavaan_mg_betag1_parameters.csv,
#   lavaan_mg_betag1_vcov.csv
#
# Assertions as 99_ (r1)-(r4) at 1e-9 relative (carried from 97_, chosen after
# measuring there, NOT re-measured here), (r5) lavaan's df = Result B4 for design
# M with nb = 1: 17, 20, 21, 25 (99_'s 16, 19, 20, 24 plus the one beta element
# no longer free), and exactly one free beta, in group 1.
#
# Run from the repository root: Rscript fixtures/generate/99b_lavaan_mg_betag1.R
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))
if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT  <- file.path("fixtures", "expected")
DOUT <- file.path("fixtures", "data")
MAN  <- file.path("fixtures", "MANIFEST_lavaan_mg_betag1.json")
written <- character(0)
TOL_REBUILD <- 1e-9

K   <- c(4L, 4L, 2L, 4L, 2L, 2L)
fac <- c(1L, 1L, 1L, 2L, 2L, 2L)
p   <- length(K)
y   <- paste0("y", seq_len(p))
NG  <- c(1200L, 2800L)
LEVELS <- c("configural", "thresholds", "metric", "scalar")
B4 <- c(configural = 17, thresholds = 20, metric = 21, scalar = 25)   # Result B4, design M, nb = 1 (g1only)

# ---- data: 99_'s committed file, read -----------------------------------------
dat <- read.csv(file.path(DOUT, "mg_beta.csv"))
stopifnot(identical(names(dat), c("g", y)), identical(as.integer(table(dat$g)), NG))
for (g in 1:2) for (j in seq_len(p))
  stopifnot(length(unique(dat[dat$g == g, j + 1L])) == K[j])     # no empty category

# ---- the level syntax (@sec-mg-beta's mgb_syntax, beta "free", psi "gen") ----
lev_syntax <- function(level) {
  poly <- vapply(1:2, function(f) any(K[fac == f] > 2L), NA)
  tied_tau <- level != "configural"; tied_lam <- level %in% c("metric", "scalar")
  out <- character()
  for (f in 1:2) {
    it <- which(fac == f)
    ld <- if (tied_lam) sprintf("c(l%d,l%d)*%s", it, it, y[it]) else sprintf("c(NA,NA)*%s", y[it])
    out <- c(out, paste0("f", f, " =~ ", paste(ld, collapse = " + ")))
  }
  if (tied_tau) for (j in seq_len(p)) {
    cc <- seq_len(K[j] - 1L)
    out <- c(out, sprintf("%s | %s", y[j], paste(sprintf("c(t%d_%d,t%d_%d)*t%d", j, cc, j, cc, cc), collapse = " + ")))
  }
  out <- c(out, "f2 ~ c(NA,0)*f1")
  psi_free <- switch(level, scalar = c(TRUE, TRUE), metric = poly, c(FALSE, FALSE))
  om_free  <- switch(level, configural = rep(FALSE, p), thresholds = K > 2L, rep(TRUE, p))
  nu_free  <- level %in% c("thresholds", "metric")
  paste(c(out,
          sprintf("%s ~*~ c(1,%s)*%s", y, ifelse(om_free, "NA", "1"), y),
          sprintf("%s ~ c(0,%s)*1", y, if (nu_free) "NA" else "0"),
          sprintf("f%d ~~ c(1,%s)*f%d", 1:2, ifelse(psi_free, "NA", "1"), 1:2),
          sprintf("f%d ~ c(0,%s)*1", 1:2, if (level == "scalar") "NA" else "0")), collapse = "\n")
}

datf <- dat; datf[y] <- lapply(datf[y], ordered)
fit_level <- function(level)
  lavaan::sem(lev_syntax(level), data = datf, group = "g", ordered = y, estimator = "WLSMV",
              parameterization = "delta", std.lv = TRUE, test = "scaled.shifted",
              auto.cov.y = FALSE, auto.cov.lv.x = FALSE,
              zero.add = c(0.5, 0.0), zero.keep.margins = TRUE, bad.marker.crit = 0)

need <- function(x, what) {
  u <- unlist(x)
  if (is.null(x) || length(x) == 0L || anyNA(u) || (is.numeric(u) && any(!is.finite(u))))
    stop("lavaan field missing or non-finite: ", what)
  x
}
rel <- function(a, b) abs(a - b) / max(abs(b), .Machine$double.xmin)

keys <- function(fit) {
  pt <- parTable(fit); fr <- pt[pt$free > 0, ]; fr <- fr[order(fr$free), ]
  fi <- function(v) as.integer(sub("f", "", v))
  vapply(seq_len(nrow(fr)), function(k) { x <- fr[k, ]; g <- as.integer(x$group)
    switch(x$op,
      "=~"  = sprintf("1:%d:%d:%d", g, match(x$rhs, y), fi(x$lhs)),
      "~~"  = if (x$lhs == x$rhs && x$lhs %in% c("f1", "f2")) sprintf("2:%d:%d:%d", g, fi(x$lhs), fi(x$lhs)) else NA_character_,
      "~"   = if (x$lhs == "f2" && x$rhs == "f1") sprintf("3:%d:2:1", g) else NA_character_,
      "|"   = sprintf("4:%d:%d:%d", g, match(x$lhs, y), as.integer(sub("t", "", x$rhs))),
      "~1"  = if (x$lhs %in% y) sprintf("5:%d:%d:1", g, match(x$lhs, y)) else sprintf("7:%d:%d:1", g, fi(x$lhs)),
      "~*~" = sprintf("6:%d:%d:1", g, match(x$lhs, y)),
      NA_character_) }, "")
}
slab <- function(fit) {          # lavaan's per-group statistic labels, ours: "y1|t1", "y1~~y2" (low index first)
  l <- rownames(lavInspect(fit, "delta")[[1]])
  vapply(l, function(s) if (grepl("~~", s, fixed = TRUE)) {
    ij <- sort(as.integer(sub("y", "", strsplit(s, "~~", fixed = TRUE)[[1]]))); sprintf("y%d~~y%d", ij[1], ij[2])
  } else s, "")
}

fit_rows <- list(); par_rows <- list(); v_rows <- list(); ss_rows <- NULL
for (lv in LEVELS) {
  fit <- fit_level(lv)
  stopifnot(lavInspect(fit, "converged"))
  ss  <- fit@test$scaled.shifted
  kk  <- keys(fit); stopifnot(!anyNA(kk), anyDuplicated(kk) == 0L)
  est <- need(coef(fit), "coef"); stopifnot(length(est) == length(kk))
  V   <- need(unname(vcov(fit)), "vcov"); stopifnot(nrow(V) == length(kk))
  N   <- sum(NG); pi <- NG / N
  Gm  <- need(lavTech(fit, "gamma"), "gamma"); Wv <- need(lavTech(fit, "WLS.V"), "WLS.V")
  wo  <- need(lavInspect(fit, "wls.obs"), "wls.obs"); we <- need(lavInspect(fit, "wls.est"), "wls.est")
  Fg  <- vapply(1:2, function(g) { r <- wo[[g]] - we[[g]]; sum(diag(Wv[[g]]) * r^2) }, 1)
  for (g in 1:2) {
    cat(sprintf("%-10s g%d  max rel |diag(WLS.V) - 1/diag(Gamma)| = %.2e; offdiag max %.2e\n", lv, g,
        max(abs(diag(Wv[[g]]) * diag(Gm[[g]]) - 1)), max(abs(Wv[[g]][row(Wv[[g]]) != col(Wv[[g]])]))))
    stopifnot(max(abs(Wv[[g]][row(Wv[[g]]) != col(Wv[[g]])])) == 0,                                   # (r2)
              max(abs(diag(Wv[[g]]) * diag(Gm[[g]]) - 1)) < TOL_REBUILD)
  }
  chisq <- need(unname(fitMeasures(fit, "chisq")), "chisq")
  stopifnot(rel(sum((NG - 1) * Fg), chisq) < TOL_REBUILD)                                              # (r1)
  D  <- do.call(rbind, lavTech(fit, "delta"))
  Bq <- lavaan:::lav_con_eq_basis(fit@Model); if (!is.null(Bq)) D <- D %*% Bq
  dg <- c(nrow(Gm[[1]]), nrow(Gm[[2]]))
  Om <- matrix(0, sum(dg), sum(dg)); W <- numeric(sum(dg)); o <- 0
  for (g in 1:2) { ix <- o + seq_len(dg[g]); Om[ix, ix] <- Gm[[g]] / pi[g]; W[ix] <- pi[g] / diag(Gm[[g]]); o <- o + dg[g] }
  WD <- W * D; U <- diag(W) - WD %*% solve(crossprod(D, WD), t(WD)); M <- U %*% Om
  tr <- c(sum(diag(M)), sum(M * t(M)))
  cat(sprintf("%-10s traces rebuilt %.10f %.10f ; lavaan %.10f %.10f\n", lv, tr[1], tr[2], ss$trace.UGamma, ss$trace.UGamma2))
  stopifnot(rel(tr[1], need(ss$trace.UGamma, "trace.UGamma")) < TOL_REBUILD,                           # (r3)
            rel(tr[2], need(ss$trace.UGamma2, "trace.UGamma2")) < TOL_REBUILD)
  stopifnot(rel(chisq / ss$scaling.factor + ss$shift.parameter, need(ss$stat, "ss$stat")) < TOL_REBUILD)   # (r4)
  stopifnot(unname(fitMeasures(fit, "df")) == B4[[lv]])                                                  # (r5)
  stopifnot(identical(kk[grepl("^3:", kk)], "3:1:2:1"))                                                 # one free beta, group 1
  pt <- parTable(fit); bb <- pt[pt$op == "~" & pt$lhs == "f2" & pt$rhs == "f1", ]
  stopifnot(nrow(bb) == 2L, bb$free[bb$group == 2] == 0L, bb$est[bb$group == 2] == 0, bb$free[bb$group == 1] > 0L)   # absent g2

  fit_rows[[lv]] <- data.frame(level = lv, N1 = NG[1], N2 = NG[2], df = B4[[lv]], npar_free = length(kk),
    chisq = chisq, F1 = Fg[1], F2 = Fg[2], T3 = ss$stat, scaling_factor = ss$scaling.factor,
    shift = ss$shift.parameter, trace_UG = ss$trace.UGamma, trace_UG2 = ss$trace.UGamma2)
  par_rows[[lv]] <- data.frame(level = lv, key = kk, est = unname(est))
  v_rows[[lv]]   <- data.frame(level = lv, key_i = rep(kk, length(kk)), key_j = rep(kk, each = length(kk)),
                               value = as.vector(V))
  sl <- slab(fit)
  sg <- do.call(rbind, lapply(1:2, function(g) data.frame(group = g, label = sl, value = unname(wo[[g]]))))
  if (is.null(ss_rows)) ss_rows <- sg else stopifnot(identical(ss_rows, sg))      # model-free: one per data set
}

fits <- do.call(rbind, fit_rows); pars <- do.call(rbind, par_rows); vc <- do.call(rbind, v_rows)
stopifnot(nrow(fits) == 4L, all(is.finite(as.matrix(fits[, -1]))), all(is.finite(pars$est)), all(is.finite(vc$value)))

f <- file.path(OUT, "lavaan_mg_betag1_fits.csv");       write_csv_lf(fits, f);    written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_betag1_parameters.csv"); write_csv_lf(pars, f);    written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_betag1_vcov.csv");       write_csv_lf(vc, f);      written <- c(written, f)

write_manifest(MAN, written, seeds = "none (data read from fixtures/data/mg_beta.csv, 99_'s)",
  extra = list(packages = c(paste("lavaan", as.character(packageVersion("lavaan"))),
    "WLSMV, delta, scaled.shifted, std.lv, auto.cov.y = FALSE, auto.cov.lv.x = FALSE, zero.add = c(0.5, 0), zero.keep.margins = TRUE, bad.marker.crit = 0 (pinned)",
    "explicit level syntax (@sec-mg-beta mgb_syntax, beta g1only (free in group 1, fixed 0 in group 2), psi gen), never group.equal",
    sprintf("rebuild tolerance %.0e relative, carried from 97_ (chosen after measuring there)", TOL_REBUILD),
    R.version.string)))
cat(sprintf("%d multigroup beta reference file(s) written; %s written.\n", length(written), MAN))
