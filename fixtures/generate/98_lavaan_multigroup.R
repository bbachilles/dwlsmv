# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan two-group fits of design E at
# the four levels (pooled fit and inference, target 0.3.0; spec 4.8 O14,
# derivation @sec-multigroup-scaling MG1-MG6 and @sec-invariance-id)
#
# Executed 2026-10-10 (R 4.6.0, lavaan 0.7.2): every assertion held on the first run.
#
# Data: design E (K = 4, 4, 4, 2, 2; one factor), two groups of UNEQUAL size
# (N1 = 1200, N2 = 2800; the user's check A, 2026-10-10: equal pi_g would hide
# a wrong weighting and makes T candidates B and C coincide). Group 1 is the
# delta reference (total latent-response variance 1); group 2 has factor mean
# .3, variance 1.44 and residual variances .4-.8, same lambda and tau as
# @sec-invariance-id's inv_sim. Written to fixtures/data/mg_E.csv (g, y1-y5,
# categories 0..K-1).
#
# Fits: the explicit level syntax of @sec-invariance-id (inv_syntax; nothing
# relies on lavaan's automatic freeing, and NOT group.equal: @sec-mg-beta
# H-L5), WLSMV, delta, scaled.shifted, std.lv, with the options pinned as in
# 20_/97_.
#
# Every free lavaan row is keyed as the parameter table keys it (spec 4.8 G6):
# key = "ptype:group:row:col", ptype 1 lambda (item, factor), 2 psi (1, 1),
# 4 tau (item, k), 5 nu (item, 1), 6 omega (item, 1), 7 alpha (1, 1). A tied
# element appears once per member group, as lavaan has it.
#
# Outputs (fixtures/expected/):
#   lavaan_mg_E_fits.csv    one row per level: level, N1, N2, df, npar_free,
#     chisq (lavaan's unscaled; candidate C), F1, F2 (each group's
#     r_g' W_g r_g at lavaan's estimate), T3 (lavaan's scaled.shifted stat),
#     scaling_factor, shift, trace_UG, trace_UG2
#   lavaan_mg_E_parameters.csv   level, key, est
#   lavaan_mg_E_vcov.csv         level, key_i, key_j, value (lavaan's vcov(),
#     over its free rows: e(V) = ACOV / (N - G), MG6)
#   lavaan_mg_E_sampstat.csv     group, label ("y1|t1", "y1~~y2"), value
#     (lavaan's per-group sample statistics: Stage 1 per group)
#
# Assertions before anything is written (each rebuilds a lavaan number from
# lavaan's own objects; tolerance 1e-9 relative, carried from 97_'s rebuild
# tolerance, which was chosen after measuring there, NOT re-measured here):
#   (r1) chisq = sum_g (N_g - 1) F_g (lavaan's chi-square is candidate C, MG3)
#   (r2) diag(WLS.V_g) = 1 / diag(Gamma_g) and WLS.V_g is diagonal (printed)
#   (r3) the pooled traces from (Delta, Gamma_g, pi_g) by explicit U and M
#        (MG3/MG4) equal lavaan's trace.UGamma, trace.UGamma2
#   (r4) T3 = chisq / scaling_factor + shift (lavaan's own composition)
#   (r5) lavaan's df = Result I4 for design E: 10, 13, 15, 19
#
# Run from the repository root: Rscript fixtures/generate/98_lavaan_multigroup.R
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))
if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT  <- file.path("fixtures", "expected")
DOUT <- file.path("fixtures", "data")
MAN  <- file.path("fixtures", "MANIFEST_lavaan_multigroup.json")
written <- character(0)
TOL_REBUILD <- 1e-9
SEED <- 20261010L

K   <- c(4L, 4L, 4L, 2L, 2L)
p   <- length(K)
y   <- paste0("y", seq_len(p))
lam <- c(.60, .65, .70, .75, .80)
tau <- lapply(seq_len(p), function(j) qnorm(seq_len(K[j] - 1L) / K[j]) + .1 * j)
NG  <- c(1200L, 2800L)
LEVELS <- c("configural", "thresholds", "loadings", "intercepts")

# ---- data -------------------------------------------------------------------
set.seed(SEED)
one <- function(n, kap, phi, th) {
  eta <- rnorm(n, kap, sqrt(phi))
  sapply(seq_len(p), function(j) findInterval(lam[j] * eta + rnorm(n, 0, sqrt(th[j])), tau[[j]]))
}
dat <- data.frame(g = rep(1:2, NG),
                  rbind(one(NG[1], 0, 1, 1 - lam^2), one(NG[2], .3, 1.44, seq(.4, .8, length.out = p))))
names(dat)[-1] <- y
for (g in 1:2) for (j in seq_len(p))
  stopifnot(length(unique(dat[dat$g == g, j + 1L])) == K[j])     # no empty category

# ---- the level table (@sec-invariance-id's inv_spec / inv_syntax) ------------
lev_spec <- function(level) {
  poly <- K > 2L
  list(tied_tau = level != "configural",
       tied_lam = level %in% c("loadings", "intercepts"),
       om  = switch(level, configural = rep("fix", p), thresholds = ifelse(poly, "free", "fix"), rep("free", p)),
       nu  = if (level %in% c("configural", "intercepts")) rep("fix", p) else rep("free", p),
       phi = level == "intercepts" || (level == "loadings" && any(poly)),
       kap = level == "intercepts")
}
lev_syntax <- function(level) {
  sp <- lev_spec(level)
  ld <- if (sp$tied_lam) sprintf("c(l%d,l%d)*%s", seq_len(p), seq_len(p), y) else sprintf("c(NA,NA)*%s", y)
  out <- paste("f =~", paste(ld, collapse = " + "))
  if (sp$tied_tau) for (j in seq_len(p)) {
    cc <- seq_len(K[j] - 1L)
    out <- c(out, sprintf("%s | %s", y[j], paste(sprintf("c(t%d_%d,t%d_%d)*t%d", j, cc, j, cc, cc), collapse = " + ")))
  }
  paste(c(out,
          sprintf("%s ~*~ c(1,%s)*%s", y, ifelse(sp$om == "free", "NA", "1"), y),
          sprintf("%s ~ c(0,%s)*1",    y, ifelse(sp$nu == "free", "NA", "0")),
          sprintf("f ~~ c(1,%s)*f", if (sp$phi) "NA" else "1"),
          sprintf("f ~ c(0,%s)*1",  if (sp$kap) "NA" else "0")), collapse = "\n")
}
I4 <- c(configural = 10, thresholds = 13, loadings = 15, intercepts = 19)   # Result I4, design E

datf <- dat; datf[y] <- lapply(datf[y], ordered)
fit_level <- function(level)
  lavaan::cfa(lev_syntax(level), data = datf, group = "g", ordered = y, estimator = "WLSMV",
              parameterization = "delta", std.lv = TRUE, test = "scaled.shifted",
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
  vapply(seq_len(nrow(fr)), function(k) { x <- fr[k, ]; g <- as.integer(x$group)
    switch(x$op,
      "=~"  = sprintf("1:%d:%d:1", g, match(x$rhs, y)),
      "~~"  = if (x$lhs == "f" && x$rhs == "f") sprintf("2:%d:1:1", g) else NA_character_,
      "|"   = sprintf("4:%d:%d:%d", g, match(x$lhs, y), as.integer(sub("t", "", x$rhs))),
      "~1"  = if (x$lhs %in% y) sprintf("5:%d:%d:1", g, match(x$lhs, y)) else sprintf("7:%d:1:1", g),
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
  kk  <- keys(fit); stopifnot(!anyNA(kk))
  est <- need(coef(fit), "coef"); stopifnot(length(est) == length(kk))
  V   <- need(unname(vcov(fit)), "vcov"); stopifnot(nrow(V) == length(kk))
  N   <- sum(NG); pi <- NG / N
  Gm  <- need(lavTech(fit, "gamma"), "gamma"); Wv <- need(lavTech(fit, "WLS.V"), "WLS.V")
  ob  <- need(lavTech(fit, "sampstat"), "sampstat")
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
  # (r3) pooled traces by explicit U and M: Delta stacked, ties through lavaan's equality basis
  D  <- do.call(rbind, lavTech(fit, "delta"))
  Bq <- lavaan:::lav_con_eq_basis(fit@Model); if (!is.null(Bq)) D <- D %*% Bq
  dg <- c(nrow(Gm[[1]]), nrow(Gm[[2]]))
  Om <- matrix(0, sum(dg), sum(dg)); W <- numeric(sum(dg)); o <- 0
  for (g in 1:2) { ix <- o + seq_len(dg[g]); Om[ix, ix] <- Gm[[g]] / pi[g]; W[ix] <- pi[g] / diag(Gm[[g]]); o <- o + dg[g] }
  WD <- W * D; U <- diag(W) - WD %*% solve(crossprod(D, WD), t(WD)); M <- U %*% Om
  tr <- c(sum(diag(M)), sum(M * t(M)))
  cat(sprintf("%-10s traces rebuilt %.10f %.10f ; lavaan %.10f %.10f\n", lv, tr[1], tr[2], ss$trace.UGamma, ss$trace.UGamma2))
  stopifnot(rel(tr[1], need(ss$trace.UGamma, "trace.UGamma")) < TOL_REBUILD,
            rel(tr[2], need(ss$trace.UGamma2, "trace.UGamma2")) < TOL_REBUILD)
  stopifnot(rel(chisq / ss$scaling.factor + ss$shift.parameter, need(ss$stat, "ss$stat")) < TOL_REBUILD)   # (r4)
  stopifnot(unname(fitMeasures(fit, "df")) == I4[[lv]])                                                  # (r5)

  fit_rows[[lv]] <- data.frame(level = lv, N1 = NG[1], N2 = NG[2], df = I4[[lv]], npar_free = length(kk),
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

f <- file.path(DOUT, "mg_E.csv");                    write_csv_lf(dat, f);     written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_E_fits.csv");         write_csv_lf(fits, f);    written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_E_parameters.csv");   write_csv_lf(pars, f);    written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_E_vcov.csv");         write_csv_lf(vc, f);      written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_E_sampstat.csv");     write_csv_lf(ss_rows, f); written <- c(written, f)

write_manifest(MAN, written, seeds = sprintf("%d (data)", SEED),
  extra = list(packages = c(paste("lavaan", as.character(packageVersion("lavaan"))),
    "WLSMV, delta, scaled.shifted, std.lv, zero.add = c(0.5, 0), zero.keep.margins = TRUE, bad.marker.crit = 0 (pinned)",
    "explicit level syntax (inv_syntax), never group.equal",
    sprintf("rebuild tolerance %.0e relative, carried from 97_ (chosen after measuring there)", TOL_REBUILD),
    R.version.string)))
cat(sprintf("%d multigroup reference file(s) written; %s written.\n", length(written), MAN))
