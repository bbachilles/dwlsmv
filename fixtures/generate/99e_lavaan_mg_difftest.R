# ---------------------------------------------------------------------------
# Drafted 2026-10-13 by a read-only agent; design D and the repo-relative output
# paths added by the main session. Executed by the main session 2026-10-13 from
# the repository root (R 4.6.0, lavaan 0.7.2): every assertion held, outputs
# identical to two scratch runs, fixtures/verify.R 88 files 0 problems.
#
# dwlsmv fixture generation -- tier 2: lavaan's scaled-shifted chi-square
# DIFFERENCE test between GENERATED two-group invariance levels (the
# generated-chain exemption gate; spec Module 12 checks 8-10, E-094,
# 2026-10-13 rows; derivation @sec-difftest, @sec-invariance-id,
# @sec-multigroup-scaling MG7, @sec-inv-nesting)
#
# Reference: lavTestLRT(method = "satorra.2000", a_method = "delta",
# scaled_shifted = TRUE) (lavaan 0.7.2 argument names, as 97_). Every pair's
# statistic is REBUILT from lavaan's own objects as mg_dt() does in the
# derivation (chunk mg-difftest): stacked, equality-constrained Delta
# (lavTech "delta" times lav_con_eq_basis), h-hat = ginv(Dst1) Dst0 (the
# unweighted Moore-Penrose on stacked rows: E-019),
# pooled P = sum_g pi_g D_g' W_g D_g with W_g = diag(Gamma_g)^-1, V1 the pooled
# ACOV (MG6), M1 = (P - P h (h'P h)^-1 h'P) V1, a_d = sqrt(d / tr(M1^2)),
# b_d = d - a_d tr(M1), T_d3 = a_d Td + b_d.
#
# Designs (derivation @sec-invariance-id inv_designs; one factor, p = 5, G = 2):
#   E  K = 4,4,4,2,2  data: the committed fixtures/data/mg_E.csv (98_, seed
#      20261010; N1 = 1200, N2 = 2800). Steps d = 3, 2, 4 (df 10, 13, 15, 19).
#   B  K = 3,3,3,3,3  data: simulated here, seed SEED_B, same DGP and group
#      sizes as 98_ (inv_sim's lambda, tau, group-2 kappa .3, phi 1.44,
#      residuals .4-.8). Result I4: df 10, 10, 14, 18, so the THRESHOLDS step
#      has d = 0 (the unverified d = 0 pattern of spec Module 12 / G8 note
#      "(i)"), and the metric and scalar steps have d = 4 with every item
#      polytomous (omega free for all items from the thresholds level on).
#   Why B: the derivation's designs are all one-factor (no two-factor design
#      exists in A-E, so none is invented here). B is the K = 3-only design;
#      it adds a d = 0 step whose lavaan behaviour must be reported, and a
#      metric step in which every item is polytomous, which E lacks.
#   D  K = 4,2,2,2,2  data: simulated here, seed SEED_D, same DGP and group
#      sizes. Result I4: df 10, 11, 11, 15, so the METRIC step has d = 0, the
#      verified Result I2 case (one polytomous item per factor). Added by the
#      main session, 2026-10-13, at the user's request (zero-df-step coverage).
#
# Levels: the explicit level syntax of 98_ / inv_syntax (never group.equal),
# derivation names configural, thresholds, loadings, intercepts; the files
# carry the surface names (spec G8 L2 (a)): loadings = metric,
# intercepts = scalar.
#
# Pairs per design: configural-thresholds, thresholds-metric, metric-scalar
# (adjacent) and configural-scalar (non-adjacent). H1 = less restricted.
#
# T conventions: lavaan's T = sum_g (N_g - 1) F_g (candidate C); dwlsmv's is
# T = (N - G) F, F = sum_g pi_g F_g (candidate B, spec 2026-10-09). Both are
# written. T_d3 is written twice: from lavaan's Td (compared with
# lavTestLRT) and from (N - G)(F0 - F1) (what dwlsmv should produce; M1 is
# the same, it does not depend on N).
#
# d = 0: lavTestLRT is called anyway and its result or error is written
# verbatim (note column, commas stripped); the rebuild records tr(M1),
# tr(M1^2) and Td; a_d, b_d, T_d3, p are NA (sqrt(0/0)); no assertion is made
# on the statistic for that row.
#
# Outputs (fixtures/expected/):
#   lavaan_mg_difftest_fits.csv        design, level, N1, N2, df, q (constrained
#     free count), npar_lav (labelled), T_lav, T_NG, F1, F2, T3, scaling_factor, shift
#   lavaan_mg_difftest_parameters.csv  design, level, key, est (98_'s keys)
#   lavaan_mg_difftest.csv             design, h1, h0, d, T0, T1, Td,
#     T0_NG, T1_NG, Td_NG, tr_M1, tr_M1sq, a_d, b_d, T_d3, p, T_d3_NG, p_NG,
#     lavaan_stat, lavaan_df, lavaan_p, note
#   fixtures/data/mg_B.csv, mg_D.csv   g, y1-y5
# No Delta or Gamma dump: dwlsmv's side refits and builds its own; lavaan's
# matrices are reachable by rerunning this file.
#
# Run from the repository root:
#   Rscript fixtures/generate/99e_lavaan_mg_difftest.R
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

t_start <- Sys.time()
source(file.path("fixtures", "generate", "00_common.R"))
if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan is required for tier-2 fixtures.")
if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS is required (ginv).")
library(lavaan)

OUT  <- file.path("fixtures", "expected")
DOUT <- file.path("fixtures", "data")
MAN  <- file.path("fixtures", "MANIFEST_lavaan_mg_difftest.json")
written <- character(0)

TOL_REBUILD <- 1e-9  # per-fit rebuilds (T, against 98_'s committed values); carried from 97_/98_
TOL_D       <- 1e-8  # T_d3 against lavTestLRT; mg_TOLD of @sec-multigroup-scaling, set BEFORE measuring there,
                     # carried here unchanged, NOT re-chosen after this run
SEED_B <- 20261013L
SEED_D <- 20261015L   # design D data (added by the main session, 2026-10-13, user)

stopifnot(all(c("method", "a_method", "scaled_shifted") %in% names(formals(lavaan::lavTestLRT))))

p   <- 5L; y <- paste0("y", seq_len(p))
lam <- c(.60, .65, .70, .75, .80)
NG  <- c(1200L, 2800L)
LEVELS  <- c("configural", "thresholds", "loadings", "intercepts")
SURFACE <- c(configural = "configural", thresholds = "thresholds", loadings = "metric", intercepts = "scalar")
DESIGNS <- list(E = c(4L, 4L, 4L, 2L, 2L), B = rep(3L, 5), D = c(4L, 2L, 2L, 2L, 2L))
I4 <- list(E = c(10, 13, 15, 19), B = c(10, 10, 14, 18), D = c(10, 11, 11, 15))      # Result I4 (inv_df), asserted below
PAIRS <- list(c("configural", "thresholds"), c("thresholds", "loadings"),
              c("loadings", "intercepts"), c("configural", "intercepts"))

# ---- data -------------------------------------------------------------------
sim_K <- function(K, seed) {                                    # 98_'s DGP, any K (B: all 3; D: 4,2,2,2,2)
  tau <- lapply(seq_len(p), function(j) qnorm(seq_len(K[j] - 1L) / K[j]) + .1 * j)
  set.seed(seed)
  one <- function(n, kap, phi, th) {
    eta <- rnorm(n, kap, sqrt(phi))
    sapply(seq_len(p), function(j) findInterval(lam[j] * eta + rnorm(n, 0, sqrt(th[j])), tau[[j]]))
  }
  d <- data.frame(g = rep(1:2, NG),
                  rbind(one(NG[1], 0, 1, 1 - lam^2), one(NG[2], .3, 1.44, seq(.4, .8, length.out = p))))
  names(d)[-1] <- y
  d
}
DATA <- list(E = read.csv(file.path("fixtures", "data", "mg_E.csv")), B = sim_K(DESIGNS$B, SEED_B), D = sim_K(DESIGNS$D, SEED_D))
for (nm in names(DATA)) {
  stopifnot(identical(names(DATA[[nm]]), c("g", y)), identical(as.integer(table(DATA[[nm]]$g)), NG))
  for (g in 1:2) for (j in seq_len(p))
    stopifnot(length(unique(DATA[[nm]][DATA[[nm]]$g == g, j + 1L])) == DESIGNS[[nm]][j])   # no empty category
}

# ---- the level syntax (98_'s lev_spec / lev_syntax with K as an argument) ------
lev_spec <- function(K, level) {
  poly <- K > 2L
  list(tied_tau = level != "configural",
       tied_lam = level %in% c("loadings", "intercepts"),
       om  = switch(level, configural = rep("fix", p), thresholds = ifelse(poly, "free", "fix"), rep("free", p)),
       nu  = if (level %in% c("configural", "intercepts")) rep("fix", p) else rep("free", p),
       phi = level == "intercepts" || (level == "loadings" && any(poly)),
       kap = level == "intercepts")
}
lev_syntax <- function(K, level) {
  sp <- lev_spec(K, level)
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

WARN <- character(0)
catch_w <- function(expr, where)
  withCallingHandlers(expr, warning = function(w) {
    WARN <<- c(WARN, sprintf("[%s] %s", where, conditionMessage(w))); invokeRestart("muffleWarning") })

fit_level <- function(dat, K, level) {
  datf <- dat; datf[y] <- lapply(datf[y], ordered)
  catch_w(lavaan::cfa(lev_syntax(K, level), data = datf, group = "g", ordered = y, estimator = "WLSMV",
                      parameterization = "delta", std.lv = TRUE, test = "scaled.shifted",
                      zero.add = c(0.5, 0.0), zero.keep.margins = TRUE, bad.marker.crit = 0),
          paste("fit", level))
}

need <- function(x, what) {
  u <- unlist(x)
  if (is.null(x) || length(x) == 0L || anyNA(u) || (is.numeric(u) && any(!is.finite(u))))
    stop("lavaan field missing or non-finite: ", what)
  x
}
rel <- function(a, b) max(abs(a - b)) / max(abs(b), .Machine$double.xmin)
keys <- function(fit) {                                         # 98_'s keys, unchanged
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

# pooled pieces of one fit (derivation mg_parts / mg_pool, trimmed to what MG7 needs)
parts <- function(fit) {
  n_g <- as.numeric(unlist(lavInspect(fit, "nobs"))); N <- sum(n_g); pg <- n_g / N
  Gam <- lapply(need(lavTech(fit, "gamma"), "gamma"), as.matrix)
  wv  <- lapply(need(lavTech(fit, "wls.v"), "wls.v"), function(x) if (is.matrix(x)) diag(x) else as.numeric(x))
  for (g in seq_along(Gam)) stopifnot(rel(wv[[g]], 1 / diag(Gam[[g]])) < 1e-12)   # W_g = diag(Gamma_g)^-1
  Bq  <- lavaan:::lav_con_eq_basis(fit@Model)
  Dg  <- lapply(need(lavTech(fit, "delta"), "delta"), function(D) if (is.null(Bq)) D else D %*% Bq)
  obs <- lapply(lavTech(fit, "wls.obs"), as.numeric)
  est <- lapply(lavaan:::lav_model_wls_est(fit@Model), as.numeric)
  Fg  <- mapply(function(x, o, e) sum(x * (o - e)^2), wv, obs, est)
  P   <- Reduce(`+`, Map(function(D, x, w) w * crossprod(D, x * D), Dg, wv, pg))
  Bm  <- Reduce(`+`, Map(function(D, x, Gm, w) { WD <- x * D; w * crossprod(WD, Gm %*% WD) }, Dg, wv, Gam, pg))
  Pi  <- solve(P)
  list(N = N, n_g = n_g, pg = pg, G = length(n_g), Gam = Gam, Dst = do.call(rbind, Dg), q = ncol(Dg[[1]]),
       Fg = Fg, F = sum(pg * Fg), P = P, V = Pi %*% Bm %*% Pi)
}

fits <- list(); fit_rows <- list(); par_rows <- list(); diff_rows <- list(); ASSERT <- character(0)
note_assert <- function(txt) { ASSERT <<- c(ASSERT, txt); cat(txt, "\n") }

for (nm in names(DESIGNS)) {
  K <- DESIGNS[[nm]]; fits[[nm]] <- list(); pp <- list()
  for (k in seq_along(LEVELS)) {
    lv <- LEVELS[k]
    fit <- fit_level(DATA[[nm]], K, lv)
    stopifnot(lavInspect(fit, "converged"))
    fits[[nm]][[lv]] <- fit; pp[[lv]] <- o <- parts(fit)
    ss <- fit@test$scaled.shifted
    Tl <- need(fit@test$standard$stat, "standard$stat")
    df <- unname(fitMeasures(fit, "df"))
    stopifnot(df == I4[[nm]][k], nrow(o$Dst) - o$q == df)                     # Result I4; D_total - q = df
    r1 <- rel(sum((o$n_g - 1) * o$Fg), Tl)                                    # lavaan's T is candidate C
    note_assert(sprintf("%s %-10s df %2d q %2d  rel|sum(N_g-1)F_g - T_lav| = %.2e (tol %.0e)", nm, lv, df, o$q, r1, TOL_REBUILD))
    stopifnot(r1 < TOL_REBUILD)
    kk <- keys(fit); stopifnot(!anyNA(kk))
    est <- need(coef(fit), "coef"); stopifnot(length(est) == length(kk))
    if (nm == "E") {                                                          # same data, same fits as 98_
      ref <- read.csv(file.path("fixtures", "expected", "lavaan_mg_E_fits.csv"))
      rp  <- read.csv(file.path("fixtures", "expected", "lavaan_mg_E_parameters.csv"), colClasses = c("character", "character", "numeric"))
      rr <- ref[ref$level == lv, ]; rq <- rp[rp$level == lv, ]
      r98 <- c(T = rel(Tl, rr$chisq), T3 = rel(ss$stat, rr$T3))
      stopifnot(identical(rq$key, kk)); r98p <- rel(unname(est), rq$est)
      note_assert(sprintf("E %-10s vs committed 98_ fixture: rel T %.2e, T3 %.2e, est %.2e (tol %.0e)", lv, r98[1], r98[2], r98p, TOL_REBUILD))
      stopifnot(r98 < TOL_REBUILD, r98p < TOL_REBUILD)
    }
    fit_rows[[length(fit_rows) + 1]] <- data.frame(design = nm, level = SURFACE[[lv]], N1 = NG[1], N2 = NG[2],
      df = df, q = o$q, npar_lav = length(kk), T_lav = Tl, T_NG = (o$N - o$G) * o$F, F1 = o$Fg[1], F2 = o$Fg[2],
      T3 = need(ss$stat, "ss$stat"), scaling_factor = need(ss$scaling.factor, "scaling.factor"),
      shift = need(ss$shift.parameter, "shift.parameter"))
    par_rows[[length(par_rows) + 1]] <- data.frame(design = nm, level = SURFACE[[lv]], key = kk, est = unname(est))
  }

  for (pr in PAIRS) {
    l1 <- pr[1]; l0 <- pr[2]; p1 <- pp[[l1]]; p0 <- pp[[l0]]
    f1 <- fits[[nm]][[l1]]; f0 <- fits[[nm]][[l0]]
    d  <- p1$q - p0$q
    stopifnot(d >= 0, d == unname(fitMeasures(f0, "df") - fitMeasures(f1, "df")))
    stopifnot(rel(p1$Gam[[1]], p0$Gam[[1]]) == 0, rel(p1$Gam[[2]], p0$Gam[[2]]) == 0)   # Gamma model-free
    hh <- MASS::ginv(p1$Dst) %*% p0$Dst
    M1 <- (p1$P - p1$P %*% hh %*% solve(t(hh) %*% p1$P %*% hh) %*% t(hh) %*% p1$P) %*% p1$V
    trM <- sum(diag(M1)); trM2 <- sum(M1 * t(M1))
    T0 <- f0@test$standard$stat; T1 <- f1@test$standard$stat; Td <- T0 - T1
    c_NG <- p1$N - p1$G; Td_NG <- c_NG * (p0$F - p1$F)
    if (d > 0) {
      a <- sqrt(d / trM2); b <- d - a * trM
      Td3 <- a * Td + b; pv <- pchisq(Td3, d, lower.tail = FALSE)
      Td3n <- a * Td_NG + b; pvn <- pchisq(Td3n, d, lower.tail = FALSE)
    } else { a <- b <- Td3 <- pv <- Td3n <- pvn <- NA_real_ }
    lv_res <- tryCatch(catch_w(lavTestLRT(f0, f1, method = "satorra.2000", a_method = "delta", scaled_shifted = TRUE),
                               sprintf("lavTestLRT %s %s-%s", nm, l1, l0)), error = function(e) e)
    if (inherits(lv_res, "error")) {
      ls <- NA_real_; ldf <- NA_real_; lp <- NA_real_
      note <- paste("lavTestLRT error:", gsub("[,\n\r]+", " ", conditionMessage(lv_res)))
    } else {
      r <- which(!is.na(lv_res[["Df diff"]]))
      if (length(r) == 1L) {
        ls <- lv_res[["Chisq diff"]][r]; ldf <- lv_res[["Df diff"]][r]; lp <- lv_res[["Pr(>Chisq)"]][r]
      } else { ls <- NA_real_; ldf <- NA_real_; lp <- NA_real_ }
      note <- paste("lavTestLRT returned; row 2:", gsub("[,\n\r]+", " ",
                    paste(names(lv_res), format(unlist(lv_res[nrow(lv_res), ]), digits = 10), sep = "=", collapse = " ")))
    }
    tag <- sprintf("%s %s-%s d=%d", nm, SURFACE[[l1]], SURFACE[[l0]], d)
    if (d > 0) {
      rd <- abs(Td3 / ls - 1)
      note_assert(sprintf("%s  T_d3 rebuilt %.12g  lavaan %.12g  rel %.2e (tol %.0e); df lavaan %s",
                          tag, Td3, ls, rd, TOL_D, format(ldf)))
      stopifnot(!is.na(ls), ldf == d, rd < TOL_D)
      cat(sprintf("   Td_lav %.10g  Td_NG %.10g  T_d3_NG %.10g  p %.6g / %.6g\n", Td, Td_NG, Td3n, pv, pvn))
    } else {
      cat(sprintf("%s  d = 0: Td_lav %.6e Td_NG %.6e tr(M1) %.3e tr(M1^2) %.3e; lavaan stat %s df %s p %s\n   %s\n",
                  tag, Td, Td_NG, trM, trM2, format(ls), format(ldf), format(lp), note))
    }
    diff_rows[[length(diff_rows) + 1]] <- data.frame(design = nm, h1 = SURFACE[[l1]], h0 = SURFACE[[l0]], d = d,
      T0 = T0, T1 = T1, Td = Td, T0_NG = c_NG * p0$F, T1_NG = c_NG * p1$F, Td_NG = Td_NG,
      tr_M1 = trM, tr_M1sq = trM2, a_d = a, b_d = b, T_d3 = Td3, p = pv, T_d3_NG = Td3n, p_NG = pvn,
      lavaan_stat = ls, lavaan_df = ldf, lavaan_p = lp, note = note)
  }
}

fitsdf <- do.call(rbind, fit_rows); pars <- do.call(rbind, par_rows); diffs <- do.call(rbind, diff_rows)
stopifnot(nrow(fitsdf) == length(LEVELS) * length(DESIGNS), nrow(diffs) == length(PAIRS) * length(DESIGNS), all(is.finite(as.matrix(fitsdf[, -(1:2)]))), all(is.finite(pars$est)))

dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(DOUT, recursive = TRUE, showWarnings = FALSE)
f <- file.path(DOUT, "mg_B.csv");                          write_csv_lf(DATA$B, f);  written <- c(written, f)
f <- file.path(DOUT, "mg_D.csv");                          write_csv_lf(DATA$D, f);  written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_difftest_fits.csv");        write_csv_lf(fitsdf, f);  written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_difftest_parameters.csv");  write_csv_lf(pars, f);    written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_difftest.csv");             write_csv_lf(diffs, f);   written <- c(written, f)

write_manifest(MAN, written, seeds = sprintf("%d (design B data), %d (design D data); design E reads mg_E.csv (98_, seed 20261010)", SEED_B, SEED_D),
  extra = list(packages = c(paste("lavaan", as.character(packageVersion("lavaan"))),
    paste("MASS", as.character(packageVersion("MASS"))),
    "WLSMV, delta, scaled.shifted, std.lv, zero.add = c(0.5, 0), zero.keep.margins = TRUE, bad.marker.crit = 0 (pinned)",
    "explicit level syntax (98_ lev_syntax), never group.equal; files use surface level names",
    "reference: lavTestLRT satorra.2000, a_method = delta, scaled_shifted; rebuilt by ginv(Dst1) Dst0 with pooled P and V1 (MG7)",
    sprintf("T_d3 tolerance %.0e relative (mg_TOLD, set before measuring in the derivation); rebuild %.0e carried from 97_/98_", TOL_D, TOL_REBUILD),
    R.version.string)))

cat("\nWarnings captured (verbatim):\n"); if (length(WARN)) cat(paste0("  ", WARN), sep = "\n") else cat("  none\n")
cat(sprintf("%d file(s) written; %s written. Elapsed %.1f s\n", length(written), MAN,
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
