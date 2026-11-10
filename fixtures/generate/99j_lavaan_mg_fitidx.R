# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan 0.7.2's scaled fit indices at
# G = 2 across generated invariance levels, and their differences
# (study E-018; spec Module 12 *Fit indices*,
# F1-F6; derivation @sec-fit-indices FI1-FI3)
#
# Written 2026-10-31 by a side session (G > 1 fit-index validation).
#
# Designs (one factor on y1-y5, two groups, N 1200/2800; levels and syntax
# exactly 99e_'s lev_spec / lev_syntax, never group.equal):
#   E  committed fixtures/data/mg_E.csv (98_)   well specified: clamps
#   D  committed fixtures/data/mg_D.csv (99e_)  K = 4,2,2,2,2: metric d = 0
#   P  simulated here, SEED_P: TWO-factor DGP (y1-y3, y4-y5, corr R_F)
#      fitted with one factor, measurement invariant across groups
#   N  simulated here, SEED_N: as P with group-2 non-invariance (y1's
#      loading, y2's thresholds)
#
# T conventions: lavaan T_lav = sum_g (N_g - 1) F_g; dwlsmv T_NG = (N - G) F,
# F = sum_g pi_g F_g (decided 2026-10-09). T3 = T / c + b with lavaan's
# scaling factor c and shift b (both free of N). Baseline: 99f_'s c1 closed
# form (Result FI2) on both conventions. Indices: FI1 with kappa = sqrt(G),
# n = N - G, F2-F4. The anchor: FI1 on lavaan's convention reproduces
# lavaan's rmsea/cfi/tli.scaled (asserted before anything is written).
#
# Steps: adjacent pairs, H0 the more restricted. dCFI, dRMSEA on the dwlsmv
# convention, H0's value minus H1's (FI3; user 2026-10-31; a worsening
# restriction gives a negative dCFI), F6: a d = 0 step's difference missing,
# the next less restricted step against the last fitted model. Rows most
# restricted first, as difftest's r(steps).
#
# Outputs (fixtures/expected/):
#   lavaan_mg_fitidx_fits.csv   design, level, N1, N2, df, T_lav, T_NG, c, b,
#     T3_lav, T3_NG, T3B_lav, T3B_NG, DB, rmsea_lav, cfi_lav, tli_lav,
#     rmsea_NG, cfi_NG, tli_NG
#   lavaan_mg_fitidx_steps.csv  design, h0, h1, d, dCFI_NG, dRMSEA_NG
#   fixtures/data/mg_P.csv, mg_N.csv   g, y1-y5
#
# Run from the repository root:
#   Rscript fixtures/generate/99j_lavaan_mg_fitidx.R
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

t_start <- Sys.time()
source(file.path("fixtures", "generate", "00_common.R"))
if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT  <- file.path("fixtures", "expected")
DOUT <- file.path("fixtures", "data")
MAN  <- file.path("fixtures", "MANIFEST_lavaan_mg_fitidx.json")
written <- character(0)

TOL_T  <- 1e-9    # sum (N_g - 1) F_g against lavaan's T (99e_'s rebuild tolerance)
TOL_T3 <- 1e-10   # T / c + b against lavaan's T3
TOL_B  <- 1e-8    # c1 closed form against lavaan's baseline (99f_'s)
TOL_A  <- 1e-10   # anchor: FI1 on lavaan's convention against fitMeasures, absolute
SEED_P <- 20261031L
SEED_N <- 20261032L
R_F    <- 0.6     # factor correlation of the P/N DGP

p   <- 5L; y <- paste0("y", seq_len(p))
lam <- c(.60, .65, .70, .75, .80)
NG  <- c(1200L, 2800L)
LEVELS  <- c("configural", "thresholds", "loadings", "intercepts")
SURFACE <- c(configural = "configural", thresholds = "thresholds", loadings = "metric", intercepts = "scalar")
DESIGNS <- list(E = c(4L, 4L, 4L, 2L, 2L), D = c(4L, 2L, 2L, 2L, 2L), P = c(4L, 4L, 4L, 2L, 2L), N = c(4L, 4L, 4L, 2L, 2L))

# ---- data -------------------------------------------------------------------
sim2 <- function(K, seed, noninv) {
  tau <- lapply(seq_len(p), function(j) qnorm(seq_len(K[j] - 1L) / K[j]) + .1 * j)
  set.seed(seed)
  one <- function(n, kap, phi, th, lm, tshift) {
    e1 <- rnorm(n); e2 <- R_F * e1 + sqrt(1 - R_F^2) * rnorm(n)
    eta <- cbind(kap + sqrt(phi) * e1, kap + sqrt(phi) * e2)
    f   <- c(1, 1, 1, 2, 2)
    sapply(seq_len(p), function(j)
      findInterval(lm[j] * eta[, f[j]] + rnorm(n, 0, sqrt(th[j])), tau[[j]] + tshift[j]))
  }
  lm2 <- lam; ts2 <- rep(0, p)
  if (noninv) { lm2[1] <- .30; ts2[2] <- .5 }
  d <- data.frame(g = rep(1:2, NG),
                  rbind(one(NG[1], 0, 1, 1 - lam^2, lam, rep(0, p)),
                        one(NG[2], .3, 1.44, seq(.4, .8, length.out = p), lm2, ts2)))
  names(d)[-1] <- y
  d
}
DATA <- list(E = read.csv(file.path(DOUT, "mg_E.csv")), D = read.csv(file.path(DOUT, "mg_D.csv")),
             P = sim2(DESIGNS$P, SEED_P, FALSE), N = sim2(DESIGNS$N, SEED_N, TRUE))
for (nm in names(DATA)) {
  stopifnot(identical(names(DATA[[nm]]), c("g", y)), identical(as.integer(table(DATA[[nm]]$g)), NG))
  for (g in 1:2) for (j in seq_len(p))
    stopifnot(length(unique(DATA[[nm]][DATA[[nm]]$g == g, j + 1L])) == DESIGNS[[nm]][j])
}

# ---- 99e_'s level syntax, unchanged ------------------------------------------
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
rel <- function(a, b) abs(a / b - 1)

# c1 baseline (99f_'s closed form) from the fit's per-group Stage 1 objects
c1 <- function(fit, wts, Tdiv) {
  obs <- lavInspect(fit, "wls.obs"); gam <- lavInspect(fit, "gamma")
  cc  <- p * (p - 1) / 2; tl <- function(k) (k - cc + 1):k
  F   <- sum(wts / sum(wts) * mapply(function(o, Gm) { k <- length(o); sum(o[tl(k)]^2 / diag(Gm)[tl(k)]) }, obs, gam))
  tr2 <- sum(sapply(gam, function(Gm) { k <- nrow(Gm); sum(cov2cor(Gm[tl(k), tl(k)])^2) }))
  Db  <- length(obs) * cc; a <- sqrt(Db / tr2); TB <- Tdiv * F
  c(T3B = a * TB + Db * (1 - a), Db = Db)
}

# FI1 with F2-F4 (spec), written from the spec text, not from dwlsmv's Mata
fi1 <- function(T3, D, n, G, T3B, DB) {
  EL <- T3 - D; EB <- T3B - DB
  rmsea <- sqrt(G) * sqrt(max(EL, 0) / (n * D))
  if (EB <= 0) return(c(rmsea = rmsea, cfi = NA, tli = NA))
  c(rmsea = rmsea, cfi = 1 - max(EL, 0) / max(EL, EB), tli = 1 - EL * DB / (EB * D))
}

rows <- list(); srows <- list()
for (nm in names(DESIGNS)) {
  K <- DESIGNS[[nm]]; X <- NULL
  for (lv in LEVELS) {
    fit <- fit_level(DATA[[nm]], K, lv)
    stopifnot(lavInspect(fit, "converged"))
    n_g <- as.numeric(unlist(lavInspect(fit, "nobs"))); N <- sum(n_g); G <- length(n_g); pg <- n_g / N
    wv  <- lapply(lavTech(fit, "wls.v"), function(x) if (is.matrix(x)) diag(x) else as.numeric(x))
    obs <- lapply(lavTech(fit, "wls.obs"), as.numeric)
    est <- lapply(lavaan:::lav_model_wls_est(fit@Model), as.numeric)
    Fg  <- mapply(function(x, o, e) sum(x * (o - e)^2), wv, obs, est)
    Tl  <- fit@test$standard$stat; ss <- fit@test$scaled.shifted
    cS  <- ss$scaling.factor; bS <- ss$shift.parameter; T3l <- ss$stat
    df  <- unname(fitMeasures(fit, "df"))
    stopifnot(rel(sum((n_g - 1) * Fg), Tl) < TOL_T, rel(Tl / cS + bS, T3l) < TOL_T3)
    TNG <- (N - G) * sum(pg * Fg); T3N <- TNG / cS + bS
    bL  <- c1(fit, n_g - 1, sum(n_g - 1)); bN <- c1(fit, n_g, N - G)
    fm  <- fitMeasures(fit, c("baseline.chisq.scaled", "baseline.df", "rmsea.scaled", "cfi.scaled", "tli.scaled"))
    stopifnot(rel(bL[["T3B"]], fm[["baseline.chisq.scaled"]]) < TOL_B, bL[["Db"]] == fm[["baseline.df"]])
    iL <- fi1(T3l, df, N - G, G, bL[["T3B"]], bL[["Db"]])
    ia <- max(abs(iL - fm[c("rmsea.scaled", "cfi.scaled", "tli.scaled")]))
    cat(sprintf("%s %-10s df %2d T3_lav %10.5f T3_NG %10.5f  anchor max|FI1 - lavaan| %.2e (tol %.0e)\n",
                nm, lv, df, T3l, T3N, ia, TOL_A))
    stopifnot(ia < TOL_A)
    iN <- fi1(T3N, df, N - G, G, bN[["T3B"]], bN[["Db"]])
    rows[[length(rows) + 1]] <- data.frame(design = nm, level = SURFACE[[lv]], N1 = n_g[1], N2 = n_g[2], df = df,
      T_lav = Tl, T_NG = TNG, c = cS, b = bS, T3_lav = T3l, T3_NG = T3N, T3B_lav = bL[["T3B"]], T3B_NG = bN[["T3B"]],
      DB = bN[["Db"]], rmsea_lav = fm[["rmsea.scaled"]], cfi_lav = fm[["cfi.scaled"]], tli_lav = fm[["tli.scaled"]],
      rmsea_NG = iN[["rmsea"]], cfi_NG = iN[["cfi"]], tli_NG = iN[["tli"]])
    X <- rbind(X, data.frame(level = SURFACE[[lv]], df = df, cfi = iN[["cfi"]], rmsea = iN[["rmsea"]], EL = T3N - df))
  }
  # X is least restricted first (F6's order): the restricted level's value
  # minus the last fitted less restricted level's (FI3; user 2026-10-31); a
  # d = 0 level is skipped. Rows are then written most restricted first, as
  # difftest's r(steps)
  last <- 1; st <- list()
  for (k in 2:nrow(X)) {
    d <- X$df[k] - X$df[k - 1]
    if (d == 0) { dc <- NA; dr <- NA } else {
      dc <- X$cfi[k] - X$cfi[last]; dr <- X$rmsea[k] - X$rmsea[last]; last <- k }
    st[[k - 1]] <- data.frame(design = nm, h0 = X$level[k], h1 = X$level[k - 1], d = d,
      dCFI_NG = dc, dRMSEA_NG = dr, EL_h0 = X$EL[k], EL_h1 = X$EL[k - 1])
  }
  srows <- c(srows, rev(st))
}
fits <- do.call(rbind, rows); steps <- do.call(rbind, srows)
print(steps, digits = 6)

# ---- coverage conditions (prereg 3a.3; lavaan values only) --------------------
cov_pos <- any(steps$dCFI_NG > 1e-4, na.rm = TRUE)
cov_neg <- any(steps$dCFI_NG < -1e-4, na.rm = TRUE)
away    <- function(h) any(steps$d > 0 & steps$h0 == h & steps$EL_h0 > 0 & steps$EL_h1 > 0)   # the step adding h
cov_met <- away("metric"); cov_sca <- away("scalar"); cov_d0 <- any(steps$d == 0)
cat(sprintf("coverage: dCFI > 1e-4 %s; dCFI < -1e-4 %s; metric step off clamp %s; scalar step off clamp %s; d = 0 step %s\n",
            cov_pos, cov_neg, cov_met, cov_sca, cov_d0))
stopifnot(cov_pos, cov_neg, cov_met, cov_sca, cov_d0)

steps$EL_h0 <- NULL; steps$EL_h1 <- NULL
f17 <- function(df) { for (j in names(df)) if (is.double(df[[j]])) df[[j]] <- sprintf("%.17g", df[[j]]); df }
f <- file.path(DOUT, "mg_P.csv");                  write_csv_lf(DATA$P, f);      written <- c(written, f)
f <- file.path(DOUT, "mg_N.csv");                  write_csv_lf(DATA$N, f);      written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_fitidx_fits.csv");  write_csv_lf(f17(fits), f);   written <- c(written, f)
f <- file.path(OUT, "lavaan_mg_fitidx_steps.csv"); write_csv_lf(f17(steps), f);  written <- c(written, f)
write_manifest(MAN, written, seeds = sprintf("%d (design P data), %d (design N data); E, D read mg_E.csv, mg_D.csv", SEED_P, SEED_N),
  extra = list(packages = c(paste("lavaan", as.character(packageVersion("lavaan"))),
    "WLSMV, delta, scaled.shifted, std.lv, zero.add = c(0.5, 0), zero.keep.margins = TRUE, bad.marker.crit = 0 (pinned)",
    "explicit level syntax (99e_ lev_syntax), never group.equal; files use surface level names",
    sprintf("asserted: T rebuild %.0e, T3 = T/c + b %.0e, c1 baseline %.0e, FI1 anchor %.0e abs", TOL_T, TOL_T3, TOL_B, TOL_A),
    R.version.string)))
cat("\nWarnings captured (verbatim):\n"); if (length(WARN)) cat(paste0("  ", WARN), sep = "\n") else cat("  none\n")
cat(sprintf("%d file(s) written; %s written. Elapsed %.1f s\n", length(written), MAN,
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
