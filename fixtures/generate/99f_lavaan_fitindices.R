# ---------------------------------------------------------------------------
# Executed 2026-10-21 from the repository root (R 4.5.2 aarch64, lavaan 0.7.2;
# the other tier-2 generators ran under R 4.6.0, see the manifest): it ran
# clean on the first attempt, every assertion held (closed form against
# lavaan's baseline <= 6.7e-16 relative).
#
# dwlsmv fixture generation -- tier 2: lavaan 0.7.2's scaled fit indices
# (family S, option A; spec Module 12 *Fit indices*, derivation @sec-fit-indices)
#
# Written 2026-10-21 by the fit-indices side session.
#
# Single group (case E, committed fixtures/data/case_E.csv, y1..y6):
#   E_2f  the two-factor correlated model (20_'s M_COR, std.lv)
#   E_1f  a one-factor model on the same data -- misspecified on purpose, so
#         T3 > df and RMSEA, CFI and TLI are away from their clamps
# For each: lavaan's fitMeasures() rows chisq.scaled, df, baseline.chisq,
# baseline.chisq.scaled, baseline.df, rmsea.scaled, cfi.scaled, tli.scaled,
# ntotal. At G = 1 dwlsmv's T = (N - 1) F is lavaan's, so every row is an
# equality target (spec 7.2's test-statistic tier).
#
# Two groups (design E, committed fixtures/data/mg_E.csv; 98_'s configural
# fit): the c1 baseline does not depend on the level, so one fit gives it.
# lavaan's T is sum_g (N_g - 1) F_g; dwlsmv's is (N - G) F with
# F = sum_g pi_g F_g (decided 2026-10-09). Both are written:
#   baseline.chisq, baseline.chisq.scaled  lavaan's own (fitMeasures)
#   TB_lav, T3B_lav   the c1 closed form (Result FI2) on lavaan's Stage 1
#                     objects with lavaan's T -- asserted equal to lavaan's
#                     baseline to 1e-8 relative before anything is written
#   TB_NG, T3B_NG     the same closed form with dwlsmv's T: what dwlsmv
#                     must post as e(T3_b) (not a lavaan quantity)
#   trM2_B            sum_g ||R_g||_F^2 (free of the T convention)
#
# Output: fixtures/expected/lavaan_fitindices.csv  (case, G, quantity, value)
# Manifest: fixtures/MANIFEST_lavaan_fitindices.json (own manifest, as 99e_)
#
# Run from the repository root:
#   Rscript fixtures/generate/99f_lavaan_fitindices.R
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))
if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_lavaan_fitindices.json")
TOL_B <- 1e-8   # closed form against lavaan's fitted baseline (set before running)
FM <- c("chisq.scaled", "df", "baseline.chisq", "baseline.chisq.scaled", "baseline.df",
        "rmsea.scaled", "cfi.scaled", "tli.scaled", "ntotal")

wlsmv_fit <- function(model, dat, ...)
  lavaan::sem(model, data = dat, ordered = setdiff(names(dat), "g"), estimator = "WLSMV",
              parameterization = "delta", test = "scaled.shifted", zero.add = c(0.5, 0.0),
              zero.keep.margins = TRUE, ...)

# c1 baseline closed form (Result FI2) from lavaan's per-group Stage 1 objects
c1 <- function(fit, wts, Tdiv) {
  G   <- lavInspect(fit, "ngroups")
  as1 <- function(x) if (G == 1) list(x) else x
  obs <- as1(lavInspect(fit, "wls.obs")); gam <- as1(lavInspect(fit, "gamma"))
  p   <- length(lavNames(fit, "ov")); cc <- p * (p - 1) / 2
  tl  <- function(k) (k - cc + 1):k
  F   <- sum(wts / sum(wts) * mapply(function(o, Gm) {
           k <- length(o); sum(o[tl(k)]^2 / diag(Gm)[tl(k)]) }, obs, gam))
  tr2 <- sum(sapply(gam, function(Gm) { k <- nrow(Gm); sum(cov2cor(Gm[tl(k), tl(k)])^2) }))
  Db  <- G * cc; a <- sqrt(Db / tr2); TB <- Tdiv * F
  c(TB = TB, T3B = a * TB + Db * (1 - a), Db = Db, trM2_B = tr2)
}

rows <- list()
add  <- function(case, G, q, v) rows[[length(rows) + 1]] <<- data.frame(case = case, G = G, quantity = q, value = v)

datE <- read.csv(file.path("fixtures", "data", "case_E.csv")); datE[] <- lapply(datE, ordered)
M_COR <- "f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5 + y6\n f1 ~~ f2"
M_ONE <- "f =~ y1 + y2 + y3 + y4 + y5 + y6"
for (cs in list(list("E_2f", M_COR), list("E_1f", M_ONE))) {
  fit <- wlsmv_fit(cs[[2]], datE, std.lv = TRUE)
  fm  <- fitMeasures(fit, FM)
  stopifnot(all(is.finite(fm)))
  N   <- lavInspect(fit, "nobs")
  bl  <- c1(fit, N - 1, N - 1)
  rb  <- abs(bl[c("TB", "T3B")] / fm[c("baseline.chisq", "baseline.chisq.scaled")] - 1)
  cat(sprintf("%s: baseline closed form rel %.2e %.2e (tol %.0e); df %g/%g\n",
              cs[[1]], rb[1], rb[2], TOL_B, bl["Db"], fm["baseline.df"]))
  stopifnot(all(rb < TOL_B), bl["Db"] == fm["baseline.df"])
  for (q in FM) add(cs[[1]], 1, q, unname(fm[q]))
}
stopifnot(fitMeasures(wlsmv_fit(M_ONE, datE, std.lv = TRUE), "chisq.scaled") >
          fitMeasures(wlsmv_fit(M_ONE, datE, std.lv = TRUE), "df"))   # E_1f is away from the clamps

datG <- read.csv(file.path("fixtures", "data", "mg_E.csv"))
datG[paste0("y", 1:5)] <- lapply(datG[paste0("y", 1:5)], ordered)
fitG <- wlsmv_fit("f =~ y1 + y2 + y3 + y4 + y5", datG, group = "g", std.lv = TRUE, bad.marker.crit = 0)
fmG  <- fitMeasures(fitG, FM)
Ng   <- lavInspect(fitG, "nobs"); G <- length(Ng)
blL  <- c1(fitG, Ng - 1, sum(Ng - 1))          # lavaan's T
blN  <- c1(fitG, Ng, sum(Ng) - G)              # dwlsmv's T
rb   <- abs(blL[c("TB", "T3B")] / fmG[c("baseline.chisq", "baseline.chisq.scaled")] - 1)
cat(sprintf("mg_E: baseline closed form (lavaan T) rel %.2e %.2e (tol %.0e); df %g/%g\n",
            rb[1], rb[2], TOL_B, blL["Db"], fmG["baseline.df"]))
stopifnot(all(rb < TOL_B), blL["Db"] == fmG["baseline.df"], blN["trM2_B"] == blL["trM2_B"])
for (q in c("baseline.chisq", "baseline.chisq.scaled", "baseline.df", "ntotal")) add("mgE", 2, q, unname(fmG[q]))
add("mgE", 2, "TB_lav", unname(blL["TB"])); add("mgE", 2, "T3B_lav", unname(blL["T3B"]))
add("mgE", 2, "TB_NG", unname(blN["TB"]));  add("mgE", 2, "T3B_NG", unname(blN["T3B"]))
add("mgE", 2, "trM2_B", unname(blN["trM2_B"]))

out <- do.call(rbind, rows)
stopifnot(all(is.finite(out$value)))
f <- file.path(OUT, "lavaan_fitindices.csv")
write_csv_lf(transform(out, value = sprintf("%.17g", value)), f)
write_manifest(MAN, f, seeds = "none (reads committed case_E.csv and mg_E.csv)",
  extra = list(packages = c(paste("lavaan", as.character(packageVersion("lavaan"))),
    "WLSMV, delta, scaled.shifted, std.lv, zero.add = c(0.5, 0), zero.keep.margins = TRUE",
    sprintf("c1 closed form vs lavaan baseline tolerance %.0e relative (set before running)", TOL_B),
    R.version.string)))
cat(sprintf("%d row(s) written to %s; %s written\n", nrow(out), f, MAN))
