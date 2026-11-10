# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan reference outputs
#
# Executed 2026-09-18 under R 4.6.0 with lavaan 0.7.2, on the development
# machine. It ran clean on the first attempt: nothing needed fixing, and both
# CHECK items below are now resolved in place. The earlier UNVERIFIED banner
# is removed on that basis.
#
# TRAP, recorded because the two files sit side by side in fixtures/expected/
# and look like the same quantity:
#
#   case_D_thresholds.csv       tier 1, POPULATION thresholds -- qnorm of the
#                               cumulative declared category probabilities.
#                               y1 and y2 are bit-identical there, and so are
#                               y3 and y4, because they are design constants.
#   lavaan_case_D_th.csv        tier 2, SAMPLE estimates from the 5000 rows.
#
# They differ by sampling error -- about 1.8 SE on y1|t1 -- and diffing them
# as though they were the same object produces a false alarm. Module 1 must
# reproduce the tier-2 numbers, not the tier-1 ones; tier 1 is a bias check,
# not an equality target.
#
# Why this tier is separate from tier 1
# ------------------------------------
# Tier 1 needs only base R and is the reference implementation of the
# mathematics. Tier 2 pulls in an external package whose behaviour can change
# between releases. Keeping them apart means:
#   * tier-1 fixtures can be regenerated on any machine, forever;
#   * a lavaan upgrade changes only tier-2 files, and the manifest says so;
#   * the Stata test suite never needs R at all -- it reads committed files.
#
# Run from the repository root:  Rscript fixtures/generate/20_lavaan_reference.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures. Tier 1 (10_fixtures.R) does not need it.")
library(lavaan)

OUT <- file.path("fixtures", "expected")
written <- character(0)

# Case D is the single-factor reference: it is the smallest case with an
# identified factor model. Cases A-C have too few indicators for a
# single-factor model to be identified, which is why Module 0 uses them for
# Gamma-level checks rather than end-to-end fits. Case E, fitted at the bottom
# of this file, is the two-factor case, and it is here for the one thing a
# single-factor design structurally cannot answer -- see its header.
dat <- read.csv(file.path("fixtures", "data", "case_D.csv"))
dat[] <- lapply(dat, ordered)

fit <- lavaan::sem(
  "f =~ y1 + y2 + y3 + y4",
  data             = dat,
  ordered          = names(dat),
  estimator        = "WLSMV",
  parameterization = "delta",
  test             = "scaled.shifted",
  zero.add         = c(0.5, 0.0),    # lavaan's documented default, stated explicitly
  # lavaan's effective default (set by lav_options_mimic), pinned 2026-09-25:
  # with it, the 2x2 zero.add moves 0.5 and keeps the margins.
  zero.keep.margins = TRUE
)

# --- the statistic vector's ordering is the authoritative one ---------------
# Gamma is by definition the asymptotic covariance of this vector, so its
# names ARE Gamma's ordering. Prefer this over dimnames(gamma), which is a
# convenience a future release could drop.
s_labels <- names(lavInspect(fit, "wls.obs"))

# RESOLVED 2026-09-18 against lavaan 0.7.2: the ordering matches Section 4.2
# exactly, all 16 labels. This is the same sequence pinned as the oracle in
# tests/unit/test_statmap.do, so Module 0-1 is now verified against a live
# lavaan rather than against the spec transcription of it. The warning below
# stays armed: it is what catches a future lavaan reordering.
expected_labels <- c(threshold_labels(c(3, 3, 2, 2)), correlation_labels(4))
if (!identical(s_labels, expected_labels))
  warning("lavaan's statistic ordering differs from the spec's Section 4.2 convention.\n",
          "  lavaan: ", paste(s_labels, collapse = " "), "\n",
          "  spec:   ", paste(expected_labels, collapse = " "), "\n",
          "  Do not silently reorder -- update Section 4.2 and this check together.")

objects <- list(
  gamma    = lavInspect(fit, "gamma"),     # N x acov of the sample statistics
  wls_obs  = lavInspect(fit, "wls.obs"),   # s
  wls_est  = lavInspect(fit, "wls.est"),   # sigma(theta-hat), NOT s
  th       = lavInspect(fit, "th"),        # thresholds
  sampstat = lavInspect(fit, "sampstat")$cov  # polychoric correlations
  # RESOLVED 2026-09-18. "th" comes back named exactly as s_labels names its
  # threshold entries, and it is on the standardized N(0,1) scale delta
  # implies. Verified precisely rather than by eye: those thresholds equal
  # qnorm(empirical cumulative proportion) computed from case_D.csv directly,
  # to 3.3e-15 across all ten. That is what Module 1 has to reproduce, and it
  # says the first stage is that and nothing more.
)

for (nm in names(objects)) {
  x <- objects[[nm]]
  f <- file.path(OUT, sprintf("lavaan_case_D_%s.csv", nm))
  if (is.matrix(x)) write.csv(cbind(label = rownames(x), as.data.frame(x)),
                              f, row.names = FALSE, quote = FALSE)
  else              write.csv(data.frame(label = names(x), value = as.numeric(x)),
                              f, row.names = FALSE, quote = FALSE)
  written <- c(written, f)
}

# --- parameter estimates and standard errors --------------------------------
# These are NOT derivable from anything above. wls.est is sigma(theta-hat), the
# model-implied sample statistics -- not theta-hat itself -- and an earlier
# version of this file saved only that, which would have left Modules 7 and 8
# with nothing to compare against.
#
# The SEs are the robust (sandwich) ones, because that is what lavaan uses by
# default under WLSMV; the se= setting is recorded below so a future default
# change is visible rather than silent. Those are Module 8's target, and the
# naive SEs are not an acceptable substitute for them.
#
# Identification is lavaan's default marker method -- first loading fixed to 1.
# fixtures/generate/30_mplus_reference.R deliberately matches it, so the two
# references are comparable element by element. Changing it here without
# changing it there makes the loadings differ by a factor and look like a bug.
se_type <- lavInspect(fit, "options")$se
cat(sprintf("lavaan se type: %s\n", se_type))

pe <- lavaan::parameterEstimates(fit)
pe$label <- paste0(pe$lhs, pe$op, pe$rhs)
f <- file.path(OUT, "lavaan_case_D_parameters.csv")
write.csv(pe[, c("label", "lhs", "op", "rhs", "est", "se", "z", "pvalue")],
          f, row.names = FALSE, quote = FALSE)
written <- c(written, f)

# Scalar fit results, for Module 9's eventual comparison.
ft <- fitMeasures(fit, c("chisq.scaled", "df.scaled", "pvalue.scaled"))
f <- file.path(OUT, "lavaan_case_D_test.csv")
write.csv(data.frame(label = names(ft), value = as.numeric(ft)),
          f, row.names = FALSE, quote = FALSE)
written <- c(written, f)

# ===========================================================================
# Case E -- the endogenous-factor metric, and the equivalent-models check
#
# Two separate things, fitted on the same two-factor data because they need
# the same ingredient: a factor that is PREDICTED by another factor.
#
# (1) THE METRIC QUESTION. Under unit-variance identification, what is fixed
#     to 1 for an endogenous factor -- its TOTAL variance, or its DISTURBANCE
#     variance? For an exogenous factor the two coincide and the question does
#     not arise. Once F2 is regressed on F1 they differ by beta^2, so the
#     answer changes every loading on F2 and is not a matter of taste.
#     Measured, not reasoned about: the checks below read it straight out of
#     lavaan's parameter table and its implied factor covariance matrix.
#
# (2) THE EQUIVALENT-MODELS CHECK. A correlated two-factor model and an
#     F2 <- F1 model are reparameterizations of each other -- both spend one
#     parameter on the link between the factors, and each maps onto the other
#     without changing sigma(theta). Fitted to the same data they must
#     therefore produce the same fitted statistics, the same df, and the same
#     T3. That identity holds regardless of how (1) comes out, which is what
#     makes it a usable internal consistency check for Module 7 rather than a
#     second way of asking the same question.
#
# Both are committed as fixtures rather than left as a note, so a lavaan
# release that changes the convention breaks a test instead of silently
# changing what spec Section 4.8 documents.
# ===========================================================================

datE <- read.csv(file.path("fixtures", "data", "case_E.csv"))
datE[] <- lapply(datE, ordered)

M_COR <- "f1 =~ y1 + y2 + y3
          f2 =~ y4 + y5 + y6
          f1 ~~ f2"
M_STR <- "f1 =~ y1 + y2 + y3
          f2 =~ y4 + y5 + y6
          f2 ~  f1"

fit_E <- function(model, ...)
  lavaan::sem(model, data = datE, ordered = names(datE),
              estimator = "WLSMV", parameterization = "delta",
              test = "scaled.shifted", zero.add = c(0.5, 0.0),
              zero.keep.margins = TRUE, ...)

fits <- list(
  cor_unitvar = fit_E(M_COR, std.lv = TRUE),
  str_unitvar = fit_E(M_STR, std.lv = TRUE),
  cor_marker  = fit_E(M_COR),
  str_marker  = fit_E(M_STR)
)

# --- (1) what is fixed, and what the fixing implies -------------------------
# psi is the (co)variance matrix of the DISTURBANCES; cov.lv is the implied
# covariance matrix of the FACTORS, Phi = (I-B)^-1 Psi (I-B)^-T. For an
# exogenous factor they agree. The whole question is which one carries the 1
# when the factor is endogenous, so both are recorded, per fit.
metric_rows <- function(tag, fit) {
  pt  <- parTable(fit)
  psi <- lavInspect(fit, "est")$psi
  phi <- lavInspect(fit, "cov.lv")
  bet <- lavInspect(fit, "est")$beta
  row <- function(q, v) data.frame(fit = tag, quantity = q, value = v)
  p22 <- pt[pt$lhs == "f2" & pt$op == "~~" & pt$rhs == "f2", ]
  b21 <- if (is.null(bet)) 0 else bet["f2", "f1"]
  rbind(
    row("psi_f2_f2_is_free",  as.numeric(p22$free != 0)),
    row("psi_f2_f2_ustart",   if (is.na(p22$ustart)) NA_real_ else p22$ustart),
    row("psi_f2_f2_est",      psi["f2", "f2"]),
    row("beta_f2_f1",         b21),
    row("var_f2_implied",     phi["f2", "f2"]),
    row("var_f1_implied",     phi["f1", "f1"]),
    row("cov_f1_f2_implied",  phi["f1", "f2"])
  )
}

metric <- do.call(rbind, Map(metric_rows, names(fits), fits))
f <- file.path(OUT, "lavaan_case_E_metric.csv")
write.csv(metric, f, row.names = FALSE, quote = FALSE)
written <- c(written, f)

get <- function(tag, q)
  metric$value[metric$fit == tag & metric$quantity == q]

# THE ANSWER, asserted rather than eyeballed. If any of these three stops
# holding, lavaan changed its convention and spec Section 4.8 is now wrong.
stopifnot(get("str_unitvar", "psi_f2_f2_is_free") == 0)   # it is fixed
stopifnot(get("str_unitvar", "psi_f2_f2_est") == 1)       # fixed to 1
stopifnot(abs(get("str_unitvar", "var_f2_implied") -
              (get("str_unitvar", "beta_f2_f1")^2 +
               get("str_unitvar", "psi_f2_f2_est"))) < 1e-12)
# ... and the total variance is consequently NOT 1. This is the line that
# distinguishes the two candidate answers; the ones above are consistent with
# either.
stopifnot(abs(get("str_unitvar", "var_f2_implied") - 1) > 0.1)
# For the EXOGENOUS factor in the same fit, and for both factors in the
# correlated fit, the two readings coincide -- which is exactly why a
# single-factor or correlated-factor test cannot settle this.
stopifnot(abs(get("str_unitvar", "var_f1_implied") - 1) < 1e-12)
stopifnot(abs(get("cor_unitvar", "var_f2_implied") - 1) < 1e-12)

cat(sprintf(paste0("\ncase E, unit-variance identification, F2 <- F1:\n",
                   "  disturbance var(zeta_2) = %.15f   (fixed: %s)\n",
                   "  structural path beta    = %.15f\n",
                   "  implied total Var(F2)   = %.15f\n",
                   "  => lavaan fixes the DISTURBANCE variance, not the",
                   " total variance.\n"),
            get("str_unitvar", "psi_f2_f2_est"),
            if (get("str_unitvar", "psi_f2_f2_is_free") == 0) "yes" else "no",
            get("str_unitvar", "beta_f2_f1"),
            get("str_unitvar", "var_f2_implied")))

# --- (2) the equivalent-models check ----------------------------------------
# Column names are written with underscores, not lavaan's dots: `import
# delimited` strips a dot from a variable name, so `chisq.scaled` and
# `chisq_scaled` would both arrive in Stata as `chisqscaled` and the two df
# columns would collide. Renaming here keeps the Stata side legible.
equiv <- do.call(rbind, Map(function(tag, fit) {
  ft <- fitMeasures(fit, c("chisq", "df", "chisq.scaled", "df.scaled",
                           "pvalue.scaled", "npar"))
  data.frame(fit = tag, chisq = ft[["chisq"]], df = ft[["df"]],
             chisq_scaled = ft[["chisq.scaled"]],
             df_scaled = ft[["df.scaled"]],
             pvalue_scaled = ft[["pvalue.scaled"]],
             npar = ft[["npar"]])
}, names(fits), fits))
f <- file.path(OUT, "lavaan_case_E_equivalence.csv")
write.csv(equiv, f, row.names = FALSE, quote = FALSE)
written <- c(written, f)

# The identity, asserted at spec Section 7.2's test-statistic tier: 1e-5
# relative + 1e-6 absolute since 2026-09-24 (it was 1e-3 relative before).
# The measured gap is many orders tighter than that; what caps it is that the
# two models are separately converged optimizations, not the
# reparameterization, which is exact.
for (idm in c("unitvar", "marker")) {
  a <- equiv[equiv$fit == paste0("cor_", idm), ]
  b <- equiv[equiv$fit == paste0("str_", idm), ]
  stopifnot(a$df_scaled == b$df_scaled, a$npar == b$npar)
  rel <- abs(a$chisq_scaled - b$chisq_scaled) / a$chisq_scaled
  stopifnot(abs(a$chisq_scaled - b$chisq_scaled) <= 1e-6 + 1e-5 * abs(a$chisq_scaled))
  cat(sprintf("case E equivalence, %-8s: T3 %.12f vs %.12f, df %g, rel gap %.2e\n",
              idm, a$chisq_scaled, b$chisq_scaled, a$df_scaled, rel))
}

# Fitted sigma(theta) is the sharper statement of the same identity: the two
# models do not merely fit equally well, they imply the SAME statistics.
for (idm in c("unitvar", "marker")) {
  g <- max(abs(lavInspect(fits[[paste0("cor_", idm)]], "wls.est") -
               lavInspect(fits[[paste0("str_", idm)]], "wls.est")))
  stopifnot(g < 1e-6)
  cat(sprintf("case E equivalence, %-8s: max |sigma(theta) gap| = %.2e\n",
              idm, g))
}

# --- parameter estimates, for Modules 7 and 8 later -------------------------
pe_E <- do.call(rbind, Map(function(tag, fit) {
  p <- lavaan::parameterEstimates(fit)
  p$label <- paste0(p$lhs, p$op, p$rhs)
  cbind(fit = tag, p[, c("label", "lhs", "op", "rhs", "est", "se", "z", "pvalue")])
}, names(fits), fits))
f <- file.path(OUT, "lavaan_case_E_parameters.csv")
write.csv(pe_E, f, row.names = FALSE, quote = FALSE)
written <- c(written, f)

# ===========================================================================
# Module 8 -- lavaan's OWN sandwich inputs and output, per fit
#
# Added and executed 2026-09-23 (R 4.6.0, lavaan 0.7.2), followed by 60_ as
# the run order requires; every pre-existing tier-2 file regenerated
# byte-identical. The SEs above are the end-to-end target; these make the
# "formula with lavaan's own inputs" test possible, the split test_gamma.do
# 5a/5b uses for Gamma. Fed lavaan's Delta (analytic, at lavaan's theta-hat)
# and lavaan's Gamma, dwlsmv_acov() and dwlsmv_ev() must reproduce lavaan's
# vcov at floating-point agreement -- no theta-hat or rho-hat tier inherited.
#
#   lavaan_acov_delta.csv   long: fit, row (Section 4.2 statistic label),
#                           col (lavaan parameter label lhs op rhs), value
#   lavaan_acov_vcov.csv    long: fit, row, col (parameter labels), value --
#                           lavInspect's vcov, i.e. ALREADY divided by N - 1
#   lavaan_case_E_gamma.csv as lavaan_case_D_gamma.csv; one Gamma serves all
#                           four case E fits, asserted identical below
#
# Parameter labels are the ones lavaan_case_*_parameters.csv carries, so a
# consumer maps them with the same routine (test_fit.do's tf_lab). Column
# order carries no meaning; it is lavaan's free order, loadings first.
#
# The formula is asserted here before anything is written -- spec Section
# 4.10 and derivation @sec-acov: nV = E^-1 Delta'W Gamma W Delta E^-1 with
# W = 1/diag(Gamma), vcov = nV / (N - 1). If lavaan changed its divisor or
# its weight, this stops rather than committing a fixture that encodes it.
# ===========================================================================
acov_fits <- c(list(D_marker = fit), fits)
rows_delta <- list(); rows_vcov <- list()
for (tag in names(acov_fits)) {
  fa   <- acov_fits[[tag]]
  pt   <- parTable(fa); fr <- pt[pt$free > 0, ]; fr <- fr[order(fr$free), ]
  plab <- paste0(fr$lhs, fr$op, fr$rhs)
  Dl   <- lavTech(fa, "delta")[[1]]
  Gm   <- lavTech(fa, "gamma")[[1]]
  V    <- lavTech(fa, "vcov")
  N    <- lavInspect(fa, "nobs")
  slab <- names(lavInspect(fa, "wls.obs"))
  stopifnot(ncol(Dl) == length(plab), nrow(Dl) == length(slab),
            all(dim(V) == length(plab)),
            max(abs(fa@SampleStats@WLS.VD[[1]] - 1 / diag(Gm))) == 0)
  WD <- (1 / diag(Gm)) * Dl
  Ei <- solve(crossprod(Dl, WD))
  nV <- Ei %*% crossprod(WD, Gm %*% WD) %*% Ei
  g  <- max(abs(nV / (N - 1) - V)) / max(abs(V))
  if (g > 1e-12)
    stop(tag, ": lavaan's vcov is not the Section 4.10 sandwich over N - 1 (gap ", g, ")")
  cat(sprintf("%-12s sandwich/(N-1) vs lavaan vcov: %.2e  (/N would be %.2e)\n",
              tag, g, max(abs(nV / N - V)) / max(abs(V))))
  rows_delta[[tag]] <- data.frame(fit = tag,
    row = rep(slab, times = length(plab)), col = rep(plab, each = length(slab)),
    value = as.numeric(Dl))
  rows_vcov[[tag]] <- data.frame(fit = tag,
    row = rep(plab, times = length(plab)), col = rep(plab, each = length(plab)),
    value = as.numeric(V))
}

gE <- lapply(fits, function(fa) lavInspect(fa, "gamma"))
stopifnot(all(vapply(gE, function(x) identical(unclass(x), unclass(gE[[1]])), logical(1))))

for (nm in c("delta", "vcov")) {
  f <- file.path(OUT, sprintf("lavaan_acov_%s.csv", nm))
  write.csv(do.call(rbind, if (nm == "delta") rows_delta else rows_vcov),
            f, row.names = FALSE, quote = FALSE)
  written <- c(written, f)
}
# ===========================================================================
# Module 9 -- lavaan's scaled.shifted components, per fit
#
# Added and executed 2026-09-24 (R 4.6.0, lavaan 0.7.2), then 60_; every
# pre-existing tier-2 file regenerated byte-identical. Pairs with the Module 8
# block above: fed lavaan's own
# Delta (lavaan_acov_delta.csv), Gamma (lavaan_case_*_gamma.csv) and T, the
# test-statistic routine must reproduce these at floating-point agreement.
#
#   lavaan_t3.csv   long: fit, quantity, value, with quantity one of
#     T       the unscaled statistic being scaled, fit@test$standard$stat,
#             which is (N-1) F (asserted in test_fit.do check 5)
#     df      the model degrees of freedom, D
#     trM     Tr(U Gamma)          fit@test$scaled.shifted$trace.UGamma
#     trM2    Tr((U Gamma)^2)      ...$trace.UGamma2
#     a       1 / scaling.factor   (fitMeasures exposes only 1/a)
#     b       shift.parameter      (stored ONLY in fit@test, not fitMeasures)
#     T3      a T + b              ...$stat, = fitMeasures chisq.scaled
#
# Asserted before writing: the trace expansion of derivation
# @sec-gate5-traces, rebuilt here from lavaan's Delta and Gamma, reproduces
# lavaan's stored traces, a, b and T3 -- so a lavaan release that changed its
# construction stops this script rather than moving a fixture.
# ===========================================================================
trace_noU <- function(Gm, Dl) {
  w  <- 1 / diag(Gm)
  WD <- w * Dl
  Ei <- solve(crossprod(Dl, WD))
  G  <- Gm %*% WD
  B  <- crossprod(WD, G)
  C  <- crossprod(G, w * G)
  EB <- Ei %*% B
  c(trM  = sum(w * diag(Gm)) - sum(Ei * B),
    trM2 = drop(crossprod(w, (Gm * Gm) %*% w)) - 2 * sum(Ei * C) + sum(EB * t(EB)))
}
rows_t3 <- list()
for (tag in names(acov_fits)) {
  fa <- acov_fits[[tag]]
  ss <- fa@test$scaled.shifted
  stopifnot(!is.null(ss), ss$df == fa@test$standard$df)
  tr <- trace_noU(lavTech(fa, "gamma")[[1]], lavTech(fa, "delta")[[1]])
  a  <- 1 / ss$scaling.factor
  b  <- ss$shift.parameter
  Tr <- fa@test$standard$stat
  g  <- c(trM  = abs(tr[["trM"]]  / ss$trace.UGamma  - 1),
          trM2 = abs(tr[["trM2"]] / ss$trace.UGamma2 - 1),
          a    = abs(sqrt(ss$df / tr[["trM2"]]) / a - 1),
          b    = abs(ss$df - sqrt(ss$df / tr[["trM2"]]) * tr[["trM"]] - b),
          T3   = abs((a * Tr + b) / ss$stat - 1))
  if (any(g > 1e-12))
    stop(tag, ": lavaan's scaled.shifted is not the @sec-gate5-traces expansion: ",
         paste(sprintf("%s %.2e", names(g), g), collapse = ", "))
  cat(sprintf("%-12s T3 components rebuilt: %s\n", tag,
              paste(sprintf("%s %.1e", names(g), g), collapse = "  ")))
  rows_t3[[tag]] <- data.frame(fit = tag,
    quantity = c("T", "df", "trM", "trM2", "a", "b", "T3"),
    value = c(Tr, ss$df, ss$trace.UGamma, ss$trace.UGamma2, a, b, ss$stat))
}
f <- file.path(OUT, "lavaan_t3.csv")
write.csv(do.call(rbind, rows_t3), f, row.names = FALSE, quote = FALSE)
written <- c(written, f)

x <- gE[[1]]
f <- file.path(OUT, "lavaan_case_E_gamma.csv")
write.csv(cbind(label = rownames(x), as.data.frame(unclass(x))),
          f, row.names = FALSE, quote = FALSE)
written <- c(written, f)

# lavaan's case E rho-hat, for test_gamma.do check 5c (Gamma at a MATCHED
# rho-hat). Same object and writer as lavaan_case_D_sampstat.csv. Like Gamma,
# one sampstat serves all four case E fits, asserted identical. Added
# 2026-09-25: case E's end-to-end Gamma gap (6.0e-8) is lavaan's rho-hat
# convergence, and showing that in Stata needs lavaan's rho-hat.
sE <- lapply(fits, function(fa) lavInspect(fa, "sampstat")$cov)
stopifnot(all(vapply(sE, function(s) identical(unclass(s), unclass(sE[[1]])), logical(1))))
x <- sE[[1]]
f <- file.path(OUT, "lavaan_case_E_sampstat.csv")
write.csv(cbind(label = rownames(x), as.data.frame(x)), f, row.names = FALSE, quote = FALSE)
written <- c(written, f)

write_manifest(file.path("fixtures", "MANIFEST_lavaan.json"), written,
               seeds = "none (deterministic given case_D.csv and case_E.csv)",
               extra = list(packages = c(
                 paste("lavaan", as.character(packageVersion("lavaan"))),
                 paste("se type:", se_type),
                 R.version.string)))

cat(sprintf("%d lavaan reference files written.\n", length(written)))
cat("Record the lavaan version in the manifest with every regeneration --\n")
cat("a reference output that moved because the package moved is not a test failure.\n")
