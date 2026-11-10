# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan standardizedSolution(), Case E
#
# Module 11 D7, step (2): the reference for the standardization map g(theta)
# and its delta-method SEs (derivation @sec-standardize, whose verification
# carries the same four fits; the choices here are spec D7's, 2026-10-05).
# Case E, F2 <- F1,
# because it is the case with an ENDOGENOUS factor: there the standardizing
# scale of F2 is sqrt(Phi_22) = sqrt(beta^2 Phi_11 + psi_22), not sqrt(psi_22),
# and a single-factor or correlated-factor case cannot tell the two apart.
# The equivalent correlated-factor model (f1 ~~ f2 free) is fitted too, as a
# cross-check: its standardized f1 ~~ f2 is the structural fit's standardized
# beta. Both identifications, options identical to 20_'s fit_E().
#
#   lavaan_case_E_standardized.csv  long: fit, label, lhs, op, rhs, est_std, se
#     fit      cor_unitvar | str_unitvar | cor_marker | str_marker  (20_'s tags)
#     label    lhs op rhs, as lavaan_case_E_parameters.csv labels it
#     est_std  standardizedSolution(type = "std.all")$est.std
#     se       its $se; 0 on rows lavaan treats as fixed after standardizing
#
# TRAP -- THE SEs ARE NOT ANALYTIC. lavaan tries a complex-step Jacobian of
# lav_standardize_all_x, that fails, and it falls back to a FORWARD difference
# with h = sqrt(eps) * max(1, |x|). The committed se therefore carries a
# truncation error of about 1.5e-8 RELATIVE against the analytic J_g of
# @sec-standardize (measured by gammaD/d7_verify.R, lavaan 0.7.2). A
# comparison against this column must use a tolerance justified by that, not
# the 1e-12 floating-point tier; est_std has no such error.
#
# Own manifest (MANIFEST_lavaan_std.json), as 70_ has: 20_ rewrites
# MANIFEST_lavaan.json from its own list, so merging this file into it would
# be dropped on the next 20_ run. fixtures/verify.R needs the matching entry.
#
# Run from the repository root:  Rscript fixtures/generate/80_lavaan_standardized.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_lavaan_std.json")
written <- character(0)

datE <- read.csv(file.path("fixtures", "data", "case_E.csv"))
datE[] <- lapply(datE, ordered)

M_COR <- "f1 =~ y1 + y2 + y3
          f2 =~ y4 + y5 + y6
          f1 ~~ f2"
M_STR <- "f1 =~ y1 + y2 + y3
          f2 =~ y4 + y5 + y6
          f2 ~  f1"

fit_E <- function(model, ...)          # identical to 20_lavaan_reference.R
  lavaan::sem(model, data = datE, ordered = names(datE),
              estimator = "WLSMV", parameterization = "delta",
              test = "scaled.shifted", zero.add = c(0.5, 0.0),
              zero.keep.margins = TRUE, ...)

fits <- list(cor_unitvar = fit_E(M_COR, std.lv = TRUE),
             str_unitvar = fit_E(M_STR, std.lv = TRUE),
             cor_marker  = fit_E(M_COR),
             str_marker  = fit_E(M_STR))
se_type <- lavInspect(fits[[1]], "options")$se

std <- do.call(rbind, Map(function(tag, fit) {
  s <- standardizedSolution(fit, type = "std.all")
  s <- s[s$op %in% c("=~", "~", "~~", "|"), ]
  data.frame(fit = tag, label = paste0(s$lhs, s$op, s$rhs),
             lhs = s$lhs, op = s$op, rhs = s$rhs,
             est_std = s$est.std, se = s$se)
}, names(fits), fits))
rownames(std) <- NULL
stopifnot(!anyNA(std$est_std), !anyNA(std$se))

# Standardized values are a property of the model, not of the metric: marker
# and unit variance must agree on est_std (and, up to the forward-difference
# error, on se). They are separately converged fits, so the gap is the
# optimizer's, not zero. Tolerances chosen after measuring (scratch run
# 2026-10-05, lavaan 0.7.2: est_std gap 4.6e-8 cor / 7.0e-8 str, se gap
# 1.4e-9 / 2.3e-9; cor-vs-str equivalence 1.1e-8, shared rows 4.4e-8); they
# are a stop-the-line check, not a comparison tier.
for (m in c("cor", "str")) {
  u <- std[std$fit == paste0(m, "_unitvar"), ]
  k <- std[std$fit == paste0(m, "_marker"), ]
  stopifnot(identical(u$label, k$label))
  d_est <- max(abs(u$est_std - k$est_std))
  d_se  <- max(abs(u$se - k$se))
  cat(sprintf("%s metric invariance: max|est_std gap| %.2e, max|se gap| %.2e\n",
              m, d_est, d_se))
  stopifnot(d_est < 1e-6, d_se < 1e-5)
}
a <- std[std$fit == "str_unitvar", ]

# The two models are equivalent (one path <-> one covariance), so with a single
# predictor the standardized f2 ~ f1 equals the standardized f1 ~~ f2, and
# every loading, threshold and residual variance agrees. Same optimizer-gap
# tolerance as above, chosen after measuring.
cu <- std[std$fit == "cor_unitvar", ]
r_cor <- cu$est_std[cu$label == "f1~~f2"]
b_str <- a$est_std[a$label == "f2~f1"]
shared <- setdiff(intersect(cu$label, a$label), c("f1~~f2", "f2~f1", "f2~~f2"))
d_eq <- max(abs(cu$est_std[match(shared, cu$label)] - a$est_std[match(shared, a$label)]))
cat(sprintf("equivalence: |cor(f1,f2) - std beta| %.2e, max shared gap %.2e\n",
            abs(r_cor - b_str), d_eq))
stopifnot(abs(r_cor - b_str) < 1e-6, d_eq < 1e-6)

# The endogenous-factor point, asserted: a standardized F2 loading is
# lambda * sqrt(Phi_22), and here that differs from lambda * sqrt(psi_22).
fu  <- fits$str_unitvar
lam <- lavInspect(fu, "est")$lambda["y4", "f2"]
phi <- lavInspect(fu, "cov.lv")["f2", "f2"]
psi <- lavInspect(fu, "est")$psi["f2", "f2"]
e44 <- a$est_std[a$label == "f2=~y4"]
stopifnot(abs(e44 - lam * sqrt(phi)) < 1e-12, abs(e44 - lam * sqrt(psi)) > 1e-2)

f <- file.path(OUT, "lavaan_case_E_standardized.csv")
# LF on every platform (decided 2026-10-05), unlike 20_/70_'s write.csv: the
# hash does not depend on the generating machine.
write_csv_lf(std, f)
written <- c(written, f)

write_manifest(MAN, written,
  seeds = "none (deterministic given case_E.csv)",
  extra = list(packages = c(
    paste("lavaan", as.character(packageVersion("lavaan"))),
    paste("se type:", se_type),
    "standardizedSolution(type = 'std.all'); se by lavaan's forward-difference Jacobian (complex step fails), ~1.5e-8 relative vs analytic",
    "zero.keep.margins = TRUE (pinned, as 20_)",
    R.version.string)))

cat(sprintf("%d standardized reference file(s) written; %s written.\n", length(written), MAN))
