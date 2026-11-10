# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: unrotated two-factor EFA, Case E
#
# ESEM groundwork (v3), step 0: does v1's pipeline already fit an UNROTATED
# m-factor EFA? An unrotated EFA is a CFA with every loading free except
# m(m-1)/2 echelon anchor zeros, Phi = I. Stage 1, Gamma, W_D, F, T3 and df
# are rotation-invariant, so no new mathematics is involved; this file is the
# reference that tests it. Rotation and rotated SEs are out of scope.
#
# Case E: six 3-category indicators, the only six-item case. m = 2 gives
# 12 - 1 = 11 free loadings against 15 correlations: df = 4.
#
# Two lavaan fits, which must agree before anything is written:
#   explicit  the echelon CFA dwlsmv is given: f2 =~ 0*y1, f1 ~~ 0*f2,
#             std.lv = TRUE. Same anchors as the Stata test, so its loadings
#             are comparable up to column signs.
#   efa       lavaan's own efa block, rotation = "none". Its internal anchors
#             are lavaan's, not ours; only rotation-invariant quantities
#             (T, T3, df, sigma, Theta) are compared with it.
# Options otherwise identical to 20_'s fit_E().
#
#   lavaan_efa_case_E.csv  long: quantity, label, value
#     test    T (= (N-1) F), df, a, b, T3    explicit fit
#     sigma   the 15 model-implied correlations, label "yi~~yj" (i < j)
#     theta   the 6 residual variances 1 - (Lambda Lambda')_jj, "yj"
#     lambda  explicit loadings, label "fk=~yj" (the fixed 0 included)
#     tau     thresholds, "yj|tk"
#
# Run from the repository root:  Rscript fixtures/generate/99g_lavaan_efa.R
# Own manifest (MANIFEST_lavaan_efa.json), for 80_'s reason.
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_lavaan_efa.json")

datE <- read.csv(file.path("fixtures", "data", "case_E.csv"))
datE[] <- lapply(datE, ordered)
yn <- names(datE)

M_EXPL <- "f1 =~ y1 + y2 + y3 + y4 + y5 + y6
           f2 =~ 0*y1 + y2 + y3 + y4 + y5 + y6
           f1 ~~ 0*f2"
M_EFA  <- 'efa("b")*f1 + efa("b")*f2 =~ y1 + y2 + y3 + y4 + y5 + y6'

fit_E <- function(model, ...)          # 20_'s fit_E(), plus std.lv
  lavaan::sem(model, data = datE, ordered = yn,
              estimator = "WLSMV", parameterization = "delta",
              test = "scaled.shifted", zero.add = c(0.5, 0.0),
              zero.keep.margins = TRUE, std.lv = TRUE, ...)

fx <- fit_E(M_EXPL)
fe <- fit_E(M_EFA, rotation = "none")
stopifnot(lavInspect(fx, "converged"), lavInspect(fe, "converged"))

sig <- function(f) { S <- lavInspect(f, "cov.ov")[yn, yn]; S[lower.tri(S)] }
st  <- function(f) f@test$scaled.shifted
ss_x <- st(fx); ss_e <- st(fe)

# Rotation invariance between the two anchorings: separately converged
# optimizations, so the gap is the optimizer's. Tolerances chosen after
# measuring (2026-10-24, lavaan 0.7.2: T 8.8e-13, T3 1.5e-10, sigma 7.8e-9);
# a stop-the-line check, not a tier.
g_T   <- abs(fx@test$standard$stat / fe@test$standard$stat - 1)
g_T3  <- abs(ss_x$stat / ss_e$stat - 1)
g_sig <- max(abs(sig(fx) - sig(fe)))
cat(sprintf("explicit vs efa(): T rel %.2e, T3 rel %.2e, sigma max %.2e, df %d vs %d\n",
            g_T, g_T3, g_sig, ss_x$df, ss_e$df))
stopifnot(ss_x$df == 4, ss_e$df == 4, g_T < 1e-5, g_T3 < 1e-5, g_sig < 1e-6)

L  <- lavInspect(fx, "est")$lambda[yn, c("f1", "f2")]
S  <- lavInspect(fx, "cov.ov")[yn, yn]
th <- lavInspect(fx, "est")$tau
ij <- which(lower.tri(S), arr.ind = TRUE)          # column-major, j < i

rows <- rbind(
  data.frame(quantity = "test", label = c("T", "df", "a", "b", "T3"),
             value = c(fx@test$standard$stat, ss_x$df,
                       1 / ss_x$scaling.factor, ss_x$shift.parameter, ss_x$stat)),
  data.frame(quantity = "sigma",
             label = paste0(yn[ij[, 2]], "~~", yn[ij[, 1]]),
             value = S[lower.tri(S)]),
  data.frame(quantity = "theta", label = yn,
             value = 1 - rowSums(L^2)),
  data.frame(quantity = "lambda",
             label = paste0(rep(c("f1", "f2"), each = 6), "=~", yn),
             value = c(L[, 1], L[, 2])),
  data.frame(quantity = "tau", label = rownames(th), value = th[, 1]))
rownames(rows) <- NULL
stopifnot(!anyNA(rows$value),
          max(abs(rows$value[rows$quantity == "theta"] -
                  diag(lavInspect(fx, "est")$theta)[yn])) < 1e-12)

f <- file.path(OUT, "lavaan_efa_case_E.csv")
write_csv_lf(rows, f)

write_manifest(MAN, f,
  seeds = "none (deterministic given case_E.csv)",
  extra = list(packages = c(
    paste("lavaan", as.character(packageVersion("lavaan"))),
    "estimator WLSMV, parameterization delta, test scaled.shifted, std.lv TRUE",
    "missing: listwise (default; case_E complete)",
    "zero.add c(0.5, 0), zero.keep.margins TRUE (as 20_)",
    "explicit echelon: f2 =~ 0*y1, f1 ~~ 0*f2; cross-checked vs efa() rotation none",
    R.version.string)))

cat("lavaan_efa_case_E.csv and", MAN, "written.\n")
