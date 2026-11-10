# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: oblique geomin EFA (0.3.0, Unit 4)
#
# lavaan's efa(rotation = "geomin") under the conventions the user settled on
# 2026-10-29 (spec section 5, Module 12; derivation @sec-geomin):
#   epsilon   Mplus's rule: .0001 at m = 2, .001 at m = 3, .01 at m >= 4
#   rows      no normalization (row_weights "none"); std_ov TRUE, the
#             identity under the delta parameterization (@sec-geomin GM1)
#   starts    lavaan's 30 random starts, seeded here so the file reproduces
#             (dwlsmv's own starts are the GM4 grid; the fixture is the
#             reference minimum, not the start rule)
#   tol       gpa_tol 1e-6, above the stall floor of about 1e-7
#   order     lavaan's reflection (positive column sums) and "index" order,
#             which is dwlsmv's convention
#
#   Case E  six 3-category items, m = 2      (fixtures/data/case_E.csv)
#   Case F  nine 4-category items, m = 3     (91_efa_data.R)
#   Case G  six binary items, m = 2          (91_efa_data.R)
#
# Fit options otherwise as 99g_ (delta, scaled.shifted, zero.add c(.5, 0)).
#
#   lavaan_efa_geomin.csv  long: case, quantity, label, value
#     test    T (= (N-1) F), df, T3
#     fit     rmsea.scaled, cfi.scaled, tli.scaled (fitMeasures; Unit 6,
#             added 2026-10-29: rotation-invariant, GM1)
#     rot     eps, Q (the geomin criterion at the rotated loadings)
#     lambda  rotated loadings, "fk=~yj"
#     phi     factor correlations, lower triangle, "fk~~fl" (k > l)
#
# Run from the repository root:  Rscript fixtures/generate/99h_lavaan_efa_geomin.R
# Own manifest (MANIFEST_lavaan_efa_geomin.json), for 80_'s reason.
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("lavaan", quietly = TRUE))
  stop("lavaan is required for tier-2 fixtures.")
library(lavaan)

OUT  <- file.path("fixtures", "expected")
MAN  <- file.path("fixtures", "MANIFEST_lavaan_efa_geomin.json")
SEED <- 20261029

eps_rule <- function(m) if (m == 2) 1e-4 else if (m == 3) 1e-3 else 1e-2
geomin_Q <- function(L, eps) sum(exp(rowSums(log(L^2 + eps)) / ncol(L)))

CASES <- list(E = list(file = "case_E.csv", m = 2),
              F = list(file = "case_F.csv", m = 3),
              G = list(file = "case_G.csv", m = 2))

fit_efa <- function(dat, m, rot, ...) {
  yn  <- names(dat)
  lhs <- paste(sprintf('efa("b")*f%d', seq_len(m)), collapse = " + ")
  lavaan::sem(paste(lhs, "=~", paste(yn, collapse = " + ")), data = dat,
              ordered = yn, estimator = "WLSMV", parameterization = "delta",
              test = "scaled.shifted", zero.add = c(0.5, 0.0),
              zero.keep.margins = TRUE, rotation = rot, ...)
}

rows <- list()
seeds <- integer(0)
for (nm in names(CASES)) {
  cs  <- CASES[[nm]]
  dat <- read.csv(file.path("fixtures", "data", cs$file))
  dat[] <- lapply(dat, ordered)
  yn  <- names(dat); m <- cs$m; eps <- eps_rule(m)

  f0 <- fit_efa(dat, m, "none")
  seeds[nm] <- SEED + utf8ToInt(nm)
  set.seed(seeds[[nm]])
  f1 <- fit_efa(dat, m, "geomin",
                rotation.args = list(geomin_epsilon = eps, rstarts = 30L,
                                     row_weights = "none", std_ov = TRUE,
                                     gpa_tol = 1e-6, order_lv_by = "index",
                                     orthogonal = FALSE))
  stopifnot(lavInspect(f0, "converged"), lavInspect(f1, "converged"))
  ra <- lavInspect(f1, "options")$rotation.args
  stopifnot(ra$geomin_epsilon == eps, ra$rstarts == 30, ra$row_weights == "none",
            isTRUE(ra$std_ov), ra$gpa_tol == 1e-6, ra$order_lv_by == "index",
            isFALSE(ra$orthogonal))

  e0 <- lavInspect(f0, "est"); e1 <- lavInspect(f1, "est")
  L  <- e1$lambda[yn, , drop = FALSE]; P <- e1$psi
  stopifnot(max(abs(diag(P) - 1)) < 1e-12,
            all(colSums(L) > 0),                    # lavaan's reflection held
            max(abs(L %*% P %*% t(L) -              # GM1: sigma unchanged
                    e0$lambda[yn, ] %*% e0$psi %*% t(e0$lambda[yn, ]))) < 1e-12)

  ss <- f1@test$scaled.shifted
  kl <- which(lower.tri(P), arr.ind = TRUE)
  fn <- colnames(L)
  rows[[nm]] <- rbind(
    data.frame(case = nm, quantity = "test", label = c("T", "df", "T3"),
               value = c(f1@test$standard$stat, ss$df, ss$stat)),
    data.frame(case = nm, quantity = "fit",
               label = c("rmsea.scaled", "cfi.scaled", "tli.scaled"),
               value = unname(fitMeasures(f1, c("rmsea.scaled", "cfi.scaled",
                                                "tli.scaled")))),
    data.frame(case = nm, quantity = "rot", label = c("eps", "Q"),
               value = c(eps, geomin_Q(L, eps))),
    data.frame(case = nm, quantity = "lambda",
               label = paste0(rep(fn, each = length(yn)), "=~", yn),
               value = as.vector(L)),
    data.frame(case = nm, quantity = "phi",
               label = paste0(fn[kl[, 1]], "~~", fn[kl[, 2]]),
               value = P[lower.tri(P)]))
  cat(sprintf("case %s: m=%d eps=%g  T3 %.4f df %d  Q %.10f  converged rotation\n",
              nm, m, eps, ss$stat, ss$df, geomin_Q(L, eps)))
}

out <- do.call(rbind, rows); rownames(out) <- NULL
stopifnot(!anyNA(out$value))
f <- file.path(OUT, "lavaan_efa_geomin.csv")
write_csv_lf(out, f)

write_manifest(MAN, f,
  seeds = seeds,
  extra = list(packages = c(
    paste("lavaan", as.character(packageVersion("lavaan"))),
    "estimator WLSMV, parameterization delta, test scaled.shifted",
    "zero.add c(0.5, 0), zero.keep.margins TRUE (as 99g_)",
    "rotation geomin, oblique; geomin_epsilon by Mplus's rule (E .0001, F .001, G .0001)",
    "rstarts 30 (seeded per case), row_weights none, std_ov TRUE, gpa_tol 1e-6, order_lv_by index",
    "inputs: case_E.csv (10_fixtures.R), case_F.csv and case_G.csv (91_efa_data.R)",
    R.version.string)))

cat("lavaan_efa_geomin.csv and", MAN, "written.\n")
