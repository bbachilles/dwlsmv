# dwlsmv -- D7 probe: lavaan's std.all equals std.lv under the delta
# parameterization with ordered indicators (Case E, F2 <- F1; marker and
# unit-variance metrics), and an endogenous factor's total variance under
# unit variance is not 1. A probe, not a fixture: nothing reads its output.
# Executed 2026-10-04 under R 4.6.0, lavaan 0.7.2 (tools/probe_d7_stdall.out).
# Run from the repository root:  Rscript tools/probe_d7_stdall.R
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek
suppressPackageStartupMessages(library(lavaan))
cat("lavaan", as.character(packageVersion("lavaan")), "R", R.version.string, "\n")
d <- read.csv("fixtures/data/case_E.csv")
ord <- paste0("y", 1:6)
models <- list(
  marker = "F1 =~ y1 + y2 + y3\nF2 =~ y4 + y5 + y6\nF2 ~ F1",
  unitvar = "F1 =~ NA*y1 + y2 + y3\nF2 =~ NA*y4 + y5 + y6\nF2 ~ F1\nF1 ~~ 1*F1\nF2 ~~ 1*F2"
)
for (m in names(models)) {
  fit <- sem(models[[m]], data = d, ordered = ord, estimator = "WLSMV",
             parameterization = "delta")
  a <- standardizedSolution(fit, type = "std.all")
  l <- standardizedSolution(fit, type = "std.lv")
  keep <- a$op %in% c("=~", "~", "~~")
  de <- max(abs(a$est.std[keep] - l$est.std[keep]))
  ds <- max(abs(a$se[keep] - l$se[keep]), na.rm = TRUE)
  cat(sprintf("STD %-8s rows %d  max|est.all - est.lv| %.3e  max|se.all - se.lv| %.3e\n",
              m, sum(keep), de, ds))
  # where do they differ, if anywhere
  w <- which(keep & abs(a$est.std - l$est.std) > 1e-12)
  if (length(w)) print(cbind(a[w, 1:3], all = a$est.std[w], lv = l$est.std[w]))
  # the endogenous factor's total variance under each metric (Phi_22)
  cat(sprintf("STD %-8s Phi_F2F2 (model-implied latent variance) %.6f\n", m,
              lavInspect(fit, "cov.lv")["F2", "F2"]))
}
