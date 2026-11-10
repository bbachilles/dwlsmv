# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 2: lavaan's full conditional fit
# (conditional.x = TRUE), one binary x, single group. Item 4b of the
# conditional-covariate work order (E-006, 2026-10-24); the
# comparison is @sec-condx-full.
#
# Written 2026-10-25 by a side session of the main session.
#
# Data: fixtures/data/condx_A.csv, committed by 99h_lavaan_condx.R. It is
# read here, never regenerated; its md5 is recorded in this manifest too.
#
# Model: f =~ y1 + y2 + y3; f ~ x. WLSMV, delta, conditional.x = TRUE
# (explicit), std.lv = TRUE (fixes f's disturbance variance psi = 1, f being
# endogenous), zero.add = c(0.5, 0), zero.keep.margins = TRUE. Written:
#   lavaan_condx_fit_est.csv   lhs, op, rhs, free, est, se  (parameterEstimates)
#   lavaan_condx_fit_test.csv  quantity, value
#     N, chisq (standard), df, chisq.scaled (scaled.shifted), its df,
#     scaling.factor (= 1/a) and shift.parameter (= b), fmin.
# The robust test's type is asserted to be "scaled.shifted" before writing.
#
# Run from the repository root:
#   Rscript fixtures/generate/99i_lavaan_condx_fit.R
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))
if (!requireNamespace("lavaan", quietly = TRUE)) stop("lavaan is required for tier-2 fixtures.")
library(lavaan)
stopifnot(packageVersion("lavaan") == "0.7.2")

DAT <- file.path("fixtures", "data", "condx_A.csv")
OUT <- file.path("fixtures", "expected")
MAN <- file.path("fixtures", "MANIFEST_lavaan_condx_fit.json")

dat <- read.csv(DAT, colClasses = "character")
dat[] <- lapply(dat, as.integer)
stopifnot(identical(names(dat), c("y1", "y2", "y3", "x")), nrow(dat) == 1500L)

fit <- lavaan::sem("f =~ y1 + y2 + y3\n f ~ x", data = dat,
                   ordered = c("y1", "y2", "y3"), estimator = "WLSMV",
                   parameterization = "delta", conditional.x = TRUE, std.lv = TRUE,
                   zero.add = c(0.5, 0.0), zero.keep.margins = TRUE)
stopifnot(isTRUE(lavInspect(fit, "options")$conditional.x),
          lavInspect(fit, "converged"))

pe <- parameterEstimates(fit, ci = FALSE, zstat = FALSE, pvalue = FALSE,
                         remove.nonfree = FALSE)
pt <- parTable(fit)
free <- pt$free[match(paste(pe$lhs, pe$op, pe$rhs), paste(pt$lhs, pt$op, pt$rhs))]
f_e <- file.path(OUT, "lavaan_condx_fit_est.csv")
write_csv_lf(data.frame(lhs = pe$lhs, op = pe$op, rhs = pe$rhs, free = free,
                        est = sprintf("%.17g", pe$est), se = sprintf("%.17g", pe$se)), f_e)

te <- lavInspect(fit, "test")
stopifnot(names(te)[2] == "scaled.shifted")
st <- te$standard; ss <- te$scaled.shifted
tq <- c(N = lavInspect(fit, "nobs"), chisq = st$stat, df = st$df,
        chisq.scaled = ss$stat, df.scaled = ss$df,
        scaling.factor = ss$scaling.factor, shift.parameter = ss$shift.parameter,
        fmin = unname(fitMeasures(fit, "fmin")))
f_t <- file.path(OUT, "lavaan_condx_fit_test.csv")
write_csv_lf(data.frame(quantity = names(tq), value = sprintf("%.17g", unname(tq))), f_t)
print(tq)

write_manifest(MAN, c(DAT, f_e, f_t), seeds = "none (reads committed condx_A.csv)",
  extra = list(packages = c(paste("lavaan", as.character(packageVersion("lavaan"))),
    "sem('f =~ y1 + y2 + y3; f ~ x'), WLSMV, delta, conditional.x = TRUE, std.lv, zero.add = c(0.5, 0), zero.keep.margins = TRUE",
    R.version.string)))
cat(sprintf("%s, %s and %s written\n", f_e, f_t, MAN))
