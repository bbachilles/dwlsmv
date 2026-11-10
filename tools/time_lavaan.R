# ===========================================================================
# Drafted 2026-09-28 by a read-only agent. Executed unchanged by the main
# session, 2026-09-28 (R 4.6.0, lavaan 0.7.2).
# Results: E-089.
# ===========================================================================
#
# dwlsmv -- spec 7.5's timing: lavaan as a TOTAL, public interface only
#
# Binding protocol: E-008, 2026-09-28, "Spec 7.5's timing
# protocol, agreed before any driver." A diagnostic, not a test.
#
#   cd <repo root>
#   "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/time_compare.do writedata   (once)
#   Rscript tools/time_lavaan.R > time_lavaan.log 2>&1
#   grep "^TIME|" time_lavaan.log
#
# RUN IT ON AN IDLE MACHINE.
#
# DATA. The CSVs tools/time_compare.do writes to fixtures/replication/timing/
# (gitignored): the 12 grid cells simulated by time_baseline.do's own code
# and seed, and replicate 1 of mc_d, mc_a1p_16, mc_b1_m1, mc_c_998, mc_B_9975
# copied from their .dta. No simulation here, so every tool sees identical data.
#
# WHAT IS TIMED. One lavaan::cfa() call per fit (grid) or lavaan::sem() call
# (MC rows, as 70_ does), wall clock (proc.time()["elapsed"]), median of 5
# repetitions (3 at p 40, N 10000), after ONE untimed warm-up fit per design.
# Excluded: package load, CSV read, the construction of the ordered
# data.frame, the warm-up, and anything after the call returns (convergence
# and iteration reads are outside the timer).
#
# OPTIONS. 70_lavaan_mc_reference.R's pinned "ref" arm, unchanged:
#   estimator "WLSMV", parameterization "delta", test "scaled.shifted",
#   zero.add c(0, 0), zero.keep.margins TRUE, bad.marker.crit 0
# plus baseline = FALSE. Verified against the installed lavaan 0.7.2
# (deparsed from the namespace, 2026-09-28):
#   - lav_options_default(): elm("baseline", TRUE) -- the default fits it.
#   - lav_step15_baseline(): the independence model is fitted only if
#     lavoptions$baseline is TRUE (and do.fit, and test != "none"); with
#     baseline = FALSE it returns list() without fitting. A toy fit confirmed
#     lavInspect(fit, "options")$baseline FALSE and an empty @baseline slot.
#   - Fit measures: lavaan() never calls fitMeasures(); they are computed only
#     by fitMeasures() / summary(fit.measures = TRUE), which are not called.
#     lav_step14_fit() only packs the already-computed test into @Fit.
#   - lav_step12_loglik(): lav_model_loglik() computes a loglik only for
#     estimator ML / MML; under WLSMV (DWLS) it returns NA without work, so
#     loglik is left at its default.
#   - h1 is left at TRUE: lav_step06_h1() -> lav_h1_implied_logl() builds the
#     saturated (polychoric) moments the WLSMV test is computed against, which
#     is work dwlsmv also does. Switching it off is not a documented way to
#     get the scaled.shifted test, so it is not touched.
# Grid model: the true structure, std.lv = TRUE (unit variance), factor
# correlations free (cfa's default), thresholds free: dwlsmv's grid spec.
# MC model: "f =~ y1 + ... + yp", marker (first loading 1), as 70_ and
# tools/paired_mc.do.
#
# WHERE THE TOTAL CANNOT BE MATCHED EXACTLY (state in the performance note)
#   - lavaan's total includes lavData(): checking and converting the
#     data.frame to its internal form. dwlsmv's starts from Y in Mata. There
#     is no public boundary to exclude it.
#   - lavaan also builds the parameter table from model syntax (dwlsmv's
#     dwls_identify() parses its canonical spec: comparable, inside both).
#   - lavaan computes vcov and the test by default (so does dwlsmv's Stage 2);
#     it also runs post-fit checks (check.post, check.vcov), which cannot be
#     turned off without changing what is checked. Left at defaults.
#   - R is largely single-threaded; matrix products use the BLAS recorded
#     below (reference BLAS unless one is swapped in).
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

suppressPackageStartupMessages(library(lavaan))

TDIR <- file.path("fixtures", "replication", "timing")
if (!dir.exists(TDIR)) stop("run tools/time_compare.do writedata first: ", TDIR)

si <- sessionInfo()
cat(sprintf("TIME|machine|%s|lavaan %s|BLAS %s|LAPACK %s (%s)|%s|processors %s logical (%s), detectCores %d|%s\n",
            R.version.string, as.character(packageVersion("lavaan")),
            if (nzchar(si$BLAS)) si$BLAS else "R default (reference)",
            La_version(), if (nzchar(La_library())) La_library() else "R default",
            si$running, Sys.getenv("NUMBER_OF_PROCESSORS"),
            Sys.getenv("PROCESSOR_IDENTIFIER"), parallel::detectCores(),
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("TIME|machine|lavaan library %s\n", find.package("lavaan")))

opts <- list(estimator = "WLSMV", parameterization = "delta",
             test = "scaled.shifted", zero.add = c(0, 0),
             zero.keep.margins = TRUE, bad.marker.crit = 0,
             baseline = FALSE)

read_ordered <- function(f, K) {
  Y <- as.matrix(read.csv(f, colClasses = "integer"))
  if (ncol(Y) != length(K)) stop(f, ": ", ncol(Y), " columns, expected ", length(K))
  d <- as.data.frame(lapply(seq_len(ncol(Y)), function(j)
    ordered(Y[, j], levels = 0:(K[j] - 1))))
  names(d) <- sprintf("y%d", seq_len(ncol(Y)))
  if (anyNA(d)) stop(f, ": a value outside 0..K-1")
  d
}

grid_model <- function(p) {
  nf <- p / 10
  paste(vapply(seq_len(nf), function(f)
    sprintf("f%d =~ %s", f, paste(sprintf("y%d", (f - 1) * 10 + 1:10), collapse = " + ")),
    ""), collapse = "\n")
}

median_time <- function(call_fit, reps) {
  fit <- call_fit()                                   # warm-up, untimed
  t <- numeric(reps)
  for (r in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    fit <- call_fit()
    t[r] <- proc.time()[["elapsed"]] - t0
  }
  list(t = median(t), fit = fit)
}

report <- function(tag, d, reps, res) {
  fit <- res$fit
  cat(sprintf("TIME|lavaan|%s|p %d K %d N %d|npar %d|reps %d|total %9.4f s|converged %d|iterations %d\n",
              tag, ncol(d), max(vapply(d, nlevels, 1L)), nrow(d),
              lavInspect(fit, "npar"), reps, res$t,
              as.integer(lavInspect(fit, "converged")),
              as.integer(lavInspect(fit, "iterations"))))
  flush.console()
}

# timer resolution: proc.time() elapsed on Windows is ms; the smallest cells
# are ~0.1 s, so reported to 4 decimals but read to 3.
for (p in c(10, 20, 40)) for (K in c(3, 5)) for (N in c(1000, 10000)) {
  d <- read_ordered(file.path(TDIR, sprintf("grid_p%d_K%d_N%d.csv", p, K, N)), rep(K, p))
  model <- grid_model(p)
  reps <- if (p == 40 && N == 10000) 3L else 5L
  res <- median_time(function() do.call(lavaan::cfa,
           c(list(model = model, data = d, ordered = names(d), std.lv = TRUE), opts)), reps)
  report("grid", d, reps, res)
}

mcK <- list(mc_d = c(4, 4, 3, 3), mc_a1p_16 = rep(3, 16), mc_b1_m1 = c(4, 4, 3, 3),
            mc_c_998 = c(4, 4, 3, 3), mc_B_9975 = c(4, 4, 3, 3))
for (nm in names(mcK)) {
  d <- read_ordered(file.path(TDIR, paste0(nm, ".csv")), mcK[[nm]])
  model <- paste("f =~", paste(names(d), collapse = " + "))
  res <- median_time(function() do.call(lavaan::sem,
           c(list(model = model, data = d, ordered = names(d)), opts)), 5L)
  pt <- parTable(res$fit)
  mk <- pt$op == "=~" & pt$rhs == "y1"
  if (!(sum(mk) == 1L && pt$free[mk] == 0L && pt$est[mk] == 1))
    stop("marker identification did not hold on ", nm)
  report(nm, d, 5L, res)
}

cat("TIME|done\n")
