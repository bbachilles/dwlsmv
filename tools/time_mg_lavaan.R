# dwlsmv -- spec 7.5's timing at G > 1: lavaan as a total (E-089)
#
# Executed 2026-11-02 (E-089). A diagnostic, not a test.
#
#   Rscript tools/time_mg_lavaan.R match    untimed fits, compared with tools/time_mg.do match
#   Rscript tools/time_mg_lavaan.R time     median of 5 after one warm-up
#   Rscript tools/time_mg_lavaan.R smoke    time, smallest dataset only
#
# Options: 70_'s pinned arm plus baseline = FALSE, as tools/time_lavaan.R. Level syntax:
# 98_lavaan_multigroup.R's table, written out for several factors and groups.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

suppressPackageStartupMessages(library(lavaan))
mode <- commandArgs(trailingOnly = TRUE)[1]
TDIR <- file.path("fixtures", "replication", "timing")
opts <- list(estimator = "WLSMV", parameterization = "delta", test = "scaled.shifted",
             zero.add = c(0, 0), zero.keep.margins = TRUE, bad.marker.crit = 0, baseline = FALSE)
cat(sprintf("TIME|machine|%s|lavaan %s|%s\n", R.version.string, packageVersion("lavaan"),
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cvec <- function(first, rest, G) sprintf("c(%s)", paste(c(first, rep(rest, G - 1L)), collapse = ","))

syntax_of <- function(G, p, K, level) {
  sc <- level == "scalar"; out <- character()
  for (f in seq_len(p / 10)) {
    v <- sprintf("y%d", (f - 1) * 10 + 1:10); j <- (f - 1) * 10 + 1:10
    ld <- if (sc) vapply(j, function(i) sprintf("%s*y%d", cvec(paste0("l", i), paste0("l", i), G), i), "")
          else sprintf("%s*%s", cvec("NA", "NA", G), v)
    out <- c(out, paste(sprintf("F%d =~", f), paste(ld, collapse = " + ")),
             sprintf("F%d ~~ %s*F%d", f, cvec("1", if (sc) "NA" else "1", G), f),
             sprintf("F%d ~ %s*1", f, cvec("0", if (sc) "NA" else "0", G)))
  }
  for (j in seq_len(p)) {
    if (sc) out <- c(out, sprintf("y%d | %s", j, paste(vapply(seq_len(K - 1L), function(k)
      sprintf("%s*t%d", cvec(sprintf("t%d_%d", j, k), sprintf("t%d_%d", j, k), G), k), ""), collapse = " + ")))
    out <- c(out, sprintf("y%d ~*~ %s*y%d", j, cvec("1", if (sc) "NA" else "1", G), j),
             sprintf("y%d ~ %s*1", j, cvec("0", "0", G)))
  }
  paste(out, collapse = "\n")
}

read_mg <- function(G, p, K, N) {
  d <- read.csv(file.path(TDIR, sprintf("mg_G%d_p%d_K%d_N%d.csv", G, p, K, N)), colClasses = "integer")
  y <- sprintf("y%d", seq_len(p))
  d[y] <- lapply(d[y], ordered, levels = 0:(K - 1L))
  d
}

fit_of <- function(d, G, p, K, level)
  do.call(lavaan::cfa, c(list(model = syntax_of(G, p, K, level), data = d, group = "g",
                              ordered = sprintf("y%d", seq_len(p)), std.lv = TRUE), opts))

if (mode == "match") {
  m  <- read.csv(file.path(TDIR, "match_dwlsmv.csv"), stringsAsFactors = FALSE)
  th <- read.csv(file.path(TDIR, "match_theta.csv"), stringsAsFactors = FALSE)
}
for (G in c(2L, 3L)) for (p in c(10L, 20L)) for (K in c(3L, 5L)) for (N in c(500L, 5000L)) {
  if (mode == "smoke" && !(G == 2L && p == 10L && K == 3L && N == 500L)) next
  d <- read_mg(G, p, K, N); cell <- sprintf("G%d_p%d_K%d_N%d", G, p, K, N)
  for (level in c("configural", "scalar")) {
    if (mode == "match") {
      f <- fit_of(d, G, p, K, level)
      r <- m[m$cell == cell & m$level == level, ]
      pe <- parameterEstimates(f); pe <- pe[pe$op %in% c("=~", "|") & pe$se > 0, ]
      ours <- th$abs_theta[th$cell == cell & th$level == level]
      gap <- max(vapply(abs(pe$est), function(x) min(abs(ours - x)), 0))
      fm <- lavInspect(f, "test")$scaled.shifted
      cat(sprintf("TIME|match|lavaan|%s|%s|npar %d vs %d|df %g vs %g|T3 %.6f vs %.6f rel %.1e|lambda,tau gap %.1e|converged %d\n",
                  cell, level, lavInspect(f, "npar"), r$npar, fm$df, r$df, fm$stat, r$T3,
                  abs(fm$stat - r$T3) / r$T3, gap, as.integer(lavInspect(f, "converged"))))
    } else {
      f <- fit_of(d, G, p, K, level)
      t <- numeric(5)
      for (i in 1:5) {
        t0 <- proc.time()[["elapsed"]]
        f <- fit_of(d, G, p, K, level)
        t[i] <- proc.time()[["elapsed"]] - t0
      }
      cat(sprintf("TIME|lavaan|%s|%s|npar %d|reps 5|total %9.4f s|converged %d|iterations %d\n", cell, level,
                  lavInspect(f, "npar"), median(t), as.integer(lavInspect(f, "converged")),
                  as.integer(lavInspect(f, "iterations"))))
    }
    flush.console()
  }
}
cat("TIME|done\n")
