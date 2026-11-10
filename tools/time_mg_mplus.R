# dwlsmv -- spec 7.5's timing at G > 1: Mplus as a total (E-089)
#
# Executed 2026-11-02 (E-089). A diagnostic, not a test.
#
#   Rscript tools/time_mg_mplus.R match    untimed runs, compared with tools/time_mg.do match
#   Rscript tools/time_mg_mplus.R time     median of 5 process runs after one, PROCESSORS 1 and 6
#   Rscript tools/time_mg_mplus.R smoke    time, smallest dataset only
#
# Each level is stated explicitly, not left to Mplus's GROUPING defaults. Mplus's statistic is
# compared with ours as a N F + b (spec 4.9), since Mplus uses N where we use N - G.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

mode <- commandArgs(trailingOnly = TRUE)[1]
TDIR <- file.path("fixtures", "replication", "timing")
mplus_exe <- "C:/Program Files/Mplus/Mplus.exe"
work <- file.path(tempdir(), "dwlsmv_time_mg_mplus")
dir.create(work, showWarnings = FALSE, recursive = TRUE)

run_once <- function(inp) {
  owd <- setwd(work); on.exit(setwd(owd))
  t0 <- proc.time()[["elapsed"]]
  processx::run(mplus_exe, inp, error_on_status = FALSE)
  proc.time()[["elapsed"]] - t0
}
wrap <- function(words, lead = "    ") {
  lines <- character(0); cur <- lead
  for (w in words) {
    if (nchar(cur) + nchar(w) + 1 > 80) { lines <- c(lines, cur); cur <- lead }
    cur <- paste0(cur, if (cur == lead) "" else " ", w)
  }
  c(lines, cur)
}

model_of <- function(G, p, K, level) {
  nf <- p / 10
  by <- unlist(lapply(seq_len(nf), function(f) {
    j <- (f - 1) * 10 + 1:10; v <- sprintf("y%d", j); v[1] <- paste0(v[1], "*")
    c(sprintf("  F%d BY", f), wrap(v), "  ;")
  }))
  fs <- sprintf("F%d", seq_len(nf)); ys <- sprintf("y%d", seq_len(p))
  overall <- c("MODEL:", by, sprintf("  %s@1;", fs), sprintf("  [%s@0];", fs),
               "  {", wrap(paste0(ys, "@1")), "  };")
  for (g in 2:G) {
    gl <- if (level == "scalar") {
      c(sprintf("  %s*;", fs), sprintf("  [%s*];", fs), "  {", wrap(paste0(ys, "*")), "  };")
    } else {
      thr <- as.vector(outer(ys, seq_len(K - 1L), function(a, b) sprintf("%s$%d", a, b)))
      c(by, sprintf("  %s@1;", fs), sprintf("  [%s@0];", fs), "  [", wrap(thr), "  ];")
    }
    overall <- c(overall, sprintf("MODEL g%d:", g), gl)
  }
  overall
}

input_of <- function(stem, G, p, K, level, np) {
  ys <- sprintf("y%d", seq_len(p))
  c(sprintf("TITLE: dwlsmv timing -- %s %s", stem, level), sprintf("DATA: FILE IS %s.dat;", stem),
    "VARIABLE:", "  NAMES ARE g", wrap(ys), "  ;", "  USEVARIABLES ARE", wrap(ys), "  ;",
    "  CATEGORICAL ARE", wrap(ys), "  ;",
    sprintf("  GROUPING IS g (%s);", paste(sprintf("%d = g%d", seq_len(G), seq_len(G)), collapse = " ")),
    "ANALYSIS:", "  ESTIMATOR = WLSMV;", "  PARAMETERIZATION = DELTA;", sprintf("  PROCESSORS = %d;", np),
    model_of(G, p, K, level))
}

read_out <- function(lines) {
  i <- grep("Chi-Square Test of Model Fit", lines)[1]
  num <- function(s) as.numeric(sub("[*]", "", tail(strsplit(trimws(s), " +")[[1]], 1)))
  chi <- num(lines[i + 2]); df <- num(lines[i + 3])
  npar <- num(grep("Number of Free Parameters", lines, value = TRUE)[1])
  st <- grep("^MODEL RESULTS", lines)[1]
  en <- st + grep("^(R-SQUARE|QUALITY OF NUMERICAL|TECHNICAL|Beginning Time)", lines[-(1:st)])[1]
  sect <- ""; est <- numeric(0)
  for (l in lines[st:en]) {
    t <- trimws(l)
    if (grepl(" BY$", t)) { sect <- "keep"; next }
    if (t == "Thresholds") { sect <- "keep"; next }
    if (grepl("^(Means|Variances|Scales| *[A-Z0-9]+ +WITH|Group )", t) || grepl(" WITH$", t)) { sect <- ""; next }
    w <- strsplit(t, " +")[[1]]
    if (sect == "keep" && length(w) == 5L && suppressWarnings(!is.na(as.numeric(w[3]))) && as.numeric(w[3]) > 0)
      est <- c(est, abs(as.numeric(w[2])))
  }
  list(chi = chi, df = df, npar = npar, est = est,
       ok = any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines, fixed = TRUE)))
}

startup <- NA_real_
if (mode != "match") {
  set.seed(1)
  write.table(matrix(round(rnorm(20), 3), 10, 2), file.path(work, "trivial.dat"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  writeLines(c("TITLE: trivial", "DATA: FILE IS trivial.dat;", "VARIABLE: NAMES ARE u1 u2;",
               "ANALYSIS: TYPE = BASIC;"), file.path(work, "trivial.inp"))
  run_once("trivial.inp")
  startup <- median(vapply(1:5, function(i) run_once("trivial.inp"), 0))
  cat(sprintf("TIME|mplus|startup|reps 5|total %9.4f s\n", startup))
} else {
  m  <- read.csv(file.path(TDIR, "match_dwlsmv.csv"), stringsAsFactors = FALSE)
  th <- read.csv(file.path(TDIR, "match_theta.csv"), stringsAsFactors = FALSE)
}

for (G in c(2L, 3L)) for (p in c(10L, 20L)) for (K in c(3L, 5L)) for (N in c(500L, 5000L)) {
  if (mode == "smoke" && !(G == 2L && p == 10L && K == 3L && N == 500L)) next
  cell <- sprintf("G%d_p%d_K%d_N%d", G, p, K, N); stem <- tolower(gsub("_", "", cell))
  Y <- as.matrix(read.csv(file.path(TDIR, sprintf("mg_%s.csv", cell)), colClasses = "integer"))
  write.table(Y, file.path(work, paste0(stem, ".dat")), row.names = FALSE, col.names = FALSE, quote = FALSE)
  for (level in c("configural", "scalar")) for (np in if (mode == "match") 1L else c(1L, 6L)) {
    inp <- sprintf("%s%s%d.inp", stem, substr(level, 1, 1), np)
    writeLines(input_of(stem, G, p, K, level, np), file.path(work, inp))
    if (mode == "match") {
      run_once(inp)
      o <- read_out(readLines(file.path(work, sub("[.]inp$", ".out", inp)), warn = FALSE))
      r <- m[m$cell == cell & m$level == level, ]
      ours <- th$abs_theta[th$cell == cell & th$level == level]
      gap <- max(vapply(o$est, function(x) min(abs(ours - x)), 0))
      cat(sprintf("TIME|match|mplus|%s|%s|npar %g vs %g|df %g vs %g|chi %.3f vs aNF+b %.3f rel %.1e|lambda,tau gap %.1e on %d|normal %d\n",
                  cell, level, o$npar, r$npar, o$df, r$df, o$chi, r$Tmplus, abs(o$chi - r$Tmplus) / r$Tmplus,
                  gap, length(o$est), as.integer(o$ok)))
    } else {
      run_once(inp)
      t <- median(vapply(1:5, function(i) run_once(inp), 0))
      o <- read_out(readLines(file.path(work, sub("[.]inp$", ".out", inp)), warn = FALSE))
      cat(sprintf("TIME|mplus|%s|%s|PROCESSORS %d|reps 5|total %9.4f s|startup (beside) %9.4f|normal %d\n",
                  cell, level, np, t, startup, as.integer(o$ok)))
    }
    flush.console()
  }
}
cat("TIME|done\n")
