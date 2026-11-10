# ===========================================================================
# Drafted 2026-09-28 by a read-only agent. Executed by the main session,
# 2026-09-28 (Mplus 8.10). The first run timed with system2(), whose wait on
# Windows resolves only to ~0.1 s; the clock is now processx::run() and the
# reported results are from the re-run.
# Results: E-089.
# ===========================================================================
#
# dwlsmv -- spec 7.5's timing: Mplus as a TOTAL, one process run per fit
#
# Binding protocol: E-008, 2026-09-28, "Spec 7.5's timing
# protocol, agreed before any driver." A diagnostic, not a test.
#
#   cd <repo root>
#   "/c/Program Files/Stata18/StataMP-64.exe" //e do tools/time_compare.do writedata   (once)
#   Rscript tools/time_mplus.R > time_mplus.log 2>&1
#   grep "^TIME|" time_mplus.log
#
# RUN IT ON AN IDLE MACHINE.
#
# DATA. The CSVs tools/time_compare.do writes to fixtures/replication/timing/:
# the 12 grid cells and replicate 1 of mc_d, mc_a1p_16, mc_b1_m1, mc_c_998,
# mc_B_9975. Each is rewritten, UNTIMED, as a headerless space-delimited .dat
# in a flat scratch directory (30_mplus_reference.R's convention: Mplus
# resolves FILE IS relative to the input and dislikes long paths).
#
# WHAT IS TIMED. Wall clock of one Mplus.exe process (processx::run, waited on),
# median over 5 process runs per design and PROCESSORS setting. Each process
# is a cold start by construction, which is why start-up is reported
# separately.
# START-UP: the median of 5 runs of a trivial input (two continuous
# variables, ten rows, TYPE = BASIC) is reported BESIDE each total and never
# subtracted.
# One untimed run of each input precedes its 5 timed runs, so the OS file
# cache holds Mplus's binaries and the data file for every timed run.
#
# MODEL. Grid: true structure, factor variances @1 and every loading free
# (Fk BY y1* ...), factor covariances free (Mplus's default for exogenous
# factors), ESTIMATOR = WLSMV, PARAMETERIZATION = DELTA. MC rows: marker
# model, f BY y1@1 y2 ...; f; exactly 30_'s setup.
# PROCESSORS = 1 and 6 in ANALYSIS. Note, from the Mplus 8 User's Guide
# (ANALYSIS command, PROCESSORS): multiple processors are documented for
# mixture, Bayes with several chains, numerical integration, continuous ML
# with missing data, and TYPE = TWOLEVEL with categorical WLSMV. Single-level
# categorical WLSMV is NOT in that list, so PROCESSORS = 6 may change nothing
# here; the two rows measure whether it does.
#
# NO TYPE = BASIC ROW. The User's Guide (Chapter 16, ANALYSIS TYPE) says only
# that BASIC "computes sample statistics and other descriptive information";
# the SAVEDATA SAMPLE option saves correlations/covariances; the only
# documented saving of a sample-statistics asymptotic covariance matrix is
# SWMATRIX, for TYPE = TWOLEVEL. Nothing documents that single-level BASIC
# computes the WLSMV weight matrix, so by the agreed rule no BASIC row appears.
# (Searched: Mplus User's Guide Version 8.pdf in C:/Program Files/Mplus/
# Documentation. The 8.1 and 8.5 Language Addenda were not searched.)
#
# CONVERGENCE AND ITERATIONS. Convergence is read from the .out ("THE MODEL
# ESTIMATION TERMINATED NORMALLY"). Mplus prints no iteration count by
# default; TECH5 would, but adds output work to the timed process, so the
# iteration count is NOT recorded. Mplus version: from the .out header.
#
# WHERE THE TOTAL CANNOT BE MATCHED EXACTLY (state in the performance note)
#   - Process start-up, licence check, parsing the .inp, and READING THE DATA
#     FILE happen inside the process and cannot be excluded (the decisions
#     row says so). Start-up is shown beside, never subtracted.
#   - Writing the .out file happens inside the process.
#   - Mplus computes its default output unconditionally: the baseline
#     (independence) model chi-square and fit measures (CFI/TLI need the
#     baseline). There is no documented switch to skip them, so Mplus does
#     work dwlsmv and the timed lavaan call do not.
#   - Mplus's sample statistics for WLSMV are its own implementation; it is
#     the same estimator, not the same code path.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

TDIR <- file.path("fixtures", "replication", "timing")
if (!dir.exists(TDIR)) stop("run tools/time_compare.do writedata first: ", TDIR)

mplus_exe <- Sys.which("Mplus")
if (!nzchar(mplus_exe)) {
  cand <- "C:/Program Files/Mplus/Mplus.exe"
  if (file.exists(cand)) mplus_exe <- cand
}
if (!nzchar(mplus_exe) || !file.exists(mplus_exe)) stop("Mplus executable not found.")

work <- file.path(tempdir(), "dwlsmv_time_mplus")
dir.create(work, showWarnings = FALSE, recursive = TRUE)

NREP  <- 5L
PROCS <- c(1L, 6L)

run_once <- function(inp) {
  owd <- setwd(work); on.exit(setwd(owd))
  # processx, not system2: on Windows system2's wait resolves only to ~0.1 s
  # (a no-op cmd.exe reads 0.11 s under system2, 0.01-0.04 s under processx;
  # measured 2026-09-28, after the first run came back in 0.11 s steps).
  t0 <- proc.time()[["elapsed"]]
  processx::run(mplus_exe, inp, error_on_status = FALSE)
  proc.time()[["elapsed"]] - t0
}

out_of <- function(inp) file.path(work, sub("[.]inp$", ".out", tolower(inp)))

time_input <- function(inp) {
  run_once(inp)                                       # untimed: warms the file cache
  t <- vapply(seq_len(NREP), function(r) run_once(inp), 0)
  out <- readLines(out_of(inp), warn = FALSE)
  list(t = median(t), out = out,
       converged = as.integer(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", out, fixed = TRUE))))
}

# --- start-up: a trivial run --------------------------------------------------
set.seed(1)
write.table(matrix(round(rnorm(20), 3), 10, 2), file.path(work, "trivial.dat"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
writeLines(c("TITLE: dwlsmv timing -- trivial run (start-up)",
             "DATA: FILE IS trivial.dat;",
             "VARIABLE: NAMES ARE u1 u2;",
             "ANALYSIS: TYPE = BASIC;"), file.path(work, "trivial.inp"))
triv <- time_input("trivial.inp")
ver  <- grep("Mplus VERSION", triv$out, value = TRUE)[1]
cat(sprintf("TIME|machine|%s|%s|processors %s logical (%s)|R %s|%s\n",
            trimws(ver), mplus_exe, Sys.getenv("NUMBER_OF_PROCESSORS"),
            Sys.getenv("PROCESSOR_IDENTIFIER"), R.version.string,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("TIME|mplus|startup|trivial TYPE = BASIC run|reps %d|total %9.4f s\n", NREP, triv$t))

# --- one design ----------------------------------------------------------------
wrap <- function(words, lead = "  ") {                 # Mplus lines must stay < 90 chars
  lines <- character(0); cur <- lead
  for (w in words) {
    if (nchar(cur) + nchar(w) + 1 > 80) { lines <- c(lines, cur); cur <- lead }
    cur <- paste0(cur, if (cur == lead) "" else " ", w)
  }
  c(lines, cur)
}

design <- function(tag, csv, model_lines) {
  Y <- as.matrix(read.csv(csv, colClasses = "integer"))
  p <- ncol(Y)
  stem <- tolower(gsub("[^A-Za-z0-9]", "", sub("[.]csv$", "", basename(csv))))
  write.table(Y, file.path(work, paste0(stem, ".dat")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  vars <- sprintf("y%d", seq_len(p))
  for (np in PROCS) {
    inp <- sprintf("%s_p%d.inp", stem, np)
    writeLines(c(
      sprintf("TITLE: dwlsmv timing -- %s", tag),
      sprintf("DATA: FILE IS %s.dat;", stem),
      "VARIABLE:", "  NAMES ARE", wrap(vars, "    "), "  ;",
      "  CATEGORICAL ARE", wrap(vars, "    "), "  ;",
      "ANALYSIS:", "  ESTIMATOR = WLSMV;", "  PARAMETERIZATION = DELTA;",
      sprintf("  PROCESSORS = %d;", np),
      "MODEL:", model_lines), file.path(work, inp))
    res <- time_input(inp)
    cat(sprintf("TIME|mplus|%s|p %d K %d N %d|PROCESSORS %d|reps %d|total %9.4f s|startup (beside, not subtracted) %9.4f s|converged %d|iterations not recorded\n",
                tag, p, max(Y) + 1L, nrow(Y), np, NREP, res$t, triv$t, res$converged))
    flush.console()
  }
}

for (p in c(10, 20, 40)) for (K in c(3, 5)) for (N in c(1000, 10000)) {
  nf <- p / 10
  model <- unlist(lapply(seq_len(nf), function(f) {
    v <- sprintf("y%d", (f - 1) * 10 + 1:10); v[1] <- paste0(v[1], "*")
    c(sprintf("  f%d BY", f), wrap(v, "    "), "  ;", sprintf("  f%d@1;", f))
  }))
  design("grid", file.path(TDIR, sprintf("grid_p%d_K%d_N%d.csv", p, K, N)), model)
}
for (nm in c("mc_d", "mc_a1p_16", "mc_b1_m1", "mc_c_998", "mc_B_9975")) {
  csv <- file.path(TDIR, paste0(nm, ".csv"))
  p <- ncol(read.csv(csv, nrows = 1))
  v <- sprintf("y%d", seq_len(p)); v[1] <- "y1@1"
  design(nm, csv, c("  f BY", wrap(v, "    "), "  ;", "  f;"))
}
cat("TIME|done\n")
