# ---------------------------------------------------------------------------
# EXECUTED 2026-11-10 by the main session from the repository root (R 4.6.0;
# log docs/notes/condx-validation/run_99m.log). Drafted by a read-only agent;
# the replicate-1 self-check was written by the main session.
#
# dwlsmv fixture generation -- gate link 3 datasets (conditional mode)
#
# Regenerates the datasets of the 100 selected replicates of the frozen
# conditional.x simulation (prereg docs/notes/condx-link3/prereg.md sections 2
# and 8, [L1]): replicates 1-25 of each of C500, C1500, D500, D1500. The draws
# are replayed, not stored: seed 20261026 + 10 * cell index, then one
# sample.int() per replicate in replicate order (condx_sim_run.R, cxs_rep).
# The population (cx_population) comes from the pinned oracle chunks, exactly
# as docs/notes/condx-link3/link3_extract.R does.
#
# Outputs (all under OUTROOT, default "fixtures"):
#   data/condx_link3.csv  cell,rep,y1,y2,y3,x -- one row per observation.
#       Row order: cells C500, C1500, D500, D1500; replicates 1..25; within a
#       replicate, the population's pattern order (expand.grid, x = 0 block
#       then x = 1 block), each pattern repeated count times. The simulation
#       tabulates the draw, so the original within-replicate observation order
#       does not exist. Row order reaches the estimators only through summation rounding.
#       ONE LONG FILE rather than 100: fixtures/data/ holds one file per
#       dataset only where a dataset is a case (condx_A.csv); here the 100 are
#       one design. 00_common.R sets no one-file-per-dataset rule; a single
#       file (about 100,000 rows) needs one manifest hash, is one
#       import delimited in Stata followed by a keep on cell and rep, and
#       cannot lose a member.
#   data/condx_D.csv      cell D1500 replicate 1 alone, y1,y2,y3,x (the format
#       of condx_A.csv). The link 2 dataset.
#   MANIFEST_condx_link3.json  (at OUTROOT, i.e. fixtures/)
#
# Checks (stop on failure): SHA-256 of the pinned oracle text (git show
# 6e37ddd:...) and of condx_sim_run.R against identities.md; replayed pattern
# count and count checksum against docs/notes/condx-link3/link3_ids.csv for all
# 100 replicates (ENABLED). The comparison against the frozen rows of
# replicate 1 always runs (the self-check section below).
#
# What this does not certify: that the replay IS the frozen draw (the frozen run
# used R 4.5.2; the RNG kinds are pinned, the rows.csv replicate-1 comparison is
# what establishes it); any estimate. Writes no theta-hat, SE or T3.
#
# Run from the repository root:  Rscript fixtures/generate/99m_condx_link3.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

# ---- switches and paths -----------------------------------------------------
# (the replicate-1 self-check against the frozen rows always runs; see below)
OUTROOT <- Sys.getenv("OUTROOT", "fixtures")   # tests point this away from the repository

ORACLE_COMMIT <- "6e37ddd"
ORACLE_PATH <- "spec/derivations/dwlsmv-derivations.qmd"
ORACLE_SHA <- "542808567127005e4529e9be6c58d5975b8a3d6df736ab7f163906171d5f67b1"
SIM_COMMIT <- "0dd8174"
SIM_PATH <- file.path("docs", "notes", "condx-sim", "condx_sim_run.R")
SIM_SHA <- "ebc26cb8fc936868d5df2cb99890a5f3bfc282088a767cfda2673de03f36a2f9"
IDS_PATH <- file.path("docs", "notes", "condx-link3", "link3_ids.csv")

DAT_DIR <- file.path(OUTROOT, "data")
f_long <- file.path(DAT_DIR, "condx_link3.csv")
f_D <- file.path(DAT_DIR, "condx_D.csv")
MAN <- file.path(OUTROOT, "MANIFEST_condx_link3.json")
dir.create(DAT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- pins -------------------------------------------------------------------
sha256 <- function(path) {
  out <- tryCatch(suppressWarnings(system2("sha256sum", shQuote(path), stdout = TRUE, stderr = FALSE)),
                  error = function(e) character(0))
  if (length(out) && grepl("^[0-9a-f]{64}", out[1])) return(substr(out[1], 1, 64))
  out <- system2("certutil", c("-hashfile", shQuote(normalizePath(path)), "SHA256"), stdout = TRUE)
  tolower(gsub("[^0-9a-fA-F]", "", out[2]))
}
tmp_q <- tempfile(fileext = ".qmd")
rc <- system2("git", c("show", shQuote(paste0(ORACLE_COMMIT, ":", ORACLE_PATH))), stdout = tmp_q)
stopifnot(identical(as.integer(rc), 0L), file.exists(tmp_q))
if (!identical(sha256(tmp_q), ORACLE_SHA))
  stop("oracle SHA-256 differs from identities.md: ", sha256(tmp_q))
if (!identical(sha256(SIM_PATH), SIM_SHA))
  stop("condx_sim_run.R SHA-256 differs from identities.md: ", sha256(SIM_PATH))

# ---- replay (as link3_extract.R, which copies condx_sim_run.R's draw) -------
qmd <- readLines(tmp_q, encoding = "UTF-8")
chunk <- function(lab) {
  s <- grep(paste0("^```\\{r ", lab, "\\}"), qmd); stopifnot(length(s) == 1)
  e <- s + which(qmd[(s + 1):length(qmd)] == "```")[1]
  qmd[(s + 1):(e - 1)]
}
for (lab in c("cx-setup", "cx-stage1", "cx-model"))
  eval(parse(text = chunk(lab)), globalenv())      # cx_true, cx_population, cx_p, ...

SEED0 <- 20261026L
tau_dif <- list(cx_true$tau[[1]] + c(0, 0.4), cx_true$tau[[2]], cx_true$tau[[3]] + 0.3)
pops <- list(C = cx_population(cx_true), D = cx_population(cx_true, tau_x1 = tau_dif))
cells <- expand.grid(N = c(500L, 1500L), dgp = c("C", "D"), stringsAsFactors = FALSE)[, 2:1]
cells$id <- seq_len(nrow(cells))                   # C500 = 1, C1500 = 2, D500 = 3, D1500 = 4
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
REPS <- 1:25                                       # [L1]

t0 <- proc.time()[["elapsed"]]
blocks <- list(); idrows <- list(); firstrep <- list()
for (k in seq_len(nrow(cells))) {
  pop <- pops[[cells$dgp[k]]]; N <- cells$N[k]; cid <- paste0(cells$dgp[k], N)
  set.seed(SEED0 + 10L * cells$id[k])
  for (rep in REPS) {                              # one sample.int per replicate, in order
    cnt <- tabulate(sample.int(length(pop$w), N, TRUE, pop$w), length(pop$w))
    stopifnot(sum(cnt) == N)
    ix <- rep(seq_along(cnt), cnt)                 # pattern order, each repeated count times
    Y <- cbind(pop$Y[ix, , drop = FALSE], pop$x[ix]); colnames(Y) <- c("y1", "y2", "y3", "x")
    blocks[[length(blocks) + 1]] <- data.frame(cell = cid, rep = rep, Y)
    idrows[[length(idrows) + 1]] <- data.frame(cell = cid, cell_index = cells$id[k],
      seed = SEED0 + 10L * cells$id[k], rep = rep, n_patterns = sum(cnt > 0),
      count_checksum = sum(cnt * seq_along(cnt)), N = N)
    if (rep == 1) firstrep[[cid]] <- list(pattern = which(cnt > 0), count = cnt[cnt > 0])
  }
}
ids <- do.call(rbind, idrows)
long <- do.call(rbind, blocks); rownames(long) <- NULL
vals <- unlist(long[c("y1", "y2", "y3", "x")])
stopifnot(nrow(ids) == 100, nrow(long) == 25 * (500 + 1500 + 500 + 1500), all(vals == round(vals)))

# ---- ENABLED check: replayed identities equal link3_ids.csv, all 100 -------
ref <- read.csv(IDS_PATH, colClasses = "character", stringsAsFactors = FALSE)
stopifnot(nrow(ref) == 100)
chk <- c("cell", "cell_index", "seed", "rep", "n_patterns", "count_checksum")
mine <- as.matrix(data.frame(lapply(ids[chk], as.character)))
theirs <- as.matrix(ref[chk]); dimnames(mine) <- dimnames(theirs) <- NULL
if (!identical(mine, theirs)) stop("replayed identities differ from ", IDS_PATH, " in ",
  sum(rowSums(mine != theirs) > 0), " of 100 replicates")
cat("link3_ids.csv identity check: all 100 replicates equal\n")

# ---- DISABLED hook: replicate 1 of each cell against the frozen rows -------
# Set SELF_CHECK_ROWS <- TRUE above to run it. It may read only replicate-1
# rows of docs/notes/condx-sim/rows.csv (prereg section 2: replicate 1 of each
# cell, refitted, must reproduce its frozen row). The body is a stub: the
# columns to compare (e.g. status, T3, th) and the tolerance belong to the
# prereg, not to this draft.
#
# Written by the main session, 2026-11-10 (the user's condition: the generator
# must prove it uses the exact frozen replicate). ALWAYS RUN; the flag above
# no longer disables it.
#
# For replicate 1 of each cell, the frozen script's own cxs_rep() is run,
# verbatim, on the rows about to be written: sample.int is masked so that the
# "draw" inside cxs_rep is the pattern index of each row of `long`. The result
# is compared with that cell's replicate-1 row of rows.csv. Only replicate-1
# lines of rows.csv are parsed; no value from it is printed, only the verdict
# and the largest relative difference.
#
# Rule, set before the first run:
#   EXACT  status and iteration count identical, and every numeric field
#          (T, a, b, T3, th1-9, seh1-9, seb1-9) equal to the frozen value after
#          rounding to the 15 significant digits write.csv kept;
#   DRAW   not EXACT, but status and iterations identical and every numeric
#          field within 1e-9 relative (+1e-12 absolute): the same data, with
#          arithmetic differing in the last digits (the frozen run used
#          R 4.5.2). A different draw would differ at about 1e-2;
#   FAIL   anything else: stop, nothing written.
SELF_CHECK_ROWS <- TRUE
for (ex in parse(SIM_PATH, keep.source = FALSE)) {          # definitions only; the study is not run
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) &&
      as.character(ex[[2]]) %in% c("GN_TOL", "GN_MAX", "th0", "nth", "cxs_gn", "cxs_start", "cxs_rep"))
    eval(ex, globalenv())
}
eval(parse(text = chunk("cxn-psi")), globalenv())            # cx_psi_n, cx_gamma, pinned text
stopifnot(exists("cxs_rep"), exists("cx_psi_n"), GN_TOL == 1e-10)
fr_lines <- readLines(file.path("docs", "notes", "condx-sim", "rows.csv"))
fr_keep <- c(TRUE, vapply(strsplit(fr_lines[-1], ","), function(z) z[3] == "1", NA))
fr1 <- read.csv(text = fr_lines[fr_keep], stringsAsFactors = FALSE); rm(fr_lines)
stopifnot(nrow(fr1) == 4, all(fr1$rep == 1))
num_cols <- c("T", "a", "b", "T3", paste0("th", 1:9), paste0("seh", 1:9), paste0("seb", 1:9))
self <- list()
for (k in seq_len(nrow(cells))) {
  pop <- pops[[cells$dgp[k]]]; N <- cells$N[k]; cid <- paste0(cells$dgp[k], N)
  rows1 <- long[long$cell == cid & long$rep == 1, c("y1", "y2", "y3", "x")]
  ix <- match(do.call(paste, rows1), paste(pop$Y[, 1], pop$Y[, 2], pop$Y[, 3], pop$x))
  stopifnot(nrow(rows1) == N, !anyNA(ix))
  sample.int <- function(...) ix                             # masks base::sample.int for cxs_rep
  o <- cxs_rep(pop, N); rm(sample.int)
  mine <- c(o$T, o$a, o$b, o$T3, o$th, o$se_h, o$se_b)
  fz <- fr1[fr1$dgp == cells$dgp[k] & fr1$N == N, ]; stopifnot(nrow(fz) == 1)
  theirs <- as.numeric(unlist(fz[num_cols]))
  same_si <- identical(o$status, fz$status) && identical(as.integer(o$it), as.integer(fz$it))
  exact <- same_si && all(as.numeric(sprintf("%.15g", mine)) == theirs)
  rel <- max(abs(mine - theirs) / abs(theirs))
  draw <- same_si && all(abs(mine - theirs) <= 1e-9 * abs(theirs) + 1e-12)
  self[[cid]] <- if (exact) "EXACT" else if (draw) "DRAW" else "FAIL"
  cat(sprintf("replicate-1 self-check %s: %s (status and iterations %s; %d of %d fields equal at 15 digits; max relative difference %.3e)\n",
              cid, self[[cid]], if (same_si) "identical" else "DIFFER",
              sum(as.numeric(sprintf("%.15g", mine)) == theirs), length(theirs), rel))
  cat("    equal at 15 digits:", paste(num_cols[as.numeric(sprintf("%.15g", mine)) == theirs], collapse = " "), "
")
  cat(sprintf("    max relative difference by class: T,a,b,T3 %.2e; theta %.2e; SE (Hessian) %.2e; SE (BHHH) %.2e
",
              max((abs(mine - theirs) / abs(theirs))[1:4]), max((abs(mine - theirs) / abs(theirs))[5:13]),
              max((abs(mine - theirs) / abs(theirs))[14:22]), max((abs(mine - theirs) / abs(theirs))[23:31])))
  cat(sprintf("    absolute difference: T %.2e; a %.2e; b %.2e; T3 %.2e; theta max %.2e; SE (Hessian) max %.2e
",
              abs(mine - theirs)[1], abs(mine - theirs)[2], abs(mine - theirs)[3], abs(mine - theirs)[4],
              max(abs(mine - theirs)[5:13]), max(abs(mine - theirs)[14:22])))
}
if (any(unlist(self) == "FAIL")) stop("replicate-1 self-check FAILED: the replay is not shown to be the frozen draw; nothing written")

# ---- write ------------------------------------------------------------------
write_csv_lf(long, f_long)
dD <- long[long$cell == "D1500" & long$rep == 1, c("y1", "y2", "y3", "x")]
stopifnot(nrow(dD) == 1500)
write_csv_lf(dD, f_D)

# ---- manifest (file keys are repo-relative, as verify.R expects) ------------
rel <- c("fixtures/data/condx_link3.csv", "fixtures/data/condx_D.csv")
md5 <- unname(tools::md5sum(c(f_long, f_D))); names(md5) <- rel
esc <- function(x) gsub('"', '\\\\"', as.character(x))
fmt <- function(x, ind = "  ") {
  if (is.list(x)) paste0("{\n", paste(sprintf('%s  "%s": %s', ind, names(x),
      vapply(x, function(v) fmt(v, paste0(ind, "  ")), character(1))), collapse = ",\n"), "\n", ind, "}")
  else if (length(x) > 1) paste0("[", paste(sprintf('"%s"', esc(x)), collapse = ", "), "]")
  else sprintf('"%s"', esc(x))
}
info <- list(
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string,
  platform = R.version$platform,
  seeds = paste0("20261026 + 10 * cell index; C500 20261036, C1500 20261046, D500 20261056, D1500 20261066; ",
                 "set.seed once per cell, then one sample.int per replicate, in replicate order"),
  rng = paste(RNGkind(), collapse = ", "),
  files = as.list(md5),
  pins = list(
    oracle = paste0(ORACLE_PATH, " at ", ORACLE_COMMIT, ", SHA-256 ", ORACLE_SHA),
    simulation = paste0(SIM_PATH, " at ", SIM_COMMIT, ", SHA-256 ", SIM_SHA)),
  selection = "L1: replicates 1-25 of C500, C1500, D500, D1500 (100 datasets); rows ordered cell, rep, pattern",
  n_rows = list(condx_link3 = nrow(long), condx_D = nrow(dD)),
  identities = sprintf("%s:%d:%d", ids$cell, ids$rep, ids$N),
  packages = c("none (base R only)",
    "frozen study ran on R 4.5.2; replay on this R verified against link3_ids.csv (100 of 100)",
    "condx_D.csv = cell D1500 replicate 1, the link 2 dataset"))
writeLines(fmt(info), MAN)
cat(sprintf("%d rows to %s, %d to %s; %s written; %.1f s.\n", nrow(long), f_long, nrow(dD), f_D, MAN,
            proc.time()[["elapsed"]] - t0))
