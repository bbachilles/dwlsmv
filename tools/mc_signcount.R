# mc_signcount.R -- Module 2 option B: the development-time multiple-root test.
#
# Drafted by a read-only agent, 2026-09-30. EXECUTED by the main session on
# 2026-09-30 over all 42 conditions, 6 cores, ~70 min: 1,086,996 pairs, 0 with two
# or more sign changes, 0 gaps, 0 minima (record: @sec-rootfind, E-056).
#
# What it is for.  Derivation @sec-rootfind ("Uniqueness", "Consequence for the
# algorithm"): no global scan in production; multiple roots are tested for during
# development by "a sign-change count on a fine grid over every pair of the paired
# Monte Carlo re-run".  The endpoint check sees an odd number of extra sign changes
# on one side; an even number goes unseen, and that is what this test is for.
#
# What it does.  For every variable pair of every replicate of the 42 paired
# Monte Carlo conditions (fixtures/replication/index.csv; the same replicates
# tools/paired_mc.do fits), with Module 1's thresholds held fixed (ec_thr()):
#   - evaluates, at every point of the grid G below, each observed cell's computed
#     pi, feasibility (min_O pi > 0), resolvability (min_O pi > delta_P), the
#     observed-cells mean score S/N and its evaluation-error bound tau_S -- the
#     same quantities, formulas and operation order as the .qmd's rf_eval();
#   - classifies each resolvable point: +1 if S/N > tau_S, -1 if S/N < -tau_S,
#     0 ("undetermined") if |S/N| <= tau_S, exact zeros included (the section's
#     strict rule: an exact 0 is not a root, and |S| <= tau_S has no certain sign).
#     Unresolvable points get no class;
#   - counts sign changes in the sequence of decisive (+1/-1) classes, ordered by
#     rho, skipping undetermined and unresolvable points.  Reported both over all
#     of R_P ∩ [-B, B] (changes_all) and over the connected run of resolvable grid
#     points containing 0 (changes_main, the piece option B searches);
#   - a "gap" is R_P ∩ G forming more than one run of consecutive grid points.
#
# Where the numbers come from.  delta_Phi, delta_P, eps, and the reference
# evaluator rf_eval() (with cellP, dP_drho, Phi2, phi2, ec_thr, rf_absD) are
# SOURCED from spec/derivations/dwlsmv-derivations.qmd: the named top-level
# assignments are parsed out of chunks setup, gate3-core, ec-core and rf-core and
# evaluated alone (nothing else in those chunks runs).  rf_eval() uses the .qmd's
# integrate()-based Phi2 (~0.1 s per evaluation), far too slow for ~1.1e6 pairs x
# 1523 points, so the scan uses fast_eval(): the same formulas vectorised over the
# grid, with corners from pbivnorm (Genz's BVND, the algorithm family of Stata's
# binormal() and mvtnorm's TVPACK).  Before the scan, a VALIDATION step compares
# fast_eval() with the .qmd's rf_eval() on a subset fixed in advance (replicate 1,
# pair (1,2), of every condition, at the 11 rho values VAL_X below) and reports
# class disagreements.  It stops only if, at a point both evaluators call
# resolvable and decisive, the two decisive signs differ.
#
# THE GRID, chosen before any result was looked at, and not to be changed after
# running.  G = { tanh(k * h_z) : |k h_z| < atanh(B) } ∪ { 0, -B, +B },
#   h_z = 0.005 in Fisher-z,  B = 0.999,  1523 points.
# Rationale:
#   1. tau_rho (<= 5.4e-12 on Level R1) is not a usable scale: a grid at it would
#      be ~4e11 points.  No root-separation bound is derived (@sec-rootfind:
#      "a scan whose spacing cannot be derived without a bound on the number of
#      roots"), so the spacing is a coverage choice, stated as such.
#   2. The paired MC's high-rho conditions put roots at .99-.999, within 1e-3 of
#      B.  A uniform-in-rho grid would give one or two points there.  Uniform in
#      z = atanh(rho) the spacing is d rho = h_z (1 - rho^2): .005 at 0, 1e-4 at
#      .99, 1e-5 at .999.  The polychoric information grows like 1/(1-rho^2)^2 near
#      the boundary, so equal steps in z are roughly equal steps in the score's own
#      scale.
#   3. d rho <= .005 everywhere, so G is at least as fine as Level R3's grid
#      (spacing .005 in rho) at every rho, and far finer near +-1.
#   4. h_z = .005 is about 1/3 of the sampling SE of z at the largest N (5000;
#      SE(z) ~ 1/sqrt(N) = .014, heuristic, Pearson-z) and smaller relative to it
#      at every smaller N.  Two extra sign changes closer together than one grid
#      step go unseen; that is the test's stated blind spot.
#   5. Cost: ~1.1e6 pairs x 1523 points is feasible in parallel (runtime below).
#
# Excluded replicates.  Replicates with an empty category (index.csv empty_cat != 0;
# mc_b1_m1/m3, mc_b2_m1/m3, 834 in all) are Module 1's deliberate abort in
# tools/paired_mc.do and are never fitted, so they are skipped and counted.  A pair
# whose table has an unobserved category in 0..max(y) is also skipped and flagged.
#
# Usage:
#   Rscript mc_signcount.R OUT.csv [CORES] [CONDITIONS] [MAXREP]
#     OUT.csv      per-pair CSV (LF, %.17g).  OUT_summary.csv: per condition.
#                  OUT_validation.csv: the fast/.qmd comparison.
#     CORES        worker processes (default: detectCores() - 2)
#     CONDITIONS   comma-separated subset, or "all" (default all)
#     MAXREP       replicates 1..MAXREP only (default all 1000)
#   Run from anywhere; the repository path is REPO below.
#
# Smoke test (grid unchanged): see the agent's report.
# -----------------------------------------------------------------------------

REPO <- normalizePath(".")
QMD  <- file.path(REPO, "spec", "derivations", "dwlsmv-derivations.qmd")
B    <- 0.999
H_Z  <- 0.005                                            # fixed before any run
VAL_X <- c(-0.999, -0.99, -0.9, -0.5, 0, 0.5, 0.9, 0.98, 0.99, 0.995, 0.999)   # fixed before any run

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript mc_signcount.R OUT.csv [CORES] [CONDITIONS] [MAXREP]")
OUT    <- args[1]
CORES  <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)
CONDS  <- if (length(args) >= 3 && args[3] != "all") strsplit(args[3], ",")[[1]] else NULL
MAXREP <- if (length(args) >= 4) as.integer(args[4]) else Inf
suppressPackageStartupMessages({ library(pbivnorm); library(haven); library(parallel) })

# ---- Source the .qmd's definitions (named assignments only) -----------------
qmd_defs <- function(qmd, want) {
  L <- readLines(qmd, warn = FALSE); env <- new.env()
  for (lab in names(want)) {
    s <- grep(sprintf("^```\\{r %s[,}]", lab), L)
    if (length(s) != 1) stop("chunk not found once: ", lab)
    e <- s + which(L[(s + 1):length(L)] == "```")[1]
    ex <- parse(text = L[(s + 1):(e - 1)], keep.source = FALSE)
    got <- character(0)
    for (x in ex) if (is.call(x) && as.character(x[[1]]) %in% c("<-", "=") &&
                      is.name(x[[2]]) && as.character(x[[2]]) %in% want[[lab]]) {
      eval(x, env); got <- c(got, as.character(x[[2]]))
    }
    miss <- setdiff(want[[lab]], got)
    if (length(miss)) stop("not found in chunk ", lab, ": ", paste(miss, collapse = ", "))
  }
  env
}
Q <- qmd_defs(QMD, list(
  "setup"      = c("Phi2", "phi2"),
  "gate3-core" = c("phi2s", "cellP", "dP_drho"),
  "ec-core"    = c("ec_thr"),
  "rf-core"    = c("rf_eps", "rf_ulp", "rf_dPhi", "rf_dP", "rf_absD", "rf_eval")))
for (f in ls(Q)) if (is.function(Q[[f]])) environment(Q[[f]]) <- Q   # resolve among themselves
EPS <- Q$rf_eps; DP <- Q$rf_dP; ec_thr <- Q$ec_thr; phi2 <- Q$phi2
stopifnot(EPS == .Machine$double.eps, abs(DP - 4.42e-15) < 1e-17)

# ---- The grid ---------------------------------------------------------------
zB  <- atanh(B)
kk  <- seq_len(ceiling(zB / H_Z) - 1)
kk  <- kk[kk * H_Z < zB]
GRID <- c(-B, rev(-tanh(kk * H_Z)), 0, tanh(kk * H_Z), B)
stopifnot(!is.unsorted(GRID, strictly = TRUE), GRID[1] == -B, GRID[length(GRID)] == B,
          sum(GRID == 0) == 1, max(diff(GRID)) <= 0.005)
cat(sprintf("grid: h_z = %g, %d points, min step %.3g, max step %.3g\n",
            H_Z, length(GRID), min(diff(GRID)), max(diff(GRID))))

# ---- Fast evaluator: rf_eval()'s quantities, vectorised over x --------------
# Same operation order as cellP (C[-1,-1] - C[-n,-1] - C[-1,-m] + C[-n,-m]),
# dP_drho (four phi2s terms, +,-,-,+) and rf_absD (w[-1,-1] + w[-n,-1] + ...).
fast_eval <- function(n, a, b, x) {
  A <- c(-Inf, a, Inf); Bv <- c(-Inf, b, Inf); nA <- length(A); nB <- length(Bv); G <- length(x)
  C <- array(0, c(G, nA, nB)); Ph <- array(0, c(G, nA, nB)); W <- array(0, c(G, nA, nB))
  fin <- which(outer(is.finite(A), is.finite(Bv), "&"), arr.ind = TRUE)
  if (nrow(fin)) {
    ai <- rep(A[fin[, 1]], each = G); bi <- rep(Bv[fin[, 2]], each = G); xi <- rep(x, nrow(fin))
    pv <- pbivnorm(ai, bi, xi)
    ph <- exp(-(ai^2 - 2 * xi * ai * bi + bi^2) / (2 * (1 - xi^2))) / (2 * pi * sqrt(1 - xi^2))  # = phi2
    wv <- ph * ((ai^2 - 2 * xi * ai * bi + bi^2) / (2 * (1 - xi^2)) + 7)
    for (k in seq_len(nrow(fin))) {
      ix <- ((k - 1) * G + 1):(k * G)
      C[, fin[k, 1], fin[k, 2]] <- pv[ix]; Ph[, fin[k, 1], fin[k, 2]] <- ph[ix]; W[, fin[k, 1], fin[k, 2]] <- wv[ix]
    }
  }
  # Phi2 at infinite corners, as the .qmd's Phi2(): -Inf -> 0; (Inf, Inf) -> 1; Inf -> pnorm(other)
  for (r in seq_len(nA)) for (s in seq_len(nB)) if (!(is.finite(A[r]) && is.finite(Bv[s]))) {
    C[, r, s] <- if (A[r] == -Inf || Bv[s] == -Inf) 0 else if (A[r] == Inf && Bv[s] == Inf) 1 else
                 if (A[r] == Inf) pnorm(Bv[s]) else pnorm(A[r])
  }
  i1 <- 2:nA; i0 <- 1:(nA - 1); j1 <- 2:nB; j0 <- 1:(nB - 1)
  P <- C[, i1, j1, drop = FALSE] - C[, i0, j1, drop = FALSE] - C[, i1, j0, drop = FALSE] + C[, i0, j0, drop = FALSE]
  D <- Ph[, i1, j1, drop = FALSE] - Ph[, i0, j1, drop = FALSE] - Ph[, i1, j0, drop = FALSE] + Ph[, i0, j0, drop = FALSE]
  aD <- W[, i1, j1, drop = FALSE] + W[, i0, j1, drop = FALSE] + W[, i1, j0, drop = FALSE] + W[, i0, j0, drop = FALSE]
  O <- which(n > 0)                                     # column-major, as n[O] in rf_eval
  P <- matrix(P, G)[, O, drop = FALSE]; D <- matrix(D, G)[, O, drop = FALSE]; aD <- matrix(aD, G)[, O, drop = FALSE]
  w <- n[O] / sum(n); Wm <- matrix(w, G, length(O), byrow = TRUE)
  mn <- apply(P, 1, min)
  u <- D / P
  list(x = x, pos = mn > 0, res = mn > DP,
       S = rowSums(Wm * u),
       tauS = rowSums(Wm * (EPS * aD + abs(u) * DP) / P) + (length(O) + 2) * EPS * rowSums(Wm * abs(u)))
}

classify <- function(e) ifelse(!e$res, NA_integer_,
                         ifelse(e$S > e$tauS, 1L, ifelse(e$S < -e$tauS, -1L, 0L)))

# ---- Per-pair table, as the .qmd's rf_pair() --------------------------------
pair_tab <- function(yi, yj) {
  n <- unclass(table(factor(yi, 0:max(yi)), factor(yj, 0:max(yj))))
  list(n = n, ok = all(rowSums(n) > 0) && all(colSums(n) > 0) && nrow(n) >= 2 && ncol(n) >= 2,
       a = ec_thr(rowSums(n)), b = ec_thr(colSums(n)))
}

fmt <- function(v) if (is.character(v)) v else if (is.logical(v)) as.character(as.integer(v)) else
  ifelse(is.na(v), "", sprintf("%.17g", v))
write_lf <- function(df, path) {
  body <- do.call(paste, c(lapply(df, fmt), sep = ","))
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(c(paste(names(df), collapse = ","), body), con, sep = "\n", useBytes = TRUE)
}

# ---- Replicates ---------------------------------------------------------------
IDX <- read.csv(file.path(REPO, "fixtures", "replication", "index.csv"), stringsAsFactors = FALSE)
IDX <- IDX[grepl("^mc_", IDX$condition), ]
if (!is.null(CONDS)) { stopifnot(all(CONDS %in% IDX$condition)); IDX <- IDX[IDX$condition %in% CONDS, ] }
IDX <- IDX[IDX$replicate <= MAXREP, ]
IDX$path <- file.path(REPO, IDX$file)
miss <- !file.exists(IDX$path)
if (any(miss)) stop(sum(miss), " replicate files missing (gitignored generated data; ",
                    "regenerate with fixtures/generate/50_replicates.R), e.g. ", IDX$path[which(miss)[1]])
cat(sprintf("%d conditions, %d replicates (%d skipped: empty category)\n",
            length(unique(IDX$condition)), nrow(IDX), sum(IDX$empty_cat != 0)))

# ---- Validation against the .qmd's rf_eval() ---------------------------------
val <- do.call(rbind, lapply(unique(IDX$condition), function(cd) {
  r <- IDX[IDX$condition == cd & IDX$empty_cat == 0, ]
  if (!nrow(r)) return(NULL); r <- r[which.min(r$replicate), ]
  Y <- as.matrix(as.data.frame(read_dta(r$path))); q <- pair_tab(Y[, 1], Y[, 2])
  if (!q$ok) return(NULL)
  fe <- fast_eval(q$n, q$a, q$b, VAL_X); fc <- classify(fe)
  do.call(rbind, lapply(seq_along(VAL_X), function(k) {
    e <- Q$rf_eval(q$n, q$a, q$b, VAL_X[k])
    qc <- if (!e$res) NA_integer_ else if (e$S > e$tauS) 1L else if (e$S < -e$tauS) -1L else 0L
    data.frame(cond = cd, rep = r$replicate, x = VAL_X[k], res_fast = fe$res[k], res_qmd = e$res,
               class_fast = fc[k], class_qmd = qc, S_fast = if (fe$res[k]) fe$S[k] else NA,
               S_qmd = if (e$res) e$S else NA, tauS_fast = if (fe$res[k]) fe$tauS[k] else NA,
               tauS_qmd = if (e$res) e$tauS else NA)
  }))
}))
val$abs_diff_over_tauS <- abs(val$S_fast - val$S_qmd) / val$tauS_fast
write_lf(val, sub("\\.csv$", "_validation.csv", OUT))
both <- !is.na(val$class_fast) & !is.na(val$class_qmd)
cat(sprintf("validation: %d points; resolvability disagrees at %d; class disagrees at %d of %d both-resolvable; max |dS|/tau_S = %.3g\n",
            nrow(val), sum(val$res_fast != val$res_qmd), sum(val$class_fast[both] != val$class_qmd[both]),
            sum(both), max(val$abs_diff_over_tauS, na.rm = TRUE)))
bad <- both & val$class_fast != 0 & val$class_qmd != 0 & val$class_fast != val$class_qmd
if (any(bad)) { print(val[bad, ]); stop("fast and .qmd evaluators give opposite decisive signs") }

# ---- The scan ------------------------------------------------------------------
scan_rep <- function(row) {
  Y <- as.matrix(as.data.frame(read_dta(row$path))); p <- ncol(Y); N <- nrow(Y)
  g <- which(upper.tri(diag(p)), arr.ind = TRUE); g <- g[order(g[, 1], g[, 2]), , drop = FALSE]
  i0 <- which(GRID == 0)
  do.call(rbind, lapply(seq_len(nrow(g)), function(k) {
    i <- g[k, 1]; j <- g[k, 2]; q <- pair_tab(Y[, i], Y[, j])
    base <- data.frame(cond = row$condition, rep = row$replicate, i = i, j = j, N = N,
                       K1 = nrow(q$n), K2 = ncol(q$n), skipped = !q$ok)
    if (!q$ok) return(cbind(base, n_res = NA, n_infeasible = NA, res_runs = NA, gap = NA, res_at_0 = NA,
                            n_undet = NA, class_at_0 = NA, first_sign = NA, last_sign = NA,
                            changes_all = NA, changes_down = NA, changes_up = NA, changes_main = NA,
                            lo_main = NA, hi_main = NA, change_locs = ""))
    e <- fast_eval(q$n, q$a, q$b, GRID); cl <- classify(e)
    rl <- rle(e$res); runs <- sum(rl$values)
    # connected run of resolvable points containing 0
    lo <- i0; hi <- i0
    if (e$res[i0]) { while (lo > 1 && e$res[lo - 1]) lo <- lo - 1; while (hi < length(GRID) && e$res[hi + 1]) hi <- hi + 1 }
    dec <- which(!is.na(cl) & cl != 0); sd <- cl[dec]; ch <- which(diff(sd) != 0)
    dm <- dec[dec >= lo & dec <= hi]; chm <- if (e$res[i0]) sum(diff(cl[dm]) != 0) else NA
    locs <- (GRID[dec[ch]] + GRID[dec[ch + 1]]) / 2
    cbind(base, n_res = sum(e$res), n_infeasible = sum(!e$pos), res_runs = runs, gap = runs > 1,
          res_at_0 = e$res[i0], n_undet = sum(cl == 0, na.rm = TRUE), class_at_0 = cl[i0],
          first_sign = if (length(sd)) sd[1] else NA, last_sign = if (length(sd)) sd[length(sd)] else NA,
          changes_all = length(ch), changes_down = sum(diff(sd) < 0), changes_up = sum(diff(sd) > 0),
          changes_main = chm, lo_main = if (e$res[i0]) GRID[lo] else NA, hi_main = if (e$res[i0]) GRID[hi] else NA,
          change_locs = paste(sprintf("%.17g", locs), collapse = ";"))
  }))
}

todo <- IDX[IDX$empty_cat == 0, ]
t0 <- proc.time()[["elapsed"]]
if (CORES > 1) {
  cl <- makeCluster(CORES); on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, suppressPackageStartupMessages({ library(pbivnorm); library(haven) }))
  clusterExport(cl, c("fast_eval", "classify", "pair_tab", "scan_rep", "GRID", "EPS", "DP", "ec_thr"))
  res <- parLapplyLB(cl, split(todo, seq_len(nrow(todo))), scan_rep)
} else res <- lapply(split(todo, seq_len(nrow(todo))), scan_rep)
R <- do.call(rbind, res); rownames(R) <- NULL
cat(sprintf("scan: %d pairs in %.1f s\n", nrow(R), proc.time()[["elapsed"]] - t0))
write_lf(R, OUT)

# ---- Per-condition summary -------------------------------------------------------
ok <- !R$skipped
S <- do.call(rbind, lapply(split(R, R$cond), function(d) {
  o <- d[!d$skipped, ]
  data.frame(cond = d$cond[1],
             reps = length(unique(d$rep)),
             reps_skipped_empty_category = sum(IDX$condition == d$cond[1] & IDX$empty_cat != 0),
             pairs = nrow(d), pairs_skipped = sum(d$skipped),
             ch0 = sum(o$changes_all == 0), ch1 = sum(o$changes_all == 1), ch2plus = sum(o$changes_all >= 2),
             ch1_up_only = sum(o$changes_all == 1 & o$changes_up == 1),   # a single - to + change: a minimum
             main_ch2plus = sum(o$changes_main >= 2, na.rm = TRUE),
             gap = sum(o$gap), not_res_at_0 = sum(!o$res_at_0),
             undet_at_0 = sum(o$class_at_0 == 0, na.rm = TRUE),
             max_undet = if (nrow(o)) max(o$n_undet) else NA)
}))
write_lf(S, sub("\\.csv$", "_summary.csv", OUT))
print(S, row.names = FALSE)
cat(sprintf("TOTAL pairs %d, skipped %d, 0 changes %d, 1 change %d, >=2 changes %d, gap %d\n",
            nrow(R), sum(R$skipped), sum(R$changes_all[ok] == 0), sum(R$changes_all[ok] == 1),
            sum(R$changes_all[ok] >= 2), sum(R$gap[ok])))
m <- ok & R$changes_all >= 2
if (any(m)) { cat("pairs with >= 2 sign changes:\n")
  print(R[m, c("cond", "rep", "i", "j", "K1", "K2", "changes_all", "changes_main", "gap", "change_locs")], row.names = FALSE) }
