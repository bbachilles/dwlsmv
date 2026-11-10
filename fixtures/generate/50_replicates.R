# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- replicate datasets for the paired validation
#
# Spec Section 7.1 steps 1 and 2: generate many replicate raw item-level
# datasets from known population models, once, in one place, under a recorded
# seed, and write each to a format both tools read identically.
#
# Output goes to fixtures/replication/, which .gitignore excludes. These files
# are regenerated, never committed: there are thousands of them, they are fully
# determined by the seed and this script, and committing them would bloat the
# repository to no purpose. What IS committed is this generator, which with
# its seeds makes the set reproducible. index.csv is written under
# fixtures/replication/ and is therefore gitignored too; an earlier version of
# this header said it was committed, and it never was.
#
# THREE POPULATION MODELS, TWO DATA-GENERATING PROCESSES.
#
# Section 7.4 asks for a single-factor CFA, a correlated two-factor CFA, and a
# small SEM with a structural path between factors. For DATA GENERATION the
# third is not a separate case: a path F1 -> F2 with a residual on F2 induces
# exactly a correlation between the two factors, so the observed data it
# produces is distributionally identical to the correlated two-factor model
# with that correlation. The models differ in what is FITTED to the data, not
# in the data. Generating a third set would produce a set indistinguishable
# from the second and invite the belief that a third process had been covered.
# The structural-path condition is therefore a fitting condition for Module 7,
# and is deliberately absent here.
#
# Edge cases from Section 7.4 ARE separate processes and are included: binary
# indicators (the tetrachoric special case) and a sparse condition combining a
# small sample with the skewed four-category design, which is where empty and
# near-empty cells actually appear.
#
# No package is needed to GENERATE the data -- the simulator below is base R,
# like 00_common.R -- so the replicate set is reproducible on any R
# installation forever. haven is needed only to write .dta.
#
# Run from the repository root:
#     Rscript fixtures/generate/50_replicates.R
#     DWLSMV_NREP=1000 Rscript fixtures/generate/50_replicates.R
#     DWLSMV_CONDITIONS=mc_c_99,mc_d Rscript fixtures/generate/50_replicates.R
#
# DWLSMV_CONDITIONS restricts the run to the named conditions. Seeds are keyed
# on a condition's position in the FULL list, so a subset run writes exactly
# the files a full run would. The mc_* conditions (Module 10's Monte Carlo
# reference, below) take DWLSMV_NREP_MC replicates, default 1000, instead of
# DWLSMV_NREP. Their .dta files are TRANSIENT: generate a condition, fit it
# (70_lavaan_mc_reference.R now, the production paired run later), then
# delete its directory. All of them at once would be ~1.5 GB; the repo is not
# a Monte Carlo warehouse. What persists is this specification and its seeds.
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

if (!requireNamespace("haven", quietly = TRUE))
  stop("haven is required to write replicate .dta files.")

SEED_BASE <- 20260919
NREP <- as.integer(Sys.getenv("DWLSMV_NREP", "200"))
if (is.na(NREP) || NREP < 1) stop("DWLSMV_NREP must be a positive integer.")
NREP_MC <- as.integer(Sys.getenv("DWLSMV_NREP_MC", "1000"))
if (is.na(NREP_MC) || NREP_MC < 1) stop("DWLSMV_NREP_MC must be a positive integer.")
SELECT <- strsplit(Sys.getenv("DWLSMV_CONDITIONS", ""), ",")[[1]]

OUT <- file.path("fixtures", "replication")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# --- the simulator ----------------------------------------------------------
# Deliberately NOT added to 00_common.R. That file is the reference
# implementation of quantities derived in the .qmd, and its header forbids
# changing a formula there without changing the derivation too. This is
# study machinery, not a derived quantity, so it lives with the study.
#
# Draw from a general confirmatory factor model on the latent response
# variables, then cut at thresholds. Residual variances are set so every y*
# has unit total variance, which is the delta parameterization the estimator
# targets -- generating on any other scale would make the population loadings
# incomparable with what lavaan and Mplus report.
simulate_factor_model <- function(n, lambda_mat, phi, thr) {
  p <- nrow(lambda_mat)
  m <- ncol(lambda_mat)
  fac <- matrix(rnorm(n * m), n, m) %*% chol(phi)
  communality <- diag(lambda_mat %*% phi %*% t(lambda_mat))
  if (any(communality >= 1))
    stop("a communality reached 1; the residual variance would be non-positive.")
  resid_sd <- sqrt(1 - communality)
  ystar <- fac %*% t(lambda_mat) +
    matrix(rnorm(n * p), n, p) %*% diag(resid_sd, p, p)
  out <- vapply(seq_len(p), function(j) findInterval(ystar[, j], thr[[j]]),
                integer(n))
  colnames(out) <- sprintf("y%d", seq_len(p))
  as.data.frame(out)
}

# --- population models ------------------------------------------------------
PI4 <- c(0.15, 0.35, 0.30, 0.20)   # skewed 4-category, as in the fixtures
PI3 <- c(0.40, 0.35, 0.25)
PI2 <- c(0.30, 0.70)
thr_of <- function(pmf) qnorm(cumsum(pmf)[-length(pmf)])

lam1 <- matrix(c(0.7, 0.7, 0.7, 0.7), ncol = 1)
lam2 <- matrix(c(0.7, 0.7, 0.7, 0, 0, 0,
                 0, 0, 0, 0.7, 0.7, 0.7), ncol = 2)

CONDITIONS <- list(
  onefactor_n500 = list(
    n = 500, lambda = lam1, phi = matrix(1, 1, 1),
    thr = list(thr_of(PI4), thr_of(PI4), thr_of(PI3), thr_of(PI3)),
    note = "single-factor CFA, four ordinal indicators"),
  onefactor_n2000 = list(
    n = 2000, lambda = lam1, phi = matrix(1, 1, 1),
    thr = list(thr_of(PI4), thr_of(PI4), thr_of(PI3), thr_of(PI3)),
    note = "single-factor CFA, larger N"),
  twofactor_n500 = list(
    n = 500, lambda = lam2, phi = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    thr = rep(list(thr_of(PI3)), 6),
    note = "correlated two-factor CFA; also the DGP for the structural-path fit"),
  twofactor_n2000 = list(
    n = 2000, lambda = lam2, phi = matrix(c(1, 0.5, 0.5, 1), 2, 2),
    thr = rep(list(thr_of(PI3)), 6),
    note = "correlated two-factor CFA, larger N"),
  binary_n500 = list(
    n = 500, lambda = lam1, phi = matrix(1, 1, 1),
    thr = rep(list(thr_of(PI2)), 4),
    note = "binary indicators; the tetrachoric special case"),
  sparse_n200 = list(
    n = 200, lambda = lam1, phi = matrix(1, 1, 1),
    thr = list(thr_of(PI4), thr_of(PI4), thr_of(PI4), thr_of(PI4)),
    note = "small N with the skewed 4-category design; empty cells expected")
)

# --- Module 10: the Monte Carlo SE reference, by partition -------------------
# E-054, Rank-deficient Gamma, step 3: the Monte
# Carlo SD of theta-hat against the mean reported SE as sigma_min -> 0. Each
# mode of derivation @sec-gamma-rank is a SEQUENCE of DGPs, not one point.
# Appended after the six conditions above, so none of their seeds move. Never
# reorder or delete an entry below: a condition's seeds are its position.
#
# partition P1  numerical ill-conditioning; calibrates the cutoff
#   mc_d        case D's DGP at its own N; the healthy control
#   mc_c_*      near-duplicate pair: lambda_1 = lambda_2 = sqrt(rho), so
#               rho_12 = rho (design c), rho = .9 ... .995
#   mc_B_*      uniformly high loadings, lambda = .8 ... .995 on every item
#               (max rho .64 ... .99). The second clean geometry: Gamma
#               deteriorates as in c, the structural bread does not
#               (derivation @sec-gamma-pop; candidate A, a three-item
#               cluster, was rejected there as Gamma-identical to c)
# partition P2  exact deficiency; informs the rank-deficient-Gamma policy
#   mc_a1p_*    design a1's items (K = 3, lambda .7), p = 8 ... 16 at N = 100:
#               n_pat - 1 crosses D between p = 12 and 14
#   mc_a2p_14   binary items (q = .3, lambda .7), p = 14 at N = 100: the
#               binary analogue of mc_a1p_* at the same N. D = 105 > N - 1,
#               so n_pat - 1 < D on EVERY draw, by counting. Not redrawn on an
#               empty 2x2 cell. One fixed point; the mc_a1p_* sequence carries
#               the trajectory.
# partition S   secondary stress; not used to calibrate anything
#   mc_b1_m*    y1's BOUNDARY category with expected count m (b1 is m = 1)
#   mc_b2_m*    y1's INTERIOR category with expected count m (b2 is m = 1)
#   mc_c_998    c's boundary point, reclassified from P1 on 2026-09-25: the
#               P1 screen found 481/1000 fits inadmissible (Heywood, residual
#               variance near 0) and its one |z| excursion (3.9) confounded
#               with that. Boundary evidence, not clean cutoff evidence.
#   mc_c_999    the bounded-estimator point: ~half of all rho-hat at 0.999
#   mc_B_9975   B's boundary point (max rho .995; 9/50 rho-hat capped in the
#               one-off screen pilot), kept as a documented stress result
#   mc_a1n_*    design a1 over N (N is confounded with rank). Its N = 100
#               point is mc_a1p_10's data, not a second draw.
#   mc_a2n_*    design a2 (binary, p = 12, D = 78) over N = 40 ... 200, not
#               redrawn on an empty 2x2 cell. An N sequence, so N is
#               confounded with rank as in mc_a1n_*: stress, not P2.
#   (Every condition with a pair of binary items is also fitted at lavaan's
#   default zero.add, by 70_; that is a second FIT, not data, and it is S.)
# partition R   reproducibility
#   mc_d_dup    mc_d's DGP on its own seed stream; a seed-independence check,
#               never a design point
#
# Every mc_* condition is fitted with marker-variable identification, the
# first indicator's loading fixed to 1; `ident` records it.
#
# b1/b2 draws with an empty category ARE written. Module 1 aborts on them and
# lavaan's ordered() silently drops the level, so both pipelines must see them
# and record the outcome; the index flags them (empty_cat) rather than this
# script dropping them, which would make the reference conditional unseen.
N_MC  <- 5000
thrD  <- list(thr_of(PI4), thr_of(PI4), thr_of(PI3), thr_of(PI3))
lamD  <- function(l12 = 0.7) matrix(c(l12, l12, 0.7, 0.7), ncol = 1)
ONE   <- matrix(1, 1, 1)
mcc <- function(n, lambda, thr, partition, note)
  list(n = n, lambda = lambda, phi = ONE, thr = thr, partition = partition,
       ident = "marker: first indicator", nrep = NREP_MC, note = note)
tag <- function(x) sub("^0[.]", "", format(x))
mc <- list(mc_d = mcc(N_MC, lamD(), thrD, "P1", "case D's DGP; the healthy control"))
# mc_c_998 keeps its place in this loop, so its seeds do not move, but it is
# labelled S: reclassified 2026-09-25 after the P1 screen (E-060).
for (r in c(0.9, 0.95, 0.98, 0.99, 0.995, 0.998))
  mc[[paste0("mc_c_", tag(r))]] <- mcc(N_MC, lamD(sqrt(r)), thrD,
    if (r == 0.998) "S" else "P1",
    sprintf("near-duplicate pair, rho_12 = %g", r))
for (l in c(0.8, 0.9, 0.95, 0.97, 0.98, 0.99, 0.995))
  mc[[paste0("mc_B_", tag(l))]] <- mcc(N_MC, matrix(l, 4, 1), thrD, "P1",
    sprintf("uniform loadings %g on case D's items", l))
for (p in c(8, 10, 12, 14, 16))
  mc[[sprintf("mc_a1p_%d", p)]] <- mcc(100, matrix(0.7, p, 1), rep(list(thr_of(PI3)), p), "P2",
    sprintf("design a1's items, p = %d at N = 100", p))
mc$mc_a2p_14 <- mcc(100, matrix(0.7, 14, 1), rep(list(qnorm(0.3)), 14), "P2",
  "binary items, p = 14 at N = 100: D = 105 > N - 1, deficient on every draw")
for (m in c(1, 3, 10, 30, 100)) {
  mc[[sprintf("mc_b1_m%d", m)]] <- mcc(N_MC, lamD(),
    c(list(thr_of(c(m / N_MC, 0.50 - m / N_MC, 0.30, 0.20))), thrD[-1]), "S",
    sprintf("y1 boundary category, expected count %d", m))
  mc[[sprintf("mc_b2_m%d", m)]] <- mcc(N_MC, lamD(),
    c(list(thr_of(c(0.15, m / N_MC, 0.65 - m / N_MC, 0.20))), thrD[-1]), "S",
    sprintf("y1 interior category, expected count %d", m))
}
mc$mc_c_999 <- mcc(N_MC, lamD(sqrt(0.999)), thrD, "S",
  "near-duplicate pair, rho_12 = .999: the bounded-estimator point")
mc$mc_B_9975 <- mcc(N_MC, matrix(0.9975, 4, 1), thrD, "S",
  "uniform loadings .9975 (max rho .995): B's boundary point")
for (n in c(40, 50, 65, 200, 500))
  mc[[sprintf("mc_a1n_%d", n)]] <- mcc(n, matrix(0.7, 10, 1), rep(list(thr_of(PI3)), 10), "S",
    sprintf("design a1 at N = %d", n))
for (n in c(40, 60, 100, 200))
  mc[[sprintf("mc_a2n_%d", n)]] <- mcc(n, matrix(0.7, 12, 1), rep(list(qnorm(0.3)), 12), "S",
    sprintf("design a2 at N = %d, not redrawn; N confounded with rank", n))
mc$mc_d_dup <- mcc(N_MC, lamD(), thrD, "R", "mc_d's DGP on its own seed stream")
CONDITIONS <- c(CONDITIONS, mc)
if (length(SELECT) && !all(SELECT %in% names(CONDITIONS)))
  stop("unknown condition(s): ", paste(setdiff(SELECT, names(CONDITIONS)), collapse = ", "))

# --- generate ---------------------------------------------------------------
index <- data.frame()
empty_cell_reports <- 0L

for (ci in seq_along(CONDITIONS)) {
  nm <- names(CONDITIONS)[ci]
  if (length(SELECT) && !(nm %in% SELECT)) next   # ci still counts: seeds unmoved
  cond <- CONDITIONS[[ci]]
  nrep <- if (is.null(cond$nrep)) NREP else cond$nrep
  cdir <- file.path(OUT, nm)
  dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
  cindex <- data.frame()

  for (r in seq_len(nrep)) {
    # One seed per (condition, replicate), derived rather than sequential, so
    # adding a condition later cannot shift the seeds of existing ones.
    seed <- SEED_BASE + ci * 1000000L + r
    set.seed(seed)
    dat <- simulate_factor_model(cond$n, cond$lambda, cond$phi, cond$thr)

    # Count pairs with an empty cell. This is not a failure -- the sparse
    # condition exists to produce them -- but it must be recorded per
    # replicate, because Section 7.1 asks for discrepancy patterns against
    # dataset characteristics, and sparsity is the characteristic that
    # matters most here.
    #
    # The table is built over every DECLARED category, not the observed ones.
    # table() on the raw codes drops an unobserved level, and with it that
    # level's whole row of empty cells, so this count used to be too low
    # exactly when a category was empty -- the replicate where it matters.
    nempty <- 0L; mincell <- .Machine$integer.max
    for (a in seq_len(ncol(dat) - 1L)) {
      for (b in seq(a + 1L, ncol(dat))) {
        tab <- table(factor(dat[[a]], levels = 0:length(cond$thr[[a]])),
                     factor(dat[[b]], levels = 0:length(cond$thr[[b]])))
        nempty <- nempty + sum(tab == 0L)
        mincell <- min(mincell, tab)
      }
    }
    empty_cell_reports <- empty_cell_reports + (nempty > 0L)

    f <- file.path(cdir, sprintf("%s_%04d.dta", nm, r))
    haven::write_dta(dat, f)

    # Module 10's dataset characteristics (spec 4.5.6): the pattern count on
    # Y as Stage 1 receives it (no missing values here, so nothing is deleted),
    # the smallest category count, and the number of declared categories that
    # are empty. empty_cat > 0 is a replicate Module 1 aborts on.
    cnt <- lapply(seq_along(cond$thr), function(j)
      tabulate(dat[[j]] + 1L, length(cond$thr[[j]]) + 1L))
    cindex <- rbind(cindex, data.frame(
      condition = nm, replicate = r, seed = seed, n = cond$n,
      p = ncol(dat), empty_cells = nempty,
      n_pat = nrow(unique(dat)),
      min_cat = min(unlist(cnt)),
      min_cell = mincell,
      empty_cat = sum(unlist(cnt) == 0L),
      file = f, stringsAsFactors = FALSE
    ))
  }
  # One index per condition, so a subset run cannot truncate the others'.
  write.csv(cindex, file.path(cdir, "index.csv"), row.names = FALSE)
  cat(sprintf("  %-18s %4d replicates, n = %4d, p = %d, %d with an empty category\n",
              nm, nrep, cond$n, nrow(cond$lambda), sum(cindex$empty_cat > 0L)))
}

# The combined index is rebuilt from every per-condition index on disk.
index <- do.call(rbind, lapply(names(CONDITIONS), function(nm) {
  f <- file.path(OUT, nm, "index.csv")
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE)
}))
empty_cell_reports <- sum(index$empty_cells > 0L)
write.csv(index, file.path(OUT, "index.csv"), row.names = FALSE)

# Population values, so bias and coverage can be computed later without
# re-deriving what generated the data.
pop <- do.call(rbind, lapply(names(CONDITIONS), function(nm) {
  cond <- CONDITIONS[[nm]]
  data.frame(condition = nm, n = cond$n,
             loading = I(list(cond$lambda)), phi = I(list(cond$phi)),
             thresholds = I(list(cond$thr)),   # needed for threshold bias/SD
             nrep = if (is.null(cond$nrep)) NREP else cond$nrep,
             partition = if (is.null(cond$partition)) NA else cond$partition,
             ident = if (is.null(cond$ident)) NA else cond$ident,
             seed_first = SEED_BASE + match(nm, names(CONDITIONS)) * 1000000L + 1L,
             note = cond$note, stringsAsFactors = FALSE)
}))
saveRDS(pop, file.path(OUT, "population.rds"))

cat(sprintf("\n%d replicate file(s) across %d condition(s) in %s\n",
            nrow(index), length(CONDITIONS), OUT))
cat(sprintf("%d replicate(s) contain at least one empty pair cell\n",
            empty_cell_reports))
cat("Not committed -- .gitignore excludes fixtures/replication/. The generator\n")
cat("and index.csv are what make this set reproducible.\n")
