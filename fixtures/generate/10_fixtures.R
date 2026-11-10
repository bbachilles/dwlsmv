# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 1: self-contained, no external packages
#
# Produces the Module 0 canonical reference cases (spec Section 5, Module 0):
#
#   Case A  two binary indicators                      -- the only fully
#                                                         closed-form case
#   Case B  three binary indicators                    -- smallest genuinely
#                                                         joint Gamma
#   Case C  two 3-category ordinal indicators          -- multiple thresholds
#                                                         per variable
#   Case D  four ordinal indicators, unequal categories -- the first case with
#                                                         p > 2 structure
#   Case E  six 3-category indicators, TWO factors       -- the first case that
#                                                         can carry a
#                                                         structural path
#
# For each case this writes:
#   fixtures/data/<case>.csv          raw ordinal data, integer-valued
#   fixtures/expected/<case>_thresholds.csv
#   fixtures/expected/<case>_gamma_tau.csv   analytic threshold block of Gamma
#
# and, from two embedded pair tables rather than a case (see the block before
# the manifest):
#   fixtures/expected/highrho_polychoric.csv counts, thresholds and rho-hat of
#                                            two near-duplicate pairs
#
# CSV rather than .dta is deliberate. These are small and integer-valued, so
# there is no precision question, and CSV diffs in git while a .dta is an
# opaque blob. The bulk replication study of spec Section 7.1 is a different
# activity with different needs and should use .dta there.
#
# Run from the repository root:  Rscript fixtures/generate/10_fixtures.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

SEED <- 20260917
N    <- 5000
LAMBDA <- 0.7

dir.create(file.path("fixtures", "data"),     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("fixtures", "expected"), showWarnings = FALSE, recursive = TRUE)

# Category probabilities used throughout. Deliberately unequal and moderately
# skewed: q = .15 on the first threshold of the 4-category items makes psi as
# skewed as it gets at realistic sparsity, which is where the noisiest Gamma
# entry lives.
PI4 <- c(0.15, 0.35, 0.30, 0.20)   # 4 categories -> 3 thresholds
PI3 <- c(0.40, 0.35, 0.25)         # 3 categories -> 2 thresholds
PI2 <- c(0.30, 0.70)               # binary        -> 1 threshold

thr_of <- function(pmf) qnorm(cumsum(pmf)[-length(pmf)])

CASES <- list(
  A = list(thr = list(thr_of(PI2), thr_of(PI2)),
           note = "two binary indicators; saturated, closed-form oracle"),
  B = list(thr = list(thr_of(PI2), thr_of(PI2), thr_of(PI2)),
           note = "three binary indicators; smallest joint Gamma"),
  C = list(thr = list(thr_of(PI3), thr_of(PI3)),
           note = "two 3-category ordinal indicators"),
  D = list(thr = list(thr_of(PI4), thr_of(PI4), thr_of(PI3), thr_of(PI3)),
           note = "four ordinal indicators, unequal category counts")
)

written <- character(0)

for (nm in names(CASES)) {
  thr <- CASES[[nm]]$thr
  p   <- length(thr)

  set.seed(SEED + utf8ToInt(nm))            # per-case seed, reproducible
  Y <- simulate_onefactor(N, LAMBDA, thr)
  colnames(Y) <- sprintf("y%d", seq_len(p))

  # population quantities -- the fixture's expected values are analytic, not
  # estimated, so they do not inherit the sample's noise
  tau    <- unlist(thr)
  q      <- pnorm(tau)
  var_id <- rep(seq_len(p), lengths(thr))
  R      <- matrix(LAMBDA^2, p, p); diag(R) <- 1
  G      <- gamma_tau(tau, q, var_id, R)

  lab <- threshold_labels(lengths(thr))
  dimnames(G) <- list(lab, lab)

  f_data <- file.path("fixtures", "data",     sprintf("case_%s.csv", nm))
  f_thr  <- file.path("fixtures", "expected", sprintf("case_%s_thresholds.csv", nm))
  f_gam  <- file.path("fixtures", "expected", sprintf("case_%s_gamma_tau.csv", nm))

  write_csv_lf(Y, f_data)
  write_csv_lf(data.frame(label = lab, tau = tau, q = q), f_thr)
  write_csv_lf(cbind(label = lab, as.data.frame(G)), f_gam)

  written <- c(written, f_data, f_thr, f_gam)

  # --- self-checks: these must hold or the fixture is not worth committing ---
  stopifnot(isTRUE(all.equal(G, t(G))))                       # symmetric
  stopifnot(min(eigen(G, only.values = TRUE)$values) > 0)     # positive definite

  # the within-variable identity of derivations 4.5.1 / spec 4.5.1: the
  # empirical outer product and the closed form are the same arithmetic
  Z  <- cumulative_indicators(Y, thr)
  qh <- colMeans(Z)
  ps <- sweep(Z, 2, qh) / rep(dnorm(qnorm(qh)), each = N)
  w  <- outer(var_id, var_id, "==")
  gap <- max(abs((crossprod(ps) / N - gamma_tau(qnorm(qh), qh, var_id, R))[w]))
  stopifnot(gap < 1e-10)

  cat(sprintf("case %s: p=%d, T=%d, within-variable identity gap %.2e  [%s]\n",
              nm, p, length(tau), gap, CASES[[nm]]$note))
}

# ---------------------------------------------------------------------------
# Case E -- two factors, six indicators, three categories each.
#
# Cases A-D are all single-factor, which is enough for every Stage 1 check but
# cannot carry a structural path: a factor has to be predicted by another
# factor before "what does unit-variance identification fix for an ENDOGENOUS
# factor" is even a question. Case E exists to make that question askable, and
# to carry the equivalent-models check (a correlated two-factor model and an
# F2 <- F1 model are reparameterizations of each other, so they must fit
# identically).
#
# The population is the `twofactor` condition of 50_replicates.R -- lambda 0.7
# throughout, phi_12 = 0.5, three equiprobable-ish categories -- at this file's
# N rather than that file's, so the two are recognizably the same design.
#
# Deliberately written as its own block rather than folded into the loop above.
# The loop's simulator is single-factor and its R is LAMBDA^2 off the diagonal;
# generalizing it would touch the RNG draw order for A-D, and A-D are committed
# files whose whole value is that they do not move. A separate block costs a
# dozen lines and cannot disturb them.
# ---------------------------------------------------------------------------

LAM_E <- matrix(c(0.7, 0.7, 0.7, 0, 0, 0,
                  0, 0, 0, 0.7, 0.7, 0.7), ncol = 2)
PHI_E <- matrix(c(1, 0.5, 0.5, 1), 2, 2)

thr_E <- rep(list(thr_of(PI3)), 6)
p_E   <- nrow(LAM_E)

R_E <- LAM_E %*% PHI_E %*% t(LAM_E)
diag(R_E) <- 1                      # delta parameterization: Var(y*) = 1

set.seed(SEED + utf8ToInt("E"))
communality <- diag(LAM_E %*% PHI_E %*% t(LAM_E))
stopifnot(all(communality < 1))
fac   <- matrix(rnorm(N * ncol(LAM_E)), N, ncol(LAM_E)) %*% chol(PHI_E)
ystar <- fac %*% t(LAM_E) +
  matrix(rnorm(N * p_E), N, p_E) %*% diag(sqrt(1 - communality), p_E, p_E)
Y_E <- vapply(seq_len(p_E), function(j) findInterval(ystar[, j], thr_E[[j]]),
              integer(N))
colnames(Y_E) <- sprintf("y%d", seq_len(p_E))

tau_E    <- unlist(thr_E)
q_E      <- pnorm(tau_E)
var_id_E <- rep(seq_len(p_E), lengths(thr_E))
G_E      <- gamma_tau(tau_E, q_E, var_id_E, R_E)
lab_E    <- threshold_labels(lengths(thr_E))
dimnames(G_E) <- list(lab_E, lab_E)

f_data <- file.path("fixtures", "data",     "case_E.csv")
f_thr  <- file.path("fixtures", "expected", "case_E_thresholds.csv")
f_gam  <- file.path("fixtures", "expected", "case_E_gamma_tau.csv")

write_csv_lf(Y_E, f_data)
write_csv_lf(data.frame(label = lab_E, tau = tau_E, q = q_E), f_thr)
write_csv_lf(cbind(label = lab_E, as.data.frame(G_E)), f_gam)

written <- c(written, f_data, f_thr, f_gam)

stopifnot(isTRUE(all.equal(G_E, t(G_E))))
stopifnot(min(eigen(G_E, only.values = TRUE)$values) > 0)

Z_E  <- cumulative_indicators(Y_E, thr_E)
qh_E <- colMeans(Z_E)
ps_E <- sweep(Z_E, 2, qh_E) / rep(dnorm(qnorm(qh_E)), each = N)
w_E  <- outer(var_id_E, var_id_E, "==")
gap_E <- max(abs((crossprod(ps_E) / N -
                  gamma_tau(qnorm(qh_E), qh_E, var_id_E, R_E))[w_E]))
stopifnot(gap_E < 1e-10)

# Every category must be occupied in all six indicators, or Module 1 aborts on
# this fixture for a reason that has nothing to do with what it is for.
stopifnot(all(apply(Y_E, 2, function(v) length(unique(v))) == length(PI3)))

cat(sprintf("case E: p=%d, T=%d, within-variable identity gap %.2e  [%s]\n",
            p_E, length(tau_E), gap_E,
            "two factors, phi=0.5; carries the structural-path fits"))

# ---------------------------------------------------------------------------
# Two near-duplicate pair tables and their tier-1 rho-hat (Module 2 near +/-1).
#
# Pair (y1, y2) of Monte Carlo replicates mc_c_98_0007 and mc_c_98_0002 (design
# c, rho = .98, N = 5000). Production aborted on the first and stalled on the
# second (E-055, "Module 2 near +/-1"). Their empty corner cells
# have true probabilities below 2^-55 near the root, so every rectangle
# difference computes them as 0 or +/-2^-55; the root is reachable only by a
# likelihood summed over OBSERVED cells (derivation @sec-emptycells).
#
# The COUNTS are embedded, not read: the replicate .dta files are gitignored.
# Thresholds are Module 1's arithmetic (dwlsmv_cumprop: cumulative count / N).
# rho-hat is the derivation's rk_rho() (chunk rank-core), restated here on the
# table rather than on the two variables; cellP and dP_drho are the
# derivation's gate3-core, restated. Do not change any of the three here
# without changing it there. Nothing asserts that the two agree: they were
# compared once, on 2026-09-26, and this file's 17 digits equalled those
# printed by the derivation's chunk ec-L7 (@sec-emptycells) exactly.
#
# Numbers are written with sprintf("%.17g") BEFORE write_csv_lf(), because
# write.csv() formats doubles to 15 significant digits.
# ---------------------------------------------------------------------------

hr_cellP <- function(a, b, rho) {
  A <- c(-Inf, a, Inf); B <- c(-Inf, b, Inf)
  C <- outer(seq_along(A), seq_along(B), Vectorize(function(r,s) Phi2(A[r],B[s],rho)))
  C[-1,-1] - C[-nrow(C),-1] - C[-1,-ncol(C)] + C[-nrow(C),-ncol(C)]
}
hr_phi2s <- function(x, y, rho) if (!is.finite(x) || !is.finite(y)) 0 else phi2(x, y, rho)
hr_dP_drho <- function(a, b, rho) {
  A <- c(-Inf, a, Inf); B <- c(-Inf, b, Inf)
  outer(seq_len(length(a)+1), seq_len(length(b)+1), Vectorize(function(r,s)
    hr_phi2s(A[r+1],B[s+1],rho) - hr_phi2s(A[r],B[s+1],rho) -
    hr_phi2s(A[r+1],B[s],  rho) + hr_phi2s(A[r],B[s],  rho)))
}
hr_rho <- function(n, a, b) {                  # rk_rho(), on the table
  S <- function(r) { v <- n * hr_dP_drho(a, b, r) / hr_cellP(a, b, r); sum(v[n > 0]) }
  lo <- -0.5; hi <- 0.5
  while (isTRUE(S(hi) > 0)) { nx <- 1 - (1 - hi) / 2; if (!is.finite(S(nx)) || nx > 1 - 1e-7) break; hi <- nx }
  while (isTRUE(S(lo) < 0)) { nx <- -1 + (1 + lo) / 2; if (!is.finite(S(nx)) || nx < 1e-7 - 1) break; lo <- nx }
  r <- uniroot(S, c(lo, hi), tol = 1e-14)$root
  for (it in 1:25) {
    P <- hr_cellP(a, b, r); Dr <- hr_dP_drho(a, b, r)
    step <- (S(r) / sum(n)) / sum((Dr^2 / P)[P > 0])
    if (abs(step) < 1e-15) break
    r <- r + step
  }
  stopifnot(abs(step) < 1e-15, abs(r) < 0.999)  # met the stopping rule, interior
  r
}

HR_TABLES <- list(
  mc_c_98_0007 = matrix(c(679,   92,    0,   0,
                           81, 1523,  152,   0,
                            0,  169, 1226, 102,
                            0,    0,   91, 885), 4, byrow = TRUE),
  mc_c_98_0002 = matrix(c(665,   88,    0,   0,
                          100, 1456,  167,   0,
                            0,  152, 1275, 101,
                            0,    0,  116, 880), 4, byrow = TRUE))

hr_thr <- function(m) qnorm(cumsum(m)[-length(m)] / sum(m))
g17    <- function(x) sprintf("%.17g", x)

hr_rows <- lapply(names(HR_TABLES), function(nm) {
  n <- HR_TABLES[[nm]]
  stopifnot(sum(n) == 5000, all(rowSums(n) > 0), all(colSums(n) > 0))
  a <- hr_thr(rowSums(n)); b <- hr_thr(colSums(n)); r <- hr_rho(n, a, b)
  cnt <- setNames(as.list(as.vector(t(n))),                     # row-major: n_<y1>_<y2>
                  sprintf("n_%d_%d", rep(0:3, each = 4), rep(0:3, 4)))
  cat(sprintf("%s: rho-hat %s\n", nm, g17(r)))
  data.frame(table = nm, N = sum(n), cnt,
             a1 = g17(a[1]), a2 = g17(a[2]), a3 = g17(a[3]),
             b1 = g17(b[1]), b2 = g17(b[2]), b3 = g17(b[3]),
             rho_hat = g17(r))
})

f_hr <- file.path("fixtures", "expected", "highrho_polychoric.csv")
write_csv_lf(do.call(rbind, hr_rows), f_hr)
written <- c(written, f_hr)

write_manifest(file.path("fixtures", "MANIFEST.json"), written,
               seeds = sprintf("%d + utf8ToInt(case)", SEED))

cat(sprintf("\n%d files written; manifest at fixtures/MANIFEST.json\n", length(written)))
