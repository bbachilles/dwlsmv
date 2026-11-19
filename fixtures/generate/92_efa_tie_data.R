# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 1: data for the A4 near-tie fixtures
#
# Drafted by a read-only agent 2026-11-17; executed by the main session 2026-11-18.
#
# Release A, A4 near-tie question (design and acceptance criterion declared
# 2026-11-17 before any data: E-109, A4 near-tie criterion): at a
# tie in the column-ordering statistic, do dwlsmv and lavaan post the same
# column order? Two cases, p = 8 four-category items, m = 2 factors:
#
#   Case T  exact tie.   f1 loads items {1,4,5,8}, f2 loads items {2,3,6,7},
#           all eight primary loadings .70, no cross-loadings, phi_12 .3.
#           Both columns have mean row index 4.5 by construction.
#   Case B  boundary.    As T, but y8 on f1 is .56 = 0.8 * .70, exactly the
#           population cutoff of the ordering statistic (0.8 * max|loading|).
#
# The simulation is 91_'s sim_case (delta parameterization, var(y*) = 1, cut at
# the population thresholds, PI4 as 10_ and 91_). SEED is used as
# SEED + utf8ToInt(case name), the 91_ convention. No seed or design is
# changed after data exist (criterion.md).
#
# Own generator and manifest (91_'s precedent), so no existing manifest moves.
#
#   fixtures/data/efatie_T.csv, efatie_B.csv   raw ordinal data, integer-valued
#   fixtures/MANIFEST_efa_tie_data.json    seeds and population values
#
# Run from the repository root:  Rscript fixtures/generate/92_efa_tie_data.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

SEED <- 20261117
N    <- 5000
MAN  <- file.path("fixtures", "MANIFEST_efa_tie_data.json")

PI4 <- c(0.15, 0.35, 0.30, 0.20)   # as 10_ and 91_
thr_of <- function(pmf) qnorm(cumsum(pmf)[-length(pmf)])

ITEMS_F1 <- c(1, 4, 5, 8)
ITEMS_F2 <- c(2, 3, 6, 7)

LAM_T <- matrix(0, 8, 2)
LAM_T[ITEMS_F1, 1] <- 0.70
LAM_T[ITEMS_F2, 2] <- 0.70
LAM_B <- LAM_T
LAM_B[8, 1] <- 0.56                # = 0.8 * .70, the population cutoff
PHI <- matrix(c(1, 0.3, 0.3, 1), 2, 2)

# The intended pattern, stated independently of the assignments above.
pat <- matrix(FALSE, 8, 2)
pat[cbind(c(1, 4, 5, 8), 1)] <- TRUE
pat[cbind(c(2, 3, 6, 7), 2)] <- TRUE
stopifnot(identical(LAM_T != 0, pat), identical(LAM_B != 0, pat),
          sum(pat) == 8,
          # Case T: all eight primaries .70, nothing else
          all(LAM_T[pat] == 0.70), all(LAM_T[!pat] == 0),
          # Case B: identical to T except y8 on f1
          identical(LAM_B[-8, ], LAM_T[-8, ]), LAM_B[8, 2] == 0,
          # y8's loading is 0.8 * the column maximum (double arithmetic: 0.8*0.7
          # differs from 0.56 by 1.1e-16, so equality is to 1e-15, not identical)
          abs(LAM_B[8, 1] - 0.8 * max(abs(LAM_B[, 1]))) < 1e-15,
          max(abs(LAM_B[, 1])) == 0.70,
          # every other primary in B is strictly above the cutoff and every
          # zero strictly below it (the cutoff admits exactly the intended set)
          all(LAM_B[setdiff(ITEMS_F1, 8), 1] > 0.8 * 0.70 + 0.1),
          all(LAM_B[ITEMS_F2, 2] == 0.70),
          # mean row index of the primaries, by construction
          mean(ITEMS_F1) == 4.5, mean(ITEMS_F2) == 4.5,
          all(diag(PHI) == 1), isSymmetric(PHI), PHI[1, 2] == 0.3,
          min(eigen(PHI, symmetric = TRUE, only.values = TRUE)$values) > 0)

sim_case <- function(nm, LAM, PHI, thr) {
  p <- nrow(LAM)
  comm <- diag(LAM %*% PHI %*% t(LAM))
  stopifnot(all(comm < 1))
  set.seed(SEED + utf8ToInt(nm))
  fac   <- matrix(rnorm(N * ncol(LAM)), N, ncol(LAM)) %*% chol(PHI)
  ystar <- fac %*% t(LAM) +
    matrix(rnorm(N * p), N, p) %*% diag(sqrt(1 - comm), p, p)
  Y <- vapply(seq_len(p), function(j) findInterval(ystar[, j], thr[[j]]),
              integer(N))
  colnames(Y) <- sprintf("y%d", seq_len(p))
  # every category observed in every item, or the case is not what it says
  stopifnot(all(vapply(seq_len(p), function(j)
    length(unique(Y[, j])) == length(thr[[j]]) + 1, logical(1))))
  cat(sprintf("case %s: p=%d, N=%d, max communality %.3f, smallest category count %d\n",
              nm, p, N, max(comm), min(vapply(seq_len(p), function(j)
                min(tabulate(Y[, j] + 1)), integer(1)))))
  Y
}

# Seeds must be distinct: "T" and "B" give different offsets.
stopifnot(utf8ToInt("T") != utf8ToInt("B"))
Y_T <- sim_case("T", LAM_T, PHI, rep(list(thr_of(PI4)), 8))
Y_B <- sim_case("B", LAM_B, PHI, rep(list(thr_of(PI4)), 8))

f_T <- file.path("fixtures", "data", "efatie_T.csv")
f_B <- file.path("fixtures", "data", "efatie_B.csv")
write_csv_lf(Y_T, f_T)
write_csv_lf(Y_B, f_B)

vec <- function(M) paste(format(as.vector(M)), collapse = " ")
write_manifest(MAN, c(f_T, f_B),
  seeds = c(T = SEED + utf8ToInt("T"), B = SEED + utf8ToInt("B")),
  extra = list(packages = c(
    "none (base R only)",
    paste("N", N),
    "design and criterion: docs/notes/efa-se/a4-neartie/criterion.md (declared 2026-11-17 before any data)",
    paste("thresholds (all items, both cases):", vec(round(thr_of(PI4), 6))),
    paste("case T Lambda (column-major 8x2):", vec(LAM_T)),
    paste("case B Lambda (column-major 8x2):", vec(LAM_B)),
    paste("Phi (column-major 2x2, both cases):", vec(PHI)),
    "f1 loads items 1 4 5 8; f2 loads items 2 3 6 7; case B y8 on f1 = .56 = 0.8 * .70")))

cat("efatie_T.csv, efatie_B.csv and", MAN, "written.\n")
