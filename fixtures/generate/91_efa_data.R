# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 1: data for the geomin EFA fixtures
#
# 0.3.0 EFA with oblique geomin rotation, Unit 4 (design approved by the user
# 2026-10-29, before any data were generated). Two cases no existing data set
# covers:
#
#   Case F  nine 4-category indicators, THREE factors, two cross-loadings:
#           the m = 3 case. Primary loadings .7 .7 .6 per factor; y3 also
#           loads .3 on F2 and y6 .3 on F3, so the rotation has real work
#           to do; factor correlations .3.
#   Case G  six BINARY indicators, two factors, loadings .7, phi_12 = .4:
#           the all-binary case (df = 4 at m = 2).
#
# Written as its own generator with its own manifest, after 95_'s precedent,
# so MANIFEST.json and Cases A-E cannot move. The simulation is 10_'s Case E
# block (delta parameterization, var(y*) = 1, cut at population thresholds).
#
#   fixtures/data/case_F.csv, case_G.csv   raw ordinal data, integer-valued
#   fixtures/MANIFEST_efa_data.json        seeds and population values
#
# Run from the repository root:  Rscript fixtures/generate/91_efa_data.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

SEED <- 20261029
N    <- 5000
MAN  <- file.path("fixtures", "MANIFEST_efa_data.json")

PI4 <- c(0.15, 0.35, 0.30, 0.20)   # as 10_
PI2 <- c(0.30, 0.70)
thr_of <- function(pmf) qnorm(cumsum(pmf)[-length(pmf)])

LAM_F <- matrix(0, 9, 3)
LAM_F[1:3, 1] <- c(0.7, 0.7, 0.6)
LAM_F[4:6, 2] <- c(0.7, 0.7, 0.6)
LAM_F[7:9, 3] <- c(0.7, 0.7, 0.6)
LAM_F[3, 2] <- 0.3
LAM_F[6, 3] <- 0.3
PHI_F <- matrix(0.3, 3, 3); diag(PHI_F) <- 1

LAM_G <- matrix(c(rep(0.7, 3), rep(0, 3),
                  rep(0, 3), rep(0.7, 3)), ncol = 2)
PHI_G <- matrix(c(1, 0.4, 0.4, 1), 2, 2)

# The intended pattern, stated independently of the assignments above: the
# primary block structure plus exactly the two cross-loadings (y3, F2) and
# (y6, F3). And Phi must be a positive definite correlation matrix.
pat_F <- matrix(FALSE, 9, 3)
pat_F[cbind(1:9, rep(1:3, each = 3))] <- TRUE
pat_F[cbind(c(3, 6), c(2, 3))] <- TRUE
stopifnot(identical(LAM_F != 0, pat_F), sum(pat_F) == 11,
          all(LAM_F[cbind(c(3, 6), c(2, 3))] == 0.3),
          all(diag(PHI_F) == 1), isSymmetric(PHI_F),
          min(eigen(PHI_F, symmetric = TRUE, only.values = TRUE)$values) > 0,
          min(eigen(PHI_G, symmetric = TRUE, only.values = TRUE)$values) > 0)

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

Y_F <- sim_case("F", LAM_F, PHI_F, rep(list(thr_of(PI4)), 9))
Y_G <- sim_case("G", LAM_G, PHI_G, rep(list(thr_of(PI2)), 6))

f_F <- file.path("fixtures", "data", "case_F.csv")
f_G <- file.path("fixtures", "data", "case_G.csv")
write_csv_lf(Y_F, f_F)
write_csv_lf(Y_G, f_G)

vec <- function(M) paste(format(as.vector(M)), collapse = " ")
write_manifest(MAN, c(f_F, f_G),
  seeds = c(F = SEED + utf8ToInt("F"), G = SEED + utf8ToInt("G")),
  extra = list(packages = c(
    "none (base R only)",
    paste("N", N),
    paste("case F thresholds (all items):", vec(round(thr_of(PI4), 6))),
    paste("case F Lambda (column-major 9x3):", vec(LAM_F)),
    paste("case F Phi (column-major 3x3):", vec(PHI_F)),
    paste("case G thresholds (all items):", vec(round(thr_of(PI2), 6))),
    paste("case G Lambda (column-major 6x2):", vec(LAM_G)),
    paste("case G Phi (column-major 2x2):", vec(PHI_G)))))

cat("case_F.csv, case_G.csv and", MAN, "written.\n")
