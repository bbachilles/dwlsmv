# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- tier 1: data for the second-order fixtures
#
# Release B, unit B4 (E-108; design in
# E-114 (the B4 design), to be approved by the owner before
# any data are generated). Three cases no existing data set covers:
#
#   Case H  one second-order factor G over F1-F4, three 4-category items per
#           F (y1-y12), and one item y13 loading on G directly. Standardized
#           second-order loadings b* = .8 .7 .6 .75; item loadings .7 .7 .6;
#           y13 loads .6 on G. Serves three models: G over F1-F4 (y1-y12; the
#           second level has 2 df), G over F1-F3 (y1-y9; the region where it
#           equals the correlated model, SO2.3), and B0's MD (y1-y13).
#           Dropping F4 or y13 from a population of this form leaves a
#           population of the same form, so each submodel is correct.
#   Case I  a CORRELATED three-factor population with factor correlations
#           r12 = r13 = .8, r23 = .4 (B0's example, E-114),
#           items .7 .7 .6. Fitted as G over F1-F3 it has no proper solution:
#           b*_1^2 = r12 r13 / r23 = 1.6, so psi_1 < 0 in metric V and metric
#           U cannot represent it (SO3.3, SO4.4). The improper case.
#   Case J  two correlated second-order factors, G1 over F1-F3 and G2 over
#           F4-F6, corr(G1, G2) = .5, b* = .8 .7 .6 under each, items .7 .7
#           .6 (y1-y18). B0's M33.
#
# Written as its own generator with its own manifest, after 91_'s precedent,
# so MANIFEST.json and Cases A-G cannot move. The simulation is 91_'s
# (delta parameterization, var(y*) = 1, cut at population thresholds); the
# factor covariance over ALL latent variables, G included, is built from the
# structural model Phi = (I - B)^-1 Psi (I - B)^-T in standardized form.
#
#   fixtures/data/case_H.csv, case_I.csv, case_J.csv   raw ordinal data
#   fixtures/MANIFEST_so_data.json                     seeds and population
#
# Run from the repository root:  Rscript fixtures/generate/93_so_data.R
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

source(file.path("fixtures", "generate", "00_common.R"))

SEED <- 20261119
N    <- 1000
MAN  <- file.path("fixtures", "MANIFEST_so_data.json")

PI4 <- c(0.15, 0.35, 0.30, 0.20)   # as 10_ and 91_
thr_of <- function(pmf) qnorm(cumsum(pmf)[-length(pmf)])
LAM3 <- c(0.7, 0.7, 0.6)           # item loadings of every first-order factor

# Standardized structural covariance over (G's, F's): G's with unit variance
# and correlation matrix RG; F_k = b*_k G_{g(k)} + disturbance with variance
# 1 - b*_k^2, so every F has unit total variance.
so_phi <- function(bstar, g_of, RG) {
  ng <- nrow(RG); K <- length(bstar); m <- ng + K
  B <- matrix(0, m, m)
  B[cbind(ng + seq_len(K), g_of)] <- bstar
  Psi <- diag(c(rep(1, ng), 1 - bstar^2), m)
  Psi[seq_len(ng), seq_len(ng)] <- RG
  A <- solve(diag(m) - B)
  P <- A %*% Psi %*% t(A)
  stopifnot(max(abs(diag(P) - 1)) < 1e-14)
  P
}
# Item loadings: three items per F (LAM3), F's after the G's; optional direct
# items on a G appended after them.
so_lam <- function(ng, K, direct = NULL) {
  L <- matrix(0, 3 * K, ng + K)
  for (k in seq_len(K)) L[3 * (k - 1) + 1:3, ng + k] <- LAM3
  if (!is.null(direct)) for (d in direct)
    L <- rbind(L, replace(numeric(ng + K), d$g, d$lambda))
  L
}

# Case H: G, F1-F4; y13 on G.
B_H   <- c(0.8, 0.7, 0.6, 0.75)
PHI_H <- so_phi(B_H, rep(1, 4), matrix(1))
LAM_H <- so_lam(1, 4, list(list(g = 1, lambda = 0.6)))

# Case I: correlated three-factor population (no G).
PHI_I <- matrix(c(1, .8, .8,
                  .8, 1, .4,
                  .8, .4, 1), 3, 3)
LAM_I <- so_lam(0, 3)
# the improper point, stated independently: b*_1^2 = r12 r13 / r23 > 1
stopifnot(PHI_I[1, 2] * PHI_I[1, 3] / PHI_I[2, 3] > 1)

# Case J: G1 over F1-F3, G2 over F4-F6, corr(G1, G2) = .5.
B_J   <- rep(c(0.8, 0.7, 0.6), 2)
RG_J  <- matrix(c(1, .5, .5, 1), 2, 2)
PHI_J <- so_phi(B_J, rep(1:2, each = 3), RG_J)
LAM_J <- so_lam(2, 6)

for (P in list(PHI_H, PHI_I, PHI_J))
  stopifnot(isSymmetric(P), all(abs(diag(P) - 1) < 1e-14),
            min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) > 0)
# the intended structure, checked against the population: in H and J each
# F-F correlation under one G is b*_k b*_l; across G's in J, b*_k b*_l r_G
stopifnot(abs(PHI_H[2, 3] - .8 * .7) < 1e-14, abs(PHI_H[1, 5] - .75) < 1e-14,
          abs(PHI_J[3, 4] - .8 * .7) < 1e-14,
          abs(PHI_J[3, 6] - .8 * .8 * .5) < 1e-14)

sim_case <- function(nm, LAM, PHI, thr) {        # 91_'s, unchanged
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
  stopifnot(all(vapply(seq_len(p), function(j)
    length(unique(Y[, j])) == length(thr[[j]]) + 1, logical(1))))
  cat(sprintf("case %s: p=%d, N=%d, max communality %.3f, smallest category count %d\n",
              nm, p, N, max(comm), min(vapply(seq_len(p), function(j)
                min(tabulate(Y[, j] + 1)), integer(1)))))
  Y
}

thr4 <- function(p) rep(list(thr_of(PI4)), p)
Y_H <- sim_case("H", LAM_H, PHI_H, thr4(nrow(LAM_H)))
Y_I <- sim_case("I", LAM_I, PHI_I, thr4(nrow(LAM_I)))
Y_J <- sim_case("J", LAM_J, PHI_J, thr4(nrow(LAM_J)))

f <- file.path("fixtures", "data", sprintf("case_%s.csv", c("H", "I", "J")))
write_csv_lf(Y_H, f[1])
write_csv_lf(Y_I, f[2])
write_csv_lf(Y_J, f[3])

vec <- function(M) paste(format(as.vector(M)), collapse = " ")
write_manifest(MAN, f,
  seeds = c(H = SEED + utf8ToInt("H"), I = SEED + utf8ToInt("I"),
            J = SEED + utf8ToInt("J")),
  extra = list(packages = c(
    "none (base R only)",
    paste("N", N),
    paste("thresholds (all items):", vec(round(thr_of(PI4), 6))),
    paste("case H b* (G -> F1..F4):", vec(B_H), "; y13 on G: 0.6; items per F:", vec(LAM3)),
    paste("case H Phi over (G, F1..F4), column-major:", vec(round(PHI_H, 6))),
    paste("case I Phi over (F1..F3), column-major:", vec(PHI_I),
          "; implied b*_1^2 =", PHI_I[1, 2] * PHI_I[1, 3] / PHI_I[2, 3]),
    paste("case J b* (G1 -> F1..F3, G2 -> F4..F6):", vec(B_J), "; corr(G1, G2) 0.5"),
    paste("case J Phi over (G1, G2, F1..F6), column-major:", vec(round(PHI_J, 6))))))

cat("case_H.csv, case_I.csv, case_J.csv and", MAN, "written.\n")
