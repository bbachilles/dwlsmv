# ---------------------------------------------------------------------------
# dwlsmv fixture generation -- the marker-metric start trap (E-016,
# "Marker-metric start trap and recovery through unit variance", 2026-10-06)
#
# Executed 2026-10-06 (haven): 200 replicates per condition, smallest
# category count 29 (P) and 31 (C). Results in E-016.
#
# Replicates for the trap-rate study. The design was fixed and approved by the
# user on 2026-10-06 BEFORE any data were generated:
#
#   A (MathAtt) ~ N(0, 1);  B (MathAb) = 0.55 A + zeta,  zeta ~ N(0, 1)
#   y* = lambda * factor + e, var(y*) = 1 (delta parameterization), cut at
#   the thresholds below. Population values are dwlsmv's unit-variance fit
#   to Stata's gsem_cfa example, rounded to two decimals; the Stata-Press
#   data themselves are not used.
#
#   P  primary:  att3, att4 load NEGATIVELY (-0.76, -0.54)
#   C  control:  the same with those two loadings' signs flipped (+0.76, +0.54),
#                thresholds unchanged, so only the correlations change sign
#
# N = 500, 200 replicates per condition. The trap definition, the rate report
# and the rule choosing the committed fixture (the first P replicate, in seed
# order, that meets the definition) are in E-016 and were fixed
# before generation. Fitting is tools/marker_trap_study.do (Stata), not here.
#
# Output goes to fixtures/replication/marker_trap_{P,C}/, which .gitignore
# excludes, with one index.csv per condition. Run from the repository root:
#     Rscript fixtures/generate/95_marker_trap.R
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek
# ---------------------------------------------------------------------------

SEED_BASE <- 20261006L
NREP <- 200L
N <- 500L
OUT <- file.path("fixtures", "replication")

lam_q   <- c(0.59, 0.27, 0.43, 0.21, 0.52, 0.43, 0.46, 0.43)
thr_q   <- c(-0.02, 0.27, -0.09, 0.19, 0.03, 0.17, -0.05, 0.02)
thr_att <- list(c(-0.52, -0.11, 0.15, 0.64),
                c(-0.62, -0.12, 0.21, 0.66),
                c(-0.41, -0.02, 0.22, 0.64),
                c(-0.55, -0.11, 0.21, 0.71),
                c(-0.70, -0.19, 0.12, 0.52))
BETA <- 0.55

CONDITIONS <- list(
  P = c(0.56, 0.25, -0.76, -0.54, 0.36),
  C = c(0.56, 0.25,  0.76,  0.54, 0.36)
)

# Structural draw, then latent responses with unit total variance. Factor A is
# exogenous with variance 1; B's total variance is BETA^2 + 1.
simulate_trap <- function(n, lam_att) {
  A <- rnorm(n)
  B <- BETA * A + rnorm(n)
  varB <- BETA^2 + 1
  comm <- c(lam_q^2 * varB, lam_att^2)
  if (any(comm >= 1)) stop("a communality reached 1")
  fac  <- cbind(matrix(B, n, 8), matrix(A, n, 5))
  lam  <- c(lam_q, lam_att)
  ystar <- sweep(fac, 2, lam, `*`) +
           sweep(matrix(rnorm(n * 13), n, 13), 2, sqrt(1 - comm), `*`)
  thr <- c(as.list(thr_q), thr_att)
  out <- vapply(seq_len(13), function(j) findInterval(ystar[, j], thr[[j]]),
                integer(n))
  colnames(out) <- c(sprintf("q%d", 1:8), sprintf("att%d", 1:5))
  as.data.frame(out)
}

# --- fixture mode ---------------------------------------------------------------
# DWLSMV_TRAP_FIXTURE=1 writes the TIER-1 FIXTURE instead of the replicates:
# P replicate 1 (seed SEED_BASE + 1000000 + 1), the replicate the pre-set rule
# chose (E-016, trap-rate study), as
#   fixtures/data/marker_trap.csv                 q1..q8, att1..att5 (0-based codes)
#   fixtures/expected/marker_trap_population.csv  label, value: the population
#       parameters in both metrics (unit variance as generated; marker, att1 and
#       q1 at 1, by c_f = lambda of the marker:
#       lambda/c_f, beta_BA * c_B / c_A, psi_ff * c_f^2)
# with fixtures/MANIFEST_marker_trap.json. If the replicate file exists it must
# hold identical codes, so the fixture is provably the studied replicate.
if (Sys.getenv("DWLSMV_TRAP_FIXTURE") == "1") {
  source(file.path("fixtures", "generate", "00_common.R"))
  seed <- SEED_BASE + 1L * 1000000L + 1L
  set.seed(seed)
  dat <- simulate_trap(N, CONDITIONS$P)
  rep1 <- file.path(OUT, "marker_trap_P", "marker_trap_P_0001.dta")
  if (file.exists(rep1)) {
    old <- as.data.frame(lapply(haven::read_dta(rep1), as.integer))
    stopifnot(identical(dim(old), dim(dat)), all(old == dat))
    cat("fixture codes identical to", rep1, "\n")
  }
  DATA <- file.path("fixtures", "data", "marker_trap.csv")
  POP  <- file.path("fixtures", "expected", "marker_trap_population.csv")
  write_csv_lf(dat, DATA)

  lam <- c(lam_q, CONDITIONS$P)
  nms <- c(sprintf("q%d", 1:8), sprintf("att%d", 1:5))
  fac <- c(rep("MathAb", 8), rep("MathAtt", 5))
  cA <- CONDITIONS$P[1]; cB <- lam_q[1]
  cc <- ifelse(fac == "MathAb", cB, cA)
  thr <- c(as.list(thr_q), thr_att)
  tl <- unlist(lapply(seq_along(thr), function(j)
    sprintf("%s|t%d", nms[j], seq_along(thr[[j]]))))
  pop <- data.frame(
    label = c(tl,
              paste0("uv:", fac, "=~", nms), "uv:MathAb~MathAtt",
              "uv:MathAtt~~MathAtt", "uv:MathAb~~MathAb",
              paste0("marker:", fac, "=~", nms), "marker:MathAb~MathAtt",
              "marker:MathAtt~~MathAtt", "marker:MathAb~~MathAb"),
    value = c(unlist(thr),
              lam, BETA, 1, 1,
              lam / cc, BETA * cB / cA, cA^2, cB^2))
  write_csv_lf(pop, POP)
  write_manifest(file.path("fixtures", "MANIFEST_marker_trap.json"), c(DATA, POP),
    seeds = paste("set.seed", seed, "(P replicate 1, chosen by the pre-set rule)"),
    extra = list(packages = c(paste("haven", as.character(packageVersion("haven"))),
                              R.version.string)))
  cat("tier-1 marker-trap fixture written:", DATA, POP, "\n")
  quit(save = "no")
}

for (ci in seq_along(CONDITIONS)) {
  nm <- names(CONDITIONS)[ci]
  cdir <- file.path(OUT, paste0("marker_trap_", nm))
  dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
  idx <- data.frame()
  for (r in seq_len(NREP)) {
    seed <- SEED_BASE + ci * 1000000L + r
    set.seed(seed)
    dat <- simulate_trap(N, CONDITIONS[[nm]])
    f <- file.path(cdir, sprintf("marker_trap_%s_%04d.dta", nm, r))
    haven::write_dta(dat, f)
    cnt <- unlist(lapply(dat, function(v) tabulate(v + 1L, max(v) + 1L)))
    idx <- rbind(idx, data.frame(condition = nm, replicate = r, seed = seed,
                                 min_cat = min(cnt),
                                 ncat_att = min(sapply(dat[9:13], function(v) length(unique(v)))),
                                 file = f, stringsAsFactors = FALSE))
  }
  write.csv(idx, file.path(cdir, "index.csv"), row.names = FALSE)
  cat(sprintf("%s: %d replicates, N = %d, smallest category count %d\n",
              nm, NREP, N, min(idx$min_cat)))
}
