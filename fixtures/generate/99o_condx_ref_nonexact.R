# ---------------------------------------------------------------------------
# Executed 2026-11-14 by the main session (R, from the repository root); drafted by a read-only agent.
#
# dwlsmv fixture generation -- tier 1: a NEW dataset (condx_U) and the analytic
# base-R reference for conditional-mode fit and inference on data that DO NOT
# FIT EXACTLY. One binary x, one group, six ordinal items. This is a development
# fixture; it is NOT a dataset of the conditional-mode preregistration
# (E-098, frozen 2026-11-10): link 2 uses condx_A
# and condx_D, link 3 uses condx_link3.csv. Nothing here reads, compares against
# or reuses those datasets or condx_ref_A_*/condx_ref_D_*; the generator reads
# only code (the derivation's functions, the frozen run script's Gauss-Newton
# loop, 00_common.R).
#
# WHY. tests/unit/test_condx_fit.do covers exact-fit populations only. At an
# exact fit F = 0 and the residual is zero, so the slope-row weights in W_D,
# T_raw = (N - 1) F, and the slope rows of the inference are invisible. Here the
# fitted model omits a true cross-loading, so F > 0 and each of them moves a
# number.
#
# DATA (seeded, base R). N = 600. K = (3, 2, 4, 3, 2, 3). Two factors, F1 ->
# y1-y3, F2 -> y4-y6, AND a cross-loading y4 <- F1 (0.20), lambda_4 = 0.65; disturbance
# covariance psi_12 = 0.30; F_k = gamma_k x + zeta_k, gamma = (0.5, -0.4); a
# direct effect of x on y3, kappa_3 = 0.4; P(x = 1) = 0.4; Var(y* | x) = 1
# (delta parameterization). The FITTED model omits the cross-loading:
#   F1 =~ y1 y2 y3, F2 =~ y4 y5 y6, unitvariance F1 F2, psi F1 F2 free,
#   covariate x, gamma F1 x free, gamma F2 x free, kappa y3 x free.
# theta (spec 4.8 order tau, lambda, beta, gamma, kappa, psi) =
#   11 thresholds | 6 loadings | gamma_1, gamma_2 | kappa_3 | psi_12   (21 free)
# s has 11 + 6 + 15 = 32 rows, so df = 11.
# Constants: N, K, thresholds, loadings, gammas, kappa, P(x=1), seed fixed before the first run;
# the cross-loading 0.20, lambda_4 0.65 and psi_12 0.30 were CHOSEN AFTER MEASURING (a first dry run
# with 0.30, 0.70, 0.35 gave lambda_4-hat above 1, an inadmissible solution; three settings were then
# tried and this one kept: all loading estimates below 0.9, T_raw above 2 df). DESIGN FLOORS asserted
# below (every category in each x group observed at least 5 times; max |rho-hat|
# < 0.9; T_raw >= 2 df) are floors on the dataset, not on any estimate of interest.
# If a floor fails the generator stops and nothing is written.
#
# REFERENCE. No mathematics is defined in the Stage 1 or Gamma path. Stage 1
# (cx_stage1), the Hessian-form influence function (cxh_psi), Gamma (cx_gamma)
# come from the CURRENT spec/derivations/dwlsmv-derivations.qmd (chunks cx-setup,
# cx-stage1, cxn-psi, cxh-forms), exactly as 99n_condx_ref.R runs them; they read
# the globals cx_K, cx_p, cx_pairs, which are reassigned below for six items
# (no chunk is edited). The model map and its Jacobian use chunk cx-anchor-setup's
# ca_mod / ca_s / ca_delta (a multi-factor, kappa-carrying reduced form, central
# differences at h = 1e-6; the map is a polynomial of degree 3 in theta, so the
# difference error is rounding only, and is checked below against h = 1e-4).
# Gauss-Newton is cxs_gn from the simulation run script of E-065 (definitions
# GN_TOL, GN_MAX, cxs_gn only), which calls the globals cx_sig_th and cx_delta:
# they are defined HERE for this model. The pipeline after that copies
# 99n_condx_ref.R's ref_one lines for sandwich, T and T3 (99n is a script, not a
# library, and is hard-wired to y1-y3, so its driver cannot be sourced).
#   W_D = diag(1 / diag(Gamma-hat)); SE = sqrt(diag(B D'W G W D B) / (N - 1));
#   T_raw = (N - 1) r'W r; a, b from U Gamma; T3 = a T + b; df = 32 - 21.
# The start of the Gauss-Newton is the reference's own (ref_start below: tau =
# s_tau; loadings from the one-factor triad inside each factor block; gamma from
# the slopes over the loadings; kappa 0; psi_12 from the cross-block
# correlations). It is independent of production's start (zeros for gamma and
# kappa, the table's starts otherwise).
#
# Outputs (OUTROOT = fixtures by default):
#   data/condx_U.csv                         y1 .. y6, x
#   expected/condx_ref_U_s.csv               label, value        s in spec 4.2 order
#   expected/condx_ref_U_gamma.csv           row, col, value     Gamma-hat, lower triangle
#                                            with diagonal, column-major (as 99n)
#   expected/condx_ref_U_theta.csv           label, est, se      theta in the order above
#   expected/condx_ref_U_scalars.csv         quantity, value     N, npar, df, F, T_raw, a,
#                                            b, T3, p_T3, gn_iter, convergence diagnostics
#   MANIFEST_condx_ref_U.json
# Numbers %.17g. Run from the repository root:
#   Rscript fixtures/generate/99o_condx_ref_nonexact.R            # OUTROOT = fixtures
#   Rscript fixtures/generate/99o_condx_ref_nonexact.R <outroot>  # elsewhere (testing)
# fixtures/verify.R needs one more entry for MANIFEST_condx_ref_U.json.
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

args <- commandArgs(TRUE)
OUTROOT <- if (length(args) >= 1) args[1] else "fixtures"
QMD <- file.path("spec", "derivations", "dwlsmv-derivations.qmd")
RUN <- file.path("docs", "notes", "condx-sim", "condx_sim_run.R")
COMMON <- file.path("fixtures", "generate", "00_common.R")
CHUNKS <- c("cx-setup", "cx-stage1", "cxn-psi", "cxh-forms", "cx-anchor-setup")
RUN_DEFS <- c("GN_TOL", "GN_MAX", "cxs_gn")
CONV_S1 <- 1e-12   # as 99n_condx_ref.R (chosen there; carried over, not re-chosen on this data)
TAG <- "U"
SEED <- 20261114L

stopifnot(dir.exists(file.path(OUTROOT, "expected")), dir.exists(file.path(OUTROOT, "data")),
          file.exists(QMD), file.exists(RUN))
DAT <- file.path(OUTROOT, "data", sprintf("condx_%s.csv", TAG))
if (file.exists(DAT)) stop("99o: ", DAT, " exists; refusing to overwrite", call. = FALSE)

t_start <- proc.time()[["elapsed"]]

source(COMMON)
rm(Phi2, phi2)   # the .qmd's own definitions are used (as 99n)

# ---- the derivation's functions ---------------------------------------------
qmd <- readLines(QMD, encoding = "UTF-8")
chunk <- function(lab) {
  s <- grep(paste0("^```\\{r ", lab, "\\}"), qmd); stopifnot(length(s) == 1)
  e <- s + which(qmd[(s + 1):length(qmd)] == "```")[1]
  qmd[(s + 1):(e - 1)]
}
for (lab in CHUNKS) eval(parse(text = chunk(lab)), globalenv())
stopifnot(exists("cx_stage1"), exists("cxh_psi"), exists("cx_gamma"), exists("ca_mod"),
          exists("ca_s"), exists("ca_delta"), exists("Phi2"))
for (ex in parse(RUN, keep.source = FALSE)) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) &&
      as.character(ex[[2]]) %in% RUN_DEFS)
    eval(ex, globalenv())
}
stopifnot(all(vapply(RUN_DEFS, exists, NA)), GN_TOL == 1e-10, GN_MAX == 200L)

# ---- the design --------------------------------------------------------------
N <- 600L
q <- 0.4
K <- c(3L, 2L, 4L, 3L, 2L, 3L)
thr <- list(c(-0.6, 0.5), 0.1, c(-0.9, 0.0, 1.0), c(-0.4, 0.7), -0.2, c(-0.7, 0.4))
stopifnot(identical(lengths(thr) + 1L, K))
Lam <- matrix(0, 6, 2); Lam[1:3, 1] <- c(0.75, 0.65, 0.70); Lam[4:6, 2] <- c(0.65, 0.60, 0.75)
Lam[4, 1] <- 0.20                                   # the cross-loading the fitted model omits
Psi <- matrix(c(1, 0.30, 0.30, 1), 2)
gam <- c(0.5, -0.4)
kap <- c(0, 0, 0.4, 0, 0, 0)
res_var <- 1 - diag(Lam %*% Psi %*% t(Lam))         # Var(y* | x) = 1
stopifnot(all(res_var > 0.1))

set.seed(SEED)
x <- rbinom(N, 1L, q)
Z <- matrix(rnorm(2L * N), N, 2) %*% chol(Psi)      # rows ~ N(0, Psi)
Fac <- outer(x, gam) + Z
Y <- sapply(1:6, function(j) {
  ystar <- drop(Fac %*% Lam[j, ]) + kap[j] * x + rnorm(N, 0, sqrt(res_var[j]))
  findInterval(ystar, thr[[j]])
})
dat <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4], y5 = Y[, 5], y6 = Y[, 6], x = x)
# design floors: every (category, x) cell observed at least 5 times
for (j in 1:6) stopifnot(all(table(factor(dat[[j]], 0:(K[j] - 1)), dat$x) >= 5))

# ---- the reference's globals: six items, K as above --------------------------
cx_K <- K
cx_p <- length(cx_K)
cx_pairs <- t(combn(cx_p, 2))
nt <- sum(cx_K - 1)
npar <- nt + 10L
s_lab <- c(threshold_labels(cx_K - 1), sprintf("y%d~x", seq_len(cx_p)), correlation_labels(cx_p))
stopifnot(identical(correlation_labels(cx_p), sprintf("y%d~~y%d", cx_pairs[, 1], cx_pairs[, 2])))
th_lab <- c(threshold_labels(cx_K - 1), sprintf("F%d=~y%d", c(1, 1, 1, 2, 2, 2), 1:6),
            "F1~x", "F2~x", "y3~x", "F1~~F2")
stopifnot(length(th_lab) == npar)
f17 <- function(v) sprintf("%.17g", v)

# the fitted model, in chunk cx-anchor-setup's terms: ca theta order is
# L (column-major), P (lower triangle), G, K = (l1..l6, psi12, g1, g2, kappa3);
# this file's (and the spec's) is (l1..l6, g1, g2, kappa3, psi12): to_ca maps one to the other.
Lfit <- matrix(0, 6, 2); Lfit[1:3, 1] <- 0.6; Lfit[4:6, 2] <- 0.6
MOD <- ca_mod(Lfit, Pfree = matrix(c(FALSE, TRUE, FALSE, FALSE), 2),
              gam = c(0, 0), gfree = c(TRUE, TRUE), kap = rep(0, 6),
              kfree = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
stopifnot(identical(MOD$kind, c(rep("L", 6), "P", "G", "G", "K")),
          identical(MOD$i, c(1:6, 2L, 1L, 2L, 3L)),
          identical(MOD$j, c(1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 1L)))
to_ca <- c(1:6, 10L, 7L, 8L, 9L)
cx_sig_th <- function(th) c(th[1:nt], ca_s(MOD, th[nt + to_ca]))
cx_delta <- function(th) {
  Dca <- ca_delta(MOD, th[nt + to_ca])
  D <- matrix(0, nt + nrow(Dca), length(th))
  D[1:nt, 1:nt] <- diag(nt)
  D[nt + seq_len(nrow(Dca)), nt + to_ca] <- Dca
  D
}

# the reference's own start (fixed rule, written before the first run)
ref_start <- function(s) {
  b <- s[nt + 1:cx_p]
  R <- diag(cx_p); R[cx_pairs] <- s[nt + cx_p + seq_len(nrow(cx_pairs))]; R <- R + t(R) - diag(cx_p)
  tri <- function(i, j, k) sqrt(min(max(abs(R[i, j] * R[i, k] / R[j, k]), 0.01), 0.95^2))
  lam <- c(tri(1, 2, 3), tri(2, 1, 3), tri(3, 1, 2), tri(4, 5, 6), tri(5, 4, 6), tri(6, 4, 5))
  psi12 <- min(max(mean(R[1:3, 4:6] / outer(lam[1:3], lam[4:6])), -0.9), 0.9)
  c(s[1:nt], lam, mean((b / lam)[1:2]), mean((b / lam)[4:6]), 0, psi12)
}

# ---- one dataset -------------------------------------------------------------
ref_one <- function(dat) {
  stopifnot(identical(names(dat), c(sprintf("y%d", 1:6), "x")), nrow(dat) > 0, !anyNA(dat))
  for (j in seq_len(cx_p)) stopifnot(all(tabulate(dat[[j]] + 1L, cx_K[j]) > 0))
  # 1. weighted patterns (as 99n)
  key <- do.call(paste, dat); tab <- table(key); first <- match(names(tab), key)
  d <- list(Y = as.matrix(dat[first, 1:6]), x = dat$x[first], w = as.numeric(tab) / N)
  dimnames(d$Y) <- NULL
  # 2. P1 refusal check (as 99n)
  for (j in seq_len(cx_p)) {
    r0 <- range(d$Y[d$x == 0, j]); r1 <- range(d$Y[d$x == 1, j])
    if (min(r0[2], r1[2]) - max(r0[1], r1[1]) + 1 <= 1) stop(sprintf("P1 refusal on item %d", j), call. = FALSE)
  }
  # 3. Stage 1
  f <- cx_stage1(d)
  stopifnot(all(is.finite(f$s)), length(f$s) == length(s_lab), max(abs(f$rho)) < 0.9)
  # 4. psi from the analytic blocks; Gamma-hat (Hessian form)
  psi <- cxh_psi(d, f)
  Gh <- cx_gamma(d, psi)
  stopifnot(all(is.finite(Gh)), !inherits(tryCatch(chol(Gh), error = function(e) e), "error"))
  # 5. W_D
  W <- diag(1 / diag(Gh))
  # 6. Gauss-Newton from the reference's own start
  st <- ref_start(f$s)
  g <- cxs_gn(f$s, W, st)
  if (!g$ok) stop(sprintf("Gauss-Newton did not converge (%s iterations)", g$it), call. = FALSE)
  th <- g$th
  # 6a. convergence, demonstrated and asserted (as 99n step 6a)
  item_sc <- vapply(seq_len(cx_p), function(j)
    max(abs(colSums(d$w * cx_uni_score(d$Y[, j], d$x, f$th[[j]])))), 0)
  pair_cv <- t(vapply(seq_len(nrow(cx_pairs)), function(m) {
    i <- cx_pairs[m, 1]; j <- cx_pairs[m, 2]
    cells <- cx_pair_cells(d$Y[, i], d$Y[, j], d$x, d$w)
    S <- function(rr) sum(cells$w * cx_pair_u(cells, f$th[[i]], f$th[[j]], rr))
    sc <- S(f$rho[m]); sl <- (S(f$rho[m] + 1e-6) - S(f$rho[m] - 1e-6)) / 2e-6
    c(score = abs(sc), slope = sl, step = abs(sc / sl))
  }, c(score = 0, slope = 0, step = 0)))
  D <- cx_delta(th); r <- f$s - cx_sig_th(th)
  grad <- drop(t(D) %*% W %*% r)
  step_next <- drop(solve(t(D) %*% W %*% D, grad))
  conv <- c(s1_item_score_max = max(item_sc), s1_pair_score_max = max(pair_cv[, "score"]),
            s1_pair_step_max = max(pair_cv[, "step"]), s1_pair_slope_max = max(pair_cv[, "slope"]),
            s1_rho_abs_max = max(abs(f$rho)), gn_next_step = max(abs(step_next)),
            gn_grad_max = max(abs(grad)))
  if (!(conv[["s1_item_score_max"]] <= CONV_S1 && conv[["s1_pair_step_max"]] <= CONV_S1 &&
        conv[["s1_pair_slope_max"]] < 0 && conv[["s1_rho_abs_max"]] < 0.99 &&
        conv[["gn_next_step"]] < GN_TOL))
    stop(sprintf("convergence not demonstrated: %s",
                 paste(names(conv), sprintf("%.3e", conv), collapse = ", ")), call. = FALSE)
  # the finite-difference Delta against a coarser step (rounding only is expected)
  D4 <- local({ h <- 1e-4; ca_delta(MOD, th[nt + to_ca], h = h) })
  fd_gap <- max(abs(D4 - D[nt + seq_len(nrow(D4)), nt + to_ca]))
  stopifnot(fd_gap < 1e-7)
  # 7. sandwich and T3 (99n's lines, verbatim)
  nth <- length(th)
  B <- solve(t(D) %*% W %*% D)
  sand <- function(G) sqrt(diag(B %*% t(D) %*% W %*% G %*% W %*% D %*% B) / (N - 1))
  U <- W - W %*% D %*% B %*% t(D) %*% W; UG <- U %*% Gh; df <- length(r) - nth
  a <- sqrt(df / sum(diag(UG %*% UG))); b <- df - a * sum(diag(UG))
  Fmin <- drop(t(r) %*% W %*% r)
  T <- (N - 1) * Fmin
  se <- sand(Gh); T3 <- a * T + b
  # design floor: the misfit is real
  stopifnot(T >= 2 * df, nth == npar)
  list(N = N, npat = length(d$w), s = f$s, G = Gh, th = th, se = se, F = Fmin, T = T, a = a, b = b,
       T3 = T3, p = pchisq(T3, df, lower.tail = FALSE), df = df, it = g$it, conv = conv,
       fd_gap = fd_gap, rho_max = max(abs(f$rho)))
}

res <- ref_one(dat)

# ---- writing (as 99n's write_one) --------------------------------------------
write_csv_lf(dat, DAT)
out <- file.path(OUTROOT, "expected", sprintf("condx_ref_%s_%s.csv", TAG, c("s", "gamma", "theta", "scalars")))
names(out) <- c("s", "gamma", "theta", "scalars")
write_csv_lf(data.frame(label = s_lab, value = f17(res$s)), out[["s"]])
ij <- which(lower.tri(res$G, diag = TRUE), arr.ind = TRUE)
ij <- ij[order(ij[, "col"], ij[, "row"]), , drop = FALSE]
write_csv_lf(data.frame(row = s_lab[ij[, "row"]], col = s_lab[ij[, "col"]], value = f17(res$G[ij])), out[["gamma"]])
write_csv_lf(data.frame(label = th_lab, est = f17(res$th), se = f17(res$se)), out[["theta"]])
write_csv_lf(data.frame(
  quantity = c("N", "npar", "df", "F", "T_raw", "a", "b", "T3", "p_T3", "gn_iter", names(res$conv)),
  value = c(as.character(res$N), as.character(npar), as.character(res$df),
            f17(c(res$F, res$T, res$a, res$b, res$T3, res$p)), as.character(res$it), f17(res$conv))),
  out[["scalars"]])

MAN <- file.path(OUTROOT, sprintf("MANIFEST_condx_ref_%s.json", TAG))
write_manifest(MAN, c(DAT, out),
  seeds = sprintf("set.seed(%dL) for condx_U.csv; RNGkind %s; order: rbinom(N, 1, 0.4), rnorm(2N) for the disturbances, then per item rnorm(N) residuals",
                  SEED, paste(RNGkind(), collapse = "/")),
  extra = list(packages = c(
    "none (base R only)", R.version.string,
    "development fixture: NOT a preregistered dataset (links 1 to 3 use condx_A, condx_D, condx_link3)",
    "design: N 600, K 3,2,4,3,2,3; two factors F1 y1-y3 and F2 y4-y6 plus cross-loading y4 <- F1 0.20 (lambda4 0.65); psi12 0.30; gamma 0.5, -0.4; kappa y3 0.4; P(x=1) 0.4; fitted model omits the cross-loading",
    paste("qmd sha256", QMD, tools::sha256sum(QMD)),
    paste("qmd chunks run whole, in order:", paste(CHUNKS, collapse = ", ")),
    paste("run-script md5", RUN, tools::md5sum(RUN)),
    paste("run-script sha256", RUN, tools::sha256sum(RUN)),
    paste("run-script definitions used:", paste(RUN_DEFS, collapse = ", ")),
    paste("common md5", COMMON, tools::md5sum(COMMON), "(write_csv_lf, write_manifest, labels; its Phi2/phi2 removed)"),
    "Stage 1, psi, Gamma-hat: cx_stage1, cxh_psi, cx_gamma with cx_K, cx_p, cx_pairs reassigned for six items; W_D = diag(1 / diag(Gamma-hat)); SE and T divisor N - 1",
    "model map: ca_mod / ca_s / ca_delta (chunk cx-anchor-setup), thresholds prepended with an identity block; Delta by central differences h 1e-6",
    sprintf("Gauss-Newton: cxs_gn, start ref_start(s) (this file), stop max|step| < %.0e, at most %d iterations", GN_TOL, GN_MAX),
    sprintf("diagnostics: Delta central-difference step 1e-6 vs 1e-4 gap %.3e; max |rho-hat| %.4f", res$fd_gap, res$rho_max),
    "Phi2: the .qmd's integrate() form, rel.tol 1e-12; phi2: cx_phi2 (plain form)",
    sprintf("condx_U: N %d, %d patterns, GN %d iterations, F %.6g, T_raw %.6f, df %d, T3 %.6f, p %.4g",
            res$N, res$npat, res$it, res$F, res$T, res$df, res$T3, res$p))))
cat(sprintf("condx_U: N %d, %d patterns, GN %d iterations, F %.6g, T_raw %.6f, df %d, T3 %.6f, p %.4g\n",
            res$N, res$npat, res$it, res$F, res$T, res$df, res$T3, res$p))
cat(sprintf("%d files and %s written; %.1f s\n", length(out) + 1L, MAN, proc.time()[["elapsed"]] - t_start))
