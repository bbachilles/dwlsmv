# ---------------------------------------------------------------------------
# EXECUTED 2026-11-10 by the main session from the repository root (R 4.6.0;
# log docs/notes/condx-validation/run_99n.log), on condx_A and condx_D.
# Drafted by a read-only agent; the convergence assertions (step 6a) were
# written by the main session.
#
# dwlsmv fixture generation -- tier 1: the analytic R reference for link 2 of
# the conditional-mode numerical validation (docs/notes/condx-validation/
# prereg.md, sections 2 and 8 [V3]). One binary x, one group, items y1, y2
# (3 categories), y3 (binary); model f =~ y1 + y2 + y3; f ~ x; psi = 1.
#
# It defines no mathematics. Every quantity comes from the derivation's own
# functions, extracted from the CURRENT spec/derivations/dwlsmv-derivations.qmd
# (SHA-256 and chunk labels in the manifest), and from the frozen simulation's
# run script docs/notes/condx-sim/condx_sim_run.R (definitions only, its MD5 in
# the manifest). The pipeline is condx_sim_run.R's cxs_rep(), with one change:
# psi comes from the analytic Hessian-form blocks (cxh_psi), not from the
# finite-difference Hessian (cx_psi_n(..., "hessian")).
#
#   1. data -> weighted patterns (w = n_k / N), as chunk cxn-lavaan builds dA
#   2. P1 refusal check, as cxs_rep (a refusal stops the generator)
#   3. Stage 1: cx_stage1 (chunk cx-stage1)
#   4. psi: cxh_psi (chunk cxh-forms); Gamma-hat = cx_gamma (chunk cxn-psi),
#      N-divided, sum_k w_k psi_k psi_k'; positive definite (chol), as cxs_rep
#   5. W_D = diag(1 / diag(Gamma-hat)), as cxs_rep
#   6. DWLS by cxs_gn (Gauss-Newton; stop when max|step| < GN_TOL = 1e-10,
#      at most GN_MAX = 200) from cxs_start(s), both from condx_sim_run.R
#   7. sandwich SEs with the N - 1 divisor; T = (N - 1) r'W r; a, b from
#      U Gamma; T3 = a T + b; df = length(s) - length(theta): the lines of
#      cxs_rep, copied below verbatim (cxs_rep itself resamples a population
#      and cannot be called on a dataset)
#   8. blocks H_j (cxh_H), U_rho and U_theta (cxh_pair_d, summed as cxh_psi
#      sums them), written as diagnostics; psi_rho rebuilt from the written
#      blocks must equal cxh_psi's column exactly
#
# Outputs, per dataset <tag> (condx_<tag>.csv), under OUTROOT/expected/:
#   condx_ref_<tag>_s.csv       label, value         s in spec 4.2 order
#   condx_ref_<tag>_gamma.csv   row, col, value      Gamma-hat (Hessian form,
#                                                   N-divided), lower triangle
#                                                   with diagonal, column-major
#   condx_ref_<tag>_blocks.csv  block, row, col, value
#                               block H_y<j> (K_j x K_j, full), Urho_y<i>~~y<j>
#                               (1 x 1), Utheta_y<i>~~y<j> (1 x (K_i + K_j);
#                               row = the pair label, col = the theta_i, theta_j
#                               labels). The blocks are the weighted sums
#                               (w = n_k / N), at s-hat.
#   condx_ref_<tag>_theta.csv   label, est, se       theta = (tau, lambda, gamma)
#   condx_ref_<tag>_scalars.csv quantity, value      N, T_raw, a, b, T3, df,
#                                                   gn_iter, gn_last_step
# and the manifest OUTROOT/MANIFEST_condx_ref.json. Numbers %.17g.
#
# Run from the repository root:
#   Rscript fixtures/generate/99n_condx_ref.R            # OUTROOT = fixtures
#   Rscript fixtures/generate/99n_condx_ref.R <outroot>  # elsewhere (testing)
# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

DATASETS <- file.path("fixtures", "data", c("condx_A.csv", "condx_D.csv"))
args <- commandArgs(TRUE)
OUTROOT <- if (length(args) >= 1) args[1] else "fixtures"
QMD <- file.path("spec", "derivations", "dwlsmv-derivations.qmd")
RUN <- file.path("docs", "notes", "condx-sim", "condx_sim_run.R")
COMMON <- file.path("fixtures", "generate", "00_common.R")
CHUNKS <- c("cx-setup", "cx-stage1", "cx-model", "cxn-psi", "cxh-forms")
RUN_DEFS <- c("GN_TOL", "GN_MAX", "cxs_gn", "cxs_start")
CONV_S1 <- 1e-12   # Stage 1 convergence limit (step 6a); chosen after measuring condx_A only

missing_in <- DATASETS[!file.exists(DATASETS)]
if (length(missing_in))
  stop("99n_condx_ref.R: input dataset(s) missing, nothing written: ",
       paste(missing_in, collapse = ", "), call. = FALSE)
stopifnot(dir.exists(file.path(OUTROOT, "expected")), file.exists(QMD), file.exists(RUN))

t_start <- proc.time()[["elapsed"]]

# write_csv_lf, write_manifest and the label helpers. 00_common.R also defines
# Phi2 and phi2; they are removed so that the .qmd's own definitions are used
# (cx-setup defines Phi2 only if none exists; the two texts are the same).
source(COMMON)
rm(Phi2, phi2)

# ---- the derivation's functions ---------------------------------------------
qmd <- readLines(QMD, encoding = "UTF-8")
chunk <- function(lab) {
  s <- grep(paste0("^```\\{r ", lab, "\\}"), qmd); stopifnot(length(s) == 1)
  e <- s + which(qmd[(s + 1):length(qmd)] == "```")[1]
  qmd[(s + 1):(e - 1)]
}
for (lab in CHUNKS) eval(parse(text = chunk(lab)), globalenv())
stopifnot(exists("cx_stage1"), exists("cxh_psi"), exists("cxh_H"), exists("cxh_pair_d"),
          exists("cx_gamma"), exists("cx_delta"), exists("cx_sig_th"), exists("Phi2"),
          identical(cx_K, c(3L, 3L, 2L)), cx_p == 3L)

# condx_sim_run.R: the named top-level assignments only (its population,
# theta* and replicate loop are not run)
for (ex in parse(RUN, keep.source = FALSE)) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) &&
      as.character(ex[[2]]) %in% RUN_DEFS)
    eval(ex, globalenv())
}
stopifnot(all(vapply(RUN_DEFS, exists, NA)), GN_TOL == 1e-10, GN_MAX == 200L)

# ---- labels -------------------------------------------------------------------
nt <- sum(cx_K - 1)
s_lab <- c(threshold_labels(cx_K - 1), sprintf("y%d~x", seq_len(cx_p)),
           correlation_labels(cx_p))
stopifnot(identical(correlation_labels(cx_p),
                    sprintf("y%d~~y%d", cx_pairs[, 1], cx_pairs[, 2])))
th_lab <- c(threshold_labels(cx_K - 1), sprintf("f=~y%d", seq_len(cx_p)), "f~x")
item_lab <- function(j) c(sprintf("y%d|t%d", j, seq_len(cx_K[j] - 1)), sprintf("y%d~x", j))
f17 <- function(x) sprintf("%.17g", x)

# ---- one dataset --------------------------------------------------------------
ref_one <- function(path) {
  raw <- read.csv(path, colClasses = "character")
  stopifnot(identical(names(raw), c("y1", "y2", "y3", "x")), nrow(raw) > 0)
  dat <- raw; dat[] <- lapply(raw, as.integer)
  stopifnot(!anyNA(dat), all(dat$x %in% 0:1),
            all(vapply(seq_len(cx_p), function(j) all(dat[[j]] %in% 0:(cx_K[j] - 1)), NA)))
  # every pooled category observed (cx_uni_fit's start needs it)
  for (j in seq_len(cx_p)) stopifnot(all(tabulate(dat[[j]] + 1L, cx_K[j]) > 0))
  N <- nrow(dat)
  # 1. weighted patterns, as chunk cxn-lavaan builds dA
  key <- do.call(paste, dat); tab <- table(key); first <- match(names(tab), key)
  d <- list(Y = as.matrix(dat[first, 1:3]), x = dat$x[first], w = as.numeric(tab) / N)
  dimnames(d$Y) <- NULL
  # 2. P1, as cxs_rep
  for (j in seq_len(cx_p)) {
    r0 <- range(d$Y[d$x == 0, j]); r1 <- range(d$Y[d$x == 1, j])
    if (min(r0[2], r1[2]) - max(r0[1], r1[1]) + 1 <= 1)
      stop(sprintf("%s: P1 refusal on item %d", path, j), call. = FALSE)
  }
  # 3. Stage 1
  f <- cx_stage1(d)
  stopifnot(all(is.finite(f$s)), length(f$s) == length(s_lab))
  # 4. psi from the analytic blocks; Gamma-hat (Hessian form)
  psi <- cxh_psi(d, f)
  Gh <- cx_gamma(d, psi)
  stopifnot(all(is.finite(Gh)), !inherits(tryCatch(chol(Gh), error = function(e) e), "error"))
  # 5. W_D
  W <- diag(1 / diag(Gh))
  # 6. Gauss-Newton from the data-based start
  st <- cxs_start(f$s)
  g <- cxs_gn(f$s, W, st)
  if (!g$ok) stop(sprintf("%s: Gauss-Newton did not converge (%s iterations)", path, g$it),
                  call. = FALSE)
  # the last step's size: re-run to it - 1 iterations (deterministic) and take
  # the same step expression as cxs_gn at that iterate
  GN_MAX_frozen <- GN_MAX
  GN_MAX <<- g$it - 1L
  th_prev <- if (g$it > 1L) cxs_gn(f$s, W, st)$th else st
  GN_MAX <<- GN_MAX_frozen
  Dp <- cx_delta(th_prev)
  step_last <- drop(solve(t(Dp) %*% W %*% Dp, t(Dp) %*% W %*% (f$s - cx_sig_th(th_prev))))
  stopifnot(identical(th_prev + step_last, g$th), max(abs(step_last)) < GN_TOL)
  # 6a. CONVERGENCE, demonstrated and asserted (main session, 2026-11-10; the
  # user's condition that the reference show convergence on every dataset and
  # not merely produce output). cx_uni_fit and cx_pair_fit return no flag, so
  # the residuals are measured here at the returned estimates:
  #   items: max |weighted score|; pairs: |score|, its slope by a central
  #   difference at 1e-6 (as cx_pair_fit's own Newton step), and the Newton
  #   step that slope implies; fit: the NEXT Gauss-Newton step at theta-hat and
  #   the gradient D'W(s - sigma).
  # Limits CONV_S1 = 1e-12 and GN_TOL, four decades inside the s tier (1e-8).
  # CONV_S1 was chosen after measuring condx_A (scores <= 2.5e-16, steps
  # <= 1.2e-16) and before any figure for condx_D existed.
  item_sc <- vapply(seq_len(cx_p), function(j)
    max(abs(colSums(d$w * cx_uni_score(d$Y[, j], d$x, f$th[[j]])))), 0)
  pair_cv <- t(vapply(seq_len(nrow(cx_pairs)), function(m) {
    i <- cx_pairs[m, 1]; j <- cx_pairs[m, 2]
    cells <- cx_pair_cells(d$Y[, i], d$Y[, j], d$x, d$w)
    S <- function(rr) sum(cells$w * cx_pair_u(cells, f$th[[i]], f$th[[j]], rr))
    sc <- S(f$rho[m]); sl <- (S(f$rho[m] + 1e-6) - S(f$rho[m] - 1e-6)) / 2e-6
    c(score = abs(sc), slope = sl, step = abs(sc / sl))
  }, c(score = 0, slope = 0, step = 0)))
  Dh <- cx_delta(g$th); rh <- f$s - cx_sig_th(g$th)
  grad <- drop(t(Dh) %*% W %*% rh)
  step_next <- drop(solve(t(Dh) %*% W %*% Dh, grad))
  conv <- c(s1_item_score_max = max(item_sc), s1_pair_score_max = max(pair_cv[, "score"]),
            s1_pair_step_max = max(pair_cv[, "step"]), s1_pair_slope_max = max(pair_cv[, "slope"]),
            s1_rho_abs_max = max(abs(f$rho)), gn_next_step = max(abs(step_next)),
            gn_grad_max = max(abs(grad)))
  if (!(conv[["s1_item_score_max"]] <= CONV_S1 && conv[["s1_pair_step_max"]] <= CONV_S1 &&
        conv[["s1_pair_slope_max"]] < 0 && conv[["s1_rho_abs_max"]] < 0.99 &&
        conv[["gn_next_step"]] < GN_TOL))
    stop(sprintf("%s: convergence not demonstrated: %s", path,
                 paste(names(conv), sprintf("%.3e", conv), collapse = ", ")), call. = FALSE)
  # 7. sandwich and T3: cxs_rep's lines, verbatim
  nth <- length(g$th)
  th <- g$th; D <- cx_delta(th); B <- solve(t(D) %*% W %*% D)
  sand <- function(G) sqrt(diag(B %*% t(D) %*% W %*% G %*% W %*% D %*% B) / (N - 1))
  r <- f$s - cx_sig_th(th)
  U <- W - W %*% D %*% B %*% t(D) %*% W; UG <- U %*% Gh; df <- length(r) - nth
  a <- sqrt(df / sum(diag(UG %*% UG))); b <- df - a * sum(diag(UG))
  T <- (N - 1) * drop(t(r) %*% W %*% r)
  se <- sand(Gh); T3 <- a * T + b
  # 8. the blocks, as cxh_psi forms them; psi_rho rebuilt from them must be
  # cxh_psi's column exactly
  H <- lapply(seq_len(cx_p), function(j) cxh_H(d$Y[, j], d$x, d$w, f$th[[j]]))
  psi1 <- lapply(seq_len(cx_p), function(j)
    cx_uni_score(d$Y[, j], d$x, f$th[[j]]) %*% t(solve(H[[j]])))
  pb <- lapply(seq_len(nrow(cx_pairs)), function(m) {
    i <- cx_pairs[m, 1]; j <- cx_pairs[m, 2]
    cells <- cx_pair_cells(d$Y[, i], d$Y[, j], d$x, d$w)
    pd <- cxh_pair_d(cells, f$th[[i]], f$th[[j]], f$rho[m])
    Ur <- sum(cells$w * pd$Ur); U1 <- colSums(cells$w * pd$Ut)
    u <- pd$u[match(cells$key, cells$ukey)]
    list(Ur = Ur, Ut = U1, psir = -(u + cbind(psi1[[i]], psi1[[j]]) %*% U1) / Ur)
  })
  stopifnot(identical(unname(psi[, nt + cx_p + seq_len(nrow(cx_pairs))]),
                      unname(sapply(pb, function(z) z$psir))))
  list(N = N, npat = length(d$w), s = f$s, G = Gh, H = H, pb = pb, th = th, se = se,
       T = T, a = a, b = b, T3 = T3, df = df, it = g$it, step = max(abs(step_last)), conv = conv,
       start = st)
}

# ---- writing ------------------------------------------------------------------
write_one <- function(res, tag) {
  out <- file.path(OUTROOT, "expected", sprintf("condx_ref_%s_%s.csv", tag,
                   c("s", "gamma", "blocks", "theta", "scalars")))
  names(out) <- c("s", "gamma", "blocks", "theta", "scalars")
  write_csv_lf(data.frame(label = s_lab, value = f17(res$s)), out[["s"]])
  ij <- which(lower.tri(res$G, diag = TRUE), arr.ind = TRUE)
  ij <- ij[order(ij[, "col"], ij[, "row"]), , drop = FALSE]
  write_csv_lf(data.frame(row = s_lab[ij[, "row"]], col = s_lab[ij[, "col"]],
                          value = f17(res$G[ij])), out[["gamma"]])
  bl <- list()
  for (j in seq_len(cx_p)) {
    lab <- item_lab(j); e <- expand.grid(r = seq_along(lab), c = seq_along(lab))
    bl[[length(bl) + 1]] <- data.frame(block = sprintf("H_y%d", j), row = lab[e$r],
                                       col = lab[e$c], value = f17(res$H[[j]][cbind(e$r, e$c)]))
  }
  for (m in seq_len(nrow(cx_pairs))) {
    i <- cx_pairs[m, 1]; j <- cx_pairs[m, 2]; pl <- sprintf("y%d~~y%d", i, j)
    bl[[length(bl) + 1]] <- data.frame(block = paste0("Urho_", pl), row = pl, col = pl,
                                       value = f17(res$pb[[m]]$Ur))
    bl[[length(bl) + 1]] <- data.frame(block = paste0("Utheta_", pl), row = pl,
                                       col = c(item_lab(i), item_lab(j)),
                                       value = f17(res$pb[[m]]$Ut))
  }
  write_csv_lf(do.call(rbind, bl), out[["blocks"]])
  write_csv_lf(data.frame(label = th_lab, est = f17(res$th), se = f17(res$se)), out[["theta"]])
  write_csv_lf(data.frame(
    quantity = c("N", "T_raw", "a", "b", "T3", "df", "gn_iter", "gn_last_step", names(res$conv)),
    value = c(as.character(res$N), f17(c(res$T, res$a, res$b, res$T3)),
              as.character(res$df), as.character(res$it), f17(res$step), f17(res$conv))), out[["scalars"]])
  out
}

files <- character(); notes <- character()
for (path in DATASETS) {
  tag <- sub("^condx_(.*)\\.csv$", "\\1", basename(path))
  stopifnot(nzchar(tag), tag != basename(path))
  t0 <- proc.time()[["elapsed"]]
  res <- ref_one(path)
  files <- c(files, path, write_one(res, tag))
  notes <- c(notes, sprintf(
    "%s: N %d, %d patterns, GN %d iterations, last step %.3e, next step %.3e, item score %.3e, pair step %.3e, T3 %.6f, %.1f s",
    tag, res$N, res$npat, res$it, res$step, res$conv[["gn_next_step"]], res$conv[["s1_item_score_max"]], res$conv[["s1_pair_step_max"]], res$T3, proc.time()[["elapsed"]] - t0))
  cat(tail(notes, 1), "\n")
}

MAN <- file.path(OUTROOT, "MANIFEST_condx_ref.json")
write_manifest(MAN, files,
  seeds = "none (deterministic; each dataset's seed is in its own generator's manifest)",
  extra = list(packages = c(
    "none (base R only)", R.version.string,
    paste("qmd sha256", QMD, tools::sha256sum(QMD)),
    paste("qmd chunks run whole, in order:", paste(CHUNKS, collapse = ", ")),
    paste("run-script md5", RUN, tools::md5sum(RUN)),
    paste("run-script sha256", RUN, tools::sha256sum(RUN)),
    paste("run-script definitions used:", paste(RUN_DEFS, collapse = ", ")),
    paste("common md5", COMMON, tools::md5sum(COMMON), "(write_csv_lf, write_manifest, labels; its Phi2/phi2 removed)"),
    "pipeline: condx_sim_run.R cxs_rep with psi = cxh_psi (analytic Hessian-form blocks) in place of cx_psi_n(, 'hessian')",
    "Gamma-hat: cx_gamma, N-divided (w = n_k / N); W_D = diag(1 / diag(Gamma-hat)); SE and T divisor N - 1",
    sprintf("Gauss-Newton: cxs_gn, start cxs_start(s), stop max|step| < %.0e, at most %d iterations", GN_TOL, GN_MAX),
    "Phi2: the .qmd's integrate() form, rel.tol 1e-12; phi2: cx_phi2 (plain form)",
    "univariate P: Stage 1 and the score g use Phi(up) - Phi(lo) (cx_uni_score); H_j uses cxh_P's upper-tail branch (prereg [V4])",
    "pair corner order in cxh_pair_d: (a1,b1), (a1,b2), (a2,b1), (a2,b2) (prereg phi2 clause)",
    notes)))
cat(sprintf("%d files and %s written; %.1f s\n", length(files), MAN,
            proc.time()[["elapsed"]] - t_start))
